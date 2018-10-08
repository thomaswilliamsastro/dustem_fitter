
MODULE MGET_TDIST

  USE CONSTANTS
  IMPLICIT NONE

  REAL (KIND=dp), PUBLIC, ALLOCATABLE  :: u_cal(:), lu_cal(:)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE  :: cal(:), lcal(:)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE  :: t_cal(:)

  PRIVATE
  PUBLIC :: GET_TDIST

CONTAINS

!----------------------------------------------------------------------------

SUBROUTINE GET_TDIST(nt, ns, a, enerabs, udist, p, t, p1, hcap, tmoy, tequi, tmax)
!========================================================================
!    GET_TDIST computes the t distribution of interstellar grains heated
!    by the radiation field ISRF.DAT following Desert et al. 1986
!
!    HEATING
!    - fully discrete (state-to-state U) from radiation field
!
!    COOLING
!    - fully continuous (Qabs*Bnu(T))
!
!    t derived from u with C(T), tmax corresponds to 2*hnumax
!    Equilibrium temperature tequi is found from cooling = heating
!
!    dP/dU derived as in Desert et al (1986)
!    (formulae 24 and 25) then dP/dT=C(T)*dP/dU
!
!========================================================================

  ! modules for global variables
  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  !  Arguments
  INTEGER,        INTENT (IN)    :: nt, ns
  REAL (KIND=dp), INTENT (IN)    :: a 
  REAL (KIND=dp), INTENT (OUT)   :: enerabs
  REAL (KIND=dp), INTENT (OUT)   :: udist(ndist)
  REAL (KIND=dp), INTENT (OUT)   :: tmoy
  REAL (KIND=dp), INTENT (OUT)   :: tequi
  REAL (KIND=dp), INTENT (OUT)   :: tmax
  REAL (KIND=dp), INTENT (OUT)   :: p(ndist)
  REAL (KIND=dp), INTENT (OUT)   :: p1(ndist)
  REAL (KIND=dp), INTENT (OUT)   :: t(ndist)
  REAL (KIND=dp), INTENT (OUT)   :: hcap(ndist)

  !  Local parameters
  REAL (KIND=dp)               :: aux, aux1
  REAL (KIND=dp)               :: e1, faux

  REAL (KIND=dp)               :: intchauf                ! total heating by photons
  INTEGER                      :: i1, i2
  INTEGER                      :: nequi                   ! index of equilibrium temperature tequi
  REAL (KIND=dp)               :: du, l0
  REAL (KIND=dp)               :: surf
  REAL (KIND=dp)               :: tlow, tmin

  REAL (KIND=dp), ALLOCATABLE  :: t2(:)
  REAL (KIND=dp)               :: nbrpho(n_qabs)
  REAL (KIND=dp)               :: nuinusigma(n_qabs), intnuinusigma(n_qabs)

  REAL (KIND=dp)               :: g(ndist)
  REAL (KIND=dp)               :: lij(ndist,ndist)

  ! JLB
  INTEGER                      :: ntmp
  REAL (KIND=dp)               :: uequi
  REAL (KIND=dp)               :: tau, t_aux1, t_aux2
  REAL (KIND=dp), ALLOCATABLE  :: cdist(:)
  INTEGER                      :: jlo

  !**********************************************************************

  ALLOCATE (fdist(1:ndist,1:ndist))
  ALLOCATE (t2(1:ndist))
  ALLOCATE (cdist(1:ndist))

  ! inits
  tlow = MINVAL(10.0_dp**temp_cal(nt,1:n_temp(nt))) ! lowest value (from C_TYPE.DAT) for tmin, tequi search
  IF (tlow > t_cmb) THEN
     tlow = t_cmb
  ENDIF
  surf = xpi * a**2

  ! get grain Qabs
  ALLOCATE (qauv(n_qabs))
  qauv(:) = qi_abs(nt,ns,:)

  ! VG : because of possible anistropic radiation field (anistropic heating)
  ! We distinguish between emission cross-section for emission (qaem, always isotropic) 
  ! and absorption cross-sections for anisotropic heating (qauv)
  ALLOCATE (qaem(n_qabs))

 ! Emission is always isotropic
  qaem(:) = qi_abs(nt,ns,:)

  IF (anisG0 .gt. 0) THEN
     ! Grains are aligned and radiation field is anisotropic with degree anisG0
     qauv(:) = anisg0 * ( f_pol(nt,ns) * qiH_abs(nt,ns,:)  &
          + (1-f_pol(nt,ns)) * qi_abs(nt,ns,:) ) &
          + (1-anisg0) * qi_abs(nt,ns,:)
   ELSE
      ! Radiation field is isotropic
      qauv(:) = qi_abs(nt,ns,:)
   ENDIF

  ! Energy absorbed per sec
  nuinusigma(1:n_qabs) = qauv(1:n_qabs) * freq_qabs(1:n_qabs) * isrfuv(1:n_qabs)
  CALL PRIMITIV2 (1, n_qabs, lfrq_qabs, nuinusigma, intnuinusigma)
  nbrpho(1:n_qabs) = qauv(1:n_qabs) * isrfuv(1:n_qabs) / xhp

  ! total heating
  intchauf = intnuinusigma(n_qabs)
  enerabs = surf * intchauf

  ! get tequi
  t_aux1 = tlow
  t_aux2 = 3000.0_dp

  CALL GET_TEQUIL (nt, ns, t_aux1, t_aux2, intchauf, a)
  tequi = t_aux2

  ntmp = n_temp(nt)
  ALLOCATE (u_cal(ntmp))
  ALLOCATE (lu_cal(ntmp))
  ALLOCATE (cal(ntmp))
  ALLOCATE (lcal(ntmp))
  ALLOCATE (t_cal(ntmp))
  CALL GET_UC (nt, a, ntmp)
  ! get continuous cooling
  CALL GET_TC (nt, a, tequi, uequi, nequi, tmin, udist(1:ndist), t(1:ndist), hcap(1:ndist))  ! returns nequi
  tmax = t(ndist)
  CALL COOLING (nt, ns, t, g, ndist, a)
  g(:) = g(:) * surf
  ! get photon energies from internal energy and J(U1,U2) (lij)
  lij(:,:) = 0.0_dp
  DO i1=1,ndist-1
     jlo = 1
     DO i2=i1+1,ndist
        e1 = udist(i2) - udist(i1)
        IF (e1 >= istiny) lij(i1,i2) = INTPOL3 (nbrpho, freq_qabs, n_qabs, e1/xhp, jlo) / e1
     ENDDO
  ENDDO
  lij(:,:) = lij(:,:) * surf

  ! get kernel fdist: needs full i2 (not only i2>i1)
  fdist(:,:) = 0.0_dp
  l0 = XINTEG2 (1, ndist, ndist, udist(1:ndist), lij(1:ndist,ndist) ) ! same as below but faster
  DO i1 = 1,ndist
     fdist(i1,ndist) = 0.0_dp
     DO i2 = ndist-1,1,-1
        du = udist(i2+1) - udist(i2)
        IF (g(i2) <= istiny) THEN
           faux = 0.0_dp
        ELSE
           tau = l0*du / g(i2)
           IF (tau <= -LOG10(istiny)) THEN 
              faux = (lij(i1,i2) * du + fdist(i1,i2+1) * g(i2+1)) * EXP(-tau) / g(i2)
           ELSE
              faux = lij(i1,i2)*du / g(i2)
           ENDIF
           IF (faux < istiny) faux = 0.0_dp
        ENDIF
        fdist(i1,i2) = faux
     ENDDO
  ENDDO

  p(:) = 0.0_dp
  p(nequi) = 1.0_dp
  IF (MAXVAL(fdist) >= istiny) THEN 
     fdist = TRANSPOSE(fdist)
     nit = 80
     DO i1=1,nit
        p = MATMUL(fdist,p)
        aux = XINTEG2 (1, ndist, ndist, LOG(udist(1:ndist)), udist(1:ndist)*p(1:ndist))
        p = p / aux
!        p = p / SQRT(DOT_PRODUCT(p,p))
     ENDDO
  ELSE ! case of Dirac type distribution, add 2nd point for integration 
     p(nequi+1) = 1.0_dp
     udist(nequi+1) = udist(nequi) * (1.0_dp + 1.0e-10_dp)
     t(nequi+1) = t(nequi) * (1.0_dp + 1.0e-10_dp)
  ENDIF
  
  ! define arrays for output and SED integration
  ! NB using definition dP/dT = C(T)*dP/dU
  p1(1:ndist) = p(1:ndist) * t(1:ndist) * hcap(1:ndist)
  t2(1:ndist) = LOG(t(1:ndist))

  ! then normalize p, p1 (2nd T-loop but fast enough)
  ! could be removed because energy renormalized
  ! SED * (enerabs/enerem) in DM_compute.f90
  aux = XINTEG2 (1, ndist, ndist, LOG(udist(1:ndist)), udist(1:ndist)*p(1:ndist))
  aux1 = XINTEG2 (1, ndist, ndist, t2(1:ndist), p1(1:ndist))
  IF( (aux*aux1) <= tiniest) THEN
     WRITE(6,*)   ' (F) GET_TDIST: ', gtype(nt),' a(nm)= ',a*1.e7_dp
     WRITE(6,'(1P,2(1X,A,E9.2))') '   pb with sum of p= ',aux,' - p1= ',aux1
     STOP
  ELSE
     p(1:ndist) = p(1:ndist) / aux
     p1(1:ndist) = p1(1:ndist) / aux1
  ENDIF

  ! mean temperature
  tmoy = XINTEG2 (1, ndist, ndist, t2(1:ndist), t(1:ndist)*p1(1:ndist))

  DEALLOCATE (fdist)
  DEALLOCATE (t2)
  DEALLOCATE (qauv)
  DEALLOCATE (qaem)
  DEALLOCATE (cdist)
  DEALLOCATE (u_cal)
  DEALLOCATE (lu_cal)
  DEALLOCATE (cal)
  DEALLOCATE (lcal)
  DEALLOCATE (t_cal)

END SUBROUTINE GET_TDIST

!---------------------------------------------------------------------

SUBROUTINE GET_TEQUIL (nt, ns, ta, tb, c, a)
! finds equilibrium temperature tb

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  INTEGER, INTENT (IN)           :: nt               ! index of grain type
  INTEGER, INTENT (IN)           :: ns               ! index of grain size
  REAL (KIND=dp), INTENT (INOUT) :: ta, tb
  REAL (KIND=dp), INTENT (IN)    :: c
  REAL (KIND=dp), INTENT (IN)    :: a                ! grain size

  REAL (KIND=dp)                 :: fa, fb
  REAL (KIND=dp)                 :: f, t
  INTEGER                        :: i

  fa = COOL1(ta, nt, ns, a) - c
  fb = COOL1(tb, nt, ns, a) - c

  i = 0
  IF (fa*fb > 0.0_dp) THEN
     PRINT *, "  (W) GET_TEQUIL:  Wrong initial guess"
     PRINT *, "      Ta = ", ta, "  fa = ", fa
     PRINT *, "      Tb = ", tb, "  fb = ", fb
  ENDIF
  DO
     t = 0.5_dp * (ta + tb)
     f = COOL1(t, nt, ns, a) - c
     IF (f*fa > 0.0_dp) THEN
       ta = t
       fa = f
     ELSE
       tb = t
       fb = f
     ENDIF
     i = i + 1
     IF ((tb-ta) < 1.0e-0_dp) EXIT
  ENDDO

  CALL GET_TEQUIL2 (nt, ns, tb, c, a)                    !  Add one pass of Newton to get more digits

END SUBROUTINE GET_TEQUIL

!---------------------------------------------------------------------

SUBROUTINE GET_TEQUIL2 (nt, ns, ta, c, a)
! Solution of COOL(ta) = c by Newton scheme.

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  INTEGER, INTENT (IN)           :: nt               ! index of grain type
  REAL (KIND=dp), INTENT (INOUT) :: ta
  REAL (KIND=dp), INTENT (IN)    :: c
  REAL (KIND=dp), INTENT (IN)    :: a                ! grain size
  INTEGER, INTENT (IN)           :: ns               ! index of grain size

  REAL (KIND=dp)                 :: dta
  INTEGER                        :: i

  i = 0
  DO
     dta = (COOL1(ta,nt,ns,a) -c) / DCOOL1(ta,nt,ns,a)
     ta = ta - dta
     i = i + 1
     IF (ABS(dta) < 1.0e-10_dp) EXIT
  ENDDO

END SUBROUTINE GET_TEQUIL2

!----------------------------------------------------------------------------

SUBROUTINE GET_TC (nt, a, tequi, uequi, nequi, tmin, udist, tdist, cdist)
! gets T-grid for dP/dT and C(T)

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt
  REAL (KIND=dp), INTENT (IN)  :: a
  REAL (KIND=dp), INTENT (IN)  :: tequi
  REAL (KIND=dp), INTENT (OUT) :: uequi
  INTEGER, INTENT (OUT)        :: nequi
  REAL (KIND=dp), INTENT (OUT) :: tmin
  REAL (KIND=dp), INTENT (OUT) :: udist(ndist)
  REAL (KIND=dp), INTENT (OUT) :: tdist(ndist)
  REAL (KIND=dp), INTENT (OUT) :: cdist(ndist)

  REAL (KIND=dp)               :: ludist(ndist)
  INTEGER                      :: i, i0, ntmp
  REAL (KIND=dp)               :: vol_gr, alp, d_lu, u_min
  REAL (KIND=dp)               :: u0, t0, c0
  REAL (KIND=dp)               :: auxil_1, auxil_2
  REAL (KIND=dp)               :: aa, bb, xx0, xx
  INTEGER                      :: F_LOGS = 1
  INTEGER                      :: jlo

  ntmp = n_temp(nt)

  vol_gr = 4.0_dp / 3.0_dp * xpi * a**3         ! grain volume

  ! get alp = dlogC/dlogT from points 2 and 6:
  alp = (lcal(6) - lcal(2)) / (temp_cal(nt,6) - temp_cal(nt,2))

  ! Fix normalisation (use second point)
  i0 = 2
  t0 = t_cal(i0)
  c0 = cal(i0)
  u0 = c0 * t0 / (1.0_dp+alp)
  auxil_1 = (alp+1.0_dp) * t0**alp / c0
  auxil_2 = 1.0_dp / (alp+1.0_dp)

  ! get Umin and Uequi
  jlo = 1
  tmin = t_cmb
  u_min = 10.0_dp**INTPOL3(lu_cal(1:ntmp), temp_cal(nt,1:ntmp), ntmp, LOG10(tmin), jlo)
  uequi = 10.0_dp**INTPOL3(lu_cal(1:ntmp), temp_cal(nt,1:ntmp), ntmp, LOG10(tequi), jlo)

  ! define U scale
  IF (uequi - hnumax > u_min) THEN ! equilibrium case
     d_lu = hnumax *   (1.0_dp + 2.5_dp * LOG(uequi / hnumax)) / DBLE(ndist/2)
     d_lu = MIN(d_lu, (uequi - u_min) / DBLE(ndist/2))
     DO i=1,ndist
        udist(i) = uequi - DBLE(ndist/2-i+1) * d_lu
     ENDDO
     nequi = ndist / 2
  ELSE
     IF (F_LOGS == 1) THEN ! fluctuation case
!        aa = 4.5_dp * (uequi+hnumax) / 3.0_dp
        aa = 1.2_dp * (uequi+hnumax)
        xx0 = 0.15_dp
        bb = 70.0_dp
        DO i=1,ndist
           xx = DBLE(i) / DBLE(ndist)
           udist(i) = aa * (LOG(COSH(bb*(xx-xx0)) / COSH(-bb*xx0)) / bb + xx)
        ENDDO
        nequi = 1
     ELSE
        ! Test with constant DELTU
        d_lu = hnumax / DBLE(ndist/2)
        DO i=1,ndist
           udist(i) = u_min + DBLE(i-1) * d_lu
        ENDDO
        nequi = 1
     ENDIF
  ENDIF

  ! get T and C
  ludist(1:ndist) = LOG10(udist(1:ndist))
  jlo = 1
  DO i=1,ndist
     IF (udist(i) < u0) THEN
        tdist(i) = (auxil_1 * udist(i))**auxil_2
        cdist(i) = (alp+1.0_dp) * udist(i) / tdist(i)
     ELSE
        tdist(i) = 10.0_dp**INTPOL3 (temp_cal(nt,1:ntmp), lu_cal(1:ntmp), ntmp, ludist(i), jlo)
        cdist(i) = 10.0_dp**INTPOL3 (lcal(1:ntmp), lu_cal(1:ntmp), ntmp, ludist(i), jlo) * vol_gr
     ENDIF
  ENDDO

END SUBROUTINE GET_TC

!----------------------------------------------------------------------------

SUBROUTINE GET_UC (nt, a, ntmp)
! gets U-grid for dP/dU 

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt, ntmp
  REAL (KIND=dp), INTENT (IN)  :: a

  INTEGER                      :: nst, i, i0
  REAL (KIND=dp)               :: vol_gr, alp
  REAL (KIND=dp)               :: u0, t0, c0
  REAL (KIND=dp)               :: sizemax(1), sizemin(1)
  INTEGER                      :: jlo

  nst  = nsize_type(nt)

  sizemax = size_type(nt,MAXLOC(size_type(nt,1:nst)))
  sizemin = size_type(nt,MINLOC(size_type(nt,1:nst)))
  IF ( (a > sizemax(1)) .OR. (a < sizemin(1)) ) THEN
     WRITE (*,*) ''
     WRITE (*,*) '(F) DM_get_tdist/GET_UC: Size out of range available in C file'
     WRITE (*,*) '                         For grain type ', TRIM(gtype(nt))
     WRITE (*,*) '                         size max  (microns)     = ', sizemax(1)
     WRITE (*,*) '                         size min  (microns)     = ', sizemin(1)
     WRITE (*,*) '                         Required size (microns) = ', a
     WRITE (*,*) ''
     STOP
  ENDIF

  vol_gr = 4.0_dp / 3.0_dp * xpi * a**3         ! grain volume
  ! interpolate Heat capacity in size and on initial T-grid (both Log and linear)
  jlo = 1
  DO i=1,ntmp
     lcal(i) = INTPOL3(calor(nt,1:nst,i), size_type(nt,1:nst), nst, a, jlo)
     cal(i) = vol_gr * 10.0_dp**lcal(i)         ! heat capacity is in erg/K/cm3
     t_cal(i) = 10.0_dp**(temp_cal(nt,i))
  ENDDO

  ! get alp = dlogC/dlogT from points 2 and 6:
  alp = (lcal(6) - lcal(2)) / (temp_cal(nt,6) - temp_cal(nt,2))

  ! Fix normalisation (use second point)
  i0 = 2
  t0 = t_cal(i0)
  c0 = cal(i0)
  u0 = c0 * t0 / (1.0_dp+alp)

  ! get Umin and Uequi
  CALL PRIMITIV2(i0, ntmp, temp_cal(nt,1:ntmp), l_10*cal(1:ntmp)*t_cal(1:ntmp), u_cal(1:ntmp))    ! l_10 is log(10.0)
  u_cal(1:ntmp) = u_cal(1:ntmp) + u0
  lu_cal(1:ntmp) = LOG10(u_cal(1:ntmp))

END SUBROUTINE GET_UC

!----------------------------------------------------------------------------

END MODULE MGET_TDIST
