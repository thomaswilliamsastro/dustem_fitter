
MODULE MGET_QEXT

  USE CONSTANTS
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: GET_QEXT, GET_QEXT_POL, GET_QEXT_CIRC, COOLING, COOL1, DCOOL1, DBETA

CONTAINS

SUBROUTINE GET_QEXT(a, nt, q_ab, q_di, g_fa)
!     gets Qabs and Qsca of grains

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  ! arguments
  INTEGER, INTENT (IN)         :: nt
  REAL (KIND=dp), INTENT (IN)  :: a
  REAL (KIND=dp), INTENT (OUT) :: q_ab(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: q_di(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: g_fa(n_qabs)

  INTEGER                      :: i
  REAL (KIND=dp)               :: sizemax(1), sizemin(1)
  REAL (KIND=dp), ALLOCATABLE  :: qabs_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: size_tmp(:)
  REAL (KIND=dp), ALLOCATABLE  :: qdiff_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: gfac_tmp(:,:)
  INTEGER                      :: nsize_tmp
  REAL (KIND=dp), ALLOCATABLE  :: tempa(:), tempd(:), tempg(:)

!============================================================================

  nsize_tmp = nsize_type(nt)
  ALLOCATE (qabs_tmp(nsize_tmp, n_qabs))
  ALLOCATE (size_tmp(nsize_tmp))
  ALLOCATE (qdiff_tmp(nsize_tmp, n_qabs))
  ALLOCATE (gfac_tmp(nsize_tmp, n_qabs))
  ALLOCATE (tempa(n_qabs))
  ALLOCATE (tempd(n_qabs))
  ALLOCATE (tempg(n_qabs))
  size_tmp       = size_type(nt, 1:nsize_tmp)
  qabs_tmp(:,:)  = q_abs(nt, 1:nsize_tmp, :)
  qdiff_tmp(:,:) = qdiff(nt, 1:nsize_tmp, :)
  gfac_tmp(:,:)  = gfac(nt, 1:nsize_tmp, :)

  ! check that size exist in Q_TYPE.DAT
  sizemax = size_tmp(MAXLOC(size_tmp))
  sizemin = size_tmp(MINLOC(size_tmp))
  IF (a > sizemax(1) .OR. a < sizemin(1) ) THEN
     WRITE (*,*) ''
     WRITE (*,*) '(F) DM_get_qext/GET_QEXT: size out of range available in Q and G files'
     WRITE (*,*) '                         for grain type ', TRIM(gtype(nt))
     WRITE (*,*) '                         min size (microns)      = ', sizemin(1)
     WRITE (*,*) '                         max size (microns)      = ', sizemax(1)
     WRITE (*,*) '                         required size (microns) = ', a
     WRITE (*,*) ''
     STOP
  ENDIF

  ! interpolate Q to size a -> q_ab(n_qabs), q_di(n_qabs)
  IF (a == size_tmp(nsize_tmp)) THEN
     q_ab(:) = qabs_tmp(nsize_tmp,:)
     q_di(:) = qdiff_tmp(nsize_tmp,:)
  ELSE
     DO i = 1, n_qabs
        q_ab(i) = INTPOL ( qabs_tmp(:,i), size_tmp(:), nsize_tmp, a )
        q_di(i) = INTPOL ( qdiff_tmp(:,i), size_tmp(:), nsize_tmp, a )
        g_fa(i) = INTPOL ( gfac_tmp(:,i), size_tmp(:), nsize_tmp, a )
     ENDDO
  ENDIF
  tempa = q_ab
  tempd = q_di
  tempg = g_fa
  DO i=1,n_qabs
     q_ab(i) = tempa(n_qabs-i+1)
     q_di(i) = tempd(n_qabs-i+1)
  ENDDO

  DEALLOCATE (qabs_tmp)
  DEALLOCATE (size_tmp)
  DEALLOCATE (qdiff_tmp)
  DEALLOCATE (gfac_tmp)
  DEALLOCATE (tempa, tempd, tempg)

END SUBROUTINE GET_QEXT

!---------------------------------------------------------------------

SUBROUTINE GET_QEXT_POL(a, nt, q1_ab, q1_di, q2_ab, q2_di,qH_ab)
!     gets Q1ext Q2ext of grains

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  ! arguments
  INTEGER, INTENT (IN)         :: nt
  REAL (KIND=dp), INTENT (IN)  :: a
  REAL (KIND=dp), INTENT (OUT) :: q1_ab(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: q2_ab(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: q1_di(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: q2_di(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: qH_ab(n_qabs)

  INTEGER                      :: i
  REAL (KIND=dp)               :: sizemax(1), sizemin(1)
  REAL (KIND=dp), ALLOCATABLE  :: q1abs_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: q2abs_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: size_tmp(:)
  REAL (KIND=dp), ALLOCATABLE  :: q1diff_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: q2diff_tmp(:,:)
  INTEGER                      :: nsize_tmp
  REAL (KIND=dp), ALLOCATABLE  :: qHabs_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: tempa1(:), tempa2(:), tempd1(:), tempd2(:), tempaH(:)

!============================================================================

  nsize_tmp = nsize_type(nt)
  ALLOCATE (q1abs_tmp(nsize_tmp, n_qabs))
  ALLOCATE (q2abs_tmp(nsize_tmp, n_qabs))
  ALLOCATE (size_tmp(nsize_tmp))
  ALLOCATE (q1diff_tmp(nsize_tmp, n_qabs))
  ALLOCATE (q2diff_tmp(nsize_tmp, n_qabs))
  ALLOCATE (qHabs_tmp(nsize_tmp, n_qabs))
  ALLOCATE (tempa1(n_qabs))
  ALLOCATE (tempa2(n_qabs))
  ALLOCATE (tempd1(n_qabs))
  ALLOCATE (tempd2(n_qabs))
  ALLOCATE (tempaH(n_qabs))
  size_tmp       = size_type(nt, 1:nsize_tmp)
  q1abs_tmp(:,:)  = q1_abs(nt, 1:nsize_tmp, :)
  q2abs_tmp(:,:)  = q2_abs(nt, 1:nsize_tmp, :)
  q1diff_tmp(:,:) = q1diff(nt, 1:nsize_tmp, :)
  q2diff_tmp(:,:) = q2diff(nt, 1:nsize_tmp, :)
  qHabs_tmp(:,:)  = qH_abs(nt, 1:nsize_tmp, :)

  ! check that size exist in Q_TYPE.DAT
  sizemax = size_tmp(MAXLOC(size_tmp))
  sizemin = size_tmp(MINLOC(size_tmp))
  IF (a > sizemax(1) .OR. a < sizemin(1) ) THEN
     WRITE (*,*) ''
     WRITE (*,*) '(F) DM_get_qext/GET_QEXT_POL: size out of range available in Q and G files'
     WRITE (*,*) '                         for grain type ', TRIM(gtype(nt))
     WRITE (*,*) '                         min size (microns)      = ', sizemin(1)
     WRITE (*,*) '                         max size (microns)      = ', sizemax(1)
     WRITE (*,*) '                         required size (microns) = ', a
     WRITE (*,*) ''
     STOP
  ENDIF

  ! interpolate Q to size a -> q_ab(n_qabs), q_di(n_qabs)
  IF (a == size_tmp(nsize_tmp)) THEN
     q1_ab(:) = q1abs_tmp(nsize_tmp,:)
     q2_ab(:) = q2abs_tmp(nsize_tmp,:)
     q1_di(:) = q1diff_tmp(nsize_tmp,:)
     q2_di(:) = q2diff_tmp(nsize_tmp,:)
     qH_ab(:) = qHabs_tmp(nsize_tmp,:)
  ELSE
     DO i = 1, n_qabs
        q1_ab(i) = INTPOL ( q1abs_tmp(:,i), size_tmp(:), nsize_tmp, a )
        q2_ab(i) = INTPOL ( q2abs_tmp(:,i), size_tmp(:), nsize_tmp, a )
        q1_di(i) = INTPOL ( q1diff_tmp(:,i), size_tmp(:), nsize_tmp, a )
        q2_di(i) = INTPOL ( q2diff_tmp(:,i), size_tmp(:), nsize_tmp, a )
        qH_ab(i) = INTPOL ( qHabs_tmp(:,i), size_tmp(:), nsize_tmp, a )
     ENDDO
  ENDIF
  tempa1 = q1_ab
  tempa2 = q2_ab
  tempd1 = q1_di
  tempd2 = q2_di
  tempaH = qH_ab
  DO i=1,n_qabs
     q1_ab(i) = tempa1(n_qabs-i+1)
     q2_ab(i) = tempa2(n_qabs-i+1)
     q1_di(i) = tempd1(n_qabs-i+1)
     q2_di(i) = tempd2(n_qabs-i+1)
     qH_ab(i) = tempaH(n_qabs-i+1)
  ENDDO

  DEALLOCATE (q1abs_tmp)
  DEALLOCATE (q2abs_tmp)
  DEALLOCATE (size_tmp)
  DEALLOCATE (q1diff_tmp)
  DEALLOCATE (q2diff_tmp)
  DEALLOCATE (qHabs_tmp)
  DEALLOCATE (tempa1,tempa2,tempd1,tempd2,tempaH)

END SUBROUTINE GET_QEXT_POL

!---------------------------------------------------------------------
SUBROUTINE GET_QEXT_CIRC(a, nt, q_ci)
  ! Get circular polarization coefficients

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  ! arguments
  INTEGER, INTENT (IN)         :: nt
  REAL (KIND=dp), INTENT (IN)  :: a
  REAL (KIND=dp), INTENT (OUT) :: q_ci(n_qabs)

  INTEGER                      :: i
  REAL (KIND=dp)               :: sizemax(1), sizemin(1)
  REAL (KIND=dp), ALLOCATABLE  :: qcirc_tmp(:,:)
  REAL (KIND=dp), ALLOCATABLE  :: size_tmp(:)
  INTEGER                      :: nsize_tmp
  REAL (KIND=dp), ALLOCATABLE  :: tempc(:)

!============================================================================

  nsize_tmp = nsize_type(nt)
  ALLOCATE (qcirc_tmp(nsize_tmp, n_qabs))
  ALLOCATE (size_tmp(nsize_tmp))
  ALLOCATE (tempc(n_qabs))
  size_tmp       = size_type(nt, 1:nsize_tmp)
  qcirc_tmp(:,:)  = qcirc(nt, 1:nsize_tmp, :)

  ! check that size exist in Q_TYPE.DAT
  sizemax = size_tmp(MAXLOC(size_tmp))
  sizemin = size_tmp(MINLOC(size_tmp))
  IF (a > sizemax(1) .OR. a < sizemin(1) ) THEN
     WRITE (*,*) ''
     WRITE (*,*) '(F) DM_get_qext/GET_QEXT_CIRC: size out of range available in Q and G files'
     WRITE (*,*) '                         for grain type ', TRIM(gtype(nt))
     WRITE (*,*) '                         min size (microns)      = ', sizemin(1)
     WRITE (*,*) '                         max size (microns)      = ', sizemax(1)
     WRITE (*,*) '                         required size (microns) = ', a
     WRITE (*,*) ''
     STOP
  ENDIF

  ! interpolate Q to size a -> q_ci(n_qabs)
  IF (a == size_tmp(nsize_tmp)) THEN
     q_ci(:) = qcirc_tmp(nsize_tmp,:)
  ELSE
     DO i = 1, n_qabs
        q_ci(i) = INTPOL ( qcirc_tmp(:,i), size_tmp(:), nsize_tmp, a )
     ENDDO
  ENDIF
  tempc = q_ci
  DO i=1,n_qabs
     q_ci(i) = tempc(n_qabs-i+1)
  ENDDO

  DEALLOCATE (qcirc_tmp)
  DEALLOCATE (size_tmp)
  DEALLOCATE (tempc)

END SUBROUTINE GET_QEXT_CIRC
!---------------------------------------------------------------------


SUBROUTINE COOLING (nt, ns, t, f, nn, a)

! cooling power (erg/s) emitted by a grain at given temperature t
! NB uses BB_lambda(T)

  USE CONSTANTS
  USE UTILITY
  USE MDTLS

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt   ! index of grain type
  INTEGER, INTENT (IN)         :: nn   ! nr of T-values
  INTEGER, INTENT (IN)         :: ns   ! index of grain size
  REAL (KIND=dp), INTENT (IN)  :: t(nn)
  REAL (KIND=dp), INTENT (IN)  :: a    ! size
  REAL (KIND=dp), INTENT (OUT) :: f(nn)

  INTEGER                      :: l, ktemp
  REAL (KIND=dp)               :: energiemise
  REAL (KIND=dp)               :: xx(n_qabs), fnut(n_qabs)
  REAL (KIND=dp)               :: sig(1)

  IF (n_beta == 0 .AND. n_dtls == 0) THEN
     ! temperature loop
     DO ktemp=1,nn
        ! frequency loop
        xx(:) = hcsurk / (t(ktemp) * lamb_qabs(:))
        DO l=1,n_qabs
           fnut(l) = F_BB (xx(l)) * qaem(n_qabs-l+1) / lamb_qabs(l)**5
        ENDDO
        energiemise = XINTEG2(1, n_qabs, n_qabs, lamb_qabs, fnut)
        f(ktemp) = cte1 * energiemise
     ENDDO
  ELSE IF (n_beta == 1 .AND. n_dtls == 0) THEN    ! apply BETA(T)
     ! temperature loop
     DO ktemp=1,nn
        ! frequency loop
        xx(:) = hcsurk / (t(ktemp) * lamb_qabs(:))
        DO l=1,n_qabs
           fnut(l) = F_BB (xx(l)) * qaem(n_qabs-l+1) * & 
           ! faster
           & EXP( DBETA(t(ktemp),nt)*f_beta(nt,l)*LOG(ltresh(nt)/lamb_qabs(l)) ) / &
           & lamb_qabs(l)**5
        ENDDO
        energiemise = XINTEG2(1, n_qabs, n_qabs, lamb_qabs, fnut)
        f(ktemp) = cte1 * energiemise
     ENDDO
  ELSE IF (n_beta == 0 .AND. n_dtls == 1) THEN    ! apply DCD/TLS
     ! temperature loop
     DO ktemp=1,nn
        ! frequency loop
        xx(:) = hcsurk / (t(ktemp) * lamb_qabs(:))
        DO l=1,n_qabs
           IF (lamb_qabs(l)*1e4_dp < ldtresh(nt)) THEN 
              fnut(l) = F_BB (xx(l)) * qaem(n_qabs-l+1) / lamb_qabs(l)**5
           ELSE 
              CALL DTLS (nt, 1, ns, n_qabs-l+1, a, t(ktemp), sig)
              fnut(l) = F_BB (xx(l)) * sig(1) / lamb_qabs(l)**5
           ENDIF
        ENDDO
        energiemise = XINTEG2(1, n_qabs, n_qabs, lamb_qabs, fnut)
        f(ktemp) = cte1 * energiemise
     ENDDO
  ENDIF

END SUBROUTINE COOLING

!---------------------------------------------------------------------

FUNCTION COOL1 (tt, nt, ns, a)

! cooling power (erg/s) emitted by a grain at given temperature t
! used to find Tequil
! NB uses BB_nu(T)

  USE CONSTANTS
  USE UTILITY
  USE MDTLS

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt      ! index of grain type
  REAL (KIND=dp), INTENT (IN)  :: tt
  REAL (KIND=dp)               :: COOL1
  REAL (KIND=dp), INTENT (IN)  :: a       ! grain size
  INTEGER, INTENT (IN)         :: ns       ! index of grain size

  INTEGER                      :: l
  REAL (KIND=dp)               :: xx(n_qabs), fnut(n_qabs)
  REAL (KIND=dp)               :: sig(1), temp(1)

  ! loop on frequencies
  xx(:) = xhp * freq_qabs(:) / (xkb * tt)
  IF (n_beta == 0 .AND. n_dtls == 0) THEN
     DO l=1,n_qabs
        fnut(l) = F_BB (xx(l)) * qaem(l) * freq_qabs(l)**3 
     ENDDO
  ELSE IF (n_beta == 1 .AND. n_dtls == 0) THEN   ! apply BETA(T)
     DO l=1,n_qabs
        fnut(l) = F_BB (xx(l)) * qaem(l) * freq_qabs(l)**3 * &
        ! faster
        & EXP( DBETA(tt,nt)*f_beta(nt,n_qabs-l+1)*LOG(freq_qabs(l)*ltresh(nt)/clight) )
     ENDDO
  ELSE IF (n_beta == 0 .AND. n_dtls == 1) THEN   ! apply DCD/TLS
     temp(1) = tt
     DO l=1,n_qabs
        IF (lamb_qabs(n_qabs-l+1)*1e4_dp < ldtresh(nt)) THEN 
           fnut(l) = F_BB (xx(l)) * freq_qabs(l)**3 * qaem(l)
        ELSE 
           CALL DTLS (nt, 1, ns, l, a, temp, sig)
           fnut(l) = F_BB (xx(l)) * freq_qabs(l)**3 * sig(1)
        ENDIF
     ENDDO
  ENDIF

  COOL1 = cte2 * XINTEG2(1, n_qabs, n_qabs, freq_qabs, fnut)

END FUNCTION COOL1

!---------------------------------------------------------------------

FUNCTION DCOOL1 (tt, nt, ns, a)

! derivative of COOL1: used to find Tequil

  USE CONSTANTS
  USE UTILITY
  USE MDTLS

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt   ! index of grain type
  REAL (KIND=dp), INTENT (IN)  :: tt
  REAL (KIND=dp), INTENT (IN)  :: a    ! grain size
  INTEGER, INTENT (IN)         :: ns   ! index of grain size
  REAL (KIND=dp)               :: DCOOL1

  INTEGER                      :: l
  REAL (KIND=dp)               :: xx(n_qabs), fnut(n_qabs)
  REAL (KIND=dp)               :: sig(1), temp(1)

  ! loop on frequencies
  xx(:) = xhp  * freq_qabs(:) / (xkb * tt)
  IF (n_beta == 0 .AND. n_dtls == 0) THEN
     DO l=1,n_qabs
        fnut(l) = G_BB (xx(l)) * qaem(l) * freq_qabs(l)**4
     ENDDO
  ELSE IF (n_beta == 1 .AND. n_dtls == 0) THEN  ! apply BETA(T)
     DO l=1,n_qabs
        fnut(l) = G_BB (xx(l)) * qaem(l) * freq_qabs(l)**4 * &
        ! faster
        & EXP( DBETA(tt,nt)*f_beta(nt,n_qabs-l+1)*LOG(freq_qabs(l)*ltresh(nt)/clight) )
     ENDDO
  ELSE IF (n_beta == 0 .AND. n_dtls == 1) THEN  ! apply DCD/TLS
     temp(1) = tt
     DO l=1,n_qabs
!VG Feb 2016: corrected error on indices
!        IF (lamb_qabs(n_qabs+l-1)*1e4_dp < ldtresh(nt)) THEN 
        IF (lamb_qabs(n_qabs-l+1)*1e4_dp < ldtresh(nt)) THEN 
           fnut(l) = G_BB (xx(l)) * freq_qabs(l)**4 * qaem(l)
        ELSE 
           CALL DTLS (nt, 1, ns, l, a, temp, sig)
           fnut(l) = G_BB (xx(l)) * freq_qabs(l)**4 * sig(1)
        ENDIF
     ENDDO
  ENDIF

  DCOOL1 = cte2 * XINTEG2(1, n_qabs, n_qabs, freq_qabs, fnut) * 0.5_dp * (xhp / xkb) / tt**2

END FUNCTION DCOOL1

!---------------------------------------------------------------------

FUNCTION DBETA (tt, nt)
! correction to generate a beta(T) behaviour of Qabs above some lambda(nu) threshold
! Q(nu)=Q0(nu)*(nu/nutresh)**(DBETA(T)*F_BETA(nu,nutresh)) 
! with BETA(T) = BETA0 + DBETA(T)

  USE CONSTANTS
  USE UTILITY

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt   ! index of grain type
  REAL (KIND=dp), INTENT (IN)  :: tt
  REAL (KIND=dp)               :: DBETA

  INTEGER                      :: jlo
  REAL (KIND=dp)               :: tmp, bm

  if (nbeta(nt) == 0) THEN
     bm = -beta0(nt) + bmax(nt)
     !  tmp  = -beta0(nt) + abeta(nt) * tt**(gbeta(nt))
     tmp  = -beta0(nt) + abeta(nt) * EXP(gbeta(nt)*LOG(tt))  ! faster
     IF (tmp > bm) tmp = bm
  ELSE 
     jlo = 1
     ! constant at edges (extrapolation)
     tmp = -beta0(nt) + INTPOL2( betav(nt,1:nbeta(nt)), tbeta(nt,1:nbeta(nt)), nbeta(nt), tt )    
  ENDIF

  DBETA = tmp

END FUNCTION DBETA

!---------------------------------------------------------------------

END MODULE MGET_QEXT
