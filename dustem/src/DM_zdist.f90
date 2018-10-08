MODULE MZDIST
  
  USE CONSTANTS
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: ZDIST



CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE ZDIST (nt, ns, a, zav, sig_z, jpel, hpel, hspel, jgg, cgg)
! 
! computes the charge distribution of dust grains
! Weingartner & Draine 01 (WD01) formalism with van Hoof et al 04 revision for Emin (also Weingartner et al 06, W06)
!

  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE
 
    ! arguments
    INTEGER,        INTENT (IN)    :: nt, ns                       ! index of grain type and size
    REAL (KIND=dp), INTENT (IN)    :: a                            ! grain radius (cm)
    REAL (KIND=dp), INTENT (OUT)   :: zav, sig_z                   ! mean charge and standard deviation
    REAL (KIND=dp), INTENT(OUT)    :: jpel, hpel, hspel, jgg, cgg  ! PE rate and heating, gas-grain rates and cooling

    ! local variables
    INTEGER                     :: i, k, imin, nzbi, nzbl
    REAL (KIND=dp)              :: uait, tmp                       ! autoionization threshold
    REAL (KIND=dp)              :: aa, zqa, zqb, zq, zc, jj, hh, hs, xqe2a
    REAL (KIND=dp)              :: z1(1), je(1), tt(1)
    REAL (KIND=dp), ALLOCATABLE :: zbl(:), fzl(:), jp(:), jm(:), cp(:), cm(:)


    aa = a/1.0e-8_dp   ! cm to Angstroems

    ! min and max charges
    uait = -(p_uait(nt,1) + p_uait(nt,2)*aa + p_uait(nt,3)/aa)
    zmin = DBLE( FLOOR(uait/14.4_dp*aa) + 1 )
    zmax = DBLE( FLOOR( ((hnumax/everg-wf(nt))/14.4_dp*aa + 0.5_dp - 0.3_dp/aa) / (1.0_dp + 0.3_dp/aa) ) )

    ! first estimate Zeq solution of Jpe+Jp=Jm
    zqa = zmin
    zqb = zmax
    CALL GET_ZEQ (nt, ns, zqa, zqb, aa)
    zeq = (zqa+zqb)/2.0_dp 
    zq = DBLE( NINT(zeq) )

    ! get nr of Z bins and Z grid
    nzbi = nz_bg + NINT( DBLE(nz_sg-nz_bg)/(ABS(zq)/ztrans + 1.0_dp) )
    nzbl = MIN( nzbi, NINT(ABS(zmax-zmin))+1)
    IF (nzbl /= nzbi) THEN  ! get the charge midpoint
       ! use the full Z range (small grains)
       zc = DBLE(NINT((zmin+zmax)/2.0_dp))
       DO i=1,nzbl 
          tmp = zmin + DBLE(i-1)
          IF (tmp <= zmax) k = i
       ENDDO
       nzbl = k
       ALLOCATE (zbl(nzbl))
       DO i=1,nzbl
          zbl(i) = zmin + DBLE(i-1)
       ENDDO
    ELSE
       ! use a range around Zeq (big grains)
       zc = zq
       ALLOCATE (zbl(nzbl))
       DO i=1,nzbl
          zbl(i) = zc + DBLE(i-nzbl/2-1)
       ENDDO
    ENDIF
    ALLOCATE(fzl(nzbl),jp(nzbl),jm(nzbl),cp(nzbl),cm(nzbl))

    ! get full charge distribution
    CALL GAS_CURRENTS (nt, nzbl, 1, nion, zbl, aa, t_gas, jp, jm, cp, cm)
    fzl = 0.0_dp
    fzl(nzbl/2+1) = 1.0_dp ! start value @ z = zc
    imin = 0
    DO i=nzbl/2+1, 2, -1
       CALL PE_CURRENT( nt, ns, zbl(i-1), aa, jpel, hpel, hspel )
       tmp = jm(i) / (jpel + jp(i-1) )
       IF ( ((tmp>istiny) .AND. (tmp<1.0_dp/istiny)).AND. ((fzl(i)>tiniest/istiny) .AND. (fzl(i)<istiny/tiniest)) )THEN 
          fzl(i-1) = fzl(i) * tmp
       ELSE 
          fzl(i-1) = 0.0_dp
       ENDIF
       IF (zbl(i-1)==zmin) imin = i-1
    ENDDO
    DO i=nzbl/2+1,nzbl-1
       CALL PE_CURRENT( nt, ns, zbl(i), aa, jpel, hpel, hspel )
       tmp = (jpel + jp(i) ) / jm(i+1)
       IF ( ((tmp>istiny) .AND. (tmp<1.0_dp/istiny)).AND. ((fzl(i)>tiniest/istiny) .AND. (fzl(i)<istiny/tiniest)) )THEN 
          fzl(i+1) = fzl(i) * tmp
       ELSE 
          fzl(i+1) = 0.0_dp
       ENDIF
    ENDDO

    ! build the final charge distribution
    fzl = fzl / SUM(fzl)
    nzb = 0
    DO i=1, nzbl
       IF ((fzl(i)>fzmin) .AND. (fzl(i)<=1.0_dp)) nzb = nzb + 1
    ENDDO
    DEALLOCATE (jp,jm,cp,cm)
    IF (nzb>0) THEN 
       ALLOCATE (zb(nzb),fz(nzb)) 
       ALLOCATE (jp(nzb),jm(nzb),cp(nzb),cm(nzb))
       k = 0
       DO i=1, nzbl
          IF ((fzl(i)>fzmin) .AND. (fzl(i)<=1.0_dp)) THEN 
             k = k + 1
             zb(k) = zbl(i)
             fz(k) = fzl(i)
          ENDIF
       ENDDO
       fz = fz / SUM(fz)

       ! compute the average charge and final rates
       zav = SUM(zb*fz)                       ! mean charge 
       sig_z = SQRT(SUM(zb**2*fz)-zav**2)     ! variance of fz 
       jpel = 0.0_dp
       hpel = 0.0_dp
       hspel = 0.0_dp
       DO i = 1, nzb 
          CALL PE_CURRENT( nt, ns, zb(i), aa, jj, hh, hs)
          jpel = jpel + jj*fz(i)
          hpel = hpel + hh*fz(i)
          hspel = hspel + hh*fz(i)
       ENDDO

       CALL GAS_CURRENTS (nt, nzb, 1, nion, zb, aa, t_gas, jp, jm, cp, cm)
       jgg = SUM( (jp + jm)*fz )
       cgg = SUM( (cp + cm)*fz )

       IF (imin > 0) THEN
          z1(1) = zmin
          CALL GAS_CURRENTS (nt, 1, 1, 1, z1, aa, t_gas, tt, je, tt, tt)
          xqe2a = xqe2/(1.e-8_dp*aa)
          cgg = cgg + fzl(imin)*je(1)*( (wf(nt)-ebg(nt,ns))*everg + &
               & (zmin-5e-1_dp)*xqe2a-(p_ea(nt,1)/(aa+p_ea(nt,2)) )*xqe2a)
       ENDIF

    ELSE 
       WRITE (*,*) ''
       WRITE (*,*) '(F) DM_zdist/ZDIST: f(Z) has no values above ',fzmin
       WRITE (*,*) ' grain type= ',gtype(nt),' a (nm)= ',aa
       WRITE (*,*) ''
       STOP
    ENDIF

  END SUBROUTINE ZDIST

!----------------------------------------------------------------------------
  SUBROUTINE GET_ZEQ (nt, ns, z_a, z_b, a)
  ! finds equilibirum charge from Jpe + Jp = Jm

    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)           :: nt               ! index of grain type
    INTEGER, INTENT (IN)           :: ns               ! index of grain size
    REAL (KIND=dp), INTENT (INOUT) :: z_a, z_b
    REAL (KIND=dp), INTENT (IN)    :: a                ! grain size

    INTEGER                        :: i
    REAL (KIND=dp)                 :: z, fa, fb, f, je, hh, hs
    REAL (KIND=dp)                 :: z1(1), jp(1), jm(1), cp(1), cm(1)

    z1(1) = z_a
    CALL GAS_CURRENTS (nt, 1, 1, nion, z1, a, t_gas, jp, jm, cp, cm)
    CALL PE_CURRENT ( nt, ns, z_a, a, je, hh, hs )
    fa = je + jp(1) - jm(1)

    z1(1) = z_b
    CALL GAS_CURRENTS (nt, 1, 1, nion, z1, a, t_gas, jp, jm, cp, cm)
    CALL PE_CURRENT ( nt, ns, z_b, a, je, hh, hs )
    fb = je + jp(1) - jm(1)
    
    i = 0
    IF (fa*fb > 0.0_dp) THEN
       PRINT *, "  (W) GET_ZEQ:  Wrong initial guess"
       PRINT *, "      Za = ", z_a, "  fa = ", fa
       PRINT *, "      Zb = ", z_b, "  fb = ", fb
    ENDIF
    DO
       z = 0.5_dp * (z_a + z_b)
       z1(1) = z
       CALL GAS_CURRENTS (nt, 1, 1, nion, z1, a, t_gas, jp, jm, cp, cm)
       CALL PE_CURRENT ( nt, ns, z, a, je, hh, hs )
       f = je + jp(1) - jm(1)
       IF (f*fa > 0.0_dp) THEN
          z_a = z
          fa = f
       ELSE
          z_b = z
          fb = f
       ENDIF
       i = i + 1
       IF ((z_b-z_a) < 5.0e-1_dp) EXIT
    ENDDO
    
END SUBROUTINE GET_ZEQ

!----------------------------------------------------------------------------

SUBROUTINE GAS_CURRENTS (nt, nz, i1, i2, zg, ag, tgas, jpi, jmi, cpi, cmi)
! computes the gas current (s-1) from ions or electrons for the charge distribution

    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)         :: nt, nz           ! index for grain type, nr of charge bins
    INTEGER, INTENT (IN)         :: i1, i2           ! min and max indices to select ions
    REAL (KIND=dp), INTENT (IN)  :: ag, zg(nz)       ! size (Angstroems) and charge grid of grain (unit of e)
    REAL (KIND=dp), INTENT (IN)  :: tgas             ! gas temperature
    REAL (KIND=dp), INTENT (OUT) :: jpi(nz), jmi(nz) ! rates for positive and negative currents
    REAL (KIND=dp), INTENT (OUT) :: cpi(nz), cmi(nz) ! cooling rates for positive and negative currents

    INTEGER                      :: k
    REAL (KIND=dp)               :: mii, tau
    REAL (KIND=dp)               :: nu(nz), theta(nz), jt(nz), ct(nz), stick(nz)

    jpi = 0.0_dp
    jmi = 0.0_dp
    cpi = 0.0_dp
    cmi = 0.0_dp

    ! ion list in GAS.DAT
    DO k = i1, i2
       tau = (ag*1e-8_dp) * xkb*tgas / xqe2/zi(k)**2
       nu = zg/zi(k)
       mii = mi(k)*amu 

       ! get Jtilde and lambda_tilde function from Draine & Sutin 87
       WHERE (nu<0.0_dp) 
          jt = (1.0_dp-nu/tau) * ( 1.0_dp + SQRT(2.0_dp/(tau-2.0_dp*nu)) )
          ct = (2.0_dp-nu/tau) * ( 1.0_dp + 1.0_dp/SQRT(tau-nu) ) 
       ENDWHERE
       WHERE (nu==0.0_dp) 
          jt = 1.0_dp + SQRT(xpi/2.0_dp/tau)
          ct = 2.0_dp + 1.5_dp*SQRT(xpi/2.0_dp/tau)
       ENDWHERE
       WHERE (nu>0.0_dp)
          theta = nu / ( 1.0_dp + 1.0_dp/SQRT(nu) )
          jt = ( 1.0_dp + 1.0_dp/SQRT(4.0_dp*tau + 3.0_dp*nu) )**2 * EXP(-theta/tau)
          ct = (2.0_dp+nu/tau)*( 1.0_dp + 1.0_dp/SQRT(3.0_dp/2.0_dp/tau + 3.0_dp*nu) )*EXP(-theta/tau) 
       ENDWHERE

       ! get sticking coefficient (default is WD01)
       IF ( (zi(k) == -1.0_dp) .AND. (ABS(mii/xme-1.0_dp)-1.0_dp < 1e-2_dp) ) THEN 
          ! electrons
          IF ((INDEX(t_opt(nt),'WD')>0) .OR. (INDEX(t_opt(nt),'WD')+(INDEX(t_opt(nt),'BT')) == 0)) THEN
             ! default is WD01
             stick = 0.5_dp * (1.0_dp - EXP(-ag/le(nt)))
             WHERE (zg <= 0.0_dp) stick = stick / (1.0_dp + EXP(2e1_dp - 4.68e-1_dp*ag**3))  ! correction Eq.28 of WD01
          ELSE IF (INDEX(t_opt(nt),'BT')>0) THEN 
             stick = 1.0_dp
          ENDIF
       ELSE
          ! ions
          stick = 1.0_dp
       ENDIF

       IF (zi(k) < 0.0_dp) THEN 
          jmi = jmi + xpi*(ag*1.0e-8_dp)**2 * stick * iden(k) * SQRT(8.0_dp*xkb*tgas/xpi/mii) * jt
          cmi = cmi + xpi*(ag*1.0e-8_dp)**2 * stick*iden(k)*SQRT(8.0_dp*xkb*tgas/xpi/mii)*xkb*tgas * ct
       ELSE IF (zi(k) > 0.0_dp) THEN
          jpi = jpi + xpi*(ag*1.0e-8_dp)**2 * stick * iden(k) * SQRT(8.0_dp*xkb*tgas/xpi/mii) * jt
          cpi = cpi + xpi*(ag*1.0e-8_dp)**2 * stick*iden(k)*SQRT(8.0_dp*xkb*tgas/xpi/mii)*xkb*tgas * ct
       ENDIF

    ENDDO

  END SUBROUTINE GAS_CURRENTS
  
!----------------------------------------------------------------------------

  SUBROUTINE PE_CURRENT( nt, ns, zg, ag, jp, hp, hsp )
    ! computes the photelectric current (s-1) for the charge distribution

    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    REAL (KIND=dp), INTENT (OUT) :: jp, hp, hsp

    INTEGER, INTENT (IN)         :: nt, ns          ! index for grain type and size
    REAL (KIND=dp), INTENT (IN)  :: ag, zg          ! size (Angstroems) and charge grid of grain (unit of e)

    INTEGER                      :: i, k, nfrq, jtresh   
    REAL (KIND=dp)               :: jp1, hp1, hsp1
    REAL (KIND=dp)               :: xqe2a, emin, theta, ipdt, ipet, n1
    REAL (KIND=dp), ALLOCATABLE  :: hnu(:), teeta(:), elow(:), ehigh(:), e2(:), alpha(:), beta(:), tt(:)
    REAL (KIND=dp), ALLOCATABLE  :: y0(:), y1(:), y2(:), yt(:)
    REAL (KIND=dp)               :: jp2, hp2, hsp2
    REAL (KIND=dp), ALLOCATABLE  :: sig_pd(:)

    xqe2a = xqe2/(1.e-8_dp*ag)/everg   ! in eV

    ! get Emin (eV)
    IF (zg < -1.0_dp) THEN 
       theta = ABS(zg+1.0_dp) / ( 1.0_dp + 1.0_dp/SQRT(ABS(zg+1.0_dp)) )
       emin = theta*xqe2a * ( 1.0_dp - 3e-1_dp/(ag/10.0_dp)**4.5e-1_dp/(ABS(zg+1.0_dp))**2.6e-1_dp )  ! W06 Eq.3 (see also van Hoof et al 04)
!       emin = ABS(zg+1.0_dp)*xqe2a / (1.0_dp+(27.0_dp/ag)**7.5e-1_dp) ! WD01
    ELSE 
       emin = 0.0_dp
    ENDIF

    ! PE threshhold and frequency grid
    ipet = emin + wf(nt) + (zg+5e-1_dp)*xqe2a + (zg+2.0_dp)*(3e-1_dp/ag)*xqe2a    ! WD01 Eqs. 2 and 6
    IF ( (1.0_dp-ipet*everg/hnumax) > 1e-4_dp ) THEN ! PE works only if ipet < hnumax
       ALLOCATE (hnu(n_qabs))
       hnu = xhp*freq_qabs/everg
       i = 1
       DO WHILE (hnu(i) <= ipet) 
          i = i+1
       ENDDO
       jtresh = i
       nfrq = jfreqmax-jtresh+1
       DEALLOCATE (hnu)
       ALLOCATE (hnu(nfrq),teeta(nfrq),elow(nfrq),ehigh(nfrq),e2(nfrq),alpha(nfrq),beta(nfrq),tt(nfrq))
       ALLOCATE(y0(nfrq),y1(nfrq),y2(nfrq),yt(nfrq))
       hnu = xhp*freq_qabs(jtresh:jfreqmax)/everg

    ! get PE yield WD01
       n1 = 4.0_dp*xpi*(ag*1e-8_dp)**3*rho(nt,ns)/3.0_dp/m1mass(nt)    ! nr of species (e.g., C, SiO4) in grain / Navogadro
       ! qi_abs sorted as freq_qabs (reverse to lamb_qabs) after CALL GET_QEXT
       beta = qi_abs(nt,ns,jtresh:jfreqmax) * xpi*(ag*1.0e-8_dp)**3 * rho(nt,ns)/m1mass(nt) / n1
       alpha = beta + ag/le(nt)
       tt = alpha**2 - 2.0_dp*alpha + 2.0_dp*(1.0_dp-EXP(-alpha))
       y1 = (beta/alpha)**2 * tt / ( beta**2 - 2.0_dp*beta + 2.0_dp*(1.0_dp-EXP(-beta)) )

       ehigh = emin + hnu - ipet
       IF (zg >= 0.0_dp) THEN 
          elow = -(zg+1.0_dp)*xqe2a
          teeta = hnu - ipet - elow 
          y2 = ehigh**2 * (ehigh-3.0_dp*elow) / (ehigh-elow)**3 
          e2 = ehigh * (ehigh-2.0_dp*elow) / (ehigh-3.0_dp*elow)
       ELSE
          elow = emin
          teeta = hnu - ipet
          y2 = 1.0_dp
          e2 = (ehigh+elow) / 2.0_dp
       ENDIF
       y0 = p_y(nt,1)*(teeta/wf(nt))**p_y(nt,3) / (1.0_dp + p_y(nt,2)*(teeta/wf(nt))**p_y(nt,3))
       DO i=jtresh,jfreqmax
          k = i-jtresh+1
          yt(k) = y2(k) * MIN( y0(k)*y1(k), 1.0_dp )
       ENDDO
       tt = yt*qi_abs(nt,ns,jtresh:jfreqmax)*isrfuv(jtresh:jfreqmax)*(e2+(zg+1.0_dp)*xqe2a)/hnu
       hsp1 = xpi*(ag*1e-8_dp)**2 * XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )
       tt = yt*qi_abs(nt,ns,jtresh:jfreqmax)*isrfuv(jtresh:jfreqmax)*e2/hnu
       hp1 = xpi*(ag*1e-8_dp)**2 * XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )
       tt = tt/e2/everg
       jp1 = xpi*(ag*1e-8_dp)**2 * XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )

    ELSE   ! no photons of energies above PE threshhold

       jp1 = 0.0_dp
       hp1 = 0.0_dp
       hsp1 = 0.0_dp

    ENDIF

    ! photodetachment part: !!!! hnu_pdt = Emin + EA(Z+1) !!!! Eq. 18 of WD01
    ipdt = emin + wf(nt) - ebg(nt,ns) + (zg+5e-1_dp)*xqe2a - ( p_ea(nt,1)/(ag+p_ea(nt,2)) )*xqe2a
    IF ((ipdt < hnumax/everg) .AND. (zg < 0.0_dp)) THEN  ! for photons of energies above PD threshhold
       IF (ALLOCATED(hnu)) DEALLOCATE (hnu)
       IF (ALLOCATED(tt)) DEALLOCATE (tt)
       ALLOCATE (hnu(n_qabs))
       hnu = xhp*freq_qabs/everg
       i = 1
       DO WHILE (hnu(i) <= ipdt) 
          i = i+1
       ENDDO
       jtresh = i
       nfrq = jfreqmax-jtresh+1
       DEALLOCATE (hnu)
       ALLOCATE (hnu(nfrq),sig_pd(nfrq),tt(nfrq))
       hnu = xhp*freq_qabs(jtresh:jfreqmax)/everg
       tt = (hnu-ipdt)/s_ea(nt,2)
       sig_pd = ABS(zg)*(2.0_dp*xpi*xqe2*xhp*s_ea(nt,1)/3.0_dp/xme/clight/s_ea(nt,2)/everg) * & 
            & tt/(1.0_dp+tt**2/3.0_dp)**2
       sig_pd = sig_pd * (ag/2.34_dp)**2 ! normalized to C6F6- for which the PD cross-section has been measured (WD01)

       tt = sig_pd*isrfuv(jtresh:jfreqmax)*(hnu-ipdt+emin+(zg+1.0_dp)*xqe2a)/hnu
       hsp2 = XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )
       tt = sig_pd*isrfuv(jtresh:jfreqmax)*(hnu-ipdt+emin)/hnu
       hp2 = XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )
       tt = tt/(hnu-ipdt+emin)/everg
       jp2 = XINTEG2( 1, nfrq, nfrq, freq_qabs(jtresh:jfreqmax), tt )

    ELSE   ! no photons of energies above PD threshhold

       jp2 = 0.0_dp
       hp2 = 0.0_dp
       hsp2 = 0.0_dp

    ENDIF
    
    jp = jp1 + jp2
    hp = hp1 + hp2
    hsp = hsp1 + hsp2

  END SUBROUTINE PE_CURRENT

!----------------------------------------------------------------------------

END MODULE MZDIST
