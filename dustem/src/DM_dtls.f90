MODULE MDTLS
  
  USE CONSTANTS
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: GET_DCD, GET_RES, GET_PHON, GET_HOP, DTLS

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE DTLS (nt, ntsed, ns, iqabs, a, t, sig)
! 
! computes the absorption cross-section for DCD-TLS effects as in Meny et al. (2007)
!

  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE
 
    INTEGER, INTENT (IN)         :: nt                   ! index of grain type
    INTEGER, INTENT (IN)         :: ntsed                ! nr of T-values
    INTEGER, INTENT (IN)         :: iqabs                ! index of frequency
    INTEGER, INTENT (IN)         :: ns                   ! index of grain size
    REAL (KIND=dp), INTENT (IN)  :: a                    ! grain size
    REAL (KIND=dp), INTENT (IN)  :: t(ntsed)             ! temperature vector
    REAL (KIND=dp), INTENT (OUT) :: sig(ntsed)           ! absorption cross-section
    REAL (KIND=dp)               :: res(ntsed), phon(ntsed), hop(ntsed), s0(ntsed)
    REAL (KIND=dp)               :: Kdcd, Ktls, dcd
    REAL (KIND=dp)               :: omega

    omega = 2.0_dp * xpi * freq_qabs(iqabs)
    CALL GET_DCD  (omega, nt, dcd)
    CALL GET_RES  (omega, ntsed, nt, t, res)
    CALL GET_PHON (omega, ntsed, nt, t, phon)
    CALL GET_HOP  (omega, ntsed, nt, t, hop)
    Kdcd = xqe2 / (16.0_dp*amu) / (3.0_dp * vt(nt)**3 * clight * rhom(nt))
    Ktls = 4.0_dp/3.0_dp * xpi**2 * Pmu(nt) / clight / rhom(nt)
    sig(:) = Kdcd*dcd + a_dtls(nt)*Ktls*(phon(:)+res(:)+hop(:) ) !DEB; j'ai multiplié les effets TLS par a_dtls plutot que le DCD
    sig(:) = sig(:) * 4.0_dp/3.0_dp * rhom(nt) * a

    ! get DTLS value @ ldtresh
    omega = 2.0_dp * xpi * (1e4_dp*clight/ldtresh(nt))
    CALL GET_DCD  (omega, nt, dcd)
    CALL GET_RES  (omega, ntsed, nt, t, res)
    CALL GET_PHON (omega, ntsed, nt, t, phon)
    CALL GET_HOP  (omega, ntsed, nt, t, hop)
    Kdcd = xqe2 / (16.0_dp*amu) / (3.0_dp * vt(nt)**3 * clight * rhom(nt))
    Ktls = 4.0_dp/3.0_dp * xpi**2 * Pmu(nt) / clight / rhom(nt)
    s0(:) = Kdcd*dcd + a_dtls(nt)*Ktls*(phon(:)+ res(:)+hop(:) )  ! voir remarque précedente
    s0(:) = s0(:) * 4.0_dp/3.0_dp * rhom(nt) * a

    ! normalize to ref value @ ldtresh
    sig(:) = Qdtls(nt,ns) * sig(:) / s0(:)            



  END SUBROUTINE DTLS


  SUBROUTINE GET_DCD (omega, nt, dcd)
!
! compute opacity for disorderd charge distribution
!
  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)        :: nt
    REAL (KIND=dp), INTENT(IN)  :: omega
    REAL (KIND=dp), INTENT(OUT) :: dcd
    REAL (KIND=dp)              :: omega_c

    IF (omega > 2.0_dp*xpi*clight/10.0e-4_dp) THEN
       dcd = 0.0_dp   ! => DCD effect valid for 10 microns < lambda
    ELSE
       omega_c = 2.0_dp*xpi*vt(nt) / (lc(nt)*1.0e-7_dp)  ! DEB: il manque le 2pi que j'ai rajouté
       dcd     = omega**2 * ( 1.0_dp - (1.0_dp + (omega/omega_c)**2)**(-2) )
    ENDIF

  END SUBROUTINE GET_DCD


  SUBROUTINE GET_RES (omega, ntsed, nt, t, res)
!
! compute opacity for resonant absorption    
!
  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)        :: ntsed, nt
    REAL (KIND=dp), INTENT(IN)  :: omega
    REAL (KIND=dp), INTENT(IN)  :: t(ntsed)
    REAL (KIND=dp), INTENT(OUT) :: res(ntsed)
    REAL (KIND=dp)              :: x, gres, g1, g2
    
    x     = omega / omega_m(nt)
    g1    = 4.0_dp/15.0_dp * (5.0_dp - 6.0_dp * x**2)
    IF (x <= 1.0_dp) THEN
       gres = 1.0_dp + g1 * x**2
    ELSE IF (x > 1.0_dp .AND. x < 340.0_dp) THEN
       g2   = 8.0_dp/15.0_dp * (1.0_dp - x**2) * (2.0_dp + 3.0_dp * x**2)
       gres = 1.0_dp + g1 * x**2 - g2 * sqrt(1.0_dp - 1.0_dp/x**2)
    ELSE
       gres = 0.0_dp
    ENDIF
 !   res(:) = omega * gres * tanh(xhbar*omega/2.0_dp/xkb/t(:))
    res(:) = omega * Pmu(nt) * tanh(xhbar*omega/2.0_dp/xkb/t(:))
  END SUBROUTINE GET_RES


  SUBROUTINE GET_PHON (omega, ntsed, nt, t, phon)
!
! compute opacity for tunneling relaxation    
!
  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)        :: ntsed, nt
    REAL (KIND=dp), INTENT(IN)  :: omega
    REAL (KIND=dp), INTENT(IN)  :: t(ntsed)
    REAL (KIND=dp), INTENT(OUT) :: phon(ntsed)
    REAL (KIND=dp)              :: aphon
    REAL (KIND=dp)              :: x(ntsed), phon_high(ntsed)

    aphon = xpi * vt(nt)**5 * rhom(nt) * xhbar**4 / (gamma_e(nt)*everg)**2 / (2.0_dp*xkb)**3
    x(:)  = aphon * omega / t(:)**3
    phon(:) = 1.4696_dp*(1.0_dp - TANH( 0.89_dp*LOG10(4.3_dp*x(:)) ) ) / 2.0_dp    ! LV speed up: good to 10% 
    phon_high(:) = 1e-3_dp*(1e3_dp/x)
    WHERE( phon > phon_high) phon = phon_high                                      ! use power law extrapolation at high x
    phon(:) = omega * phon(:) / 4.0_dp / xpi**2

  END SUBROUTINE GET_PHON

  SUBROUTINE GET_HOP (omega, ntsed, nt, t, hop)
!
! compute opacity for hopping relaxation   
!
  ! modules with global variables
    USE CONSTANTS
    USE UTILITY

    IMPLICIT NONE

    INTEGER, INTENT (IN)        :: ntsed, nt
    REAL (KIND=dp), INTENT(IN)  :: omega
    REAL (KIND=dp), INTENT(IN)  :: t(ntsed)
    REAL (KIND=dp), INTENT(OUT) :: hop(ntsed)
    REAL (KIND=dp)              :: x, Cv, sw, sig
    INTEGER                     :: i, nV
    REAL (KIND=dp), ALLOCATABLE :: V(:), Pv(:), tau(:), f(:)

    nV = 10
    ALLOCATE (V(nV))
    ALLOCATE (Pv(nV))
    ALLOCATE (tau(nV))
    ALLOCATE (f(nV))
    x     = (Vm(nt) - Vmin(nt)) / V0(nt)
    Cv    = ( 0.5_dp * (XERF(x) + 1.0_dp) )**(-1) / V0(nt)/sqpi
    sig   = V0(nt) / SQRT(2.0_dp)
    sw = 3.0_dp   ! integrate Pv out to Vmax=sw*sigma 
    V(1)  = Vmin(nt)
    Pv(1) = Cv * exp(- (V(1)-Vm(nt))**2 / V0(nt)**2)
    DO i=2,nV
       V(i)  = V(i-1) + (sw*sig-Vmin(nt))/(REAL(nv,dp)-1.0_dp)
       Pv(i) = Cv * exp(- (V(i)-Vm(nt))**2 / V0(nt)**2)
    ENDDO
    DO i=1,ntsed
       IF (MAXVAL(V)/t(i) > 100.0_dp) THEN  !DEB: j'ai enlevé le /xkb
          hop(i) = 0.0_dp
       ELSE 
          tau(:) = tau_0(nt) * exp(V(:) / t(i))  !DEB: j'ai enlevé le /xkb
          f(:)   = omega * tau(:) * Pv(:) / ( 1.0_dp + omega**2 * tau(:)**2 )
          hop(i) = 2.0_dp/xpi * omega * (5.8_dp+c_delta(nt) + log(t(i))) * XINTEG2(1, nV, nV, V, f) ! DEB: j'ai rajouté le 5.8
       ENDIF
    ENDDO
!write(*,*)  XINTEG2(1, nV, nV, V, f)
  END SUBROUTINE GET_HOP

END MODULE MDTLS
