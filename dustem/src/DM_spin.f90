MODULE MSPIN
  
  USE CONSTANTS
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: SPIN, SP_RATE

CONTAINS

!----------------------------------------------------------------------------

  SUBROUTINE SPIN (nt, ns, a, t_eq, nufnu, sp_nufnu, trot)
! 
! computes the emission of spinning dust grains
!

  ! modules with global variables
    USE CONSTANTS
    USE UTILITY
    USE MGET_QEXT

    IMPLICIT NONE
 
    ! arguments
    INTEGER,        INTENT (IN)  :: nt                     ! current index of grain type
    INTEGER,        INTENT (IN)  :: ns                     ! index of grain size
    REAL (KIND=dp), INTENT (IN)  :: a                      ! grain radius (cm)
    REAL (KIND=dp), INTENT (IN)  :: t_eq                   ! equilibrium temperature (K)
    REAL (KIND=dp), INTENT (IN)  :: nufnu(n_qabs)          ! IR SED (erg/s/sr/Hz)
    REAL (KIND=dp), INTENT (OUT) :: sp_nufnu(n_qabs)       ! spinning emission nu*Inu (erg/s/H)
    REAL (KIND=dp), INTENT (OUT) :: trot                   ! rotational temperature (K)

    ! local variables
    INTEGER        :: i, j
    REAL (KIND=dp) :: omega_quad             ! mean quadratic angular frequency (rad/s)
    REAL (KIND=dp) :: mu                     ! electric dipole moment (Debye)
    REAL (KIND=dp) :: omega(n_qabs)          ! rotational frequency (rad/s)
    REAL (KIND=dp) :: xi                     ! moment of inertia (g cm^2)

    CALL SP_RATE(nt, ns, a, t_eq, mu, nufnu, xi, omega_quad)

    omega = 2.0_dp * xpi * freq_qabs
    DO i = 1,n_qabs
       j           = n_qabs-i+1
       sp_nufnu(i) = omega(j)**6 * exp(-3.0_dp/2.0_dp * omega(j)**2/omega_quad) * freq_qabs(j)
    ENDDO
    sp_nufnu = sp_nufnu / omega_quad**1.5_dp * sqrt(8.0_dp/3.0_dp/xpi) * 1.0_dp/clight**3 * (dcgs*mu)**2 * 4.0_dp*xpi   ! erg/s

    trot = xi * omega_quad / xkb / 3.0_dp

  END SUBROUTINE SPIN

!----------------------------------------------------------------------------

  SUBROUTINE SP_RATE(nt, ns, a, t_eq, mu, nufnu, xi, omega_quad)

    USE CONSTANTS
    USE UTILITY
    USE MGET_QEXT
    USE MGET_TDIST

    IMPLICIT NONE

    ! arguments
    INTEGER,        INTENT (IN)    :: nt                      ! current index of grain type
    INTEGER,        INTENT (IN)    :: ns                      ! index of grain size
    REAL (KIND=dp), INTENT (IN)    :: a                       ! grain radius
    REAL (KIND=dp), INTENT (IN)    :: t_eq                    ! equilibrium temperature (K)
    REAL (KIND=dp), INTENT (INOUT) :: mu                      ! electric dipole moment (Debye)
    REAL (KIND=dp), INTENT (IN)    :: nufnu(n_qabs)           ! IR SED in erg/s/sr/Hz
    REAL (KIND=dp), INTENT (INOUT) :: omega_quad              ! mean quadratic angular frequency (rad/s)
    REAL (KIND=dp), INTENT (OUT)   :: xi                      ! moment of inertia (g cm^2)

    ! local arguments
    INTEGER                     :: j, k
    REAL (KIND=dp)              :: cs, cx, zeta, tmp
    REAL (KIND=dp)              :: as                            ! equivalent radius for surface reactions
    REAL (KIND=dp)              :: ax                            ! equivalent radius for spinning excitation 
    REAL (KIND=dp)              :: vol_gr                        ! volume of grain
    REAL (KIND=dp)              :: mass_gr                       ! mass of grain
    REAL (KIND=dp)              :: tau_h                         ! DL98 normalization factors for spin rates
    REAL (KIND=dp)              :: tau_ed                        ! characteristic damping time
    REAL (KIND=dp)              :: f_ir, g_ir                    ! IR rates
    REAL (KIND=dp)              :: f_plasma                      ! plasma drag rates (f_plasma = g_plasma)
    REAL (KIND=dp)              :: g_H2                          ! H2 formation rate
    REAL (KIND=dp)              :: f_n, g_nin, g_nev             ! neutral gas particule/grain collision rates
    REAL (KIND=dp)              :: f_i, g_iin, g_iev             ! ions/grain collision rates
    REAL (KIND=dp)              :: f_tot                         ! total damping rate
    REAL (KIND=dp)              :: g_tot                         ! total exciting rate
    REAL (KIND=dp), ALLOCATABLE :: f(:)
    REAL (KIND=dp)              :: gamma, Ef, J2, y              ! constants for H2 excitation
    REAL (KIND=dp)              :: ld, bom, bq                   ! characteristic lenghts and angle, used for plasma drag
    REAL (KIND=dp)              :: cosi                          ! electronic density
    REAL (KIND=dp)              :: mu_tild, psi, eps_i, eps_e    ! variables for gas/grain collisions
    REAL (KIND=dp)              :: eps_n, phi, uo, t_ev, tq
    REAL (KIND=dp)              :: g1, g2, h1, h2, pol, nn, mn
    REAL (KIND=dp)              :: cst1, cst2
    REAL (KIND=dp)              :: Rabs, Pabs, Eabs              ! ISRF photons absorption rate, absorbed power,
                                                                 ! mean absorbed energy (for Tev)
    REAL (KIND=dp), ALLOCATABLE :: temperature_cal(:), calo(:)   ! heat capacity (for Tev calculation)
    REAL (KIND=dp), ALLOCATABLE :: energie_cal(:)
    INTEGER                     :: jlo, ntmp, nst

    ! defining grain shape and radius
    cs      = 1.0_dp    ! 1 for sphere, 1.0_dp/2.0_dp for thin disk 
    cx      = 1.0_dp    ! 1 for sphere, (3.0_dp/8.0_dp)**0.25_dp for thin disk 
    zeta    = 1.0_dp    ! 1 for sphere, 5.0_dp/4.0_dp for thin disk
    as      = cs * a
    ax      = cx * a 
    vol_gr  = 4.0_dp / 3.0_dp * xpi * a**3
    mass_gr = rhom(nt) * vol_gr
    xi      = zeta * 0.4_dp * mass_gr * a**2

    ! electric dipole moment
    mu = m0(nt)*sqrt(5.45e23_dp*a**3)

    ! timescale for collisions with H atoms
    tau_h = 1.0_dp / ( hden * xmh * SQRT(2.0_dp*xkb*t_gas / xpi/xmh) * ax**4 * xpi * 4.0_dp/3.0_dp / xi )
   
    ! ions polarizability
    pol_i = pol_i*(1.0e-8_dp)**3 ! from A^3 to cm^3

    ! characteristic damping time
    tau_ed = 3.0_dp * xi**2 * clight**3 / 4.0_dp / xkb / t_gas / (dcgs*mu)**2

    ! get damping and excitational rates for rotation
    f_tot = 0.0_dp
    g_tot = 0.0_dp

    ! evaporation temperature after gas/grain collision
    ALLOCATE (f(1:n_qabs))
    f    = sect_eff(nt,ns) * qi_abs(nt,ns,:) * isrfuv(:)
    Pabs = XINTEG2(1, n_qabs, n_qabs, freq_qabs, f)
    f    = f / xhp / freq_qabs
    Rabs = XINTEG2(1, n_qabs, n_qabs, freq_qabs, f)
    DEALLOCATE (f)
    Eabs = Pabs / Rabs
    ntmp = n_temp(nt)
    nst  = nsize_type(nt)
    jlo  = 1
    ALLOCATE (calo(1:ntmp))
    ALLOCATE (temperature_cal(1:ntmp))
    ALLOCATE (energie_cal(1:ntmp))
    DO j=1,ntmp
       calo(j)            = INTPOL3(calor(nt,1:nst,j), size_type(nt,1:nst), nst, a, jlo)
       calo(j)            = vol_gr * 10.0_dp**calo(j)
       temperature_cal(j) = 10.0_dp**(temp_cal(nt,j))
       energie_cal(j)     = XINTEG2(1, j, ntmp, temperature_cal, calo)
    ENDDO
    tq   = INTPOL3(temperature_cal, energie_cal, ntmp, Eabs, jlo)
    t_ev = max(t_eq, tq)
    DEALLOCATE (calo)
    DEALLOCATE (temperature_cal)
    DEALLOCATE (energie_cal)

    ! get IR rates
    f_ir  = 0.0_dp
    g_ir  = 0.0_dp
    f_ir  = 2.0_dp * tau_h / xpi / xi
    g_ir  = xhp / 6.0_dp/xpi/xi * tau_h / xkb/t_gas
    f_ir  = f_ir * XINTEG2(1, n_qabs, n_qabs, freq_qabs, nufnu/freq_qabs**2)
    g_ir  = g_ir * XINTEG2(1, n_qabs, n_qabs, freq_qabs, nufnu/freq_qabs)
    f_tot = f_tot + f_ir
    g_tot = g_tot + g_ir

    ! get plasma drag rates
    f_plasma = 0.0_dp
    cosi     = 1.0_dp / 3.0_dp
    ld       = sqrt(xkb*t_gas / 4.0_dp/xpi/eden/xqe2)
    bom      = sqrt(xi/mass_gr)
    DO j=1,nion
       bq       = sqrt(2.0_dp*xkb*t_gas/mi(j)/amu) * xi/xhbar
       tmp      = ( log(bom/as) + cosi*log(min(bq,ld)/bom) ) * (dcgs*mu)**2 * iden(j)/hden * sqrt(mi(j)*amu/xmh)
       tmp      = tmp * 2.0_dp/3.0_dp * ( Zi(j)*xqe / xkb/t_gas/ax**2 )**2
       f_plasma = f_plasma+tmp
    ENDDO
    f_tot    = f_tot + f_plasma
    g_tot    = g_tot + f_plasma

! get photoelectric rate (Weingartner et Draine 2001)
    cst1    = xme/xmh / ( 2.0_dp*xpi*as**2 * hden * sqrt(2.0_dp*xkb*t_gas/xpi/xmh) )
    cst2    = xme / 4.0_dp/ hden / sqrt(8.0_dp*xpi*xmh*xkb*t_gas) / ax**2 / xkb / t_gas
    f_tot = f_tot + jpe(ns) * cst1
    g_tot = g_tot + hspe(ns) * cst2
   
    ! get H2 formation rate
    g_H2  = 0.0_dp
    gamma = 0.1_dp / 4.0_dp
    Ef    = 0.2_dp*everg
    J2    = 100.0_dp
    y     = 2.0_dp * h2den/hden
    g_H2  = gamma * (1.0_dp-y) * Ef/xkb/t_gas
    g_H2  = g_H2 * (1.0_dp + J2 * xhbar**2 / 2.0_dp /xmh /Ef /ax**2)
    g_tot = g_tot + g_H2

    ! get gas-grain collision rates
    ! collisions with neutrals
    f_n   = 0.0_dp
    g_nin = 0.0_dp
    g_nev = 0.0_dp
    DO j=1,nZb
       ! collisions with H atom
       tmp   = 0.0_dp
       nn    = hden
       mn    = xmh
       pol   = 0.67_dp*(1.0e-8_dp)**3                                  ! from A^3 to cm^3
       eps_n = sqrt( pol*(Zb(j)*xqe/a**2)**2 / 2.0_dp/xkb/t_gas )
       eps_e = sqrt( pol*(Zb(j)*xqe/a**2)**2 / 2.0_dp/xkb/t_ev  )
       tmp   = nn/hden * sqrt(mn/xmh) * (exp(-eps_e**2) + 2.0_dp*eps_e**2)
       tmp   = tmp * (exp(-eps_n**2)+sqrt(xpi*eps_n)*xerf(eps_n)) / (exp(-eps_e**2)+sqrt(xpi*eps_e)*xerf(eps_e))
       f_n   = f_n + tmp*fZ(j)
       g_nin = g_nin + fZ(j) * nn/2.0_dp/hden * sqrt(mn/xmh) * (exp(-eps_n**2) + 2.0_dp*eps_n**2)
       g_nev = g_nev + tmp*fZ(j) * t_ev/2.0_dp/t_gas
       ! collision with H2 molecule
       tmp   = 0.0_dp
       nn    = h2den
       mn    = 2.0_dp*xmh
       pol   = 0.79_dp*(1.0e-8_dp)**3                                  ! from A^3 to cm^3
       eps_n = sqrt( pol*(Zb(j)*xqe/a**2)**2 / 2.0_dp/xkb/t_gas)
       eps_e = sqrt( pol*(Zb(j)*xqe/a**2)**2 / 2.0_dp/xkb/t_ev)
       tmp   = nn/hden * sqrt(mn/xmh) * (exp(-eps_e**2) + 2.0_dp*eps_e**2)
       tmp   = tmp * (exp(-eps_n**2)+sqrt(xpi*eps_n)*xerf(eps_n)) / (exp(-eps_e**2)+sqrt(xpi*eps_e)*xerf(eps_e))
       f_n   = f_n + tmp*fZ(j)
       g_nin = g_nin + fZ(j) * nn/2.0_dp/hden * sqrt(mn/xmh) * (exp(-eps_n**2) + 2.0_dp*eps_n**2)
       g_nev = g_nev + tmp*fZ(j) * t_ev/2.0_dp/t_gas
       ! other species...
       ! ...
    ENDDO
    f_tot = f_tot + f_n
    g_tot = g_tot + g_nin + g_nev
    ! collisions with ions
    f_i   = 0.0_dp
    g_iin = 0.0_dp
    g_iev = 0.0_dp
    DO j=1,nion 
       IF (zi(j)>0.0_dp) THEN ! only treat positive ions
          phi     = sqrt( 2.0_dp * (Zi(j)*xqe)**2 / ax/xkb/t_gas )
          mu_tild = Zi(j) * xqe * dcgs*mu / ax**2 / xkb / t_gas
          uo      = ( -phi + sqrt(phi**2+4.0_dp*mu_tild) ) / 2.0_dp
          DO k=1,nZb
             h1    = 0.0_dp
             h2    = 0.0_dp
             g1    = 0.0_dp
             g2    = 0.0_dp
             psi   = Zb(k)*Zi(j)*xqe**2 / ax/xkb/t_gas
             eps_i = sqrt( pol_i(j) * (Zb(k)*xqe/a**2)**2 / 2.0_dp/xkb/t_ev )
             IF (Zb(k) == 0.0_dp) THEN
                h1    = 0.5_dp + mu_tild/4.0_dp + (2.0_dp+phi**2)*(1.0_dp-exp(-uo**2)) / 4.0_dp/mu_tild
                h1    = h1 - phi*uo*exp(-uo**2) / 4.0_dp/mu_tild
                h1    = h1 + sqpi*phi/2.0_dp * (1.0_dp+(3.0_dp-2.0_dp*mu_tild)*xerf(uo)/4.0_dp/mu_tild)
                h2    = 0.5_dp + 3.0_dp*sqpi*phi/4.0_dp + phi**2/4.0_dp + mu_tild**2/12.0_dp + mu_tild/4.0_dp
                h2    = h2 + (1.0_dp+phi**2)*(1.0_dp-exp(-uo**2))/2.0_dp/mu_tild
                h2    = h2 + (2.0_dp*mu_tild*phi**2+phi*(2.0_dp*mu_tild-7.0_dp)*uo)*exp(-uo**2) / 16.0_dp/mu_tild
                h2    = h2 + sqpi*phi*xerf(uo)*(4.0_dp*mu_tild**2-12.0_dp*mu_tild+15.0_dp+2.0_dp*phi**2)/32.0_dp/mu_tild
                f_i   = f_i + fZ(k) * iden(j)/hden * sqrt(mi(j)*amu/xmh) * h1
                g_iin = g_iin + fZ(k) * iden(j)/2.0_dp/hden * sqrt(mi(j)*amu/xmh) * h2
                g_iev = g_iev + fZ(k) * iden(j)/hden * sqrt(mi(j)*amu/xmh) * h1 * t_ev/2.0_dp/t_gas
             ELSE
                IF (mu_tild <= abs(psi)) THEN
                   IF (psi < 0.0_dp) THEN 
                      g1 = 1.0_dp - psi
                      g2 = 1.0_dp - psi + psi**2/2.0_dp + mu_tild**2/6.0_dp
                   ELSE
                      g1 = exp(-psi) * sinh(1.0_dp)
                      g2 = g1
                   ENDIF
                ELSE
                   g1 = ( 1.0_dp - exp(-psi-mu_tild) + mu_tild - psi + 0.5_dp*(mu_tild-psi)**2 ) / 2.0_dp/mu_tild
                   g2 = g1 + (mu_tild-psi)**3 / 12.0_dp/mu_tild
                ENDIF
                tmp   = iden(j)/hden * sqrt(mi(j)*amu/xmh) * g1
                tmp   = tmp * (exp(-eps_i**2) + 2.0_dp*eps_i**2) / (exp(-eps_i**2) + sqpi*eps_i*xerf(eps_i))
                f_i   = f_i + tmp*fZ(k)
                g_iin = g_iin + fZ(k) * iden(j)/2.0_dp/hden * sqrt(mi(j)*amu/xmh) * g2
                g_iev = g_iev + tmp*fZ(k) * t_ev/2.0_dp/t_gas
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    f_tot = f_tot + f_i
    g_tot = g_tot + g_iin + g_iev

    omega_quad = 20.0_dp/3.0_dp * g_tot/f_tot**2 * tau_h/tau_ed
    omega_quad = 2.0_dp / ( 1.0_dp + sqrt(1.0_dp + omega_quad) ) * g_tot/f_tot * 3.0_dp*xkb*t_gas/xi

  END SUBROUTINE SP_RATE

END MODULE MSPIN
