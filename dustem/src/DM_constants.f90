
MODULE CONSTANTS

! JLB - Sept 2009
! Gather all constantes of DUSTEM here.
! This module should be used everywhere.

! "kind PARAMETER"  for double precision

  INTEGER, PUBLIC, PARAMETER        :: dp      = SELECTED_REAL_KIND(P=15)

! usual numbers

  REAL (KIND=dp), PUBLIC, PARAMETER :: istiny  = 1.0e-60_dp
  REAL (KIND=dp), PUBLIC, PARAMETER :: hugest  = HUGE(0.0_dp)
  REAL (KIND=dp), PUBLIC, PARAMETER :: l_hu    = 300.0_dp
  REAL (KIND=dp), PUBLIC, PARAMETER :: tiniest = TINY(0.0_dp)
  REAL (KIND=dp), PUBLIC, PARAMETER :: l_ty    = -300.0_dp
  REAL (KIND=dp), PUBLIC, PARAMETER :: l_10    = 2.30258509299405_dp

! Physical constants (Partly from  rev. mod. phys. 59 , 1121 (1987) )
! Revised from physics.nist.gov/constants

  REAL (KIND=dp), PUBLIC, PARAMETER :: amu     = 1.66053873e-24_dp        ! Atomic mass unit (g)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xmh     = 1.6735e-24_dp            ! H atom mass (g)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xmp     = 1.67262158e-24_dp        ! proton mass (g)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xme     = 9.10938188e-28_dp        ! electron mass (g)
!  REAL (KIND=dp), PUBLIC, PARAMETER :: Na      = 6.0221e23_dp             ! Avogadro constant
  REAL (KIND=dp), PUBLIC, PARAMETER :: xqe     = 4.80325e-10_dp           ! electric charge (esu)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xqe2    = 2.30712e-19_dp           ! xqe * xqe
  REAL (KIND=dp), PUBLIC, PARAMETER :: e2smc2  = 2.817939e-13_dp          ! r0 = e2 / (m*c**2) (cm) - Lang
  REAL (KIND=dp), PUBLIC, PARAMETER :: xkb     = 1.3806503e-16_dp         ! Boltzman constant k (erg deg-1)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xhp     = 6.62606876e-27_dp        ! Planck constant h (erg s)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xhbar   = 1.054571628e-27_dp       ! h/2/pi
  REAL (KIND=dp), PUBLIC, PARAMETER :: dcgs    = 1.0e-18_dp               ! Debye*dcgs -> cgs
  REAL (KIND=dp), PUBLIC, PARAMETER :: sigma   = 5.6705e-5_dp             ! Stefan constant (cgs)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xr      = 8.314472e7_dp            ! Perfect gaz constant R (erg deg-1 mole-1)
  REAL (KIND=dp), PUBLIC, PARAMETER :: clight  = 2.99792458e10_dp         ! Celerity of light c (cm s-1)
  REAL (KIND=dp), PUBLIC, PARAMETER :: xpi     = 3.1415926535897932384_dp ! pi (!)
  REAL (KIND=dp), PUBLIC, PARAMETER :: sqpi    = 1.77245385090551602_dp   ! Sqrt(pi)
  REAL (KIND=dp), PUBLIC, PARAMETER :: g_euler = 0.5772156649015328606065120900824024310422_dp
  REAL (KIND=dp), PUBLIC, PARAMETER :: t_cmb   = 2.728_dp

! Conversion factors

  REAL (KIND=dp), PUBLIC, PARAMETER :: hcsurk  = 1.43883442658354_dp      ! From Wave number (cm-1) to Temperature (Kelvin)
  REAL (KIND=dp), PUBLIC, PARAMETER :: everg   = 1.6022e-12_dp            ! From Electron-Volt to Erg
  REAL (KIND=dp), PUBLIC, PARAMETER :: evtoK   = 11604.676_dp             ! eV to K
  REAL (KIND=dp), PUBLIC, PARAMETER :: calev   = 4.3363e-2_dp             ! From kilocal mole-1 to Ev
  REAL (KIND=dp), PUBLIC, PARAMETER :: pccm    = 3.086e18_dp              ! From parsec to cm
  REAL (KIND=dp), PUBLIC, PARAMETER :: ryd2ang = 911.267101235_dp         ! Rydberg to Angstrom conversion
  REAL (KIND=dp), PUBLIC, PARAMETER :: deksamu = 166289442.703935_dp      ! 2 * k / amu
  REAL (KIND=dp), PUBLIC, PARAMETER :: dehpcde = 1.19104272254325e-5_dp   ! 2 * h c**2
  REAL (KIND=dp), PUBLIC, PARAMETER :: cte1    = 0.000149670842690138_dp  ! 4 pi * dehpcde
  REAL (KIND=dp), PUBLIC, PARAMETER :: cte2    = 1.85291028578918e-46_dp  ! 4 pi * 2 * h / c**2
  REAL (KIND=dp), PUBLIC, PARAMETER :: hpc     = 1.98644544043741e-16_dp  ! h * c

! sizes of array for size, charge and temperature distributions

  INTEGER, PUBLIC, PARAMETER        :: ncar  = 40                         ! maximum length of CHARACTER defining grain type
  INTEGER, PUBLIC, PARAMETER        :: nsize_max_qabs = 100               ! max nr of size bins in data files (Q, C and G)
  INTEGER, PUBLIC, PARAMETER        :: nsize_max = 100                    ! max nr of size bins requested
  INTEGER, PUBLIC, PARAMETER        :: ntempmax = 100                     ! max nr of T bins in C files 
  INTEGER, PUBLIC, PARAMETER        :: ndist = 200                        ! nr of T bins in dP/dT
  INTEGER, PUBLIC, PARAMETER        :: nbeta_max = 50                     ! max nr of bins in BETA(T) values
  INTEGER, PUBLIC, PARAMETER        :: nz_sg = 20                         ! max nr of charge states in f(Z) (small grains)
  INTEGER, PUBLIC, PARAMETER        :: nz_bg = 200                        ! max nr of charge states in f(Z) (large grains)
  REAL (KIND=dp), PUBLIC, PARAMETER :: ztrans = 50.0_dp                   ! factor for transition between small and large grain charge (to be raised with nz_bg)
  REAL (KIND=dp), PUBLIC, PARAMETER :: fzmin = istiny                     ! min value in charge distribution

! Input and Output directories

  CHARACTER (len=100)        :: data_path='/Users/lverstra/Desktop/dustem4.2/' !!!!! TO BE CHANGED
  CHARACTER (len=108)        :: dir_DAT ='data/'
  CHARACTER (len=109)        :: dir_QABS='oprop/'
  CHARACTER (len=109)        :: dir_CAPA='hcap/'
  CHARACTER (len=108)        :: dir_RES ='out/'
  CHARACTER (len=100)        :: dir_PDR ='/Users/lverstra/Desktop/dustem4.2/' !!!!! TO BE CHANGED
  INTEGER, PUBLIC, PARAMETER :: max_len = 200

! Definition of derived types:

  TYPE :: FICH    ! structure CHARACTERizing a file
     CHARACTER(LEN=max_len) :: nom          ! file name
     INTEGER                :: unit         ! file unit
  END TYPE FICH

END MODULE CONSTANTS

