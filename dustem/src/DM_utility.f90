
MODULE UTILITY

  USE CONSTANTS
  IMPLICIT NONE

! variables for radiation field
  INTEGER, PUBLIC                            :: nisrf            ! nr of wave points in ISRF.DAT
  REAL (KIND=dp), PUBLIC                     :: hnumin, hnumax   ! Max and min energy of incident photons
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: isrfuv(:)        ! radiation field interpolated on lambda_qabs grid (from LAMBDA.DAT)

! variables for grain type
  CHARACTER (LEN=20),PUBLIC, ALLOCATABLE     :: sbulk(:)         ! string for bulk materials
  INTEGER, PUBLIC, ALLOCATABLE               :: nsize(:)         ! nr of size bins for each grain type
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: masstot(:)       ! Mass in grains of type i

! size dependent quantities for each grain type
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: size_ava(:,:)    ! sizes for each grain type
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: si_ava_l(:,:)    ! Log of sizes for each grain type
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: sect_eff(:,:)    ! grain cross section
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: ava(:,:)         ! grain volume distribution a^4*dn/da
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: dloga(:,:)       ! size step in log scale
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: rho(:,:)         ! mass density (g/cm3) per species and per size
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: f_mix(:,:)       ! mixing factor for SED
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: f_pol(:,:)       ! alignment efficiency for polarization (type,size)

! variables for GRAIN.DAT
  INTEGER, PUBLIC                            :: ntype            ! nr of grain type (for allocatables below)

! global keywords
  INTEGER, PUBLIC                            :: n_ftemp          ! for temperature distribution
  INTEGER, PUBLIC                            :: n_fsize          ! to use SIZE_*.DAT files
  INTEGER, PUBLIC                            :: n_quiet          ! verbose off
  INTEGER, PUBLIC                            :: n_res_a          ! for size resolved SED output
  INTEGER, PUBLIC                            :: n_sdist          ! for size distribution output
  INTEGER, PUBLIC                            :: n_zdist          ! for charge distribution output
  INTEGER, PUBLIC                            :: n_pdr_o          ! Write specific outputs for Meudon PDR code

! type keywords
  INTEGER, PUBLIC                            :: n_chrg           ! charge distribution
  INTEGER, PUBLIC                            :: n_zm             ! mix & charge distribution (PAH)
  INTEGER, PUBLIC                            :: n_spin           ! spinning dust
  INTEGER, PUBLIC                            :: n_beta           ! beta-correction
  INTEGER, PUBLIC                            :: n_pol            ! polarization 
  INTEGER, PUBLIC                            :: n_dtls           ! DCD/TLS effects

  ! Polarization
  INTEGER, PUBLIC                            :: n_lin            ! Linear polarization
  INTEGER, PUBLIC                            :: n_univ           ! Universal alignment law
  INTEGER, PUBLIC                            :: n_circ           ! Circular polarization
  INTEGER, PUBLIC                            :: n_anis           ! Take into account anisotropy (i.e. dust alignment) in emission and extinction
  INTEGER, PUBLIC                            :: n_rrf            ! The alignment function is the RRF : interpolate Q coefficients on RRF files
  INTEGER, PUBLIC                            :: n_falig          ! to force to use POL_*.DAT files

  REAL (KIND=dp), PUBLIC                     :: g0               ! radiation field scaling factor
  REAL (KIND=dp), PUBLIC                     :: cr_rate          ! cosmic ray rate
  REAL (KIND=dp), PUBLIC                     :: anisg0           ! anisotropy factor for the radiation field 
  CHARACTER (LEN=NCAR), PUBLIC, ALLOCATABLE  :: gtype(:)         ! string for grain type
  CHARACTER (LEN=20), PUBLIC, ALLOCATABLE    :: t_opt(:)         ! options for grain TYPE (PLAW, LOGN, SIZE, MIX, POL, SPIN, DTLS)
  CHARACTER (LEN=20), PUBLIC, ALLOCATABLE    :: p_opt(:)         ! options for alignment function (IDG, RAT, PAR)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: mprop(:)         ! dust-to-gas mass ratio (grain mass)/(nH*mH)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: rhom(:)          ! mean (over size) density of grain (g/cm3)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: as1(:)           ! for power law size dist: min sizes
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: as2(:)           !                        : max sizes
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: sd_par(:)        ! other parameters for size dist: alpha-a0,sig and 3 ED or CV parameters
                                                                 ! sd_par(1)  : power law index or logn centroid
                                                                 ! sd_par(2)  : logn width
                                                                 ! then parameters IN ORDER for exp. decay (ED) THEN for curvature (CV) (3 param each)
                                                                 ! sd_par(3:5): at, ac and gam for ED    EXP(-((a-at)/ac)**gam)
                                                                 ! sd_par(6:8): at, delta and gam for CV  (1+ABS(beta)*(a/at)**gam)**SGN(delta)

! variables for GAS.DAT file
  REAL (KIND=dp), PUBLIC                     :: hden             ! gas proton density (cm-3)
  REAL (KIND=dp), PUBLIC                     :: t_gas            ! gas temperature (K)
  REAL (KIND=dp), PUBLIC                     :: h2den            ! gas H2 density (cm-3)
  INTEGER, PUBLIC                            :: nion             ! nr of ion type involved in gas/grain interaction
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: iden(:)          ! ion density (cm-3)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: mi(:)            ! ion mass (amu)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: pol_i(:)         ! ion polarizability (angstrom**3)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: Zi(:)            ! ion charge
  REAL (KIND=dp), PUBLIC                     :: eden             ! electronic density (cm-3)

! variables for CHRG_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: wf(:)               ! work function for PE effect - array(ntype)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: ebg(:,:)            ! band gap - array(ntype, nsize_max)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: p_ea(:,:)           ! for EA correction - array(ntype,2)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: s_ea(:,:)           ! for EA cross-section - array(ntype,2)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: p_y(:,:),le(:)      ! for PE yield - array(ntype)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: p_uait(:,:)         ! for Uait - array(ntype)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: m1mass(:)           ! molecular mass of grain - array(ntype)

! variables for SPIN_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: m0(:)                 ! array(ntype) factor for the grain electric dipole moment (Debye)

! variables for BETA_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: beta0(:)              ! standard value of submm beta 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: abeta(:)              ! parameters for beta function, DBETA(T)=-beta0+abeta*T**gbeta
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: gbeta(:)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: bmax(:)               ! max value of beta
  INTEGER,        PUBLIC, ALLOCATABLE        :: nbeta(:)              ! nr of read beta values
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: tbeta(:,:),betav(:,:) ! T, BETA(T) read in BETA_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: f_beta(:,:)           ! threshold array(ntype,n_qabs) for beta-correction 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: ltresh(:)             ! threshold (microns) for beta-correction 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: lstiff(:)             ! stiffness (relative to ltresh) of threshold for beta-correction
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: f_dtls(:,:)           ! threshold array(ntype,n_qabs) for DTLS correction 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: ldtresh(:)            ! threshold (microns) for DTLS correction
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: ldstiff(:)            ! stiffness (relative to ltresh) of threshold for DTLS correction
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qdtls(:,:)            ! ref value @ ldtresh (type,size) to apply DTLS correction

! variables for POL_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: athresh(:)             ! threshold (microns) for polarization
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: pstiff(:)             ! stiffness (relative to athresh) of threshold
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: plev(:)               ! level after threshold

! variables for DTLS_*.DAT file
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: vt(:)                 ! sound transversal speed
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: lc(:)                 ! correlation length (nm)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: a_dtls(:)             ! relative weight of DCD to TLS in Qabs
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: Pmu(:)                ! 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: gamma_e(:)            ! elastic dipole moment (eV)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: omega_m(:)            ! threshold frequency for resonant absorption
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: c_delta(:)            ! 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: tau_0(:)              ! relaxation time for tunneling relaxation
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: V0(:)                 ! parameters to define the potential barrier
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: Vmin(:)               ! (hopping term)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: Vm(:)                 !

! variables for dP/dT
  INTEGER, PUBLIC                            :: nit              ! nr of iterations to compute dP/dT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: fdist(:,:)       ! kernel to get dP/dT

! Variables for Q_ext and Q_abs
  INTEGER, PUBLIC, ALLOCATABLE               :: nsize_type(:)    ! array(ntype): for each grain type nr of sizes read in QCG files 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: size_type(:,:)   ! array(ntype,nsize_max_qabs): sizes read in QCG files
  INTEGER, PUBLIC, ALLOCATABLE               :: nsz1(:)          ! same as nsize_type for reading and checking
  INTEGER, PUBLIC                            :: n_qabs           ! nr of wavelength values read in LAMBDA.DAT
  INTEGER, PUBLIC                            :: jfreqmax         ! in radiation field, index of max frequency where flux is not 0
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: lamb_qabs(:)     ! common wavelength grid for all Q's. Read in LAMBDA.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: freq_qabs(:)     ! frequencies corresponding to lambda_qabs
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: lfrq_qabs(:)     ! Log(freq_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q_abs(:,:,:)     ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in Q_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qdiff(:,:,:)     ! array(ntype,nsize_max_qabs,n_qabs) of Qsca as read in Q_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q1_abs(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in Q1_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qH_abs(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in QH_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qH1_abs(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in QH_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qH2_abs(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in QH_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q1diff(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qsca as read in Q1_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q2_abs(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qabs as read in Q2_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q2diff(:,:,:)    ! array(ntype,nsize_max_qabs,n_qabs) of Qsca as read in Q2_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qcirc(:,:,:)     ! array(ntype,nsize_max_qabs,n_qabs) of Qsca as read in Qc_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: gfac(:,:,:)      ! array(ntype,nsize_max_qabs,n_qabs) of g factors as read in G_*.DAT
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qi_abs(:,:,:)    ! Qabs interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qidiff(:,:,:)    ! Qsca interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qiH_abs(:,:,:)   ! Qabs interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q1i_abs(:,:,:)   ! Q1abs interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q1idiff(:,:,:)   ! Q1sca interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q2i_abs(:,:,:)   ! Q2abs interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: q2idiff(:,:,:)   ! Q2sca interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qicirc(:,:,:)    ! Qpha interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: gifac(:,:,:)     ! g factor interpolated on size_ava: array(ntype,max(nsize),n_qabs)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qauv(:)          ! Qabs used in GET_TDIST to calculate heating
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: qaem(:)          ! Qabs used in GET_TDIST and COOLING to calculate cooling

! heat capacities
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: calor(:,:,:)     ! array(ntype, NSIZE_MAX_QABS, ntempmax) of heat capacities (erg/K/cm3)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: temp_cal(:,:)    ! array(nt, ntempmax) temperature grid 
  INTEGER, PUBLIC, ALLOCATABLE               :: n_temp(:)        ! nr of temperature points for each type

! Emitted spectrum
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: nuinuem(:,:)     ! emitted spectrum nu*Inu per type 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: nuinuemtot(:)    ! total emitted spectrum (erg/s/H)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: spnuinuem(:,:)   ! emitted spinning spectrum nu*Inu per type 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: spnuinuemtot(:)  ! total emitted spinning spectrum (erg/s/H)
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: nuinuemp(:,:)    ! polarized emitted spectrum nu*Inu per type 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE        :: nuinuemptot(:)   ! total polarized emitted spectrum (erg/s/H)

! for charge distribution
  INTEGER, PUBLIC                     :: nzb                      ! nr of charge bins for zb and fz
  REAL (KIND=dp), PUBLIC              :: zmin,zmax,zeq            ! min, max and equilibrium grain charge
  REAL (KIND=dp), PUBLIC              :: sigp, sigm, sigpe        ! cross-section for charging processes
  REAL (KIND=dp), PUBLIC, ALLOCATABLE :: zmean(:), sd_z(:)        ! mean grain charge and standard deviation
  REAL (KIND=dp), PUBLIC, ALLOCATABLE :: zb(:),fz(:)              ! grain charge and distribution per size and type 
  REAL (KIND=dp), PUBLIC, ALLOCATABLE :: jpe(:), hpe(:), hspe(:)  ! PE rate, heating and G coeff for spin as a function of size
  REAL (KIND=dp), PUBLIC, ALLOCATABLE :: jgas(:), cgas(:)         ! gas rate and cooling (electrons + ions) as a function of size

  PRIVATE

  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE
  INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8 

  PUBLIC :: TRIMCAT, READ_COM, PARSE, UPCASE, PRIMITIV2, XINTEG2, INTPOL, INTPOL2, INTPOL3, F_BB, G_BB, XERF

CONTAINS

!------------------------ USEFUL SUBROUTINES and FUNCTIONS ---------------------

FUNCTION TRIMCAT(ch1, ch2)
! concatenate character strings ch1 & ch2

  USE CONSTANTS

  IMPLICIT NONE

  CHARACTER (LEN=max_len)        :: TRIMCAT

  CHARACTER (LEN=*), INTENT (IN) :: ch1
  CHARACTER (LEN=*), INTENT (IN) :: ch2

  INTEGER                        :: i1, i2
  CHARACTER (LEN=max_len)        :: lch1, lch2, ch3

  lch1 = ADJUSTL(ch1)
  lch2 = ADJUSTL(ch2)
  i1 = LEN_TRIM(lch1)
  i2 = LEN_TRIM(lch2)

  ch3(1:i1) = TRIM(lch1)
  ch3(i1+1:i1+i2) = TRIM(lch2)
  ch3(i1+i2+1:max_len) = REPEAT(' ',max_len-i1-i2)
  TRIMCAT = TRIM(ch3)

END FUNCTION TRIMCAT

!----------------------------------------------------------------

 FUNCTION READ_COM(file_input) RESULT(line)

! reads comment lines beginning with #
! returns 1st line that is not a comment as
! CHARACTER

  USE CONSTANTS

  IMPLICIT NONE

  INTEGER, INTENT (IN)    :: file_input
  LOGICAL                 :: comment_found
  CHARACTER (LEN=max_len) :: line

  ! -----------------------------------------
  ! Skip comments
  ! -----------------------------------------
  comment_found = .TRUE.
  DO WHILE (comment_found)
     READ (file_input, FMT='(A200)') line
     IF (line(1:1) /= '#') comment_found = .FALSE.
  ENDDO

END FUNCTION READ_COM

!----------------------------------------------------------------

SUBROUTINE COMPACT(str)

! Converts multiple spaces and tabs to single spaces; deletes control characters;
! removes initial spaces.

  character(len=*):: str
  character(len=1):: ch
  character(len=len_trim(str)):: outstr
  integer         :: i,ich,isp,k,lenstr

  str=adjustl(str)
  lenstr=len_trim(str)
  outstr=' '
  isp=0
  k=0

  do i=1,lenstr
     ch=str(i:i)
     ich=iachar(ch)
  
     select case(ich)
  
     case(9,32)     ! space or tab character
        if(isp==0) then
           k=k+1
           outstr(k:k)=' '
        end if
        isp=1
      
     case(33:)      ! not a space, quote, or control character
        k=k+1
        outstr(k:k)=ch
        isp=0
        
     end select
  
  end do

  str=adjustl(outstr)

END SUBROUTINE COMPACT

!----------------------------------------------------------------

SUBROUTINE SPLIT(str,delims,before,sep)

! Routine finds the first instance of a character from 'delims' in the
! the string 'str'. The characters before the found delimiter are
! output in 'before'. The characters after the found delimiter are
! output in 'str'. The optional output character 'sep' contains the 
! found delimiter. A delimiter in 'str' is treated like an ordinary 
! character if it is preceded by a backslash (\). If the backslash 
! character is desired in 'str', then precede it with another backslash.

  character(len=*) :: str,delims,before
  character,optional :: sep
  logical :: pres
  character :: ch,cha
  integer :: i,ibsl,iposa,ipos,k,lenstr
  
  pres=present(sep)
  str=adjustl(str)
  CALL COMPACT(str)
  lenstr=len_trim(str)
  if(lenstr == 0) return        ! string str is empty
  k=0
  ibsl=0                        ! backslash initially inactive
  before=' '
  do i=1,lenstr
     ch=str(i:i)
     if(ibsl == 1) then          ! backslash active
        k=k+1
        before(k:k)=ch
        ibsl=0
        cycle
     end if
     if(ch == '\') then          ! backslash with backslash inactive
        k=k+1
        before(k:k)=ch
        ibsl=1
        cycle
     end if
     ipos=index(delims,ch)         
     if(ipos == 0) then          ! character is not a delimiter
        k=k+1
        before(k:k)=ch
        cycle
     end if
     if(ch /= ' ') then          ! character is a delimiter that is not a space
        str=str(i+1:)
        if(pres) sep=ch
        exit
     end if
     cha=str(i+1:i+1)            ! character is a space delimiter
     iposa=index(delims,cha)
     if(iposa > 0) then          ! next character is a delimiter
        str=str(i+2:)
        if(pres) sep=cha
        exit
     else
        str=str(i+1:)
        if(pres) sep=ch
        exit
     end if
  end do
  if(i >= lenstr) str=''
  str=adjustl(str)              ! remove initial spaces
  return
  
END SUBROUTINE SPLIT

!----------------------------------------------------------------

SUBROUTINE REMOVEBKSL(str)

! Removes backslash (\) characters. Double backslashes (\\) are replaced
! by a single backslash.

  character (len=*)             :: str
  character (len=1)             :: ch
  character (len=len_trim(str)) :: outstr
  integer                       :: i,ibsl,k,lenstr

  str=adjustl(str)
  lenstr=len_trim(str)
  outstr=' '
  k=0
  ibsl=0                        ! backslash initially inactive

  do i=1,lenstr
     ch=str(i:i)
     if(ibsl == 1) then          ! backslash active
        k=k+1
        outstr(k:k)=ch
        ibsl=0
        cycle
     end if
     if(ch == '\') then          ! backslash with backslash inactive
       ibsl=1
       cycle
     end if
     k=k+1
     outstr(k:k)=ch              ! non-backslash with backslash inactive
  end do

  str=adjustl(outstr)

END SUBROUTINE REMOVEBKSL

!----------------------------------------------------------------

SUBROUTINE PARSE(str,delims,args,nargs)

! Parses the string 'str' into arguments args(1), ..., args(nargs) based on
! the delimiters contained in the string 'delims'. Preceding a delimiter in
! 'str' by a backslash (\) makes this particular instance not a delimiter.
! The integer output variable nargs contains the number of arguments found.

  character(len=*)              :: str,delims
  character(len=len_trim(str))  :: strsav
  character(len=*),dimension(:) :: args
  integer                       :: i, k, lenstr, na, nargs
  
  strsav=str
  CALL COMPACT(str)
  na=size(args)
  do i=1,na
     args(i)=' '
  end do
  nargs=0
  lenstr=len_trim(str)
  if(lenstr==0) return
  k=0
  
  do
     if(len_trim(str) == 0) exit
     nargs=nargs+1
     CALL SPLIT(str,delims,args(nargs))
     CALL REMOVEBKSL(args(nargs))
  end do
  str=strsav
  
END SUBROUTINE PARSE

!----------------------------------------------------------------
 
FUNCTION UPCASE(string) RESULT(upper)
  CHARACTER(LEN=*), INTENT(IN) :: string
  CHARACTER(LEN=len(string))   :: upper
  INTEGER                      :: j

  DO j = 1,len(string)
     IF(string(j:j) >= "a" .AND. string(j:j) <= "z") THEN
        upper(j:j) = achar(iachar(string(j:j)) - 32)
     ELSE
        upper(j:j) = string(j:j)
     ENDIF
  ENDDO
 END FUNCTION UPCASE

!----------------------------------------------------------------

SUBROUTINE PRIMITIV2 (init, n, xcoor, fsub, primit)
! computes primitive of fsub (step not constant)

  USE CONSTANTS

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: init
  INTEGER, INTENT (IN)         :: n
  REAL (KIND=dp), INTENT (IN)  :: xcoor(n)
  REAL (KIND=dp), INTENT (IN)  :: fsub(n)
  REAL (KIND=dp), INTENT (OUT) :: primit(n)
  REAL (KIND=dp)               :: primitaux
  REAL (KIND=dp)               :: xa, xb, ya, yb
  INTEGER                      :: i

  xa = xcoor(1)
  ya = fsub(1)
  primitaux = 0.0_dp
  primit(1) = primitaux

  DO i=2,n
     xb = xcoor(i)
     yb = fsub(i)
     primitaux = primitaux + (xb - xa) * (yb + ya) * 0.5_dp     ! Trapeze
     primit(i) = primitaux
     xa = xb
     ya = yb
  ENDDO

  ! Recentrage sur init
  IF (init > 0) THEN
     primitaux = primit(init)
     primit(:) = primit(:) - primitaux
  ENDIF

END SUBROUTINE PRIMITIV2

!----------------------------------------------------------------

FUNCTION XINTEG2(imin, imax, n, xin, yin)
! computes integral of yin, variable step 
! make sure you have at least two points of integration

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp) :: XINTEG2

  INTEGER, INTENT (IN)        :: imin
  INTEGER, INTENT (IN)        :: imax
  INTEGER, INTENT (IN)        :: n
  REAL (KIND=dp), INTENT (IN) :: xin(n)
  REAL (KIND=dp), INTENT (IN) :: yin(n)
  REAL (KIND=dp)              :: xa, ya, xb, yb
  REAL (KIND=dp)              :: primitaux
  INTEGER                     :: i

  xa = xin(imin)
  ya = yin(imin)
  primitaux = 0.0_dp
  DO i=imin+1,imax
     xb = xin(i)
     yb = yin(i)
     primitaux = primitaux + (xb - xa) * (ya + yb)
     xa = xb
     ya = yb
  ENDDO
  XINTEG2 = primitaux * 0.5_dp

END FUNCTION XINTEG2

!----------------------------------------------------------------

FUNCTION INTPOL(fint, xint, ni, xess)
! linear interoplation of fint @ xess
! xint assumed increasing, extrapolation set to 0
! make sure xess belongs to [xint(1),xint(ni)]

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp) :: INTPOL

  INTEGER, INTENT (IN)        :: ni
  REAL (KIND=dp), INTENT (IN) :: fint(ni)
  REAL (KIND=dp), INTENT (IN) :: xint(ni)
  REAL (KIND=dp), INTENT (IN) :: xess
  INTEGER                     :: i

  i = 1
  DO WHILE (xint(i) <= xess .AND. i /= ni)
     i = i + 1
  ENDDO
  i = i - 1

! extrapolation set to 0
  IF( (xess < xint(1)) .OR. (xess > xint(ni)) ) THEN
     INTPOL = 0.0_dp
  ELSE
     INTPOL = fint(i) * (xint(i+1) - xess) + fint(i+1) * (xess-xint(i))
     INTPOL = INTPOL / (xint(i+1) - xint(i))
  ENDIF

END FUNCTION INTPOL

!----------------------------------------------------------------

FUNCTION INTPOL2(fint, xint, ni, xess)
! linear interoplation of fint @ xess
! xint assumed increasing
! extrapolation is constant at value at edges
! make sure xess belongs to [xint(1),xint(ni)]

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp) :: INTPOL2

  INTEGER, INTENT (IN)        :: ni
  REAL (KIND=dp), INTENT (IN) :: fint(ni)
  REAL (KIND=dp), INTENT (IN) :: xint(ni)
  REAL (KIND=dp), INTENT (IN) :: xess
  INTEGER                     :: i, j

  IF(xess <= xint(1)) THEN
     INTPOL2 = fint(1)
  ELSE IF(xess >= xint(ni)) THEN
     INTPOL2 = fint(ni)
  ELSE
     DO j=2,ni
        IF (xint(j) > xess) THEN
           i = j-1
           EXIT
        ENDIF
     ENDDO
     INTPOL2 = fint(i) * (xint(i+1) - xess) + fint(i+1) * (xess-xint(i))
     INTPOL2 = INTPOL2 / (xint(i+1) - xint(i))
  ENDIF

END FUNCTION INTPOL2

!----------------------------------------------------------------

FUNCTION INTPOL3(fint, xint, ni, xess, jlo)
! linear interoplation of fint @ xess
! xint assumed increasing, uses HUNT
! extrapolation is constant at value at edges

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp) :: INTPOL3

  INTEGER, INTENT (IN)        :: ni
  REAL (KIND=dp), INTENT (IN) :: fint(ni)
  REAL (KIND=dp), INTENT (IN) :: xint(ni)
  REAL (KIND=dp), INTENT (IN) :: xess
  INTEGER, INTENT (INOUT)     :: jlo

  CALL HUNT(xint, ni, xess, jlo)

  IF(jlo == 0) THEN
     INTPOL3 = fint(1)
  ELSE IF(jlo >= ni) THEN
     INTPOL3 = fint(ni)
  ELSE
     INTPOL3 = (fint(jlo) * (xint(jlo+1) - xess) + fint(jlo+1) * (xess-xint(jlo))) &
             / (xint(jlo+1) - xint(jlo))
  ENDIF

END FUNCTION INTPOL3

!----------------------------------------------------------------

FUNCTION F_BB(x)
! for Planck function

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp)              :: F_BB
  REAL (KIND=dp), INTENT (IN) :: x

  IF (x > l_hu) THEN   ! Double precision (use 70.0 for single)
    F_BB = 0.0_dp
  ELSE
    F_BB = 1.0_dp / (EXP(x) - 1.0_dp)
  ENDIF

END FUNCTION F_BB

!----------------------------------------------------------------

FUNCTION G_BB(x)
! for derivative of Planck function

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp)              :: G_BB
  REAL (KIND=dp), INTENT (IN) :: x

  IF (x > l_hu) THEN   ! Double precision (use 70.0 for single)
    G_BB = 0.0_dp
  ELSE
    G_BB = 1.0_dp / (COSH(x) - 1.0_dp)
  ENDIF

END FUNCTION G_BB

!----------------------------------------------------------------

! Adapted from "Numerical Recipes" - JLB IX 09
SUBROUTINE HUNT(xx, n, x, jlo)

  USE CONSTANTS

  IMPLICIT NONE

  REAL (KIND=dp), INTENT (IN) :: xx(:)            !  xx must be in ascending order
  INTEGER, INTENT (IN)        :: n
  REAL (KIND=dp), INTENT (IN) :: x
  INTEGER, INTENT (INOUT)     :: jlo

  INTEGER                     :: inc, jhi, jm
  LOGICAL                     :: dicho

  dicho = .FALSE.
  IF (jlo <= 0 .OR. jlo > n) THEN
     jlo = 0
     jhi = n + 1
     dicho = .TRUE.
  ENDIF

  inc = 1
  IF (dicho .EQV. .FALSE.) THEN
     IF (x >= xx(jlo)) THEN
        DO
           jhi = jlo + inc
           IF (jhi > n) THEN
              jhi = n + 1
              EXIT
           ELSE IF (x >= xx(jhi)) THEN
              jlo = jhi
              inc = inc + inc
           ELSE
              EXIT
           ENDIF
        ENDDO
     ELSE
        jhi = jlo
        DO
           jlo = jhi - inc
           IF (jlo < 1) THEN
              jlo = 0
              EXIT
           ELSE IF (x < xx(jlo)) THEN
              jhi = jlo
              inc = inc + inc
           ELSE
              EXIT
           ENDIF
        ENDDO
     ENDIF
  ENDIF

  IF (jlo /= 0 .AND. jhi /= n+1) THEN
     IF (x < xx(jlo) .OR. x >= xx(jhi)) THEN
        print *, xx(jlo), x, xx(jhi)
     ENDIF
  ELSE
!    print *, jlo, jhi, n
  ENDIF

  DO
     IF (jhi - jlo == 1) THEN
        EXIT
     ENDIF
     jm = (jhi + jlo) / 2
     IF (x > xx(jm)) THEN
        jlo = jm
     ELSE
        jhi = jm
     ENDIF
  ENDDO

END SUBROUTINE HUNT

!----------------------------------------------------------------

! Gauss error function (and associated functions) - N. Ysard Feb 2010
! adapted from Numerical Recipes
! NB depending on compiler, ERF is not always defined as intrinsic function

FUNCTION XERF(x)
        
  IMPLICIT None
  
  REAL (KIND=dp)             :: XERF
  REAL (KIND=dp), INTENT(in) :: x

  XERF = gammp_s(0.5_dp,x**2.0_dp)
  IF (x < 0.0_dp) XERF = -XERF  

END FUNCTION XERF

FUNCTION gammp_s(a,x)
  
  IMPLICIT None
  
  REAL (KIND=dp), INTENT(in) :: a,x
  REAL (KIND=dp)             :: gammp_s
  
  IF (x < a+1.0_dp) THEN 
     gammp_s = gser_s(a,x)
  ELSE
     gammp_s = 1.0_dp-gcf_s(a,x)
  END IF
  
END FUNCTION gammp_s

FUNCTION gser_s(a,x,gln)
  
  IMPLICIT None
  
  REAL (KIND=dp), INTENT(in)            :: a,x
  REAL (KIND=dp), OPTIONAL, INTENT(out) :: gln
  REAL (KIND=dp)                        :: gser_s
  INTEGER, PARAMETER                    :: itmax=100
  REAL (KIND=dp)                        :: EPS=epsilon(x)
  INTEGER                               :: n
  REAL (KIND=dp)                        :: ap,del,summ
  
  IF (x == 0.0_dp) THEN
     gser_s = 0.0_dp
     RETURN
  END IF
  
  ap   = a
  summ = 1.0_dp/a
  del  = summ
  
  DO n = 1,itmax
     ap   = ap+1.0_dp
     del  = del*x/ap
     summ = summ+del
     IF (abs(del) < abs(summ)*EPS) EXIT
  END DO
  
  IF (n > itmax) WRITE(*,*) 'Problem in gser_s, bad result'
  
  IF (present(gln)) THEN
     gln    = gammln_s(a)
     gser_s = summ*exp(-x+a*log(x)-gln)
  ELSE
     gser_s = summ*exp(-x+a*log(x)-gammln_s(a))
  END IF
  
END FUNCTION gser_s

FUNCTION gcf_s(a,x,gln)
  
  IMPLICIT None
  
  REAL (KIND=dp), INTENT(in)            :: a,x
  REAL (KIND=dp), OPTIONAL, INTENT(out) :: gln
  REAL (KIND=dp)                        :: gcf_s
  INTEGER, PARAMETER                    :: ITMAX=100
  REAL (KIND=dp), PARAMETER             :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER                               :: i
  REAL (KIND=dp)                        :: an,b,c,d,del,h
   
  IF (x == 0.0_dp) THEN
     gcf_s = 1.0_dp
     RETURN
  END IF
  
  b = x+1.0_dp-a
  c = 1.0_dp/FPMIN
  d = 1.0_dp/b
  h = d
  
  DO i = 1,ITMAX
     an = -DBLE(i)*(DBLE(i)-a)
     b  = b+2.0_dp
     d  = an*d+b
     IF (ABS(d) < FPMIN) d = FPMIN
     c = b+an/c
     IF (ABS(c) < FPMIN) c = FPMIN
     d   = 1.0_dp/d
     del = d*c
     h   = h*del
     IF (ABS(del-1.0_dp) <= EPS) EXIT
  END DO
  
  IF (i > ITMAX) WRITE(*,*) 'Problem in gcf_s, bad result'

  IF (present(gln)) THEN
     gln   = gammln_s(a)
     gcf_s = exp(-x+a*log(x)-gln)*h
  ELSE
     gcf_s = exp(-x+a*log(x)-gammln_s(a))*h
  ENDIF
  
END FUNCTION gcf_s

FUNCTION gammln_s(xx)
  
  IMPLICIT None
  
  REAL (KIND=dp), INTENT(in)   :: xx
  REAL (KIND=dp)               :: gammln_s
  REAL (KIND=dp)               :: tmp,x
  REAL (KIND=dp)               :: stp=2.5066282746310005_dp
  REAL (KIND=dp), DIMENSION(6) :: coeff=(/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,               &
       -1.231739572450155_dp,0.1208650973866179e-2_dp,           &
       -0.5395239384953e-5_dp/)
  
  x        = xx
  tmp      = x+5.5_dp
  tmp      = (x+0.5_dp)*LOG(tmp)-tmp
  gammln_s = tmp+LOG(stp*(1.000000000190015_dp+SUM(coeff(:)/arth(x+1.0_dp,1.0_dp,SIZE(coeff)) ))/x)
  
END FUNCTION gammln_s

FUNCTION arth_r(first,increment,n)
  REAL (KIND=4), INTENT(IN)   :: first, increment
  INTEGER, INTENT(IN)          :: n
  REAL (KIND=4), DIMENSION(n) :: arth_r
  INTEGER                      :: k, k2
  REAL (KIND=4)               :: temp
  IF (n > 0) arth_r(1) = first
  IF (n <= NPAR_ARTH) THEN
     DO k = 2,n
        arth_r(k) = arth_r(k-1)+increment
     END DO
  ELSE
     DO k = 2,NPAR2_ARTH
        arth_r(k) = arth_r(k-1)+increment
     END DO
     temp = increment*REAL(NPAR2_ARTH)
     k    = NPAR2_ARTH
     DO
        IF (k >= n) EXIT
        k2                    = k+k
        arth_r(k+1:min(k2,n)) = temp+arth_r(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
     END DO
  END IF
END FUNCTION arth_r

FUNCTION arth_d(first,increment,n)
  REAL (KIND=dp), INTENT(IN)   :: first, increment
  INTEGER, INTENT(IN)          :: n
  REAL (KIND=dp), DIMENSION(n) :: arth_d
  INTEGER                      :: k, k2
  REAL (KIND=dp)               :: temp
  IF (n > 0) arth_d(1) = first
  IF (n <= NPAR_ARTH) THEN
     DO k = 2,n
        arth_d(k) = arth_d(k-1)+increment
     END DO
  ELSE
     DO k = 2,NPAR2_ARTH
        arth_d(k) = arth_d(k-1)+increment
     END DO
     temp = increment*DBLE(NPAR2_ARTH)
     k    = NPAR2_ARTH
     DO
        IF (k >= n) EXIT
        k2                    = k+k
        arth_d(k+1:min(k2,n)) = temp+arth_d(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
     END DO
  END IF
END FUNCTION arth_d

FUNCTION arth_i(first,increment,n)
  INTEGER, INTENT(IN)   :: first, increment, n
  INTEGER, DIMENSION(n) :: arth_i
  INTEGER               :: k, k2, temp
  IF (n > 0) arth_i(1) = first
  IF (n <= NPAR_ARTH) THEN
     DO k = 2,n
        arth_i(k) = arth_i(k-1)+increment
     END DO
  ELSE
     DO k = 2,NPAR2_ARTH
        arth_i(k) = arth_i(k-1)+increment
     END DO
     temp = increment*NPAR2_ARTH
     k    = NPAR2_ARTH
     DO
        IF (k >= n) EXIT
        k2                    = k+k
        arth_i(k+1:min(k2,n)) = temp+arth_i(1:min(k,n-k))
        temp                  = temp+temp
        k                     = k2
     END DO
  END IF
END FUNCTION arth_i

END MODULE UTILITY
