
MODULE IN_OUT

  USE CONSTANTS
!  USE EXTENSION
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: READ_DATA, WRITE_DATA, EXTINCTION

CONTAINS

SUBROUTINE READ_DATA(fgrain)
! reads the .DAT files

  ! global variables modules
  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT
  USE MGET_TDIST

  IMPLICIT NONE

  ! local variables
  CHARACTER (LEN=max_len)     :: filename_tmp
  CHARACTER (LEN=max_len)     :: c_tmp
  CHARACTER (LEN=max_len)     :: the_char
  CHARACTER (LEN=10)          :: sub_char(20)

  ! data files
  TYPE (FICH)                 :: fgrain     != FICH ('GRAIN.DAT', 11)
  TYPE (FICH)                 :: falign     = FICH ('ALIGN.DAT', 21)
  TYPE (FICH)                 :: fmix       = FICH ('MIX_', 12)
  TYPE (FICH)                 :: fpol       = FICH ('POL_', 13)
  TYPE (FICH)                 :: fsize      = FICH ('SIZE_', 14)
  TYPE (FICH)                 :: fchrg      = FICH ('CHRG_', 14)
  TYPE (FICH)                 :: fspin      = FICH ('SPIN_', 15)
  TYPE (FICH)                 :: fisrf      = FICH ('ISRF.DAT', 16)
  TYPE (FICH)                 :: fgas       = FICH ('GAS.DAT', 17)
  TYPE (FICH)                 :: fcalor     = FICH ('C_', 19)
  TYPE (FICH)                 :: fbeta      = FICH ('BETA_', 20)
  TYPE (FICH)                 :: fdtls      = FICH ('DTLS_', 20)

  ! wavelength dependent data files
  TYPE (FICH)                 :: flamb_qabs = FICH ('LAMBDA.DAT', 50)
  TYPE (FICH)                 :: fqext      = FICH ('Q_', 51)
  TYPE (FICH)                 :: fgfac      = FICH ('G_', 52)
  TYPE (FICH)                 :: fq1ext     = FICH ('Q1_', 53)
  TYPE (FICH)                 :: fq2ext     = FICH ('Q2_', 54)
  TYPE (FICH)                 :: fqcirc     = FICH ('Qc_', 55)
  TYPE (FICH)                 :: fqHabs     = FICH ('QH', 56)

  INTEGER                     :: nmat, eof        ! nr of grain material in a composite grain
  INTEGER                     :: i, j, k, u       ! loop indices
  INTEGER                     :: nbg              ! nr of size values for bandgap Ebg
  INTEGER                     :: n_gas            ! check nr to read GAS.DAT only first time
  INTEGER                     :: n1_sopt          ! nr of param other than amin, amax for PLAW (1) and LOGN (2)
  INTEGER                     :: n2_sopt, n_sopt  ! nr of param (other than PLAW and LOGN) for size distribution, total
  REAL (KIND=dp)              :: au, zeta, zxp, gama
  REAL (KIND=dp)              :: da, aux, argu
  REAL (KIND=dp)              :: xx(nsize_max), yy(nsize_max)
  INTEGER                     :: ns_i, nsub
  REAL (KIND=dp), ALLOCATABLE :: lambisrf(:)      ! wavelengths in microns for radiation field
  REAL (KIND=dp), ALLOCATABLE :: isrf(:)          ! radiation field
  INTEGER                     :: flag_cut
  REAL (KIND=dp), ALLOCATABLE :: tmp1(:),tmp2(:), tmp3(:), tmp4(:), tmp5(:)
 
  !=====================================================
  ! init run keywords 
  n_ftemp = 0
  n_pdr_o = 0
  n_quiet = 0
  n_res_a = 0
  
  ! init type counters
  n_chrg = 0 
  n_zm   = 0 
  n_spin = 0
  n_beta = 0
  n_pol  = 0
  n_dtls = 0

  ! reading GRAIN.DAT, set up size dist and run parameters
  !-----------------------------------------------------------------
!  CALL GETENVIRONMENT("DUSTEM_DATA_PATH",data_path) ! for paths set up in your bash
  filename_tmp = TRIMCAT (data_path, dir_dat)
  filename_tmp = TRIMCAT (filename_tmp, fgrain%nom)
  OPEN (UNIT = fgrain%unit, FILE = filename_tmp, STATUS = 'old')
  the_char = READ_COM (fgrain%unit)
  the_char = UPCASE ( TRIM (the_char) )

  ! get global keywords (type insensitive)
  n_res_a = INDEX (the_char, 'RES_A')
  n_ftemp = INDEX (the_char, 'TEMP')
  n_quiet = INDEX (the_char, 'QUIET')
  n_pdr_o = INDEX (the_char, 'PDR')
  n_sdist = INDEX (the_char, 'SDIST')
  n_zdist = INDEX (the_char, 'ZDIST')

  IF (n_quiet == 0) THEN
     WRITE (*,*) '==========================================='
     WRITE (*,*) '             Running DUSTEM'
     WRITE (*,*) ''
     WRITE (*, FMT='(A18)')' >> read GRAIN.DAT'
  ENDIF

  ! get G0 factor for radiation field
  READ (UNIT=fgrain%unit, FMT=*) g0
  IF (n_quiet == 0) THEN
     WRITE (*, FMT='(A12,1X,A20,1x,A3,1x,1PE9.2)') 'Run Keywords', the_char, 'G0=', g0
     WRITE (*, FMT='(45x,A93)') 'nsize        t_opt             M_dust/M_H  rho(g/cm3)      a-range(cm)      alpha-a0    sigma'
  ENDIF

  ! get number of grain types
  u = 0
  DO
     READ (UNIT=fgrain%unit, FMT=*, END=999) the_char
     u = u + 1
  ENDDO
999 ntype = u
  REWIND(fgrain%unit)
  the_char = READ_COM (fgrain%unit)
  READ (UNIT=fgrain%unit, FMT=*) the_char

  ALLOCATE (gtype(ntype))
  ALLOCATE (mprop(ntype))
  ALLOCATE (rhom(ntype))
  ALLOCATE (as1(ntype))
  ALLOCATE (as2(ntype))
  ALLOCATE (nsize(ntype))
  ALLOCATE (t_opt(ntype))
  ALLOCATE (p_opt(ntype))
  ALLOCATE (size_ava(ntype,nsize_max))
  ALLOCATE (si_ava_l(ntype,nsize_max))
  ALLOCATE (sect_eff(ntype,nsize_max))
  ALLOCATE (ava(ntype,nsize_max))
  ALLOCATE (dloga(ntype,nsize_max))
  ALLOCATE (rho(ntype,nsize_max))
  ALLOCATE (f_mix(ntype,nsize_max))
  ALLOCATE (beta0(ntype),abeta(ntype),gbeta(ntype),bmax(ntype))
  ALLOCATE (nbeta(ntype),tbeta(ntype,nbeta_max),betav(ntype,nbeta_max))
  ALLOCATE (ltresh(ntype),lstiff(ntype))
  ALLOCATE (ldtresh(ntype),qdtls(ntype,nsize_max))
  ALLOCATE (f_pol(ntype,nsize_max))
  ALLOCATE (athresh(ntype),pstiff(ntype),plev(ntype))
  ALLOCATE (a_dtls(ntype),lc(ntype),c_delta(ntype))
  ALLOCATE (vt(ntype),Pmu(ntype),gamma_e(ntype),omega_m(ntype))
  ALLOCATE (tau_0(ntype),V0(ntype),Vmin(ntype),Vm(ntype))

  n_gas = 0
  DO i=1,ntype  !!! type loop
     READ (UNIT=fgrain%unit,FMT=*) gtype(i), nsize(i), c_tmp, mprop(i)
     t_opt(i) = UPCASE(TRIM(c_tmp))
     n2_sopt = 0
     BACKSPACE(fgrain%unit)
     READ (UNIT=fgrain%unit, FMT='(A)') the_char
     CALL PARSE(the_char,' ', sub_char, nsub) ! get nsub number of size distribution parameters
     
     ! get the size distribution from GRAIN.DAT parameters
     IF (INDEX(t_opt(i),'SIZE') == 0) THEN
        n_fsize = INDEX(t_opt(i),'PLAW') + INDEX(t_opt(i),'LOGN')
        IF (n_fsize == 0) THEN
           WRITE(*,*)'(F) DM_inout/READ_DATA: ',TRIM(gtype(i)),' undefined size distribution '
           STOP
        ENDIF
        IF (INDEX(t_opt(i),'LOGN') > 0) n1_sopt = 2
        IF (INDEX(t_opt(i),'PLAW') > 0) THEN
           n1_sopt = 1
           IF ( ((INDEX(t_opt(i),'-ED') > 0 ) .AND. (INDEX(t_opt(i),'-CV') > 0 ) .AND. (nsub /= 14)) .AND. &
                ((INDEX(t_opt(i),'-ED') > 0 ) .OR. (INDEX(t_opt(i),'-CV') > 0 ) .AND. (nsub /= 11)) ) THEN
              WRITE(*,*)'(F) DM_inout/READ_DATA: ',TRIM(gtype(i)),' missing parameters for PLAW-...'
              STOP
           ENDIF
           IF (INDEX(t_opt(i),'-ED') > 0) n2_sopt = 3
           IF (INDEX(t_opt(i),'-CV') > 0) n2_sopt = n2_sopt + 3
        ENDIF
        n_sopt = n1_sopt + n2_sopt
        ALLOCATE ( sd_par(n_sopt) )
        
        BACKSPACE(fgrain%unit)
        READ (UNIT=fgrain%unit,FMT=*) gtype(i), nsize(i), c_tmp, mprop(i), rhom(i),as1(i), as2(i), &
             & (sd_par(u), u=1,n_sopt)
        IF (n_quiet == 0) THEN
           WRITE (*, FMT='(A40,5x,I3,7x,A20,2x,1P,E9.2,3x,E9.2,2x,2(E9.2,1x),8(E9.2,1x))') &
                & gtype(i), nsize(i), t_opt(i), mprop(i), rhom(i), as1(i), as2(i), sd_par(:)
        ENDIF

        IF (INDEX(t_opt(i),'PLAW') > 0 ) THEN
           IF (nsize(i) /= 1) THEN
              da = ( LOG(as2(i)) - LOG(as1(i)) ) / DBLE(nsize(i)-1)
           ELSE
              da = LOG(as2(i)) - LOG(as1(i))
           ENDIF
           DO j=1,nsize(i)
              aux = LOG(as1(i)) + DBLE(j-1) * da
              argu = (4.0_dp + sd_par(1)) * aux
              IF (argu > -350.0_dp) THEN
                 ava(i,j) = EXP(argu)
              ELSE
                 ava(i,j) = 0.0_dp
              ENDIF
              si_ava_l(i,j) = aux
              size_ava(i,j) = EXP(aux)
              sect_eff(i,j) = xpi * size_ava(i,j)**2
              dloga(i,j) = da
              rho(i,j) = rhom(i)
           ENDDO
           IF (size_ava(i,nsize(i)) /= as2(i)) THEN ! correct rounding
              aux = LOG(as2(i))
              argu = (4.0_dp + sd_par(1)) * aux
              IF (argu > -350.0_dp) THEN
                 ava(i,nsize(i)) = EXP(argu)
              ELSE
                 ava(i,nsize(i)) = 0.0_dp
              ENDIF
              si_ava_l(i,nsize(i)) = aux
              size_ava(i,nsize(i)) = EXP(aux)
              sect_eff(i,nsize(i)) = xpi * size_ava(i,nsize(i))**2
           ENDIF
           
           ! apply exponential decay 
           IF (INDEX(t_opt(i),'-ED') > 0) THEN 
              WHERE (size_ava(i,:) >= sd_par(n1_sopt+1))
                 ava(i,1:nsize(i)) = ava(i,1:nsize(i)) * &
                      & EXP(-( (size_ava(i,1:nsize(i))-sd_par(n1_sopt+1))/sd_par(n1_sopt+2) )**sd_par(n1_sopt+3))
              ENDWHERE
           ENDIF
           
           ! apply curvature
           IF (INDEX(t_opt(i),'-CV') > 0) THEN
              IF (INDEX(t_opt(i),'-ED') > 0) THEN ! get the CV parameters
                 au = sd_par(n1_sopt+4)
                 zeta = ABS(sd_par(n1_sopt+5))
                 zxp = SIGN(1.0_dp,sd_par(n1_sopt+5))
                 gama = sd_par(n1_sopt+6)
              ELSE
                 au = sd_par(n1_sopt+1)
                 zeta = ABS(sd_par(n1_sopt+2))
                 zxp = SIGN(1.0_dp,sd_par(n1_sopt+2))
                 gama = sd_par(n1_sopt+3)
              ENDIF
              ava(i,1:nsize(i)) = ava(i,1:nsize(i)) * ( 1.0_dp+zeta*(size_ava(i,1:nsize(i))/au)**gama )**zxp
           ENDIF
           
        ELSE IF (INDEX(t_opt(i),'LOGN') > 0 ) THEN
           IF ((sd_par(1) == 0.0_dp) .OR. (sd_par(2) == 0.0_dp)) THEN
              WRITE(*,*)'(F) DM_inout/READ_DATA:',TRIM(gtype(i)),' centroid or sigma of log-normal cannot be 0'
              STOP
           ENDIF
           IF (nsize(i) /= 1) THEN
              da = ( LOG(as2(i)) - LOG(as1(i)) ) / DBLE(nsize(i)-1)
           ELSE
              da = LOG(as2(i)) - LOG(as1(i))
           ENDIF
           DO j=1,nsize(i)
              aux = LOG(as1(i)) + DBLE(j-1) * da
              argu = 3.0_dp*aux - 0.5_dp * ( (aux - LOG(sd_par(1))) / sd_par(2) )**2
              IF (argu > -350.0_dp) THEN
                 ava(i,j) = EXP(argu)
              ELSE
                 ava(i,j) = 0.0_dp
              ENDIF
              si_ava_l(i,j) = aux
              size_ava(i,j) = EXP(aux)
              sect_eff(i,j) = xpi * size_ava(i,j)**2
              dloga(i,j) = da
              rho(i,j) = rhom(i)
           ENDDO
           IF (size_ava(i,nsize(i)) /= as2(i)) THEN ! correct rounding
              aux = LOG(as2(i))
              argu = 3.0_dp*aux - 0.5_dp * ( (aux - LOG(sd_par(1))) / sd_par(2) )**2
              IF (argu > -350.0_dp) THEN
                 ava(i,nsize(i)) = EXP(argu)
              ELSE
                 ava(i,nsize(i)) = 0.0_dp
              ENDIF
              si_ava_l(i,nsize(i)) = aux
              size_ava(i,nsize(i)) = EXP(aux)
              sect_eff(i,nsize(i)) = xpi * size_ava(i,nsize(i))**2
           ENDIF
        ENDIF
        DEALLOCATE (sd_par)

     ELSE IF (INDEX(t_opt(i),'SIZE') > 0) THEN  ! use size dist. from SIZE_*.DAT file

       ! open SIZE_TYPE.DAT
        filename_tmp = TRIMCAT (data_path, dir_dat)
        filename_tmp = TRIMCAT (filename_tmp, fsize%nom)
        filename_tmp = TRIMCAT (filename_tmp, gtype(i))
        filename_tmp = TRIMCAT (filename_tmp, '.DAT')
        c_tmp = TRIMCAT (fsize%nom, gtype(i))
        c_tmp = TRIMCAT (c_tmp, '.DAT')
        OPEN (UNIT = fsize%unit, FILE = filename_tmp, STATUS = 'old')

        ! read doc
        the_char = READ_COM (fsize%unit)
        READ (the_char, FMT=*) nmat
        IF (.NOT. ALLOCATED(sbulk)) ALLOCATE (sbulk(nmat))
        READ (UNIT=fsize%unit, FMT=*) (sbulk(u), u=1,nmat)
        READ (UNIT=fsize%unit, FMT=*) the_char
        READ (UNIT=fsize%unit, FMT=*) the_char

        ! read data
        READ (UNIT=fsize%unit, FMT=*) nsize(i)
        IF (n_quiet == 0) THEN
           WRITE (*, FMT='(A7,2x,A30,A40)') '>>read ',TRIM(c_tmp), &
                & '    a-range (cm)     nsize'
        ENDIF
        DO j=1,nsize(i)
           READ (UNIT=fsize%unit, FMT=*) size_ava(i,j), dloga(i,j), ava(i,j), rho(i,j)
        ENDDO
        si_ava_l(i,1:nsize(i)) = LOG(size_ava(i,1:nsize(i)))
        sect_eff(i,1:nsize(i)) = xpi * size_ava(i,1:nsize(i))**2
        IF (n_quiet == 0) THEN
           c_tmp = 'bulk'
           DO u=1,nmat
              c_tmp = TRIMCAT(c_tmp, '-'//sbulk(u))
           ENDDO
           !WRITE (*, FMT='(A39,1P,2(E9.2,1x),1x,I3,8x,i1,2x,1PE9.2)') c_tmp, size_ava(i,1), size_ava(i,nsize(i)), nsize(i)
           WRITE (*, FMT='(A39,1P,2(E9.2,1x),1x,I3)') c_tmp, size_ava(i,1), size_ava(i,nsize(i)), nsize(i)
        ENDIF
        CLOSE (UNIT=fsize%unit)
        
     ENDIF ! size distribution

     ! treat the other type options
     IF (INDEX(t_opt(i),'MIX') > 0) THEN
        filename_tmp = TRIMCAT (data_path, dir_dat)
        filename_tmp = TRIMCAT (filename_tmp, fmix%nom)
        filename_tmp = TRIMCAT (filename_tmp, gtype(i))
        filename_tmp = TRIMCAT (filename_tmp, '.DAT')
        the_char = TRIMCAT (fmix%nom, gtype(i))
        the_char = TRIMCAT (the_char, '.DAT')
        c_tmp = the_char
        IF (n_quiet == 0) THEN
           WRITE (*, FMT='(A7,2X,A40)') '>>read ', TRIM(the_char)
        ENDIF
        OPEN (UNIT=fmix%unit, FILE=filename_tmp, STATUS='old')
        the_char = READ_COM (fmix%unit)
        READ (the_char, FMT=*, iostat=eof) f_mix(i,1)
        j = 1
        DO WHILE (eof == 0)
           j = j + 1
           READ (UNIT=fmix%unit, FMT=*, iostat=eof) f_mix(i,j)
        ENDDO
        IF (j-1 /= nsize(i)) THEN
           WRITE (*,*) '(F) DM_INOUT/MIX:',j-1,' points in ',TRIM(c_tmp),' not equal to ', &
                & nsize(i),' nr of sizes'
           STOP
        ENDIF
        CLOSE (UNIT=fmix%unit)
     ELSE IF (INDEX(t_opt(i),'MIX') == 0) THEN
        f_mix(i,:) = 1.0_dp
     ENDIF

     ! charge distribution and spinning dust emission
     IF ((INDEX(t_opt(i),'CHRG') > 0) .OR. (INDEX(t_opt(i),'SPIN') > 0)) THEN
        IF (n_gas == 0) THEN
           ! first time get the gas parameters
           filename_tmp = TRIMCAT (data_path, dir_dat)
           filename_tmp = TRIMCAT (filename_tmp, fgas%nom)
           the_char = TRIMCAT (fgas%nom, '.DAT')
           IF (n_quiet == 0) THEN
              WRITE (*, FMT='(A7,2X,A40)') '>>read ', TRIM(the_char)
           ENDIF
           OPEN (UNIT=fgas%unit, FILE=filename_tmp, STATUS='old')
           the_char = READ_COM (fgas%unit)
           READ (the_char, FMT=*) t_gas,hden,h2den,cr_rate,aux
           IF (aux > 0.0_dp) THEN 
              g0 = aux ! superseding GRAIN.DAT value
              IF (n_quiet == 0) WRITE (*, FMT='(A23,1X,1PE9.2)') '        superseding G0=', g0
           ENDIF
           READ(UNIT=fgas%unit, FMT=*) nion
           ALLOCATE (iden(nion), mi(nion), zi(nion), pol_i(nion))
           DO j=1,nion
              READ(UNIT=fgas%unit, FMT=*) iden(j), mi(j), zi(j), pol_i(j)
           ENDDO
           eden = iden(1) ! electrons on 1st line
           CLOSE (UNIT=fgas%unit)
           n_gas = n_gas + 1
        ENDIF

        ! get the charge parameters
        if (n_chrg == 0) THEN 
           ALLOCATE (wf(ntype), ebg(ntype,nsize_max), p_ea(ntype,2), s_ea(ntype,2))
           ALLOCATE (le(ntype), p_y(ntype,3))
           ALLOCATE (p_uait(ntype,3))
           ALLOCATE (m1mass(ntype))
        ENDIF
        filename_tmp = TRIMCAT (data_path, dir_dat)
        filename_tmp = TRIMCAT (filename_tmp, fchrg%nom)
        filename_tmp = TRIMCAT (filename_tmp, gtype(i))
        filename_tmp = TRIMCAT (filename_tmp, '.DAT')
        the_char = TRIMCAT (fchrg%nom, gtype(i))
        the_char = TRIMCAT (the_char, '.DAT')
        IF (n_quiet == 0) THEN
           WRITE (*, FMT='(A7,2X,A40)') '>>read ', TRIM(the_char)
        ENDIF
        OPEN (UNIT=fchrg%unit, FILE=filename_tmp, STATUS='old')
        the_char = READ_COM (fchrg%unit)
        READ (the_char, FMT=*) wf(i), p_ea(i,1), p_ea(i,2)
        READ (UNIT=fchrg%unit, FMT=*) s_ea(i,1), s_ea(i,2)
        READ (UNIT=fchrg%unit, FMT=*) le(i), p_y(i,1), p_y(i,2), p_y(i,3)
        READ (UNIT=fchrg%unit, FMT=*) p_uait(i,1), p_uait(i,2), p_uait(i,3)
        READ (UNIT=fchrg%unit, FMT=*) m1mass(i)
        READ (UNIT=fchrg%unit, FMT=*) nbg
        DO u=1,nbg 
           READ (UNIT=fchrg%unit, FMT=*) ebg(i,u)
        ENDDO
        IF ((nbg /= nsize(i)) .AND. (nbg /= 1)) THEN
           WRITE (*,*) '(F) DM_INOUT/CHRG:',nbg,' bandgap points, not equal to ', &
                & nsize(i),' nr of sizes'
           STOP
        ELSE IF (nbg==1) THEN 
           ebg(i,:) = ebg(i,1) 
        ENDIF
        CLOSE (UNIT=fchrg%unit)
        n_chrg = n_chrg + 1

        ! get the spin parameters
        IF (INDEX(t_opt(i),'SPIN') > 0) THEN 
           if (n_spin == 0) ALLOCATE (m0(ntype))
           filename_tmp = TRIMCAT (data_path, dir_dat)
           filename_tmp = TRIMCAT (filename_tmp, fspin%nom)
           filename_tmp = TRIMCAT (filename_tmp, gtype(i))
           filename_tmp = TRIMCAT (filename_tmp, '.DAT')
           the_char = TRIMCAT (fspin%nom, gtype(i))
           the_char = TRIMCAT (the_char, '.DAT')
           IF (n_quiet == 0) THEN
              WRITE (*, FMT='(A7,2X,A40)') '>>read ', TRIM(the_char)
           ENDIF
           OPEN (UNIT=fspin%unit, FILE=filename_tmp, STATUS='old')
           the_char = READ_COM (fspin%unit)
           READ (the_char, FMT=*) m0(i)
           CLOSE (UNIT=fspin%unit)
           n_spin = n_spin + 1
        ENDIF
     ENDIF
 
     ! beta(T)-correction
     IF (INDEX(t_opt(i),'BETA') > 0) THEN
        n_beta = n_beta + 1
        filename_tmp = TRIMCAT (data_path, dir_dat)
        filename_tmp = TRIMCAT (filename_tmp, fbeta%nom)
        filename_tmp = TRIMCAT (filename_tmp, gtype(i))
        filename_tmp = TRIMCAT (filename_tmp, '.DAT')
        the_char = TRIMCAT (fbeta%nom, gtype(i))
        the_char = TRIMCAT (the_char, '.DAT')
        IF (n_quiet == 0) THEN
           WRITE (*,FMT='(A7,2X,A40)')'>>read ', TRIM(the_char)
        ENDIF
        OPEN (UNIT=fbeta%unit, FILE=filename_tmp, STATUS='old')
        the_char = READ_COM (fbeta%unit)
        READ (the_char,FMT=*) beta0(i), abeta(i), gbeta(i), bmax(i)
        READ (UNIT=fbeta%unit, FMT=*) ltresh(i), lstiff(i)
        READ (UNIT=fbeta%unit, FMT=*) nbeta(i)
        IF (nbeta(i) > 0) THEN
           DO k=1, nbeta(i) 
              READ (UNIT=fbeta%unit, FMT=*) tbeta(i,k), betav(i,k)
           ENDDO
        ENDIF
        CLOSE (UNIT=fbeta%unit)
        ltresh(i) = ltresh(i) * 1.00e-4_dp   ! microns to cm
     ENDIF

     ! DCD/TLS effects
     IF (INDEX(t_opt(i),'DTLS') > 0) THEN
        n_dtls = n_dtls + 1
        filename_tmp = TRIMCAT (data_path, dir_dat)
        filename_tmp = TRIMCAT (filename_tmp, fdtls%nom)
        filename_tmp = TRIMCAT (filename_tmp, gtype(i))
        filename_tmp = TRIMCAT (filename_tmp, '.DAT')
        the_char = TRIMCAT (fdtls%nom, gtype(i))
        the_char = TRIMCAT (the_char, '.DAT')
        IF (n_quiet == 0) THEN
           WRITE (*,FMT='(A7,2X,A40)')'>>read ', TRIM(the_char)
        ENDIF
        OPEN (UNIT=fdtls%unit, FILE=filename_tmp, STATUS='old')
        the_char = READ_COM (fdtls%unit)
        READ (the_char,FMT=*) a_dtls(i), lc(i), c_delta(i)
        READ (UNIT=fdtls%unit, FMT=*) vt(i), Pmu(i), gamma_e(i) 
        READ (UNIT=fdtls%unit, FMT=*) omega_m(i), tau_0(i), V0(i), Vmin(i), Vm(i)
        READ (UNIT=fdtls%unit, FMT=*) ldtresh(i)
        CLOSE (UNIT=fdtls%unit)
     ENDIF

     ! nr of grain types for polarisation
     IF (INDEX(t_opt(i),'POL') > 0) n_pol = n_pol + 1
     
  ENDDO   !!! type loop
  CLOSE (UNIT = fgrain%UNIT)

  ! polarisation
  IF (n_pol > 0) THEN

     ! Reading ALIGN file
     filename_tmp = TRIMCAT (data_path, dir_dat)
     filename_tmp = TRIMCAT (filename_tmp, falign%nom)
     OPEN (UNIT = falign%unit, FILE = filename_tmp, STATUS = 'old')
     the_char = READ_COM (falign%unit)
     the_char = UPCASE ( TRIM (the_char) )

     ! get global keywords (type insensitive)
     ! same alignment law for all types if univ set
     n_univ = INDEX(the_char,'UNIV')
     ! polarization model
     n_lin = INDEX(the_char,'LIN')
     n_circ = INDEX (the_char, 'CIRC')
     n_anis = INDEX (the_char, 'ANIS')
     IF (n_lin + n_circ + n_anis .eq. 0) THEN
        WRITE(*,*)'(F) DM_inout/READ_DATA ALIGN.DAT: lacking keywords LIN or CIRC or ANIS'
        STOP
     ENDIF     

     IF (n_quiet == 0) THEN
        WRITE (*,*)
        WRITE (*, FMT='(A17)')'>> read ALIGN.DAT'
     ENDIF

     ! get anisotropy factor for the radiation field: anisG0
     READ (UNIT=falign%unit, FMT=*) anisg0

     ! get alignment law
     i = 1
     READ (UNIT=falign%unit,FMT=*) p_opt(i), c_tmp
     p_opt(i) = UPCASE ( TRIM (p_opt(i)) )
     n_falig = INDEX(p_opt(i),'IDG') + INDEX(p_opt(i),'RAT') + INDEX(p_opt(i),'PAR')
     IF (n_falig == 0) THEN
        WRITE(*,*)'(F) DM_inout/READ_DATA: in ALIGN.DAT ',TRIM(p_opt(i)),' undefined alignment function'
        STOP
     ENDIF
     READ (UNIT=falign%unit,FMT=*,iostat=eof) c_tmp ! checks end of file
     BACKSPACE(falign%unit)
     BACKSPACE(falign%unit)

    IF (n_univ .NE. 0) THEN ! universal alignment
        READ (UNIT=falign%unit,FMT=*) p_opt(i), athresh(i), pstiff(i), plev(i)
        p_opt(:) = UPCASE(TRIM(p_opt(i)))
        athresh(:) = athresh(i)
        pstiff(:) = pstiff(i)
        plev(:) = plev(i)
     ELSE                   ! type sensitive alignment
        IF ((eof < 0) .AND. (n_pol >1)) then
           WRITE(*,*)'(F) DM_inout/READ_DATA: in ALIGN.DAT while UNIV not set not enough lines'
           STOP 
        ENDIF
        DO i=1,ntype
           IF (INDEX(t_opt(i),'POL') > 0) then
              READ (UNIT=falign%unit,FMT=*) p_opt(i), athresh(i), pstiff(i), plev(i)
              p_opt(i) = UPCASE(TRIM(p_opt(i)))
           ENDIF
        ENDDO
     ENDIF

     DO i=1,ntype
        IF (INDEX(t_opt(i),'POL') > 0) THEN
           IF (INDEX(p_opt(i),'PAR') > 0) THEN           ! Parametric model
              f_pol(i,1:nsize(i))= 0.5_dp*plev(i) * (1.0_dp + TANH(LOG(size_ava(i,1:nsize(i))*1d4/athresh(i)) / pstiff(i) ) )
           endif

           IF (INDEX(p_opt(i),'IDG') > 0) THEN           ! Imperfect Davis-Greenstein
           endif

           IF (INDEX(p_opt(i),'RAT') > 0) THEN           ! RATs
           endif
        ENDIF
     ENDDO
  ENDIF ! polarisation
  
  ! reading heat capacities (C_*.DAT)
  !-----------------------------------------------------------------
  ALLOCATE (nsize_type(ntype),nsz1(ntype))
  ALLOCATE (size_type(ntype,nsize_max_qabs))
  ALLOCATE (calor(ntype,nsize_max_qabs,ntempmax))
  ALLOCATE (temp_cal(ntype,ntempmax))
  ALLOCATE (n_temp(ntype))

  ! initialize
  temp_cal(:,:) = 0.0_dp
  calor(:,:,:) = 0.0_dp

  IF (n_quiet == 0) THEN
     WRITE (*,*)
     WRITE (*, FMT='(A17,32x,A51)') '>>read C_TYPE.DAT','a-range (cm)     nsize       T-range (K)      NTEMP'
  ENDIF

  DO i=1,ntype

     ! open file
     filename_tmp = TRIMCAT (data_path, dir_capa)
     filename_tmp = TRIMCAT (filename_tmp, fcalor%nom)
     filename_tmp = TRIMCAT (filename_tmp, gtype(i))
     filename_tmp = TRIMCAT (filename_tmp, '.DAT')
     OPEN ( UNIT=fcalor%unit,FILE =filename_tmp,STATUS='old' )

     ! read doc
     the_char = READ_COM(fcalor%unit)
     READ (the_char, FMT=*) nsize_type(i)

     ! checking C file for nr of sizes
     IF (nsize_type(i) > nsize_max_qabs) THEN
        WRITE (*,*) ''
        WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', filename_tmp
        WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
        WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
        WRITE (*,*) ''
        STOP
     ENDIF

     ! get sizes
     READ (UNIT = fcalor%unit,FMT = *) (size_type(i,u), u = 1,nsize_type(i))

     ! read data
     READ (UNIT = fcalor%unit,FMT = *) n_temp(i)
     DO k=1,n_temp(i)
        READ (unit = fcalor%unit, FMT = *) temp_cal(i,k), (calor(i,j,k), j=1,nsize_type(i))
     ENDDO
     IF (n_quiet == 0) THEN
        WRITE (*, FMT='((A40,5x,1P,2( 2(E9.2,1x),1x,I3,5X)))') &
             gtype(i), 1.0e-4_dp*size_type(i,1), 1.0e-4_dp*size_type(i,nsize_type(i)), nsize_type(i), &
             10.0_dp**temp_cal(i,1), 10.0_dp**temp_cal(i,n_temp(i)), n_temp(i)
     ENDIF
     CLOSE (UNIT=fcalor%UNIT)
  ENDDO

  ! reading lambda grid for all Q files (LAMBDA.DAT)
  !---------------------------------------------------------------------------------
  filename_tmp = TRIMCAT (data_path, dir_qabs)
  filename_tmp = TRIMCAT (filename_tmp, flamb_qabs%nom)
  OPEN (UNIT=flamb_qabs%unit, FILE=filename_tmp, STATUS='old')
  the_char = READ_COM (flamb_qabs%unit)
  READ (the_char, FMT=*) n_qabs
  ALLOCATE (lamb_qabs(n_qabs))
  ALLOCATE (freq_qabs(n_qabs))
  ALLOCATE (lfrq_qabs(n_qabs))
  ALLOCATE (tmp1(n_qabs))
  ALLOCATE (f_beta(ntype,n_qabs))
  DO k=1,n_qabs
    READ (UNIT=flamb_qabs%unit, FMT=*) lamb_qabs(k)
  ENDDO
  IF (n_quiet == 0) THEN
     WRITE (*,*)
     WRITE (*, FMT='(A51)') '>>read LAMBDA.DAT          w-range (microns)  nwave'
     WRITE (*, FMT='(A25,1P,2(E9.2,1x),1x,I4)') '', lamb_qabs(1), lamb_qabs(n_qabs), n_qabs
  ENDIF
  CLOSE (UNIT=flamb_qabs%unit)
  lamb_qabs = lamb_qabs * 1.0e-4_dp              ! Convert from micron to cm
  freq_qabs = clight / lamb_qabs
  ! We sort frequencies in reverse order
  tmp1 = freq_qabs
  DO k=1,n_qabs
     freq_qabs(k) = tmp1(n_qabs-k+1)
  ENDDO
  lfrq_qabs = LOG(freq_qabs)


  ! reading Qabs, Qsca (Q_*.DAT) and G_*.DAT 
  !---------------------------------------------------------------
  ALLOCATE (q_abs(ntype,  nsize_max_qabs, n_qabs))
  ALLOCATE (qdiff(ntype,  nsize_max_qabs, n_qabs))
  ALLOCATE (gfac(ntype,   nsize_max_qabs, n_qabs))

  ALLOCATE (qi_abs(ntype, nsize_max, n_qabs))
  ALLOCATE (qidiff(ntype, nsize_max, n_qabs))
  ALLOCATE (gifac(ntype,  nsize_max, n_qabs))
  ALLOCATE (qicirc(ntype, nsize_max, n_qabs))
  ALLOCATE (qiH_abs(ntype,nsize_max, n_qabs))

  ALLOCATE (tmp2(n_qabs))
  ALLOCATE (tmp3(n_qabs))
  ALLOCATE (masstot(ntype))

  ! initialise q_abs, g & size_type
  q_abs(:,:,:)   = 0.0_dp
  qdiff(:,:,:)   = 0.0_dp
  gfac(:,:,:)    = 0.0_dp
  qi_abs(:,:,:)  = 0.0_dp
  qidiff(:,:,:)  = 0.0_dp
  qicirc(:,:,:)  = 0.0_dp
  gifac(:,:,:)   = 0.0_dp
  size_type(:,:) = 0.0_dp

  IF (n_quiet == 0) THEN
     WRITE (*,*)
     WRITE (*, FMT='(A17,32x,A22)') '>>read Q_TYPE.DAT','a-range (cm)     nsize'          
     IF (n_pdr_o /= 0) THEN
        WRITE (*, FMT='(A17,32x,A22)') '>>read G_TYPE.DAT','a-range (cm)     nsize'
     ENDIF
  ENDIF

  DO i=1,ntype
    filename_tmp = TRIMCAT (data_path, dir_qabs)
    filename_tmp = TRIMCAT (filename_tmp, fqext%nom)
    filename_tmp = TRIMCAT (filename_tmp, gtype(i))
    filename_tmp = TRIMCAT (filename_tmp, '.DAT')
    OPEN (UNIT=fqext%unit, FILE=filename_tmp, STATUS='old')

    ! read doc
    the_char = READ_COM (fqext%unit)
    READ (the_char, FMT=*) nsz1(i)

    ! checking Q file
    IF (nsz1(i) > nsize_max_qabs) THEN
       WRITE (*,*) ''
       WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', TRIM(filename_tmp)
       WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
       WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
       WRITE (*,*) ''
       STOP
    ENDIF
    IF (nsz1(i) /= nsize_type(i)) THEN
       WRITE (*,*) ''
       WRITE (*,*) '(F) DM_inout/READ_DATA: odd number of sizes in ', TRIM(filename_tmp)
       WRITE (*,*) '                        found', nsz1(i), ' sizes but ', nsize_type(i), ' sizes expected'
       WRITE (*,*) '                        Q, C and G files must have same number of sizes'
       WRITE (*,*) ''
       STOP
    ENDIF
    nsize_type(i) = nsz1(i)

    ! get sizes
    READ (UNIT=fqext%unit, FMT=*) (size_type(i,u), u=1,nsize_type(i))
    IF (n_quiet == 0) THEN
       WRITE(*, FMT='(A40,5x,1P,2(E9.2,1x),1x,I3)') gtype(i), 1.0e-4_dp*size_type(i,1), &
             1.0e-4_dp*size_type(i,(nsize_type(i))), nsize_type(i)
    ENDIF

    ! get Qabs
    the_char = READ_COM (fqext%unit)
    BACKSPACE fqext%unit   ! READ_COM only reads 1st column
    DO k=1,n_qabs
       READ (UNIT=fqext%unit, FMT=*) (q_abs(i,u,k), u=1,nsize_type(i))
    ENDDO

    ! get Qsca
    the_char = READ_COM (fqext%unit)
    BACKSPACE fqext%unit   ! READ_COM only reads 1st column
    DO k=1,n_qabs
       READ (UNIT=fqext%unit, FMT=*) (qdiff(i,u,k), u=1,nsize_type(i))
    ENDDO

    CLOSE (UNIT=fqext%unit)

    IF (n_pdr_o /= 0) THEN
       filename_tmp = TRIMCAT (data_path, dir_qabs)
       filename_tmp = TRIMCAT (filename_tmp, fgfac%nom)
       filename_tmp = TRIMCAT (filename_tmp, gtype(i))
       filename_tmp = TRIMCAT (filename_tmp, '.DAT')
       OPEN (UNIT=fgfac%unit, FILE=filename_tmp, STATUS='old')

       ! read doc
       the_char = READ_COM (fgfac%unit)
       READ (the_char, FMT=*) nsz1(i)

       ! checking G file
       IF (nsz1(i) > nsize_max_qabs) THEN
          WRITE (*,*) ''
          WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', TRIM(filename_tmp)
          WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
          WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
          WRITE (*,*) ''
          STOP
       ENDIF
       IF (nsz1(i) /= nsize_type(i)) THEN
          WRITE (*,*) ''
          WRITE (*,*) '(F) DM_inout/READ_DATA: odd number of sizes in ', TRIM(filename_tmp)
          WRITE (*,*) '                        found', nsz1(i), ' sizes but ', nsize_type(i), ' sizes expected'
          WRITE (*,*) '                        Q, C and G files must have same number of sizes'
          WRITE (*,*) ''
          STOP
       ENDIF
       nsize_type(i) = nsz1(i)

       ! get sizes
       READ (UNIT=fgfac%unit, FMT=*) (size_type(i,u), u=1,nsize_type(i))
       IF (n_quiet == 0) THEN
          WRITE(*, FMT='(A40,5x,1P,2(E9.2,1x),1x,I3)') gtype(i), 1.0e-4_dp*size_type(i,1), &
               1.0e-4_dp*size_type(i,(nsize_type(i))), nsize_type(i)
       ENDIF

       ! get g factors
       the_char = READ_COM (fgfac%unit)
       BACKSPACE fgfac%unit   ! READ_COM only reads 1st column
       DO k=1,n_qabs
          READ (UNIT=fgfac%unit, FMT=*) (gfac(i,u,k), u=1,nsize_type(i))
       ENDDO

       CLOSE (UNIT=fgfac%unit)
    ENDIF

    ! size interpolation (once for all)
    ! mass and g-normalization for each grain type
    DO j=1,nsize(i)
       aux  = 1.e4_dp*size_ava(i,j)
       CALL GET_QEXT (aux, i, tmp1, tmp2, tmp3)
       qi_abs(i,j,:) = tmp1(:)
       qidiff(i,j,:) = tmp2(:)
       gifac(i,j,:)  = tmp3(:)
       ns_i = nsize(i)
       xx(1:ns_i) = si_ava_l(i,1:ns_i)
       yy(1:ns_i) = ava(i,1:ns_i) * rho(i,1:ns_i) * 4.0_dp * xpi / 3.0_dp
       masstot(i) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
    ENDDO

    ! get the beta threshold
    IF (INDEX(t_opt(i),'BETA') > 0) THEN
!     f_beta(i,:) = 0.5_dp * ( 1.0_dp + TANH( 4.0_dp*lstiff(i)*(lamb_qabs(:)/ltresh(i)-1.0_dp)) )
       f_beta(i,:) = 0.5_dp * ( 1.0_dp + TANH( 4.0_dp*(LOG10(lamb_qabs(:)) - LOG10(ltresh(i)))/lstiff(i) ) )
    ENDIF

    ! get the DTLS threshold and normalization
    IF (INDEX(t_opt(i),'DTLS') > 0) THEN
       DO j=1, nsize(i) 
          Qdtls(i,j)   = INTPOL ( qi_abs(i,j,:), freq_qabs(:), n_qabs, 1e4_dp*clight/ldtresh(i) )
       ENDDO          
    ENDIF

  ENDDO ! TYPE loop on i

  ! convert sizes microns --> cm
  size_type(:,:) = size_type(:,:) * 1.0e-4_dp


  ! reading Q1_*.DAT and Q2_*.DAT files
  !-----------------------------------------------------------------
  IF (n_pol > 0) THEN

     ALLOCATE (q1_abs(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (q2_abs(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (q1diff(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (q2diff(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (qcirc(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (qH1_abs(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (qH2_abs(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (qH_abs(ntype, nsize_max_qabs, n_qabs))
     ALLOCATE (q1i_abs(ntype, nsize_max, n_qabs))
     ALLOCATE (q2i_abs(ntype, nsize_max, n_qabs))
     ALLOCATE (q1idiff(ntype, nsize_max, n_qabs))
     ALLOCATE (q2idiff(ntype, nsize_max, n_qabs))
     ALLOCATE (tmp4(n_qabs),tmp5(n_qabs))

     q1_abs(:,:,:)   = 0.0_dp
     q2_abs(:,:,:)   = 0.0_dp
     q1diff(:,:,:)   = 0.0_dp
     q2diff(:,:,:)   = 0.0_dp
     q1i_abs(:,:,:)  = 0.0_dp
     q2i_abs(:,:,:)  = 0.0_dp
     q1idiff(:,:,:)  = 0.0_dp
     q2idiff(:,:,:)  = 0.0_dp

     DO i=1,ntype

        IF (INDEX(t_opt(i),'POL') > 0) THEN

           ! reading Q1_*.DAT file
           IF (n_quiet == 0) THEN
              WRITE (*,*)
              WRITE (*, FMT='(A18,31x,A22)') '>>read Q1_TYPE.DAT','a-range (cm)     nsize'
           ENDIF
           filename_tmp = TRIMCAT (data_path, dir_qabs)
           filename_tmp = TRIMCAT (filename_tmp, fq1ext%nom)
           filename_tmp = TRIMCAT (filename_tmp, gtype(i))
           filename_tmp = TRIMCAT (filename_tmp, '.DAT')
           OPEN (UNIT=fq1ext%unit, FILE=filename_tmp, STATUS='old')

           ! read doc
           the_char = READ_COM (fq1ext%unit)
           READ (the_char, FMT=*) nsz1(i)

           ! checking Q file
           IF (nsz1(i) > nsize_max_qabs) THEN
              WRITE (*,*) ''
              WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', TRIM(filename_tmp)
              WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
              WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
              WRITE (*,*) ''
              STOP
           ENDIF
           IF (nsz1(i) /= nsize_type(i)) THEN
              WRITE (*,*) ''
              WRITE (*,*) '(F) DM_inout/READ_DATA: odd number of sizes in ', TRIM(filename_tmp)
              WRITE (*,*) '                        found', nsz1(i), ' sizes but ', nsize_type(i), ' sizes expected'
              WRITE (*,*) '                        Q, C and G files must have same number of sizes'
              WRITE (*,*) ''
              STOP
           ENDIF
           nsize_type(i) = nsz1(i)

           ! get sizes
           READ (UNIT=fq1ext%unit, FMT=*) (size_type(i,u), u=1,nsize_type(i))
           IF (n_quiet == 0) THEN
              WRITE(*, FMT='(A40,5x,1P,2(E9.2,1x),1x,I3)') gtype(i), 1.0e-4_dp*size_type(i,1), &
                   1.0e-4_dp*size_type(i,(nsize_type(i))), nsize_type(i)
           ENDIF

           ! get Q1abs
           the_char = READ_COM (fq1ext%unit)
           BACKSPACE fq1ext%unit   ! READ_COM only reads 1st column
           DO k=1,n_qabs
              READ (UNIT=fq1ext%unit, FMT=*) (q1_abs(i,u,k), u=1,nsize_type(i))
           ENDDO

           ! get Q1sca
           the_char = READ_COM (fq1ext%unit)
           BACKSPACE fq1ext%unit   ! READ_COM only reads 1st column
           DO k=1,n_qabs
              READ (UNIT=fq1ext%unit, FMT=*) (q1diff(i,u,k), u=1,nsize_type(i))
           ENDDO
           CLOSE (UNIT=fq1ext%unit)

           ! reading Q2_*.DAT file
           IF (n_quiet == 0) THEN
              WRITE (*,*)
              WRITE (*, FMT='(A18,31x,A22)') '>>read Q2_TYPE.DAT', 'a-range (cm)     nsize'
           ENDIF
           filename_tmp = TRIMCAT (data_path, dir_qabs)
           filename_tmp = TRIMCAT (filename_tmp, fq2ext%nom)
           filename_tmp = TRIMCAT (filename_tmp, gtype(i))
           filename_tmp = TRIMCAT (filename_tmp, '.DAT')
           OPEN (UNIT=fq2ext%unit, FILE=filename_tmp, STATUS='old')

           ! read doc
           the_char = READ_COM (fq2ext%unit)
           READ (the_char, FMT=*) nsz1(i)

           ! checking Q file
           IF (nsz1(i) > nsize_max_qabs) THEN
              WRITE (*,*) ''
              WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', TRIM(filename_tmp)
              WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
              WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
              WRITE (*,*) ''
              STOP
           ENDIF
           IF (nsz1(i) /= nsize_type(i)) THEN
              WRITE (*,*) ''
              WRITE (*,*) '(F) DM_inout/READ_DATA: odd number of sizes in ', TRIM(filename_tmp)
              WRITE (*,*) '                        found', nsz1(i), ' sizes but ', nsize_type(i), ' sizes expected'
              WRITE (*,*) '                        Q, C and G files must have same number of sizes'
              WRITE (*,*) ''
              STOP
           ENDIF
           nsize_type(i) = nsz1(i)

           ! get sizes
           READ (UNIT=fq2ext%unit, FMT=*) (size_type(i,u), u=1,nsize_type(i))
           IF (n_quiet == 0) THEN
              WRITE(*, FMT='(A40,5x,1P,2(E9.2,1x),1x,I3)') gtype(i), 1.0e-4_dp*size_type(i,1), &
                   1.0e-4_dp*size_type(i,(nsize_type(i))), nsize_type(i)
           ENDIF

           ! get Q2abs
           the_char = READ_COM (fq2ext%unit)
           BACKSPACE fq2ext%unit   ! READ_COM only reads 1st column
           DO k=1,n_qabs
              READ (UNIT=fq2ext%unit, FMT=*) (q2_abs(i,u,k), u=1,nsize_type(i))
           ENDDO

           ! get Q2sca
           the_char = READ_COM (fq2ext%unit)
           BACKSPACE fq2ext%unit   ! READ_COM only reads 1st column
           DO k=1,n_qabs
              READ (UNIT=fq2ext%unit, FMT=*) (q2diff(i,u,k), u=1,nsize_type(i))
           ENDDO
           CLOSE (UNIT=fq2ext%unit)

           DO j=1,nsize(i)
              aux  = 1.e4_dp*size_ava(i,j)
              CALL GET_QEXT_POL(aux, i, tmp1,tmp2,tmp3,tmp4,tmp5)
              q1i_abs(i,j,:) = tmp1(:)
              q1idiff(i,j,:) = tmp2(:)
              q2i_abs(i,j,:) = tmp3(:)
              q2idiff(i,j,:) = tmp4(:)
              qiH_abs(i,j,:) = tmp5(:)
           ENDDO

	   ! Read circular polarization
	   if (n_circ > 0) then 

              IF (n_quiet == 0) THEN
              	 WRITE (*,*)
                 WRITE (*, FMT='(A18,31x,A22)') '>>read Qc_TYPE.DAT', 'a-range (cm)     nsize'
              ENDIF
              filename_tmp = TRIMCAT (data_path, dir_qabs)
              filename_tmp = TRIMCAT (filename_tmp, fqcirc%nom)
              filename_tmp = TRIMCAT (filename_tmp, gtype(i))
              filename_tmp = TRIMCAT (filename_tmp, '.DAT')
              OPEN (UNIT=fqcirc%unit, FILE=filename_tmp, STATUS='old')

              ! read doc
              the_char = READ_COM (fqcirc%unit)
              READ (the_char, FMT=*) nsz1(i)

              ! checking Q file
              IF (nsz1(i) > nsize_max_qabs) THEN
                 WRITE (*,*) ''
                 WRITE (*,*) '(F) DM_inout/READ_DATA: too many sizes required in ', TRIM(filename_tmp)
                 WRITE (*,*) '                        nsize_max_qabs=', nsize_max_qabs, ' while nsize_type=', nsize_type(i)
                 WRITE (*,*) '                        nsize_max_qabs PARAMETER can be changed in DM_utility.f90'
                 WRITE (*,*) ''
                 STOP
              ENDIF
              IF (nsz1(i) /= nsize_type(i)) THEN
                 WRITE (*,*) ''
                 WRITE (*,*) '(F) DM_inout/READ_DATA: odd number of sizes in ', TRIM(filename_tmp)
                 WRITE (*,*) '                        found', nsz1(i), ' sizes but ', nsize_type(i), ' sizes expected'
                 WRITE (*,*) '                        Q, C and G files must have same number of sizes'
                 WRITE (*,*) ''
                 STOP
              ENDIF
              nsize_type(i) = nsz1(i)

              ! get sizes
              READ (UNIT=fqcirc%unit, FMT=*) (size_type(i,u), u=1,nsize_type(i))
              IF (n_quiet == 0) THEN
                 WRITE(*, FMT='(A40,5x,1P,2(E9.2,1x),1x,I3)') gtype(i), 1.0e-4_dp*size_type(i,1), &
                      1.0e-4_dp*size_type(i,(nsize_type(i))), nsize_type(i)
              ENDIF
   
              ! get Qcirc
              the_char = READ_COM (fqcirc%unit)
              BACKSPACE fqcirc%unit   ! READ_COM only reads 1st column
              DO k=1,n_qabs
                 READ (UNIT=fqcirc%unit, FMT=*) (qcirc(i,u,k), u=1,nsize_type(i))
              ENDDO
              CLOSE (UNIT=fqcirc%unit)

              DO j=1,nsize(i)
                 aux  = 1.e4_dp*size_ava(i,j)
                 CALL GET_QEXT_CIRC(aux, i, tmp1)
                 qicirc(i,j,:) = tmp1(:)
              ENDDO

	   endif

           ! convert sizes microns --> cm after POL read
           size_type(i,:) = size_type(i,:) * 1.0e-4_dp ! microns --> cm

        ENDIF

     ENDDO ! type loop

     DEALLOCATE(qH1_abs)
     DEALLOCATE(qH2_abs)

  ENDIF

  ! reading radiation field (ISRF.DAT)
  !-----------------------------------------------------------------
  filename_tmp = TRIMCAT (data_path, dir_dat)
  filename_tmp = TRIMCAT (filename_tmp, fisrf%nom)
  OPEN (UNIT=fisrf%unit, FILE=filename_tmp, STATUS='old')
  the_char = READ_COM (fisrf%unit)
  READ (the_char, FMT=*) nisrf
  ALLOCATE (lambisrf(nisrf))
  ALLOCATE (isrf(nisrf))
  DO k=1,nisrf
     READ (unit=fisrf%unit, FMT=*) lambisrf(k), isrf(k)
  ENDDO
  isrf = g0 * isrf
  CLOSE (UNIT=fisrf%unit)
  lambisrf = lambisrf * 1.0e-4_dp

  ! ISRF interpolation - extrapolation outside lambisrf bounds is set to 0
  ! beware of flux=0 above lyman limit
  ALLOCATE (isrfuv(n_qabs))
  flag_cut = 0
  jfreqmax = 1
  DO i=1,n_qabs
     isrfuv(i) = INTPOL (isrf, lambisrf, nisrf, lamb_qabs(i))
     IF (flag_cut == 0 .AND. isrfuv(i) <= istiny) THEN
       jfreqmax = i
     ELSE IF (flag_cut == 0 .AND. isrfuv(i) > istiny) THEN
       flag_cut = 1
     ENDIF
  ENDDO
  jfreqmax = jfreqmax + 1
  tmp1 = isrfuv
  DO i=1,n_qabs
     isrfuv(i) = tmp1(n_qabs-i+1)
  ENDDO
  jfreqmax = n_qabs -jfreqmax +1
  hnumin = xhp * freq_qabs(1)
  hnumax = xhp * freq_qabs(jfreqmax)
  IF (n_quiet == 0) THEN
     WRITE (*,*)
     WRITE (*, FMT='(A64)') '>>read ISRF.DAT            w-range (microns)  hnu_max(eV)  nwave'
     WRITE (*, FMT='(A25,1P,3(E9.2,1x),3x,I4)') '', lambisrf(1)*1.0e4_dp, lambisrf(nisrf)*1.0e4_dp, hnumax/everg, nisrf
  ENDIF

  ! allocate arrays for grain emission
  !-----------------------------------------------------------------
  ALLOCATE (nuinuemtot(n_qabs))
  ALLOCATE (nuinuem(ntype,n_qabs))
  nuinuem = 0.0_dp
  nuinuemtot = 0.0_dp

  IF (n_spin > 0) THEN
     ALLOCATE (spnuinuem(ntype,n_qabs))
     ALLOCATE (spnuinuemtot(n_qabs))
     spnuinuem = 0.0_dp
     spnuinuemtot = 0.0_dp
  ENDIF

  IF (n_pol > 0 .and. n_lin > 0) THEN
     ALLOCATE (nuinuemp(ntype,n_qabs))
     ALLOCATE (nuinuemptot(n_qabs))
     nuinuemp = 0.0_dp
     nuinuemptot = 0.0_dp
  ENDIF

  DEALLOCATE (lambisrf)
  DEALLOCATE (isrf)
  DEALLOCATE (tmp1,tmp2,tmp3)
  IF (n_pol > 0) DEALLOCATE (tmp4,tmp5)

END SUBROUTINE READ_DATA

!----------------------------------------------------------------

SUBROUTINE WRITE_DATA(fsed)
! WRITEs DUSTEM outputs integrated over all sizes

  !global variables modules
  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT
  USE MGET_TDIST
  USE MCOMPUTE

  IMPLICIT NONE

  TYPE(FICH)                  :: fsed  ! = FICH ('SED.RES', 18)
  TYPE(FICH)                  :: fpsed  = FICH ('SED_POL.RES', 19)
  TYPE(FICH)                  :: fext   = FICH ('EXT.RES', 22)
  TYPE(FICH)                  :: fpext  = FICH ('EXT_POL.RES', 23)
  TYPE(FICH)                  :: fsdist = FICH ('SDIST.RES', 24)
  TYPE (FICH)                 :: fspd   = FICH ('SPIN.RES', 25)
  TYPE(FICH)                  :: fpcirc = FICH ('EXT_CIRC.RES', 26)
  TYPE(FICH)                  :: fHabs  = FICH ('ABS_HEATING.RES', 27)
  TYPE(FICH)                  :: fpdr_e = FICH ('EXTINCTION_DUSTEM.RES', 32)

  CHARACTER (LEN=max_len)     :: filename_tmp

  INTEGER                     :: i,j,k
  REAL (KIND=dp), ALLOCATABLE :: abscsuv(:,:)               ! absorption cross-section cm2/gram integrated over sizes
  REAL (KIND=dp), ALLOCATABLE :: diffcsuv(:,:)              ! scattering cross-section cm2/gram integrated over sizes
  REAL (KIND=dp), ALLOCATABLE :: circcsuv(:,:)              ! Circular polarization cross-section cm2/gram integrated over sizes
  REAL (KIND=dp), ALLOCATABLE :: absH(:,:)                  ! absorption cross-section cm2/gram for grain heating
  REAL (KIND=dp), ALLOCATABLE :: ggfac(:,:)                 ! g-factor for scattering integrated over sizes
  REAL (KIND=dp), ALLOCATABLE :: sdiff(:,:)                 ! weight for scattering integrated over sizes
  REAL (KIND=dp), ALLOCATABLE :: tmp1(:), tmp2(:), tmp3(:), tmp4(:), tmp5(:)
  REAL(KIND=DP)               :: fact

 !---------- WRITE THE OUTPUT FILES -----------

  ALLOCATE (abscsuv(ntype,n_qabs))
  ALLOCATE (diffcsuv(ntype,n_qabs))
  ALLOCATE (circcsuv(ntype,n_qabs))
  ALLOCATE (absH(ntype,n_qabs))
  ALLOCATE (ggfac(ntype,n_qabs))
  ALLOCATE (sdiff(ntype,n_qabs))
  ALLOCATE (tmp1(n_qabs),tmp2(n_qabs),tmp3(n_qabs),tmp4(n_qabs),tmp5(n_qabs))

  filename_tmp = TRIMCAT(data_path,dir_res)
  filename_tmp = TRIMCAT(filename_tmp,fsed%nom)
  OPEN  (UNIT=fsed%unit, FILE=filename_tmp, STATUS='unknown')
  WRITE (UNIT=fsed%unit, FMT='(A40)')            '# DUSTEM SED:  4*pi*nu*I_nu/NH (erg/s/H)'
  WRITE (UNIT=fsed%unit, FMT='(A1)')             '#'
  WRITE (UNIT=fsed%unit, FMT='(A14)')            '# Grain types '
  WRITE (UNIT=fsed%unit, FMT='(A34)')            '# nr of grain types   nr of lambda'
  WRITE (UNIT=fsed%unit, FMT='(A52)')            '# lambda (microns)   SED(1)...SED(ntype)   SED total'
  WRITE (UNIT=fsed%unit, FMT='(A1)')             '#'
  WRITE (UNIT=fsed%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
  WRITE (UNIT=fsed%unit, FMT='(I2,2x,I4)')       ntype, n_qabs
  DO k=1,n_qabs
     WRITE (UNIT=fsed%unit, FMT='(1P,25E16.6E3)') lamb_qabs(k)*1.0e4_dp, (nuinuem(i,k), i= 1,ntype), &
                                                      nuinuemtot(k)
  ENDDO
  CLOSE (UNIT=fsed%unit)

! get and write extinction
  DO i=1,ntype
     CALL EXTINCTION (i, tmp1, tmp2, tmp3, tmp4, tmp5)
     abscsuv(i,:)  = tmp1(:)
     diffcsuv(i,:) = tmp2(:)
     ggfac(i,:)    = tmp3(:)
     sdiff(i,:)    = tmp4(:)
     absH(i,:)     = tmp5(:)
  ENDDO

!VG Conversion factor from optical depth per gram to optical depth per NH (for NH=1e21 H/cm2)
  fact = xmp * 1.0e21_dp

  filename_tmp = TRIMCAT(data_path,dir_res)
  filename_tmp = TRIMCAT(filename_tmp,fext%nom)
  OPEN  (UNIT=fext%unit, FILE=filename_tmp, STATUS='unknown')
  ! VG : extinction per H, not per gram, like for SED, to avoid reading GRAIN.DAT for plotting
  ! WRITE (UNIT=fext%unit, FMT='(A43)')          '# DUSTEM extinction per mass: sigma (cm2/g)'
  WRITE (UNIT=fext%unit, FMT='(A78)')            '# DUSTEM extinction cross-section for NH=10^21 H/cm2 : sigma (x 1e-21) [cm2/H]'
  WRITE (UNIT=fext%unit, FMT='(A2)')             '# '
  WRITE (UNIT=fext%unit, FMT='(A14)')            '# Grain types '
  WRITE (UNIT=fext%unit, FMT='(A34)')            '# nr of grain types - nr of lambda'
  WRITE (UNIT=fext%unit, FMT='(A77)')            '# lambda (microns)   ABS(1) ... ABS(ntype)  SCA(1) ... SCA(ntype)  EXT(Total)'
  WRITE (UNIT=fext%unit, FMT='(A2)')             '# '
  WRITE (UNIT=fext%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
  WRITE (UNIT=fext%unit, FMT='(I2,2X,I4)')       ntype, n_qabs
  DO k = 1,n_qabs
     !VG: per H, not per gram
     !WRITE (UNIT=fext%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, (abscsuv(i,n_qabs-k+1), i=1,ntype), &
     !                                               (diffcsuv(i,n_qabs-k+1), i=1,ntype)
     WRITE (UNIT=fext%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, (fact*mprop(i)* abscsuv(i,n_qabs-k+1), i=1,ntype), &
                                                                         (fact*mprop(i)*diffcsuv(i,n_qabs-k+1), i=1,ntype), &
!VG : add total extinction
						  fact*sum(mprop(:)*(abscsuv(:,n_qabs-k+1)+diffcsuv(:,n_qabs-k+1)))
  ENDDO
  CLOSE (UNIT=fext%unit)

  IF (anisG0 .gt. 0) THEN

     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = TRIMCAT(filename_tmp,fHabs%nom)
     OPEN  (UNIT=fHabs%unit, FILE=filename_tmp, STATUS='unknown')
     WRITE (UNIT=fHabs%unit, FMT='(A28)')            '# DUSTEM absorption for heating: '
     WRITE (UNIT=fHabs%unit, FMT='(I3)') ntype
     WRITE (UNIT=fHabs%unit, FMT='(A2)')             '# '
     WRITE (UNIT=fHabs%unit, FMT='(A14)')            '# Grain types '
     WRITE (UNIT=fHabs%unit, FMT='(A34)')            '# nr of grain types - nr of lambda'
     WRITE (UNIT=fHabs%unit, FMT='(A53)')            '# lambda (microns)   ABS(1) ... ABS(ntype) ABS(Total)'
     WRITE (UNIT=fHabs%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
     WRITE (UNIT=fHabs%unit, FMT='(I2,2X,I4)')       ntype, n_qabs
     DO k = 1,n_qabs
        WRITE (UNIT=fHabs%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, (fact*mprop(i)* absH(i,n_qabs-k+1), i=1,ntype), &
						  fact*sum(mprop(:)*(absH(:,n_qabs-k+1)))
     ENDDO
     CLOSE (UNIT=fHabs%unit)

  ENDIF

  IF (n_sdist /= 0) THEN
     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = TRIMCAT(filename_tmp,fsdist%nom)
     OPEN  (UNIT=fsdist%unit, FILE=filename_tmp, STATUS='unknown')
     WRITE (UNIT=fsdist%unit, FMT='(A28)')            '# DUSTEM size distribution: '
     WRITE (UNIT=fsdist%unit, FMT='(A2)')             '# '
     WRITE (UNIT=fsdist%unit, FMT='(A19)')            '# nr of grain types'
     WRITE (UNIT=fsdist%unit, FMT='(A26)')            '# grain type - nr of sizes'
     WRITE (UNIT=fsdist%unit, FMT='(A36)')            '# size (cm) - a**3 * dn/dlna (cm3/H)'
     WRITE (UNIT=fsdist%unit, FMT='(A2)')             '# '
     WRITE (UNIT=fsdist%unit, FMT='(I3)') ntype
     DO i=1,ntype
             WRITE (UNIT=fsdist%unit, FMT='(A40,1x,I3)') gtype(i), nsize(i)
             DO j=1, nsize(i)
                WRITE (UNIT=fsdist%unit, FMT='(2(1PE15.6E3,1x))') size_ava(i,j), &
                     & ava(i,j)*mprop(i)*xmp/masstot(i)
             ENDDO
     ENDDO
     CLOSE (UNIT=fsdist%unit)
  ENDIF

  IF (n_pdr_o /= 0) THEN
  ! write transfer stuff: cross-sections, SED per H atom and mean g-factor
     filename_tmp = TRIMCAT(dir_PDR,fpdr_e%nom)
     OPEN (UNIT=fpdr_e%unit, FILE=filename_tmp, STATUS='unknown')
     WRITE(UNIT=fpdr_e%unit,FMT='(1i5)') n_qabs
     ! 6 columns: Wavelength, Absorption, Scattering, Emission, Albedo, g (Mean of cos(theta))
     DO k = 1,n_qabs
        WRITE (UNIT=fpdr_e%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, &
              SUM( (/ (abscsuv(i,n_qabs-k+1)*mprop(i)*xmp, i=1,ntype) /) ), &
              SUM( (/ (diffcsuv(i,n_qabs-k+1)*mprop(i)*xmp, i=1,ntype) /) ), &
              nuinuemtot(k), &
              SUM( (/ ((diffcsuv(i,n_qabs-k+1))*mprop(i)*xmp, i=1,ntype) /) )/ &
              & SUM( (/ ((abscsuv(i,n_qabs-k+1)+diffcsuv(i,n_qabs-k+1))*mprop(i)*xmp, i=1,ntype) /) ), &
              SUM( (/ (ggfac(i,k), i=1,ntype) /) ) / SUM( (/ (sdiff(i,k), i=1,ntype) /) )
     ENDDO
     CLOSE (UNIT=fpdr_e%unit)
  ENDIF

  IF (n_spin > 0) THEN
     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = TRIMCAT(filename_tmp,fspd%nom)
     OPEN  (UNIT=fspd%unit, FILE=filename_tmp, STATUS='unknown')
     WRITE (UNIT=fspd%unit, FMT='(A50)')            '# DUSTEM spinning SED:  4*pi*nu*I_nu/NH (erg/s/H)'
     WRITE (UNIT=fspd%unit, FMT='(A1)')             '#'
     WRITE (UNIT=fspd%unit, FMT='(A14)')            '# Grain types '
     WRITE (UNIT=fspd%unit, FMT='(A34)')            '# nr of grain types   nr of lambda'
     WRITE (UNIT=fspd%unit, FMT='(A52)')            '# lambda (microns)   SED(1)...SED(ntype)   SED total'
     WRITE (UNIT=fspd%unit, FMT='(A1)')             '#'
     WRITE (UNIT=fspd%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
     WRITE (UNIT=fspd%unit, FMT='(I2,2x,I4)')       ntype, n_qabs
     DO k=1,n_qabs
        WRITE (UNIT=fspd%unit, FMT='(1P,25E16.6E3)') lamb_qabs(k)*1.0e4_dp, (spnuinuem(i,k), i= 1,ntype), &
                                                      spnuinuemtot(k)
     ENDDO
     CLOSE (UNIT=fspd%unit)
  ENDIF

! get and write polarized SED and extinction 
  IF (n_pol > 0 .and. n_lin > 0) THEN

     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = TRIMCAT(filename_tmp,fpsed%nom)
     OPEN  (UNIT=fpsed%unit, FILE=filename_tmp, STATUS='unknown')
     WRITE (UNIT=fpsed%unit, FMT='(A50)')            '# DUSTEM polarized SED:  4*pi*nu*I_nu/NH (erg/s/H)'
     WRITE (UNIT=fpsed%unit, FMT='(A1)')             '#'
     WRITE (UNIT=fpsed%unit, FMT='(A14)')            '# Grain types '
     WRITE (UNIT=fpsed%unit, FMT='(A34)')            '# nr of grain types   nr of lambda'
     WRITE (UNIT=fpsed%unit, FMT='(A52)')            '# lambda (microns)   SED(1)...SED(ntype)   SED total'
     WRITE (UNIT=fpsed%unit, FMT='(A1)')             '#'
     WRITE (UNIT=fpsed%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
     WRITE (UNIT=fpsed%unit, FMT='(I2,2x,I4)')       ntype, n_qabs
     DO k=1,n_qabs
        WRITE (UNIT=fpsed%unit, FMT='(1P,25E16.6E3)') lamb_qabs(k)*1.0e4_dp, (nuinuemp(i,k), i= 1,ntype), &
                                                      nuinuemptot(k)
     ENDDO
     CLOSE (UNIT=fpsed%unit)

     DO i=1,ntype
        CALL EXTINCTION_POL (i, tmp1, tmp2)
        abscsuv(i,:)  = tmp1(:)
        diffcsuv(i,:) = tmp2(:)
     ENDDO
     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = TRIMCAT(filename_tmp,fpext%nom)
     OPEN  (UNIT=fpext%unit, FILE=filename_tmp, STATUS='unknown')
     !VG : per H, not per gram
     !WRITE (UNIT=fpext%unit, FMT='(A57)') '# DUSTEM polarized extinction per mass: sigma_pol (cm2/g)'
     !WRITE (UNIT=fpext%unit, FMT='(A61)') '# DUSTEM polarized extinction per hydrogen: sigma_pol (cm2/H)'
     WRITE (UNIT=fpext%unit, FMT='(A98)') &
            adjustl('# DUSTEM polarization cross-section in extinction for NH=10^21 H/cm2 : sigma_pol (x 1e-21) [cm2/H]')
     WRITE (UNIT=fpext%unit, FMT='(A2)')   '# '
     WRITE (UNIT=fpext%unit, FMT='(A14)')  '# Grain types '
     WRITE (UNIT=fpext%unit, FMT='(A34)')  '# nr of grain types - nr of lambda'
     WRITE (UNIT=fpext%unit, FMT='(A77)')  '# lambda (microns)   ABS(1) ... ABS(ntype)  SCA(1) ... SCA(ntype)  PEXT Total'
     WRITE (UNIT=fpext%unit, FMT='(A2)')   '# '
     WRITE (UNIT=fpext%unit, FMT='(A2,10(A,X))')     '# ', (trim(gtype(i)), i=1,ntype)
     WRITE (UNIT=fpext%unit, FMT='(I2,2X,I4)')       ntype, n_qabs
     DO k = 1,n_qabs
        WRITE (UNIT=fpext%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, (fact*mprop(i)* abscsuv(i,n_qabs-k+1), i=1,ntype), &
                                                                             (fact*mprop(i)*diffcsuv(i,n_qabs-k+1), i=1,ntype), &
						  fact*sum(mprop(:)*(abscsuv(:,n_qabs-k+1)+diffcsuv(:,n_qabs-k+1)))
     ENDDO
     CLOSE (UNIT=fpext%unit)

  ENDIF

  IF (n_pol > 0 .and. n_circ > 0) THEN
     	DO i=1,ntype
        	CALL EXTINCTION_CIRC(i, tmp1)
        	circcsuv(i,:)  = tmp1(:)
     	ENDDO
     	filename_tmp = TRIMCAT(data_path,dir_res)
     	filename_tmp = TRIMCAT(filename_tmp,fpcirc%nom)
     	OPEN  (UNIT=fpcirc%unit, FILE=filename_tmp, STATUS='unknown')
     	!VG : per H, not per gram
     	WRITE (UNIT=fpcirc%unit, FMT='(A63)')        '# DUSTEM circular polarization cross-section for NH=10^21 H/cm2'
     	WRITE (UNIT=fpcirc%unit, FMT='(A2)')         '# '
     	WRITE (UNIT=fpcirc%unit, FMT='(A14)')        '# Grain types '
     	WRITE (UNIT=fpcirc%unit, FMT='(A34)')        '# nr of grain types - nr of lambda'
     	WRITE (UNIT=fpcirc%unit, FMT='(A50)')        '# lambda (microns)   CIRC(1) ... CIRC(ntype) TOTAL'
     	WRITE (UNIT=fpcirc%unit, FMT='(A2)')         '# '
     	WRITE (UNIT=fpcirc%unit, FMT='(A19,50(A14,1X))') '# lambda (microns) ', ("C"//gtype(i), i=1,ntype),'PCIRC Total        '
     	WRITE (UNIT=fpcirc%unit, FMT='(I2,2X,I4)')       ntype, n_qabs
     	DO k = 1,n_qabs
        	WRITE (UNIT=fpcirc%unit, FMT='(1P,50E15.6E3)') lamb_qabs(k)*1.0e4_dp, (fact*mprop(i)* circcsuv(i,n_qabs-k+1), i=1,ntype), &
						  	fact*sum(mprop(:)*(circcsuv(:,n_qabs-k+1)))
     	ENDDO
     	CLOSE (UNIT=fpcirc%unit)

  ENDIF

  DEALLOCATE (nuinuemtot, nuinuem, isrfuv)
  DEALLOCATE (calor, temp_cal, n_temp, lamb_qabs, freq_qabs, lfrq_qabs, q_abs, qdiff)
  DEALLOCATE (ggfac, sdiff)
  DEALLOCATE (tmp1, tmp2, tmp3, tmp4, tmp5)
  IF (n_spin > 0) DEALLOCATE (spnuinuemtot, spnuinuem)
  IF (n_pol > 0 .and. n_lin > 0) THEN
     DEALLOCATE (q1_abs, q1diff,q2_abs, q2diff)
     DEALLOCATE (nuinuemptot, nuinuemp)
  ENDIF
  DEALLOCATE (rho, rhom, gtype, mprop, as2, as1, t_opt, nsize)
  DEALLOCATE (size_ava, ava, f_mix, f_pol, dloga, sect_eff, si_ava_l, masstot)
  DEALLOCATE (nsize_type, size_type)

END SUBROUTINE WRITE_DATA

!----------------------------------------------------------------

SUBROUTINE EXTINCTION (nt, absuv, diffuv, ggfac, sdiff, absH)
! for each grain type NT, computes absorption and scattering cross-section 
! integrated over size distribution. Result in cm2/gram

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  ! in-out arguments
  INTEGER, INTENT (IN)         :: nt                 ! index for grain type
  REAL (KIND=dp), INTENT (OUT) :: absuv(n_qabs)      ! absorption cross-section in cm2/g
  REAL (KIND=dp), INTENT (OUT) :: diffuv(n_qabs)     ! scattering cross-section in cm2/g
  REAL (KIND=dp), INTENT (OUT) :: ggfac(n_qabs)      ! g-factor for scattering
  REAL (KIND=dp), INTENT (OUT) :: sdiff(n_qabs)      ! weight for scattering
  REAL (KIND=dp), INTENT (OUT) :: absH(n_qabs)       ! absorption for Heating

  ! local
  INTEGER                      :: i, ns_i
  REAL (KIND=dp)               :: yy(nsize_max), xx(nsize_max),zz(nsize_max)

  ! inits
  absuv = 0.0_dp
  diffuv= 0.0_dp
  ggfac = 0.0_dp
  absH  = 0.0_dp

  ns_i = nsize(nt)

  ! sum over sizes
    DO i=1,n_qabs
    xx(1:ns_i) = si_ava_l(nt,1:ns_i)
    IF (n_anis .eq. 0 .or. INDEX(t_opt(nt),'POL') == 0) THEN
	! Polarization is not calculated : grains are not aligned
    	yy(1:ns_i) = sect_eff(nt,1:ns_i) * qi_abs(nt,1:ns_i,i) * ava(nt,1:ns_i) * f_mix(nt,1:ns_i) &
                   / size_ava(nt,1:ns_i)**3 / masstot(nt)
    else
	!VG : extinction depends on alignment efficiency 
    	yy(1:ns_i) = sect_eff(nt,1:ns_i) &
               * ( (q1i_abs(nt,1:ns_i,i) + q2i_abs(nt,1:ns_i,i))/2.0_dp * f_pol(nt,1:ns_i) &
                 + qi_abs(nt,1:ns_i,i) * (1._dp - f_pol(nt,1:ns_i) ) ) &
               * ava(nt,1:ns_i) * f_mix(nt,1:ns_i)   &
               / size_ava(nt,1:ns_i)**3 / masstot(nt) 
    endif
    absuv(i)  = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))

    IF (anisG0 .gt. 0) THEN
        ! Grains are aligned and radiation field is anisotropic with degree anisG0
        zz(1:ns_i) = sect_eff(nt,1:ns_i) &
                   * (anisg0 * (    f_pol(nt,1:ns_i)  * qiH_abs(nt,1:ns_i,i)  &
                               + (1-f_pol(nt,1:ns_i)) * qi_abs(nt,1:ns_i,i) ) &
                     + (1-anisg0) *                     qi_abs(nt,1:ns_i,i) ) &
                   * ava(nt,1:ns_i) * f_mix(nt,1:ns_i)                        &
                   / size_ava(nt,1:ns_i)**3 / masstot(nt)
        absH(i)  = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), zz(1:ns_i))
    ELSE
        ! Radiation field is isotropic
        absH(i) = absuv(i)
     ENDIF

    IF (n_anis .eq. 0 .or. INDEX(t_opt(nt),'POL') == 0) THEN
    !if (.TRUE.) then
	! Polarization is not calculated : grains are not aligned
    	yy(1:ns_i) = sect_eff(nt,1:ns_i) * qidiff(nt,1:ns_i,i) * ava(nt,1:ns_i) * f_mix(nt,1:ns_i) / &
                & size_ava(nt,1:ns_i)**3 / masstot(nt)
    else
	!VG : extinction depends on alignment efficiency 
    	yy(1:ns_i) = sect_eff(nt,1:ns_i) &
               * ( (q1idiff(nt,1:ns_i,i) + q2idiff(nt,1:ns_i,i))/2.0_dp * f_pol(nt,1:ns_i) &
                 + qidiff(nt,1:ns_i,i) * (1._dp - f_pol(nt,1:ns_i) ) ) &
               * ava(nt,1:ns_i) * f_mix(nt,1:ns_i)   &
               / size_ava(nt,1:ns_i)**3 / masstot(nt)
    endif
    diffuv(i) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))

! NB: qi_abs and qi_diff are reversed wrt lamb_qabs (CALL GET_QEXT) but NOT gifac !!!
    yy(1:ns_i) = gifac(nt,1:ns_i,i)*sect_eff(nt,1:ns_i)*qidiff(nt,1:ns_i,n_qabs-i+1)*f_mix(nt,1:ns_i)* &
                & ava(nt,1:ns_i) / size_ava(nt,1:ns_i)**3
    ggfac(i)  = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
    yy(1:ns_i) = sect_eff(nt,1:ns_i)*qidiff(nt,1:ns_i,n_qabs-i+1)*ava(nt,1:ns_i)*f_mix(nt,1:ns_i)/ &
                & size_ava(nt,1:ns_i)**3
    sdiff(i) =  XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))

  ENDDO

END SUBROUTINE EXTINCTION


SUBROUTINE EXTINCTION_CIRC (nt, circuv)

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  ! in-out arguments
  INTEGER, INTENT (IN)         :: nt                 ! index for grain type
  REAL (KIND=dp), INTENT (OUT) :: circuv(n_qabs)      ! circular polarization cross-section in cm2/g

  ! local
  INTEGER                      :: i, ns_i
  REAL (KIND=dp)               :: yy(nsize_max), xx(nsize_max)

  ! inits
  circuv = 0.0_dp

! sum over sizes
  DO i=1,n_qabs
     ns_i = nsize(nt)
     xx(1:ns_i) = si_ava_l(nt,1:ns_i)
     yy(1:ns_i) = sect_eff(nt,1:ns_i) * qicirc(nt,1:ns_i,i) * f_pol(nt,1:ns_i) * & 
          & ava(nt,1:ns_i) * f_mix(nt,1:ns_i) / size_ava(nt,1:ns_i)**3 / masstot(nt)
     circuv(i)  = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
  ENDDO

END SUBROUTINE EXTINCTION_CIRC 

SUBROUTINE EXTINCTION_POL (nt, absuv, diffuv)

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  ! in-out arguments
  INTEGER, INTENT (IN)         :: nt                 ! index for grain type
  REAL (KIND=dp), INTENT (OUT) :: absuv(n_qabs)      ! absorption cross-section in cm2/g
  REAL (KIND=dp), INTENT (OUT) :: diffuv(n_qabs)     ! scattering cross-section in cm2/g

  ! local
  INTEGER                      :: i, ns_i
  REAL (KIND=dp)               :: yy(nsize_max), xx(nsize_max)

  ! inits
  absuv = 0.0_dp
  diffuv= 0.0_dp

! sum over sizes
  DO i=1,n_qabs
     ns_i = nsize(nt)
     xx(1:ns_i) = si_ava_l(nt,1:ns_i)
     yy(1:ns_i) = sect_eff(nt,1:ns_i) * ((q2i_abs(nt,1:ns_i,i)-q1i_abs(nt,1:ns_i,i))/2.0_dp)* f_pol(nt,1:ns_i) * & 
          & ava(nt,1:ns_i) * f_mix(nt,1:ns_i) / size_ava(nt,1:ns_i)**3 / masstot(nt)
     absuv(i)  = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
     yy(1:ns_i) = sect_eff(nt,1:ns_i) * ((q2idiff(nt,1:ns_i,i)-q1idiff(nt,1:ns_i,i))/2.0_dp)* f_pol(nt,1:ns_i) * &
          & ava(nt,1:ns_i) * f_mix(nt,1:ns_i) / size_ava(nt,1:ns_i)**3 / masstot(nt)
     diffuv(i) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
  ENDDO

END SUBROUTINE EXTINCTION_POL


END MODULE IN_OUT
