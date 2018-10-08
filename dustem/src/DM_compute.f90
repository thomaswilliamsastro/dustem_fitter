MODULE MCOMPUTE

  USE CONSTANTS
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: COMPUTE

CONTAINS

! ------------------------------------------------------------------------------------

SUBROUTINE COMPUTE

  ! global variables modules

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT
  USE MGET_TDIST
  USE MZDIST
  USE MSPIN
  USE MDTLS

  ! local variables

  IMPLICIT NONE

  TYPE (FICH)                 :: ftemp   = FICH('TEMP.RES',21)
  TYPE (FICH)                 :: fzdist  = FICH('ZDIST.RES',26)
  TYPE (FICH)                 :: fsed_a  = FICH('SED_A.RES',31)
  TYPE (FICH)                 :: fpsed_a = FICH('SED_POL_A.RES',32)
  TYPE (FICH)                 :: fspd_a  = FICH('SPIN_A.RES',33)
  CHARACTER (LEN=max_len)     :: filename_tmp

  ! loop indices
  INTEGER                     :: i             ! grain type index
  INTEGER                     :: j             ! grain size index
  INTEGER                     :: kt

  REAL (KIND=dp)              :: enerabs(nsize_max)  ! power absorbed by grain
  REAL (KIND=dp)              :: enerem(nsize_max)   ! power emitted by grain
  REAL (KIND=dp)              :: tempmoy(nsize_max)  ! grain mean temperature 
  REAL (KIND=dp)              :: tempequi(nsize_max) ! grain equilibrium temperature 
  REAL (KIND=dp)              :: tempmax(nsize_max)  ! grain max temperature 
  REAL (KIND=dp)              :: t_rot(nsize_max)    ! rotational temperature

  REAL (KIND=dp)              :: temp(ndist)   ! temperature grid
  REAL (KIND=dp)              :: t2(ndist)     ! LOG(T)
  REAL (KIND=dp)              :: hcap(ndist)   ! heat capacity
  REAL (KIND=dp)              :: dpt(ndist)    ! dP/dln T
  REAL (KIND=dp)              :: uint(ndist)   ! internal energy
  REAL (KIND=dp)              :: dpu(ndist)    ! dP/dU

  ! for dust IR (vibrational) SED
  REAL (KIND=dp), ALLOCATABLE :: nuflux_j(:)   ! flux*frequency for 1 type and 1 size
  REAL (KIND=dp), ALLOCATABLE :: nuflux(:,:)   ! flux*frequency for 1 type and all sizes

  ! for dust polarized emission
  REAL (KIND=dp), ALLOCATABLE :: nufluxp_j(:)  ! flux*frequency for 1 type and 1 size
  REAL (KIND=dp), ALLOCATABLE :: nufluxp(:,:)  ! flux*frequency for 1 type and all sizes

  ! for charge and mix
  REAL (KIND=dp), ALLOCATABLE :: fz0(:)
  REAL (KIND=dp)              :: f0(nsize_max), f1(nsize_max)        ! fractions of neutrals and charged 

  ! for spinning dust SED
  REAL (KIND=dp), ALLOCATABLE :: nuf(:)
  REAL (KIND=dp), ALLOCATABLE :: spnuflux_j(:)  ! flux*frequency for 1 type and 1 size
  REAL (KIND=dp), ALLOCATABLE :: spnuflux(:,:)  ! flux*frequency for 1 type and all sizes

  ! for integration overs sizes
  INTEGER                     :: ns_i
  REAL (KIND=dp)              :: xx(nsize_max), yy(nsize_max)       ! Integrand

  !-------------------------------------------------------------------------

  ! allocate IR SED arrays 
  ALLOCATE (nuflux_j(n_qabs))
  ALLOCATE (nuflux(nsize_max,n_qabs))

  ! allocate polarized SED arrays
  IF (n_pol > 0 .and. n_lin > 0) THEN 
     ALLOCATE (nufluxp_j(n_qabs))
     ALLOCATE (nufluxp(nsize_max,n_qabs))
  ENDIF

  ! allocate charge distribution arrays
  IF (n_chrg > 0) THEN 
     ALLOCATE (zmean(nsize_max),sd_z(nsize_max))
     ALLOCATE (jpe(nsize_max),hpe(nsize_max),hspe(nsize_max),jgas(nsize_max),cgas(nsize_max))
  ENDIF
  f0 = 1.0_dp
  f1 = 1.0_dp

  ! allocate spinning arrays
  IF (n_spin > 0) THEN 
     ALLOCATE (spnuflux_j(n_qabs))
     ALLOCATE (spnuflux(nsize_max, n_qabs))
  ENDIF

  IF (n_ftemp > 0) THEN
     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = trimcat(filename_tmp,ftemp%nom)
     OPEN (UNIT=ftemp%unit, FILE = filename_tmp, STATUS = 'unknown')
     WRITE (UNIT=ftemp%unit, FMT='(A40)')'# DUSTEM: grain temperature distribution'
     WRITE (UNIT=ftemp%unit, FMT='(A1)') '#'
     WRITE (UNIT=ftemp%unit, FMT='(A119)') &
     & '# TYPE                   size(cm)          Teq(K)              Tmoy(K)            Tmax(K)        Trot(K)          ndist'
     WRITE (UNIT=ftemp%unit, FMT='(A1)') '#'
     WRITE (UNIT=ftemp%unit, FMT='(A87)') &
     & '#       T(K)             dP/dlnT           C(T)(erg/K)         U(erg)             dP/dU'
     WRITE (UNIT=ftemp%unit, FMT='(A1)') '#'
  ENDIF

  IF ((n_zdist > 0) .AND. (n_chrg > 0)) THEN
     filename_tmp = TRIMCAT(data_path,dir_res)
     filename_tmp = trimcat(filename_tmp,fzdist%nom)
     OPEN (UNIT=fzdist%unit, FILE = filename_tmp, STATUS = 'unknown')
     WRITE (UNIT=fzdist%unit, FMT='(A35)')'# DUSTEM: grain charge distribution'
     WRITE (UNIT=fzdist%unit, FMT='(A1)') '#'
     WRITE (UNIT=fzdist%unit, FMT='(A76,A103)') &
     & '# TYPE                     size(cm)        Zmin           Zmax           Zeq', &
     & '            Zmean          sd_Z         Jpe(s-1)      Hpe(erg/s)     Jgas(s-1)      Cgas(erg/s)     nzb'
     WRITE (UNIT=fzdist%unit, FMT='(A1)') '#'
     WRITE (UNIT=fzdist%unit, FMT='(A31)') &
     & '#         Z                f(Z)'
     WRITE (UNIT=fzdist%unit, FMT='(A1)') '#'
  ENDIF

  IF (n_res_a > 0) THEN
     filename_tmp = TRIMCAT(data_path, dir_res)
     filename_tmp = trimcat(filename_tmp,fsed_a%nom)
     OPEN (UNIT = fsed_a%unit, FILE = filename_tmp, STATUS = 'unknown')
     WRITE (UNIT = fsed_a%unit, FMT='(A58)') '# DUSTEM: grain SED per size 4*pi*nu*I_nu/Ng (erg/s/grain)'
     WRITE (UNIT = fsed_a%unit, FMT='(A1)')  '#'
     WRITE (UNIT = fsed_a%unit, FMT='(A33)') '# TYPE  nr of sizes   nr of waves'
     WRITE (UNIT = fsed_a%unit, FMT='(A28)') '# size(1)...size(nsize) (cm)'
     WRITE (UNIT = fsed_a%unit, FMT='(A38)') '# lambda(microns) SED(1)....SED(nsize)'
     WRITE (UNIT = fsed_a%unit, FMT='(A1)')  '#'

     IF (n_pol > 0 .and. n_lin > 0) THEN 
        filename_tmp = TRIMCAT(data_path, dir_res)
        filename_tmp = trimcat(filename_tmp,fpsed_a%nom)
        OPEN (UNIT = fpsed_a%unit, FILE = filename_tmp, STATUS = 'unknown')
        WRITE (UNIT = fpsed_a%unit, FMT='(A68)') '# DUSTEM: grain polarized SED per size 4*pi*nu*I_nu/Ng (erg/s/grain)'
        WRITE (UNIT = fpsed_a%unit, FMT='(A1)')  '#'
        WRITE (UNIT = fpsed_a%unit, FMT='(A33)') '# TYPE  nr of sizes   nr of waves'
        WRITE (UNIT = fpsed_a%unit, FMT='(A28)') '# size(1)...size(nsize) (cm)'
        WRITE (UNIT = fpsed_a%unit, FMT='(A38)') '# lambda(microns) SED(1)....SED(nsize)'
        WRITE (UNIT = fpsed_a%unit, FMT='(A1)')  '#'
     ENDIF
  ENDIF

  IF (n_spin > 0) THEN
     filename_tmp = TRIMCAT(data_path, dir_res)
     filename_tmp = trimcat(filename_tmp,fspd_a%nom)
     OPEN (UNIT = fspd_a%unit, FILE = filename_tmp, STATUS = 'unknown')
     WRITE (UNIT = fspd_a%unit, FMT='(A67)') '# DUSTEM: grain spinning SED per type 4*pi*nu*I_nu/Ng (erg/s/grain)'
     WRITE (UNIT = fspd_a%unit, FMT='(A1)')  '#'
     WRITE (UNIT = fspd_a%unit, FMT='(A33)') '# TYPE  nr of types   nr of waves'
     WRITE (UNIT = fspd_a%unit, FMT='(A28)') '# size(1)...size(nsize) (cm)'
     WRITE (UNIT = fspd_a%unit, FMT='(A38)') '# lambda(microns) SED(1)....SED(ntype)'
     WRITE (UNIT = fspd_a%unit, FMT='(A1)')  '#'
  ENDIF

  ! compute and write
  DO i=1,ntype                                     ! loop on grain types

     t_rot(:) = 0.0_dp
     ns_i     = nsize(i)
     n_beta   = 0
     n_dtls = 0
     IF (INDEX(t_opt(i),'BETA') > 0)   n_beta   = 1 
     IF (INDEX(t_opt(i),'DTLS') > 0)   n_dtls   = 1 

     DO j=1,nsize(i)                               ! loop on sizes

        CALL GET_TDIST( &
             & i,             &                    ! (I): index of grain TYPE
             & j,             &                    ! (I): index of grain size
             & size_ava(i,j), &                    ! (I): grain size
             & enerabs(j),    &                    ! (O): power absorbed by grain
             & uint,          &                    ! (O): internal energy grid
             & dpu,           &                    ! (O): dP/dU
             & temp,          &                    ! (O): grain temperature T
             & dpt,           &                    ! (O): dP/dln T
             & hcap,          &                    ! (O): grain heat capacity vs temp
             & tempmoy(j) ,   &                    ! (O): average temperature 
             & tempequi(j),   &                    ! (O): equilibrium temperature 
             & tempmax(j) )                        ! (O): max temperature 

        t2(:) = LOG(temp(:))

        CALL GET_ASED (i,j,size_ava(i,j),ndist,temp(1:ndist),t2(1:ndist),dpt(1:ndist),nuflux_j(:),enerem(j))
        nuflux(j,:) = nuflux_j(:)

        IF (n_pol > 0 .and. n_lin > 0 .and. INDEX(t_opt(i),'POL') > 0) THEN
           CALL GET_ASED_POL (i,j,size_ava(i,j),ndist,temp(1:ndist),t2(1:ndist),dpt(1:ndist),nufluxp_j(:))
           nufluxp(j,:) = nufluxp_j(:)
        ENDIF

        IF ((INDEX(t_opt(i),'CHRG') > 0 ) .OR. (INDEX(t_opt(i),'SPIN') > 0 )) THEN
          CALL ZDIST(i, j, size_ava(i,j), zmean(j), sd_z(j), jpe(j), hpe(j), hspe(j), jgas(j), cgas(j))
          IF (n_zdist > 0) THEN 
             WRITE (UNIT = fzdist%unit, FMT='(/,A2,A20,1x,10(1PE12.4,3X),i4)') '# ',gtype(i),size_ava(i,j),zmin,zmax,zeq,zmean(j), &
                  & sd_z(j),jpe(j),hpe(j),jgas(j),cgas(j),nzb
             DO kt = 1, nzb
                WRITE (UNIT = fzdist%unit, FMT='(1P,2(E18.10E3,1X))') zb(kt),fz(kt)
             ENDDO
          ENDIF
          IF (INDEX(t_opt(i),'ZM') > 0 ) THEN                     ! zm: 1st pass
             ALLOCATE (fz0(nzb))
             fz0 = fz
             WHERE(zb/=0.0_dp) fz0=0.0_dp  
             f0(j) = SUM(fz0)
             f1(j) = 1.0_dp - f0(j)
             f_mix(i,j) = f0(j)
          ENDIF
        ENDIF

        IF ((INDEX(t_opt(i),'ZM') > 0 ) .AND. (n_zm == 1)) THEN    ! zm: 2nd pass
           f_mix(i,j) = f1(j)
        ENDIF

        IF (INDEX(t_opt(i),'SPIN') > 0 ) THEN
           ALLOCATE (nuf(n_qabs))
           nuf(:) = nuflux_j(:) * enerabs(j)/enerem(j) / 4.0_dp/xpi
           DO kt=1,n_qabs 
              nuf(kt) = nuflux_j(n_qabs-kt+1) * enerabs(j)/enerem(j) / 4.0_dp/xpi / freq_qabs(kt)
           ENDDO
           CALL SPIN(i, j, size_ava(i,j), tempequi(j), nuf, spnuflux_j, t_rot(j))
           DEALLOCATE (nuf)
           spnuflux(j,:) = spnuflux_j(:)
           nuflux(j,:) = nuflux(j,:) + spnuflux(j,:)
        ENDIF

        IF (n_ftemp > 0) THEN
           WRITE (UNIT = ftemp%unit, FMT='(/,A2,A20,1x,5(1PE12.4,6X),i4)') '# ',gtype(i),size_ava(i,j),tempequi(j), &
                & tempmoy(j),tempmax(j),t_rot(j),ndist
           DO kt = 1, ndist
              WRITE (UNIT = ftemp%unit, FMT='(1P,5(E18.10E3,1X))') temp(kt), dpt(kt), hcap(kt), uint(kt), dpu(kt)
           ENDDO
        ENDIF

        IF ((INDEX(t_opt(i),'CHRG')>0) .OR. (INDEX(t_opt(i),'SPIN')>0)) THEN 
           DEALLOCATE(zb,fz)
           IF (INDEX(t_opt(i),'ZM')>0) DEALLOCATE(fz0)
        ENDIF

     ENDDO    ! end of size loop (j)

     xx(1:ns_i) = si_ava_l(i,1:ns_i)
     DO kt=1,n_qabs
       yy(1:ns_i) = ( nuflux(1:ns_i,kt)/ size_ava(i,1:ns_i)**3 ) &
                  * (enerabs(1:ns_i) / enerem(1:ns_i)) &
                  * f_mix(i,1:ns_i) * ava(i,1:ns_i)
       nuinuem(i,kt) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
     ENDDO

     ! get SED in erg/s/H (4*pi*nu*inu/NH)
     nuinuem(i,:) = nuinuem(i,:) * mprop(i) * xmp / masstot(i)
     nuinuemtot = nuinuemtot + nuinuem(i,:)

     IF (n_res_a > 0) THEN
        WRITE (UNIT=fsed_a%unit, FMT='(A2,a20,1x,i3,1X,I4)') '# ', gtype(i), nsize(i), n_qabs
        WRITE (UNIT=fsed_a%unit, FMT='(100(1PE12.4,1X))')  (size_ava(i,j), j=1,nsize(i))
        DO kt=1,n_qabs
           WRITE (UNIT=fsed_a%unit, FMT='(101(1PE14.6E3,1X))') lamb_qabs(kt)*1.e4_dp, &
                (nuflux(j,kt)*enerabs(j)/enerem(j), j=1,nsize(i))
        ENDDO
     ENDIF

     IF ( ((INDEX(t_opt(i),'CHRG') > 0) .AND. (INDEX(t_opt(i),'ZM') > 0 ) .AND. (n_zm < 1)) &
          & .OR. ((INDEX(t_opt(i),'SPIN') > 0) .AND. (INDEX(t_opt(i),'ZM') > 0 ) .AND. (n_zm < 1)) ) THEN
        n_zm = n_zm + 1
     ELSE 
        n_zm = 0
     ENDIF

     IF (INDEX(t_opt(i),'SPIN') > 0 ) THEN
        DO kt=1,n_qabs
           yy(1:ns_i) = ( spnuflux(1:ns_i,kt)/ size_ava(i,1:ns_i)**3 ) &
                * (enerabs(1:ns_i) / enerem(1:ns_i)) &
                * f_mix(i,1:ns_i) * ava(i,1:ns_i)
           spnuinuem(i,kt) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
        ENDDO
        spnuinuem(i,:) = spnuinuem(i,:) * mprop(i) * xmp / masstot(i)
        spnuinuemtot = spnuinuemtot + spnuinuem(i,:)
        IF (n_res_a > 0) THEN
           WRITE(UNIT = fspd_a%unit, FMT='(A2,a20,1x,i3,1X,I4)') '# ', gtype(i), nsize(i), n_qabs
           WRITE(UNIT = fspd_a%unit, FMT='(100(1PE12.4,1X))')  (size_ava(i,j), j=1,nsize(i))
           DO kt=1,n_qabs
              WRITE(UNIT = fspd_a%unit, FMT='(101(1PE14.6E3,1X))') lamb_qabs(kt)*1.e4_dp, (spnuflux(j,kt), j=1,nsize(i))
           ENDDO
        ENDIF        
     ENDIF

     IF (n_pol > 0 .and. n_lin > 0 .and. INDEX(t_opt(i),'POL') > 0) THEN 
        DO kt=1,n_qabs
           yy(1:ns_i) = ( nufluxp(1:ns_i,kt)/ size_ava(i,1:ns_i)**3 ) &
                * (enerabs(1:ns_i) / enerem(1:ns_i)) &
                * f_mix(i,1:ns_i) * ava(i,1:ns_i)
           nuinuemp(i,kt) = XINTEG2 (1, ns_i, ns_i, xx(1:ns_i), yy(1:ns_i))
        ENDDO
        nuinuemp(i,:) = nuinuemp(i,:) * mprop(i) * xmp / masstot(i)
        nuinuemptot = nuinuemptot + nuinuemp(i,:)
        IF (n_res_a > 0) THEN
           WRITE(UNIT = fpsed_a%unit, FMT='(A2,a20,1x,i3,1X,I4)') '# ', gtype(i), nsize(i), n_qabs
           WRITE(UNIT = fpsed_a%unit, FMT='(100(1PE12.4,1X))')  (size_ava(i,j), j=1,nsize(i))
           DO kt=1,n_qabs
              WRITE(UNIT = fpsed_a%unit, FMT='(101(1PE14.6E3,1X))') lamb_qabs(kt)*1.e4_dp, (nufluxp(j,kt), j=1,nsize(i))
           ENDDO
        ENDIF
     ENDIF

  ENDDO  ! end type loop (i)

  ! deallocate arrays
  DEALLOCATE(nuflux_j)
  DEALLOCATE(nuflux)
  IF (n_pol > 0 .and. n_lin > 0) THEN 
     DEALLOCATE (nufluxp_j)
     DEALLOCATE (nufluxp)
  ENDIF
  IF (n_spin > 0) THEN 
     DEALLOCATE(spnuflux_j)
     DEALLOCATE(spnuflux)
  ENDIF

  ! close files
  IF ((n_zdist > 0) .AND. (n_chrg > 0)) CLOSE (UNIT = fzdist%unit)
  IF (n_ftemp > 0) CLOSE (UNIT = ftemp%unit)
  IF (n_res_a > 0) THEN 
     CLOSE (UNIT = fsed_a%unit)
     IF (n_pol > 0 .and. n_lin > 0) CLOSE (UNIT = fpsed_a%unit)
  ENDIF

END SUBROUTINE COMPUTE

! -------------------------------------------------------------------------------

SUBROUTINE GET_ASED (nt, ns, a, ntsed, t, ti, p, xnufnu_los, energiemise)

! computes the SED (erg/s) emitted by a grain of given size and with  
! temperature distribution p (usually dp/dlnT). 
! Temperature grid for integration is ti (usually LOG(T))
! uses lambda*B_lambda(T)

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT
  USE MDTLS

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt               ! index of grain type
  INTEGER, INTENT (IN)         :: ns               ! index of grain size
  INTEGER, INTENT (IN)         :: ntsed            ! nr of T-values
  REAL (KIND=dp), INTENT (IN)  :: a                ! grain size
  REAL (KIND=dp), INTENT (IN)  :: t(ntsed)
  REAL (KIND=dp), INTENT (IN)  :: ti(ntsed)
  REAL (KIND=dp), INTENT (IN)  :: p(ntsed)
  ! Anisotropic emission : energy emitted toward the observer
  REAL (KIND=dp), INTENT (OUT) :: xnufnu_los(n_qabs)
  ! Isotropic emission : energy emitted in all directions Only useful to calculate energiemise
  REAL (KIND=dp)               :: xnufnu(n_qabs)
  REAL (KIND=dp), INTENT (OUT) :: energiemise

  INTEGER        :: l, ktemp
  REAL (KIND=dp) :: al1
  REAL (KIND=dp) :: fnusom
  REAL (KIND=dp) :: fnut(ntsed), xx(ntsed)
  REAL (KIND=dp) :: yy(n_qabs)
  REAL (KIND=dp) :: coef
  REAL (KIND=dp) :: sig(ntsed)

  coef = xpi * a**2 * cte1

  ! main loop on wavelengths
  energiemise = 0.0_dp
  IF (n_beta == 0 .AND. n_dtls == 0) THEN
     DO l=1,n_qabs
        al1 = lamb_qabs(l)
        ! temperature loop
        xx(:) = hcsurk / (al1 * t(:))
        DO ktemp=1,ntsed
           fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp))
        ENDDO
        fnusom = XINTEG2 (1, ntsed, ntsed, ti, fnut)
        ! qauv is backward...
        xnufnu(l) = coef * fnusom * qi_abs(nt,ns,n_qabs-l+1) / al1**4    

        IF (n_anis .eq. 0 .or. INDEX(t_opt(nt),'POL') == 0) THEN
            ! We ignore alignment for total emission
            xnufnu_los(l) = xnufnu(l)
        ELSE 
            ! We take alignment into account for total emission
            ! Contribution of not aligned grains (1-f) + that of aligned grains (f)
           xnufnu_los(l) = coef * fnusom / al1**4 &
                * ( (1. - f_pol(nt,ns)) * qi_abs(nt,ns,n_qabs-l+1) &
                + f_pol(nt,ns) * (q1i_abs(nt,ns,n_qabs-l+1)+q2i_abs(nt,ns,n_qabs-l+1))/2 )
        ENDIF
     ENDDO
  ELSE IF (n_beta == 1 .AND. n_dtls == 0) THEN    ! apply BETA(T)
     DO l=1,n_qabs
        al1 = lamb_qabs(l)
        ! temperature loop
        xx(:) = hcsurk / (al1 * t(:))
        DO ktemp=1,ntsed
           fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp)) * &
           ! faster
           & EXP( DBETA(t(ktemp),nt)*f_beta(nt,l)*LOG(ltresh(nt)/al1) )
        ENDDO
        fnusom = XINTEG2 (1, ntsed, ntsed, ti, fnut)
        ! qauv is backward...
        xnufnu(l) = coef * fnusom * qi_abs(nt,ns,n_qabs-l+1) / al1**4 
     ENDDO
     xnufnu_los = xnufnu
  ELSE IF (n_beta == 0 .AND. n_dtls == 1) THEN  ! apply DCD-TLS
     DO l=1,n_qabs
        al1 = lamb_qabs(l)
        ! temperature loop
        xx(:) = hcsurk / (al1 * t(:))
        IF (al1*1e4_dp < ldtresh(nt)) THEN 
           DO ktemp=1,ntsed
              fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp)) * qi_abs(nt,ns,n_qabs-l+1)
           ENDDO
        ELSE 
           CALL DTLS (nt, ntsed, ns, n_qabs-l+1, a, t, sig)
           DO ktemp=1,ntsed
              fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp)) * sig(ktemp)
           ENDDO
        ENDIF
        fnusom = XINTEG2 (1, ntsed, ntsed, ti, fnut)
        xnufnu(l) = coef * fnusom / al1**4 
     ENDDO
     xnufnu_los = xnufnu
  ENDIF

  yy(:) = xnufnu(:) / lamb_qabs(:)
  l = 1
  energiemise = XINTEG2 (l, n_qabs, n_qabs, lamb_qabs, yy)

END SUBROUTINE GET_ASED

! -------------------------------------------------------------------------------

SUBROUTINE GET_ASED_POL (nt, ns, a, ntsed, t, ti, p, xnufnu)

! computes the SED (erg/s) emitted by a grain of given size and with  
! temperature distribution p (usually dp/dlnT). 
! Temperature grid for integration is ti (usually LOG(T))
! uses lambda*B_lambda(T)

  USE CONSTANTS
  USE UTILITY
  USE MGET_QEXT

  IMPLICIT NONE

  INTEGER, INTENT (IN)         :: nt               ! index of grain type
  INTEGER, INTENT (IN)         :: ns               ! index of grain size
  INTEGER, INTENT (IN)         :: ntsed            ! nr of T-values
  REAL (KIND=dp), INTENT (IN)  :: a                ! grain size
  REAL (KIND=dp), INTENT (IN)  :: t(ntsed)
  REAL (KIND=dp), INTENT (IN)  :: ti(ntsed)
  REAL (KIND=dp), INTENT (IN)  :: p(ntsed)
  REAL (KIND=dp), INTENT (OUT) :: xnufnu(n_qabs)

  INTEGER        :: l, ktemp
  REAL (KIND=dp) :: al1
  REAL (KIND=dp) :: fnusom
  REAL (KIND=dp) :: fnut(ntsed), xx(ntsed)
  REAL (KIND=dp) :: coef

  coef = xpi * a**2 * cte1

  ! main loop on wavelengths
  IF (n_beta == 0) THEN
     DO l=1,n_qabs
        al1 = lamb_qabs(l)
        ! temperature loop
        xx(:) = hcsurk / (al1 * t(:))
        DO ktemp=1,ntsed
           fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp))
        ENDDO
        fnusom = XINTEG2 (1, ntsed, ntsed, ti, fnut)
        ! qauv is backward...
        xnufnu(l) = coef * fnusom * (q2i_abs(nt,ns,n_qabs-l+1)-q1i_abs(nt,ns,n_qabs-l+1))/2 / &
             & al1**4    
     ENDDO
  ELSE                  ! apply BETA(T)
     DO l=1,n_qabs
        al1 = lamb_qabs(l)
        ! temperature loop
        xx(:) = hcsurk / (al1 * t(:))
        DO ktemp=1,ntsed
           fnut(ktemp) = p(ktemp) * F_BB (xx(ktemp)) * &
           ! faster
           & EXP( DBETA(t(ktemp),nt)*f_beta(nt,l)*LOG(ltresh(nt)/al1) )
        ENDDO
        fnusom = XINTEG2 (1, ntsed, ntsed, ti, fnut)
        ! qauv is backward...
        xnufnu(l) = coef * fnusom * (q2i_abs(nt,ns,n_qabs-l+1)-q1i_abs(nt,ns,n_qabs-l+1))/2 / &
             & al1**4 
     ENDDO
  ENDIF

  xnufnu(:) = xnufnu(:) * f_pol(nt,ns) 

END SUBROUTINE GET_ASED_POL

END MODULE MCOMPUTE
