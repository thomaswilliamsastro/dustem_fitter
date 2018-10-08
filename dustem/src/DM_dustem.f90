
! DUSTEM - Dust Emissivity
! See: Desert et al., 1986, A&A 160, 295 for temperature distribution
! 2008-09:  unique spectral range and new dust properties: M Compiegne
! Aut 2009: code fully rewritten: L Verstraete & J Le Bourlot
! 2009-10:  adding size dist parameters, beta function & output: LV
! Spr 2010: adding polarization: LV, L Masson & V Guillet
!           adding DCD-TLS and spinning: N Ysard & LV

PROGRAM DUSTEM

  USE CONSTANTS
  USE UTILITY
  USE MCOMPUTE
  USE IN_OUT

  IMPLICIT NONE

  REAL (KIND=dp) t1
  TYPE (FICH)                 :: input     != FICH ('GRAIN.DAT', 11)
  TYPE(FICH)                  :: output   != FICH ('SED_test.RES', 18)

  integer::cptArg,n
  character(len=99)::name

  n = command_argument_count()
  !loop across options
  do cptArg=1,n
    call get_command_argument(cptArg,name)
      if (cptArg == 1) then
        input = FICH(name,11)
      else if (cptArg == 2) then
        output = FICH(name,18)
      end if
  end do


  ! WRITE (*,*) '--> Reading...'
  CALL READ_DATA(input)

  ! WRITE (*,*) '--> WORKING...'
  CALL COMPUTE

  ! WRITE (*,*) '--> writing...'
  CALL WRITE_DATA(output)

  CALL CPU_TIME (t1)

  IF (n_quiet == 0) THEN
     WRITE (*,*) ''
     WRITE (*,FMT='(A20,1PE8.2,A4)')'Exit OK - CPU Time ', t1,' sec'
     WRITE (*,*) '==========================================='
  ENDIF

END PROGRAM DUSTEM
