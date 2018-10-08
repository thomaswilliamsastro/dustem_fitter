
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

  ! WRITE (*,*) '--> Reading...'
  CALL READ_DATA

  ! WRITE (*,*) '--> WORKING...'
  CALL COMPUTE

  ! WRITE (*,*) '--> writing...'
  CALL WRITE_DATA

  CALL CPU_TIME (t1)

  IF (n_quiet == 0) THEN
     WRITE (*,*) ''
     WRITE (*,FMT='(A20,1PE8.2,A4)')'Exit OK - CPU Time ', t1,' sec'
     WRITE (*,*) '==========================================='
  ENDIF

END PROGRAM DUSTEM
