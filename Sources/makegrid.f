!******************************************************************************
!  File makegrid.f
!*******************************************************************************
!**********************************************************************
!**                                                                  **
!**   MAIN PROGRAM:  Make Grid                                       **
!**                                                                  **
!**                                                                  **
!**     PROGRAM DESCRIPTION:                                         **
!**        MAKEGRID reads in a coils-dot file, and generates         **
!**        an MGRID file.                                            **
!**                                                                  **
!**     Author: Steve Hirshman                                       **
!**             James Hanson                                         **
!**             Jonathan Hebert                                      **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**********************************************************************
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!  This file uses the following modules:
!    stel_kinds
!       located in LIBSTELL/Sources/Modules
!
!    stel_constants
!       located in LIBSTELL/Sources/Modules
!
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  Executing 'xgrid -h' will printout a help message
!
!    INPUT FILES
!
!
!   OUTPUT FILES
!      
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!
!*******************************************************************************
!  CODE MAKEGRID
!    
! SECTION I.	Main Program
!
! SECTION II.	Initialization Subroutines
!    interactive_input
!    namelist_input
!
! SECTION III.	TASK SUBROUTINES
!    task_mgrid
!    task_mgrid_rs
!    task_circ_tor_grid
!
! SECTION IV. AUXILIARY SUBROUTINES
!    coil_group_report
!    fft_local
!
! SECTION V.	SUBROUTINES FOR TESTNG
!
! SECTION VI.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
! SECTION I. MAIN PROGRAM
!*******************************************************************************

      PROGRAM makegrid
!
!     THIS CODE (MAKEGRID) GENERATES B-FIELD COMPONENTS ON R,Z, PHI GRID
!     AND COMPUTES THE POLOIDAL FLUX AT NOBS OBSERVATION POINTS
!
!     NOTE TO USER: EXPERIENCE SHOWS THAT A GRID SPACING OF
!     DEL R = DEL Z <= 1/50 WORKS WELL.  HERE, DEL R = rmax-rmin.
!     GRID SPACING MUCH LARGER THAN THIS MAY ADVERSELY EFFECT CONVERGENCE
!     OF VACUUM CODE.
!
!     BOX DIMENSIONS: rmin <= R <= rmax,  zmin <= Z <= zmax
!                     KP = NO. TOROIDAL PLANES/FIELD PERIOD (MUST EQUAL VMEC VALUE)
!-------------------------------------------------------------------------------
!
!     grid dimensions:
!          ir  = no. radial (r) points in box
!          jz  = no. z points in box
!
!          suggest choosing hr == (rmax-rmin)/(ir-1) equal to
!                           hz == (zmax-zmin)/(jz-1)
!
!
!     THIS CAN BE RUN EITHER FROM THE COMMAND LINE (INTERACTIVELY) OR FROM A COMMAND FILE:
!
!     xgrid < filename.cmd
!
!     WHERE the command file (filename.cmd) has entries:
!     coils file extension
!     stell_sym (T/F)
!     rmin value
!     rmax value
!     zmin value
!     zmax value
!

!
!     SET UP GRID DIMENSIONS. USE EITHER INTERACTIVE (OR DRIVER FILE)
!     OR COMMAND LINE ARGUMENT LIST
!
!-------------------------------------------------------------------------------
!  CHANGE LOG
!-------------------------------------------------------------------------------
!     1-18-2010  JDHe 
!     Added task structure to makegrid (does not change functionality of old
!        code)
!     Added namelist input support (usable with all tasks)
!     Added help message
!     Added task to calculate vector potential on inside of vacuum vessel
!        for use with NIMROD (task CIRC_TOR_GRID)
!     Previous functionality available as called before (no command line
!        inputs prompts for mgrid data, 10 command line inputs runs mgrid with
!        those inputs) as well as through the namelist (task MGRID)
!     Added 2 arguments to command line input to reflect interactive input
!---------------------------------------------------------------------------
!  Use statements are followed by variables and subroutines used in this file
!  (subroutines are followed by parentheses)
!---------------------------------------------------------------------------

      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &                       rmin, rmax, zmin, zmax, kp, ir, jz

      USE makegrid_global, only: task

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numargs, i
      CHARACTER(LEN=100) :: arg1
      CHARACTER(LEN=20) :: task_use
      
      CHARACTER(LEN=57), DIMENSION(23), PARAMETER  :: help_message=
     & (/' Program makegrid command line help                      ',
     &   '   -h    Display this message.                           ',
     &   ' Number of Command Line Arguments:                       ',
     &   '  1 (not -h):                                            ',
     &   '   Run makegrid using Namelist input file supplied in    ',
     &   '   argument                                              ',
     &   '  10:                                                    ',
     &   '   Use command line arguments as input to makegrid in    ',
     &   '   order (only for task MGRID):                          ',
     &   '     mgrid_ext                                           ',
     &   '     mgrid_mode                                          ',
     &   '     lstell_sym                                          ',
     &   '     rmin                                                ',
     &   '     rmax                                                ',
     &   '     zmin                                                ',
     &   '     zmax                                                ',
     &   '     kp                                                  ',
     &   '     ir  (optional,default=121)                          ',
     &   '     jz  (optional,default=121)                          ',
     &   '  0:                                                     ',
     &   '   Get command line input interactively (only for task   ',
     &   '   MGRID)                                                ',
     &   ' END HELP MESSAGE                                        '/)

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      CALL getcarg(1, arg1, numargs)

      SELECT CASE(numargs)

         CASE (1)
!  INTERACTIVE (USE REDIRECTED DRIVER FILE ALSO, XGRID < DRIVER)
            IF (arg1 .EQ. '-h' .OR. arg1 .EQ. '--h' .OR. arg1 .EQ. '-H'
     &          .OR. arg1 .EQ. '--H') THEN
               DO i = 1, SIZE(help_message)
                  WRITE(6,*) help_message(i)
               END DO
               STOP
            ELSE
               CALL namelist_input(arg1)
            END IF

         CASE(0)
!  Only used for task 'MGRID'
            CALL interactive_input() !sets task='MGRID'

         CASE(8:10)
!  Only used for task 'MGRID'
            arg1 = ADJUSTL(arg1)
            mgrid_ext = TRIM(arg1)
            CALL getcarg(2, arg1, numargs)
            IF (arg1(1:1) == 'R' .or. arg1(1:1) == 'r') THEN
               mgrid_mode = 'R'
            END IF
            CALL getcarg(3, arg1, numargs)
            lstell_sym = arg1(1:1) == 'Y' .or. arg1(1:1) == 'y' .or.
     &                   arg1(1:1) == 'T' .or. arg1(1:1) == 't'
            CALL getcarg(4, arg1, numargs)
            READ (arg1, *) rmin
            CALL getcarg(5, arg1, numargs)
            READ (arg1, *) rmax
            CALL getcarg(6, arg1, numargs)
            READ (arg1, *) zmin
            CALL getcarg(7, arg1, numargs)
            READ (arg1, *) zmax
            CALL getcarg(8, arg1, numargs)
            READ (arg1, *) kp
            IF (numargs .GE. 9) THEN
               CALL getcarg(9, arg1, numargs)
               READ (arg1, *) ir
            END IF
            IF (numargs .EQ. 10) THEN
               CALL getcarg(10, arg1, numargs)
               READ (arg1, *) jz
            END IF
            task='MGRID'

         CASE DEFAULT
            STOP 'Unknown number of arguments, type "xgrid -h" ' //            &
     &           'for help.'
      
      END SELECT

!-----------------------------------------------------------------------
!  TASKS:
!  'makegrid' can now run two separate tasks:
!    'B_GRID':         find the magnetic field on a grid of points inside the
!                      plasma volume (for VMEC, V3FIT, etc.)
!    'CIRC_TOR_GRID':  find the vector potential on the outer-most surface
!                      (for NIMROD)
!-----------------------------------------------------------------------

!  Convert to lower case, avoid simple misinterpretations.
      task_use = task
      CALL tolower(task_use)
      
      SELECT CASE(TRIM(ADJUSTL(task_use)))
      
         CASE('mgrid')
            CALL task_mgrid()

         CASE('mgrid_rs')
            CALL task_mgrid_rs()

         CASE('circ_tor_mgrid')
            CALL task_circ_tor_grid()
          
         CASE DEFAULT
            WRITE(*,*) 'Unknown task ', TRIM(ADJUSTL(task_use))
      
      END SELECT

      END PROGRAM makegrid

!*******************************************************************************
! SECTION II.	Initialization Subroutines
!*******************************************************************************
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE namelist_input(arg1)
!  namelist_input
!    Read a namelist input file for either task.

      USE stel_constants
      
      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &                       rmin, rmax, zmin, zmax, kp, ir, jz,               &
     &                       use_eddy

      USE makegrid_global, only: task, rmajor, aminor, nphi, ntheta,           &
     &                           extcur_mgrid, cg_shift_1, cg_shift_2,         &
     &                           cg_rot_xcent, cg_rot_theta, cg_rot_phi,       &
     &                           cg_rot_angle, l_rot_coil_center
     
      USE sym_check, ONLY: sym_ir, sym_jz, sym_kp, sym_perform_tests

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      CHARACTER(LEN=100), INTENT(IN) :: arg1
      INTEGER                        :: ferror
      INTEGER, PARAMETER             :: iou_nli=12

!-----------------------------------------------
!   **  Namelist Variables  **
!
!  Used in all tasks:
!    task        String to select which mode 'makegrid' to run.  See comment
!                     in program 'makegrid' above.               (makegrid_global)
!    mgrid_ext   String to concatenate to all file names.  Must be the extension
!                     of the 'coils-dot' file to be used.        (write_mgrid)
!
!  Used in task 'mgrid' and 'mgrid_rs'
!    mgrid_mode   Character to choose wheter to run in raw or scaled mode   (write_mgrid)
!    lstell_sym   Logical of whether or not to assume stellarator symmetry  (write_mgrid)
!    rmin         Minimum radial position on grid                           (write_mgrid)
!    rmax         Maximum radial position on grid                           (write_mgrid)
!    zmin         Minimum vertical position on grid                         (write_mgrid)
!    zmax         Maximum vertical position on grid                         (write_mgrid)
!    kp           Number of toroidal planes per period                      (write_mgrid)
!    ir           Number of radial mesh points                              (write_mgrid)
!    jz           Number of vertical mesh points                            (write_mgrid)
!    use_eddy     Add an aditional blank coil group to store vacuum eddy
!                 current contribution. Eddy currents are not coputed by
!                 makegrid
!
!  Used in task 'mgrid_rs'. 
!        (All variables are declared in module makegrid_global)
!        (  : indicates dimension nextcur_dim - parameter in makegrid_global)
!    cg_shift_1(:,3)        Vector to shift all the coils. (Before rotation)
!    cg_rot_theta(:)        Spherical polar angle to specify axis of rotation
!    cg_rot_phi(:)          Spherical azimuthal angle to specify axis of rotation
!    cg_rot_angle(:)        Angle to rotate about axis of rotation. 
!                             NB - LEFT HAND convention. Put left thumb along 
!                             axis of rotation, fingers indicate direction of 
!                             positive rotation.
!    cg_rot_xcent(:,3)      Position of center of rotation
!    l_rot_coil_center(:)   Logical. True - use current-averaged center of
!                             coil-group for center of rotation
!                             False - use position specified in cg_rot_xcent
!                             for center of rotation
!    cg_shift_2(:,3)        Vector to shift all the coils. (After rotation)    
!
!  Used in task 'circ_tor_grid':
!    rmajor          Major radius of circular torus            (makegrid_global)
!    aminor          Minor radius of circular torus            (makegrid_global)
!    nphi            Number of toroidal mesh points            (makegrid_global)
!    ntheta          Number of poloidal mesh points            (makegrid_global)
!    extcur_mgrid    Current in each coil of 'coils-dot' file  (makegrid_global)
!-----------------------------------------------

      NAMELIST /mgrid_nli/ task,                                               &
     &    mgrid_ext, mgrid_mode, lstell_sym,                                   &
     &    rmin, rmax, zmin, zmax, kp, ir, jz,                                  &
     &    rmajor, aminor, nphi, ntheta, extcur_mgrid,                          &
     &    cg_shift_1, cg_shift_2, cg_rot_xcent, cg_rot_theta,                  &
     &    cg_rot_phi, cg_rot_angle, l_rot_coil_center,                         &
     &    sym_ir, sym_jz, sym_kp, sym_perform_tests, use_eddy
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  Initialize namelist

      task = 'MGRID'
      mgrid_ext = 'dummy.cth.m.d12c.f5ss'
      mgrid_mode = 'R'
      lstell_sym = .TRUE.
      rmin = 0.45
      rmax = 1.05
      zmin = -0.3
      zmax = 0.3
      kp = 5
!  Commented out below, so that ir and jz have default values as specified in
!    module write_mgrid
!      ir=20          
!      jz=20
      rmajor = 1
      aminor = 0.5
      nphi = 3
      ntheta = 3
      extcur_mgrid = 1.
      sym_ir = 5
      sym_jz = 5
      sym_kp = 4
      sym_perform_tests = .FALSE.
      
      cg_shift_1 = zero
      cg_shift_2 = zero
      cg_rot_xcent = zero
      cg_rot_theta = zero
      cg_rot_phi = zero
      cg_rot_angle = zero
      l_rot_coil_center = .true.

      use_eddy = .false.

      WRITE(*,*) ' Running makegrid using NLI file ', arg1
      OPEN(iou_nli, FILE=TRIM(ADJUSTL(arg1)),STATUS='OLD',IOSTAT=ferror)
      IF(ferror .NE. 0) THEN
         WRITE(*,*) 'Could not open NLI file ', arg1
         WRITE(*,*) 'Open failed with error status ', ferror
      END IF
            
      READ(iou_nli, mgrid_nli)
      
      CLOSE(iou_nli)

      END SUBROUTINE namelist_input
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE interactive_input()
!  interactive_input
!    For task MGRID, input parameters when prompted. This is for
!    backwards compatibility with the earlier version of makegrid.

      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                & 
     &                       rmin, rmax, zmin, zmax, kp, ir, jz,               &
     &                       mgrid_file

      USE makegrid_global, only: task

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER             :: numargs
      INTEGER             :: i
      CHARACTER (LEN=100) :: arg1
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      WRITE (6, 220, advance='no')                                             &
     &   ' Enter extension of "coils" file     : '
      READ (*, *) mgrid_ext

      WRITE (6, '(a,/,a)', advance='no')                                       &
     &   ' Scale (S) bfield to unit current/turn OR',                          &
     &   ' use raw (R) currents from coils file: '
      READ (*, *) mgrid_file
      IF (mgrid_file(1:1) == 'R' .or. mgrid_file(1:1) == 'r') THEN
         mgrid_mode = 'R'
      END IF

      WRITE (6, 220, advance='no')                                             & 
     &   ' Assume stellarator symmetry (Y/N)?  : '
      READ (*, *) mgrid_file
      lstell_sym = mgrid_file(1:1) == 'Y' .or. mgrid_file(1:1) == 'y'

      WRITE (6, 220, advance='no')                                             & 
     &   ' Enter rmin (min radial grid dimension)  : '
      READ (*, *) rmin

      WRITE (6, 220, advance='no')                                             &  
     &   ' Enter rmax (max radial grid dimension)  : '
      READ (*, *) rmax

      WRITE (6, 220, advance='no')                                             & 
     &   ' Enter zmin (min vertical grid dimension): '
      READ (*, *) zmin

      WRITE (6, 220, advance='no')                                             &  
     &   ' Enter zmax (max vertical grid dimension): '
      READ (*, *) zmax

      WRITE (6, 220, advance='no')                                             & 
     &   ' Enter number of toroidal planes/period  : '
      READ (*, *) kp

      WRITE (6, 220, advance='no')                                             & 
     &   ' Enter number of r (radial) mesh points  : '
      READ (*, *) ir

      WRITE (6, 220, advance='no')                                             & 
     &   ' Enter number of z mesh points  : '
      READ (*, *) jz
      task = 'MGRID'
      
 220  FORMAT(a)
      
      END SUBROUTINE interactive_input
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!

!*******************************************************************************
! SECTION III.	TASK SUBROUTINES
!*******************************************************************************
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_mgrid()
!  task_mgrid
!    By this point, the parameters (mgrid_ext, mgrid_mode, lstell_sym, rmin,
!    rmax, zmin, zmax, kp, ir, jz) have been input.  This subroutine evaluates 
!    the magnetic field in the plasma volume.

      USE stel_kinds
      
      USE write_mgrid
       !  , only: mgrid_ext, mgrid_mode, lstell_sym,    
       !   rmin, rmax, zmin, zmax, kp, ir, jz,   
       !   br, bz, bp, kp2, jz2, kp_odd, jz_odd, coil_file, mgrid_file,     
       !   iextc, nfp, nextcur, extcur, fperiod, delr, delp, delz
       !  subroutines: write_mgrid_nc, cleanup_biotsavart
       !  access to bsc_b in module bsc.
     
      USE biotsavart !  variables nfp_bs, coil_group
                     !  subroutines parse_coils_file

      USE makegrid_global, only: task

      USE safe_open_mod !safe_open()
      
      USE sym_check, ONLY: init_symmetry, cleanup_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                                 :: istat
      CHARACTER (LEN=100)                     :: extcur_file
      REAL (rprec)                            :: time_start
      REAL (rprec)                            :: time_end1
      REAL (rprec)                            :: time_end2

      REAL (rprec)                            :: r_ave
      REAL (rprec)                            :: z_ave
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: phi_array
      INTEGER                                 :: kmax
      INTEGER                                 :: k
      INTEGER                                 :: icoll
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      
      
      WRITE(*,*) ' Running makegrid with the following parameters:'
      WRITE(*,*) '   task       = ', task
      WRITE(*,*) '   mgrid_ext  = ', mgrid_ext
      WRITE(*,*) '   mgrid_mode = ', mgrid_mode
      WRITE(*,*) '   lstell_sym = ', lstell_sym
      WRITE(*,*) '   rmin       = ', rmin
      WRITE(*,*) '   rmax       = ', rmax
      WRITE(*,*) '   zmin       = ', zmin
      WRITE(*,*) '   zmax       = ', zmax
      WRITE(*,*) '   kp         = ', kp
      WRITE(*,*) '   ir         = ', ir
      WRITE(*,*) '   jz         = ', jz
      WRITE(*,*) '   use_eddy   = ', use_eddy
      WRITE(*,*)

      IF (rmin .lt. 0.) STOP ' rmin must be > 0 in xgrid'
      IF (rmax .le. rmin) STOP ' rmax must be > rmin in xgrid'
      IF (zmax .le. zmin) STOP ' zmax must be > zmin in xgrid'
      IF (kp .le. 0) STOP 'kp must be > 0 in xgrid'

      ALLOCATE (br(ir,jz,kp), bz(ir,jz,kp), bp(ir,jz,kp), stat=istat)
      IF (istat .ne. 0) THEN
         STOP ' allocation error in xgrid'
      END IF
      ALLOCATE (ar(ir,jz,kp), az(ir,jz,kp), ap(ir,jz,kp), stat=istat)
      IF (istat .ne. 0) THEN
         STOP ' allocation error in xgrid'
      END IF

      IF (lstell_sym) THEN
         kp2 = kp/2
         jz2 = jz/2
         kp_odd = MOD(kp,2)
         jz_odd = MOD(jz,2)
!
!        Must be sure zmax = -zmin
!
         IF (ABS(zmax) > ABS(zmin)) THEN
            zmax = ABS(zmax)
            zmin = -zmax
         ELSE
            zmin = -ABS(zmin)
            zmax = -zmin
         END IF
      ELSE
         kp2 = kp
         jz2 = jz
         kp_odd = 0
         jz_odd = 0
      END IF

      coil_file = 'coils.' // TRIM(mgrid_ext)
      mgrid_file = 'mgrid_' // TRIM(mgrid_ext)
      extcur_file = 'extcur.' // TRIM(mgrid_ext)
      IF (lstell_sym) THEN
         WRITE (6,*) 'Stellarator symmetry IS assumed'
      ELSE
         WRITE (6,*) 'Stellarator symmetry IS NOT assumed'
      END IF
      WRITE (6, *) 'rmin = ', rmin,' rmax = ', rmax
      WRITE (6, *) 'zmin = ', zmin,' zmax = ', zmax
      WRITE (6, *) 'kp = ',  kp, ' ir = ',ir,' jz = ',jz
      PRINT *
      WRITE (6, *) 'Input  file: ',TRIM(coil_file)
      WRITE (6, *) 'Mgrid  file: ',TRIM(mgrid_file)
      WRITE (6, *) 'Extcur file: ',TRIM(extcur_file)


      iextc = 100
      CALL safe_open(iextc, istat, TRIM(extcur_file),
     &   'replace', 'formatted')
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(extcur_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', iextc
         STOP 25
      END IF

!-----------------------------------------------
!
!     PARSE FILAMENT FILE FOR NUMBER OF FIELD PERIODS
!     SPLIT INTO COIL GROUPS. DETERMINE NEXTCUR
!     COMING OUT, IGROUP+100=UNIT NO IS OPENED AND READY TO READ
!     AFTER REWINDING
!
      CALL second0(time_start)

!  parse_coils_file is in module biotsavart, made available through the
!  USE write_mgrid
      CALL parse_coils_file(coil_file)
      nfp = nfp_bs
      nextcur = SIZE(coil_group)
      IF (use_eddy) THEN
         nextcur = nextcur + 1
      END IF

      ALLOCATE (extcur(nextcur))

      CALL second0(time_end1)

      fperiod = (8*ATAN(one))/nfp
      delr = (rmax-rmin)/(ir-1)
      delz = (zmax-zmin)/(jz-1)                    
      delp = fperiod/kp

#if defined(NETCDF)
      CALL write_mgrid_nc
#else
      CALL write_mgrid_bin
#endif

      CALL second0(time_end2)

      WRITE (*, '(2(/,a,f8.3,a))')                                             &
     &   ' TIME IN PARSER = ', time_end1 - time_start,' SECONDS',              &
     &   ' TIME IN BFIELD = ', time_end2 - time_end1,' SECONDS'
      
!-----------------------------------------------
!  Print out information about each coil group
!    NB coil_group is an allocatable array of bsc_coilcoll declared in module
!    biotsavart
      WRITE(*,'(/a/)') "Extra Information about coil groups:"
      kmax = MAX(1, kp2 + kp_odd)                    ! logic from write_mgrid
      IF ((kp_odd .eq. 0) .and. lstell_sym) THEN
         kmax = MAX(kmax,kp2 + 1)
      END IF
      ALLOCATE (phi_array(kmax))
      r_ave = (rmax + rmin) / 2.
      z_ave = (zmax + zmin) / 2.
      DO k = 1, kmax
         phi_array(k) = (k - 1) * delp
      END DO
      CALL init_symmetry
      DO icoll = 1, SIZE(coil_group)
         CALL coil_group_report(coil_group(icoll), icoll, r_ave, z_ave,        &
     &                          kmax, phi_array(1:kmax))
      END DO
      DEALLOCATE (phi_array)
      WRITE(*,'(/80("*"))')

!-----------------------------------------------
!  clean up and finish
      CALL cleanup_symmetry
      CALL cleanup_biotsavart

      DEALLOCATE (br, bp, bz)
      DEALLOCATE (ar, ap, az)

      IF (mgrid_mode .eq. 'R') THEN
         WRITE (6, 205)
         WRITE (iextc, 205)
      ELSE
         WRITE (6, 200) extcur_file
      END IF

      CLOSE (iextc)

      DEALLOCATE (extcur)
      
 200  FORMAT(/,
     &       ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN SCALED',
     &       ' MODE. THE EXTERNAL',/,' CURRENTS CORRESPONDING TO THOSE',
     &       ' IN THE COILS-DOT FILE',/,' ARE GIVEN IN THE EXTCUR ARRAY',
     &       ' IN THE FILE',1x,a,'. THEY SHOULD BE ENTERED INTO THE',
     &       ' VMEC INPUT (INDATA) FILE.'/)
 205  FORMAT(/,
     &       ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN RAW',
     &       /' MODE. THE USER MUST PROVIDE THE EXTCUR ARRAY VALUES FOR',
     &       ' THE VMEC INPUT (INDATA) FILE.',
     &       /' CHOOSING EXTCUR=1 CORRESPONDS TO THE CURRENTS PRESENT',
     &       ' IN THE COILS FILE.')

      
      END SUBROUTINE task_mgrid
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_mgrid_rs()
!  Task: mgrid_rs, generate an MGRID file with Rotated, Shifted coil groups

      USE stel_kinds
      
      USE bsc_T
       !  Derived type bsc_rs
       !  subroutines  bsc_spill_coil, bsc_mean_r, bsc_construct_rs,
       !    bsc_rot_shift
      
      USE write_mgrid
       !  , only: mgrid_ext, mgrid_mode, lstell_sym,    
       !   rmin, rmax, zmin, zmax, kp, ir, jz,   
       !   br, bz, bp, kp2, jz2, kp_odd, jz_odd, coil_file, mgrid_file,     
       !   iextc, nfp, nextcur, extcur, fperiod, delr, delp, delz
       !  subroutines: write_mgrid_nc, cleanup_biotsavart
     
      USE biotsavart !  variables nfp_bs, coil_group
                     !  subroutines parse_coils_file

      USE makegrid_global, only: task, nextcur_dim, cg_shift_1,                &
     &                           cg_shift_2, cg_rot_xcent,                     &
     &                           cg_rot_theta, cg_rot_phi,                     &
     &                           cg_rot_angle, l_rot_coil_center

      USE safe_open_mod !safe_open()
      
      USE sym_check, ONLY: init_symmetry, cleanup_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                                 :: istat
      CHARACTER (LEN=100)                     :: extcur_file
      REAL (rprec)                            :: time_start
      REAL (rprec)                            :: time_end1
      REAL (rprec)                            :: time_end2
 
      INTEGER                                 :: iextcur
      INTEGER                                 :: icoil
      REAL (rprec), DIMENSION(3)              :: cgxcent
      REAL (rprec), DIMENSION(3)              :: cg_rot_xcent_use
      REAL (rprec), DIMENSION(3)              :: mean_r
      REAL (rprec)                            :: cur_coil
      REAL (rprec)                            :: cur_total

      REAL (rprec)                            :: r_ave
      REAL (rprec)                            :: z_ave
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: phi_array
      INTEGER                                 :: kmax
      INTEGER                                 :: k
      INTEGER                                 :: icoll
      
      TYPE (bsc_rs)                           :: this_bsc_rs      !  Derived type from bsc, for rotation and shift

      REAL (rprec), DIMENSION(3), PARAMETER ::                                 &
     &   zero_a3 = (/zero,zero,zero/)

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------      
      
      WRITE (*,*) ' Running makegrid with the following parameters:'
      WRITE (*,*) '   task       = ', task
      WRITE (*,*) '   mgrid_ext  = ', mgrid_ext
      WRITE (*,*) '   mgrid_mode = ', mgrid_mode
      WRITE (*,*) '   lstell_sym = ', lstell_sym
      WRITE (*,*) '   rmin       = ', rmin
      WRITE (*,*) '   rmax       = ', rmax
      WRITE (*,*) '   zmin       = ', zmin
      WRITE (*,*) '   zmax       = ', zmax
      WRITE (*,*) '   kp         = ', kp
      WRITE (*,*) '   ir         = ', ir
      WRITE (*,*) '   jz         = ', jz
      WRITE (*,*)

      IF (rmin .lt. 0.) THEN
         STOP ' rmin must be > 0 in xgrid'
      END IF
      IF (rmax .le. rmin) THEN
         STOP ' rmax must be > rmin in xgrid'
      END IF
      IF (zmax .le. zmin) THEN
         STOP ' zmax must be > zmin in xgrid'
      END IF
      IF (kp .le. 0) THEN
         STOP 'kp must be > 0 in xgrid'
      END IF

      ALLOCATE (br(ir,jz,kp), bz(ir,jz,kp), bp(ir,jz,kp), stat=istat)
      IF (istat .ne. 0) THEN
         STOP ' allocation error in xgrid'
      END IF
      ALLOCATE (ar(ir,jz,kp), az(ir,jz,kp), ap(ir,jz,kp), stat=istat)
      IF (istat .ne. 0) THEN
         STOP ' allocation error in xgrid'
      END IF

      IF (lstell_sym) THEN
         kp2 = kp/2
         jz2 = jz/2
         kp_odd = MOD(kp,2)
         jz_odd = MOD(jz,2)
!
!        Must be sure zmax = -zmin
!
         IF (ABS(zmax) > ABS(zmin)) THEN
            zmax = ABS(zmax)
            zmin = -zmax
         ELSE
            zmin = -ABS(zmin)
            zmax = -zmin
         END IF
      ELSE
         kp2 = kp
         jz2 = jz
         kp_odd = 0
         jz_odd = 0
      END IF

      coil_file = 'coils.' // TRIM(mgrid_ext)
      mgrid_file = 'mgrid_' // TRIM(mgrid_ext)
      extcur_file = 'extcur.' // TRIM(mgrid_ext)
      IF (lstell_sym) THEN
         WRITE (6,*) 'Stellarator symmetry IS assumed'
      ELSE
         WRITE (6,*) 'Stellarator symmetry IS NOT assumed'
      END IF
      WRITE (6, *) 'rmin = ', rmin,' rmax = ', rmax
      WRITE (6, *) 'zmin = ', zmin,' zmax = ', zmax
      WRITE (6, *) 'kp = ',  kp, ' ir = ',ir,' jz = ',jz
      PRINT *
      WRITE (6, *) 'Input  file: ',TRIM(coil_file)
      WRITE (6, *) 'Mgrid  file: ',TRIM(mgrid_file)
      WRITE (6, *) 'Extcur file: ',TRIM(extcur_file)

      iextc = 100
      CALL safe_open(iextc, istat, TRIM(extcur_file),
     &               'replace', 'formatted')
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(extcur_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', iextc
         STOP 25
      END IF
!
!     PARSE FILAMENT FILE FOR NUMBER OF FIELD PERIODS
!     SPLIT INTO COIL GROUPS. DETERMINE NEXTCUR
!     COMING OUT, IGROUP+100=UNIT NO IS OPENED AND READY TO READ
!     AFTER REWINDING
!
      CALL second0(time_start)

!  parse_coils_file is in module biotsavart, made available through the
!  USE write_mgrid
      CALL parse_coils_file(coil_file)
      nfp = nfp_bs
      nextcur = SIZE(coil_group)
      IF (use_eddy) THEN
         nextcur = nextcur + 1
      END IF

! Test to make sure that nextcur <= nextcur_dim
      IF (nextcur .gt. nextcur_dim) THEN
         WRITE(*,*) 'Number of coils greater than default number of ',         &
     &              'currents.'
         STOP
      END IF

      ALLOCATE (extcur(nextcur))

!-----------------------------------------------
!  Rotate and shift the coil groups
!-----------------------------------------------
!  Coils are stored in array of bsc_coilcoll named coil_group, 
!  declared in module biotsavart
!  Loop over coil groups
      WRITE(*,*)
      WRITE(*,*) ' Rotate and Shift of the Coil Groups'
      DO iextcur = 1, nextcur
         WRITE(*,*)
         IF (iextcur .eq. nextcur .and. use_eddy) THEN
            EXIT
         ELSE
            WRITE(*,*) ' Coil Group ', iextcur,' with s_name ',
     &                 coil_group(iextcur)%s_name
         END IF

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) before rs:')
!         WRITE(*,*)

!    Compute current-averaged center of coil group (cgxcent)
         cgxcent(1:3) = zero
         cur_total = zero
         DO icoil = 1, coil_group(iextcur)%ncoil
            cur_coil = coil_group(iextcur)%coils(icoil)%current
            cur_total = cur_total + cur_coil
            CALL bsc_mean_r(coil_group(iextcur)%coils(icoil),                  &
     &                      mean_r)
            cgxcent(1:3) = cgxcent(1:3) + mean_r(1:3)*cur_coil
         END DO
         IF (cur_total .ne. 0) THEN
            cgxcent(1:3) = cgxcent(1:3)/cur_total
         END IF
         IF (l_rot_coil_center(iextcur)) THEN
            cg_rot_xcent_use = cgxcent
         ELSE
            cg_rot_xcent_use = cg_rot_xcent(iextcur,1:3)
         ENDIF
 
!    Generate bsc_rs for first shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs, zero, zero, zero, zero_a3,         &
     &                         cg_shift_1(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur), this_bsc_rs)

!    Generate bsc_rs for rotation and second shift, and apply it
         CALL bsc_construct_rs(this_bsc_rs, cg_rot_theta(iextcur),             &
     &                         cg_rot_phi(iextcur),                            &
     &                         cg_rot_angle(iextcur),                          &
     &                         cg_rot_xcent_use(1:3),                          &
     &                         cg_shift_2(iextcur,1:3))
         CALL bsc_rot_shift(coil_group(iextcur), this_bsc_rs)


         WRITE(*,1000) '   Current-Averaged center of cg = ',                  &
     &      cgxcent(1:3)
         WRITE(*,1000) '   First Shift = ', cg_shift_1(iextcur,1:3)
         WRITE(*,1000) '   Center of Rotation Used =  ',                       &
     &      cg_rot_xcent_use
         WRITE(*,1000) '   Rotation theta, phi, angle = ',                     &
     &      cg_rot_theta(iextcur), cg_rot_phi(iextcur),                        &
     &      cg_rot_angle(iextcur)
         WRITE(*,1000) '   Second Shift = ', cg_shift_2(iextcur,1:3)
1000  FORMAT(a34,3(2x,es12.5))

!  Debug/test - spill the first coil in the coil group. Comment out if unneeded.
!         CALL bsc_spill_coil(coil_group(iextcur) % coils(1),                   &
!     &      'coil(1) after rs:')

      END DO
!  End loop over coil groups

      CALL second0(time_end1)

      fperiod = (8*ATAN(one))/nfp
      delr = (rmax-rmin)/(ir-1)
      delz = (zmax-zmin)/(jz-1)                    
      delp = fperiod/kp

#if defined(NETCDF)
      CALL write_mgrid_nc
#else
      CALL write_mgrid_bin
#endif

      CALL second0(time_end2)

      WRITE (*, '(2(/,a,f8.3,a))') 
     &     ' TIME IN PARSER = ', time_end1-time_start,' SECONDS',
     &     ' TIME IN BFIELD = ', time_end2-time_end1,' SECONDS'
 
!-----------------------------------------------
!  Print out information about each coil group
!    NB coil_group is an allocatable array of bsc_coilcoll declared in module
!    biotsavart
      WRITE (*,'(/a/)') "Extra Information about coil groups:"
      kmax = MAX(1,kp2 + kp_odd)                    ! logic from write_mgrid
      IF ((kp_odd .eq. 0) .and. lstell_sym) THEN
         kmax = MAX(kmax,kp2 + 1)
      ENDIF
      ALLOCATE (phi_array(kmax))
      r_ave = (rmax + rmin) / 2.
      z_ave = (zmax + zmin) / 2.
      DO k = 1, kmax
         phi_array(k) = (k - 1)*delp
      END DO
      CALL init_symmetry
      DO icoll = 1, SIZE(coil_group)
         CALL coil_group_report(coil_group(icoll), icoll, r_ave, z_ave,        &
     &                          kmax, phi_array(1:kmax))
      END DO
      DEALLOCATE (phi_array)
      WRITE (*,'(/80("*"))')

!-----------------------------------------------
!  Clean up and finish
      CALL cleanup_symmetry
      CALL cleanup_biotsavart

      DEALLOCATE (br, bp, bz)
      DEALLOCATE (ar, ap, az)

      IF (mgrid_mode .eq. 'R') THEN
         WRITE (6, 205)
         WRITE (iextc, 205)
      ELSE
         WRITE (6, 200) extcur_file
      END IF

      CLOSE (iextc)

      DEALLOCATE (extcur)
      
 200  FORMAT(/,
     &       ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN SCALED',
     &       ' MODE. THE EXTERNAL',/,' CURRENTS CORRESPONDING TO THOSE',
     &       ' IN THE COILS-DOT FILE',/,' ARE GIVEN IN THE EXTCUR ARRAY',
     &       ' IN THE FILE',1x,a,'. THEY SHOULD BE ENTERED INTO THE',
     &       ' VMEC INPUT (INDATA) FILE.'/)
 205  FORMAT(/,
     &       ' THE BFIELDS HAVE BEEN STORED IN THE MGRID FILE IN RAW',
     &       /' MODE. THE USER MUST PROVIDE THE EXTCUR ARRAY VALUES FOR',
     &       ' THE VMEC INPUT (INDATA) FILE.',
     &       /' CHOOSING EXTCUR=1 CORRESPONDS TO THE CURRENTS PRESENT',
     &       ' IN THE COILS FILE.')

      
      END SUBROUTINE task_mgrid_rs
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_circ_tor_grid()
!  task_circ_tor_grid
!    This subroutine evaluates the vector potential and 
!    magnetic field on the surface of the (circular) torus described by rmajor 
!    and aminor at nphi toroidal slices, ntheta times poloidally per slice.

!    By this point, the parameters (mgrid_ext, rmajor, aminor, nphi, ntheta) 
!    have been input.  

      USE stel_kinds
      
      USE stel_constants
      
      USE write_mgrid, only: mgrid_ext, mgrid_mode, lstell_sym,                &
     &                       rmin, rmax, zmin, zmax, kp, ir, jz,               &
     &                       coil_file, nextcur

      USE makegrid_global, only: task, rmajor, aminor, nphi, ntheta,           &
     &                           extcur_mgrid, nextcur_dim

      USE biotsavart
       !  , only: coil_group
       !  subroutines: parse_coils_file
       !  access to bsc_b in module bsc.

      IMPLICIT NONE
      
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      REAL (rprec), DIMENSION(3)             :: x
      REAL (rprec), DIMENSION(3)             :: btot
      REAL (rprec), DIMENSION(3)             :: bdum
      REAL (rprec)                           :: tfrac
      REAL (rprec)                           :: pfrac
      REAL (rprec)                           :: brr
      REAL (rprec)                           :: bpp
      REAL (rprec)                           :: brout
      REAL (rprec)                           :: btheta
      REAL (rprec)                           :: cost
      REAL (rprec)                           :: sint
      REAL (rprec)                           :: cosp
      REAL (rprec)                           :: sinp
      REAL (rprec)                           :: norm
      REAL (rprec)                           :: imagr
      REAL (rprec)                           :: imagp
      REAL (rprec)                           :: imagt
      REAL (rprec)                           :: realr
      REAL (rprec)                           :: realp
      REAL (rprec)                           :: realt
      REAL (rprec), DIMENSION(nphi,ntheta,2) :: brf
      REAL (rprec), DIMENSION(nphi,ntheta,2) :: bpf
      REAL (rprec), DIMENSION(nphi,ntheta,2) :: bthetaf
      CHARACTER (LEN=80)                     :: ctg_file
      INTEGER                                :: i
      INTEGER                                :: j
      INTEGER                                :: k
      INTEGER                                :: nmode
      INTEGER                                :: mmode
      INTEGER, PARAMETER                     :: ctg_iou = 17
      INTEGER, PARAMETER                     :: fp = 5
      !fp=5, do one field period and multiply phi mode numbers by 5
      !fp=1, do all field periods
      
!-----------------------------------------------
!  Subroutine Interface
!-----------------------------------------------

      INTERFACE
         SUBROUTINE fft_local(dat)
            USE stel_kinds
            REAL(rprec), DIMENSION(:,:,:) :: dat
         END SUBROUTINE fft_local
      END INTERFACE
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      WRITE (*,*) ' Running makegrid with the following parameters:'
      WRITE (*,*) '   task      = ', task
      WRITE (*,*) '   mgrid_ext = ', mgrid_ext
      WRITE (*,*) '   rmajor    = ', rmajor
      WRITE (*,*) '   aminor    = ', aminor
      WRITE (*,*) '   nphi      = ', nphi
      WRITE (*,*) '   ntheta    = ', ntheta
      WRITE (*,*)
      
      coil_file = 'coils.' // TRIM(mgrid_ext)
      WRITE (6, *) 'Input  file: ',TRIM(coil_file)
      
      ctg_file = 'ctg.'// TRIM(mgrid_ext)
      WRITE (*,*) 'CTG file: ', TRIM(ADJUSTL(ctg_file))
      OPEN (UNIT=ctg_iou, FILE=TRIM(ADJUSTL(ctg_file)))
      
      CALL parse_coils_file(coil_file)
      nextcur = SIZE(coil_group)
      
      x = (/0.0,0.0,0.0/)
      
      WRITE (ctg_iou,*) 'Fourier coefficients of magnetic field.'
      WRITE (ctg_iou,*) 'This magnetic geometry is 5 fold periodic in ',       &
     &                  'phi.  Only modes which preserve this ',               &
     &                  'periodicity are kept.'
      WRITE (ctg_iou,*) 'If the field is stellarator symmetric, the ',         &
     &                  'Fourier transform can be expressed in terms ',        &
     &                  'of only sines or only cosines.'
      WRITE (ctg_iou,*)
      WRITE (ctg_iou,*) 'Major radius: R0 = ', rmajor
      WRITE (ctg_iou,*) 'Minor radius: a = ', aminor
      WRITE (ctg_iou,*) 'nphi (number of points in phi) = ', nphi
      WRITE (ctg_iou,*) 'ntheta (number of points in theta) = ', ntheta
      WRITE (ctg_iou,*)
      WRITE (ctg_iou,*) 'Coordinate system:'
      WRITE (ctg_iou,*) 'r = SQRT((SQRT(x**2+y**2)-R0**2)**2+z**2)'
      WRITE (ctg_iou,*) 'THETA = ARCSIN(z/(SQRT(x**2+y**2)-R0)'
      WRITE (ctg_iou,*) 'PHI = ARCTAN(y/x)'
      WRITE (ctg_iou,*) 'Br = B.r_hat'
      WRITE (ctg_iou,*) 'Btheta = B.theta_hat'
      WRITE (ctg_iou,*) 'Bphi = B.phi_hat'
      WRITE (ctg_iou,*)
      WRITE (ctg_iou,*) 'The magnetic field is given by'
      WRITE (ctg_iou,*) 'Br(theta,phi) = 1/SQRT(ntheta*nphi) Sum_over_',       &
     &                  'm_and_n (BR_smn * (SIN(m*theta+n*phi))'
      WRITE (ctg_iou,*) 'Bphi(theta,phi) = 1/SQRT(ntheta*nphi)' //             &
     &                  'Sum_over_m_and_n(BPHI_cmn*COS(m*theta+n*phi))'
      WRITE (ctg_iou,*) 'Btheta(theta,phi) = 1/SQRT(ntheta*nphi) Sum_',        &
     &                  'over_m_and_n(BTHETA_cmn * (COS(m*theta+n*phi))'
      WRITE (ctg_iou,*) 'Where the subscripts "cmn" and"smn" is ' //           &
     &                  'shorthand for the trig function involved ' //         &
     &                  '("c" or "s") and the mode number ("m", "n)'

      WRITE (ctg_iou,*)
      WRITE (ctg_iou,245) "N","M","BR_smn","BPHI_cmn","BTHETA_cmn"

! Test to make sure that nextcur <= nextcur_dim
      
      IF (nextcur .gt. nextcur_dim) THEN
         WRITE(*,*) 'Number of coils greater than default number of ',         &
     &              'currents.'
         STOP
      END IF
! nphi, toroidal
! ntheta, poloidal
      DO i = 1, nphi
         pfrac = (i - 1)*1.0/nphi
         DO j = 1, ntheta
            tfrac = (j - 1)*1.0/ntheta
            
! Define spatial position to calculate field
            cost = COS(twopi*tfrac)
            sint = SIN(twopi*tfrac)
            cosp = COS(twopi*pfrac/fp)
            sinp = SIN(-twopi*pfrac/fp)
            x(1) = (rmajor + aminor*cost)*cosp
            x(2) = (rmajor + aminor*cost)*sinp
            x(3) = aminor*sint
            
! Calculate field
            btot = 0
            DO k = 1, nextcur
               CALL bsc_b(coil_group(k), x, bdum)
               btot = btot + bdum*extcur_mgrid(k)
            END DO
            
! Transform field to local normal coordinates
            brr = btot(1)*cosp + btot(2)*sinp
            brout = cost*brr + sint*btot(3)
            btheta = -sint*brr + cost*btot(3)
            bpp = -btot(1)*sinp + btot(2)*cosp
            
! Arrange fields for Fourier transform
            brf(i,j,1) = brout
            brf(i,j,2) = 0
            bpf(i,j,1) = bpp
            bpf(i,j,2) = 0
            bthetaf(i,j,1) = btheta
            bthetaf(i,j,2) = 0

         END DO
      END DO
      
!Fourier transform the data in theta and phi
      CALL fft_local(brf)
      CALL fft_local(bpf)
      CALL fft_local(bthetaf)
      
!Arrange the data so that it runs from negative modes to positive modes and
!write to output file
      DO i = nphi/2 + 1, NPHI
         DO j = NTHETA/2 + 1, NTHETA
            
            norm = SQRT(nphi*one)*SQRT(ntheta*one)
            realr = brf(i,j,1)/norm
            imagr = brf(i,j,2)/norm
            realp = bpf(i,j,1)/norm
            imagp = bpf(i,j,2)/norm
            realt = bthetaf(i,j,1)/norm
            imagt = bthetaf(i,j,2)/norm
            
            mmode = j - NTHETA - 1
            nmode = fp*(i - nphi - 1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
         DO j = 1, NTHETA/2
            
            norm = SQRT(nphi*one)*SQRT(ntheta*one)
            realr = brf(i,j,1)/norm
            imagr = brf(i,j,2)/norm
            realp = bpf(i,j,1)/norm
            imagp = bpf(i,j,2)/norm
            realt = bthetaf(i,j,1)/norm
            imagt = bthetaf(i,j,2)/norm
            
            mmode = j - 1
            nmode = fp*(i - nphi - 1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
      END DO

      DO i = 1, NPHI/2
         DO j = NTHETA/2 + 1, NTHETA
            
            norm = SQRT(nphi*one)*SQRT(ntheta*one)
            realr = brf(i,j,1)/norm
            imagr = brf(i,j,2)/norm
            realp = bpf(i,j,1)/norm
            imagp = bpf(i,j,2)/norm
            realt = bthetaf(i,j,1)/norm
            imagt = bthetaf(i,j,2)/norm
            
            mmode = j - NTHETA - 1
            nmode = fp*(i - 1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
         DO j = 1, NTHETA/2
            
            norm = SQRT(nphi*one)*SQRT(ntheta*one)
            realr = brf(i,j,1)/norm
            imagr = brf(i,j,2)/norm
            realp = bpf(i,j,1)/norm
            imagp = bpf(i,j,2)/norm
            realt = bthetaf(i,j,1)/norm
            imagt = bthetaf(i,j,2)/norm
            
            mmode = j - 1
            nmode = fp*(i - 1)

            WRITE(ctg_iou,240) nmode, mmode, -imagr, realp, realt
         END DO
      END DO
      
 240  FORMAT(2I6,3ES20.12)
 245  FORMAT(2A6,3A20)
      END SUBROUTINE task_circ_tor_grid

!*******************************************************************************
! SECTION IV. AUXILIARY SUBROUTINES
!*******************************************************************************
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE coil_group_report(cg, igroup, rR, zZ, nphi, phi_array)
!  Subroutine to print out information about a coil group
!  (Stored in a bsc_coilcoll)

!  JDH 2011-08-15 - First Version

      USE stel_kinds
      USE stel_constants
      
      USE bsc_T
       !  Access to derived types, magnetic field computation 
       
      USE sym_check, ONLY: check_symmetry

      IMPLICIT NONE
!-----------------------------------------------
!   A R G U M E N T   V a r i a b l e s
!-----------------------------------------------
      
      TYPE (bsc_coilcoll),  INTENT(inout)         :: cg
      INTEGER, INTENT(in)                         :: igroup
      INTEGER, INTENT(in)                         :: nphi
      REAL (rprec), INTENT(in)                    :: rR
      REAL (rprec), INTENT(in)                    :: zZ
      REAL (rprec), DIMENSION(1:nphi), INTENT(in) :: phi_array

!  cg          A bsc_coilcoll (hodling a coil group) to report on
!  igroup      Coil group number (merely reported here)
!  rR          Cylindrical R at which to report B
!  zZ          Cylindrical Z at which to report B
!  nphi        Number of phi value at which to report B (length of phi_array)
!  phi_array   Array of cylindrical phi values at which to report B
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL (rprec), PARAMETER                     :: big = 1.E40_rprec
      INTEGER, PARAMETER                          :: int_big = 987654321

      INTEGER                                     :: n_coil_fl
      INTEGER                                     :: n_coil_fc
      INTEGER                                     :: n_coil_other
      INTEGER                                     :: nfil_max_fl
      INTEGER                                     :: nfil_min_fl
      INTEGER                                     :: nfil_sum_fl
      INTEGER                                     :: nfil_sumsq_fl

      INTEGER                                     :: ic
      INTEGER                                     :: ncoil
      INTEGER                                     :: nfil
      INTEGER                                     :: iphi
      INTEGER                                     :: i

      REAL (rprec)                                :: cur_max_fl
      REAL (rprec)                                :: cur_min_fl
      REAL (rprec)                                :: cur_sum_fl
      REAL (rprec)                                :: cur_sumsq_fl
      REAL (rprec)                                :: cur_max_fc
      REAL (rprec)                                :: cur_min_fc
      REAL (rprec)                                :: cur_sum_fc
      REAL (rprec)                                :: cur_sumsq_fc
      REAL (rprec)                                :: ave_nfil_fl
      REAL (rprec)                                :: sd_nfil_fl
      REAL (rprec)                                :: cur_ave_fl
      REAL (rprec)                                :: cur_sd_fl
      REAL (rprec)                                :: cur_ave_fc
      REAL (rprec)                                :: cur_sd_fc
     
      TYPE(bsc_coil), POINTER                     :: coil
      
      REAL (rprec), DIMENSION(3)                  :: x
      REAL (rprec), DIMENSION(3)                  :: b
      
      REAL (rprec)                                :: c
      REAL (rprec)                                :: s
      REAL (rprec)                                :: br
      REAL (rprec)                                :: bphi
      REAL (rprec)                                :: bz
      REAL (rprec)                                :: cur
      REAL (rprec)                                :: phi

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------  
      
!  Find out how many coils
      ncoil = cg%ncoil
      IF (ncoil .le. 0) THEN
         WRITE(*,4000)  igroup, cg % s_name
4000  FORMAT(/80("*")/"No coils in coil group # ",i4,", group ID ",a30)
         RETURN
      ELSE
         WRITE(*,4100) igroup, cg % s_name
4100  FORMAT(/80("*")/"Coil group # ",i4,", with group ID ",a30)
      ENDIF

!  Initialize variables before loop over coils
      n_coil_fl = 0
      n_coil_fc = 0
      n_coil_other = 0
      nfil_max_fl = 0
      nfil_min_fl = int_big
      nfil_sum_fl = 0
      nfil_sumsq_fl = 0
      cur_max_fl = - big
      cur_min_fl = big
      cur_sum_fl = zero
      cur_sumsq_fl = zero
      cur_max_fc = - big
      cur_min_fc = big
      cur_sum_fc = zero
      cur_sumsq_fc = zero

!  Loop over coils
      DO ic = 1, ncoil
         coil => cg%coils(ic)
         cur = coil%current

!  Different coding, depending on c_type
         SELECT CASE (coil%c_type)
            CASE ('fil_loop','floop')
               n_coil_fl = n_coil_fl + 1
            
               nfil = SIZE(coil%xnod,2) - 1
               nfil_max_fl = MAX(nfil_max_fl, nfil)
               nfil_min_fl = MIN(nfil_min_fl,nfil)
               nfil_sum_fl = nfil_sum_fl + nfil
               nfil_sumsq_fl = nfil_sumsq_fl + nfil*nfil
            
               cur_max_fl = MAX(cur_max_fl,cur)
               cur_min_fl = MIN(cur_min_fl,cur)
               cur_sum_fl = cur_sum_fl + cur
               cur_sumsq_fl = cur_sumsq_fl + cur*cur

            CASE ('fil_circ','fcirc')
               n_coil_fc = n_coil_fc + 1
               cur_max_fc = MAX(cur_max_fc,cur)
               cur_min_fc = MIN(cur_min_fc,cur)
               cur_sum_fc = cur_sum_fc + cur
               cur_sumsq_fc = cur_sumsq_fc + cur*cur
         
            CASE DEFAULT
               n_coil_other = n_coil_other + 1
         
         END SELECT
      END DO
      
      ave_nfil_fl = nfil_sum_fl*one/MAX(n_coil_fl,1)
      sd_nfil_fl = SQRT(nfil_sumsq_fl*one/MAX(n_coil_fl,1) -                   &
     &                  ave_nfil_fl**2)
      cur_ave_fl = cur_sum_fl/MAX(n_coil_fl,1)
      cur_sd_fl =  SQRT(cur_sumsq_fl/MAX(n_coil_fl,1) -                        &
     &                  cur_ave_fl**2)
      cur_ave_fc = cur_sum_fc/MAX(n_coil_fc,1)
      cur_sd_fc = SQRT(cur_sumsq_fc/MAX(n_coil_fc,1) -                         &
     &                 cur_ave_fc**2)

!  Write out results for all coils
      WRITE (*,5000)
      WRITE (*,5100) n_coil_fl, n_coil_fc, n_coil_other, ncoil
      
      IF (n_coil_fl .gt. 0) THEN
         WRITE (*,5200) nfil_min_fl, nfil_max_fl, ave_nfil_fl,                 &
     &                  sd_nfil_fl
         WRITE (*,5300) cur_min_fl, cur_max_fl, cur_ave_fl, cur_sd_fl
      ENDIF
      
      IF (n_coil_fc .gt. 0) THEN
         WRITE (*,5400) cur_min_fc, cur_max_fc, cur_ave_fc, cur_sd_fc
      ENDIF

5000  FORMAT(/"               Loops     Circles   Other   Total")
5100  FORMAT(" # coils ",4(3x,i7))
5200  FORMAT(/" Filamentary loops, number of filaments:"/                      &
     &       "      Min     Max      Average          SD "/                    &
     &       i8,3x,i8,2(2x,f11.1))
5300  FORMAT(/" Filamentary loops, current:"/                                  &
     &       "      Min          Max         Average          SD "/            &
     &       3(2x,es12.5),2x,es9.2)
5400  FORMAT(/" Filamentary circles, current:"/                                &
     &       "      Min          Max         Average          SD "/            &
     &       3(2x,es12.5),2x,es9.2)
     
!  Loop over phi values
      WRITE(*,6000) rR, zZ
      DO iphi = 1, nphi
         phi = phi_array(iphi)
         c = COS(phi)
         s = SIN(phi)
         x(1) = rR*c
         x(2) = rR*s
         x(3) = zZ
         CALL bsc_b(cg, x, b)
         br = b(1)*c + b(2)*s
         bphi = -b(1)*s + b(2)*c
         bz = b(3)
         WRITE(*,6100) phi, br, bphi, bz
      END DO
      
      CALL check_symmetry(igroup)

6000  FORMAT(/"Magnetic field with unit current multiplier, at R,Z=",          &
     &       2(2x,es12.5)/t7,"phi",t24,"B . rhat",t38,"B . phihat",            &
     &       t52,"B . zhat")
6100  FORMAT(4(3x,es12.5))
      
      END SUBROUTINE coil_group_report
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE fft_local(dat)
!  fft_local
!    Subroutine to interface with 'cftfax_g' and 'cfft99', Fourier transform
!    subroutines found in LIBSTELL/FFTpack.  'cftfax_g' conditions the vectors 'trigs' and
!    'ifax', while 'cfft99' does the actual transform.  See those subroutines for more
!    information.

!    Performs a 2D (two 1D) Fourier transform on the input data.  The data
!    should have the form dat(N1, N2, 2) where N1 and N2 are the number of
!    data points in the first and second direction, respectively, and
!    the final index separates the real (1) and imaginary (2) parts of the
!    data.

      USE stel_kinds, only: rprec
      
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      
      REAL (rprec), DIMENSION(:,:,:)          :: dat
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: a
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: trigs
      REAL (rprec), DIMENSION(:), ALLOCATABLE :: work
      INTEGER                                 :: N
      INTEGER                                 :: i
      INTEGER                                 :: inc
      INTEGER                                 :: jump
      INTEGER                                 :: lot
      INTEGER                                 :: isign
      INTEGER                                 :: j
      INTEGER, DIMENSION(13)                  :: ifax
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      N = SIZE(dat,1)
      isign = -1  !Sign of exponent in transform
      lot = SIZE(dat,2)
      jump = N  !Number to skip between data vectors.  Since we only do
                !one transform at a time, it shouldn't matter, but it
                !seems to have problems for values larger than N
      inc = 1   !Number of array spaces between value pairs
      
      ALLOCATE(a(2*N))
      ALLOCATE(trigs(2*N))
      ALLOCATE(work(lot*N))
      
      DO j = 1, lot
         DO i = 1, N
!  The values are stored in a 1-d array with alternating real and imaginary
!  parts.
            a(2*i - 1) = dat(i,j,1)
            a(2*i) = dat(i,j,2)
         END DO

!  Do the transform
         CALL cftfax_g(N, ifax, trigs)
         CALL cfft99(a, work, trigs, ifax, inc, jump, N, 1, isign)

!  Put the new values back into the input/output array
         DO i = 1, N
            dat(i,j,1)=a(2*i - 1)
            dat(i,j,2)=a(2*i)
         END DO
      END DO

!  Reallocate to do the transform in the other index
      IF (ALLOCATED(a)) THEN
         DEALLOCATE(a)
      END IF
      ALLOCATE(a(2*lot))
      IF (ALLOCATED(trigs)) THEN
         DEALLOCATE(trigs)
      END IF
      ALLOCATE(trigs(2*lot))
      IF (ALLOCATED(work)) THEN
         DEALLOCATE(work)
      END IF
      ALLOCATE(work(lot*N))

      DO i = 1, N
         DO j = 1, lot
            a(2*j - 1) = dat(i,j,1)
            a(2*j) = dat(i,j,2)
         END DO

         CALL cftfax_g(lot, ifax, trigs)
         CALL cfft99(a, work, trigs, ifax, inc, jump, lot, 1, isign)

         DO j=1, lot
            dat(i,j,1) = a(2*j - 1)
            dat(i,j,2) = a(2*j)
         END DO
      END DO
      
      END SUBROUTINE fft_local

!*******************************************************************************
! SECTION VI. COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2011-08-15
!    Added coil_group_report subroutine, and coding to call it.
