!*******************************************************************************
!>  @file check_fields.f
!>  @brief Utility to check fields written from mgrid.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Checks mgrid fields against analytic values.
!*******************************************************************************
      PROGRAM check_fields
      USE ezcdf
      USE stel_kinds

      IMPLICIT NONE

!  Local Variables
      INTEGER :: ncid
      INTEGER :: local_error
      INTEGER :: test_error
      INTEGER :: temp_int

!  Start of executable code.
      local_error = 0
      test_error = 0

!  Open the netcdf file. File name hardcoded for now.
      CALL cdf_open(ncid, 'mgrid_test.nc', 'R', local_error)
      IF (local_error .ne. 0) THEN
         WRITE (*,*) 'Open NetCDF file fialed.'
         test_error = local_error
      END IF

!  For now the tests are hard coded. Tests coils will produce fields on a 2x3
!  grid ranging from r = [0,0.5] z =[-0.5,0.5]. First check the dimensions to
!  make sure mgrid file was created as expected.
      CALL check_int(ncid, 'ir',      2, test_error)
      CALL check_int(ncid, 'jz',      3, test_error)
      CALL check_int(ncid, 'kp',      1, test_error)
      CALL check_int(ncid, 'nfp',     1, test_error)
      CALL check_int(ncid, 'nextcur', 9, test_error)

      CALL check_real(ncid, 'rmax',  0.5_dp, 0.0001_dp, test_error)
      CALL check_real(ncid, 'rmin',  0.0_dp, 0.0001_dp, test_error)
      CALL check_real(ncid, 'zmax',  0.5_dp, 0.0001_dp, test_error)
      CALL check_real(ncid, 'zmin', -0.5_dp, 0.0001_dp, test_error)

!  All the coils with canceling currents should be zero at all points in every
!  direction.
      CALL check_real_2d_all(ncid, 'ar_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ap_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'br_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bz_003',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)

      CALL check_real_2d_all(ncid, 'ar_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ap_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'br_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bz_006',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)

      CALL check_real_2d_all(ncid, 'ar_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ap_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'br_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bz_009',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)

!  These coils art all cirular so the poloidal magnetic field should be zero.
!  Note we have already tested coils 3, 6 and 9 above.
      CALL check_real_2d_all(ncid, 'bp_001',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_002',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_004',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_005',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_007',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'bp_008',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)

!  The R and Z vector potentials should be zero as well.
      CALL check_real_2d_all(ncid, 'ar_001',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ar_002',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ar_004',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ar_005',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ar_007',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'ar_008',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_001',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_002',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_004',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_005',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_007',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)
      CALL check_real_2d_all(ncid, 'az_008',  0.0_dp, 1.0E-20_dp,              &
     &                       test_error)

!  Check center of the coil along the Z axis at R. The radial components should
!  cancel along the axis.
      CALL check_real_2d(ncid, 'br_001', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_001', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_001', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_002', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_002', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_002', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_004', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_004', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_004', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_005', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_005', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_005', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_007', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_007', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_007', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_008', 0.0_dp, 1.0E-18_dp, 1, 1,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_008', 0.0_dp, 1.0E-18_dp, 1, 2,             &
     &                   test_error)
      CALL check_real_2d(ncid, 'br_008', 0.0_dp, 1.0E-18_dp, 1, 3,             &
     &                   test_error)

!  Close the netcdf file.
      CALL cdf_close(ncid, local_error)
      IF (local_error .ne. 0) THEN
         WRITE (*,*) 'Close NetCDF file fialed.'
         test_error = local_error
      END IF

      CALL EXIT(test_error)

      END PROGRAM

!*******************************************************************************
!>  @brief check an expected int value.
!>
!>  Checks an integer value for the expected result.
!>
!>  @param[in]    ncid  NetCDF file reference.
!>  @param[in]    name  Name of the
!>  @param[in]    value Expected value.
!>  @param[inout] error Current error condition.
!*******************************************************************************
      SUBROUTINE check_int(ncid, name, value, error)
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)           :: ncid
      CHARACTER (len=*), INTENT(in) :: name
      INTEGER, INTENT(in)           :: value
      INTEGER, INTENT(inout)        :: error

!  Local Variables
      INTEGER                       :: temp_int

!  Start of executable code.
      CALL cdf_read(ncid, name, temp_int)
      IF (temp_int .ne. value) THEN
         WRITE (*,1000) name, value, temp_int
         error = 666
      END IF

1000  FORMAT('Expected ',a,' value of ',i3,' recieved ',i3,'.')

      END SUBROUTINE

!*******************************************************************************
!>  @brief check an expected real value.
!>
!>  Checks an real value for the expected result with in the tolarance.
!>
!>  @param[in]    ncid      NetCDF file reference.
!>  @param[in]    name      Name of the
!>  @param[in]    value     Expected value.
!>  @param[in]    tolarance Expected value.
!>  @param[inout] error     Current error condition.
!*******************************************************************************
      SUBROUTINE check_real(ncid, name, value, tolarance, error)
      USE stel_kinds
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)           :: ncid
      CHARACTER (len=*), INTENT(in) :: name
      REAL (dp), INTENT(in)         :: value
      REAL (dp), INTENT(in)         :: tolarance
      INTEGER, INTENT(inout)        :: error

!  Local Variables
      REAL (dp)                     :: temp_real

!  Start of executable code.
      CALL cdf_read(ncid, name, temp_real)
      IF (temp_real .lt. value - tolarance .or.                                &
     &    temp_real .gt. value + tolarance) THEN
         WRITE (*,1000) name, value, tolarance, temp_real
         error = 666
      END IF

1000  FORMAT('Expected ',a,' value of ',e12.5,' with in ',e12.5,               &
     &       ' recieved ',e12.5,'.')

      END SUBROUTINE

!*******************************************************************************
!>  @brief check an expected real value.
!>
!>  Checks an real value for the expected result with in the tolarance.
!>
!>  @param[in]    ncid      NetCDF file reference.
!>  @param[in]    name      Name of the
!>  @param[in]    value     Expected value.
!>  @param[in]    tolarance Expected value.
!>  @param[in]    i         R index value.
!>  @param[in]    j         Z index value.
!>  @param[inout] error     Current error condition.
!*******************************************************************************
      SUBROUTINE check_real_2d(ncid, name, value, tolarance, i, j,             &
     &                         error)
      USE stel_kinds
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)           :: ncid
      CHARACTER (len=*), INTENT(in) :: name
      REAL (dp), INTENT(in)         :: value
      REAL (dp), INTENT(in)         :: tolarance
      INTEGER, INTENT(in)           :: i
      INTEGER, INTENT(in)           :: j
      INTEGER, INTENT(inout)        :: error

!  Local Variables
      REAL (dp), DIMENSION(2,3,1)   :: temp_real

!  Start of executable code.
      CALL cdf_read(ncid, name, temp_real)
      IF (temp_real(i,j,1) .lt. value - tolarance .or.                         &
     &    temp_real(i,j,1) .gt. value + tolarance) THEN
         WRITE (*,1000) name, value, tolarance, temp_real(i,j,1), i, j
         error = 666
      END IF

1000  FORMAT('Expected ',a,' value of ',e12.5,' with in ',e12.5,               &
     &       ' recieved ',e12.5,' at ',i3',',i3'.')

      END SUBROUTINE

!*******************************************************************************
!>  @brief check an expected real value.
!>
!>  Checks an real value for the expected result with in the tolarance at all
!>  index points.
!>
!>  @param[in]    ncid      NetCDF file reference.
!>  @param[in]    name      Name of the
!>  @param[in]    value     Expected value.
!>  @param[in]    tolarance Expected value.
!>  @param[in]    i         R index value.
!>  @param[in]    j         Z index value.
!>  @param[inout] error     Current error condition.
!*******************************************************************************
      SUBROUTINE check_real_2d_all(ncid, name, value, tolarance, error)
      USE stel_kinds
      USE ezcdf

      IMPLICIT NONE

!  Declare Arguments
      INTEGER, INTENT(in)           :: ncid
      CHARACTER (len=*), INTENT(in) :: name
      REAL (dp), INTENT(in)         :: value
      REAL (dp), INTENT(in)         :: tolarance
      INTEGER, INTENT(inout)        :: error

!  Local Variables
      REAL (dp), DIMENSION(2,3,1)   :: temp_real
      INTEGER                       :: i
      INTEGER                       :: j

!  Start of executable code.
      CALL cdf_read(ncid, name, temp_real)
      DO j = 1, 3
         DO i = 1,2
            IF (temp_real(i,j,1) .lt. value - tolarance .or.                   &
     &          temp_real(i,j,1) .gt. value + tolarance) THEN
               WRITE (*,1000) name, value, tolarance, temp_real(i,j,1),        &
     &                        i, j
               error = 666
            END IF
         END DO
      END DO

1000  FORMAT('Expected ',a,' value of ',e12.5,' with in ',e12.5,               &
     &       ' recieved ',e12.5,' at ',i3',',i3'.')

      END SUBROUTINE

