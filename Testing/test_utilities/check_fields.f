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

      IMPLICIT NONE

!  Local Variables
      INTEGER :: ncid
      INTEGER :: local_error
      INTEGER :: test_error

!  Start of executable code.
      local_error = 0
      test_error = 0

!  Open the netcdf file. File name hardcoded for now.
      CALL cdf_open(ncid, mgrid_test.nc, 'R', local_error)
      IF (local_error .ne. 0) THEN
         WRITE (*,*) 'Open NetCDF file fialed.'
         test_error = local_error
      END IF

!  Close
      CALL cdf_close(ncid, local_error)
      IF (local_error .ne. 0) THEN
         WRITE (*,*) 'Close NetCDF file fialed.'
         test_error = local_error
      END IF

      CALL EXIT(test_error)

      END PROGRAM
