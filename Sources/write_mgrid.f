      MODULE write_mgrid
      USE biotsavart
      USE mgrid_mod, ONLY: vn_nextcur, vn_mgmode, vn_ir,                       &
     &                     vn_jz, vn_kp, vn_nfp, vn_rmin, vn_rmax,             &
     &                     vn_zmin, vn_zmax, vn_coilgrp, vn_coilcur,           &
     &                     vn_br0, vn_bz0, vn_bp0, vn_ar0, vn_az0,             &
     &                     vn_ap0
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL (rprec), PARAMETER      :: one = 1
      CHARACTER (LEN=*), PARAMETER ::                                          &
     &   coildim(2) = (/'stringsize          ', 'external_coil_groups'/)
      CHARACTER (LEN=*), PARAMETER :: groupdim(1)= (/'external_coils'/)
      CHARACTER (LEN=*), PARAMETER :: cylcoord(3)= (/'rad','zee','phi'/)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER                                       :: ir = 121
      INTEGER                                       :: jz = 121
      INTEGER                                       :: iextc
      INTEGER                                       :: nfp
      INTEGER                                       :: nextcur
      INTEGER                                       :: istat
      INTEGER                                       :: kp
      INTEGER                                       :: kp2
      INTEGER                                       :: jz2
      INTEGER                                       :: kp_odd
      INTEGER                                       :: jz_odd
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: br
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: bz
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: bp
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: ar
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: az
      REAL (rprec), ALLOCATABLE, DIMENSION(:, :, :) :: ap
      REAL (rprec), ALLOCATABLE, DIMENSION(:)       :: extcur
      REAL (rprec)                                  :: rmin
      REAL (rprec)                                  :: rmax
      REAL (rprec)                                  :: zmin
      REAL (rprec)                                  :: zmax
      REAL (rprec)                                  :: fperiod
      REAL (rprec)                                  :: delr
      REAL (rprec)                                  :: delz
      REAL (rprec)                                  :: delp
      LOGICAL                                       ::                         &
     &   lstell_sym = .false.
      CHARACTER (LEN=1)                             :: mgrid_mode = 'S'
      CHARACTER (LEN=70)                            :: mgrid_file
      CHARACTER (LEN=70)                            :: coil_file
      CHARACTER (LEN=60)                            :: mgrid_ext
      LOGICAL                                       :: use_eddy
      CHARACTER (len=30), DIMENSION(:), ALLOCATABLE :: temp_sname

      PRIVATE                                       :: istat
C-----------------------------------------------

      CONTAINS
#if defined(NETCDF)
      SUBROUTINE write_mgrid_nc
      USE ezcdf
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER                                        :: ngrid
      INTEGER                                        :: ig
      CHARACTER (LEN=100)                            :: temp
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_br
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_bp
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_bz
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_ar
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_ap
      CHARACTER (LEN=100), ALLOCATABLE, DIMENSION(:) :: vn_az
C-----------------------------------------------
      mgrid_file = TRIM(mgrid_file) // '.nc'

      CALL cdf_open(ngrid,mgrid_file,'w',istat)
      IF (istat .ne. 0) THEN
         STOP 'Error opening WOUT.nc file VMEC WROUT'
      END IF

!
!     DEFINE DATA VARIABLES, DIMENSION NAMES
!
      CALL cdf_define(ngrid, vn_ir, ir)
      CALL cdf_define(ngrid, vn_jz, jz)
      CALL cdf_define(ngrid, vn_kp, kp)
      CALL cdf_define(ngrid, vn_nfp, nfp)
      CALL cdf_define(ngrid, vn_nextcur, nextcur)
      CALL cdf_define(ngrid, vn_rmin, rmin)
      CALL cdf_define(ngrid, vn_zmin, zmin)
      CALL cdf_define(ngrid, vn_rmax, rmax)
      CALL cdf_define(ngrid, vn_zmax, zmax)
      IF (nextcur .eq. 1) THEN
         CALL cdf_define(ngrid, vn_coilgrp,coil_group(1)%s_name)
      ELSE IF (use_eddy) THEN
         ALLOCATE(temp_sname(nextcur))
         temp_sname(1:nextcur - 1) = coil_group(1:nextcur - 1)%s_name
         temp_sname(nextcur) = 'Eddy currents'
         CALL cdf_define(ngrid, vn_coilgrp, temp_sname, dimname=coildim)
      ELSE
         CALL cdf_define(ngrid, vn_coilgrp,coil_group(1:nextcur)%s_name,
     &                   dimname=coildim)
      END IF
      CALL cdf_define(ngrid, vn_mgmode, mgrid_mode)
      CALL cdf_define(ngrid, vn_coilcur, extcur(1:nextcur),
     &                dimname=groupdim)
!
!     STORED AS 3D ARRAYS (ACTUALLY 4D, BUT CUT THROUGH IG)
!
      ALLOCATE (vn_br(nextcur), vn_bz(nextcur), vn_bp(nextcur))
      ALLOCATE (vn_ar(nextcur), vn_az(nextcur), vn_ap(nextcur))

      DO ig = 1, nextcur
         write (temp, '(a,i3.3)') "_",ig
         vn_br(ig) = vn_br0 // temp
         vn_bp(ig) = vn_bp0 // temp
         vn_bz(ig) = vn_bz0 // temp
         CALL cdf_define(ngrid, vn_br(ig), br, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_bp(ig), bp, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_bz(ig), bz, dimname=cylcoord)

         vn_ar(ig) = vn_ar0 // temp
         vn_ap(ig) = vn_ap0 // temp
         vn_az(ig) = vn_az0 // temp
         CALL cdf_define(ngrid, vn_ar(ig), ar, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_ap(ig), ap, dimname=cylcoord)
         CALL cdf_define(ngrid, vn_az(ig), az, dimname=cylcoord)
      END DO


!
!     WRITE OUT DATA
!
      CALL cdf_write(ngrid, vn_ir, ir)
      CALL cdf_write(ngrid, vn_jz, jz)
      CALL cdf_write(ngrid, vn_kp, kp)
      CALL cdf_write(ngrid, vn_nfp, nfp)
      CALL cdf_write(ngrid, vn_nextcur, nextcur)
      CALL cdf_write(ngrid, vn_rmin, rmin)
      CALL cdf_write(ngrid, vn_zmin, zmin)
      CALL cdf_write(ngrid, vn_rmax, rmax)
      CALL cdf_write(ngrid, vn_zmax, zmax)
      IF (nextcur .eq. 1) THEN
         CALL cdf_write(ngrid, vn_coilgrp, coil_group(1)%s_name)
      ELSE IF (use_eddy) THEN
         CALL cdf_write(ngrid, vn_coilgrp, temp_sname)
         DEALLOCATE(temp_sname)
      ELSE
         CALL cdf_write(ngrid, vn_coilgrp, coil_group(1:nextcur)%s_name)
      END IF

!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "	" TO COMPUTE THESE VALUES
!
      GROUPS: DO ig = 1, nextcur
         IF (ig .eq. nextcur .and. use_eddy) THEN
            br = 0
            bp = 0
            bz = 0
            ar = 0
            ap = 0
            az = 0
            extcur(nextcur) = 1
         ELSE
            CALL compute_bfield(ig)
         END IF

         CALL cdf_write(ngrid, vn_br(ig), br)
         CALL cdf_write(ngrid, vn_bp(ig), bp)
         CALL cdf_write(ngrid, vn_bz(ig), bz)

         CALL cdf_write(ngrid, vn_ar(ig), ar)
         CALL cdf_write(ngrid, vn_ap(ig), ap)
         CALL cdf_write(ngrid, vn_az(ig), az)

      END DO GROUPS

      CALL cdf_write(ngrid, vn_mgmode, mgrid_mode)
      CALL cdf_write(ngrid, vn_coilcur, extcur(1:nextcur))

      CALL cdf_close(ngrid)

      DEALLOCATE (vn_br, vn_bz, vn_bp)
      DEALLOCATE (vn_ar, vn_az, vn_ap)

      END SUBROUTINE write_mgrid_nc
#endif

      SUBROUTINE write_mgrid_bin
      USE safe_open_mod
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: igrid0 = 10
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER            :: igrid
      INTEGER            :: ig
      INTEGER            :: i
C-----------------------------------------------

      igrid = igrid0
      mgrid_file = TRIM(mgrid_file) // '.bin'
      CALL safe_open(igrid, istat, TRIM(mgrid_file),                           &
     &               'replace', 'unformatted')
      
      IF (istat .ne. 0) THEN
         WRITE (6,*) 'XGRID could not create ', TRIM(mgrid_file)
         WRITE (6,*) 'IOSTAT = ', istat,' IUNIT = ', igrid
         STOP 20
      END IF

!
!     MARK NEW-STYLE CODE WITH NEXTCUR < 0 SO READ_MGRID CAN DISTINGUISH
!     FROM OLD-STYLE MGRID FILE
!
      WRITE (igrid) ir, jz, kp, nfp, -nextcur
      WRITE (igrid) rmin, zmin, rmax, zmax
      WRITE (igrid) (coil_group(i)%s_name, i = 1, nextcur)

!
!     SET UP CYLINDRICAL COMPONENTS OF MAGNETIC FIELD ON GRID
!     SUM OVER CURRENT GROUPS IG = 1,NEXTCUR
!     NOTE: USER MUST SUPPLY SUBROUTINE "BFIELD" TO COMPUTE THESE VALUES
!
      GROUPS: DO ig = 1, nextcur

         CALL compute_bfield(ig)


         WRITE (igrid, iostat=istat) br, bp, bz
         WRITE (igrid, iostat=istat) ar, ap, az
!        WRITE(igrid,'(1p,9e12.4)') (((br(i,j,k), bz(i,j,k),              !OLD STYLE
!    1         bp(i,j,k), i=1,ir), j=1,jz),k=1,kp)
         IF (istat .ne. 0) THEN
            STOP ' Error writing bfield components'
         END IF

      END DO GROUPS 

      WRITE (igrid) mgrid_mode
      WRITE (igrid) extcur(1:nextcur)

!
!     ADD RECONSTRUCTION/OBSERVATION STUFF HERE
!

      CLOSE (igrid)

      END SUBROUTINE write_mgrid_bin

       
      SUBROUTINE compute_bfield(ig)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ig
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER             :: icoil
      INTEGER             :: numcoils
      INTEGER             :: k
      INTEGER             :: j
      INTEGER             :: i
      INTEGER             :: numfils
      INTEGER             :: curindex(1)
      REAL (rprec)        :: rad
      REAL (rprec)        :: zee
      REAL (rprec)        :: phi
      REAL (rprec)        :: extcur_ig
      REAL (rprec)        :: extcur_ig_inv
C-----------------------------------------------
      numcoils = coil_group(ig)%ncoil
      curindex = MAXLOC(abs(coil_group(ig)%coils(1:numcoils)%current))
      extcur_ig = coil_group(ig)%coils(curindex(1))%current
      
      extcur(ig) = extcur_ig
      IF (extcur_ig .eq. zero) THEN
         WRITE (6, 100) ig
      END IF
 100  FORMAT ('In COILS file, current(s) vanished for coil group ',i4)

      numfils = 0
      DO icoil = 1, numcoils
         numfils = numfils +                                                   &
     &      MAX(SIZE(coil_group(ig)%coils(icoil)%xnod, 2), 2) - 1
!
!        Compute field for unit current
!
!  JDH 2011-08-16. Comment out below - don't want to change bsc_coil data
!    Adjust B values late on with Raw test
!    Will behave differently if extcur_ig == 0.
!         IF (extcur_ig .ne. zero) THEN
!            coil_group(ig) % coils(icoil) % current =
!     1      coil_group(ig) % coils(icoil) % current /extcur_ig
!         ELSE
!            coil_group(ig) % coils(icoil) % current = 1
!         END IF
      END DO

      WRITE (6, '(/,2a)')      ' COIL GROUP          :  ',                     &
     &                         TRIM(coil_group(ig) % s_name)
      WRITE (6, '(a,i6,a,i6)') ' TOTAL COILS IN GROUP: ', numcoils,            &
     &                         ' TOTAL FILAMENTS: ', numfils

      k = 1                     ! this is always a symmetry plane
      phi = (k - 1)*delp

!$omp parallel do private(i, j, zee, rad)
      DO j = 1, jz2 + jz_odd
         zee = zmin + (j - 1)*delz
         DO i = 1, ir
            rad = rmin + (i - 1)*delr
            CALL bfield(rad, phi, zee, br(i,j,k),
     &                  bp(i,j,k), bz(i,j,k), ig)

            IF (lstell_sym) THEN
               br(i,jz + 1 - j,k) = -br(i,j,k)
               bz(i,jz + 1 - j,k) =  bz(i,j,k)
               bp(i,jz + 1 - j,k) =  bp(i,j,k)
            END IF

            CALL afield (rad, phi, zee, ar(i,j,k),                             &
     &                   ap(i,j,k), az(i,j,k), ig)

            IF (lstell_sym) THEN
               ar(i,jz + 1 - j,k) = -ar(i,j,k)
               az(i,jz + 1 - j,k) =  az(i,j,k)
               ap(i,jz + 1 - j,k) =  ap(i,j,k)
            END IF

         END DO
      END DO
!$omp end parallel do
      IF (kp .ne. 1) THEN
         WRITE (6, '(a,i4,a,i4,a)') ' K = ',k,' (OUT OF ',KP,')'
      END IF

!$omp parallel do private(i, j, k, zee, rad, phi)
      DO k = 2, kp2 + kp_odd
         phi = (k - 1)*delp
         DO j = 1, jz
            zee = zmin + (j - 1)*delz
            DO i = 1, ir
               rad = rmin + (i - 1)*delr
               CALL bfield(rad, phi, zee, br(i,j,k),                           &
     &                     bp(i,j,k), bz(i,j,k), ig)
               IF (lstell_sym) THEN
                  br(i,jz + 1 - j,kp + 2 - k) = -br(i,j,k)
                  bz(i,jz + 1 - j,kp + 2 - k) =  bz(i,j,k)
                  bp(i,jz + 1 - j,kp + 2 - k) =  bp(i,j,k)
               END IF

               CALL afield(rad, phi, zee, ar(i,j,k),                           &
     &                     ap(i,j,k), az(i,j,k), ig)
               IF (lstell_sym) THEN
                  ar(i,jz + 1 - j,kp + 2 - k) = -ar(i,j,k)
                  az(i,jz + 1 - j,kp + 2 - k) =  az(i,j,k)
                  ap(i,jz + 1 - j,kp + 2 - k) =  ap(i,j,k)
               END IF
            END DO
         END DO
         WRITE (6,'(a,i4)') ' K = ',k
      END DO
!$omp end parallel do

      IF ((kp_odd .eq. 0) .and. lstell_sym) THEN       ! another symmetry plane
         k = kp2 + 1
         phi = (k - 1)*delp
!$omp parallel do private(i, j, zee, rad)
         DO j = 1, jz2 + jz_odd
            zee = zmin + (j - 1)*delz
            DO i = 1, ir
               rad = rmin + (i - 1)*delr
               CALL bfield(rad, phi, zee, br(i,j,k),                           &
     &                     bp(i,j,k), bz(i,j,k), ig)
               br(i,jz + 1 - j,k) = -br(i,j,k)
               bz(i,jz + 1 - j,k) =  bz(i,j,k)
               bp(i,jz + 1 - j,k) =  bp(i,j,k)

               CALL afield (rad, phi, zee, ar(i,j,k),                          &
     &                      ap(i,j,k), az(i,j,k), ig)
               ar(i,jz + 1 - j,k) = -ar(i,j,k)
               az(i,jz + 1 - j,k) =  az(i,j,k)
               ap(i,jz + 1 - j,k) =  ap(i,j,k)
            END DO
         END DO
!$omp end parallel do
         WRITE (6,'(a,i4)') ' K = ',k
      END IF

      IF (mgrid_mode .eq. 'R') THEN
!  JDH 2011-080-16. Comment out below (1 line)
!         br = br*extcur_ig;  bp = bp*extcur_ig; bz = bz*extcur_ig    ! "Raw" fields
         IF (ig .lt. 10) THEN
            WRITE (iextc, 220) ig, extcur_ig, extcur_ig      ! currents for "raw" mode
         END IF
         IF (ig .ge. 10) THEN
            WRITE (iextc, 225) ig, extcur_ig, extcur_ig
         END IF
      ELSE
!  JDH 2011-08-16. Add below, scale b. 2011-08-24 - clean up extcur_ig_inv 
         IF (ABS(extcur_ig) .gt. 1.E-100_rprec) THEN
            extcur_ig_inv = 1./extcur_ig
         ELSE
            extcur_ig_inv = 1.E100_rprec
         ENDIF
!         extcur_ig_inv = 1. / MAX(extcur_ig, 1.E-100)
         br = br*extcur_ig_inv
         bp = bp*extcur_ig_inv
         bz = bz*extcur_ig_inv                                 ! "Scaled" fields

         ar = ar*extcur_ig_inv
         ap = ap*extcur_ig_inv
         az = az*extcur_ig_inv                                 ! "Scaled" fields
!  JDH 2011-08-16. end addition
         IF (ig .lt. 10) THEN
            WRITE (iextc, 210) ig, extcur_ig      ! currents for "scaled" mode
         END IF
         IF (ig .ge. 10) THEN
            WRITE (iextc, 215) ig, extcur_ig
         END IF
      END IF


 210  FORMAT('EXTCUR(', i1,')  = ', 1p,e22.14)
 215  FORMAT('EXTCUR(', i2,') = ', 1p,e22.14)
 220  FORMAT('EXTCUR(', i1,')  = ', 1p,e22.14,'/',e22.14)
 225  FORMAT('EXTCUR(', i2,') = ', 1p,e22.14,'/',e22.14)

      END SUBROUTINE compute_bfield

      END MODULE write_mgrid
