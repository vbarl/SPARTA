!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! An example for creating the input files needed for SPARTA                                                                     !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Dr. Vasileios Barlakas                                                                                                        !
! Leibniz Institute for Tropospheric Research (TROPOS), Permosestrasse 15, 04318 Leipzig, Germany                               !
! Email: barlakas@tropos.de                                                                                                     !
!        vasileios.barlakas@uni-leipzig.de                                                                                      !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! References:                                                                                                                   !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
! 1. Barlakas, V., Macke, A., and Wendisch, M.: SPARTA – Solver for Polarized Atmospheric Radiative Transfer Applications:      !
! Introduction and application to Saharan dust fields, J. Quant. Spectrosc. Radiat. Transfer, 178, 77 – 92,                     !
! doi:10.1016/j.jqsrt.2016.02.019, 2016a.                                                                                       !
!-------------------------------------------------------------------------------------------------------------------------------!
! 2. Barlakas, V.: A New Three–Dimensional Vector Radiative Transfer Model and Applications to Saharan Dust Fields, Ph.D.       !
! thesis, University of Leipzig University of Leipzig, Faculty of Physics and Earth Sciences,                                   !
! http://www.qucosa.de/recherche/frontdoor/?tx_slubopus4frontend[id]=20746, 2016b.                                              !
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Notes:                                                                                                                        !
!                                                                                                                               !
! 1. Main rule: the input files should be created in the same way that are read within SPARTA (see subroutine INPUT).           !
!                                                                                                                               !
! 2. Here is just an example; the user might find his own way for creating the input files (always following the main rule).    !
!                                                                                                                               !
! 3. Purpose: creation of input files for Case B3: Aerosol profile, IPRT - phase A. For more information the reader is referred !
!             to the webpage http://www.meteo.physik.uni-muenchen.de/~iprt/doku.php?id=intercomparisons:b3_aerosol              !
!-------------------------------------------------------------------------------------------------------------------------------!
program input
implicit none
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare main variables
!
integer(8), parameter    :: n_z = 30                                                ! # of layers in z-direction
integer(8)               :: iz,iy,ix,i                                              ! indices
integer(8), parameter    :: n_x = 1                                                 ! # of grid cells in x-direction
integer(8), parameter    :: n_y = 1                                                 ! # of grid cells in y-direction
real(8)                  :: pi, rad, deg                                            ! pi, rads, degrees
real(8)                  :: albedo_spheroid                                         ! single scattering albedo of spheroidal particles
character(18), parameter :: ray_abs = 'tau_molabs_350.dat'                          ! scattering optical thicknesses (file)
character(20), parameter :: ray_sca = 'tau_rayleigh_350.dat'                        ! absorption optical thickness (file)
character(15), parameter :: aer_ext = 'tau_aerosol.dat'                             ! aerosol optical thickness profile (file)
character(10)            :: dummy_string                                            ! just to read string from files that are not needed
integer(4)               :: n_angle                                                 ! # of angles
real(8),allocatable      :: p11(:,:),p12(:,:),p22(:,:),p33(:,:),p34(:,:),p44(:,:)   ! total scattering phase matrix
real(8),allocatable      :: p11h(:),p12h(:),p22h(:),p33h(:),p34h(:),p44h(:)         ! spheroidal scattering phase matrix
real(8),allocatable      :: p11r(:),p12r(:),p22r(:),p33r(:),p34r(:),p44r(:)         ! Rayleigh scattering phase matrix
real(8),allocatable      :: angle(:),theta(:)                                       ! angle arrays
real(8)                  :: depol,del,del_dot                                       ! for the needs of the Rayleigh phase matrix (depolarization)
real(8)                  :: xmesh(n_x+1),ymesh(n_y+1),zmesh(n_z+1)                  ! x,y, and z coordinanes
real(8)                  :: d_z,sca_mol_1(n_z+1),abs_mol_1(n_z+1)                   ! vertical extend, molecular scattering and absorption coeff.
real(8)                  :: aer_sca_1(n_z+1),aer_abs_1(n_z+1),aer_ext_1(n_z+1)      ! aerosol scattering, absorption, and extinction coeff.
real(8)                  :: scattering(n_x,n_y,n_z+1)                               ! total volumetric scattering coefficient
real(8)                  :: absorption(n_x,n_y,n_z+1)                               ! total volumetric absorption coefficient
real(8)                  :: extinction(n_x,n_y,n_z+1)                               ! total volumetric extinction coefficient
real(8)                  :: albedo(n_x,n_y,n_z)                                     ! total volumetric scattering coefficient
integer(4)               :: pm_index(n_x, n_y, n_z)                                 ! phase matrix index
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
!
! PART 1: Read the given files and conduct the needed calculations
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! initialize phase matrix index array
pm_index = 0d0

! open the angle file; scattering angle
open(11,file='info/angles',action = 'read')

! read the total # of scattering angles
read(11,*) n_angle

! initialize dimension of angle arrays
allocate(angle(1:n_angle));angle=0d0
allocate(theta(1:n_angle));theta=0d0

! read the scatteing angle
do i = 1,n_angle
 read(11,*) angle(i)
end do

! define pi, degress, and radians
pi = 4d0*atan(1d0)
deg = 180d0/pi
rad = 1d0/deg

! convert degrees to radians
theta = angle*rad

! initialize dimension of rayleigh scattering phase matrix
allocate(p11r(1:n_angle))
allocate(p12r(1:n_angle))
allocate(p22r(1:n_angle))
allocate(p33r(1:n_angle))
allocate(p34r(1:n_angle))
allocate(p44r(1:n_angle))

! initialize dimension of spheroidal scattering phase matrix
allocate(p11h(1:n_angle))
allocate(p12h(1:n_angle))
allocate(p22h(1:n_angle))
allocate(p33h(1:n_angle))
allocate(p34h(1:n_angle))
allocate(p44h(1:n_angle))

! open the files where the scattering phase matrix of the spheroidal particles are included
open(17,file='info/spheroidal.p11',action='read')
open(18,file='info/spheroidal.p12',action='read')
open(19,file='info/spheroidal.p22',action='read')
open(20,file='info/spheroidal.p33',action='read')
open(21,file='info/spheroidal.p44',action='read')
open(22,file='info/spheroidal.p34',action='read')

! read the the scattering phase matrix
do i = 1,n_angle
 read(17,*) p11h(i)
 read(18,*) p12h(i)
 read(19,*) p22h(i)
 read(20,*) p33h(i)
 read(21,*) p44h(i)
 read(22,*) p34h(i)
end do

! depol is the depolarization factor of the Rayleigh scattering phase matrix
! the reader is referred to equations 4.2 and 4.3 in Chapter 4, Barlakas, 2016b
depol = 0.03
del = (1d0 - depol)/(1d0 + depol/2d0)
del_dot = (1d0 - 2d0*depol)/(1d0 - depol)

! create the Rayleigh scattering phase matrix
do i = 1,n_angle
 p11r(i) = (3d0/4d0) * del * (1d0+(cos(theta(i)))**2d0) + (1d0 - del)
 p12r(i) = -(3d0/4d0) * del * (sin(theta(i))**2d0)
 p22r(i) = (3d0/4d0) * del * (1d0+(cos(theta(i)))**2d0)
 p33r(i) = (3d0/2d0) * del * cos(theta(i))
 p44r(i) = (3d0/2d0) * del * del_dot * cos(theta(i)) + (1d0 - del)
end do

p44r = 0d0
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Open the files where the scattering information is found
! read the string information out
open(10, file = 'info/'//trim(ray_sca), status = 'old')
read(10,*) dummy_string

open(11, file = 'info/'//trim(ray_abs), status = 'old')
read(11,*) dummy_string

open(12, file = 'info/'//trim(aer_ext), status = 'old')
read(12,*) dummy_string

! read the altitude and the properties information
do iz = n_z+1, 1, -1
  read(10,*) zmesh(iz), sca_mol_1(iz)
  read(11,*) zmesh(iz), abs_mol_1(iz)
  read(12,*) zmesh(iz), aer_ext_1(iz)
end do
close(10);close(11);close(12)

! this is the single scattering albedo of the spheroidal particles (is found in sizedistr_spheroid.dat file)
albedo_spheroid = 0.787581d0

! derive the scattering and absorption coefficients from the extinction coeff and the single scattering albedo
aer_sca_1 = aer_ext_1 * albedo_spheroid
aer_abs_1 = aer_ext_1 * (1d0 - albedo_spheroid)

! derive the total volumetric properties
do iz = 1,n_z+1
  do ix = 1, n_x
    do iy = 1, n_y
     extinction(ix,iy,iz) = sca_mol_1(iz) + abs_mol_1(iz) + aer_ext_1(iz)
     scattering(ix,iy,iz) = sca_mol_1(iz) + aer_sca_1(iz)
     absorption(ix,iy,iz) = abs_mol_1(iz) + aer_abs_1(iz)
    end do
  end do
end do

!albedo = scattering/extinction
do iz = 1,n_z
  do ix = 1, n_x
    do iy = 1, n_y
      albedo(ix,iy,iz) = scattering(ix,iy,iz)/extinction(ix,iy,iz)
    end do
  end do
end do

! initialize dimension of the total scattering phase matrix
allocate(p11(1:n_z,1:n_angle));p11 = 0.d0
allocate(p12(1:n_z,1:n_angle));p12 = 0.d0
allocate(p22(1:n_z,1:n_angle));p22 = 0.d0
allocate(p33(1:n_z,1:n_angle));p33 = 0.d0
allocate(p34(1:n_z,1:n_angle));p34 = 0.d0
allocate(p44(1:n_z,1:n_angle));p44 = 0.d0

! calculate the total scattering phase matrix (see equation 5.1 in Barlakas, 2016b)
do i = 1,n_angle
 do iz = 1, n_z

  p11(iz,i) = (sca_mol_1(iz) * p11r(i) + aer_sca_1(iz) * p11h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

  p12(iz,i) = (sca_mol_1(iz) * p12r(i) + aer_sca_1(iz) * p12h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

  p22(iz,i) = (sca_mol_1(iz) * p22r(i) + aer_sca_1(iz) * p22h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

  p33(iz,i) = (sca_mol_1(iz) * p33r(i) + aer_sca_1(iz) * p33h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

  p34(iz,i) = (sca_mol_1(iz) * p34r(i) + aer_sca_1(iz) * p34h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

  p44(iz,i) = (sca_mol_1(iz) * p44r(i) + aer_sca_1(iz) * p44h(i))&
  /(sca_mol_1(iz)+aer_sca_1(iz))

 end do
end do
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! PART 2: Create the input files of SPARTA
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! create x, y, and z grids
open(10, file = 'input/z_grid')
write(10,*), n_z

do iz = 1, n_z+1
 write(10,*) zmesh(iz)
end do
close(10)

open(10, file = 'input/x_grid')
xmesh(1) = 0.d0
xmesh(2) = 100.d0
write(10,*) n_x

do ix = 1, n_x+1
 write(10,*) xmesh(ix)
end do

close(10)

open(10, file = 'input/y_grid')

! In the desctiption of the case, there is no information about the dimensions of the x and y. Since we are reffering to a 1D
! case, we could define the extend by ourselves. Here, we artificially defined 100 km. However, another suggestion might be
! the d_z, the vertical extend in the z-direction (here, 30 km).
d_z = zmesh(n_z) - zmesh(1)

ymesh(1) = 0.d0
ymesh(2) = 100.d0

write(10,*) n_y
do iy = 1, n_y+1
 write(10,*) ymesh(iy)
end do
close(10)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Cloud information for each grid box (binary files)
open(10, file = 'input/cube_files.ext', form = 'unformatted')
open(11, file = 'input/cube_files.w0', form = 'unformatted')
open(12, file = 'input/cube_files.abs', form = 'unformatted')
open(13, file = 'input/cube_files.sca', form = 'unformatted')

do ix = 1, n_x
  do iy = 1, n_y
    write(10) (extinction(ix,iy,iz), iz = 1, n_z+1)
    write(11) (albedo(ix,iy,iz), iz = 1, n_z)
    write(12) (absorption(ix,iy,iz), iz = 1, n_z+1)
    write(13) (scattering(ix,iy,iz), iz = 1, n_z+1)
  end do
end do
close(10); close(11); close(12); close(13)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! scattering angle  (binary file)
open(10,file='input/angle',action='write',form = 'unformatted')
write(10) n_angle

do i=1,n_angle
  write(10) angle(i)
end do
close(10)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Total scattering phase matrix (binary files)
open(23, file = 'input/cloud.p11',status = 'replace', form = 'unformatted')
open(24, file = 'input/cloud.p12',status = 'replace', form = 'unformatted')
open(25, file = 'input/cloud.p22',status = 'replace', form = 'unformatted')
open(26, file = 'input/cloud.p33',status = 'replace', form = 'unformatted')
open(27, file = 'input/cloud.p34',status = 'replace', form = 'unformatted')
open(28, file = 'input/cloud.p44',status = 'replace', form = 'unformatted')


do iz = 1, n_z
 write(23) (p11(iz,i), i = 1, n_angle)
 write(24) (p12(iz,i), i = 1, n_angle)
 write(25) (p22(iz,i), i = 1, n_angle)
 write(26) (p33(iz,i), i = 1, n_angle)
 write(27) (p34(iz,i), i = 1, n_angle)
 write(28) (p44(iz,i), i= 1, n_angle)
end do
close(23);close(24);close(25);close(26);close(27);close(28)

! Cloud file: information for the number of the total normalized scattering phase matrices and the name list of the total normalized 
! phase matrix elements
open(29, file = 'input/cloud', status = 'replace', form = 'formatted')

write(29, "(I2)" ) 30             !# of phase functions
write(29, "(A9)") "cloud.p11"
write(29, "(A9)") "cloud.p12"
write(29, "(A9)") "cloud.p22"
write(29, "(A9)") "cloud.p33"
write(29, "(A9)") "cloud.p34"
write(29, "(A9)") "cloud.p44"

close(29)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Phase matrix index calculation
do iz = 1, n_z
 do ix = 1, n_x
  do iy = 1, n_y
   pm_index(ix,iy,iz) = iz
  end do ! iy-loop
 end do ! ix-loop
end do ! iz-loop

open(30, file = 'input/cube_files.pm-index', status = 'replace', form = 'unformatted')

! Create the phase matrix index file (binary file)
do ix = 1, n_x
 do iy = 1, n_y
  write(30) (pm_index(ix,iy,iz), iz = 1, n_z)
 end do
end do
close(30)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Contains the name list of the cloud files.
open(6000,file='input/cube_files',status='replace',form='formatted')

write(6000, "(A14)") "cube_files.ext"
write(6000, "(A13)") "cube_files.w0"
write(6000, "(A14)") "cube_files.sca"
write(6000, "(A14)") "cube_files.abs"
write(6000, "(A19)") "cube_files.pm-index"

close(6000)
!__________________________________________________________________________
!
end program input
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!