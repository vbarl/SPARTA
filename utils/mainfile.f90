!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! An example for creating the main input file                                                                                   !
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
! thesis, University of Leipzig, Faculty of Physics and Earth Sciences,                                                         !
! http://www.qucosa.de/recherche/frontdoor/?tx_slubopus4frontend[id]=20746, 2016b.                                              !
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Notes:                                                                                                                        !
!                                                                                                                               !
! 1. This is just an example                                                                                                    !
!                                                                                                                               !
! 2. Purpose: creation of the main input file for Case B3: Aerosol profile, IPRT - phase A. This is example is just an          !
!             and illustration of how to create the example main file, 1 (meets the needs of Case B3).                          !
!                                                                                                                               !
! 3. The user should modify the code according to his needs                                                                     !
!-------------------------------------------------------------------------------------------------------------------------------!
program main
implicit none
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare main variables
!
real(8)        :: albedo_lambert, x, theta_min, phi                     ! single scattering albedo, angles variables
integer(8)     :: theta_0,phi_0,i                                       ! sun zenith and azimuth angles, index
integer(8)     :: n_radiance, n_angle, n_photon, n_rand                 ! explanation bellow
character(100) :: mc_file, cloudfile, angle, xmesh, ymesh, zmesh, cloud ! explanation bellow
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
write(*,*)
write(*,*) '------------------------------------------'
write(*,*) '          SPARTA main input file          '
write(*,*) '------------------------------------------'
write(*,*)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
print *, 'Create the main input file: an example'
print *
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Give the needed information
print *, 'give the name of the main input file (type "1"): '; read *, mc_file
print *
print *, 'cloud filename (type "cube_files"): '; read *, cloudfile
print *
print *, 'phase matrix filename (type "cloud"): '; read *, cloud
print *
print *, 'angle filename (type "angle"): '; read *, angle
print *
print *, 'x grid filename (type "x_grid"): '; read *, xmesh
print *
print *, 'y grid filename (type "y_grid"): '; read *, ymesh
print *
print *, 'z grid filename(type "z_grid"): '; read *, zmesh
print *
print *, 'solar zenith angle [deg] (type "30"): '; read *, theta_0
print *
print *, 'solar angles [deg] (type "0"): '; read *, phi_0
print *
print *, 'number of photons (type "100000000"):'; read *, n_photon
print *
print *, 'number of look-up-tables bins (type "10000")'; read *, n_rand
print *
print *, 'surface albedo [0,1] (type "0.0"): '; read *, albedo_lambert
print *
print *, 'number of directions (type "18"):'; read *, n_radiance
print *
print *, 'minimum value of the detector zenithal angle (type "180"):'; read *,theta_min
print *
print *, 'detector azimuth angle(type "45"):'; read *,phi
print *
!
! Write the given information in the main input file. The user should modify the lines 107 - 118 according to their needs.
! For example, if less photons are needed, i.e., 1000000, then the corresponding "write format" should be "(i7)".
open(10, file = 'input/'//trim(mc_file))
write(10,'(a50)') cloudfile
write(10,'(a50)') cloud
write(10,'(a50)') angle
write(10,'(a50)') xmesh
write(10,'(a50)') ymesh
write(10,'(a50)') zmesh
write(10,"(i9)") n_photon
write(10,"(i2,1x,i1)") theta_0,phi_0
write(10,'(f3.1)') albedo_lambert
write(10,'(i5)') n_rand
write(10,'(i2)') n_radiance
write(10,"(f5.1,1x,f4.1)") theta_min,phi

x = theta_min
do i=1,n_radiance-1
 x = x - 5
 write(10,"(f5.1,1x,f4.1)") x,phi
end do

close(10)

end program main
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!