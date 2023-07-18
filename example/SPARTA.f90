!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Monte Carlo Radiative Transfer Solver SPARTA                                                                                  !
! SPARTA: Solver for Polarized Atmospheric Radiative Transfer Applications                                                      !
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
!###############################################################################################################################!
!########################################################### SPARTA ############################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
program SPARTA
implicit none
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare main variables found in INPUT subroutine
!
real(8), allocatable   :: theta_hauf(:)                             ! scattering phase function check
real(8), allocatable   :: angle(:)                                  ! scattering angles
real(8), allocatable   :: xmesh(:), ymesh(:), zmesh(:)              ! grid coordinates [m or km]
real(8), allocatable   :: d_zmesh(:)                                ! variable grid size in z [m or km]
real(8)                :: d_xmesh, d_ymesh                          ! constant grid sizes in x, y [m or km]
real(8), allocatable   :: beta_x(:,:,:)                             ! 3D volumme extinction coeff. [m^-1 or km^-1]
real(8), allocatable   :: beta_abs(:,:,:), beta_sca(:,:,:)          ! 3D volume absorption and scattering coeff. [m^-1 or km^-1]
real(8), allocatable   :: w_0(:,:,:)                                ! single scattering albedo, co-albedo
real(8), allocatable   :: co_w_0(:,:,:)                             ! single scattering albedo, co-albedo
real(8), allocatable   :: tau_field(:,:)                            ! optical thickness
real(8)                :: tau                                       ! temporary variable for optical thickness caluculation
real(8)                :: theta_0, phi_0                            ! solar zenith and azimuth angle
real(8)                :: costheta_0, sintheta_0                    ! cosine and sine of the solar zenith angle
real(8)                :: cosphi_0, sinphi_0                        ! cosine and sine of the solar azimuth angle
real(8)                :: cloud_cover_true                          ! cloud cover
real(8)                :: help1, help2, test1                       ! temporary variables: helping calculations
integer(4)             :: n_angle                                   ! # of angles
integer(4)             :: xmin, xmax, ymin, ymax, zmin, zmax        ! cloud boundaries in grid
integer(8)             :: photon_total, photon                      ! # of photons, for photon loop
integer(8)             :: ix, iy, iz                                ! grid coordinates
character(100)         :: mainfile, cubefile, pm_file, anglefile    ! string variables for input files
character(100)         :: xmeshfile, ymeshfile, zmeshfile           ! string variables for input files
integer(8)             :: i, j, k, i_angle                          ! loop indices, scattering angle index,
integer(8)             :: N_scat, N_refl                            ! counters for scattering and reflecting events
integer(4)             :: n_pm, g                                   ! # of scattering phase matrices
real(8)                :: tau_mean                                  ! optical thickness, model aread mean
integer(4),allocatable :: pm_index(:,:,:)                           ! scattering phase matrix index (ix,iy,iz)
character(100)         :: s11, s12, s22, s33, s34, s44              ! string variables for reading main cloud file
real(8),allocatable    :: p11(:), p12(:), p22(:)                    ! total scattering phase matrix elements
real(8),allocatable    :: p33(:), p34(:), p44(:)                    ! total scattering phase matrix elements
real(8)                :: theta_det                                 ! zenithal angle of detector in rad
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare main variables found in INIT subroutine
!
integer(4)             :: iseed, iflag                              ! variables for the radnom number generator ran(idum)
integer(4)             :: i_rand, n_rand, n_rand_p1                 ! # of look-up-table bins
real(8)                :: d_n_rand, d_n_rand_p1                     ! dble(n_rand_p1)
integer(8)             :: start_time, end_time, i_s                 ! start and end time of simulations, index
integer(8)             :: store_every                               ! indicates every when to store intermediate results
integer(8)             :: theta_i                                   ! scattering angle index
character(5)           :: zone                                      ! the time difference with respect to UTC, in the form hhmm
character(8)           :: date                                      ! date, in form CCYYMMDD, where CCYY is the four-digit year
character(10)          :: time                                      ! time, in the form hhmmss.sss, ss.sss secs and millisecs
integer                :: values(8)                                 ! integer array of 8 elements describing CCYYMMDD, hhmmss.sss
real(8), allocatable   :: angle_u(:), angle_o(:), d_angle(:)        ! scatt. angle bin boundaries
real(8), allocatable   :: arg(:)                                    ! solid angles
real(8)                :: theta_0_true                              ! interpolated zenith scattering angle
real(8), allocatable   :: theta_true(:,:)                           ! LUT for interpolated zenith scattering angle
real(8)                :: phi_true                                  ! azimuth scattering angle
real(8), allocatable   :: costheta_scat(:,:), sintheta_scat(:,:)    ! look-up tables (LUT) for zenith scattering angles
real(8), allocatable   :: cosphi_scat(:), sinphi_scat(:)            ! LUT for azimuth scattering angles
real(8)                :: t, eps, r                                 ! helping vartiable, 'cloud-epsilon' = 1d-9, random namber
real(8), allocatable   :: pintegral(:)                              ! normalization factor for scattering phase function
real(8)                :: pi, pi2, twopi, rad, deg                  ! pi, pi/2, 2pi, rads, degrees
real(8)                :: albedo_lambert, co_albedo_lambert         ! single scattering albedo and co-albedo
real(8)                :: theta_lambert_max                         ! mimimum theta to avoid endless horizontal flights below cloud
real(8)                :: angle_rand                                ! helping variable used in LUT
real(8), allocatable   :: k_x_choice(:),k_y_choice(:),k_z_choice(:) ! direction cosines of selected rad. (detector)
real(8), allocatable   :: theta_choice(:), phi_choice(:)            ! zenith and az. of selected rad. (detector)
real(8)                :: costheta_0_in, sintheta_0_in              ! cosines and sines of solar zenith angle
real(8)                :: cosphi_0_in, sinphi_0_in                  ! cosines and sines of solar azimuth angle
real(8)                :: cloudtop, cloudbottom, cloudheight        ! cloud dimensions in z
real(8)                :: cloudfront, cloudback, clouddepth         ! cloud dimensions in y
real(8)                :: cloudleft, cloudright, cloudwidth         ! cloud dimensions in x
real(8)                :: cloudtop_in, cloudbottom_in               ! cloud dimensions just inside cloud in z
real(8)                :: cloudfront_in, cloudback_in               ! cloud dimensions just inside cloud in y
real(8)                :: cloudleft_in, cloudright_in               ! cloud dimensions just inside cloud in x
real(8)                :: cloudtop_out, cloudbottom_out             ! above, below cloud
integer, parameter     :: d_theta_i = 10                            ! step width in calculating cumulative scattering phase function
real(8), allocatable   :: p11_cum(:,:)                              ! cum. scattering phase function LUT
real(8), allocatable   :: p11_all(:,:), p12_all(:,:), p22_all(:,:)  ! set of scattering phase matrix elements
real(8), allocatable   :: p33_all(:,:), p34_all(:,:), p44_all(:,:)  ! set of scattering phase matrix elements
real(8), allocatable   :: p11_arg_int_sum(:)                        ! scattering  phase function * sold angle
real(8)                :: d_p11_arg_int_sum                         ! helping variable to construct p11_arg_int_sum
real(8), allocatable   :: p11_scat(:,:),p12_scat(:,:),p22_scat(:,:) ! scattering phase matrix elements LUT
real(8), allocatable   :: p33_scat(:,:),p34_scat(:,:),p44_scat(:,:) ! scattering phase matrix elements LUT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare variables for net horizontal flux H = 1 - R - T - A
!
real(8), allocatable   :: flux_in(:,:), flux_S(:,:)                 ! Incoming flux scalar (_in) and vector (_S)
real(8), allocatable   :: transmittance_diffuse(:,:)                ! Diffuse downward flux (transmittance) - related to count
real(8), allocatable   :: transmittance_direct(:,:)                 ! Direct downward flux (transmittance) - related to count
real(8), allocatable   :: transmittance_total(:,:)                  ! Total transmittance
integer(8)             :: count                                     ! indicates whether a scattering take place (before surface)
real(8), allocatable   :: albedo(:,:)                               ! Upward flux (albedo)
real(8), allocatable   :: net_horizontal_transport(:,:)             ! Net horizontal photon transport
real(8), allocatable   :: absorptance(:,:), absorptance_ground(:,:) ! Absorbed flux, flux absorbed from the ground
!
! Declare scalar and Stokes weight
!
real(8), dimension(1:4):: S_in                                      ! Stokes weight vector (unpolarized extraterrestrial irradiance)
real(8)                :: weight, radiance_weight                   ! photon/radiance weight (scalar)
real(8), dimension(1:4):: S_out                                     ! corresponds to the output (scattered) Stokes weight vector
!
! Declare variables related to photon position
!
real(8)                :: x, y, z, zpos                             ! photon position
real(8)                :: k_x, k_y, k_z, k_x_new, k_y_new, k_z_new  ! incident and scattering directions (direction cosines)
real(8)                :: slen                                      ! for normalization purposes (direction cosines)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare variables related to the rotations of the Stokes vector (i.e., rotation angles). Subroutines polarization and
! LocalEstimateMethod (LEM)
!
real(8)                :: cosi, cos2i, sin2i                        ! first rotation angle (cosines and sines), see polarization
real(8)                :: cosj, cos2j, sin2j                        ! second rotation angle (cosines and sines), see polarization
real(8)                :: cosk, cos2k, sin2k                        ! first rotation angle for (cosines and sines), see LEM
real(8)                :: cosh, cos2h, sin2h                        ! second rotation angle for (cosines and sines), see LEM
real(8)                :: help3_1, help3_2                          ! temporary variables: helping rotation calculations (polarization)
real(8)                :: help4_1, help4_2, help4a, help4b, help4c  ! temporary variables: helping rotation calculations (LEM)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare the variables related to radiances (scalar and vector)
!
real(8), allocatable   :: select_radiance(:,:,:)                    ! scalar selected radiances for each pixel
real(8), allocatable   :: stokes_select_radiance(:,:,:,:)           ! vector (Stokes) selected radiances for each pixel
real(8), allocatable   :: DOP(:,:,:)                                ! degree of polarization for each pixel
!
! Declare variables related to errors
!
real(8)                :: mean_pixel_error                          ! this is for the intercomparison, I3RC (phase 2)
real(8), allocatable   :: error_to_mean(:)                          ! this is for the intercomparison, I3RC (phase 2)
real(8), allocatable   :: stokes_error_to_mean(:,:)                 ! the same but for the Stokes vector
real(8), allocatable   :: error_select_radiance(:,:,:)              ! statistical Monte Carlo error (scalar)
real(8), allocatable   :: error_stokes_select_radiance(:,:,:,:)     ! statistical Monte Carlo error (vector)
real(8), allocatable   :: data(:)                                   ! a vector used for statistics
integer(4)             :: i_data                                    ! statistics
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare variables found in PhotonPath subroutine: Check the subroutine for more information
!
real(8)                :: t_ray(3), t_ray_min                       ! ray parameter, minimum of t_ray (tracking the photon position)
integer(4)             :: i_ray_min                                 ! position of the minimum ray
real(8)                :: local_beta_x                              ! local volume extinction coefficient
real(8)                :: local_beta_sca, local_beta_abs            ! local volume scattering and absorption coefficients
real(8)                :: enough                                    ! randomly chosen volume scattering coefficient
real(8)                :: step_back                                 ! correct photon position ("overjumping")
real(8)                :: tau_sum, sca_sum, abs_sum                 ! extinction, scattering, and absorption tau (along the path)
integer(4)             :: up_or_down_x, up_or_down_y, up_or_down_z  ! poynting to the x, y, z axis (positive, up or negative, down)
integer(4)             :: delta_i(-1:1)                             ! indicating up or down depending on up_or_down_*
integer(8)             :: ihelp1                                    ! helping variable
real(8)                :: theta_inc, theta_sca, sintheta_s          ! incident and scattered zenithal angles, sine of theta_sca
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare variables found in PhotonPathDetector subroutine
!
integer(8)             :: ix2detector,iy2detector,iz2detector       ! grid box for the position of the detector
real(8)                :: x2detector, y2detector, z2detector        ! posotion of detector
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare variables found in subroutine lambert_reflection
!
real(8)                :: zenith                                    ! random zenithal (global zenith angle, see subroutine reflection
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Declare the rest of variables
!
integer(8)             :: i_radiance                                ! loop index
integer(4)             :: n_radiance                                ! number of selected radiances
integer(4)             :: i_pm                                      ! index of scattering phase matrices
real(8)                :: phi_inc, phi_sca, d_phi                   ! azimuth of the incident and scattered directions, and their difference
real(8)                :: cosnew                                    ! cosine of the randomly chosen scattering anlge
real(8)                :: p11_inter,p12_inter,p22_inter             ! interpolated components of the scattering phase matrix from the LUT
real(8)                :: p33_inter,p34_inter,p44_inter             ! interpolated components of the scattering phase matrix from the LUT
! 
! Found in LEM subroutine
!
real(8)                :: phi_det                                   ! azimuth zenithal angle of detector in radians
real(8)                :: d_phi_LE                                  ! difference between the incident and the detector direction for the LEM
real(8)                :: sintheta_LE                               ! sine of zenithal angle of detector
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
!
!
!
!
!
!
!
!
!-------------------------------------------------------------------------------------------------------------------------------!
!######################################################### MAIN PROGRAM ########################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Read the main ino
!
call INPUT   ! Read the main input files
call INIT    ! Basic calculations and LUT constraction
!
N_scat = 0   ! Counter for scattering events
N_refl = 0   ! Counter for reflection events
!
!-------------------------------------------------------------------------------------------------------------------------------!
!################################################### START MONTE CARLO LOOP ####################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
do photon = 1, photon_total
  call init_photon                                    ! Initiate photons
  call PhotonPath                                     ! Tracke the photon path
   do while (z .lt. cloudtop .and. weight .gt. eps)   ! Trace photon back to cloud top when reflecting surface present
    if (z .gt. cloudbottom) then
        call polarization                             ! Scattering event and change of polarization
        call PhotonPath                               ! Tracke the photon path
    else
        if (count .gt. 0) then                        ! Diffuse radiation
            transmittance_diffuse(ix,iy) = transmittance_diffuse(ix,iy) + weight
        else                                          ! Direct transmitted radiation
            transmittance_direct(ix,iy) = transmittance_direct(ix,iy) + weight
        end if
        call reflection                               ! Reflection event
        call PhotonPath                               ! Tracke the photon path
    end if
    end do ! while
    albedo(ix, iy) = albedo(ix, iy) + weight          ! Albedo calculation
    if (mod(photon,store_every) .eq. 0) call STATUS
end do! photon-loop
call OUTPUT
!-------------------------------------------------------------------------------------------------------------------------------!
!##################################################### END MONTE CARLO LOOP ####################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
! List of subroutines
!
!-------------------------------------------------------------------------------------------------------------------------------!
contains
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine INPUT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Initialization: reads main input file and assigns input to variables                                                          !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
write(*,*)
write(*,*) '------------------------------------------'
write(*,*) 'provide the name of the main input file = '
write(*,*)
read (*,*) mainfile                                   ! main input file
write(*,*)
write(*,*) '------------------------------------------'
write(*,*) '                  SPARTA                  '
write(*,*) '------------------------------------------'
write(*,*)

open(5, file = 'input/'//trim(mainfile),status='old')
read(5,'(a80)') cubefile
read(5,'(a80)') pm_file
read(5,'(a80)') anglefile
read(5,'(a80)') xmeshfile
read(5,'(a80)') ymeshfile
read(5,'(a80)') zmeshfile
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Scattering angles at exact forward and backward directions and the coresponding values of the scattering phase matrix are not !
! usually given. If given, the code has to be slightly modified, since there is no need for interpolations. Do the following:   !
!                                                                                                                               !
! (1) uncomment line 319       : in order not to inflence the computation of the cumulative scattering phase phunction          !
! (2) comment lines 379 - 418  : no extrapolation is needed                                                                     !
! (3) uncomment lines 426 - 433: read the given phase matrix correctly                                                          !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
open(9, file = 'input/'//trim(anglefile),action='read', form = 'unformatted')
read(9) n_angle                                       ! read total number of scattering angles
!n_angle = n_angle - 2                                ! if angles at exact forward and backward directions are given, comment

allocate( angle(0:n_angle+1) )                        ! scatt. angle
allocate( angle_u(0:n_angle+1), angle_o(0:n_angle+1) )! lower and upper bound
allocate( d_angle(1: n_angle) )                       ! angular interval size
allocate( arg(0:n_angle+1) )                          ! solid angle interval

allocate( p11_arg_int_sum(0:n_angle+1 + d_theta_i) );p11_arg_int_sum= 0d0

! read scattering angles
do i = 1, n_angle
 read(9) angle(i)
end do
close(9)

angle(0) = 0d0                                        ! scatt. angle at exact forward direction
angle(n_angle + 1) = 180d0                            ! scatt. angle at backward forward direction

read(5,*) photon_total                                ! # of incoming photons
read(5,*) theta_0, phi_0                              ! solar zenith and solar azimuth angle
read(5,*) albedo_lambert                              ! isotropic surface albedo
read(5,*) n_rand                                      ! # of bins in pm lookup table
read(5,*) n_radiance                                  ! directions of selected radiances

allocate( theta_choice(1:n_radiance), phi_choice(1:n_radiance) )
do i = 1, n_radiance
  read(5,*) theta_choice(i), phi_choice(i)            ! read the viewing zenithal and azimuth angles
end do
close(5)

open(55,file = 'input/'//trim(pm_file),action='read')             !read components of scattering phase matrix
 read(55,*) n_pm
 read(55,*) s11
 read(55,*) s12
 read(55,*) s22
 read(55,*) s33
 read(55,*) s34
 read(55,*) s44
close(55)

open(12,file='input/'//trim(s11),action='read',status='old', form = 'unformatted')
open(13,file='input/'//trim(s12),action='read',status='old', form = 'unformatted')
open(14,file='input/'//trim(s22),action='read',status='old', form = 'unformatted')
open(15,file='input/'//trim(s33),action='read',status='old', form = 'unformatted')
open(16,file='input/'//trim(s34),action='read',status='old', form = 'unformatted')
open(17,file='input/'//trim(s44),action='read',status='old', form = 'unformatted')

allocate( p11_all(1:n_pm, 0:n_angle + 1) ); p11_all  = 0d0
allocate( p12_all(1:n_pm, 0:n_angle + 1) ); p12_all  = 0d0
allocate( p22_all(1:n_pm, 0:n_angle + 1) ); p22_all  = 0d0
allocate( p33_all(1:n_pm, 0:n_angle + 1) ); p33_all  = 0d0
allocate( p34_all(1:n_pm, 0:n_angle + 1) ); p34_all  = 0d0
allocate( p44_all(1:n_pm, 0:n_angle + 1) ); p44_all  = 0d0
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
!   Lines 379 - 418: in case the scattering angles at exact forward and backward directions are NOT given (if given comment out)!
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
do i=1,n_pm
 read(12) (p11_all(i,i_angle), i_angle = 1, n_angle)
 read(13) (p12_all(i,i_angle), i_angle = 1, n_angle)
 read(14) (p22_all(i,i_angle), i_angle = 1, n_angle)
 read(15) (p33_all(i,i_angle), i_angle = 1, n_angle)
 read(16) (p34_all(i,i_angle), i_angle = 1, n_angle)
 read(17) (p44_all(i,i_angle), i_angle = 1, n_angle)

 ! extrapolate phase function at exact forward scattering angle:
 p11_all(i,0) = p11_all(i,1) - (p11_all(i,2) - p11_all(i,1))/(angle(2) - angle(1))*angle(1)

 p12_all(i,0) = p12_all(i,1) - (p12_all(i,2) - p12_all(i,1))/(angle(2) - angle(1))*angle(1)

 p22_all(i,0) = p22_all(i,1) - (p22_all(i,2) - p22_all(i,1))/(angle(2) - angle(1))*angle(1)

 p33_all(i,0) = p33_all(i,1) - (p33_all(i,2) - p33_all(i,1))/(angle(2) - angle(1))*angle(1)

 p34_all(i,0) = p34_all(i,1) - (p34_all(i,2) - p34_all(i,1))/(angle(2) - angle(1))*angle(1)

 p44_all(i,0) = p44_all(i,1) - (p44_all(i,2) - p44_all(i,1))/(angle(2) - angle(1))*angle(1)

 ! ... same for exact backscattering:
 p11_all(i,n_angle+1) = p11_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p11_all(i,n_angle)&
 - p11_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))

 p12_all(i,n_angle+1) = p12_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p12_all(i,n_angle)&
 - p12_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))

 p22_all(i,n_angle+1) = p22_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p22_all(i,n_angle)&
 - p22_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))

 p33_all(i,n_angle+1) = p33_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p33_all(i,n_angle)&
 - p33_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))

 p34_all(i,n_angle+1) = p11_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p34_all(i,n_angle)&
 - p34_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))

 p44_all(i,n_angle+1) = p44_all(i,n_angle-1) + (180d0 - angle(n_angle-1))*(p44_all(i,n_angle)&
 - p44_all(i,n_angle-1))/(angle(n_angle) - angle(n_angle-1))
end do
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Lines 426 - 433: in case the scattering angles at exact forward and backward directions are given (if not comment out)        !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
!do i=1,n_pm
! read(12) (p11_all(i,i_angle), i_angle = 0, n_angle+1)
! read(13) (p12_all(i,i_angle), i_angle = 0, n_angle+1)
! read(14) (p22_all(i,i_angle), i_angle = 0, n_angle+1)
! read(15) (p33_all(i,i_angle), i_angle = 0, n_angle+1)
! read(16) (p34_all(i,i_angle), i_angle = 0, n_angle+1)
! read(17) (p44_all(i,i_angle), i_angle = 0, n_angle+1)
!end do
!
!-------------------------------------------------------------------------------------------------------------------------------!
!
close(12);close(13);close(14)
close(15);close(16);close(17)

xmin = 1; ymin = 1; zmin = 1                          ! minimum cloud boundaries

open(7, file = 'input/'//trim(xmeshfile), status = 'old')
read(7,*) xmax

allocate( xmesh(1:xmax+1) )
do i = xmin, xmax + 1
 read(7,*) xmesh(i)
end do

close(7)

open(7, file = 'input/'//trim(ymeshfile), status = 'old')
read(7,*) ymax

allocate( ymesh(1:ymax+1) )
do i = ymin, ymax + 1
 read(7,*) ymesh(i)
end do
close(7)

open(7, file = 'input/'//trim(zmeshfile), status = 'old')
read(7,*) zmax

allocate( zmesh(0:zmax+2) )                           ! zmesh needs 0'th element in PhotonPath
do i = zmin, zmax + 1
 read(7,*) zmesh(i)
end do

close(7)

zmesh(0) = zmesh(1)

allocate( d_zmesh(zmin:zmax) )
do i = zmin, zmax
 d_zmesh(i) = zmesh(i+1)-zmesh(i)
end do

zmesh(zmax + 1) = zmesh(zmax) + d_zmesh(zmax)
zmesh(zmax + 2) = 1d100                               ! beyond cloudtop for subroutine getgridnumber

d_xmesh = xmesh(xmax+1) - xmesh(xmax)                 ! constant grid size in x
d_ymesh = ymesh(ymax+1) - ymesh(ymax)                 ! constant grid size in y

open(10, file = 'input/'//trim(cubefile)//'.ext', status = 'old', form = 'unformatted')  ! total volume extinction coefficient [m^-1 or km^-1]
open(11, file = 'input/'//trim(cubefile)//'.w0', status = 'old', form = 'unformatted')   ! total single scattering albedo
open(46, file = 'input/'//trim(cubefile)//'.sca', status = 'old', form = 'unformatted')  ! total volume scattering coefficient
open(47, file = 'input/'//trim(cubefile)//'.abs', status = 'old', form = 'unformatted')  ! total volume absorption coefficient

allocate( w_0(xmin:xmax, ymin:ymax, zmin:zmax) ); w_0  = 0d0
allocate( co_w_0(xmin:xmax, ymin:ymax, zmin:zmax) ); co_w_0 = 0d0
allocate( beta_x(0:xmax+2,0:ymax+2,0:zmax+2) )                                 ! + 2 or + 1 (it is not an issue)
allocate( tau_field(xmin:xmax, ymin:ymax) )                                    ! optical thickness
allocate( beta_abs(xmin:xmax, ymin:ymax, zmin:zmax) )
allocate( beta_sca(0:xmax+2,0:ymax+2,0:zmax+2) )                               ! + 2 or + 1 (it is not an issue)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! SOS: If beta_abs and beta_sca are not available COMMENT OUT lines 486, 487, 514, and 515. BUT, since these variables          !
! are in need, they need to be calculated. Then, go to INIT subroutine and uncomment lines 609 - 616                            !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
beta_x = 1d10   ! Termination criterion
beta_sca = 1d10
tau_mean = 0.
cloud_cover_true = 0.

! read the 3D scattering properties and calculate the z-integrated tau (optical thickness) and cloud cover
do i = xmin, xmax
 do j = ymin, ymax
  help1 = 0d0
  read(10) (beta_x(i,j,k), k = zmin, zmax )
  read(11) (w_0(i,j,k), k = zmin, zmax)
  read(46) (beta_sca(i,j,k), k = zmin, zmax)
  read(47) (beta_abs(i,j,k), k = zmin, zmax)

  do k = zmin, zmax
   tau = beta_x(i,j,k)*d_zmesh(k)
   tau_mean = tau_mean + tau
   help1 = help1 + tau
  end do

  tau_field(i,j) = help1
  if (help1 .gt. 0.01) then
    cloud_cover_true = cloud_cover_true + 1
  end if
 end do
end do

close(10);close(11)
close(46);close(47)

! open and read scattering phase matrix index
allocate( pm_index(xmin:xmax, ymin:ymax, zmin:zmax) ); pm_index = 0

open(18, file = 'input/'//trim(cubefile)//'.pm-index', status = 'old', form = 'unformatted')

do ix = xmin, xmax
 do iy = ymin, ymax
   read(18) (pm_index(ix,iy,iz), iz = zmin, zmax)
 end do
end do
close(18)

! write z-integrated tau, test reasons
open(18, file = 'output/'//trim(mainfile)//'.cloud')

do ix = xmin, xmax
 do iy = ymin, ymax
   write(18,*) tau_field(ix,iy)
 end do
end do

close(18)
deallocate(tau_field)

help1 = dble((xmax - xmin + 1))*dble(ymax - ymin + 1)
tau_mean = tau_mean/help1
cloud_cover_true = cloud_cover_true/help1
!
end subroutine INPUT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine INIT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Construct LUT: allocation of all arrays. Sets startvalues and dimensions of some variables. In other words, it calulates      !
! what is needed for for the simulations and stores them in a form (LUT)                                                        !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
call system_clock(start_time)                 ! it calculates simulation time
call date_and_time(date, time, zone, values)  ! returns data from the real-time clock and the date
                                              ! it returns 8 values that are used to create random numbers (see the loop)
allocate( theta_hauf(0:n_angle) )
theta_hauf = 0d0                              ! scattering phase function check

iflag = 0
do i = 1, 8                                   ! here I create the values used to create the random numbers
 iflag = iflag - values(i)
end do

iseed = iflag                                 ! store iflag -> reproducability
help1 = dble(ran(iflag))                      ! seed NR random number generator 'ran' (p. 1142)

! calculate pi, pi/2, 2pi, and degrees <-> radians conversions
pi = 4d0*atan(1d0); pi2 = pi/2d0; twopi = 2d0*pi
deg = 180d0/pi; rad = 1d0/deg

do ix = xmin, xmax                            ! here i calculate the absorption
 do iy = ymin, ymax
  do iz = zmin, zmax
    co_w_0(ix,iy,iz) = 1d0 - w_0(ix,iy,iz)
  end do
 end do
end do
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! In case beta_abs and beta_sca are not given, they should be calculated. Uncomment the following lines 609 - 616               !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
!do ix = xmin, xmax
!  do iy = ymin, ymax
!    do iz = zmin, zmax
!      beta_abs(ix,iy,iz) = co_w_0(ix,iy,iz) * beta_x(ix,iy,iz) ! the abrosption coeff
!      beta_sca(ix,iy,iz) = beta_x(ix,iy,iz) * w_0(ix,iy,iz)    ! the scattering coeff
!    end do
!  end do
!end do

co_albedo_lambert = 1d0 - albedo_lambert
theta_lambert_max = pi2 - 1d-2                ! mimimum theta to avoid endless horizontal flights below cloud


store_every = 100000                          ! save intermediate results every 'store_every' photons
eps = 1d-9                                    ! single prec.

delta_i(-1) = 0                               ! up or down in PhotonPath
delta_i( 1) = 1
!
! (1) preparations for cumulative scattering phase function
!
i_s = 0
angle_o(i_s) = 0.0

! new approach to create upper and lower bounds integral
do while ((angle_o(i_s) .lt. 180.0).and.(i_s.lt.n_angle))
   i_s = i_s + 1
   angle_u(i_s) = angle_o(i_s - 1)
   angle_o(i_s) = ( angle(i_s) + angle(i_s + 1) ) / 2
end do
!
! (2) solid angle interval:
!
do i = 1, n_angle
 angle_u(i) = angle_u(i)*rad; angle_o(i) = angle_o(i)*rad
 d_angle(i) = angle_o(i) - angle_u(i)                         ! in rads
 arg(i) = cos(angle_u(i)) - cos(angle_o(i))                   ! d_Omega
end do
!
! (3) cumulative p11
!
allocate( p11_cum(1:n_pm, 0:n_angle) );p11_cum=0d0
allocate(pintegral(1:n_pm));pintegral=0d0

do g = 1, n_pm
 help1 = 0d0
 p11_arg_int_sum(0) = 0d0

 do i_angle = 1, n_angle                                      ! cum. sum
   d_p11_arg_int_sum = p11_all(g,i_angle)*arg(i_angle)
   help1 = help1 + d_p11_arg_int_sum
   p11_arg_int_sum(i_angle) = help1
 end do

 do i_s = 1, d_theta_i + 1                                    ! extend pflux_int_sum d_theta_i beyond n_angle
   p11_arg_int_sum(n_angle + i_s) = p11_arg_int_sum(n_angle)
 end do

 pintegral(g) = p11_arg_int_sum(n_angle)                      ! normalization of the scattering phase function

 do i_angle = 1, n_angle
   p11_cum(g,i_angle) = p11_arg_int_sum(i_angle)/pintegral(g) ! cum. scattering phase function sum
 end do
 p11_cum(g,0) = 0d0
end do
!
! (4) construct LUT for cos and sin(scattering angle) and p11, p12, scattering angle:
!
allocate( costheta_scat(1:n_pm,0:n_rand), sintheta_scat(1:n_pm,0:n_rand))
allocate( cosphi_scat(0:n_rand), sinphi_scat(0:n_rand) )
allocate( p11_scat(1:n_pm,0:n_rand+1) )
allocate( p12_scat(1:n_pm,0:n_rand+1) )
allocate( p22_scat(1:n_pm,0:n_rand+1) )
allocate( p33_scat(1:n_pm,0:n_rand+1) )
allocate( p34_scat(1:n_pm,0:n_rand+1) )
allocate( p44_scat(1:n_pm,0:n_rand+1) )
allocate( theta_true(1:n_pm,0:n_rand) )

n_rand_p1 = n_rand + 1
d_n_rand_p1 = dble(n_rand_p1)
d_n_rand = dble(n_rand)

do g = 1, n_pm
 do i_rand = 0, n_rand
   r = dble(i_rand)/dble(n_rand) ! r in [0,1]
   theta_i = 1
   do while (p11_cum(g,theta_i) .lt. r)
     theta_i = theta_i + 1
   end do

   ! cumulative scattering phase function at "scattering angle" irand
   d_p11_arg_int_sum = p11_cum(g,theta_i) - p11_cum(g,theta_i-1)
   t = (r - p11_cum(g,theta_i - 1))/d_p11_arg_int_sum
   theta_0_true = angle_u(theta_i) + t*d_angle(theta_i) ! rad
   theta_true(g,i_rand) = theta_0_true
   costheta_scat(g,i_rand) = cos(theta_0_true)
   sintheta_scat(g,i_rand) = sin(theta_0_true)
 end do
end do

do g = 1, n_pm
 p11_scat(g,0) = p11_all(g,0)               ! no interplation here
 p11_scat(g, n_rand) = p11_all(g, n_angle)  ! "

 p12_scat(g,0) = p12_all(g,0)               ! no interplation here
 p12_scat(g, n_rand) = p12_all(g, n_angle)  ! "

 p22_scat(g,0) = p22_all(g,0)               ! no interplation here
 p22_scat(g, n_rand) = p22_all(g, n_angle)  ! "

 p33_scat(g,0) = p33_all(g,0)               ! no interplation here
 p33_scat(g, n_rand) = p33_all(g, n_angle)  ! "

 p34_scat(g,0) = p34_all(g,0)               ! no interplation here
 p34_scat(g, n_rand) = p34_all(g, n_angle)  ! "

 p44_scat(g,0) = p44_all(g,0)               ! no interplation here
 p44_scat(g, n_rand) = p44_all(g, n_angle)  ! "

 do i_rand = 1, n_rand - 1
   angle_rand = dble(i_rand)/dble(n_rand)*180d0 ! in [0,180]
   i_angle = 0                                  ! search pf scatt. angle

   do while (angle(i_angle) < angle_rand)
      i_angle = i_angle + 1
   end do

   i_angle = i_angle - 1
   t = (angle_rand - angle(i_angle))/(angle(i_angle+1) - angle(i_angle))
   p11_scat(g, i_rand) = p11_all(g,i_angle) + t*(p11_all(g,i_angle+1) - p11_all(g,i_angle))
   p12_scat(g, i_rand) = p12_all(g,i_angle) + t*(p12_all(g,i_angle+1) - p12_all(g,i_angle))
   p22_scat(g, i_rand) = p22_all(g,i_angle) + t*(p22_all(g,i_angle+1) - p22_all(g,i_angle))
   p33_scat(g, i_rand) = p33_all(g,i_angle) + t*(p33_all(g,i_angle+1) - p33_all(g,i_angle))
   p34_scat(g, i_rand) = p34_all(g,i_angle) + t*(p34_all(g,i_angle+1) - p34_all(g,i_angle))
   p44_scat(g, i_rand) = p44_all(g,i_angle) + t*(p44_all(g,i_angle+1) - p44_all(g,i_angle))

   ! here the scattering phase matrix is already normalized to 4pi
 end do
end do
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! The integral of the total scattering phase function over the unit sphere is normalized to 4pi. Actually, the calculated       !
! normalization factor (pintegral) is 2 (see equation 2.33 in Barlakas, 2016b). To avoid any problems with the given total      !
! scattering phase matrices, i.e., not correctly normalized (pintegral might slightly differ from 2, i.e., 2.0001 or so), one   !
! should renormalize the scattering phase matrix elements. How? Each phase matrix element should be multilpied by 2 and divided !
! by the pintegral (lines 761 - 770).                                                                                           !
!                                                                                                                               !
! To writer only: see Notebook5, page 116                                                                                       !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
do g = 1, n_pm
    do i_rand = 0, n_rand
        p11_scat(g, i_rand) = p11_scat(g, i_rand)*2.d0/pintegral(g)
        p12_scat(g, i_rand) = p12_scat(g, i_rand)*2.d0/pintegral(g)
        p22_scat(g, i_rand) = p22_scat(g, i_rand)*2.d0/pintegral(g)
        p33_scat(g, i_rand) = p33_scat(g, i_rand)*2.d0/pintegral(g)
        p34_scat(g, i_rand) = p34_scat(g, i_rand)*2.d0/pintegral(g)
        p44_scat(g, i_rand) = p44_scat(g, i_rand)*2.d0/pintegral(g)
    end do
end do

deallocate( p11_all,p12_all, angle_u, angle_o, p11_cum)   ! total scattering phase matrix is replaced by LUT
deallocate( p22_all,p33_all,p34_all,p44_all )

! create LUT for the scattering azimuth angle. Randomly chosen within [0,2pi]
do i_rand = 0, n_rand
 r = dble(i_rand)/dble(n_rand)
 phi_true = twopi*r
 cosphi_scat(i_rand) = cos(phi_true)
 sinphi_scat(i_rand) = sin(phi_true)
end do

! define the dimensions of arrays
allocate( flux_S(xmin:xmax, ymin:ymax) ); flux_S = 0d0
allocate( albedo(xmin:xmax, ymin:ymax) ); albedo = 0d0
allocate( flux_in(xmin:xmax, ymin:ymax) ); flux_in = 0d0
allocate( select_radiance(xmin:xmax, ymin:ymax, 1:n_radiance) ); select_radiance = 0d0
allocate( DOP(xmin:xmax, ymin:ymax, 1:n_radiance) ); DOP = 0d0
allocate( error_select_radiance(xmin:xmax, ymin:ymax, 1:n_radiance) ); error_select_radiance = 0d0
allocate( stokes_select_radiance(xmin:xmax, ymin:ymax,1:4 ,1:n_radiance) );stokes_select_radiance = 0d0
allocate( error_stokes_select_radiance(xmin:xmax, ymin:ymax, 1:4,1:n_radiance) ); error_stokes_select_radiance = 0d0
allocate( k_x_choice(1:n_radiance), k_y_choice(1:n_radiance), k_z_choice(1:n_radiance) )
allocate( absorptance(xmin:xmax, ymin:ymax) ); absorptance = 0d0
allocate( absorptance_ground(xmin:xmax, ymin:ymax) ); absorptance_ground = 0d0
allocate( transmittance_total(xmin:xmax, ymin:ymax) ); transmittance_total = 0d0
allocate( transmittance_direct(xmin:xmax, ymin:ymax) ); transmittance_direct = 0d0
allocate( transmittance_diffuse(xmin:xmax, ymin:ymax) ); transmittance_diffuse = 0d0
allocate( net_horizontal_transport(xmin:xmax, ymin:ymax) ); net_horizontal_transport = 0d0

! calc. detector direction vector/cosines (k_x,k_y,k_z)=(0,0,1)
do i = 1, n_radiance
 k_x_choice(i) = cos(phi_choice(i)*rad)*sin(theta_choice(i)*rad)
 k_y_choice(i) = sin(phi_choice(i)*rad)*sin(theta_choice(i)*rad)
 k_z_choice(i) = cos(theta_choice(i)*rad)
end do

theta_0 = theta_0*rad
costheta_0_in = cos(theta_0)
sintheta_0_in = sin(theta_0)
phi_0 = phi_0*rad
cosphi_0_in = cos(phi_0)
sinphi_0_in = sin(phi_0)

! geometrical cloud size = Number of cells -> Cell width = 1

cloudtop    = zmesh(zmax+1)
cloudbottom = zmesh(zmin)
cloudheight = cloudtop - cloudbottom
cloudleft   = xmesh(xmin)
cloudright  = xmesh(xmax + 1)
cloudwidth  = cloudright - cloudleft
cloudfront  = ymesh(ymin)
cloudback   = ymesh(ymax + 1)
clouddepth  = cloudback - cloudfront

help1 = max(cloudheight, cloudwidth, clouddepth) ! get largest cloud dimension
help2 = help1*eps                                ! ajust 'cloud-epsilon' to cloud dimension

cloudtop_in = cloudtop - help2                   ! numerically more stable than cloudtop - eps
cloudbottom_in = cloudbottom + help2
cloudleft_in = cloudleft + help2
cloudright_in = cloudright - help2
cloudfront_in = cloudfront + help2
cloudback_in = cloudback - help2
cloudtop_out = cloudtop + help2
cloudbottom_out = cloudbottom - help2
!
end subroutine INIT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine init_photon
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Photon initiation: Directions are specified by solar position. Set a scalar and a vector weight (Stokes weight) to each photon!
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! scalar weight
weight = 1.d0

! It is related to transmission. It counts the number of scattering events before the photon reaches the ground. If non direct
count = 0

! stokes weight
S_in = 0.d0
S_in(1) = 1.d0

! define the direction cosine of the photon.
k_x_new = sintheta_0_in*cosphi_0_in
k_y_new = sintheta_0_in*sinphi_0_in
k_z_new = -costheta_0_in

! for the needs of the polarization subroutine
phi_sca = phi_0
theta_sca = theta_0

! here you define the position of the photon at the top of the medium (x,y,z) and the position on the specific grid box (ix,iy,iz)
x = cloudleft + dble(ran(iflag))*cloudwidth
ix = xmin + int( (x - cloudleft)/d_xmesh)
y = cloudfront + dble(ran(iflag))*clouddepth
iy = ymin + int( (y - cloudfront)/d_ymesh)
z = cloudtop_in
iz = zmax

! scalar and vector fluxes
flux_in(ix,iy) = flux_in(ix,iy) + weight
flux_S(ix,iy) = flux_S(ix,iy) + S_in(1)
!
end subroutine init_photon
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine PhotonPath
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Photon tracing and absorption: Photons are traced throughout their propagation in a scattering and absorbing medium until a   !
! scattering event takes place; absorption calculation (see subsection 3.2.3, in Barlakas, 2016b)                               !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
sca_sum = 0d0
abs_sum = 0d0
tau_sum = 0d0

help1 = dble(ran(iflag)) ! init jump

! random path/optical thickness
if (help1 .gt. eps) then
  enough = -log(help1)
else
  enough = 1d10
end if

! until the cumulated scattering optical thickness exceeds the randomly chosen (exponentially distributed) optical thickness
do while (sca_sum .lt. enough)
 !
 ! (1) get local volume extinction coeff.
 !
 local_beta_x = beta_x(ix,iy,iz)
 local_beta_sca = beta_sca(ix,iy,iz)
 local_beta_abs = beta_abs(ix,iy,iz)
 !
 ! (2) Trace the photon from current grid box surface to the intersection with the nearest neighbor grid–box surface
 !
 t_ray = 1d10
 if (abs(k_x_new) .gt. eps) then
  up_or_down_x = int(sign(1.0d0, k_x_new))
  t_ray(1) = ( xmesh(ix + delta_i(up_or_down_x)) - x )/k_x_new
 end if
 if (abs(k_y_new) .gt. eps) then
  up_or_down_y = int(sign(1.0d0, k_y_new))
  t_ray(2) = ( ymesh(iy + delta_i(up_or_down_y)) - y )/k_y_new
 end if
 if (abs(k_z_new) .gt. eps) then
  up_or_down_z = int(sign(1.0d0, k_z_new))
  t_ray(3) = ( zmesh(iz + delta_i(up_or_down_z)) - z )/k_z_new
 end if

 t_ray_min = minval(t_ray)              ! gives the minimum
 i_ray_min = minloc(t_ray, dim = 1)     ! gives the place of the minimum

 ! shortest distance is realistic. Intersection with the smallest distance
 tau_sum = tau_sum + t_ray_min*local_beta_x
 sca_sum = sca_sum + t_ray_min*local_beta_sca
 abs_sum = abs_sum + t_ray_min*local_beta_abs

 ! make sure new photon is inside next box
 help1 = (t_ray_min + eps)
 x = x + help1*k_x_new
 y = y + help1*k_y_new
 z = z + help1*k_z_new

 ! periodic bounds only if tau not exceeded
 if (sca_sum .lt. enough) then

  wasnun:   select case (i_ray_min)

  case(1)
   ihelp1 = delta_i(up_or_down_x) ! (-1,+1) -> (0,1)
   ix = ix + up_or_down_x
   if (ix .lt. xmin) then
    ix = xmax; x = cloudright_in
   end if
   if (ix .gt. xmax) then
    ix = xmin; x = cloudleft_in
   end if

  case(2)
   ihelp1 = delta_i(up_or_down_y) ! (-1,+1) -> (0,1)
   iy = iy + up_or_down_y
   if (iy .gt. ymax) then
    iy = ymin; y = cloudfront_in
   end if
   if (iy .lt. ymin) then
    iy = ymax; y = cloudback_in
   end if

  case(3)
   ihelp1 = delta_i(up_or_down_z) ! (-1,+1) -> (0,1)
   iz = iz + up_or_down_z

   if (iz .lt. zmin) then
    sca_sum = enough; z = cloudbottom_out
   end if
   if (iz .gt. zmax) then
    sca_sum = enough; z = cloudtop_out
   end if

  end select wasnun
 end if ! sca_sum < enough
end do ! while
!
! (3) Subsequently, the photon steps backward to ensure that the total photon path exactly matches randomly chosen one.
!     Finally, calculate the absorbing optical thickness (correction for overjumps)
!
if (local_beta_sca .gt. 0d0) then
  step_back = (sca_sum - enough)/local_beta_sca
  z = z - step_back*k_z_new
  y = y - step_back*k_y_new
  x = x - step_back*k_x_new
  !
  !-----------------------------------------------------------------------------------------------------------------------------!
  !                                                                                                                             !
  ! BUG report: For very high tau at the edges of the domain (and a large number of photons), it might occur that the above     !
  ! corrections for overjumps is not sufficient, meaning the photon might be slightly outside the boundaries. In such case, we  !
  ! make sure that the photon is inside the closest grid box                                                                    !
  !                                                                                                                             !
  !-----------------------------------------------------------------------------------------------------------------------------!
  !
  if (x.gt.xmesh(xmax+1)) then
    x = cloudright_in
  end if
  if (y.gt.ymesh(ymax+1)) then
    y = cloudback_in
  end if

  ! find the new grid box that the photon is located
  call getgridnumber
  abs_sum = abs_sum - step_back*local_beta_abs
end if

! (4) Finally, update the Stokes and the scalar weights: attenuated by exp(-abs_sum) along the photon path

S_in(1) = S_in(1) * exp(-abs_sum)
S_in(2) = S_in(2) * exp(-abs_sum)
S_in(3) = S_in(3) * exp(-abs_sum)
S_in(4) = S_in(4) * exp(-abs_sum)
weight = weight * exp(-abs_sum)
!
end subroutine PhotonPath ! photon path through inhom. cloud
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine getgridnumber
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Gives the grid number of the specific grid box that the photon is located.                                                    !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
ix = xmin + int( (x - cloudleft)/d_xmesh)
iy = ymin + int( (y - cloudfront)/d_ymesh)
iz = zmin
zpos = cloudbottom

do while(zpos .lt. z)
 iz = iz + 1
 zpos = zmesh(iz)
end do

iz = iz - 1

end subroutine getgridnumber
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine polarization
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Anisotropic scattering and Rotation of the Stokes Vector: A scattering event takes place. The importance sampling method has  !
! been employed to sample the new direction (see Barlakas et al., 2016a and subsection 3.2.4 in Barlakas, 2016b).               !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! update the counters and absorptance
N_scat = N_scat + 1
count = count + 1
absorptance(ix,iy) = absorptance(ix,iy) + weight*co_w_0(ix,iy,iz)    ! ix from PhotonPath
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! As explained in init_photon subroutine, we assigned as k_*_new the initial direction cosines (the same stands for the         !
! zenithal and azimuth angles); this is done only once before simulation. Whenever a scattering event takes place, the          !
! scattering direction of the previous scattering event is the incident direction for the next scattering event                 !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
k_x = k_x_new
k_y = k_y_new
k_z = k_z_new
phi_inc = phi_sca
theta_inc = theta_sca
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Mishchenko Michael suggestions:                                                                                               !
! 1. If theta_inc < 1D-7, make theta_inc = 1D-7.                                                                                !
! 2. If pi - theta_inc < 1D-7, make theta_inc = pi - 1D-7.                                                                      !
! 3. If theta_sca < 1D-7, make theta_sca = 1D-7.                                                                                !
! 4. If pi - theta_sca < 1D-7, make theta_sca = pi - 1D-7.                                                                      !
! Steps 1-4 help avoid cases of theta angles being too close to zero or pi. The use of double-precision arithmetic should take  !
! care of near 0/0 cases                                                                                                        !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Suggestions 1. and 2.
if (theta_inc.lt.1.d-7) theta_inc = 1.d-7
if ((pi - theta_inc).lt.1.d-7) theta_inc = pi - 1.d-7

! scattering phase function index
i_pm = pm_index(ix,iy,iz)

! for each scattering event we calculate its contribution to the radiance calculation (check the corresponding subroutine)
call LocalEstimateMethod
!
! (1) Scattering event: determine first the scattering zenithal angle
!
! throw the dice (random number [0,1])
r =   dble(ran(iflag))

! convert r to the size of the LUT [0,10000]
help1 = d_n_rand*r

! go to the LUT and get the corresponding scattering angle (interpolation is needed)
i_rand = int(help1)
t = help1 - dble(i_rand)
costheta_0 = costheta_scat(i_pm,i_rand) + t*(costheta_scat(i_pm,i_rand+1) - costheta_scat(i_pm,i_rand))
sintheta_0 = sintheta_scat(i_pm,i_rand) + t*(sintheta_scat(i_pm,i_rand+1) - sintheta_scat(i_pm,i_rand))

!###############################################
!                                              !
! MAJOR BUG: obtain the correct i_rand
!
help1 = d_n_rand*acos(costheta_0)/pi
i_rand = int(help1)
t = help1 - dble(i_rand)
!                                              !
!###############################################
!
! (2) interpolation for scattering phase matrix
!
p11_inter = p11_scat(i_pm,i_rand)+ t*(p11_scat(i_pm,i_rand+1) - p11_scat(i_pm,i_rand))
p12_inter = p12_scat(i_pm,i_rand)+ t*(p12_scat(i_pm,i_rand+1) - p12_scat(i_pm,i_rand))
p22_inter = p22_scat(i_pm,i_rand)+ t*(p22_scat(i_pm,i_rand+1) - p22_scat(i_pm,i_rand))
p33_inter = p33_scat(i_pm,i_rand)+ t*(p33_scat(i_pm,i_rand+1) - p33_scat(i_pm,i_rand))
p34_inter = p34_scat(i_pm,i_rand)+ t*(p34_scat(i_pm,i_rand+1) - p34_scat(i_pm,i_rand))
p44_inter = p44_scat(i_pm,i_rand)+ t*(p44_scat(i_pm,i_rand+1) - p44_scat(i_pm,i_rand))
!
! (2.1) Importance sampling method: division by P11 corrects for using P11 as a PDF to obtain scattering theta
!
p12_inter = p12_inter/p11_inter
p22_inter = p22_inter/p11_inter
p33_inter = p33_inter/p11_inter
p34_inter = p34_inter/p11_inter
p44_inter = p44_inter/p11_inter
p11_inter = p11_inter/p11_inter
!
! (3) Determine the scattering azimuth angle
!
r =  dble(ran(iflag)) ! through the dice
help1 = d_n_rand*r    ! from [0,1] to [0,10000]
i_rand = int(help1)   ! i_rand in [0,n_rand-1]

! interpolation:
t = help1 - dble(i_rand)
cosphi_0 = cosphi_scat(i_rand) + t*(cosphi_scat(i_rand+1) - cosphi_scat(i_rand))
sinphi_0 = sinphi_scat(i_rand) + t*(sinphi_scat(i_rand+1) - sinphi_scat(i_rand))

! recreate P11 check weahter everything is under control (the user might want to comment this out)
theta_0_true = acos(costheta_0)
ihelp1 = idnint(theta_0_true*deg*real(n_angle)/180.)
theta_hauf(ihelp1) = theta_hauf(ihelp1) + 1.0
!
! (4) determine the new direction cosines
!
test1 = 1d0 - k_z*k_z                                      ! sin^2(theta_old)

! having defined the scattering event (zenithal, azimuth) the new direction of the photon is calculated

if (test1 .gt. 0) then                                     ! old photon direction not exactly along z-axis
  help1 = sqrt((1d0 - costheta_0*costheta_0)/test1)
  help2 = k_z*cosphi_0

  k_x_new = k_x*costheta_0 - (k_y*sinphi_0 + k_x*help2)*help1
  k_y_new = k_y*costheta_0 + (k_x*sinphi_0 - k_y*help2)*help1
  k_z_new = k_z*costheta_0 + test1*cosphi_0*help1
else
  k_x_new = sintheta_0*cosphi_0
  k_y_new = sintheta_0*sinphi_0
  k_z_new = k_z*costheta_0
end if

! normalizations
slen = k_x_new**2 + k_y_new**2 + k_z_new**2
k_x_new = k_x_new/dsqrt(slen)
k_y_new = k_y_new/dsqrt(slen)
k_z_new = k_z_new/dsqrt(slen)
!
! (5) calculate the zenithal and azimuth angles of the new direction
!
phi_sca = 0.5*twopi-atan2(k_y_new,-k_x_new)
sintheta_s = k_x_new/cos(phi_sca)
theta_sca = 0.5*twopi-atan2(sintheta_s,-k_z_new)
!
! (5.1) get rid of problematic cases (Suggestions 3. and 4.)
!
if (theta_sca.lt.1.d-7) then
    theta_sca = 1.d-7
    k_z_new = cos(theta_sca)
end if
if ((pi - theta_sca).lt.1.d-7) then
    theta_sca = pi - 1.d-7
    k_z_new = cos(theta_sca)
end if
!
! (6) define the relative azimuth angle and make sure it is within the interval [0,2pi]
!
d_phi = phi_sca - phi_inc

if (d_phi.lt.0) d_phi = d_phi + twopi
!
! (7) define the rotation angles (spherical trigonometry) and make sure there is no arithmetic problem
!
help3_1 = sin(theta_inc)*sintheta_0

if (dabs(dabs(costheta_0)-1.d0).lt.1.0E-10) then
  cosi = costheta_0/dabs(costheta_0)
else
  cosi = (k_z_new - k_z*costheta_0)/help3_1
end if

if (cosi.gt.1.d0) cosi=1.d0
if (cosi.lt.(-1.d0)) cosi=-1.d0

cos2i = 2.d0*cosi*cosi - 1.d0
sin2i = 2.d0*dsqrt(1.d0 - cosi*cosi)*cosi

if (sin2i.eq.0) sin2i = 0.d0
if (dabs(cos2i).gt.1.d0) then
if (dabs(cos2i).gt.1.000001) write(*,*) 'cos2i>1',cos2i
   cos2i = cos2i/dabs(cos2i)
   sin2i = 0.d0
end if

help3_2 = sintheta_s*sintheta_0

if (dabs(dabs(costheta_0)-1.d0).lt.1.0E-10) then
  cosj = costheta_0/dabs(costheta_0)
else
  cosj = (k_z - k_z_new*costheta_0)/help3_2
end if

if (cosj.gt.1.d0) cosj=1.d0
if (cosj.lt.(-1.d0)) cosj=-1.d0

cos2j = 2.d0*cosj*cosj-1.d0
sin2j = 2.d0*dsqrt(1 - cosj*cosj)*cosj

if (sin2j.eq.0) sin2j = 0.d0
if (dabs(cos2j).gt.1.d0) then
if (dabs(cos2j).gt.1.000001) write(*,*) 'cos2j>1',cos2j
   cos2j = cos2j/dabs(cos2j)
   sin2j = 0.d0
end if
!
! (8) perform the rotation of the Stokes vector depending on the relative azimuth angle
!
if ((d_phi.gt.0).and.(d_phi.lt.pi)) then  ! d_phi > 0

 S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) + p12_inter*sin2i*S_in(3)

 S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&
 + (p22_inter*sin2i*cos2j + p33_inter*cos2i*sin2j)*S_in(3) + p34_inter*sin2j*S_in(4)

 S_out(3)= - p12_inter*sin2j*S_in(1) + ( - p22_inter*cos2i*sin2j - p33_inter*sin2i*cos2j)*S_in(2)&
 + ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)

 S_out(4) = + p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)

else if ((d_phi.gt.pi).and.(d_phi.lt.twopi)) then

  S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) - p12_inter*sin2i*S_in(3)

  S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&
  + ( - p22_inter*sin2i*cos2j - p33_inter*cos2i*sin2j)*S_in(3) - p34_inter*sin2j*S_in(4)

  S_out(3) = p12_inter*sin2j*S_in(1) + (p22_inter*cos2i*sin2j + p33_inter*sin2i*cos2j)*S_in(2)&
  + ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)

  S_out(4) = - p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)

else

 S_out(1) =   p11_inter * S_in(1) + p12_inter * S_in(2)
 S_out(2) =   p12_inter * S_in(1) + p22_inter * S_in(2)
 S_out(3) =   p33_inter * S_in(3) + p34_inter * S_in(4)
 S_out(4) = - p34_inter * S_in(3) + p44_inter * S_in(4)

end if
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Alternatively: follow Mishchenko's suggestion: always perform rotations; You still need to make the second rotation since it  !
! affects the 12 and 13 elements of the resulting scattering phase matrix.                                                      !
!                                                                                                                               !
! Both approaches have been tested: insignificant differences, the user should decide                                           !
!                                                                                                                               !
! Mishchenko's suggestion:                                                                                                      !
!                                                                                                                               !
! if ((d_phi.ge.0).and.(d_phi.le.pi)) then                                                                                      !
!                                                                                                                               !
! S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) + p12_inter*sin2i*S_in(3)                                              !
! S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&                                 !
! + (p22_inter*sin2i*cos2j + p33_inter*cos2i*sin2j)*S_in(3) + p34_inter*sin2j*S_in(4)                                           !
! S_out(3)= - p12_inter*sin2j*S_in(1) + ( - p22_inter*cos2i*sin2j - p33_inter*sin2i*cos2j)*S_in(2)&                             !
! + ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)                                        !
! S_out(4) = + p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)                                            !
!                                                                                                                               !
! else if ((d_phi.gt.pi).and.(d_phi.lt.twopi)) then                                                                             !
!                                                                                                                               !
! S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) - p12_inter*sin2i*S_in(3)                                              !
! S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&                                 !
! + ( - p22_inter*sin2i*cos2j - p33_inter*cos2i*sin2j)*S_in(3) - p34_inter*sin2j*S_in(4)                                        !
! S_out(3) = p12_inter*sin2j*S_in(1) + (p22_inter*cos2i*sin2j + p33_inter*sin2i*cos2j)*S_in(2)&                                 !
! + ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)                                        !
! S_out(4) = - p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)                                            !
!                                                                                                                               !
! end if                                                                                                                        !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! (9) Update the Stokes weight: the scattered Stokes vector is the incident for the next scattering event
!
S_in(1) = S_out(1)
S_in(2) = S_out(2)
S_in(3) = S_out(3)
S_in(4) = S_out(4)
!
end subroutine polarization
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
function ran(idum)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Random number generator: from Numerical Recipes Fortran 90 : slightly modified                                                !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
implicit none
integer, parameter :: k4b = selected_int_kind(9)
integer(k4b), intent(inout) :: idum
real :: ran
integer(k4b), parameter :: ia = 16807, im = 2147483647, iq = 127773, ir = 2836
real, save :: am
integer(k4b), save :: ix = -1, iy = -1, k

if (idum <= 0 .or. iy < 0) then
 am = nearest(1.0, -1.0)/im
 iy = ior(ieor(888889999, abs(idum)),1)
 ix = ieor(777755555, abs(idum))
 idum = abs(idum) + 1
end if
ix = ieor(ix, ishft(ix, 13))
ix = ieor(ix, ishft(ix, -17))
ix = ieor(ix, ishft(ix, 5))
k = iy/iq
iy = ia*(iy - k*iq) - ir*k

if (iy < 0) iy = iy + im
ran = am*ior(iand(im, ieor(ix,iy)),1)
!
end function ran
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
!subroutine rotation(S,i,F,Q,U,V)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Rotation of the Stokes vector: it applies the first operation, meaning rotate the stokes vector so that is expressed with     !
! respect to the scattering plane (not used anymore since rotations are made on the "flow")                                     !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
!
!implicit none
!real(8),dimension(1:4), intent(in)  :: S
!real(8), intent(in)                 :: i
!real(8), intent(out)                :: F,Q,U,V
!real(8)                             :: cos2i,sin2i
!
!cos2i=cos(2*i)
!sin2i=sin(2*i)
!
!if (dabs(cos2i).gt.1.d0) then
!if (dabs(cos2i).gt.1.000001) write(*,*) 'cos2i>1',cos2i
!  cos2i = cos2i/dabs(cos2i)
!  sin2i = 0
!end if
!
!F = S(1)                         ! here is R(-i)
!Q = S(2)*cos2i-S(3)*sin2i
!U = S(2)*sin2i+S(3)*cos2i
!V = S(4)
!
!end subroutine rotation
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine LocalEstimateMethod
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Local Estimate Method for efficient radiance calculation: it accounts for the probability that the photon is scattered into   !
! the direction of the sensor at each scattering process, always considering the attenuation tau along the photon path (is      !
! performed in PhotonPathDetector subroutine)                                                                                   !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
do i_radiance = 1, n_radiance

 phi_det = phi_choice(i_radiance)*rad
 sintheta_LE = k_x_choice(i_radiance)/cos(phi_det)
 theta_det = theta_choice(i_radiance)*rad

 ! Mishchenko's suggestions (the same are implemented for the detector angles)
 if (theta_det.lt.1.d-7) then
   theta_det = 1.d-7
   k_x_choice(i_radiance) = cos(phi_det)*sin(theta_det)
   k_y_choice(i_radiance) = sin(phi_det)*sin(theta_det)
   k_z_choice(i_radiance) = cos(theta_det)
 end if
 if ((pi - theta_det).lt.1.d-7) then
   theta_det = pi - 1.d-7
   k_x_choice(i_radiance) = cos(phi_det)*sin(theta_det)
   k_y_choice(i_radiance) = sin(phi_det)*sin(theta_det)
   k_z_choice(i_radiance) = cos(theta_det)
 end if

 !
 ! (1) define the relative azimuth angle and make sure it is within the interval [0,2pi]
 !
 d_phi_LE = phi_det - phi_sca
 if (d_phi_LE.lt.0) d_phi_LE = d_phi_LE + twopi
 !
 ! (2) define the angle ("scattering") between the incident direction and the direction of observation
 !
 r = k_x_new*k_x_choice(i_radiance) + k_y_new*k_y_choice(i_radiance) + k_z_new*k_z_choice(i_radiance) ! cos(theta_detector)
 cosnew = cos(acos(r))
 help1 = d_n_rand*acos(r)/pi ! in [0, d_n_rand) from LUT
 i_rand = int(help1)
 t = help1 - dble(i_rand)
 !
 ! (3) scattering phase matrix interpolation
 !
 p11_inter = p11_scat(i_pm,i_rand)+ t*(p11_scat(i_pm,i_rand+1) - p11_scat(i_pm,i_rand))
 p12_inter = p12_scat(i_pm,i_rand)+ t*(p12_scat(i_pm,i_rand+1) - p12_scat(i_pm,i_rand))
 p22_inter = p22_scat(i_pm,i_rand)+ t*(p22_scat(i_pm,i_rand+1) - p22_scat(i_pm,i_rand))
 p33_inter = p33_scat(i_pm,i_rand)+ t*(p33_scat(i_pm,i_rand+1) - p33_scat(i_pm,i_rand))
 p34_inter = p34_scat(i_pm,i_rand)+ t*(p34_scat(i_pm,i_rand+1) - p34_scat(i_pm,i_rand))
 p44_inter = p44_scat(i_pm,i_rand)+ t*(p44_scat(i_pm,i_rand+1) - p44_scat(i_pm,i_rand))
 !
 ! (4) division by 4pi accounts for the normalization within the solid angle (scattering probability density)
 !
 p11_inter = p11_inter/(4d0*pi)
 p12_inter = p12_inter/(4d0*pi)
 p22_inter = p22_inter/(4d0*pi)
 p33_inter = p33_inter/(4d0*pi)
 p34_inter = p34_inter/(4d0*pi)
 p44_inter = p44_inter/(4d0*pi)
 !
 ! (5) define the rotation angles (spherical trigonometry) and make sure there is no arithmetic problem
 !
 help4a = dsqrt(1-cosnew*cosnew)
 help4b = dsqrt(1-k_z_new*k_z_new)
 help4c = dsqrt(1-k_z_choice(i_radiance)*k_z_choice(i_radiance))

 help4_1 = help4a*help4b
 help4_2 = help4a*help4c

 help4_1 = sin(theta_inc)*dsqrt(1-cosnew*cosnew) ! its exact the same as 2 lines above (checking reasons: you can comment it out)

 if (dabs(dabs(cosnew)-1.d0).lt.1.0E-10) then
   cosk = cosnew/dabs(cosnew)
 else
   cosk = (k_z_choice(i_radiance) - k_z_new*cosnew)/help4_1
 end if

 if (cosk.gt.1.d0) cosk=1.d0
 if (cosk.lt.(-1.d0)) cosk=-1.d0

 cos2k = 2.d0*cosk*cosk - 1.d0
 sin2k = 2.d0*dsqrt(1.d0 - cosk*cosk)*cosk

 if (sin2k.eq.0) sin2k = 0.d0
 if (dabs(cos2k).gt.1.d0) then
 if (dabs(cos2k).gt.1.000001) write(*,*) 'cos2k>1',cos2k
   cos2k = cos2k/dabs(cos2k)
   sin2k = 0.d0
 end if

 help4_2 = sintheta_LE*dsqrt(1-cosnew*cosnew)

 if (dabs(dabs(cosnew)-1.d0).lt.1.0E-10) then
   cosh = cosnew/dabs(cosnew)
 else
   cosh = (k_z_new - k_z_choice(i_radiance)*cosnew)/help4_2
 end if

 if (cosh.gt.1.d0) cosh=1.d0
 if (cosh.lt.(-1.d0)) cosh=-1.d0

 cos2h = 2.d0*cosh*cosh - 1.d0
 sin2h = 2.d0*dsqrt(1.d0 - cosh*cosh)*cosh

 if (sin2h.eq.0.d0) sin2h = 0.d0 ! sometimes i get - 0 ... and with this its fixed
 if (dabs(cos2h).gt.1.d0) then
 if (dabs(cos2h).gt.1.000001) write(*,*) 'cos2h>1',cos2h
   cos2h = cos2h/dabs(cos2h)
   sin2h = 0.d0
 end if
 !
 ! (6) perform the rotation of the Stokes vector depending on the relative azimuth angle
 !
 if ((d_phi_LE.gt.0).and.(d_phi_LE.lt.pi)) then  ! d_phi > 0

  S_out(1) = p11_inter*S_in(1) + p12_inter*cos2k*S_in(2) + p12_inter*sin2k*S_in(3)

  S_out(2) = p12_inter*cos2h*S_in(1) + (p22_inter*cos2k*cos2h - p33_inter*sin2k*sin2h)*S_in(2)&
  + (p22_inter*sin2k*cos2h + p33_inter*cos2k*sin2h)*S_in(3) + p34_inter*sin2h*S_in(4)

  S_out(3)= - p12_inter*sin2h*S_in(1) + ( - p22_inter*cos2k*sin2h - p33_inter*sin2k*cos2h)*S_in(2)&
  + ( - p22_inter*sin2k*sin2h + p33_inter*cos2k*cos2h)*S_in(3) + p34_inter*cos2h*S_in(4)

  S_out(4) = + p34_inter*sin2k*S_in(2) - p34_inter*cos2k*S_in(3) + p44_inter*S_in(4)

 else if ((d_phi_LE.gt.pi).and.(d_phi_LE.lt.twopi)) then

  S_out(1) = p11_inter*S_in(1) + p12_inter*cos2k*S_in(2) - p12_inter*sin2k*S_in(3)

  S_out(2) = p12_inter*cos2h*S_in(1) + (p22_inter*cos2k*cos2h - p33_inter*sin2k*sin2h)*S_in(2)&
  + ( - p22_inter*sin2k*cos2h - p33_inter*cos2k*sin2h)*S_in(3) - p34_inter*sin2h*S_in(4)

  S_out(3) = p12_inter*sin2h*S_in(1) + (p22_inter*cos2k*sin2h + p33_inter*sin2k*cos2h)*S_in(2)&
  + ( - p22_inter*sin2k*sin2h + p33_inter*cos2k*cos2h)*S_in(3) + p34_inter*cos2h*S_in(4)

  S_out(4) = - p34_inter*sin2k*S_in(2) - p34_inter*cos2k*S_in(3) + p44_inter*S_in(4)

else

 S_out(1) =   p11_inter * S_in(1) + p12_inter * S_in(2)
 S_out(2) =   p12_inter * S_in(1) + p22_inter * S_in(2)
 S_out(3) =   p33_inter * S_in(3) + p34_inter * S_in(4)
 S_out(4) = - p34_inter * S_in(3) + p44_inter * S_in(4)

end if
 !
 !-------------------------------------------------------------------------------------------------------------------------------!
 !                                                                                                                               !
 ! Alternatively: follow Mishchenko's suggestion: always perform rotations; the same as in polarization subroutine               !
 !                                                                                                                               !
 !if ((d_phi.ge.0).and.(d_phi.le.pi)) then                                                                                       !
 !-------------------------------------------------------------------------------------------------------------------------------!
 !S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) + p12_inter*sin2i*S_in(3)                                               !
 !S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&                                  !
 !+ (p22_inter*sin2i*cos2j + p33_inter*cos2i*sin2j)*S_in(3) + p34_inter*sin2j*S_in(4)                                            !
 !S_out(3)= - p12_inter*sin2j*S_in(1) + ( - p22_inter*cos2i*sin2j - p33_inter*sin2i*cos2j)*S_in(2)&                              !
 !+ ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)                                         !
 !S_out(4) = + p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)                                             !
 !                                                                                                                               !
 !else if ((d_phi.gt.pi).and.(d_phi.lt.twopi)) then                                                                              !
 !                                                                                                                               !
 !S_out(1) = p11_inter*S_in(1) + p12_inter*cos2i*S_in(2) - p12_inter*sin2i*S_in(3)                                               !
 !S_out(2) = p12_inter*cos2j*S_in(1) + (p22_inter*cos2i*cos2j - p33_inter*sin2i*sin2j)*S_in(2)&                                  !
 !+ ( - p22_inter*sin2i*cos2j - p33_inter*cos2i*sin2j)*S_in(3) - p34_inter*sin2j*S_in(4)                                         !
 !S_out(3) = p12_inter*sin2j*S_in(1) + (p22_inter*cos2i*sin2j + p33_inter*sin2i*cos2j)*S_in(2)&                                  !
 !+ ( - p22_inter*sin2i*sin2j + p33_inter*cos2i*cos2j)*S_in(3) + p34_inter*cos2j*S_in(4)                                         !
 !S_out(4) = - p34_inter*sin2i*S_in(2) - p34_inter*cos2i*S_in(3) + p44_inter*S_in(4)                                             !
 !                                                                                                                               !
 !end if                                                                                                                         !
 !                                                                                                                               !
 !-------------------------------------------------------------------------------------------------------------------------------!
 !
 ! scalar weight
 radiance_weight = weight*p11_inter

 call PhotonPathDetector ! track the photon along its path form the scattering event to the detector position

end do ! do i_radiance
!
end subroutine LocalEstimateMethod
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine PhotonPathDetector
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Same as in PhotonPath: now the photon is traced along its path towards the detector (it is called within LocalEstimateMethod).!
! Here, the attenuation along the aforementioned path is calculated and taken into account in the radiance calculation. An      !
! additional difference to the PhotonPath is that the calculations are based on the extinction and not on the scattering coeff. !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! Assign detector's position and grid indices
tau_sum = 0d0
x2detector = x
y2detector = y
z2detector = z
ix2detector = ix
iy2detector = iy
iz2detector = iz

do while (z2detector .lt. cloudtop .and. z2detector .gt. cloudbottom) ! trace out of model domain
 !
 ! (1) get local volume extinction coeff.
 !
 local_beta_x = beta_x(ix2detector,iy2detector,iz2detector)
 !
 ! (2) Trace the photon from current grid box surface to the intersection with the nearest neighbor grid–box surface
 !
 t_ray = 1d10
 if (abs(k_x_choice(i_radiance)) .gt. eps) then
  up_or_down_x = int(sign(1.0d0, k_x_choice(i_radiance)))
  t_ray(1) = ( xmesh(ix2detector + delta_i(up_or_down_x)) - x2detector )/k_x_choice(i_radiance)
 end if
 if (abs(k_y_choice(i_radiance)) .gt. eps) then
  up_or_down_y = int(sign(1.0d0, k_y_choice(i_radiance)))
  t_ray(2) = ( ymesh(iy2detector + delta_i(up_or_down_y)) - y2detector )/k_y_choice(i_radiance)
 end if
 if (abs(k_z_choice(i_radiance)) .gt. eps) then
  up_or_down_z = int(sign(1.0d0, k_z_choice(i_radiance)))
  t_ray(3) = ( zmesh(iz2detector + delta_i(up_or_down_z)) - z2detector )/k_z_choice(i_radiance)
 end if

 t_ray_min = minval(t_ray)                              ! gives the minimum
 i_ray_min = minloc(t_ray, dim = 1)                     ! gives the place of the minimum
 tau_sum = tau_sum + t_ray_min*local_beta_x             ! shortest distance is realistic. Intersection with the smallest distance

 help1 = (t_ray_min + eps)                              ! make sure new photon is inside next box
 x2detector = x2detector + help1*k_x_choice(i_radiance)
 y2detector = y2detector + help1*k_y_choice(i_radiance)
 z2detector = z2detector + help1*k_z_choice(i_radiance)

 wasnun:   select case (i_ray_min)
 case(1)
 ihelp1 = delta_i(up_or_down_x) ! (-1,+1) -> (0,1)
 ix2detector = ix2detector + up_or_down_x
 if (ix2detector .lt. xmin) then
  ix2detector = xmax; x2detector = cloudright_in
 end if
 if (ix2detector .gt. xmax) then
  ix2detector = xmin; x2detector = cloudleft_in
 end if

 case(2)
 ihelp1 = delta_i(up_or_down_y) ! (-1,+1) -> (0,1)
 iy2detector = iy2detector + up_or_down_y
 if (iy2detector .gt. ymax) then
  iy2detector = ymin; y2detector = cloudfront_in
 end if
 if (iy2detector .lt. ymin) then
  iy2detector = ymax; y2detector = cloudback_in
 end if

 case(3)
 ihelp1 = delta_i(up_or_down_z) ! (-1,+1) -> (0,1)
 iz2detector = iz2detector + up_or_down_z
 end select wasnun
end do ! while
!
! (3) Radiance contribution of each scattering event
!
! vector calculations
stokes_select_radiance(ix2detector,iy2detector, 1,i_radiance) = &
stokes_select_radiance(ix2detector,iy2detector, 1, i_radiance) + S_out(1)*exp(-tau_sum)

stokes_select_radiance(ix2detector,iy2detector, 2,i_radiance) = &
stokes_select_radiance(ix2detector,iy2detector, 2, i_radiance) + S_out(2)*exp(-tau_sum)

stokes_select_radiance(ix2detector,iy2detector, 3,i_radiance) = &
stokes_select_radiance(ix2detector,iy2detector, 3, i_radiance) + S_out(3)*exp(-tau_sum)

stokes_select_radiance(ix2detector,iy2detector, 4,i_radiance) = &
stokes_select_radiance(ix2detector,iy2detector, 4, i_radiance) + S_out(4)*exp(-tau_sum)

! scalar calculations
select_radiance(ix2detector,iy2detector,i_radiance) = &
select_radiance(ix2detector,iy2detector,i_radiance) + radiance_weight*exp(-tau_sum)
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! In the following lines the statistical Monte Carlo error is calculated. Note that this error for both the scalar and vector   !
! (all the Stokes components) calculations is identical for the same pixel (for different pixels, of course, is different).     !
! Bottom line, the user might need to keep only one calculation (i.e., scalar) and comment out the rest.                        !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! calculate the statistical Monte Carlo error
!
! vector
error_stokes_select_radiance(ix2detector,iy2detector,1,i_radiance) = &
error_stokes_select_radiance(ix2detector,iy2detector,1,i_radiance) + 1d0

error_stokes_select_radiance(ix2detector,iy2detector,2,i_radiance) = &
error_stokes_select_radiance(ix2detector,iy2detector,2,i_radiance) + 1d0

error_stokes_select_radiance(ix2detector,iy2detector,3,i_radiance) = &
error_stokes_select_radiance(ix2detector,iy2detector,3,i_radiance) + 1d0

error_stokes_select_radiance(ix2detector,iy2detector,4,i_radiance) = &
error_stokes_select_radiance(ix2detector,iy2detector,4,i_radiance) + 1d0

! scalar
error_select_radiance(ix2detector,iy2detector,i_radiance) =&
error_select_radiance(ix2detector,iy2detector,i_radiance) + 1d0
!
end subroutine PhotonPathDetector
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine reflection
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Surface contribution (isotropic surface - Lambertian surface): it is driven by the Lambert's law, diffusely reflected light   !
! is isotropic and unpolarized, independently of the state of polarization of the incident radiation and the angle of           !
! polarization. In other words, reflection at the surface fully depolarizes the incident radiation. How it is interpreted? First!
! a random zenithal angle is created (global zenith). This has to be very small, smaller than the theta_lambert_max, which is   !
! the maximum zenithal angle to avoid enless horizontal flights below the cloud. Having defined theta and phi, the direction    !
! cosine of the photon is calculated; the contribution of the reflected photon to radiance calculations is accounted            !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
! update the counter
N_refl = N_refl + 1

! throw the dice to obtain the "global zenith"
zenith = acos(sqrt(dble(ran(iflag))))

! avoid horizontal flights below cloud
if (zenith .gt. theta_lambert_max) zenith = theta_lambert_max

costheta_0 = cos(zenith)
sintheta_0 = dsqrt(1d0 - costheta_0*costheta_0)

! define the azimuth angle of reflection, randomly chosen from the LUT
r =  dble(ran(iflag))
help1 = d_n_rand*r
i_rand = int(help1)

! interpolation
t = help1 - dble(i_rand)
cosphi_0 = cosphi_scat(i_rand) + t*(cosphi_scat(i_rand+1) - cosphi_scat(i_rand))
sinphi_0 = sinphi_scat(i_rand) + t*(sinphi_scat(i_rand+1) - sinphi_scat(i_rand))

! defien the new direction cosines of the reflected radiation
k_x_new = cosphi_0*sintheta_0
k_y_new = sinphi_0*sintheta_0
k_z_new = costheta_0

! normalization to make sure that the direction cosines are unit vectors
slen = k_x_new**2 + k_y_new**2 + k_z_new**2
k_x_new = k_x_new/dsqrt(slen)
k_y_new = k_y_new/dsqrt(slen)
k_z_new = k_z_new/dsqrt(slen)

! derive the zenith and azimuth angles of the reflected radiation
phi_sca = 0.5*twopi-atan2(k_y_new,-k_x_new)
sintheta_s = k_x_new/cos(phi_sca)
theta_sca = 0.5*twopi-atan2(sintheta_s,-k_z_new)

! start reflection at cloudbottom
z = cloudbottom_in
iz = 1

! contributio to absorbed flux
absorptance_ground(ix,iy) = absorptance_ground(ix,iy) + weight*co_albedo_lambert ! ix, iy from previous PhotonPath

! up date the scalar weight
weight = weight*albedo_lambert

! update the Stokes weight: as already mentioned, the update is applied only to the first Stokes element; photon loses polarization state
S_in(1) = S_in(1) * albedo_lambert
S_in(2) = 0.d0
S_in(3) = 0.d0
S_in(4) = 0.d0

! this is the mini Local Estimate Method: calculated the contribution of the reflected photon to the radiance calculation
do i_radiance = 1, n_radiance
  if (k_z_choice(i_radiance) .gt. 0.0) then
  ! vector

  S_out(1) = S_in(1) * k_z_choice(i_radiance)/pi
  S_out(2) = S_in(2)
  S_out(3) = S_in(3)
  S_out(4) = S_in(4)

  ! scalar
  radiance_weight = weight * k_z_choice(i_radiance)/pi

  call PhotonPathDetector

  end if
end do

! this might be commneted out, it was just for initial checks
if(weight.gt.1.d0) then
  write(*,*) 'ERROR : sth wrong have happened'
stop 2
endif
!
end subroutine reflection
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine STATUS
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! This subroutine stores in a file the simulation time and some general information.                                            !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
call date_and_time(date, time, zone, values)                        ! it calculates simulation time
call system_clock(end_time)                                         ! returns data from the real-time clock and the date

open(21, file = trim(mainfile)//'.status', form = 'formatted')      ! open the file
write(21,*) ' Inputfile: ', trim(mainfile)                          ! mame of the main input file
write(21,*) ' Date & time: ', date, ' ', time ! IBM                 ! date and time the simulation started
write(21,*) ' Elapsed time: ',(end_time - start_time)/1.e6, 'sec.'  ! time the simulation lasted
write(21,*) ' '
write(21,*) photon,' out of ', photon_total,' photons'              ! the number of photons used
close(21)
!
end subroutine STATUS
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
subroutine OUTPUT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!                                                                                                                               !
! Stores the so far sampled scattering data.                                                                                    !
!                                                                                                                               !
!-------------------------------------------------------------------------------------------------------------------------------!
!
call date_and_time(date, time, zone, values)
call system_clock(end_time)


! scattering phase function check (the user might need to comment out this)
help1 = 0
do i = 1, n_angle
  help1 = help1 + theta_hauf(i)
end do
open(11, file = 'output/'//trim(mainfile)//'.p11', form = 'formatted')
do i = 0, n_angle - 1
  write(11,*) angle(i+1), theta_hauf(i)/help1/arg(i+1)*2
end do
close(11)

allocate( error_to_mean(1:n_radiance) ); error_to_mean = 0d0
allocate( stokes_error_to_mean(1:4,1:n_radiance) ); stokes_error_to_mean = 0d0


! Calculate H, R, T, and A
do ix = xmin, xmax                  ! loop along x-grid boxes
  do iy = ymin, ymax                ! loop along y-grid boxes

     help1 = (transmittance_diffuse(ix,iy) + transmittance_direct(ix,iy))*(1d0 - albedo_lambert)
     net_horizontal_transport(ix,iy) = flux_in(ix,iy) - albedo(ix,iy) - help1 - absorptance(ix,iy)
     if (flux_in(ix,iy) > 0d0) then
        net_horizontal_transport(ix,iy) = net_horizontal_transport(ix,iy)/flux_in(ix,iy)
        !
        ! this way numerically more stable than H = 1 - R - T(1-R(ground)) - A
        !
        albedo(ix,iy) = albedo(ix,iy)/flux_in(ix,iy)
        transmittance_diffuse(ix,iy) = transmittance_diffuse(ix,iy)/flux_in(ix,iy)
        transmittance_direct(ix,iy) = transmittance_direct(ix,iy)/flux_in(ix,iy)
        transmittance_total(ix,iy) = transmittance_direct(ix,iy) + transmittance_diffuse(ix,iy)
        absorptance(ix,iy) = absorptance(ix,iy)/flux_in(ix,iy)
        absorptance_ground(ix,iy) = absorptance_ground(ix,iy)/flux_in(ix,iy)
     else
        net_horizontal_transport(ix,iy) = -99.99d0
        albedo(ix,iy) =  -99.99d0
        transmittance_diffuse(ix,iy) =  -99.99d0
        transmittance_direct(ix,iy) =  -99.99d0
        transmittance_total(ix,iy) =  -99.99d0
        absorptance(ix,iy) =  -99.99d0
        absorptance_ground(ix,iy) =  -99.99d0
     end if

     ! Calculate spatial scalar radiance field (+ errors)
     do i_radiance = 1, n_radiance  ! loop along i_radiance
        if (flux_in(ix,iy) > 0d0) then
           !
           !----------------------------------------------------------------------------------------------------------------!
           ! Be aware: since this code has been used for the IPRT purposes, the output is stored in a different way.        !
           ! The normalized radiance is not given by equation 3.3.8 (Barlakas, 2016b), but by equation 4.4 (Barlakas, 2016b)!
           !                                                                                                                !
           ! Bottom line, uncomment lines 1924 and 1925 and comment out 1929 and 1930                                       !
           !                                                                                                                !
           ! select_radiance(ix,iy,i_radiance) = pi*select_radiance(ix,iy,i_radiance)/flux_in(ix,iy)&                       !
           ! /dabs(k_z_choice(i_radiance))                                                                                  !
           !                                                                                                                !
           !----------------------------------------------------------------------------------------------------------------!
           !
           select_radiance(ix,iy,i_radiance) = costheta_0_in*select_radiance(ix,iy,i_radiance)/flux_in(ix,iy)&
           /dabs(k_z_choice(i_radiance))
        else
           select_radiance(ix,iy,i_radiance) =  -99.99d0
        end if
        if (error_select_radiance(ix,iy,i_radiance) .gt. 0.) then
           error_to_mean(i_radiance) = error_to_mean(i_radiance) + error_select_radiance(ix,iy,i_radiance) ! phase 2
           error_select_radiance(ix,iy,i_radiance) = 1d0/sqrt(error_select_radiance(ix,iy,i_radiance))*100d0
        else
           error_select_radiance(ix,iy,i_radiance) = -99.99d0
        end if
     end do                         ! i_radiance-loop
  end do                            ! iy-loop
end do                              ! ix-loop



! write scalar spatial radiance field
do i = 1, n_radiance
  open(10, file = 'output/'//trim(mainfile)//'.char_radiance')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(10, '(f12.5)') select_radiance(ix,iy,i)
    end do
  end do
end do
close(10)

! Calculate spatial vector radiance field (Stokes vector) (+ errors)
do ix = xmin, xmax
  do iy = ymin, ymax
    do i = 1, n_radiance

      if (flux_S(ix,iy) > 0d0) then
        !
        !----------------------------------------------------------------------------------------------------------------!
        ! Be aware: The same as for the scalar radiances,                                                                !
        ! The normalized radiance is not given by equation 3.3.8 (Barlakas, 2016b), but by equation 4.4 (Barlakas, 2016b)!
        !                                                                                                                !
        ! Bottom line, uncomment lines 1970 - 1980 and comment out 1985 - 1995                                           !
        !                                                                                                                !
        ! stokes_select_radiance(ix,iy,1,i) = pi*stokes_select_radiance(ix,iy,1,i)/flux_in(ix,iy)&                       !
        ! /dabs(k_z_choice(i))                                                                                           !
        !                                                                                                                !
        ! stokes_select_radiance(ix,iy,2,i) = pi*stokes_select_radiance(ix,iy,2,i)/flux_in(ix,iy)&                       !
        ! /dabs(k_z_choice(i))                                                                                           !
        !                                                                                                                !
        ! stokes_select_radiance(ix,iy,3,i) = pi*stokes_select_radiance(ix,iy,3,i)/flux_in(ix,iy)&                       !
        ! /dabs(k_z_choice(i))                                                                                           !
        !                                                                                                                !
        ! stokes_select_radiance(ix,iy,4,i) = pi*stokes_select_radiance(ix,iy,4,i)/flux_in(ix,iy)&                       !
        ! /dabs(k_z_choice(i))                                                                                           !
        !                                                                                                                !
        !----------------------------------------------------------------------------------------------------------------!
        !
        ! F
        stokes_select_radiance(ix,iy,1,i) = costheta_0_in*stokes_select_radiance(ix,iy,1,i)/flux_in(ix,iy)&
        /dabs(k_z_choice(i))
        ! Q
        stokes_select_radiance(ix,iy,2,i) = costheta_0_in*stokes_select_radiance(ix,iy,2,i)/flux_in(ix,iy)&
        /dabs(k_z_choice(i))
        ! U
        stokes_select_radiance(ix,iy,3,i) = costheta_0_in*stokes_select_radiance(ix,iy,3,i)/flux_in(ix,iy)&
        /dabs(k_z_choice(i))
        ! V
        stokes_select_radiance(ix,iy,4,i) = costheta_0_in*stokes_select_radiance(ix,iy,4,i)/flux_in(ix,iy)&
        /dabs(k_z_choice(i))
      else
         stokes_select_radiance(ix,iy,1,i) =  -99.99d0
         stokes_select_radiance(ix,iy,2,i) =  -99.99d0
         stokes_select_radiance(ix,iy,3,i) =  -99.99d0
         stokes_select_radiance(ix,iy,4,i) =  -99.99d0
      end if

          ! F
      if (error_stokes_select_radiance(ix,iy,1,i) .gt. 0.) then
         stokes_error_to_mean(1,i) = stokes_error_to_mean(1,i)&
         + error_stokes_select_radiance(ix,iy,1,i)
         error_stokes_select_radiance(ix,iy,1,i) = 1d0/sqrt(error_stokes_select_radiance(ix,iy,1,i))*100d0
      else
         error_stokes_select_radiance(ix,iy,1,i) =  -99.99d0
      end if

          ! Q
      if (error_stokes_select_radiance(ix,iy,2,i) .gt. 0.) then
         stokes_error_to_mean(2,i) = stokes_error_to_mean(2,i)&
         + error_stokes_select_radiance(ix,iy,2,i)
         error_stokes_select_radiance(ix,iy,2,i) = 1d0/sqrt(error_stokes_select_radiance(ix,iy,2,i))*100d0
      else
         error_stokes_select_radiance(ix,iy,2,i) =  -99.99d0
      end if

          ! U
      if (error_stokes_select_radiance(ix,iy,3,i) .gt. 0.) then
         stokes_error_to_mean(3,i) = stokes_error_to_mean(3,i)&
         + error_stokes_select_radiance(ix,iy,3,i)
         error_stokes_select_radiance(ix,iy,3,i) = 1d0/sqrt(error_stokes_select_radiance(ix,iy,3,i))*100d0
      else
         error_stokes_select_radiance(ix,iy,3,i) =  -99.99d0
      end if

          ! V
      if (error_stokes_select_radiance(ix,iy,4,i) .gt. 0.) then
         stokes_error_to_mean(4,i) = stokes_error_to_mean(4,i)&
         + error_stokes_select_radiance(ix,iy,4,i)
         error_stokes_select_radiance(ix,iy,4,i) = 1d0/sqrt(error_stokes_select_radiance(ix,iy,4,i))*100d0
      else
         error_stokes_select_radiance(ix,iy,4,i) =  -99.99d0
      end if

    end do
  end do
end do
!
! Write the spatial vector radiance field (Stokes vector) (+ errors)
!
! Write the first Stokes vector and the corresponding Monte Carlo error (the user might need to comment this out, see 1704 - 1706)
do i = 1, n_radiance
  open(41, file = 'output/'//trim(mainfile)//'.I', form = 'formatted')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(41,*) stokes_select_radiance(ix,iy,1,i), error_stokes_select_radiance(ix,iy,1,i)
    end do
  end do
end do ! n_radiance loop
close(41)

! Write the second Stokes vector and the corresponding Monte Carlo error (the user might need to comment this out, see 1704 - 1706)
do i = 1, n_radiance ! write spatial radiance field
  open(42, file = 'output/'//trim(mainfile)//'.Q', form = 'formatted')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(42,*) stokes_select_radiance(ix,iy,2,i), error_stokes_select_radiance(ix,iy,2,i)
    end do
  end do
end do ! n_radiance loop
close(42)

! Write the third Stokes vector and the corresponding Monte Carlo error (the user might need to comment this out, see 1704 - 1706)
do i = 1, n_radiance ! write spatial radiance field
  open(43, file = 'output/'//trim(mainfile)//'.U', form = 'formatted')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(43,*) stokes_select_radiance(ix,iy,3,i), error_stokes_select_radiance(ix,iy,3,i)
    end do
  end do
end do ! n_radiance loop
close(43)

! Write the forth Stokes vector and the corresponding Monte Carlo error (the user might need to comment this out, see 1704 - 1706)
do i = 1, n_radiance
  open(44, file = 'output/'//trim(mainfile)//'.V', form = 'formatted')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(44,*) stokes_select_radiance(ix,iy,4,i),error_stokes_select_radiance(ix,iy,4,i)
    end do
  end do
end do ! n_radiance loop
close(44)

! Calculate and write the degree of polarization (if not needed, the user might need to comment this out)
do i = 1, n_radiance
  open(45, file = 'output/'//trim(mainfile)//'.DegreeOfPolarization', form = 'formatted')
  do ix = xmin, xmax
    do iy = ymin, ymax
      DOP(ix,iy,i) = (sqrt(stokes_select_radiance(ix,iy,2,i)**2+stokes_select_radiance(ix,iy,3,i)**2&
      +stokes_select_radiance(ix,iy,4,i)**2))/stokes_select_radiance(ix,iy,1,i)
      write(45,'(f27.22)') DOP(ix,iy,i)
    end do
  end do
end do
close(45)
!                                                                                                                !
!----------------------------------------------------------------------------------------------------------------!
!                                                                                                                !
! The following output was stored as outlined from IPRT (ascii file); the user might want to comment this out    !
! Note that the different models might have a different geometry convention. The zenithal angle might be         !
! measured from the upward or downward normal and the azimuth angle clockwise or anticlockwise while looking     !
! upwards. To meet the geometry convention requested. Viewing zenith angles:                                     !
! if 0 - 80 is requested, run 180 - 100 (use the supplementary angles)                                           !
! if 100 - 180 is requested, run 80 - 0 (use the supplementary angles)                                           !
!                                                                                                                !
!----------------------------------------------------------------------------------------------------------------!
!
open(10, file = 'output/'//trim(mainfile)//'.B3.dat', form = 'formatted')
!
! This is the format for IPRT - phase A
!
write(10,*)'# IPRT case B3 - Standard atmosphere with aerosol profile.'
write(10,*)'# RT model: SPARTA'
write(10,*)'# depol altitude sza saa vza vaa I Q U V'

do i = 1, n_radiance
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(10,"(f4.2,2x,i2,2x,i2,2x,i3,2x,i3,2x,i3,2x,f12.8,2x,f12.8,2x,f12.8,2x,f12.8)") &
      0.03,0,30,0,180,45,stokes_select_radiance(ix,iy,1,i),stokes_select_radiance(ix,iy,2,i),&
      stokes_select_radiance(ix,iy,3,i),stokes_select_radiance(ix,iy,4,i)
    end do
  end do
end do
close(10)


! some bulk informations:
open(12, file = 'output/'//trim(mainfile)//'.status', form = 'formatted')
write(12,*) ' Inputfile: ', trim(mainfile)
write(12,*) ' RNG seed : ', iseed
write(12,*) ' Date & time: ', date, ' ', time ! IBM
write(12,*) ' Elapsed time: ',(end_time - start_time)/1.e6, 'sec.' ! time
write(12,*) ' '
write(12,*) photon - 1,' out of ', photon_total,' photons'
write(12,*) ' '
write(12,*) ' Cloud cover: ', cloud_cover_true

write(12,*) cloudleft, cloudright,cloudback, cloudfront, cloudbottom, cloudtop
close(12)

! The rest are stored as outlined from the I3RC
!
! I3RC convention for total reflectance R
open(10, file = 'output/'//trim(mainfile)//'.R') ! I3RC convention
do ix = xmin, xmax
  do iy = ymin, ymax
    write(10, '(f15.7)') albedo(ix,iy) ! I3RC convention
  end do
end do
close(10)

! I3RC convention for total transmittance T
open(10, file = 'output/'//trim(mainfile)//'.T')
do ix = xmin, xmax
  do iy = ymin, ymax
    write(10, '(f8.5)') transmittance_total(ix,iy) ! I3RC convention
  end do
end do
close(10)

! I3RC convention for absorptance A
open(10, file = 'output/'//trim(mainfile)//'.A')
do ix = xmin, xmax
  do iy = ymin, ymax
    write(10,'(f8.5)') absorptance(ix,iy) ! I3RC convention
  end do
end do
close(10)

! I3RC convention for net horizontal flux
open(10, file = 'output/'//trim(mainfile)//'.H')
  do ix = xmin, xmax
    do iy = ymin, ymax
      write(10, '(f8.5)') net_horizontal_transport(ix,iy) ! I3RC convention
  end do
end do
close(10)

! A few statistics (not stored); the user might need to comment this out
ihelp1 = (xmax - xmin + 1)*(ymax - ymin + 1)
allocate( data(1:ihelp1) )

do i = 1, n_radiance
  mean_pixel_error = 0d0
  i_data = 0
  do ix = xmin, xmax
    do iy = ymin, ymax
      help1 = error_select_radiance(ix,iy,i)
      if (help1 .gt. 0.) then
         mean_pixel_error = mean_pixel_error + help1
         i_data = i_data + 1
         data(i_data) = help1
      end if
    end do
  end do
end do ! n_radiance-loop

! Just print in the terminal the number of scattering and reflecting events (the user might want to commnet this out)
write(*,*)  'N_scat =', N_scat
write(*,*)  'N_refl =', N_refl
write(*,*)
write(*,*)  '------------------------------------------'
write(*,*)  'SPARTA simulations terminated successfully'
write(*,*)  '------------------------------------------'
write(*,*)
end subroutine OUTPUT
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!
!
end program SPARTA
!
!-------------------------------------------------------------------------------------------------------------------------------!
!###############################################################################################################################!
!###############################################################################################################################!
!-------------------------------------------------------------------------------------------------------------------------------!