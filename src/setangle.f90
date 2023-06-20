module mod_setangle
#if defined(_USE_IBM)
use mpi
#if defined(_USE_CONTACTANGLE_DYNAMIC)
use mod_contactangle, only: getTheta
#endif
use mod_param, only: dynamic_angle,theta1,theta2,zmin_ibm,d1,pi,small
use mod_common_mpi, only: myid
use mod_types
!@cuf use cudafor
!@cuf use mod_common_mpi, only: mydev
implicit none
#if defined(_SINGLE_PRECISION)
real(rp), parameter :: eps = 10._rp**(-16)
#else
real(rp), parameter :: eps = 10._rp**(-8)
#endif
!
private
public  :: setangle
!
contains
!
subroutine setangle(n,nh_v,nh_b,dl,dli,vof_h,nabs_surf,nx_surf,ny_surf,nz_surf,deltan, &
                    i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                    WP1,WP2,dir,vel_x,vel_y,vel_z,zc,zf,dzc)
!
! Routine for imposing static/dynamic contact angle using ghost cell IBM for imposing Neumann BC
! Contact line velocity calculation adapted from Patel et al. 2017 Chem. Eng. Sci. 
!
implicit none
integer,  dimension(3), intent(in)  :: n
integer,  intent(in)                :: nh_v,nh_b
real(rp), dimension(3) , intent(in) :: dl,dli
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),     intent(inout) :: vof_h
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),     intent(in)    :: nabs_surf,nx_surf,ny_surf,nz_surf,deltan
integer,  dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),     intent(in)    :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v,1:7), intent(in)    :: WP1,WP2
real(rp), allocatable, dimension(:,:,:), save                                              :: vof
real(rp), dimension(0:,0:,0:), intent(in) :: vel_x, vel_y, vel_z
real(rp), dimension(0:), intent(in) :: zc,zf,dzc
real(rp), dimension(3) :: tang, vel_cl_xyz
integer,  intent(in) :: dir
integer  :: i,j,k,ii1,jj1,kk1,ii2,jj2,kk2,n1,n2,n3,nv,nb,do_cox
real(rp) :: temp, chem_m, dphidx, dphidy, dphidz, dPhidN, dPhidTAbs2, d, d_mirror, thetad, vel_cl, intdir, z, z_bottom
real(rp) :: dl_l, theta_1, theta_2, pi_value
real(rp) :: dx,dy,dz,dxi,dyi,dzi
real(rp) :: q_1,BB1,WW1_1,WW2_1,WW3_1,WW4_1,WW5_1,WW6_1,WW7_1,B1_1,B2_1,B3_1,B4_1,B5_1,B6_1,B7_1
real(rp) :: q_2,BB2,WW1_2,WW2_2,WW3_2,WW4_2,WW5_2,WW6_2,WW7_2,B1_2,B2_2,B3_2,B4_2,B5_2,B6_2,B7_2
real(rp) :: phi_ip_11,phi_im_11,phi_jp_11,phi_jm_11,phi_kp_11,phi_km_11
real(rp) :: phi_ip_21,phi_im_21,phi_jp_21,phi_jm_21,phi_kp_21,phi_km_21
real(rp) :: phi_ip_31,phi_im_31,phi_jp_31,phi_jm_31,phi_kp_31,phi_km_31
real(rp) :: phi_ip_41,phi_im_41,phi_jp_41,phi_jm_41,phi_kp_41,phi_km_41
real(rp) :: phi_ip_51,phi_im_51,phi_jp_51,phi_jm_51,phi_kp_51,phi_km_51
real(rp) :: phi_ip_61,phi_im_61,phi_jp_61,phi_jm_61,phi_kp_61,phi_km_61
real(rp) :: phi_ip_71,phi_im_71,phi_jp_71,phi_jm_71,phi_kp_71,phi_km_71
real(rp) :: phi_ip_12,phi_im_12,phi_jp_12,phi_jm_12,phi_kp_12,phi_km_12
real(rp) :: phi_ip_22,phi_im_22,phi_jp_22,phi_jm_22,phi_kp_22,phi_km_22
real(rp) :: phi_ip_32,phi_im_32,phi_jp_32,phi_jm_32,phi_kp_32,phi_km_32
real(rp) :: phi_ip_42,phi_im_42,phi_jp_42,phi_jm_42,phi_kp_42,phi_km_42
real(rp) :: phi_ip_52,phi_im_52,phi_jp_52,phi_jm_52,phi_kp_52,phi_km_52
real(rp) :: phi_ip_62,phi_im_62,phi_jp_62,phi_jm_62,phi_kp_62,phi_km_62
real(rp) :: phi_ip_72,phi_im_72,phi_jp_72,phi_jm_72,phi_kp_72,phi_km_72
real(rp) :: phi_x_11,phi_y_11,phi_z_11
real(rp) :: phi_x_21,phi_y_21,phi_z_21
real(rp) :: phi_x_31,phi_y_31,phi_z_31
real(rp) :: phi_x_41,phi_y_41,phi_z_41
real(rp) :: phi_x_51,phi_y_51,phi_z_51
real(rp) :: phi_x_61,phi_y_61,phi_z_61
real(rp) :: phi_x_71,phi_y_71,phi_z_71
real(rp) :: phi_x_12,phi_y_12,phi_z_12
real(rp) :: phi_x_22,phi_y_22,phi_z_22
real(rp) :: phi_x_32,phi_y_32,phi_z_32
real(rp) :: phi_x_42,phi_y_42,phi_z_42
real(rp) :: phi_x_52,phi_y_52,phi_z_52
real(rp) :: phi_x_62,phi_y_62,phi_z_62
real(rp) :: phi_x_72,phi_y_72,phi_z_72
integer  :: dphidx1,dphidy1,dphidz1,dphidx2,dphidy2,dphidz2
!@cuf integer :: istat
!@cuf attributes(managed) :: vof,vof_h,vel_x,vel_y,vel_z,zc,zf,dzc
!@cuf attributes(managed) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2
!@cuf attributes(managed) :: nabs_surf,nx_surf,ny_surf,nz_surf,deltan
!@cuf attributes(managed) :: tang,vel_cl_xyz
!
if(.not.allocated(vof)) then
 allocate(vof(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v))
 !@cuf istat = cudaMemPrefetchAsync(vof, size(vof), mydev, 0)
endif
!
!$acc kernels
vof(:,:,:) = vof_h(:,:,:)
!$acc end kernels
!
z_bottom = zmin_ibm
n1=n(1)
n2=n(2)
n3=n(3)
nv = nh_v
nb = nh_b
dl_l = dl(dir)
dx = dl(1)
dy = dl(2)
dz = dl(3)
dxi = dli(1)
dyi = dli(2)
dzi = dli(3)
theta_1 = theta1
theta_2 = theta2
d_mirror = d1
pi_value = pi
!
!$acc kernels
 do k=1,n3
  do j=1,n2
   do i=1,n1
    if (nabs_surf(i,j,k).gt.small) then
       dzi = dzc(k)**(-1)
       d = d_mirror*sqrt(dx**2 + dy**2 + dzc(k)**2)

       ii1   = i_IP1(i,j,k) ! Coordinates for cell of interpolation point 1 (IP1)
       jj1   = j_IP1(i,j,k)
       kk1   = k_IP1(i,j,k)
       ii2   = i_IP2(i,j,k) ! Coordinates for cell of interpolation point 2  (IP2)
       jj2   = j_IP2(i,j,k)
       kk2   = k_IP2(i,j,k)

       WW1_1 = WP1(i,j,k,1)  ! Weights for interpolating to IP1
       WW2_1 = WP1(i,j,k,2)
       WW3_1 = WP1(i,j,k,3)
       WW4_1 = WP1(i,j,k,4)
       WW5_1 = WP1(i,j,k,5)
       WW6_1 = WP1(i,j,k,6)
       WW7_1 = WP1(i,j,k,7)
       WW1_2 = WP2(i,j,k,1)  ! Weights for interpolating to IP2
       WW2_2 = WP2(i,j,k,2)
       WW3_2 = WP2(i,j,k,3)
       WW4_2 = WP2(i,j,k,4)
       WW5_2 = WP2(i,j,k,5)
       WW6_2 = WP2(i,j,k,6)
       WW7_2 = WP2(i,j,k,7)

       B1_1  = vof(ii1,jj1+1,kk1)  ! Values for interpolating to IP2
       B2_1  = vof(ii1,jj1,kk1+1)
       B3_1  = vof(ii1,jj1-1,kk1)
       B4_1  = vof(ii1,jj1,kk1-1)
       B5_1  = vof(ii1,jj1,kk1)
       B6_1  = vof(ii1-1,jj1,kk1)
       B7_1  = vof(ii1+1,jj1,kk1)
       B1_2  = vof(ii2,jj2+1,kk2)  ! Values for interpolating to IP2
       B2_2  = vof(ii2,jj2,kk2+1)
       B3_2  = vof(ii2,jj2-1,kk2)
       B4_2  = vof(ii2,jj2,kk2-1)
       B5_2  = vof(ii2,jj2,kk2)
       B6_2  = vof(ii2-1,jj2,kk2)
       B7_2  = vof(ii2+1,jj2,kk2)

       ! Sum weights
       q_1 = WW1_1 + WW2_1 + WW3_1 + WW4_1 + WW5_1 + WW6_1 + WW7_1
       q_2 = WW1_2 + WW2_2 + WW3_2 + WW4_2 + WW5_2 + WW6_2 + WW7_2

       ! Values at interpolation points
       BB1 = (1._rp/(q_1+eps))*(WW1_1*B1_1 + WW2_1*B2_1 + WW3_1*B3_1 + WW4_1*B4_1 + WW5_1*B5_1 + WW6_1*B6_1 + WW7_1*B7_1)
       BB2 = (1._rp/(q_2+eps))*(WW1_2*B1_2 + WW2_2*B2_2 + WW3_2*B3_2 + WW4_2*B4_2 + WW5_2*B5_2 + WW6_2*B6_2 + WW7_2*B7_2)

       chem_m = 2.0_rp*BB1-BB2 ! Value at mirror point

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Cell edge values for AP for IP1
       phi_ip_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1+1,jj1+1,kk1))
       phi_im_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1-1,jj1+1,kk1))
       phi_jp_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1,jj1+2,kk1))
       phi_jm_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1,jj1,kk1))
       phi_kp_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1,jj1+1,kk1+1))
       phi_km_11 = 0.5_rp*(vof(ii1,jj1+1,kk1)+vof(ii1,jj1+1,kk1-1))

       phi_ip_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1+1,jj1,kk1+1))
       phi_im_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1-1,jj1-1,kk1+1))
       phi_jp_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1,jj1+1,kk1+1))
       phi_jm_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1,jj1-1,kk1+1))
       phi_kp_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1,jj1,kk1+2))
       phi_km_21 = 0.5_rp*(vof(ii1,jj1,kk1+1)+vof(ii1,jj1,kk1))

       phi_ip_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1+1,jj1-1,kk1))
       phi_im_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1-1,jj1-1,kk1))
       phi_jp_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1,jj1,kk1))
       phi_jm_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1,jj1-2,kk1))
       phi_kp_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1,jj1-1,kk1+1))
       phi_km_31 = 0.5_rp*(vof(ii1,jj1-1,kk1)+vof(ii1,jj1-1,kk1-1))

       phi_ip_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1+1,jj1,kk1-1))
       phi_im_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1-1,jj1-1,kk1-1))
       phi_jp_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1,jj1+1,kk1-1))
       phi_jm_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1,jj1-1,kk1-1))
       phi_kp_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1,jj1,kk1))
       phi_km_41 = 0.5_rp*(vof(ii1,jj1,kk1-1)+vof(ii1,jj1,kk1-2))

       phi_ip_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1+1,jj1,kk1))
       phi_im_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1-1,jj1,kk1))
       phi_jp_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1,jj1+1,kk1))
       phi_jm_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1,jj1-1,kk1))
       phi_kp_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1,jj1,kk1+1))
       phi_km_51 = 0.5_rp*(vof(ii1,jj1,kk1)+vof(ii1,jj1,kk1-1))

       phi_ip_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1,jj1,kk1))
       phi_im_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1-2,jj1,kk1))
       phi_jp_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1-1,jj1+1,kk1))
       phi_jm_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1-1,jj1-1,kk1))
       phi_kp_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1-1,jj1,kk1+1))
       phi_km_61 = 0.5_rp*(vof(ii1-1,jj1,kk1)+vof(ii1-1,jj1,kk1-1))

       phi_ip_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1+2,jj1,kk1))
       phi_im_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1,jj1,kk1))
       phi_jp_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1+1,jj1+1,kk1))
       phi_jm_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1+1,jj1-1,kk1))
       phi_kp_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1+1,jj1,kk1+1))
       phi_km_71 = 0.5_rp*(vof(ii1+1,jj1,kk1)+vof(ii1+1,jj1,kk1-1))

       ! Cell edge values for AP for IP1
       phi_ip_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2+1,jj2+1,kk1))
       phi_im_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2-1,jj2+1,kk1))
       phi_jp_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2,jj2+2,kk1))
       phi_jm_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2,jj2,kk1))
       phi_kp_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2,jj2+1,kk1+1))
       phi_km_12 = 0.5_rp*(vof(ii2,jj2+1,kk2)+vof(ii2,jj2+1,kk1-1))

       phi_ip_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2+1,jj2,kk1+1))
       phi_im_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2-1,jj2-1,kk1+1))
       phi_jp_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2,jj2+1,kk1+1))
       phi_jm_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2,jj2-1,kk1+1))
       phi_kp_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2,jj2,kk1+2))
       phi_km_22 = 0.5_rp*(vof(ii2,jj2,kk2+1)+vof(ii2,jj2,kk1))

       phi_ip_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2+1,jj2-1,kk2))
       phi_im_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2-1,jj2-1,kk2))
       phi_jp_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2,jj2,kk2))
       phi_jm_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2,jj2-2,kk2))
       phi_kp_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2,jj2-1,kk2+1))
       phi_km_32 = 0.5_rp*(vof(ii2,jj2-1,kk2)+vof(ii2,jj2-1,kk2-1))

       phi_ip_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2+1,jj2,kk2-1))
       phi_im_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2-1,jj2-1,kk2-1))
       phi_jp_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2,jj2+1,kk2-1))
       phi_jm_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2,jj2-1,kk2-1))
       phi_kp_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2,jj2,kk2))
       phi_km_42 = 0.5_rp*(vof(ii2,jj2,kk2-1)+vof(ii2,jj2,kk2-2))

       phi_ip_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2+1,jj2,kk2))
       phi_im_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2-1,jj2,kk2))
       phi_jp_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2,jj2+1,kk2))
       phi_jm_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2,jj2-1,kk2))
       phi_kp_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2,jj2,kk2+1))
       phi_km_52 = 0.5_rp*(vof(ii2,jj2,kk2)+vof(ii2,jj2,kk2-1))

       phi_ip_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2,jj2,kk2))
       phi_im_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2-2,jj2,kk2))
       phi_jp_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2-1,jj2+1,kk2))
       phi_jm_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2-1,jj2-1,kk2))
       phi_kp_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2-1,jj2,kk2+1))
       phi_km_62 = 0.5_rp*(vof(ii2-1,jj2,kk2)+vof(ii2-1,jj2,kk2-1))

       phi_ip_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2+2,jj2,kk2))
       phi_im_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2,jj2,kk2))
       phi_jp_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2+1,jj2+1,kk2))
       phi_jm_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2+1,jj2-1,kk2))
       phi_kp_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2+1,jj2,kk2+1))
       phi_km_72 = 0.5_rp*(vof(ii2+1,jj2,kk2)+vof(ii2+1,jj2,kk2-1))

       ! Derivatives at AP
       phi_x_11 = (phi_ip_11-phi_im_11)*dxi
       phi_y_11 = (phi_jp_11-phi_jm_11)*dyi
       phi_z_11 = (phi_kp_11-phi_km_11)*dzi
       phi_x_21 = (phi_ip_21-phi_im_21)*dxi
       phi_y_21 = (phi_jp_21-phi_jm_21)*dyi
       phi_z_21 = (phi_kp_21-phi_km_21)*dzi
       phi_x_31 = (phi_ip_31-phi_im_31)*dxi
       phi_y_31 = (phi_jp_31-phi_jm_31)*dyi
       phi_z_31 = (phi_kp_31-phi_km_31)*dzi
       phi_x_41 = (phi_ip_41-phi_im_41)*dxi
       phi_y_41 = (phi_jp_41-phi_jm_41)*dyi
       phi_z_41 = (phi_kp_41-phi_km_41)*dzi
       phi_x_51 = (phi_ip_51-phi_im_51)*dxi
       phi_y_51 = (phi_jp_51-phi_jm_51)*dyi
       phi_z_51 = (phi_kp_51-phi_km_51)*dzi
       phi_x_61 = (phi_ip_61-phi_im_61)*dxi
       phi_y_61 = (phi_jp_61-phi_jm_61)*dyi
       phi_z_61 = (phi_kp_61-phi_km_61)*dzi
       phi_x_71 = (phi_ip_71-phi_im_71)*dxi
       phi_y_71 = (phi_jp_71-phi_jm_71)*dyi
       phi_z_71 = (phi_kp_71-phi_km_71)*dzi
       phi_x_12 = (phi_ip_12-phi_im_12)*dxi
       phi_y_12 = (phi_jp_12-phi_jm_12)*dyi
       phi_z_12 = (phi_kp_12-phi_km_12)*dzi
       phi_x_22 = (phi_ip_22-phi_im_22)*dxi
       phi_y_22 = (phi_jp_22-phi_jm_22)*dyi
       phi_z_22 = (phi_kp_22-phi_km_22)*dzi
       phi_x_32 = (phi_ip_32-phi_im_32)*dxi
       phi_y_32 = (phi_jp_32-phi_jm_32)*dyi
       phi_z_32 = (phi_kp_32-phi_km_32)*dzi
       phi_x_42 = (phi_ip_42-phi_im_42)*dxi
       phi_y_42 = (phi_jp_42-phi_jm_42)*dyi
       phi_z_42 = (phi_kp_42-phi_km_42)*dzi
       phi_x_52 = (phi_ip_52-phi_im_52)*dxi
       phi_y_52 = (phi_jp_52-phi_jm_52)*dyi
       phi_z_52 = (phi_kp_52-phi_km_52)*dzi
       phi_x_62 = (phi_ip_62-phi_im_62)*dxi
       phi_y_62 = (phi_jp_62-phi_jm_62)*dyi
       phi_z_62 = (phi_kp_62-phi_km_62)*dzi
       phi_x_72 = (phi_ip_72-phi_im_72)*dxi
       phi_y_72 = (phi_jp_72-phi_jm_72)*dyi
       phi_z_72 = (phi_kp_72-phi_km_72)*dzi

       ! Derivatives at mirror point
       dphidx1 = (1._rp/(q_1+eps))*(WW1_1*phi_x_11 + WW2_1*phi_x_21 + WW3_1*phi_x_31 + WW4_1*phi_x_41 + WW5_1*phi_x_51 + WW6_1*phi_x_61 + WW7_1*phi_x_71)
       dphidy1 = (1._rp/(q_1+eps))*(WW1_1*phi_y_11 + WW2_1*phi_y_21 + WW3_1*phi_y_31 + WW4_1*phi_y_41 + WW5_1*phi_y_51 + WW6_1*phi_y_61 + WW7_1*phi_y_71)
       dphidz1 = (1._rp/(q_1+eps))*(WW1_1*phi_z_11 + WW2_1*phi_z_21 + WW3_1*phi_z_31 + WW4_1*phi_z_41 + WW5_1*phi_z_51 + WW6_1*phi_z_61 + WW7_1*phi_z_71)
       dphidx2 = (1._rp/(q_2+eps))*(WW1_2*phi_x_12 + WW2_2*phi_x_22 + WW3_2*phi_x_32 + WW4_2*phi_x_42 + WW5_2*phi_x_52 + WW6_2*phi_x_62 + WW7_2*phi_x_72)
       dphidy2 = (1._rp/(q_2+eps))*(WW1_2*phi_y_12 + WW2_2*phi_y_22 + WW3_2*phi_y_32 + WW4_2*phi_y_42 + WW5_2*phi_y_52 + WW6_2*phi_y_62 + WW7_2*phi_y_72)
       dphidz2 = (1._rp/(q_2+eps))*(WW1_2*phi_z_12 + WW2_2*phi_z_22 + WW3_2*phi_z_32 + WW4_2*phi_z_42 + WW5_2*phi_z_52 + WW6_2*phi_z_62 + WW7_2*phi_z_72)

       ! Extrapolate to mirror point
       dphidx = 2.0_rp*dphidx1-dphidx2
       dphidy = 2.0_rp*dphidy1-dphidy2
       dphidz = 2.0_rp*dphidz1-dphidz2

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       dPhidN = dphidx*nx_surf(i,j,k) + dphidy*ny_surf(i,j,k) + dphidz*nz_surf(i,j,k)  ! Project gradient onto surface normal
       !
       dPhidTAbs2 = (dphidx - dPhidN*nx_surf(i,j,k))**2 + &
                    (dphidy - dPhidN*ny_surf(i,j,k))**2 + &
                    (dphidz - dPhidN*nz_surf(i,j,k))**2   ! Find absolute value squared of surface tangent
#if defined(_USE_CONTACTANGLE_DYNAMIC)
        vel_cl = 0._rp
        vel_cl_xyz(:) = 0._rp

        tang(1) = (dphidx - dPhidN*nx_surf(i,j,k))/sqrt(dPhidTAbs2) ! Components of unit surface tangent vector
        tang(2) = (dphidy - dPhidN*ny_surf(i,j,k))/sqrt(dPhidTAbs2)
        tang(3) = (dphidz - dPhidN*nz_surf(i,j,k))/sqrt(dPhidTAbs2)

        temp = (0.5_rp*(vel_x(i,j,k) + vel_x(i-1,j,k)))*tang(1) + &
               (0.5_rp*(vel_y(i,j,k) + vel_y(i,j-1,k)))*tang(2) + &
               (0.5_rp*(vel_z(i,j,k) + vel_z(i,j,k-1)))*tang(3)     ! Projected velocity along surface

        vel_cl_xyz(1) = temp*tang(1)/sqrt(dPhidTAbs2)
        vel_cl_xyz(2) = temp*tang(2)/sqrt(dPhidTAbs2)
        vel_cl_xyz(3) = temp*tang(3)/sqrt(dPhidTAbs2)

        vel_cl = sqrt(vel_cl_xyz(1)**2 + vel_cl_xyz(2)**2 + vel_cl_xyz(3)**2)

        intdir = sign(1._rp,( vel_cl_xyz(1)*dphidx/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) + &
                              vel_cl_xyz(2)*dphidy/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) + &
                              vel_cl_xyz(3)*dphidz/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) ) ) !Interface movement direction

        vel_cl = intdir*vel_cl ! Contact line velocity

        if( (z-z_bottom).le.dzc(k)) call getTheta(theta_1*pi_value/180._rp, dl_l, thetad, vel_cl)  ! Get dynamic contact angle (note: thetad is in radians)
        if( (z-z_bottom).gt.dzc(k)) call getTheta(theta_2*pi_value/180._rp, dl_l, thetad, vel_cl)
#else
       ! if((zc(k)-z_bottom).le.dzc(k))then
        ! thetad = theta_1*pi_value/180.0_rp
       ! elseif((zc(k)-z_bottom).gt.dzc(k)) then
        ! thetad = theta_2*pi_value/180.0_rp
       ! endif
       thetad = theta_1!*pi_value/180.0_rp
#endif
       dPhidN = cosd(thetad)*sqrt((1.0_rp/(1.0_rp - (cosd(thetad))**2))*dPhidTAbs2)    ! Adjust VoF gradient
       !
       vof_h(i,j,k) = chem_m + dPhidN*(d + deltan(i,j,k))
    endif
   enddo
  enddo
 enddo
!$acc end kernels
!
!@cuf istat=cudaDeviceSynchronize()
!
end subroutine setangle
#endif
end module mod_setangle