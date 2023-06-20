module mod_initIBM
#if defined(_USE_IBM)
use mpi
use mod_IBM,        only: IBM_Mask,normal_vectors,intersect,mirrorpoints, &
                          mirrorpoints_ijk,interpolation_mirror,InterpolationWeights
use mod_bound,      only: set_bc, updthalo
use mod_load,       only: load,load_int,load_weight
use mod_param,      only: l,ng,cbcpre,bcpre,cbcvel,bcvel,is_outflow
use mod_output    , only: write_visu_3d
use mod_common_mpi, only: myid,ierr,left,right,front,back,top,bottom
use mod_types
!@cuf use cudafor
!@cuf use mod_common_mpi, only: mydev
implicit none
private
public initIBM
contains
subroutine initIBM(dims,n,ng,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag,Level_set, &
                   nx_surf,ny_surf,nz_surf,nabs_surf, &
                   i_mirror,j_mirror,k_mirror, &
                   i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                   WP1,WP2,deltan, &
                   nh_d,nh_b,halo, &
                   zc,zf,dzc,dzf,dl,dli)
implicit none
character(len=100) :: data_dir,restart_dir
integer,  intent(in)                         :: nh_d,nh_b
integer,  dimension(3), intent(in)           :: dims,n,ng,halo
real(rp), dimension(1-nh_d:), intent(in)     :: zc,zf,dzc,dzf
real(rp), dimension(3), intent(in)           :: dl,dli
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(out) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(out) :: Level_set
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(out) :: i_mirror,j_mirror,k_mirror
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(out) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(out) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7), intent(out) :: WP1,WP2
real(rp), dimension(:,:,:), allocatable, save                                            :: x_intersect,y_intersect,z_intersect
real(rp), dimension(:,:,:), allocatable, save                                            :: x_IP1,y_IP1,z_IP1
real(rp), dimension(:,:,:), allocatable, save                                            :: x_IP2,y_IP2,z_IP2
real(rp), dimension(:,:,:), allocatable, save                                            :: x_mirror,y_mirror,z_mirror
real(rp), dimension(:,:,:), allocatable, save                                            :: temp
logical :: is_data
integer :: i,j,k,q,n1,n2,n3
integer :: ind1,ind2
integer :: istep
character(len=9) :: fldnum
!
!@cuf integer :: istat
!@cuf attributes(managed) :: dzc,dzf,zc,zf
!@cuf attributes(managed) :: Level_set,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
!@cuf attributes(managed) :: i_mirror,j_mirror,k_mirror
!@cuf attributes(managed) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
!@cuf attributes(managed) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
!@cuf attributes(managed) :: WP1,WP2
!@cuf attributes(managed) :: x_intersect,y_intersect,z_intersect
!@cuf attributes(managed) :: x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2
!@cuf attributes(managed) :: x_mirror,y_mirror,z_mirror
!@cuf attributes(managed) :: temp
!
#if defined(_DECOMP_X)
ind1 = 2
ind2 = 3
#elif _DECOMP_Y
ind1 = 1
ind2 = 3
#else /*_DECOMP_Z*/
ind1 = 1
ind2 = 2
#endif
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
data_dir = 'data/'
restart_dir = 'data/restart_dir/'
!
inquire(file=trim(restart_dir)//'stag.bin',exist=is_data)
!
if(.not.allocated(temp)) allocate(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
!
if (.not.is_data) then
  if(.not.allocated(x_mirror))    allocate(x_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(y_mirror))    allocate(y_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(z_mirror))    allocate(z_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(x_IP1))       allocate(x_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(y_IP1))       allocate(y_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(z_IP1))       allocate(z_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(x_IP2))       allocate(x_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(y_IP2))       allocate(y_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(z_IP2))       allocate(z_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(x_intersect)) allocate(x_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(y_intersect)) allocate(y_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  if(.not.allocated(z_intersect)) allocate(z_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))

  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set),       mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf),       mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan),             mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_mirror, size(i_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_mirror, size(j_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_mirror, size(k_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1, size(WP1),                   mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2, size(WP2),                   mydev, 0)
  !@cuf istat = cudaMemAdvise(x_intersect, size(x_intersect),   cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(y_intersect, size(y_intersect),   cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(z_intersect, size(z_intersect),   cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(x_IP1, size(x_IP1),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(y_IP1, size(y_IP1),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(z_IP1, size(z_IP1),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(x_IP2, size(x_IP2),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(y_IP2, size(y_IP2),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(z_IP2, size(z_IP2),               cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(x_mirror, size(x_mirror),         cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(y_mirror, size(y_mirror),         cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(z_mirror, size(z_mirror),         cudaMemAdviseSetPreferredLocation, mydev)
  !@cuf istat = cudaMemAdvise(temp, size(temp),                 cudaMemAdviseSetPreferredLocation, mydev)
  !
  !$acc kernels
  cell_phi_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 1.0_rp
  Level_set(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 1
  !$acc end kernels
  !
  call IBM_mask(n,ng,nh_d,nh_b,cell_phi_tag,Level_set,zc,zf,dl,dzc)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_phi_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_phi_tag)
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(Level_set(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  Level_set(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    cell_phi_tag( 0,:,:)  = cell_phi_tag(1,:,:)
    cell_phi_tag(-1,:,:)  = cell_phi_tag(1,:,:)
    cell_phi_tag(-2,:,:)  = cell_phi_tag(1,:,:)
    cell_phi_tag(-3,:,:)  = cell_phi_tag(1,:,:)
    cell_phi_tag(-4,:,:)  = cell_phi_tag(1,:,:)
    cell_phi_tag(-5,:,:)  = cell_phi_tag(1,:,:)
    !
    Level_set( 0,:,:)     = Level_set(1,:,:)
    Level_set(-1,:,:)     = Level_set(1,:,:)
    Level_set(-2,:,:)     = Level_set(1,:,:)
    Level_set(-3,:,:)     = Level_set(1,:,:)
    Level_set(-4,:,:)     = Level_set(1,:,:)
    Level_set(-5,:,:)     = Level_set(1,:,:)
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    cell_phi_tag(n1+1,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+2,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+3,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+4,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+5,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+6,:,:) = cell_phi_tag(n1,:,:)
    !
    Level_set(n1+1,:,:)    = Level_set(n1,:,:)
    Level_set(n1+2,:,:)    = Level_set(n1,:,:)
    Level_set(n1+3,:,:)    = Level_set(n1,:,:)
    Level_set(n1+4,:,:)    = Level_set(n1,:,:)
    Level_set(n1+5,:,:)    = Level_set(n1,:,:)
    Level_set(n1+6,:,:)    = Level_set(n1,:,:)
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    cell_phi_tag(:, 0,:)  = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-1,:)  = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-2,:)  = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-3,:)  = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-4,:)  = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-5,:)  = cell_phi_tag(:,1,:)
    !
    Level_set(:, 0,:)     = Level_set(:,1,:)
    Level_set(:,-1,:)     = Level_set(:,1,:)
    Level_set(:,-2,:)     = Level_set(:,1,:)
    Level_set(:,-3,:)     = Level_set(:,1,:)
    Level_set(:,-4,:)     = Level_set(:,1,:)
    Level_set(:,-5,:)     = Level_set(:,1,:)
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    cell_phi_tag(:,n1+1,:) = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+2,:) = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+3,:) = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+4,:) = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+5,:) = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+6,:) = cell_phi_tag(:,n1,:)
    !
    Level_set(:,n1+1,:)    = Level_set(:,n1,:)
    Level_set(:,n1+2,:)    = Level_set(:,n1,:)
    Level_set(:,n1+3,:)    = Level_set(:,n1,:)
    Level_set(:,n1+4,:)    = Level_set(:,n1,:)
    Level_set(:,n1+5,:)    = Level_set(:,n1,:)
    Level_set(:,n1+6,:)    = Level_set(:,n1,:)
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    cell_phi_tag(:,:, 0)  = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-1)  = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-2)  = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-3)  = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-4)  = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-5)  = cell_phi_tag(:,:,1)
    !
    Level_set(:,:, 0)     = Level_set(:,:,1)
    Level_set(:,:,-1)     = Level_set(:,:,1)
    Level_set(:,:,-2)     = Level_set(:,:,1)
    Level_set(:,:,-3)     = Level_set(:,:,1)
    Level_set(:,:,-4)     = Level_set(:,:,1)
    Level_set(:,:,-5)     = Level_set(:,:,1)
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    cell_phi_tag(:,:,n3+1) = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+2) = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+3) = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+4) = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+5) = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+6) = cell_phi_tag(:,:,n3)
    !
    Level_set(:,:,n3+1)   = Level_set(:,:,n3)
    Level_set(:,:,n3+2)   = Level_set(:,:,n3)
    Level_set(:,:,n3+3)   = Level_set(:,:,n3)
    Level_set(:,:,n3+4)   = Level_set(:,:,n3)
    Level_set(:,:,n3+5)   = Level_set(:,:,n3)
    Level_set(:,:,n3+6)   = Level_set(:,:,n3)
  endif
  !

  if (myid.eq.0)  print*, '**** Solid marker set ****'
  !---------------------------------------------------------------------
   !
   cell_u_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 1.0_rp
   cell_v_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 1.0_rp
   cell_w_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 1.0_rp
   !
   do k=0,n3
    do j=0,n2
     do i=0,n1
      cell_u_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j,k)+cell_phi_tag(i+1,j,k))
      cell_v_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j,k)+cell_phi_tag(i,j+1,k))
      cell_w_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j,k)+cell_phi_tag(i,j,k+1))
     enddo
    enddo
   enddo
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_u_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_u_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_v_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_v_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_w_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_w_tag)
  !
  ! along x
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    cell_u_tag( 0,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-1,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-2,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-3,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-4,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-5,:,:)     = cell_u_tag(1,:,:)
    !
    cell_v_tag( 0,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-1,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-2,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-3,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-4,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-5,:,:)     = cell_v_tag(1,:,:)
    !
    cell_w_tag( 0,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-1,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-2,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-3,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-4,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-5,:,:)     = cell_w_tag(1,:,:)
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    cell_u_tag(n1+1,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+2,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+3,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+4,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+5,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+6,:,:)   = cell_u_tag(n1,:,:)
    !
    cell_v_tag(n1+1,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+2,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+3,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+4,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+5,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+6,:,:)   = cell_v_tag(n1,:,:)
    !
    cell_w_tag(n1+1,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+2,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+3,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+4,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+5,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+6,:,:)   = cell_w_tag(n1,:,:)
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    cell_u_tag(:, 0,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-1,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-2,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-3,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-4,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-5,:)     = cell_u_tag(:,1,:)
    !
    cell_v_tag(:, 0,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-1,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-2,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-3,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-4,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-5,:)     = cell_v_tag(:,1,:)
    !
    cell_w_tag(:, 0,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-1,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-2,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-3,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-4,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-5,:)     = cell_w_tag(:,1,:)
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    cell_u_tag(:,n1+1,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+2,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+3,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+4,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+5,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+6,:)   = cell_u_tag(:,n1,:)
    !
    cell_v_tag(:,n1+1,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+2,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+3,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+4,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+5,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+6,:)   = cell_v_tag(:,n1,:)
    !
    cell_w_tag(:,n1+1,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+2,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+3,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+4,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+5,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+6,:)   = cell_w_tag(:,n1,:)
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    cell_u_tag(:,:, 0)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-1)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-2)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-3)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-4)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-5)     = cell_u_tag(:,:,1)
    !
    cell_v_tag(:,:, 0)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-1)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-2)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-3)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-4)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-5)     = cell_v_tag(:,:,1)
    !
    cell_w_tag(:,:, 0)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-1)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-2)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-3)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-4)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-5)     = cell_w_tag(:,:,1)
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    cell_u_tag(:,:,n3+1)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+2)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+3)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+4)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+5)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+6)   = cell_u_tag(:,:,n3)
    !
    cell_v_tag(:,:,n3+1)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+2)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+3)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+4)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+5)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+6)   = cell_v_tag(:,:,n3)
    !
    cell_w_tag(:,:,n3+1)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+2)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+3)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+4)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+5)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+6)   = cell_w_tag(:,:,n3)
  endif
  !

  if (myid.eq.0)  print*, '*** Volume fractions calculated ***'
  !---------------------------------------------------------------------
  nx_surf(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  ny_surf(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  nz_surf(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  nabs_surf(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)  = 0.0_rp
  !
  call normal_vectors(n,ng,nh_d,nh_b,Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,dl,dli,zc,dzc)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nabs_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nabs_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nx_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nx_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,ny_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,ny_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nz_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nz_surf)
  !
  ! along x
  !
  if(left.eq.MPI_PROC_NULL) then ! x - bottom part
    nabs_surf( 0,:,:)   = nabs_surf(1,:,:)
    nabs_surf(-1,:,:)   = nabs_surf(1,:,:)
    nabs_surf(-2,:,:)   = nabs_surf(1,:,:)
    nabs_surf(-3,:,:)   = nabs_surf(1,:,:)
    nabs_surf(-4,:,:)   = nabs_surf(1,:,:)
    nabs_surf(-5,:,:)   = nabs_surf(1,:,:)
    !
    nx_surf( 0,:,:)     = nx_surf(1,:,:)
    nx_surf(-1,:,:)     = nx_surf(1,:,:)
    nx_surf(-2,:,:)     = nx_surf(1,:,:)
    nx_surf(-3,:,:)     = nx_surf(1,:,:)
    nx_surf(-4,:,:)     = nx_surf(1,:,:)
    nx_surf(-5,:,:)     = nx_surf(1,:,:)
    !
    ny_surf( 0,:,:)     = ny_surf(1,:,:)
    ny_surf(-1,:,:)     = ny_surf(1,:,:)
    ny_surf(-2,:,:)     = ny_surf(1,:,:)
    ny_surf(-3,:,:)     = ny_surf(1,:,:)
    ny_surf(-4,:,:)     = ny_surf(1,:,:)
    ny_surf(-5,:,:)     = ny_surf(1,:,:)
    !
    nz_surf( 0,:,:)     = nz_surf(1,:,:)
    nz_surf(-1,:,:)     = nz_surf(1,:,:)
    nz_surf(-2,:,:)     = nz_surf(1,:,:)
    nz_surf(-3,:,:)     = nz_surf(1,:,:)
    nz_surf(-4,:,:)     = nz_surf(1,:,:)
    nz_surf(-5,:,:)     = nz_surf(1,:,:)
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    nabs_surf(n1+1,:,:) = nabs_surf(n1,:,:)
    nabs_surf(n1+2,:,:) = nabs_surf(n1,:,:)
    nabs_surf(n1+3,:,:) = nabs_surf(n1,:,:)
    nabs_surf(n1+4,:,:) = nabs_surf(n1,:,:)
    nabs_surf(n1+5,:,:) = nabs_surf(n1,:,:)
    nabs_surf(n1+6,:,:) = nabs_surf(n1,:,:)
    !
    nx_surf(n1+1,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+2,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+3,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+4,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+5,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+6,:,:)   = nx_surf(n1,:,:)
    !
    ny_surf(n1+1,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+2,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+3,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+4,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+5,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+6,:,:)   = ny_surf(n1,:,:)
    !
    nz_surf(n1+1,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+2,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+3,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+4,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+5,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+6,:,:)   = nz_surf(n1,:,:)
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    nabs_surf(:, 0,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-1,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-2,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-3,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-4,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-5,:)   = nabs_surf(:,1,:)
    !
    nx_surf(:, 0,:)     = nx_surf(:,1,:)
    nx_surf(:,-1,:)     = nx_surf(:,1,:)
    nx_surf(:,-2,:)     = nx_surf(:,1,:)
    nx_surf(:,-3,:)     = nx_surf(:,1,:)
    nx_surf(:,-4,:)     = nx_surf(:,1,:)
    nx_surf(:,-5,:)     = nx_surf(:,1,:)
    !
    ny_surf(:, 0,:)     = ny_surf(:,1,:)
    ny_surf(:,-1,:)     = ny_surf(:,1,:)
    ny_surf(:,-2,:)     = ny_surf(:,1,:)
    ny_surf(:,-3,:)     = ny_surf(:,1,:)
    ny_surf(:,-4,:)     = ny_surf(:,1,:)
    ny_surf(:,-5,:)     = ny_surf(:,1,:)
    !
    nz_surf(:, 0,:)     = nz_surf(:,1,:)
    nz_surf(:,-1,:)     = nz_surf(:,1,:)
    nz_surf(:,-2,:)     = nz_surf(:,1,:)
    nz_surf(:,-3,:)     = nz_surf(:,1,:)
    nz_surf(:,-4,:)     = nz_surf(:,1,:)
    nz_surf(:,-5,:)     = nz_surf(:,1,:)
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    nabs_surf(:,n1+1,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+2,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+3,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+4,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+5,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+6,:) = nabs_surf(:,n1,:)
    !
    nx_surf(:,n1+1,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+2,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+3,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+4,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+5,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+6,:)   = nx_surf(:,n1,:)
    !
    ny_surf(:,n1+1,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+2,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+3,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+4,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+5,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+6,:)   = ny_surf(:,n1,:)
    !
    nz_surf(:,n1+1,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+2,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+3,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+4,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+5,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+6,:)   = nz_surf(:,n1,:)
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    nabs_surf(:,:, 0)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-1)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-2)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-3)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-4)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-5)   = nabs_surf(:,:,1)
    !
    nx_surf(:,:, 0)     = nx_surf(:,:,1)
    nx_surf(:,:,-1)     = nx_surf(:,:,1)
    nx_surf(:,:,-2)     = nx_surf(:,:,1)
    nx_surf(:,:,-3)     = nx_surf(:,:,1)
    nx_surf(:,:,-4)     = nx_surf(:,:,1)
    nx_surf(:,:,-5)     = nx_surf(:,:,1)
    !
    ny_surf(:,:, 0)     = ny_surf(:,:,1)
    ny_surf(:,:,-1)     = ny_surf(:,:,1)
    ny_surf(:,:,-2)     = ny_surf(:,:,1)
    ny_surf(:,:,-3)     = ny_surf(:,:,1)
    ny_surf(:,:,-4)     = ny_surf(:,:,1)
    ny_surf(:,:,-5)     = ny_surf(:,:,1)
    !
    nz_surf(:,:, 0)     = nz_surf(:,:,1)
    nz_surf(:,:,-1)     = nz_surf(:,:,1)
    nz_surf(:,:,-2)     = nz_surf(:,:,1)
    nz_surf(:,:,-3)     = nz_surf(:,:,1)
    nz_surf(:,:,-4)     = nz_surf(:,:,1)
    nz_surf(:,:,-5)     = nz_surf(:,:,1)
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    nabs_surf(:,:,n3+1) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+2) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+3) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+4) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+5) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+6) = nabs_surf(:,:,n3)
    !
    nx_surf(:,:,n3+1)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+2)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+3)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+4)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+5)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+6)   = nx_surf(:,:,n3)
    !
    ny_surf(:,:,n3+1)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+2)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+3)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+4)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+5)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+6)   = ny_surf(:,:,n3)
    !
    nz_surf(:,:,n3+1)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+2)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+3)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+4)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+5)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+6)   = nz_surf(:,:,n3)
  endif
  !

  if (myid.eq.0)  print*, '*** Normal vectors calculated ***'
  !*********************************************************************
  x_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  y_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  z_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  !
  call intersect(n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dl,dzc,zc,zf)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,x_intersect)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,x_intersect)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,y_intersect)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,y_intersect)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,z_intersect)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,z_intersect)
  !
  ! along x
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    x_intersect( 0,:,:)     = x_intersect(1,:,:)
    x_intersect(-1,:,:)     = x_intersect(1,:,:)
    x_intersect(-2,:,:)     = x_intersect(1,:,:)
    x_intersect(-3,:,:)     = x_intersect(1,:,:)
    x_intersect(-4,:,:)     = x_intersect(1,:,:)
    x_intersect(-5,:,:)     = x_intersect(1,:,:)
    !
    y_intersect( 0,:,:)     = y_intersect(1,:,:)
    y_intersect(-1,:,:)     = y_intersect(1,:,:)
    y_intersect(-2,:,:)     = y_intersect(1,:,:)
    y_intersect(-3,:,:)     = y_intersect(1,:,:)
    y_intersect(-4,:,:)     = y_intersect(1,:,:)
    y_intersect(-5,:,:)     = y_intersect(1,:,:)
    !
    z_intersect( 0,:,:)     = z_intersect(1,:,:)
    z_intersect(-1,:,:)     = z_intersect(1,:,:)
    z_intersect(-2,:,:)     = z_intersect(1,:,:)
    z_intersect(-3,:,:)     = z_intersect(1,:,:)
    z_intersect(-4,:,:)     = z_intersect(1,:,:)
    z_intersect(-5,:,:)     = z_intersect(1,:,:)
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    x_intersect(n1+1,:,:)   = x_intersect(n1,:,:)
    x_intersect(n1+2,:,:)   = x_intersect(n1,:,:)
    x_intersect(n1+3,:,:)   = x_intersect(n1,:,:)
    x_intersect(n1+4,:,:)   = x_intersect(n1,:,:)
    x_intersect(n1+5,:,:)   = x_intersect(n1,:,:)
    x_intersect(n1+6,:,:)   = x_intersect(n1,:,:)
    !
    y_intersect(n1+1,:,:)   = y_intersect(n1,:,:)
    y_intersect(n1+2,:,:)   = y_intersect(n1,:,:)
    y_intersect(n1+3,:,:)   = y_intersect(n1,:,:)
    y_intersect(n1+4,:,:)   = y_intersect(n1,:,:)
    y_intersect(n1+5,:,:)   = y_intersect(n1,:,:)
    y_intersect(n1+6,:,:)   = y_intersect(n1,:,:)
    !
    z_intersect(n1+1,:,:)   = z_intersect(n1,:,:)
    z_intersect(n1+2,:,:)   = z_intersect(n1,:,:)
    z_intersect(n1+3,:,:)   = z_intersect(n1,:,:)
    z_intersect(n1+4,:,:)   = z_intersect(n1,:,:)
    z_intersect(n1+5,:,:)   = z_intersect(n1,:,:)
    z_intersect(n1+6,:,:)   = z_intersect(n1,:,:)
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    x_intersect(:, 0,:)     = x_intersect(:,1,:)
    x_intersect(:,-1,:)     = x_intersect(:,1,:)
    x_intersect(:,-2,:)     = x_intersect(:,1,:)
    x_intersect(:,-3,:)     = x_intersect(:,1,:)
    x_intersect(:,-4,:)     = x_intersect(:,1,:)
    x_intersect(:,-5,:)     = x_intersect(:,1,:)
    !
    y_intersect(:, 0,:)     = y_intersect(:,1,:)
    y_intersect(:,-1,:)     = y_intersect(:,1,:)
    y_intersect(:,-2,:)     = y_intersect(:,1,:)
    y_intersect(:,-3,:)     = y_intersect(:,1,:)
    y_intersect(:,-4,:)     = y_intersect(:,1,:)
    y_intersect(:,-5,:)     = y_intersect(:,1,:)
    !
    z_intersect(:, 0,:)     = z_intersect(:,1,:)
    z_intersect(:,-1,:)     = z_intersect(:,1,:)
    z_intersect(:,-2,:)     = z_intersect(:,1,:)
    z_intersect(:,-3,:)     = z_intersect(:,1,:)
    z_intersect(:,-4,:)     = z_intersect(:,1,:)
    z_intersect(:,-5,:)     = z_intersect(:,1,:)
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    x_intersect(:,n1+1,:)   = x_intersect(:,n1,:)
    x_intersect(:,n1+2,:)   = x_intersect(:,n1,:)
    x_intersect(:,n1+3,:)   = x_intersect(:,n1,:)
    x_intersect(:,n1+4,:)   = x_intersect(:,n1,:)
    x_intersect(:,n1+5,:)   = x_intersect(:,n1,:)
    x_intersect(:,n1+6,:)   = x_intersect(:,n1,:)
    !               
    y_intersect(:,n1+1,:)   = y_intersect(:,n1,:)
    y_intersect(:,n1+2,:)   = y_intersect(:,n1,:)
    y_intersect(:,n1+3,:)   = y_intersect(:,n1,:)
    y_intersect(:,n1+4,:)   = y_intersect(:,n1,:)
    y_intersect(:,n1+5,:)   = y_intersect(:,n1,:)
    y_intersect(:,n1+6,:)   = y_intersect(:,n1,:)
    !               
    z_intersect(:,n1+1,:)   = z_intersect(:,n1,:)
    z_intersect(:,n1+2,:)   = z_intersect(:,n1,:)
    z_intersect(:,n1+3,:)   = z_intersect(:,n1,:)
    z_intersect(:,n1+4,:)   = z_intersect(:,n1,:)
    z_intersect(:,n1+5,:)   = z_intersect(:,n1,:)
    z_intersect(:,n1+6,:)   = z_intersect(:,n1,:)
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    x_intersect(:,:, 0)     = x_intersect(:,:,1)
    x_intersect(:,:,-1)     = x_intersect(:,:,1)
    x_intersect(:,:,-2)     = x_intersect(:,:,1)
    x_intersect(:,:,-3)     = x_intersect(:,:,1)
    x_intersect(:,:,-4)     = x_intersect(:,:,1)
    x_intersect(:,:,-5)     = x_intersect(:,:,1)
    !                         
    y_intersect(:,:, 0)     = y_intersect(:,:,1)
    y_intersect(:,:,-1)     = y_intersect(:,:,1)
    y_intersect(:,:,-2)     = y_intersect(:,:,1)
    y_intersect(:,:,-3)     = y_intersect(:,:,1)
    y_intersect(:,:,-4)     = y_intersect(:,:,1)
    y_intersect(:,:,-5)     = y_intersect(:,:,1)
    !                         
    z_intersect(:,:, 0)     = z_intersect(:,:,1)
    z_intersect(:,:,-1)     = z_intersect(:,:,1)
    z_intersect(:,:,-2)     = z_intersect(:,:,1)
    z_intersect(:,:,-3)     = z_intersect(:,:,1)
    z_intersect(:,:,-4)     = z_intersect(:,:,1)
    z_intersect(:,:,-5)     = z_intersect(:,:,1)
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    x_intersect(:,:,n3+1)   = x_intersect(:,:,n3)
    x_intersect(:,:,n3+2)   = x_intersect(:,:,n3)
    x_intersect(:,:,n3+3)   = x_intersect(:,:,n3)
    x_intersect(:,:,n3+4)   = x_intersect(:,:,n3)
    x_intersect(:,:,n3+5)   = x_intersect(:,:,n3)
    x_intersect(:,:,n3+6)   = x_intersect(:,:,n3)
    !                         
    y_intersect(:,:,n3+1)   = y_intersect(:,:,n3)
    y_intersect(:,:,n3+2)   = y_intersect(:,:,n3)
    y_intersect(:,:,n3+3)   = y_intersect(:,:,n3)
    y_intersect(:,:,n3+4)   = y_intersect(:,:,n3)
    y_intersect(:,:,n3+5)   = y_intersect(:,:,n3)
    y_intersect(:,:,n3+6)   = y_intersect(:,:,n3)
    !                         
    z_intersect(:,:,n3+1)   = z_intersect(:,:,n3)
    z_intersect(:,:,n3+2)   = z_intersect(:,:,n3)
    z_intersect(:,:,n3+3)   = z_intersect(:,:,n3)
    z_intersect(:,:,n3+4)   = z_intersect(:,:,n3)
    z_intersect(:,:,n3+5)   = z_intersect(:,:,n3)
    z_intersect(:,:,n3+6)   = z_intersect(:,:,n3)
  endif
  !
  
  if (myid.eq.0)  print*, '*** Intersection points calculated ***'
  !*********************************************************************
  x_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = -1000.0_rp
  y_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = -1000.0_rp
  z_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = -1000.0_rp
  x_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  y_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  z_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  x_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  y_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  z_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = -1000.0_rp
  deltan(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)   = -1000.0_rp
  !
  call mirrorpoints(n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf, &
                    x_intersect,y_intersect,z_intersect, &
                    x_mirror,y_mirror,z_mirror, &
                    deltan, &
                    x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                    dl,dzc,zc)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,deltan)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,deltan)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,x_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,x_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,y_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,y_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,z_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,z_mirror)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,x_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,x_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,y_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,y_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,z_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,z_IP1)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,x_IP2)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,x_IP2)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,y_IP2)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,y_IP2)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,z_IP2)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,z_IP2)
  !
  ! along x
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    deltan( 0,:,:)     = deltan(1,:,:)
    deltan(-1,:,:)     = deltan(1,:,:)
    deltan(-2,:,:)     = deltan(1,:,:)
    deltan(-3,:,:)     = deltan(1,:,:)
    deltan(-4,:,:)     = deltan(1,:,:)
    deltan(-5,:,:)     = deltan(1,:,:)
    !
    x_mirror( 0,:,:)   = x_mirror(1,:,:)
    x_mirror(-1,:,:)   = x_mirror(1,:,:)
    x_mirror(-2,:,:)   = x_mirror(1,:,:)
    x_mirror(-3,:,:)   = x_mirror(1,:,:)
    x_mirror(-4,:,:)   = x_mirror(1,:,:)
    x_mirror(-5,:,:)   = x_mirror(1,:,:)
    !
    y_mirror( 0,:,:)   = y_mirror(1,:,:)
    y_mirror(-1,:,:)   = y_mirror(1,:,:)
    y_mirror(-2,:,:)   = y_mirror(1,:,:)
    y_mirror(-3,:,:)   = y_mirror(1,:,:)
    y_mirror(-4,:,:)   = y_mirror(1,:,:)
    y_mirror(-5,:,:)   = y_mirror(1,:,:)
    !
    z_mirror( 0,:,:)   = z_mirror(1,:,:)
    z_mirror(-1,:,:)   = z_mirror(1,:,:)
    z_mirror(-2,:,:)   = z_mirror(1,:,:)
    z_mirror(-3,:,:)   = z_mirror(1,:,:)
    z_mirror(-4,:,:)   = z_mirror(1,:,:)
    z_mirror(-5,:,:)   = z_mirror(1,:,:)
    !
    x_IP1( 0,:,:)      = x_IP1(1,:,:)
    x_IP1(-1,:,:)      = x_IP1(1,:,:)
    x_IP1(-2,:,:)      = x_IP1(1,:,:)
    x_IP1(-3,:,:)      = x_IP1(1,:,:)
    x_IP1(-4,:,:)      = x_IP1(1,:,:)
    x_IP1(-5,:,:)      = x_IP1(1,:,:)
    !                  
    y_IP1( 0,:,:)      = y_IP1(1,:,:)
    y_IP1(-1,:,:)      = y_IP1(1,:,:)
    y_IP1(-2,:,:)      = y_IP1(1,:,:)
    y_IP1(-3,:,:)      = y_IP1(1,:,:)
    y_IP1(-4,:,:)      = y_IP1(1,:,:)
    y_IP1(-5,:,:)      = y_IP1(1,:,:)
    !                  
    z_IP1( 0,:,:)      = z_IP1(1,:,:)
    z_IP1(-1,:,:)      = z_IP1(1,:,:)
    z_IP1(-2,:,:)      = z_IP1(1,:,:)
    z_IP1(-3,:,:)      = z_IP1(1,:,:)
    z_IP1(-4,:,:)      = z_IP1(1,:,:)
    z_IP1(-5,:,:)      = z_IP1(1,:,:)
    !                  
    x_IP2( 0,:,:)      = x_IP2(1,:,:)
    x_IP2(-1,:,:)      = x_IP2(1,:,:)
    x_IP2(-2,:,:)      = x_IP2(1,:,:)
    x_IP2(-3,:,:)      = x_IP2(1,:,:)
    x_IP2(-4,:,:)      = x_IP2(1,:,:)
    x_IP2(-5,:,:)      = x_IP2(1,:,:)
    !                  
    y_IP2( 0,:,:)      = y_IP2(1,:,:)
    y_IP2(-1,:,:)      = y_IP2(1,:,:)
    y_IP2(-2,:,:)      = y_IP2(1,:,:)
    y_IP2(-3,:,:)      = y_IP2(1,:,:)
    y_IP2(-4,:,:)      = y_IP2(1,:,:)
    y_IP2(-5,:,:)      = y_IP2(1,:,:)
    !                  
    z_IP2( 0,:,:)      = z_IP2(1,:,:)
    z_IP2(-1,:,:)      = z_IP2(1,:,:)
    z_IP2(-2,:,:)      = z_IP2(1,:,:)
    z_IP2(-3,:,:)      = z_IP2(1,:,:)
    z_IP2(-4,:,:)      = z_IP2(1,:,:)
    z_IP2(-5,:,:)      = z_IP2(1,:,:)
    !
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    deltan(n1+1,:,:)   = deltan(n1,:,:)
    deltan(n1+2,:,:)   = deltan(n1,:,:)
    deltan(n1+3,:,:)   = deltan(n1,:,:)
    deltan(n1+4,:,:)   = deltan(n1,:,:)
    deltan(n1+5,:,:)   = deltan(n1,:,:)
    deltan(n1+6,:,:)   = deltan(n1,:,:)
    !
    x_mirror(n1+1,:,:)   = x_mirror(n1,:,:)
    x_mirror(n1+2,:,:)   = x_mirror(n1,:,:)
    x_mirror(n1+3,:,:)   = x_mirror(n1,:,:)
    x_mirror(n1+4,:,:)   = x_mirror(n1,:,:)
    x_mirror(n1+5,:,:)   = x_mirror(n1,:,:)
    x_mirror(n1+6,:,:)   = x_mirror(n1,:,:)
    !
    y_mirror(n1+1,:,:)   = y_mirror(n1,:,:)
    y_mirror(n1+2,:,:)   = y_mirror(n1,:,:)
    y_mirror(n1+3,:,:)   = y_mirror(n1,:,:)
    y_mirror(n1+4,:,:)   = y_mirror(n1,:,:)
    y_mirror(n1+5,:,:)   = y_mirror(n1,:,:)
    y_mirror(n1+6,:,:)   = y_mirror(n1,:,:)
    !
    z_mirror(n1+1,:,:)   = z_mirror(n1,:,:)
    z_mirror(n1+2,:,:)   = z_mirror(n1,:,:)
    z_mirror(n1+3,:,:)   = z_mirror(n1,:,:)
    z_mirror(n1+4,:,:)   = z_mirror(n1,:,:)
    z_mirror(n1+5,:,:)   = z_mirror(n1,:,:)
    z_mirror(n1+6,:,:)   = z_mirror(n1,:,:)
    !
    x_IP1(n1+1,:,:)   = x_IP1(n1,:,:)
    x_IP1(n1+2,:,:)   = x_IP1(n1,:,:)
    x_IP1(n1+3,:,:)   = x_IP1(n1,:,:)
    x_IP1(n1+4,:,:)   = x_IP1(n1,:,:)
    x_IP1(n1+5,:,:)   = x_IP1(n1,:,:)
    x_IP1(n1+6,:,:)   = x_IP1(n1,:,:)
    !
    y_IP1(n1+1,:,:)   = y_IP1(n1,:,:)
    y_IP1(n1+2,:,:)   = y_IP1(n1,:,:)
    y_IP1(n1+3,:,:)   = y_IP1(n1,:,:)
    y_IP1(n1+4,:,:)   = y_IP1(n1,:,:)
    y_IP1(n1+5,:,:)   = y_IP1(n1,:,:)
    y_IP1(n1+6,:,:)   = y_IP1(n1,:,:)
    !
    z_IP1(n1+1,:,:)   = z_IP1(n1,:,:)
    z_IP1(n1+2,:,:)   = z_IP1(n1,:,:)
    z_IP1(n1+3,:,:)   = z_IP1(n1,:,:)
    z_IP1(n1+4,:,:)   = z_IP1(n1,:,:)
    z_IP1(n1+5,:,:)   = z_IP1(n1,:,:)
    z_IP1(n1+6,:,:)   = z_IP1(n1,:,:)
    !
    x_IP2(n1+1,:,:)   = x_IP2(n1,:,:)
    x_IP2(n1+2,:,:)   = x_IP2(n1,:,:)
    x_IP2(n1+3,:,:)   = x_IP2(n1,:,:)
    x_IP2(n1+4,:,:)   = x_IP2(n1,:,:)
    x_IP2(n1+5,:,:)   = x_IP2(n1,:,:)
    x_IP2(n1+6,:,:)   = x_IP2(n1,:,:)
    !
    y_IP2(n1+1,:,:)   = y_IP2(n1,:,:)
    y_IP2(n1+2,:,:)   = y_IP2(n1,:,:)
    y_IP2(n1+3,:,:)   = y_IP2(n1,:,:)
    y_IP2(n1+4,:,:)   = y_IP2(n1,:,:)
    y_IP2(n1+5,:,:)   = y_IP2(n1,:,:)
    y_IP2(n1+6,:,:)   = y_IP2(n1,:,:)
    !
    z_IP2(n1+1,:,:)   = z_IP2(n1,:,:)
    z_IP2(n1+2,:,:)   = z_IP2(n1,:,:)
    z_IP2(n1+3,:,:)   = z_IP2(n1,:,:)
    z_IP2(n1+4,:,:)   = z_IP2(n1,:,:)
    z_IP2(n1+5,:,:)   = z_IP2(n1,:,:)
    z_IP2(n1+6,:,:)   = z_IP2(n1,:,:)
    !
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    deltan(:, 0,:)     = deltan(:,1,:)
    deltan(:,-1,:)     = deltan(:,1,:)
    deltan(:,-2,:)     = deltan(:,1,:)
    deltan(:,-3,:)     = deltan(:,1,:)
    deltan(:,-4,:)     = deltan(:,1,:)
    deltan(:,-5,:)     = deltan(:,1,:)
    !
    x_mirror(:, 0,:)     = x_mirror(:,1,:)
    x_mirror(:,-1,:)     = x_mirror(:,1,:)
    x_mirror(:,-2,:)     = x_mirror(:,1,:)
    x_mirror(:,-3,:)     = x_mirror(:,1,:)
    x_mirror(:,-4,:)     = x_mirror(:,1,:)
    x_mirror(:,-5,:)     = x_mirror(:,1,:)
    !
    y_mirror(:, 0,:)     = y_mirror(:,1,:)
    y_mirror(:,-1,:)     = y_mirror(:,1,:)
    y_mirror(:,-2,:)     = y_mirror(:,1,:)
    y_mirror(:,-3,:)     = y_mirror(:,1,:)
    y_mirror(:,-4,:)     = y_mirror(:,1,:)
    y_mirror(:,-5,:)     = y_mirror(:,1,:)
    !
    z_mirror(:, 0,:)     = z_mirror(:,1,:)
    z_mirror(:,-1,:)     = z_mirror(:,1,:)
    z_mirror(:,-2,:)     = z_mirror(:,1,:)
    z_mirror(:,-3,:)     = z_mirror(:,1,:)
    z_mirror(:,-4,:)     = z_mirror(:,1,:)
    z_mirror(:,-5,:)     = z_mirror(:,1,:)
    !
    x_IP1(:, 0,:)     = x_IP1(:,1,:)
    x_IP1(:,-1,:)     = x_IP1(:,1,:)
    x_IP1(:,-2,:)     = x_IP1(:,1,:)
    x_IP1(:,-3,:)     = x_IP1(:,1,:)
    x_IP1(:,-4,:)     = x_IP1(:,1,:)
    x_IP1(:,-5,:)     = x_IP1(:,1,:)
    !
    y_IP1(:, 0,:)     = y_IP1(:,1,:)
    y_IP1(:,-1,:)     = y_IP1(:,1,:)
    y_IP1(:,-2,:)     = y_IP1(:,1,:)
    y_IP1(:,-3,:)     = y_IP1(:,1,:)
    y_IP1(:,-4,:)     = y_IP1(:,1,:)
    y_IP1(:,-5,:)     = y_IP1(:,1,:)
    !
    z_IP1(:, 0,:)     = z_IP1(:,1,:)
    z_IP1(:,-1,:)     = z_IP1(:,1,:)
    z_IP1(:,-2,:)     = z_IP1(:,1,:)
    z_IP1(:,-3,:)     = z_IP1(:,1,:)
    z_IP1(:,-4,:)     = z_IP1(:,1,:)
    z_IP1(:,-5,:)     = z_IP1(:,1,:)
    !
    x_IP2(:, 0,:)     = x_IP2(:,1,:)
    x_IP2(:,-1,:)     = x_IP2(:,1,:)
    x_IP2(:,-2,:)     = x_IP2(:,1,:)
    x_IP2(:,-3,:)     = x_IP2(:,1,:)
    x_IP2(:,-4,:)     = x_IP2(:,1,:)
    x_IP2(:,-5,:)     = x_IP2(:,1,:)
    !
    y_IP2(:, 0,:)     = y_IP2(:,1,:)
    y_IP2(:,-1,:)     = y_IP2(:,1,:)
    y_IP2(:,-2,:)     = y_IP2(:,1,:)
    y_IP2(:,-3,:)     = y_IP2(:,1,:)
    y_IP2(:,-4,:)     = y_IP2(:,1,:)
    y_IP2(:,-5,:)     = y_IP2(:,1,:)
    !
    z_IP2(:, 0,:)     = z_IP2(:,1,:)
    z_IP2(:,-1,:)     = z_IP2(:,1,:)
    z_IP2(:,-2,:)     = z_IP2(:,1,:)
    z_IP2(:,-3,:)     = z_IP2(:,1,:)
    z_IP2(:,-4,:)     = z_IP2(:,1,:)
    z_IP2(:,-5,:)     = z_IP2(:,1,:)
    !
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    deltan(:,n1+1,:)   = deltan(:,n1,:)
    deltan(:,n1+2,:)   = deltan(:,n1,:)
    deltan(:,n1+3,:)   = deltan(:,n1,:)
    deltan(:,n1+4,:)   = deltan(:,n1,:)
    deltan(:,n1+5,:)   = deltan(:,n1,:)
    deltan(:,n1+6,:)   = deltan(:,n1,:)
    !               
    x_mirror(:,n1+1,:)   = x_mirror(:,n1,:)
    x_mirror(:,n1+2,:)   = x_mirror(:,n1,:)
    x_mirror(:,n1+3,:)   = x_mirror(:,n1,:)
    x_mirror(:,n1+4,:)   = x_mirror(:,n1,:)
    x_mirror(:,n1+5,:)   = x_mirror(:,n1,:)
    x_mirror(:,n1+6,:)   = x_mirror(:,n1,:)
    !               
    y_mirror(:,n1+1,:)   = y_mirror(:,n1,:)
    y_mirror(:,n1+2,:)   = y_mirror(:,n1,:)
    y_mirror(:,n1+3,:)   = y_mirror(:,n1,:)
    y_mirror(:,n1+4,:)   = y_mirror(:,n1,:)
    y_mirror(:,n1+5,:)   = y_mirror(:,n1,:)
    y_mirror(:,n1+6,:)   = y_mirror(:,n1,:)
    !
    z_mirror(:,n1+1,:)   = z_mirror(:,n1,:)
    z_mirror(:,n1+2,:)   = z_mirror(:,n1,:)
    z_mirror(:,n1+3,:)   = z_mirror(:,n1,:)
    z_mirror(:,n1+4,:)   = z_mirror(:,n1,:)
    z_mirror(:,n1+5,:)   = z_mirror(:,n1,:)
    z_mirror(:,n1+6,:)   = z_mirror(:,n1,:)
    !               
    x_IP1(:,n1+1,:)   = x_IP1(:,n1,:)
    x_IP1(:,n1+2,:)   = x_IP1(:,n1,:)
    x_IP1(:,n1+3,:)   = x_IP1(:,n1,:)
    x_IP1(:,n1+4,:)   = x_IP1(:,n1,:)
    x_IP1(:,n1+5,:)   = x_IP1(:,n1,:)
    x_IP1(:,n1+6,:)   = x_IP1(:,n1,:)
    !               
    y_IP1(:,n1+1,:)   = y_IP1(:,n1,:)
    y_IP1(:,n1+2,:)   = y_IP1(:,n1,:)
    y_IP1(:,n1+3,:)   = y_IP1(:,n1,:)
    y_IP1(:,n1+4,:)   = y_IP1(:,n1,:)
    y_IP1(:,n1+5,:)   = y_IP1(:,n1,:)
    y_IP1(:,n1+6,:)   = y_IP1(:,n1,:)
    !
    z_IP1(:,n1+1,:)   = z_IP1(:,n1,:)
    z_IP1(:,n1+2,:)   = z_IP1(:,n1,:)
    z_IP1(:,n1+3,:)   = z_IP1(:,n1,:)
    z_IP1(:,n1+4,:)   = z_IP1(:,n1,:)
    z_IP1(:,n1+5,:)   = z_IP1(:,n1,:)
    z_IP1(:,n1+6,:)   = z_IP1(:,n1,:)
    !               
    x_IP2(:,n1+1,:)   = x_IP2(:,n1,:)
    x_IP2(:,n1+2,:)   = x_IP2(:,n1,:)
    x_IP2(:,n1+3,:)   = x_IP2(:,n1,:)
    x_IP2(:,n1+4,:)   = x_IP2(:,n1,:)
    x_IP2(:,n1+5,:)   = x_IP2(:,n1,:)
    x_IP2(:,n1+6,:)   = x_IP2(:,n1,:)
    !
    y_IP2(:,n1+1,:)   = y_IP2(:,n1,:)
    y_IP2(:,n1+2,:)   = y_IP2(:,n1,:)
    y_IP2(:,n1+3,:)   = y_IP2(:,n1,:)
    y_IP2(:,n1+4,:)   = y_IP2(:,n1,:)
    y_IP2(:,n1+5,:)   = y_IP2(:,n1,:)
    y_IP2(:,n1+6,:)   = y_IP2(:,n1,:)
    !
    z_IP2(:,n1+1,:)   = z_IP2(:,n1,:)
    z_IP2(:,n1+2,:)   = z_IP2(:,n1,:)
    z_IP2(:,n1+3,:)   = z_IP2(:,n1,:)
    z_IP2(:,n1+4,:)   = z_IP2(:,n1,:)
    z_IP2(:,n1+5,:)   = z_IP2(:,n1,:)
    z_IP2(:,n1+6,:)   = z_IP2(:,n1,:)
    ! 
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    deltan(:,:, 0)     = deltan(:,:,1)
    deltan(:,:,-1)     = deltan(:,:,1)
    deltan(:,:,-2)     = deltan(:,:,1)
    deltan(:,:,-3)     = deltan(:,:,1)
    deltan(:,:,-4)     = deltan(:,:,1)
    deltan(:,:,-5)     = deltan(:,:,1)
    !                         
    x_mirror(:,:, 0)     = x_mirror(:,:,1)
    x_mirror(:,:,-1)     = x_mirror(:,:,1)
    x_mirror(:,:,-2)     = x_mirror(:,:,1)
    x_mirror(:,:,-3)     = x_mirror(:,:,1)
    x_mirror(:,:,-4)     = x_mirror(:,:,1)
    x_mirror(:,:,-5)     = x_mirror(:,:,1)
    !                         
    y_mirror(:,:, 0)     = y_mirror(:,:,1)
    y_mirror(:,:,-1)     = y_mirror(:,:,1)
    y_mirror(:,:,-2)     = y_mirror(:,:,1)
    y_mirror(:,:,-3)     = y_mirror(:,:,1)
    y_mirror(:,:,-4)     = y_mirror(:,:,1)
    y_mirror(:,:,-5)     = y_mirror(:,:,1)
    !
    z_mirror(:,:, 0)     = z_mirror(:,:,1)
    z_mirror(:,:,-1)     = z_mirror(:,:,1)
    z_mirror(:,:,-2)     = z_mirror(:,:,1)
    z_mirror(:,:,-3)     = z_mirror(:,:,1)
    z_mirror(:,:,-4)     = z_mirror(:,:,1)
    z_mirror(:,:,-5)     = z_mirror(:,:,1)
    !                         
    x_IP1(:,:, 0)     = x_IP1(:,:,1)
    x_IP1(:,:,-1)     = x_IP1(:,:,1)
    x_IP1(:,:,-2)     = x_IP1(:,:,1)
    x_IP1(:,:,-3)     = x_IP1(:,:,1)
    x_IP1(:,:,-4)     = x_IP1(:,:,1)
    x_IP1(:,:,-5)     = x_IP1(:,:,1)
    !                         
    y_IP1(:,:, 0)     = y_IP1(:,:,1)
    y_IP1(:,:,-1)     = y_IP1(:,:,1)
    y_IP1(:,:,-2)     = y_IP1(:,:,1)
    y_IP1(:,:,-3)     = y_IP1(:,:,1)
    y_IP1(:,:,-4)     = y_IP1(:,:,1)
    y_IP1(:,:,-5)     = y_IP1(:,:,1)
    !
    z_IP1(:,:, 0)     = z_IP1(:,:,1)
    z_IP1(:,:,-1)     = z_IP1(:,:,1)
    z_IP1(:,:,-2)     = z_IP1(:,:,1)
    z_IP1(:,:,-3)     = z_IP1(:,:,1)
    z_IP1(:,:,-4)     = z_IP1(:,:,1)
    z_IP1(:,:,-5)     = z_IP1(:,:,1)
    !                         
    x_IP2(:,:, 0)     = x_IP2(:,:,1)
    x_IP2(:,:,-1)     = x_IP2(:,:,1)
    x_IP2(:,:,-2)     = x_IP2(:,:,1)
    x_IP2(:,:,-3)     = x_IP2(:,:,1)
    x_IP2(:,:,-4)     = x_IP2(:,:,1)
    x_IP2(:,:,-5)     = x_IP2(:,:,1)
    !
    y_IP2(:,:, 0)     = y_IP2(:,:,1)
    y_IP2(:,:,-1)     = y_IP2(:,:,1)
    y_IP2(:,:,-2)     = y_IP2(:,:,1)
    y_IP2(:,:,-3)     = y_IP2(:,:,1)
    y_IP2(:,:,-4)     = y_IP2(:,:,1)
    y_IP2(:,:,-5)     = y_IP2(:,:,1)
    !
    z_IP2(:,:, 0)     = z_IP2(:,:,1)
    z_IP2(:,:,-1)     = z_IP2(:,:,1)
    z_IP2(:,:,-2)     = z_IP2(:,:,1)
    z_IP2(:,:,-3)     = z_IP2(:,:,1)
    z_IP2(:,:,-4)     = z_IP2(:,:,1)
    z_IP2(:,:,-5)     = z_IP2(:,:,1)
    ! 
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    deltan(:,:,n3+1)   = deltan(:,:,n3)
    deltan(:,:,n3+2)   = deltan(:,:,n3)
    deltan(:,:,n3+3)   = deltan(:,:,n3)
    deltan(:,:,n3+4)   = deltan(:,:,n3)
    deltan(:,:,n3+5)   = deltan(:,:,n3)
    deltan(:,:,n3+6)   = deltan(:,:,n3)
    !                         
    x_mirror(:,:,n3+1)   = x_mirror(:,:,n3)
    x_mirror(:,:,n3+2)   = x_mirror(:,:,n3)
    x_mirror(:,:,n3+3)   = x_mirror(:,:,n3)
    x_mirror(:,:,n3+4)   = x_mirror(:,:,n3)
    x_mirror(:,:,n3+5)   = x_mirror(:,:,n3)
    x_mirror(:,:,n3+6)   = x_mirror(:,:,n3)
    !                         
    y_mirror(:,:,n3+1)   = y_mirror(:,:,n3)
    y_mirror(:,:,n3+2)   = y_mirror(:,:,n3)
    y_mirror(:,:,n3+3)   = y_mirror(:,:,n3)
    y_mirror(:,:,n3+4)   = y_mirror(:,:,n3)
    y_mirror(:,:,n3+5)   = y_mirror(:,:,n3)
    y_mirror(:,:,n3+6)   = y_mirror(:,:,n3)
    !
    z_mirror(:,:,n3+1)   = z_mirror(:,:,n3)
    z_mirror(:,:,n3+2)   = z_mirror(:,:,n3)
    z_mirror(:,:,n3+3)   = z_mirror(:,:,n3)
    z_mirror(:,:,n3+4)   = z_mirror(:,:,n3)
    z_mirror(:,:,n3+5)   = z_mirror(:,:,n3)
    z_mirror(:,:,n3+6)   = z_mirror(:,:,n3)
    !                      
    x_IP1(:,:,n3+1)   = x_IP1(:,:,n3)
    x_IP1(:,:,n3+2)   = x_IP1(:,:,n3)
    x_IP1(:,:,n3+3)   = x_IP1(:,:,n3)
    x_IP1(:,:,n3+4)   = x_IP1(:,:,n3)
    x_IP1(:,:,n3+5)   = x_IP1(:,:,n3)
    x_IP1(:,:,n3+6)   = x_IP1(:,:,n3)
    !                         
    y_IP1(:,:,n3+1)   = y_IP1(:,:,n3)
    y_IP1(:,:,n3+2)   = y_IP1(:,:,n3)
    y_IP1(:,:,n3+3)   = y_IP1(:,:,n3)
    y_IP1(:,:,n3+4)   = y_IP1(:,:,n3)
    y_IP1(:,:,n3+5)   = y_IP1(:,:,n3)
    y_IP1(:,:,n3+6)   = y_IP1(:,:,n3)
    !
    z_IP1(:,:,n3+1)   = z_IP1(:,:,n3)
    z_IP1(:,:,n3+2)   = z_IP1(:,:,n3)
    z_IP1(:,:,n3+3)   = z_IP1(:,:,n3)
    z_IP1(:,:,n3+4)   = z_IP1(:,:,n3)
    z_IP1(:,:,n3+5)   = z_IP1(:,:,n3)
    z_IP1(:,:,n3+6)   = z_IP1(:,:,n3)
    !                      
    x_IP2(:,:,n3+1)   = x_IP2(:,:,n3)
    x_IP2(:,:,n3+2)   = x_IP2(:,:,n3)
    x_IP2(:,:,n3+3)   = x_IP2(:,:,n3)
    x_IP2(:,:,n3+4)   = x_IP2(:,:,n3)
    x_IP2(:,:,n3+5)   = x_IP2(:,:,n3)
    x_IP2(:,:,n3+6)   = x_IP2(:,:,n3)
    !
    y_IP2(:,:,n3+1)   = y_IP2(:,:,n3)
    y_IP2(:,:,n3+2)   = y_IP2(:,:,n3)
    y_IP2(:,:,n3+3)   = y_IP2(:,:,n3)
    y_IP2(:,:,n3+4)   = y_IP2(:,:,n3)
    y_IP2(:,:,n3+5)   = y_IP2(:,:,n3)
    y_IP2(:,:,n3+6)   = y_IP2(:,:,n3)
    !
    z_IP2(:,:,n3+1)   = z_IP2(:,:,n3)
    z_IP2(:,:,n3+2)   = z_IP2(:,:,n3)
    z_IP2(:,:,n3+3)   = z_IP2(:,:,n3)
    z_IP2(:,:,n3+4)   = z_IP2(:,:,n3)
    z_IP2(:,:,n3+5)   = z_IP2(:,:,n3)
    z_IP2(:,:,n3+6)   = z_IP2(:,:,n3)
    !
  endif
  !
  i_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
  j_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
  k_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
  i_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  j_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  k_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  i_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  j_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  k_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)       =   -1000
  !
  if (myid.eq.0)  print*, '*** Mirror points set ***'
  !*********************************************************************
  call mirrorpoints_ijk(n,ng,nh_d,nh_b,nabs_surf,x_mirror,y_mirror,z_mirror, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        i_mirror,j_mirror,k_mirror, &
                        i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                        dl,zf)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    i_mirror( 0,:,:)  = i_mirror(1,:,:)
    i_mirror(-1,:,:)  = i_mirror(1,:,:)
    i_mirror(-2,:,:)  = i_mirror(1,:,:)
    i_mirror(-3,:,:)  = i_mirror(1,:,:)
    i_mirror(-4,:,:)  = i_mirror(1,:,:)
    i_mirror(-5,:,:)  = i_mirror(1,:,:)
    !
    j_mirror( 0,:,:)  = j_mirror(1,:,:)
    j_mirror(-1,:,:)  = j_mirror(1,:,:)
    j_mirror(-2,:,:)  = j_mirror(1,:,:)
    j_mirror(-3,:,:)  = j_mirror(1,:,:)
    j_mirror(-4,:,:)  = j_mirror(1,:,:)
    j_mirror(-5,:,:)  = j_mirror(1,:,:)
    !
    k_mirror( 0,:,:)  = k_mirror(1,:,:)
    k_mirror(-1,:,:)  = k_mirror(1,:,:)
    k_mirror(-2,:,:)  = k_mirror(1,:,:)
    k_mirror(-3,:,:)  = k_mirror(1,:,:)
    k_mirror(-4,:,:)  = k_mirror(1,:,:)
    k_mirror(-5,:,:)  = k_mirror(1,:,:)
    !
    i_IP1( 0,:,:)     = i_IP1(1,:,:)
    i_IP1(-1,:,:)     = i_IP1(1,:,:)
    i_IP1(-2,:,:)     = i_IP1(1,:,:)
    i_IP1(-3,:,:)     = i_IP1(1,:,:)
    i_IP1(-4,:,:)     = i_IP1(1,:,:)
    i_IP1(-5,:,:)     = i_IP1(1,:,:)
    !
    j_IP1( 0,:,:)     = j_IP1(1,:,:)
    j_IP1(-1,:,:)     = j_IP1(1,:,:)
    j_IP1(-2,:,:)     = j_IP1(1,:,:)
    j_IP1(-3,:,:)     = j_IP1(1,:,:)
    j_IP1(-4,:,:)     = j_IP1(1,:,:)
    j_IP1(-5,:,:)     = j_IP1(1,:,:)
    !
    k_IP1( 0,:,:)     = k_IP1(1,:,:)
    k_IP1(-1,:,:)     = k_IP1(1,:,:)
    k_IP1(-2,:,:)     = k_IP1(1,:,:)
    k_IP1(-3,:,:)     = k_IP1(1,:,:)
    k_IP1(-4,:,:)     = k_IP1(1,:,:)
    k_IP1(-5,:,:)     = k_IP1(1,:,:)
    !
    i_IP2( 0,:,:)     = i_IP2(1,:,:)
    i_IP2(-1,:,:)     = i_IP2(1,:,:)
    i_IP2(-2,:,:)     = i_IP2(1,:,:)
    i_IP2(-3,:,:)     = i_IP2(1,:,:)
    i_IP2(-4,:,:)     = i_IP2(1,:,:)
    i_IP2(-5,:,:)     = i_IP2(1,:,:)
    !
    j_IP2( 0,:,:)     = j_IP2(1,:,:)
    j_IP2(-1,:,:)     = j_IP2(1,:,:)
    j_IP2(-2,:,:)     = j_IP2(1,:,:)
    j_IP2(-3,:,:)     = j_IP2(1,:,:)
    j_IP2(-4,:,:)     = j_IP2(1,:,:)
    j_IP2(-5,:,:)     = j_IP2(1,:,:)
    !
    k_IP2( 0,:,:)     = k_IP2(1,:,:)
    k_IP2(-1,:,:)     = k_IP2(1,:,:)
    k_IP2(-2,:,:)     = k_IP2(1,:,:)
    k_IP2(-3,:,:)     = k_IP2(1,:,:)
    k_IP2(-4,:,:)     = k_IP2(1,:,:)
    k_IP2(-5,:,:)     = k_IP2(1,:,:)
    !
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    i_mirror(n1+1,:,:) = i_mirror(n1,:,:)
    i_mirror(n1+2,:,:) = i_mirror(n1,:,:)
    i_mirror(n1+3,:,:) = i_mirror(n1,:,:)
    i_mirror(n1+4,:,:) = i_mirror(n1,:,:)
    i_mirror(n1+5,:,:) = i_mirror(n1,:,:)
    i_mirror(n1+6,:,:) = i_mirror(n1,:,:)
    !                    
    j_mirror(n1+1,:,:) = j_mirror(n1,:,:)
    j_mirror(n1+2,:,:) = j_mirror(n1,:,:)
    j_mirror(n1+3,:,:) = j_mirror(n1,:,:)
    j_mirror(n1+4,:,:) = j_mirror(n1,:,:)
    j_mirror(n1+5,:,:) = j_mirror(n1,:,:)
    j_mirror(n1+6,:,:) = j_mirror(n1,:,:)
    !                    
    k_mirror(n1+1,:,:) = k_mirror(n1,:,:)
    k_mirror(n1+2,:,:) = k_mirror(n1,:,:)
    k_mirror(n1+3,:,:) = k_mirror(n1,:,:)
    k_mirror(n1+4,:,:) = k_mirror(n1,:,:)
    k_mirror(n1+5,:,:) = k_mirror(n1,:,:)
    k_mirror(n1+6,:,:) = k_mirror(n1,:,:)
    !
    i_IP1(n1+1,:,:)   = i_IP1(n1,:,:)
    i_IP1(n1+2,:,:)   = i_IP1(n1,:,:)
    i_IP1(n1+3,:,:)   = i_IP1(n1,:,:)
    i_IP1(n1+4,:,:)   = i_IP1(n1,:,:)
    i_IP1(n1+5,:,:)   = i_IP1(n1,:,:)
    i_IP1(n1+6,:,:)   = i_IP1(n1,:,:)
    !                   
    j_IP1(n1+1,:,:)   = j_IP1(n1,:,:)
    j_IP1(n1+2,:,:)   = j_IP1(n1,:,:)
    j_IP1(n1+3,:,:)   = j_IP1(n1,:,:)
    j_IP1(n1+4,:,:)   = j_IP1(n1,:,:)
    j_IP1(n1+5,:,:)   = j_IP1(n1,:,:)
    j_IP1(n1+6,:,:)   = j_IP1(n1,:,:)
    !                   
    k_IP1(n1+1,:,:)   = k_IP1(n1,:,:)
    k_IP1(n1+2,:,:)   = k_IP1(n1,:,:)
    k_IP1(n1+3,:,:)   = k_IP1(n1,:,:)
    k_IP1(n1+4,:,:)   = k_IP1(n1,:,:)
    k_IP1(n1+5,:,:)   = k_IP1(n1,:,:)
    k_IP1(n1+6,:,:)   = k_IP1(n1,:,:)
    !
    i_IP2(n1+1,:,:)   = i_IP2(n1,:,:)
    i_IP2(n1+2,:,:)   = i_IP2(n1,:,:)
    i_IP2(n1+3,:,:)   = i_IP2(n1,:,:)
    i_IP2(n1+4,:,:)   = i_IP2(n1,:,:)
    i_IP2(n1+5,:,:)   = i_IP2(n1,:,:)
    i_IP2(n1+6,:,:)   = i_IP2(n1,:,:)
    !                   
    j_IP2(n1+1,:,:)   = j_IP2(n1,:,:)
    j_IP2(n1+2,:,:)   = j_IP2(n1,:,:)
    j_IP2(n1+3,:,:)   = j_IP2(n1,:,:)
    j_IP2(n1+4,:,:)   = j_IP2(n1,:,:)
    j_IP2(n1+5,:,:)   = j_IP2(n1,:,:)
    j_IP2(n1+6,:,:)   = j_IP2(n1,:,:)
    !                   
    k_IP2(n1+1,:,:)   = k_IP2(n1,:,:)
    k_IP2(n1+2,:,:)   = k_IP2(n1,:,:)
    k_IP2(n1+3,:,:)   = k_IP2(n1,:,:)
    k_IP2(n1+4,:,:)   = k_IP2(n1,:,:)
    k_IP2(n1+5,:,:)   = k_IP2(n1,:,:)
    k_IP2(n1+6,:,:)   = k_IP2(n1,:,:)
    !
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    i_mirror(:, 0,:)  = i_mirror(:,1,:)
    i_mirror(:,-1,:)  = i_mirror(:,1,:)
    i_mirror(:,-2,:)  = i_mirror(:,1,:)
    i_mirror(:,-3,:)  = i_mirror(:,1,:)
    i_mirror(:,-4,:)  = i_mirror(:,1,:)
    i_mirror(:,-5,:)  = i_mirror(:,1,:)
    !                   
    j_mirror(:, 0,:)  = j_mirror(:,1,:)
    j_mirror(:,-1,:)  = j_mirror(:,1,:)
    j_mirror(:,-2,:)  = j_mirror(:,1,:)
    j_mirror(:,-3,:)  = j_mirror(:,1,:)
    j_mirror(:,-4,:)  = j_mirror(:,1,:)
    j_mirror(:,-5,:)  = j_mirror(:,1,:)
    !                   
    k_mirror(:, 0,:)  = k_mirror(:,1,:)
    k_mirror(:,-1,:)  = k_mirror(:,1,:)
    k_mirror(:,-2,:)  = k_mirror(:,1,:)
    k_mirror(:,-3,:)  = k_mirror(:,1,:)
    k_mirror(:,-4,:)  = k_mirror(:,1,:)
    k_mirror(:,-5,:)  = k_mirror(:,1,:)
    !
    i_IP1(:, 0,:)     = i_IP1(:,1,:)
    i_IP1(:,-1,:)     = i_IP1(:,1,:)
    i_IP1(:,-2,:)     = i_IP1(:,1,:)
    i_IP1(:,-3,:)     = i_IP1(:,1,:)
    i_IP1(:,-4,:)     = i_IP1(:,1,:)
    i_IP1(:,-5,:)     = i_IP1(:,1,:)
    !                   
    j_IP1(:, 0,:)     = j_IP1(:,1,:)
    j_IP1(:,-1,:)     = j_IP1(:,1,:)
    j_IP1(:,-2,:)     = j_IP1(:,1,:)
    j_IP1(:,-3,:)     = j_IP1(:,1,:)
    j_IP1(:,-4,:)     = j_IP1(:,1,:)
    j_IP1(:,-5,:)     = j_IP1(:,1,:)
    !                   
    k_IP1(:, 0,:)     = k_IP1(:,1,:)
    k_IP1(:,-1,:)     = k_IP1(:,1,:)
    k_IP1(:,-2,:)     = k_IP1(:,1,:)
    k_IP1(:,-3,:)     = k_IP1(:,1,:)
    k_IP1(:,-4,:)     = k_IP1(:,1,:)
    k_IP1(:,-5,:)     = k_IP1(:,1,:)
    !
    i_IP2(:, 0,:)     = i_IP2(:,1,:)
    i_IP2(:,-1,:)     = i_IP2(:,1,:)
    i_IP2(:,-2,:)     = i_IP2(:,1,:)
    i_IP2(:,-3,:)     = i_IP2(:,1,:)
    i_IP2(:,-4,:)     = i_IP2(:,1,:)
    i_IP2(:,-5,:)     = i_IP2(:,1,:)
    !                   
    j_IP2(:, 0,:)     = j_IP2(:,1,:)
    j_IP2(:,-1,:)     = j_IP2(:,1,:)
    j_IP2(:,-2,:)     = j_IP2(:,1,:)
    j_IP2(:,-3,:)     = j_IP2(:,1,:)
    j_IP2(:,-4,:)     = j_IP2(:,1,:)
    j_IP2(:,-5,:)     = j_IP2(:,1,:)
    !                   
    k_IP2(:, 0,:)     = k_IP2(:,1,:)
    k_IP2(:,-1,:)     = k_IP2(:,1,:)
    k_IP2(:,-2,:)     = k_IP2(:,1,:)
    k_IP2(:,-3,:)     = k_IP2(:,1,:)
    k_IP2(:,-4,:)     = k_IP2(:,1,:)
    k_IP2(:,-5,:)     = k_IP2(:,1,:)
    !
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part               
    i_mirror(:,n1+1,:) = i_mirror(:,n1,:)
    i_mirror(:,n1+2,:) = i_mirror(:,n1,:)
    i_mirror(:,n1+3,:) = i_mirror(:,n1,:)
    i_mirror(:,n1+4,:) = i_mirror(:,n1,:)
    i_mirror(:,n1+5,:) = i_mirror(:,n1,:)
    i_mirror(:,n1+6,:) = i_mirror(:,n1,:)
    !                    
    j_mirror(:,n1+1,:) = j_mirror(:,n1,:)
    j_mirror(:,n1+2,:) = j_mirror(:,n1,:)
    j_mirror(:,n1+3,:) = j_mirror(:,n1,:)
    j_mirror(:,n1+4,:) = j_mirror(:,n1,:)
    j_mirror(:,n1+5,:) = j_mirror(:,n1,:)
    j_mirror(:,n1+6,:) = j_mirror(:,n1,:)
    !                    
    k_mirror(:,n1+1,:) = k_mirror(:,n1,:)
    k_mirror(:,n1+2,:) = k_mirror(:,n1,:)
    k_mirror(:,n1+3,:) = k_mirror(:,n1,:)
    k_mirror(:,n1+4,:) = k_mirror(:,n1,:)
    k_mirror(:,n1+5,:) = k_mirror(:,n1,:)
    k_mirror(:,n1+6,:) = k_mirror(:,n1,:)
    !               
    i_IP1(:,n1+1,:)   = i_IP1(:,n1,:)
    i_IP1(:,n1+2,:)   = i_IP1(:,n1,:)
    i_IP1(:,n1+3,:)   = i_IP1(:,n1,:)
    i_IP1(:,n1+4,:)   = i_IP1(:,n1,:)
    i_IP1(:,n1+5,:)   = i_IP1(:,n1,:)
    i_IP1(:,n1+6,:)   = i_IP1(:,n1,:)
    !                   
    j_IP1(:,n1+1,:)   = j_IP1(:,n1,:)
    j_IP1(:,n1+2,:)   = j_IP1(:,n1,:)
    j_IP1(:,n1+3,:)   = j_IP1(:,n1,:)
    j_IP1(:,n1+4,:)   = j_IP1(:,n1,:)
    j_IP1(:,n1+5,:)   = j_IP1(:,n1,:)
    j_IP1(:,n1+6,:)   = j_IP1(:,n1,:)
    !                   
    k_IP1(:,n1+1,:)   = k_IP1(:,n1,:)
    k_IP1(:,n1+2,:)   = k_IP1(:,n1,:)
    k_IP1(:,n1+3,:)   = k_IP1(:,n1,:)
    k_IP1(:,n1+4,:)   = k_IP1(:,n1,:)
    k_IP1(:,n1+5,:)   = k_IP1(:,n1,:)
    k_IP1(:,n1+6,:)   = k_IP1(:,n1,:)
    !               
    i_IP2(:,n1+1,:)   = i_IP2(:,n1,:)
    i_IP2(:,n1+2,:)   = i_IP2(:,n1,:)
    i_IP2(:,n1+3,:)   = i_IP2(:,n1,:)
    i_IP2(:,n1+4,:)   = i_IP2(:,n1,:)
    i_IP2(:,n1+5,:)   = i_IP2(:,n1,:)
    i_IP2(:,n1+6,:)   = i_IP2(:,n1,:)
    !                   
    j_IP2(:,n1+1,:)   = j_IP2(:,n1,:)
    j_IP2(:,n1+2,:)   = j_IP2(:,n1,:)
    j_IP2(:,n1+3,:)   = j_IP2(:,n1,:)
    j_IP2(:,n1+4,:)   = j_IP2(:,n1,:)
    j_IP2(:,n1+5,:)   = j_IP2(:,n1,:)
    j_IP2(:,n1+6,:)   = j_IP2(:,n1,:)
    !                   
    k_IP2(:,n1+1,:)   = k_IP2(:,n1,:)
    k_IP2(:,n1+2,:)   = k_IP2(:,n1,:)
    k_IP2(:,n1+3,:)   = k_IP2(:,n1,:)
    k_IP2(:,n1+4,:)   = k_IP2(:,n1,:)
    k_IP2(:,n1+5,:)   = k_IP2(:,n1,:)
    k_IP2(:,n1+6,:)   = k_IP2(:,n1,:)
    ! 
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary                         
    i_mirror(:,:, 0)  = i_mirror(:,:,1)
    i_mirror(:,:,-1)  = i_mirror(:,:,1)
    i_mirror(:,:,-2)  = i_mirror(:,:,1)
    i_mirror(:,:,-3)  = i_mirror(:,:,1)
    i_mirror(:,:,-4)  = i_mirror(:,:,1)
    i_mirror(:,:,-5)  = i_mirror(:,:,1)
    !                     
    j_mirror(:,:, 0)  = j_mirror(:,:,1)
    j_mirror(:,:,-1)  = j_mirror(:,:,1)
    j_mirror(:,:,-2)  = j_mirror(:,:,1)
    j_mirror(:,:,-3)  = j_mirror(:,:,1)
    j_mirror(:,:,-4)  = j_mirror(:,:,1)
    j_mirror(:,:,-5)  = j_mirror(:,:,1)
    !                   
    k_mirror(:,:, 0)  = k_mirror(:,:,1)
    k_mirror(:,:,-1)  = k_mirror(:,:,1)
    k_mirror(:,:,-2)  = k_mirror(:,:,1)
    k_mirror(:,:,-3)  = k_mirror(:,:,1)
    k_mirror(:,:,-4)  = k_mirror(:,:,1)
    k_mirror(:,:,-5)  = k_mirror(:,:,1)
    !                         
    i_IP1(:,:, 0)     = i_IP1(:,:,1)
    i_IP1(:,:,-1)     = i_IP1(:,:,1)
    i_IP1(:,:,-2)     = i_IP1(:,:,1)
    i_IP1(:,:,-3)     = i_IP1(:,:,1)
    i_IP1(:,:,-4)     = i_IP1(:,:,1)
    i_IP1(:,:,-5)     = i_IP1(:,:,1)
    !                        
    j_IP1(:,:, 0)     = j_IP1(:,:,1)
    j_IP1(:,:,-1)     = j_IP1(:,:,1)
    j_IP1(:,:,-2)     = j_IP1(:,:,1)
    j_IP1(:,:,-3)     = j_IP1(:,:,1)
    j_IP1(:,:,-4)     = j_IP1(:,:,1)
    j_IP1(:,:,-5)     = j_IP1(:,:,1)
    !                   
    k_IP1(:,:, 0)     = k_IP1(:,:,1)
    k_IP1(:,:,-1)     = k_IP1(:,:,1)
    k_IP1(:,:,-2)     = k_IP1(:,:,1)
    k_IP1(:,:,-3)     = k_IP1(:,:,1)
    k_IP1(:,:,-4)     = k_IP1(:,:,1)
    k_IP1(:,:,-5)     = k_IP1(:,:,1)
    !                         
    i_IP2(:,:, 0)     = i_IP2(:,:,1)
    i_IP2(:,:,-1)     = i_IP2(:,:,1)
    i_IP2(:,:,-2)     = i_IP2(:,:,1)
    i_IP2(:,:,-3)     = i_IP2(:,:,1)
    i_IP2(:,:,-4)     = i_IP2(:,:,1)
    i_IP2(:,:,-5)     = i_IP2(:,:,1)
    !                   
    j_IP2(:,:, 0)     = j_IP2(:,:,1)
    j_IP2(:,:,-1)     = j_IP2(:,:,1)
    j_IP2(:,:,-2)     = j_IP2(:,:,1)
    j_IP2(:,:,-3)     = j_IP2(:,:,1)
    j_IP2(:,:,-4)     = j_IP2(:,:,1)
    j_IP2(:,:,-5)     = j_IP2(:,:,1)
    !                   
    k_IP2(:,:, 0)     = k_IP2(:,:,1)
    k_IP2(:,:,-1)     = k_IP2(:,:,1)
    k_IP2(:,:,-2)     = k_IP2(:,:,1)
    k_IP2(:,:,-3)     = k_IP2(:,:,1)
    k_IP2(:,:,-4)     = k_IP2(:,:,1)
    k_IP2(:,:,-5)     = k_IP2(:,:,1)
    ! 
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary                        
    i_mirror(:,:,n3+1) = i_mirror(:,:,n3)
    i_mirror(:,:,n3+2) = i_mirror(:,:,n3)
    i_mirror(:,:,n3+3) = i_mirror(:,:,n3)
    i_mirror(:,:,n3+4) = i_mirror(:,:,n3)
    i_mirror(:,:,n3+5) = i_mirror(:,:,n3)
    i_mirror(:,:,n3+6) = i_mirror(:,:,n3)
    !                      
    j_mirror(:,:,n3+1) = j_mirror(:,:,n3)
    j_mirror(:,:,n3+2) = j_mirror(:,:,n3)
    j_mirror(:,:,n3+3) = j_mirror(:,:,n3)
    j_mirror(:,:,n3+4) = j_mirror(:,:,n3)
    j_mirror(:,:,n3+5) = j_mirror(:,:,n3)
    j_mirror(:,:,n3+6) = j_mirror(:,:,n3)
    !                    
    k_mirror(:,:,n3+1) = k_mirror(:,:,n3)
    k_mirror(:,:,n3+2) = k_mirror(:,:,n3)
    k_mirror(:,:,n3+3) = k_mirror(:,:,n3)
    k_mirror(:,:,n3+4) = k_mirror(:,:,n3)
    k_mirror(:,:,n3+5) = k_mirror(:,:,n3)
    k_mirror(:,:,n3+6) = k_mirror(:,:,n3)
    !                      
    i_IP1(:,:,n3+1)   = i_IP1(:,:,n3)
    i_IP1(:,:,n3+2)   = i_IP1(:,:,n3)
    i_IP1(:,:,n3+3)   = i_IP1(:,:,n3)
    i_IP1(:,:,n3+4)   = i_IP1(:,:,n3)
    i_IP1(:,:,n3+5)   = i_IP1(:,:,n3)
    i_IP1(:,:,n3+6)   = i_IP1(:,:,n3)
    !                        
    j_IP1(:,:,n3+1)   = j_IP1(:,:,n3)
    j_IP1(:,:,n3+2)   = j_IP1(:,:,n3)
    j_IP1(:,:,n3+3)   = j_IP1(:,:,n3)
    j_IP1(:,:,n3+4)   = j_IP1(:,:,n3)
    j_IP1(:,:,n3+5)   = j_IP1(:,:,n3)
    j_IP1(:,:,n3+6)   = j_IP1(:,:,n3)
    !                   
    k_IP1(:,:,n3+1)   = k_IP1(:,:,n3)
    k_IP1(:,:,n3+2)   = k_IP1(:,:,n3)
    k_IP1(:,:,n3+3)   = k_IP1(:,:,n3)
    k_IP1(:,:,n3+4)   = k_IP1(:,:,n3)
    k_IP1(:,:,n3+5)   = k_IP1(:,:,n3)
    k_IP1(:,:,n3+6)   = k_IP1(:,:,n3)
    !                      
    i_IP2(:,:,n3+1)   = i_IP2(:,:,n3)
    i_IP2(:,:,n3+2)   = i_IP2(:,:,n3)
    i_IP2(:,:,n3+3)   = i_IP2(:,:,n3)
    i_IP2(:,:,n3+4)   = i_IP2(:,:,n3)
    i_IP2(:,:,n3+5)   = i_IP2(:,:,n3)
    i_IP2(:,:,n3+6)   = i_IP2(:,:,n3)
    !                   
    j_IP2(:,:,n3+1)   = j_IP2(:,:,n3)
    j_IP2(:,:,n3+2)   = j_IP2(:,:,n3)
    j_IP2(:,:,n3+3)   = j_IP2(:,:,n3)
    j_IP2(:,:,n3+4)   = j_IP2(:,:,n3)
    j_IP2(:,:,n3+5)   = j_IP2(:,:,n3)
    j_IP2(:,:,n3+6)   = j_IP2(:,:,n3)
    !                   
    k_IP2(:,:,n3+1)   = k_IP2(:,:,n3)
    k_IP2(:,:,n3+2)   = k_IP2(:,:,n3)
    k_IP2(:,:,n3+3)   = k_IP2(:,:,n3)
    k_IP2(:,:,n3+4)   = k_IP2(:,:,n3)
    k_IP2(:,:,n3+5)   = k_IP2(:,:,n3)
    k_IP2(:,:,n3+6)   = k_IP2(:,:,n3)
    !
  endif
  !
  if (myid.eq.0)  print*, '**** Mirror points calculated ****'
  !------------------------------------------------------------
  !
  WP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7) = 0.0_rp
  WP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7) = 0.0_rp
  !
  call InterpolationWeights(n,ng,nh_d,nh_b,nabs_surf,Level_set, &
                            x_IP1,y_IP1,z_IP1, &
                            x_IP2,y_IP2,z_IP2, &
                            i_IP1,j_IP1,k_IP1, &
                            i_IP2,j_IP2,k_IP2, &
                            WP1,WP2, &
                            dl,dzc,zc)
  !
  !@cuf istat=cudaDeviceSynchronize()
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  !
  do q = 1,7
   call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,WP1(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,WP1(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,WP2(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,WP2(:,:,:,q))
  enddo
  !
  do q = 1,7
   if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    !
    WP1( 0,:,:,q)     = WP1(1,:,:,q)
    WP1(-1,:,:,q)     = WP1(1,:,:,q)
    WP1(-2,:,:,q)     = WP1(1,:,:,q)
    WP1(-3,:,:,q)     = WP1(1,:,:,q)
    WP1(-4,:,:,q)     = WP1(1,:,:,q)
    WP1(-5,:,:,q)     = WP1(1,:,:,q)
    !                      
    WP2( 0,:,:,q)     = WP2(1,:,:,q)
    WP2(-1,:,:,q)     = WP2(1,:,:,q)
    WP2(-2,:,:,q)     = WP2(1,:,:,q)
    WP2(-3,:,:,q)     = WP2(1,:,:,q)
    WP2(-4,:,:,q)     = WP2(1,:,:,q)
    WP2(-5,:,:,q)     = WP2(1,:,:,q)
    !
   endif
   if(right.eq.MPI_PROC_NULL) then ! x - upper part
    !
    WP1(n1+1,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+2,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+3,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+4,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+5,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+6,:,:,q)   = WP1(n1,:,:,q)
    !                      
    WP2(n1+1,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+2,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+3,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+4,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+5,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+6,:,:,q)   = WP2(n1,:,:,q)
    !
   endif
   !
   ! along y
   !
   if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    !
    WP1(:, 0,:,q)     = WP1(:,1,:,q)
    WP1(:,-1,:,q)     = WP1(:,1,:,q)
    WP1(:,-2,:,q)     = WP1(:,1,:,q)
    WP1(:,-3,:,q)     = WP1(:,1,:,q)
    WP1(:,-4,:,q)     = WP1(:,1,:,q)
    WP1(:,-5,:,q)     = WP1(:,1,:,q)
    !                      
    WP2(:, 0,:,q)     = WP2(:,1,:,q)
    WP2(:,-1,:,q)     = WP2(:,1,:,q)
    WP2(:,-2,:,q)     = WP2(:,1,:,q)
    WP2(:,-3,:,q)     = WP2(:,1,:,q)
    WP2(:,-4,:,q)     = WP2(:,1,:,q)
    WP2(:,-5,:,q)     = WP2(:,1,:,q)
    !
    endif
   if(back .eq.MPI_PROC_NULL) then ! y - upper part               
    !               
    WP1(:,n1+1,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+2,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+3,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+4,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+5,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+6,:,q)   = WP1(:,n1,:,q)
    !                      
    WP2(:,n1+1,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+2,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+3,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+4,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+5,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+6,:,q)   = WP2(:,n1,:,q)
    ! 
   endif
   !
   ! along z
   !
   if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    !  
    WP1(:,:, 0,q)     = WP1(:,:,1,q)
    WP1(:,:,-1,q)     = WP1(:,:,1,q)
    WP1(:,:,-2,q)     = WP1(:,:,1,q)
    WP1(:,:,-3,q)     = WP1(:,:,1,q)
    WP1(:,:,-4,q)     = WP1(:,:,1,q)
    WP1(:,:,-5,q)     = WP1(:,:,1,q)
    !                      
    WP2(:,:, 0,q)     = WP2(:,:,1,q)
    WP2(:,:,-1,q)     = WP2(:,:,1,q)
    WP2(:,:,-2,q)     = WP2(:,:,1,q)
    WP2(:,:,-3,q)     = WP2(:,:,1,q)
    WP2(:,:,-4,q)     = WP2(:,:,1,q)
    WP2(:,:,-5,q)     = WP2(:,:,1,q)
    !
   endif
   if(top.eq.MPI_PROC_NULL) then ! z - top boundary                        
    !                      
    WP1(:,:,n3+1,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+2,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+3,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+4,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+5,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+6,q)   = WP1(:,:,n3,q)
    !                        
    WP2(:,:,n3+1,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+2,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+3,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+4,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+5,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+6,q)   = WP2(:,:,n3,q)
    !
   endif
  enddo
  !
  if (myid.eq.0)  print*, '*** Interpolation Weights calculated ***'
  !
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf , size(nx_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf , size(ny_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf , size(nz_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf , size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_mirror , size(i_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_mirror , size(j_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_mirror , size(k_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1 , size(i_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1 , size(j_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1 , size(k_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2 , size(i_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2 , size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2 , size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1 , size(WP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2 , size(WP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan , size(deltan), cudaCpuDeviceId, 0)
  !
  call load('w',trim(restart_dir)//'utag.bin',n,         cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'vtag.bin',n,         cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'wtag.bin',n,         cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'stag.bin',n,       cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'lvl.bin',n,       Level_set(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nx.bin',n,              nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ny.bin',n,              ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nz.bin',n,              nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ns.bin',n,            nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'deltan.bin',n,           deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'imirror.bin',n,    i_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jmirror.bin',n,    j_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kmirror.bin',n,    k_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP1.bin',n,          i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP1.bin',n,          j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP1.bin',n,          k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP2.bin',n,          i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP2.bin',n,          j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP2.bin',n,          k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('w',trim(restart_dir)//'WP1.bin',n,          WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('w',trim(restart_dir)//'WP2.bin',n,          WP2(1:n(1),1:n(2),1:n(3),1:7))
  !
  if(allocated(x_mirror))    deallocate(x_mirror)
  if(allocated(y_mirror))    deallocate(y_mirror)
  if(allocated(z_mirror))    deallocate(z_mirror)
  if(allocated(x_IP1))       deallocate(x_IP1)
  if(allocated(y_IP1))       deallocate(y_IP1)
  if(allocated(z_IP1))       deallocate(z_IP1)
  if(allocated(x_IP2))       deallocate(x_IP2)
  if(allocated(y_IP2))       deallocate(y_IP2)
  if(allocated(z_IP2))       deallocate(z_IP2)
  if(allocated(x_intersect)) deallocate(x_intersect)
  if(allocated(y_intersect)) deallocate(y_intersect)
  if(allocated(z_intersect)) deallocate(z_intersect)
  !
else
  !
  call load('r',trim(restart_dir)//'utag.bin',n,        cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'vtag.bin',n,        cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'wtag.bin',n,        cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'stag.bin',n,      cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'lvl.bin',n,      Level_set(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nx.bin',n,             nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ny.bin',n,             ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nz.bin',n,             nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ns.bin',n,           nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'deltan.bin',n,          deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'imirror.bin',n,   i_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jmirror.bin',n,   j_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kmirror.bin',n,   k_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP1.bin',n,         i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP1.bin',n,         j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP1.bin',n,         k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP2.bin',n,         i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP2.bin',n,         j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP2.bin',n,         k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('r',trim(restart_dir)//'WP1.bin',n,         WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('r',trim(restart_dir)//'WP2.bin',n,         WP2(1:n(1),1:n(2),1:n(3),1:7))
  !
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag),     mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set),       mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf),           mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf),       mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan),             mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_mirror, size(i_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_mirror, size(j_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_mirror, size(k_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1, size(WP1),                   mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2, size(WP2),                   mydev, 0)
  !
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(Level_set(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  Level_set(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = &
  int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_phi_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_phi_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_u_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_u_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_v_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_v_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,cell_w_tag)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,cell_w_tag)
  !
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nabs_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nabs_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nx_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nx_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,ny_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,ny_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,nz_surf)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,nz_surf)
  !
  !
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_mirror(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_IP1(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(i_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  i_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(j_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  j_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !$acc wait
  !$acc kernels
  temp(1:n(1),1:n(2),1:n(3)) = real(k_IP2(1:n(1),1:n(2),1:n(3)),rp)
  !$acc end kernels
  call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,temp)
  call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,temp)
  !$acc kernels
  k_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = int(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  !$acc end kernels
  !
  do q = 1,7
   call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,WP1(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,WP1(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind1),ind1,WP2(:,:,:,q))
   call updthalo(n1,n2,n3,nh_b,halo(ind2),ind2,WP2(:,:,:,q))
  enddo
  !
  ! along x
  !
  if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    cell_u_tag( 0,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-1,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-2,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-3,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-4,:,:)     = cell_u_tag(1,:,:)
    cell_u_tag(-5,:,:)     = cell_u_tag(1,:,:)
    !
    cell_v_tag( 0,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-1,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-2,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-3,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-4,:,:)     = cell_v_tag(1,:,:)
    cell_v_tag(-5,:,:)     = cell_v_tag(1,:,:)
    !
    cell_w_tag( 0,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-1,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-2,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-3,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-4,:,:)     = cell_w_tag(1,:,:)
    cell_w_tag(-5,:,:)     = cell_w_tag(1,:,:)
    !
    cell_phi_tag( 0,:,:)   = cell_phi_tag(1,:,:)
    cell_phi_tag(-1,:,:)   = cell_phi_tag(1,:,:)
    cell_phi_tag(-2,:,:)   = cell_phi_tag(1,:,:)
    cell_phi_tag(-3,:,:)   = cell_phi_tag(1,:,:)
    cell_phi_tag(-4,:,:)   = cell_phi_tag(1,:,:)
    cell_phi_tag(-5,:,:)   = cell_phi_tag(1,:,:)
    !
    Level_set( 0,:,:)     = Level_set(1,:,:)
    Level_set(-1,:,:)     = Level_set(1,:,:)
    Level_set(-2,:,:)     = Level_set(1,:,:)
    Level_set(-3,:,:)     = Level_set(1,:,:)
    Level_set(-4,:,:)     = Level_set(1,:,:)
    Level_set(-5,:,:)     = Level_set(1,:,:)
    !
    nabs_surf( 0,:,:)     = nabs_surf(1,:,:)
    nabs_surf(-1,:,:)     = nabs_surf(1,:,:)
    nabs_surf(-2,:,:)     = nabs_surf(1,:,:)
    nabs_surf(-3,:,:)     = nabs_surf(1,:,:)
    nabs_surf(-4,:,:)     = nabs_surf(1,:,:)
    nabs_surf(-5,:,:)     = nabs_surf(1,:,:)
    !
    nx_surf( 0,:,:)     = nx_surf(1,:,:)
    nx_surf(-1,:,:)     = nx_surf(1,:,:)
    nx_surf(-2,:,:)     = nx_surf(1,:,:)
    nx_surf(-3,:,:)     = nx_surf(1,:,:)
    nx_surf(-4,:,:)     = nx_surf(1,:,:)
    nx_surf(-5,:,:)     = nx_surf(1,:,:)
    !
    ny_surf( 0,:,:)     = ny_surf(1,:,:)
    ny_surf(-1,:,:)     = ny_surf(1,:,:)
    ny_surf(-2,:,:)     = ny_surf(1,:,:)
    ny_surf(-3,:,:)     = ny_surf(1,:,:)
    ny_surf(-4,:,:)     = ny_surf(1,:,:)
    ny_surf(-5,:,:)     = ny_surf(1,:,:)
    !
    nz_surf( 0,:,:)     = nz_surf(1,:,:)
    nz_surf(-1,:,:)     = nz_surf(1,:,:)
    nz_surf(-2,:,:)     = nz_surf(1,:,:)
    nz_surf(-3,:,:)     = nz_surf(1,:,:)
    nz_surf(-4,:,:)     = nz_surf(1,:,:)
    nz_surf(-5,:,:)     = nz_surf(1,:,:)
  endif
  if(right.eq.MPI_PROC_NULL) then ! x - upper part
    cell_u_tag(n1+1,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+2,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+3,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+4,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+5,:,:)   = cell_u_tag(n1,:,:)
    cell_u_tag(n1+6,:,:)   = cell_u_tag(n1,:,:)
    !
    cell_v_tag(n1+1,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+2,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+3,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+4,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+5,:,:)   = cell_v_tag(n1,:,:)
    cell_v_tag(n1+6,:,:)   = cell_v_tag(n1,:,:)
    !
    cell_w_tag(n1+1,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+2,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+3,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+4,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+5,:,:)   = cell_w_tag(n1,:,:)
    cell_w_tag(n1+6,:,:)   = cell_w_tag(n1,:,:)
    !
    cell_phi_tag(n1+1,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+2,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+3,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+4,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+5,:,:) = cell_phi_tag(n1,:,:)
    cell_phi_tag(n1+6,:,:) = cell_phi_tag(n1,:,:)
    !
    Level_set(n1+1,:,:)   = Level_set(n1,:,:)
    Level_set(n1+2,:,:)   = Level_set(n1,:,:)
    Level_set(n1+3,:,:)   = Level_set(n1,:,:)
    Level_set(n1+4,:,:)   = Level_set(n1,:,:)
    Level_set(n1+5,:,:)   = Level_set(n1,:,:)
    Level_set(n1+6,:,:)   = Level_set(n1,:,:)
    !
    nabs_surf(n1+1,:,:)   = nabs_surf(n1,:,:)
    nabs_surf(n1+2,:,:)   = nabs_surf(n1,:,:)
    nabs_surf(n1+3,:,:)   = nabs_surf(n1,:,:)
    nabs_surf(n1+4,:,:)   = nabs_surf(n1,:,:)
    nabs_surf(n1+5,:,:)   = nabs_surf(n1,:,:)
    nabs_surf(n1+6,:,:)   = nabs_surf(n1,:,:)
    !
    nx_surf(n1+1,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+2,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+3,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+4,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+5,:,:)   = nx_surf(n1,:,:)
    nx_surf(n1+6,:,:)   = nx_surf(n1,:,:)
    !
    ny_surf(n1+1,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+2,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+3,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+4,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+5,:,:)   = ny_surf(n1,:,:)
    ny_surf(n1+6,:,:)   = ny_surf(n1,:,:)
    !
    nz_surf(n1+1,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+2,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+3,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+4,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+5,:,:)   = nz_surf(n1,:,:)
    nz_surf(n1+6,:,:)   = nz_surf(n1,:,:)
  endif
  !
  ! along y
  !
  if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    cell_u_tag(:, 0,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-1,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-2,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-3,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-4,:)     = cell_u_tag(:,1,:)
    cell_u_tag(:,-5,:)     = cell_u_tag(:,1,:)
    !
    cell_v_tag(:, 0,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-1,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-2,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-3,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-4,:)     = cell_v_tag(:,1,:)
    cell_v_tag(:,-5,:)     = cell_v_tag(:,1,:)
    !
    cell_w_tag(:, 0,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-1,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-2,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-3,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-4,:)     = cell_w_tag(:,1,:)
    cell_w_tag(:,-5,:)     = cell_w_tag(:,1,:)
    !
    cell_phi_tag(:, 0,:)     = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-1,:)     = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-2,:)     = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-3,:)     = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-4,:)     = cell_phi_tag(:,1,:)
    cell_phi_tag(:,-5,:)     = cell_phi_tag(:,1,:)
    !
    Level_set(:, 0,:)     = Level_set(:,1,:)
    Level_set(:,-1,:)     = Level_set(:,1,:)
    Level_set(:,-2,:)     = Level_set(:,1,:)
    Level_set(:,-3,:)     = Level_set(:,1,:)
    Level_set(:,-4,:)     = Level_set(:,1,:)
    Level_set(:,-5,:)     = Level_set(:,1,:)
    !
    nabs_surf(:, 0,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-1,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-2,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-3,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-4,:)   = nabs_surf(:,1,:)
    nabs_surf(:,-5,:)   = nabs_surf(:,1,:)
    !
    nx_surf(:, 0,:)     = nx_surf(:,1,:)
    nx_surf(:,-1,:)     = nx_surf(:,1,:)
    nx_surf(:,-2,:)     = nx_surf(:,1,:)
    nx_surf(:,-3,:)     = nx_surf(:,1,:)
    nx_surf(:,-4,:)     = nx_surf(:,1,:)
    nx_surf(:,-5,:)     = nx_surf(:,1,:)
    !
    ny_surf(:, 0,:)     = ny_surf(:,1,:)
    ny_surf(:,-1,:)     = ny_surf(:,1,:)
    ny_surf(:,-2,:)     = ny_surf(:,1,:)
    ny_surf(:,-3,:)     = ny_surf(:,1,:)
    ny_surf(:,-4,:)     = ny_surf(:,1,:)
    ny_surf(:,-5,:)     = ny_surf(:,1,:)
    !
    nz_surf(:, 0,:)     = nz_surf(:,1,:)
    nz_surf(:,-1,:)     = nz_surf(:,1,:)
    nz_surf(:,-2,:)     = nz_surf(:,1,:)
    nz_surf(:,-3,:)     = nz_surf(:,1,:)
    nz_surf(:,-4,:)     = nz_surf(:,1,:)
    nz_surf(:,-5,:)     = nz_surf(:,1,:)
   endif
  if(back .eq.MPI_PROC_NULL) then ! y - upper part
    cell_u_tag(:,n1+1,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+2,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+3,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+4,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+5,:)   = cell_u_tag(:,n1,:)
    cell_u_tag(:,n1+6,:)   = cell_u_tag(:,n1,:)
    !
    cell_v_tag(:,n1+1,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+2,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+3,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+4,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+5,:)   = cell_v_tag(:,n1,:)
    cell_v_tag(:,n1+6,:)   = cell_v_tag(:,n1,:)
    !
    cell_w_tag(:,n1+1,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+2,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+3,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+4,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+5,:)   = cell_w_tag(:,n1,:)
    cell_w_tag(:,n1+6,:)   = cell_w_tag(:,n1,:)
    !
    cell_phi_tag(:,n1+1,:)   = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+2,:)   = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+3,:)   = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+4,:)   = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+5,:)   = cell_phi_tag(:,n1,:)
    cell_phi_tag(:,n1+6,:)   = cell_phi_tag(:,n1,:)
    !
    Level_set(:,n1+1,:)   = Level_set(:,n1,:)
    Level_set(:,n1+2,:)   = Level_set(:,n1,:)
    Level_set(:,n1+3,:)   = Level_set(:,n1,:)
    Level_set(:,n1+4,:)   = Level_set(:,n1,:)
    Level_set(:,n1+5,:)   = Level_set(:,n1,:)
    Level_set(:,n1+6,:)   = Level_set(:,n1,:)
    !
    nabs_surf(:,n1+1,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+2,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+3,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+4,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+5,:) = nabs_surf(:,n1,:)
    nabs_surf(:,n1+6,:) = nabs_surf(:,n1,:)
    !
    nx_surf(:,n1+1,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+2,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+3,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+4,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+5,:)   = nx_surf(:,n1,:)
    nx_surf(:,n1+6,:)   = nx_surf(:,n1,:)
    !
    ny_surf(:,n1+1,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+2,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+3,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+4,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+5,:)   = ny_surf(:,n1,:)
    ny_surf(:,n1+6,:)   = ny_surf(:,n1,:)
    !
    nz_surf(:,n1+1,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+2,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+3,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+4,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+5,:)   = nz_surf(:,n1,:)
    nz_surf(:,n1+6,:)   = nz_surf(:,n1,:)
  endif
  !
  ! along z
  !
  if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    cell_u_tag(:,:, 0)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-1)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-2)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-3)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-4)     = cell_u_tag(:,:,1)
    cell_u_tag(:,:,-5)     = cell_u_tag(:,:,1)
    !
    cell_v_tag(:,:, 0)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-1)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-2)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-3)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-4)     = cell_v_tag(:,:,1)
    cell_v_tag(:,:,-5)     = cell_v_tag(:,:,1)
    !
    cell_w_tag(:,:, 0)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-1)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-2)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-3)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-4)     = cell_w_tag(:,:,1)
    cell_w_tag(:,:,-5)     = cell_w_tag(:,:,1)
    !
    cell_phi_tag(:,:, 0)     = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-1)     = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-2)     = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-3)     = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-4)     = cell_phi_tag(:,:,1)
    cell_phi_tag(:,:,-5)     = cell_phi_tag(:,:,1)
    !
    Level_set(:,:, 0)     = Level_set(:,:,1)
    Level_set(:,:,-1)     = Level_set(:,:,1)
    Level_set(:,:,-2)     = Level_set(:,:,1)
    Level_set(:,:,-3)     = Level_set(:,:,1)
    Level_set(:,:,-4)     = Level_set(:,:,1)
    Level_set(:,:,-5)     = Level_set(:,:,1)
    !
    nabs_surf(:,:, 0)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-1)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-2)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-3)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-4)   = nabs_surf(:,:,1)
    nabs_surf(:,:,-5)   = nabs_surf(:,:,1)
    !
    nx_surf(:,:, 0)     = nx_surf(:,:,1)
    nx_surf(:,:,-1)     = nx_surf(:,:,1)
    nx_surf(:,:,-2)     = nx_surf(:,:,1)
    nx_surf(:,:,-3)     = nx_surf(:,:,1)
    nx_surf(:,:,-4)     = nx_surf(:,:,1)
    nx_surf(:,:,-5)     = nx_surf(:,:,1)
    !
    ny_surf(:,:, 0)     = ny_surf(:,:,1)
    ny_surf(:,:,-1)     = ny_surf(:,:,1)
    ny_surf(:,:,-2)     = ny_surf(:,:,1)
    ny_surf(:,:,-3)     = ny_surf(:,:,1)
    ny_surf(:,:,-4)     = ny_surf(:,:,1)
    ny_surf(:,:,-5)     = ny_surf(:,:,1)
    !
    nz_surf(:,:, 0)     = nz_surf(:,:,1)
    nz_surf(:,:,-1)     = nz_surf(:,:,1)
    nz_surf(:,:,-2)     = nz_surf(:,:,1)
    nz_surf(:,:,-3)     = nz_surf(:,:,1)
    nz_surf(:,:,-4)     = nz_surf(:,:,1)
    nz_surf(:,:,-5)     = nz_surf(:,:,1)
  endif
  if(top.eq.MPI_PROC_NULL) then ! z - top boundary
    cell_u_tag(:,:,n3+1)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+2)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+3)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+4)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+5)   = cell_u_tag(:,:,n3)
    cell_u_tag(:,:,n3+6)   = cell_u_tag(:,:,n3)
    !
    cell_v_tag(:,:,n3+1)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+2)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+3)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+4)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+5)   = cell_v_tag(:,:,n3)
    cell_v_tag(:,:,n3+6)   = cell_v_tag(:,:,n3)
    !
    cell_w_tag(:,:,n3+1)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+2)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+3)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+4)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+5)   = cell_w_tag(:,:,n3)
    cell_w_tag(:,:,n3+6)   = cell_w_tag(:,:,n3)
    !
    cell_phi_tag(:,:,n3+1)   = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+2)   = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+3)   = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+4)   = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+5)   = cell_phi_tag(:,:,n3)
    cell_phi_tag(:,:,n3+6)   = cell_phi_tag(:,:,n3)
    !
    Level_set(:,:,n3+1)   = Level_set(:,:,n3)
    Level_set(:,:,n3+2)   = Level_set(:,:,n3)
    Level_set(:,:,n3+3)   = Level_set(:,:,n3)
    Level_set(:,:,n3+4)   = Level_set(:,:,n3)
    Level_set(:,:,n3+5)   = Level_set(:,:,n3)
    Level_set(:,:,n3+6)   = Level_set(:,:,n3)
    !
    nabs_surf(:,:,n3+1) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+2) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+3) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+4) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+5) = nabs_surf(:,:,n3)
    nabs_surf(:,:,n3+6) = nabs_surf(:,:,n3)
    !
    nx_surf(:,:,n3+1)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+2)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+3)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+4)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+5)   = nx_surf(:,:,n3)
    nx_surf(:,:,n3+6)   = nx_surf(:,:,n3)
    !
    ny_surf(:,:,n3+1)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+2)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+3)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+4)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+5)   = ny_surf(:,:,n3)
    ny_surf(:,:,n3+6)   = ny_surf(:,:,n3)
    !
    nz_surf(:,:,n3+1)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+2)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+3)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+4)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+5)   = nz_surf(:,:,n3)
    nz_surf(:,:,n3+6)   = nz_surf(:,:,n3)
  endif
  do q = 1,7
   if(left .eq.MPI_PROC_NULL) then ! x - bottom part
    !
    WP1( 0,:,:,q)     = WP1(1,:,:,q)
    WP1(-1,:,:,q)     = WP1(1,:,:,q)
    WP1(-2,:,:,q)     = WP1(1,:,:,q)
    WP1(-3,:,:,q)     = WP1(1,:,:,q)
    WP1(-4,:,:,q)     = WP1(1,:,:,q)
    WP1(-5,:,:,q)     = WP1(1,:,:,q)
    !                      
    WP2( 0,:,:,q)     = WP2(1,:,:,q)
    WP2(-1,:,:,q)     = WP2(1,:,:,q)
    WP2(-2,:,:,q)     = WP2(1,:,:,q)
    WP2(-3,:,:,q)     = WP2(1,:,:,q)
    WP2(-4,:,:,q)     = WP2(1,:,:,q)
    WP2(-5,:,:,q)     = WP2(1,:,:,q)
    !
   endif
   if(right.eq.MPI_PROC_NULL) then ! x - upper part
    !
    WP1(n1+1,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+2,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+3,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+4,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+5,:,:,q)   = WP1(n1,:,:,q)
    WP1(n1+6,:,:,q)   = WP1(n1,:,:,q)
    !                      
    WP2(n1+1,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+2,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+3,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+4,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+5,:,:,q)   = WP2(n1,:,:,q)
    WP2(n1+6,:,:,q)   = WP2(n1,:,:,q)
    !
   endif
   !
   ! along y
   !
   if(front.eq.MPI_PROC_NULL) then ! y - bottom part
    !
    WP1(:, 0,:,q)     = WP1(:,1,:,q)
    WP1(:,-1,:,q)     = WP1(:,1,:,q)
    WP1(:,-2,:,q)     = WP1(:,1,:,q)
    WP1(:,-3,:,q)     = WP1(:,1,:,q)
    WP1(:,-4,:,q)     = WP1(:,1,:,q)
    WP1(:,-5,:,q)     = WP1(:,1,:,q)
    !                      
    WP2(:, 0,:,q)     = WP2(:,1,:,q)
    WP2(:,-1,:,q)     = WP2(:,1,:,q)
    WP2(:,-2,:,q)     = WP2(:,1,:,q)
    WP2(:,-3,:,q)     = WP2(:,1,:,q)
    WP2(:,-4,:,q)     = WP2(:,1,:,q)
    WP2(:,-5,:,q)     = WP2(:,1,:,q)
    !
    endif
   if(back .eq.MPI_PROC_NULL) then ! y - upper part               
    !               
    WP1(:,n1+1,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+2,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+3,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+4,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+5,:,q)   = WP1(:,n1,:,q)
    WP1(:,n1+6,:,q)   = WP1(:,n1,:,q)
    !                      
    WP2(:,n1+1,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+2,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+3,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+4,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+5,:,q)   = WP2(:,n1,:,q)
    WP2(:,n1+6,:,q)   = WP2(:,n1,:,q)
    ! 
   endif
   !
   ! along z
   !
   if(bottom.eq.MPI_PROC_NULL) then ! z - bottom boundary
    !  
    WP1(:,:, 0,q)     = WP1(:,:,1,q)
    WP1(:,:,-1,q)     = WP1(:,:,1,q)
    WP1(:,:,-2,q)     = WP1(:,:,1,q)
    WP1(:,:,-3,q)     = WP1(:,:,1,q)
    WP1(:,:,-4,q)     = WP1(:,:,1,q)
    WP1(:,:,-5,q)     = WP1(:,:,1,q)
    !                      
    WP2(:,:, 0,q)     = WP2(:,:,1,q)
    WP2(:,:,-1,q)     = WP2(:,:,1,q)
    WP2(:,:,-2,q)     = WP2(:,:,1,q)
    WP2(:,:,-3,q)     = WP2(:,:,1,q)
    WP2(:,:,-4,q)     = WP2(:,:,1,q)
    WP2(:,:,-5,q)     = WP2(:,:,1,q)
    !
   endif
   if(top.eq.MPI_PROC_NULL) then ! z - top boundary                        
    !                      
    WP1(:,:,n3+1,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+2,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+3,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+4,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+5,q)   = WP1(:,:,n3,q)
    WP1(:,:,n3+6,q)   = WP1(:,:,n3,q)
    !                        
    WP2(:,:,n3+1,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+2,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+3,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+4,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+5,q)   = WP2(:,:,n3,q)
    WP2(:,:,n3+6,q)   = WP2(:,:,n3,q)
    !
   endif
  enddo
  !
  !
  !*********************************************************************

  if (myid.eq.0)  print*, '****** Existing IBM data loaded ******'

endif
!
!@cuf istat=cudaDeviceSynchronize()
call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!
!@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(nx_surf , size(nx_surf), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(ny_surf , size(ny_surf), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(nz_surf , size(nz_surf), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(nabs_surf , size(nabs_surf), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(i_mirror , size(i_mirror), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(j_mirror , size(j_mirror), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(k_mirror , size(k_mirror), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(i_IP1 , size(i_IP1), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(j_IP1 , size(j_IP1), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(k_IP1 , size(k_IP1), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(i_IP2 , size(i_IP2), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(j_IP2 , size(j_IP2), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(k_IP2 , size(k_IP2), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(WP1 , size(WP1), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(WP2 , size(WP2), cudaCpuDeviceId, 0)
!@cuf istat = cudaMemPrefetchAsync(deltan , size(deltan), cudaCpuDeviceId, 0)
!
istep = 0
write(fldnum,'(i9.9)') istep
call write_visu_3d(trim(data_dir),'sol_fld_'//fldnum//'.bin','log_visu_3d.out','sol', &
                   (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),0.0_rp,0, &
                   (1.0_rp - cell_phi_tag(1:n(1),1:n(2),1:n(3))))
call write_visu_3d(trim(data_dir),'nab_fld_'//fldnum//'.bin','log_visu_3d.out','nab', &
                   (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),0.0_rp,0, &
                   nabs_surf(1:n(1),1:n(2),1:n(3)))
call write_visu_3d(trim(data_dir),'n_x_fld_'//fldnum//'.bin','log_visu_3d.out','nx', &
                   (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),0.0_rp,0, &
                   nx_surf(1:n(1),1:n(2),1:n(3)))
call write_visu_3d(trim(data_dir),'n_y_fld_'//fldnum//'.bin','log_visu_3d.out','ny', &
                   (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),0.0_rp,0, &
                   ny_surf(1:n(1),1:n(2),1:n(3)))
call write_visu_3d(trim(data_dir),'n_z_fld_'//fldnum//'.bin','log_visu_3d.out','nz', &
                   (/1,1,1/),(/ng(1),ng(2),ng(3)/),(/1,1,1/),0.0_rp,0, &
                   nz_surf(1:n(1),1:n(2),1:n(3)))
end subroutine initIBM
#endif
end module mod_initIBM