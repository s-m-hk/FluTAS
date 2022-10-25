module mod_initIBM
#if defined(_USE_IBM)
use mpi
use mod_IBM
use mod_bound
use mod_load
use mod_param
use mod_common_mpi, only: myid,ierr
use mod_types
!@cuf use cudafor
!@cuf use mod_common_mpi, only: mydev
implicit none
private
public initIBM
contains
subroutine initIBM(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag,Level_set, &
                   nx_surf,ny_surf,nz_surf,nabs_surf, &
                   i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2,deltan, &
                   nh_d,nh_b,halo,&
                   zc,zf,dzc,dzf,dl,dli)
implicit none
!
character(len=100) :: datadir,restart_dir
integer, intent(in)                           :: nh_d,nh_b
integer, intent(in), dimension(3)             :: halo
real(rp), intent(in ), dimension(1-nh_d:)     :: zc,zf,dzc,dzf
real(rp), dimension(3), intent(in )           :: dl,dli
!
real(rp), intent(out),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
integer,  intent(out),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) :: Level_set,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), intent(out),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
real(rp), intent(out),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7) :: WP1,WP2
real(rp), allocatable :: temp(:,:,:)
real(rp), allocatable :: z_intersect(:,:,:),y_intersect(:,:,:), x_intersect(:,:,:)
integer,  allocatable :: i_mirror(:,:,:), j_mirror(:,:,:), k_mirror(:,:,:)
real(rp), allocatable :: x_mirror(:,:,:), y_mirror(:,:,:), z_mirror(:,:,:)
real(rp), allocatable :: x_IP1(:,:,:), y_IP1(:,:,:), z_IP1(:,:,:) !can be deallocated later
real(rp), allocatable :: x_IP2(:,:,:), y_IP2(:,:,:), z_IP2(:,:,:) !can be deallocated later 
logical:: is_data
integer::i,j,k,i1,j1,k1
!@cuf integer :: istat
!@cuf attributes(managed) :: dzc,dzf,zc,zf
!@cuf attributes(managed) :: temp,Level_set,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
!@cuf attributes(managed) :: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
!@cuf attributes(managed) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
!@cuf attributes(managed) :: i_mirror, j_mirror, k_mirror
!@cuf attributes(managed) :: x_mirror, y_mirror, z_mirror, x_IP1, y_IP1, z_IP1, x_IP2, y_IP2, z_IP2
!@cuf attributes(managed) :: z_intersect, y_intersect, x_intersect
!@cuf attributes(managed) :: WP1,WP2
!
!@cuf istat = cudaMemAdvise(z_intersect, size(z_intersect), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(y_intersect, size(y_intersect), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(x_intersect, size(x_intersect), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(i_mirror, size(i_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(j_mirror, size(j_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(k_mirror, size(k_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(x_mirror, size(x_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(y_mirror, size(y_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(z_mirror, size(z_mirror), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(deltan, size(deltan), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(x_IP1, size(x_IP1), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(y_IP1, size(y_IP1), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(z_IP1, size(z_IP1), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(x_IP2, size(x_IP2), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(y_IP2, size(y_IP2), cudaMemAdviseSetReadMostly, mydev)
!@cuf istat = cudaMemAdvise(z_IP2, size(z_IP2), cudaMemAdviseSetReadMostly, mydev)
!
restart_dir = 'data/restart_dir/'
!
i1=n(1)+1
j1=n(2)+1
k1=n(3)+1
!
inquire(file=trim(restart_dir)//'stag.bin',exist=is_data)

  allocate(z_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           x_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(i_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           j_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           k_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(temp(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))

  cell_u_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  cell_v_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  cell_w_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  cell_phi_tag(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)  = 0.0_rp
  Level_set(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)     = 0

if (.not.is_data) then
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), cudaCpuDeviceId, 0)
#if !defined(_OPENACC)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), cudaCpuDeviceId, 0)
#endif
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set), cudaCpuDeviceId, 0)
  call IBM_mask(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag,Level_set,zc,zf,dl,dzc)
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), mydev, 0)
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_phi_tag)
   !$acc kernels 
   do k=0,n(3)
    do j=0,n(2)
     do i=0,n(1)
      cell_u_tag(i,j,k) = 0.5*(cell_phi_tag(i+1,j,k)+cell_phi_tag(i,j,k))
      cell_v_tag(i,j,k) = 0.5*(cell_phi_tag(i,j+1,k)+cell_phi_tag(i,j,k))
      cell_w_tag(i,j,k) = 0.5*(cell_phi_tag(i,j,k+1)+cell_phi_tag(i,j,k))
     enddo
    enddo
   enddo
   !$acc end kernels
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_b,halo,is_outflow,dl,dzc,dzf,cell_u_tag,cell_v_tag,cell_w_tag)

  if (myid.eq.0)  print*, '*** Volume fractions calculated ***'
  !---------------------------------------------------------------------
      !@cuf istat = cudaMemPrefetchAsync(temp,        size(temp),        mydev, 0)
      temp(:,:,:) = real(Level_set(:,:,:),rp)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,temp)
      Level_set(:,:,:) = int(temp(:,:,:))
  
  if (myid.eq.0)  print*, '*** Marker function calculated ***'
  !---------------------------------------------------------------------

  nx_surf(  1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  ny_surf(  1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  nz_surf(  1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  nabs_surf(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    = 0.0_rp
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  call normal_vectors(Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,dl,dli,zc,dzc)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), mydev, 0)
  call bounduvw(cbcvel,n,bcvel,nh_d,nh_b,halo,is_outflow,dl,dzc,dzf,nx_surf,ny_surf,nz_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nabs_surf)

  if (myid.eq.0)  print*, '*** Normal vectors calculated ***'
  !*********************************************************************
  x_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  y_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  z_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b) = 0.0_rp
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_intersect, size(x_intersect), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_intersect, size(y_intersect), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_intersect, size(z_intersect), cudaCpuDeviceId, 0)
  call intersect(nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dl,dzc,zc)
  !@cuf istat = cudaMemPrefetchAsync(x_intersect, size(x_intersect), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_intersect, size(y_intersect), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_intersect, size(z_intersect), mydev, 0)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_intersect)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_intersect)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_intersect)
  
  if (myid.eq.0)  print*, '*** Intersect points have been calculated! ***'
  !*********************************************************************
  x_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  y_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  z_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  x_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  x_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  y_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  y_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  z_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  z_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  deltan(  1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000.0_rp
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_intersect, size(x_intersect), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_intersect, size(y_intersect), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_intersect, size(z_intersect), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_mirror, size(x_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_mirror, size(y_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_mirror, size(z_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP1, size(x_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP1, size(y_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP1, size(z_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP2, size(x_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP2, size(y_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP2, size(z_IP2), cudaCpuDeviceId, 0)
      call mirrorpoints(nx_surf,ny_surf,nz_surf,nabs_surf, &
                        x_intersect,y_intersect,z_intersect, &
                        x_mirror,y_mirror,z_mirror, &
                        deltan, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        dl,dzc,zc)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_mirror, size(x_mirror), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_mirror, size(y_mirror), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_mirror, size(z_mirror), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP1, size(x_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP1, size(y_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP1, size(z_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP2, size(x_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP2, size(y_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP2, size(z_IP2), mydev, 0)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,deltan)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_mirror)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_mirror)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_mirror)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_IP1)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_IP2)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_IP1)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_IP2)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_IP1)
      call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_IP2)

      i_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      j_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      k_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      i_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      j_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      k_IP1(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      i_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      j_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
      k_IP2(   1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    =   -1000
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_mirror, size(x_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_mirror, size(y_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_mirror, size(z_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP1, size(x_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP1, size(y_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP1, size(z_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP2, size(x_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP2, size(y_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP2, size(z_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_mirror, size(i_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_mirror, size(j_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_mirror, size(k_mirror), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
      call mirrorpoints_ijk(nabs_surf,x_mirror,y_mirror,z_mirror, &
                            x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                            i_mirror,j_mirror,k_mirror, &
                            i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                            dl,zf)

      i_mirror(:,:,0)        = i_mirror(:,:,1)
      i_mirror(:,:,-1)       = i_mirror(:,:,1)
      i_mirror(:,:,-2)       = i_mirror(:,:,1)
      i_mirror(:,:,-3)       = i_mirror(:,:,1)
      i_mirror(:,:,-4)       = i_mirror(:,:,1)
      i_mirror(:,:,-5)       = i_mirror(:,:,1)
      i_mirror(:,:,k1)       = i_mirror(:,:,n(3))
      i_mirror(:,:,k1+1)     = i_mirror(:,:,n(3))
      i_mirror(:,:,k1+2)     = i_mirror(:,:,n(3))
      i_mirror(:,:,k1+3)     = i_mirror(:,:,n(3))
      i_mirror(:,:,k1+4)     = i_mirror(:,:,n(3))
      i_mirror(:,:,k1+5)     = i_mirror(:,:,n(3))
  
      i_IP1(:,:,0)           = i_IP1(:,:,1)
      i_IP1(:,:,-1)          = i_IP1(:,:,1)
      i_IP1(:,:,-2)          = i_IP1(:,:,1)
      i_IP1(:,:,-3)          = i_IP1(:,:,1)
      i_IP1(:,:,-4)          = i_IP1(:,:,1)
      i_IP1(:,:,-5)          = i_IP1(:,:,1)
      i_IP1(:,:,k1)          = i_IP1(:,:,n(3))
      i_IP1(:,:,k1+1)        = i_IP1(:,:,n(3))
      i_IP1(:,:,k1+2)        = i_IP1(:,:,n(3)) 
      i_IP1(:,:,k1+3)        = i_IP1(:,:,n(3))
      i_IP1(:,:,k1+4)        = i_IP1(:,:,n(3))
      i_IP1(:,:,k1+5)        = i_IP1(:,:,n(3))
  
      i_IP2(:,:,0)           = i_IP2(:,:,1)
      i_IP2(:,:,-1)          = i_IP2(:,:,1)
      i_IP2(:,:,-2)          = i_IP2(:,:,1)
      i_IP2(:,:,-3)          = i_IP2(:,:,1)
      i_IP2(:,:,-4)          = i_IP2(:,:,1)
      i_IP2(:,:,-5)          = i_IP2(:,:,1)
      i_IP2(:,:,k1)          = i_IP2(:,:,n(3))
      i_IP2(:,:,k1+1)        = i_IP2(:,:,n(3))
      i_IP2(:,:,k1+2)        = i_IP2(:,:,n(3)) 
      i_IP2(:,:,k1+3)        = i_IP2(:,:,n(3))
      i_IP2(:,:,k1+4)        = i_IP2(:,:,n(3))
      i_IP2(:,:,k1+5)        = i_IP2(:,:,n(3))
  
      j_mirror(:,:,0)        = j_mirror(:,:,1)
      j_mirror(:,:,-1)       = j_mirror(:,:,1)
      j_mirror(:,:,-2)       = j_mirror(:,:,1)
      j_mirror(:,:,-3)       = j_mirror(:,:,1)
      j_mirror(:,:,-4)       = j_mirror(:,:,1)
      j_mirror(:,:,-5)       = j_mirror(:,:,1)
      j_mirror(:,:,k1)       = j_mirror(:,:,n(3))
      j_mirror(:,:,k1+1)     = j_mirror(:,:,n(3))
      j_mirror(:,:,k1+2)     = j_mirror(:,:,n(3))
      j_mirror(:,:,k1+3)     = j_mirror(:,:,n(3))
      j_mirror(:,:,k1+4)     = j_mirror(:,:,n(3))
      j_mirror(:,:,k1+5)     = j_mirror(:,:,n(3))
  
  
      j_IP1(:,:,0)           = j_IP1(:,:,1)
      j_IP1(:,:,-1)          = j_IP1(:,:,1)
      j_IP1(:,:,-2)          = j_IP1(:,:,1)
      j_IP1(:,:,-3)          = j_IP1(:,:,1)
      j_IP1(:,:,-4)          = j_IP1(:,:,1)
      j_IP1(:,:,-5)          = j_IP1(:,:,1)
      j_IP1(:,:,k1)          = j_IP1(:,:,n(3))
      j_IP1(:,:,k1+1)        = j_IP1(:,:,n(3))
      j_IP1(:,:,k1+2)        = j_IP1(:,:,n(3))    
      j_IP1(:,:,k1+3)        = j_IP1(:,:,n(3))
      j_IP1(:,:,k1+4)        = j_IP1(:,:,n(3))
      j_IP1(:,:,k1+5)        = j_IP1(:,:,n(3))
  
      j_IP2(:,:,0)           = j_IP2(:,:,1)
      j_IP2(:,:,-1)          = j_IP2(:,:,1)
      j_IP2(:,:,-2)          = j_IP2(:,:,1)
      j_IP2(:,:,-3)          = j_IP2(:,:,1)
      j_IP2(:,:,-4)          = j_IP2(:,:,1)
      j_IP2(:,:,-5)          = j_IP2(:,:,1)
      j_IP2(:,:,k1)          = j_IP2(:,:,n(3))
      j_IP2(:,:,k1+1)        = j_IP2(:,:,n(3))
      j_IP2(:,:,k1+2)        = j_IP2(:,:,n(3))   
      j_IP2(:,:,k1+3)        = j_IP2(:,:,n(3))
      j_IP2(:,:,k1+4)        = j_IP2(:,:,n(3))
      j_IP2(:,:,k1+5)        = j_IP2(:,:,n(3))
  
      k_mirror(:,:,0)        = k_mirror(:,:,1)
      k_mirror(:,:,-1)       = k_mirror(:,:,1)
      k_mirror(:,:,-2)       = k_mirror(:,:,1)
      k_mirror(:,:,-3)       = k_mirror(:,:,1)
      k_mirror(:,:,-4)       = k_mirror(:,:,1)
      k_mirror(:,:,-5)       = k_mirror(:,:,1)
      k_mirror(:,:,k1)       = k_mirror(:,:,n(3))
      k_mirror(:,:,k1+1)     = k_mirror(:,:,n(3))
      k_mirror(:,:,k1+2)     = k_mirror(:,:,n(3))
      k_mirror(:,:,k1+3)     = k_mirror(:,:,n(3))
      k_mirror(:,:,k1+4)     = k_mirror(:,:,n(3))
      k_mirror(:,:,k1+5)     = k_mirror(:,:,n(3))
  
      k_IP1(:,:,0)           = k_IP1(:,:,1)
      k_IP1(:,:,-1)          = k_IP1(:,:,1)
      k_IP1(:,:,-2)          = k_IP1(:,:,1)
      k_IP1(:,:,-3)          = k_IP1(:,:,1)
      k_IP1(:,:,-4)          = k_IP1(:,:,1)
      k_IP1(:,:,-5)          = k_IP1(:,:,1)
      k_IP1(:,:,k1)          = k_IP1(:,:,n(3))
      k_IP1(:,:,k1+1)        = k_IP1(:,:,n(3))
      k_IP1(:,:,k1+2)        = k_IP1(:,:,n(3))    
      k_IP1(:,:,k1+3)        = k_IP1(:,:,n(3))
      k_IP1(:,:,k1+4)        = k_IP1(:,:,n(3))
      k_IP1(:,:,k1+5)        = k_IP1(:,:,n(3))
  
      k_IP2(:,:,0)           = k_IP2(:,:,1)
      k_IP2(:,:,-1)          = k_IP2(:,:,1)
      k_IP2(:,:,-2)          = k_IP2(:,:,1)
      k_IP2(:,:,-3)          = k_IP2(:,:,1)
      k_IP2(:,:,-4)          = k_IP2(:,:,1)
      k_IP2(:,:,-5)          = k_IP2(:,:,1)
      k_IP2(:,:,k1)          = k_IP2(:,:,n(3))
      k_IP2(:,:,k1+1)        = k_IP2(:,:,n(3))
      k_IP2(:,:,k1+2)        = k_IP2(:,:,n(3))   
      k_IP2(:,:,k1+3)        = k_IP2(:,:,n(3))
      k_IP2(:,:,k1+4)        = k_IP2(:,:,n(3))
      k_IP2(:,:,k1+5)        = k_IP2(:,:,n(3))

  if (myid.eq.0)  print*, '*** Mirror points calculated ***'
  !------------------------------------------------------------
  WP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,:) = 0.0_rp
  WP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,:) = 0.0_rp
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP1, size(x_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP1, size(y_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP1, size(z_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP2, size(x_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP2, size(y_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP2, size(z_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1, size(WP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2, size(WP2), cudaCpuDeviceId, 0)
  call InterpolationWeights(nabs_surf,Level_set, &
                            x_IP1,y_IP1,z_IP1, &
                            x_IP2,y_IP2,z_IP2, &
                            i_IP1,j_IP1,k_IP1, &
                            i_IP2,j_IP2,k_IP2, &
                            WP1,WP2, &
                            dl,dzc)

  WP1(:,:,0,1)           = WP1(:,:,1,1)
  WP1(:,:,-1,1)          = WP1(:,:,1,1)
  WP1(:,:,-2,1)          = WP1(:,:,1,1)
  WP1(:,:,-3,1)          = WP1(:,:,1,1)
  WP1(:,:,-4,1)          = WP1(:,:,1,1)
  WP1(:,:,-5,1)          = WP1(:,:,1,1)
  WP1(:,:,k1,1)          = WP1(:,:,n(3),1)
  WP1(:,:,k1+1,1)        = WP1(:,:,n(3),1)
  WP1(:,:,k1+2,1)        = WP1(:,:,n(3),1)
  WP1(:,:,k1+3,1)        = WP1(:,:,n(3),1)
  WP1(:,:,k1+4,1)        = WP1(:,:,n(3),1)
  WP1(:,:,k1+5,1)        = WP1(:,:,n(3),1)
  
  WP1(:,:,0,2)           = WP1(:,:,1,2)
  WP1(:,:,-1,2)          = WP1(:,:,1,2)
  WP1(:,:,-2,2)          = WP1(:,:,1,2)
  WP1(:,:,-3,2)          = WP1(:,:,1,2)
  WP1(:,:,-4,2)          = WP1(:,:,1,2)
  WP1(:,:,-5,2)          = WP1(:,:,1,2)
  WP1(:,:,k1,2)          = WP1(:,:,n(3),2)
  WP1(:,:,k1+1,2)        = WP1(:,:,n(3),2)
  WP1(:,:,k1+2,2)        = WP1(:,:,n(3),2)
  WP1(:,:,k1+3,2)        = WP1(:,:,n(3),2)
  WP1(:,:,k1+4,2)        = WP1(:,:,n(3),2)
  WP1(:,:,k1+5,2)        = WP1(:,:,n(3),2)
  
  WP1(:,:,0,3)           = WP1(:,:,1,3)
  WP1(:,:,-1,3)          = WP1(:,:,1,3)
  WP1(:,:,-2,3)          = WP1(:,:,1,3)
  WP1(:,:,-3,3)          = WP1(:,:,1,3)
  WP1(:,:,-4,3)          = WP1(:,:,1,3)
  WP1(:,:,-5,3)          = WP1(:,:,1,3)
  WP1(:,:,k1,3)          = WP1(:,:,n(3),3)
  WP1(:,:,k1+1,3)        = WP1(:,:,n(3),3)
  WP1(:,:,k1+2,3)        = WP1(:,:,n(3),3)
  WP1(:,:,k1+3,3)        = WP1(:,:,n(3),3)
  WP1(:,:,k1+4,3)        = WP1(:,:,n(3),3)
  WP1(:,:,k1+5,3)        = WP1(:,:,n(3),3)

  WP1(:,:,0,4)           = WP1(:,:,1,4)
  WP1(:,:,-1,4)          = WP1(:,:,1,4)
  WP1(:,:,-2,4)          = WP1(:,:,1,4)
  WP1(:,:,-3,4)          = WP1(:,:,1,4)
  WP1(:,:,-4,4)          = WP1(:,:,1,4)
  WP1(:,:,-5,4)          = WP1(:,:,1,4)
  WP1(:,:,k1,4)          = WP1(:,:,n(3),4)
  WP1(:,:,k1+1,4)        = WP1(:,:,n(3),4)
  WP1(:,:,k1+2,4)        = WP1(:,:,n(3),4)
  WP1(:,:,k1+3,4)        = WP1(:,:,n(3),4)
  WP1(:,:,k1+4,4)        = WP1(:,:,n(3),4)
  WP1(:,:,k1+5,4)        = WP1(:,:,n(3),4)

  WP1(:,:,0,5)           = WP1(:,:,1,5)
  WP1(:,:,-1,5)          = WP1(:,:,1,5)
  WP1(:,:,-2,5)          = WP1(:,:,1,5)
  WP1(:,:,-3,5)          = WP1(:,:,1,5)
  WP1(:,:,-4,5)          = WP1(:,:,1,5)
  WP1(:,:,-5,5)          = WP1(:,:,1,5)
  WP1(:,:,k1,5)          = WP1(:,:,n(3),5)
  WP1(:,:,k1+1,5)        = WP1(:,:,n(3),5)
  WP1(:,:,k1+2,5)        = WP1(:,:,n(3),5)
  WP1(:,:,k1+3,5)        = WP1(:,:,n(3),5)
  WP1(:,:,k1+4,5)        = WP1(:,:,n(3),5)
  WP1(:,:,k1+5,5)        = WP1(:,:,n(3),5)
  
  WP1(:,:,0,6)           = WP1(:,:,1,6)
  WP1(:,:,-1,6)          = WP1(:,:,1,6)
  WP1(:,:,-2,6)          = WP1(:,:,1,6)
  WP1(:,:,-3,6)          = WP1(:,:,1,6)
  WP1(:,:,-4,6)          = WP1(:,:,1,6)
  WP1(:,:,-5,6)          = WP1(:,:,1,6)
  WP1(:,:,k1,6)          = WP1(:,:,n(3),6)
  WP1(:,:,k1+1,6)        = WP1(:,:,n(3),6)
  WP1(:,:,k1+2,6)        = WP1(:,:,n(3),6)
  WP1(:,:,k1+3,6)        = WP1(:,:,n(3),6)
  WP1(:,:,k1+4,6)        = WP1(:,:,n(3),6)
  WP1(:,:,k1+5,6)        = WP1(:,:,n(3),6)

  WP1(:,:,0,7)           = WP1(:,:,1,7)
  WP1(:,:,-1,7)          = WP1(:,:,1,7)
  WP1(:,:,-2,7)          = WP1(:,:,1,7)
  WP1(:,:,-3,7)          = WP1(:,:,1,7)
  WP1(:,:,-4,7)          = WP1(:,:,1,7)
  WP1(:,:,-5,7)          = WP1(:,:,1,7)
  WP1(:,:,k1,7)          = WP1(:,:,n(3),7)
  WP1(:,:,k1+1,7)        = WP1(:,:,n(3),7)
  WP1(:,:,k1+2,7)        = WP1(:,:,n(3),7)
  WP1(:,:,k1+3,7)        = WP1(:,:,n(3),7)
  WP1(:,:,k1+4,7)        = WP1(:,:,n(3),7)
  WP1(:,:,k1+5,7)        = WP1(:,:,n(3),7)

  WP2(:,:,0,1)           = WP2(:,:,1,2)
  WP2(:,:,-1,1)          = WP2(:,:,1,2)
  WP2(:,:,-2,1)          = WP2(:,:,1,2)
  WP2(:,:,-3,1)          = WP2(:,:,1,2)
  WP2(:,:,-4,1)          = WP2(:,:,1,2)
  WP2(:,:,-5,1)          = WP2(:,:,1,2)
  WP2(:,:,k1,1)          = WP2(:,:,n(3),2)
  WP2(:,:,k1+1,1)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+2,1)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+3,1)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+4,1)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+5,1)        = WP2(:,:,n(3),2)
  
  WP2(:,:,0,2)           = WP2(:,:,1,2)
  WP2(:,:,-1,2)          = WP2(:,:,1,2)
  WP2(:,:,-2,2)          = WP2(:,:,1,2)
  WP2(:,:,-3,2)          = WP2(:,:,1,2)
  WP2(:,:,-4,2)          = WP2(:,:,1,2)
  WP2(:,:,-5,2)          = WP2(:,:,1,2)
  WP2(:,:,k1,2)          = WP2(:,:,n(3),2)
  WP2(:,:,k1+1,2)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+2,2)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+3,2)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+4,2)        = WP2(:,:,n(3),2)
  WP2(:,:,k1+5,2)        = WP2(:,:,n(3),2)
  
  WP2(:,:,0,3)           = WP2(:,:,1,3)
  WP2(:,:,-1,3)          = WP2(:,:,1,3)
  WP2(:,:,-2,3)          = WP2(:,:,1,3)
  WP2(:,:,-3,3)          = WP2(:,:,1,3)
  WP2(:,:,-4,3)          = WP2(:,:,1,3)
  WP2(:,:,-5,3)          = WP2(:,:,1,3)
  WP2(:,:,k1,3)          = WP2(:,:,n(3),3)
  WP2(:,:,k1+1,3)        = WP2(:,:,n(3),3)
  WP2(:,:,k1+2,3)        = WP2(:,:,n(3),3)
  WP2(:,:,k1+3,3)        = WP2(:,:,n(3),3)
  WP2(:,:,k1+4,3)        = WP2(:,:,n(3),3)
  WP2(:,:,k1+5,3)        = WP2(:,:,n(3),3)

  WP2(:,:,0,4)           = WP2(:,:,1,4)
  WP2(:,:,-1,4)          = WP2(:,:,1,4)
  WP2(:,:,-2,4)          = WP2(:,:,1,4)
  WP2(:,:,-3,4)          = WP2(:,:,1,4)
  WP2(:,:,-4,4)          = WP2(:,:,1,4)
  WP2(:,:,-5,4)          = WP2(:,:,1,4)
  WP2(:,:,k1,4)          = WP2(:,:,n(3),4)
  WP2(:,:,k1+1,4)        = WP2(:,:,n(3),4)
  WP2(:,:,k1+2,4)        = WP2(:,:,n(3),4)
  WP2(:,:,k1+3,4)        = WP2(:,:,n(3),4)
  WP2(:,:,k1+4,4)        = WP2(:,:,n(3),4)
  WP2(:,:,k1+5,4)        = WP2(:,:,n(3),4)

  WP2(:,:,0,5)           = WP2(:,:,1,5)
  WP2(:,:,-1,5)          = WP2(:,:,1,5)
  WP2(:,:,-2,5)          = WP2(:,:,1,5)
  WP2(:,:,-3,5)          = WP2(:,:,1,5)
  WP2(:,:,-4,5)          = WP2(:,:,1,5)
  WP2(:,:,-5,5)          = WP2(:,:,1,5)
  WP2(:,:,k1,5)          = WP2(:,:,n(3),5)
  WP2(:,:,k1+1,5)        = WP2(:,:,n(3),5)
  WP2(:,:,k1+2,5)        = WP2(:,:,n(3),5)
  WP2(:,:,k1+3,5)        = WP2(:,:,n(3),5)
  WP2(:,:,k1+4,5)        = WP2(:,:,n(3),5)
  WP2(:,:,k1+5,5)        = WP2(:,:,n(3),5)
  
  WP2(:,:,0,6)           = WP2(:,:,1,6)
  WP2(:,:,-1,6)          = WP2(:,:,1,6)
  WP2(:,:,-2,6)          = WP2(:,:,1,6)
  WP2(:,:,-3,6)          = WP2(:,:,1,6)
  WP2(:,:,-4,6)          = WP2(:,:,1,6)
  WP2(:,:,-5,6)          = WP2(:,:,1,6)
  WP2(:,:,k1,6)          = WP2(:,:,n(3),6)
  WP2(:,:,k1+1,6)        = WP2(:,:,n(3),6)
  WP2(:,:,k1+2,6)        = WP2(:,:,n(3),6)
  WP2(:,:,k1+3,6)        = WP2(:,:,n(3),6)
  WP2(:,:,k1+4,6)        = WP2(:,:,n(3),6)
  WP2(:,:,k1+5,6)        = WP2(:,:,n(3),6)

  WP2(:,:,0,7)           = WP2(:,:,1,7)
  WP2(:,:,-1,7)          = WP2(:,:,1,7)
  WP2(:,:,-2,7)          = WP2(:,:,1,7)
  WP2(:,:,-3,7)          = WP2(:,:,1,7)
  WP2(:,:,-4,7)          = WP2(:,:,1,7)
  WP2(:,:,-5,7)          = WP2(:,:,1,7)
  WP2(:,:,k1,7)          = WP2(:,:,n(3),7)
  WP2(:,:,k1+1,7)        = WP2(:,:,n(3),7)
  WP2(:,:,k1+2,7)        = WP2(:,:,n(3),7)
  WP2(:,:,k1+3,7)        = WP2(:,:,n(3),7)
  WP2(:,:,k1+4,7)        = WP2(:,:,n(3),7)
  WP2(:,:,k1+5,7)        = WP2(:,:,n(3),7)

  if (myid.eq.0)  print*, '*** Interpolation Weights calculated ***'
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1, size(WP1), cudaCpuDeviceId, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2, size(WP2), cudaCpuDeviceId, 0)
  call load('w',trim(restart_dir)//'utag.bin',n,      cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'vtag.bin',n,      cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'wtag.bin',n,      cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'stag.bin',n,      cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'lvl.bin',n,      Level_set(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nx.bin',n,      nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ny.bin',n,      ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nz.bin',n,      nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ns.bin',n,      nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'deltan.bin',n,      deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP1.bin',n,      i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP1.bin',n,      j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP1.bin',n,      k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP2.bin',n,      i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP2.bin',n,      j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP2.bin',n,      k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('w',trim(restart_dir)//'WP1.bin',n,      WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('w',trim(restart_dir)//'WP2.bin',n,      WP2(1:n(1),1:n(2),1:n(3),1:7))

else

  call load('r',trim(restart_dir)//'utag.bin',n,      cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'vtag.bin',n,      cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'wtag.bin',n,      cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'stag.bin',n,      cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'lvl.bin',n,      Level_set(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nx.bin',n,      nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ny.bin',n,      ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nz.bin',n,      nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ns.bin',n,      nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'deltan.bin',n,      deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP1.bin',n,      i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP1.bin',n,      j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP1.bin',n,      k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP2.bin',n,      i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP2.bin',n,      j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP2.bin',n,      k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('r',trim(restart_dir)//'WP1.bin',n,      WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('r',trim(restart_dir)//'WP2.bin',n,      WP2(1:n(1),1:n(2),1:n(3),1:7))
  !@cuf istat = cudaMemPrefetchAsync(cell_u_tag, size(cell_u_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_v_tag, size(cell_v_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_w_tag, size(cell_w_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(cell_phi_tag, size(cell_phi_tag), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(Level_set, size(Level_set), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nx_surf, size(nx_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(ny_surf, size(ny_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nz_surf, size(nz_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(nabs_surf, size(nabs_surf), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(deltan, size(deltan), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP1, size(i_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP1, size(j_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP1, size(k_IP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(i_IP2, size(i_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(j_IP2, size(j_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(k_IP2, size(k_IP2), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP1, size(WP1), mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(WP2, size(WP2), mydev, 0)
   call bounduvw(cbcvel,n,bcvel,nh_d,nh_b,halo,is_outflow,dl,dzc,dzf,cell_u_tag,cell_v_tag,cell_w_tag)
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_phi_tag)
   temp(:,:,:) = real(Level_set(:,:,:),rp)
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,temp)
   Level_set(:,:,:) = int(temp(:,:,:))
   call bounduvw(cbcvel,n,bcvel,nh_d,nh_b,halo,is_outflow,dl,dzc,dzf,nx_surf,ny_surf,nz_surf)
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nabs_surf)

   i_mirror(:,:,0)        = i_mirror(:,:,1)
   i_mirror(:,:,-1)       = i_mirror(:,:,1)
   i_mirror(:,:,-2)       = i_mirror(:,:,1)
   i_mirror(:,:,-3)       = i_mirror(:,:,1)
   i_mirror(:,:,-4)       = i_mirror(:,:,1)
   i_mirror(:,:,-5)       = i_mirror(:,:,1)
   i_mirror(:,:,k1)       = i_mirror(:,:,n(3))
   i_mirror(:,:,k1+1)     = i_mirror(:,:,n(3))
   i_mirror(:,:,k1+2)     = i_mirror(:,:,n(3))
   i_mirror(:,:,k1+3)     = i_mirror(:,:,n(3))
   i_mirror(:,:,k1+4)     = i_mirror(:,:,n(3))
   i_mirror(:,:,k1+5)     = i_mirror(:,:,n(3))
  
   j_mirror(:,:,0)        = j_mirror(:,:,1)
   j_mirror(:,:,-1)       = j_mirror(:,:,1)
   j_mirror(:,:,-2)       = j_mirror(:,:,1)
   j_mirror(:,:,-3)       = j_mirror(:,:,1)
   j_mirror(:,:,-4)       = j_mirror(:,:,1)
   j_mirror(:,:,-5)       = j_mirror(:,:,1)
   j_mirror(:,:,k1)       = j_mirror(:,:,n(3))
   j_mirror(:,:,k1+1)     = j_mirror(:,:,n(3))
   j_mirror(:,:,k1+2)     = j_mirror(:,:,n(3))
   j_mirror(:,:,k1+3)     = j_mirror(:,:,n(3))
   j_mirror(:,:,k1+4)     = j_mirror(:,:,n(3))
   j_mirror(:,:,k1+5)     = j_mirror(:,:,n(3))
  
   k_mirror(:,:,0)        = k_mirror(:,:,1)
   k_mirror(:,:,-1)       = k_mirror(:,:,1)
   k_mirror(:,:,-2)       = k_mirror(:,:,1)
   k_mirror(:,:,-3)       = k_mirror(:,:,1)
   k_mirror(:,:,-4)       = k_mirror(:,:,1)
   k_mirror(:,:,-5)       = k_mirror(:,:,1)
   k_mirror(:,:,k1)       = k_mirror(:,:,n(3))
   k_mirror(:,:,k1+1)     = k_mirror(:,:,n(3))
   k_mirror(:,:,k1+2)     = k_mirror(:,:,n(3))
   k_mirror(:,:,k1+3)     = k_mirror(:,:,n(3))
   k_mirror(:,:,k1+4)     = k_mirror(:,:,n(3))
   k_mirror(:,:,k1+5)     = k_mirror(:,:,n(3))

   i_IP1(:,:,0)           = i_IP1(:,:,1)
   i_IP1(:,:,-1)          = i_IP1(:,:,1)
   i_IP1(:,:,-2)          = i_IP1(:,:,1)
   i_IP1(:,:,-3)          = i_IP1(:,:,1)
   i_IP1(:,:,-4)          = i_IP1(:,:,1)
   i_IP1(:,:,-5)          = i_IP1(:,:,1)
   i_IP1(:,:,k1)          = i_IP1(:,:,n(3))
   i_IP1(:,:,k1+1)        = i_IP1(:,:,n(3))
   i_IP1(:,:,k1+2)        = i_IP1(:,:,n(3))
   i_IP1(:,:,k1+3)        = i_IP1(:,:,n(3))
   i_IP1(:,:,k1+4)        = i_IP1(:,:,n(3))
   i_IP1(:,:,k1+5)        = i_IP1(:,:,n(3))

   i_IP2(:,:,0)           = i_IP2(:,:,1)
   i_IP2(:,:,-1)          = i_IP2(:,:,1)
   i_IP2(:,:,-2)          = i_IP2(:,:,1)
   i_IP2(:,:,-3)          = i_IP2(:,:,1)
   i_IP2(:,:,-4)          = i_IP2(:,:,1)
   i_IP2(:,:,-5)          = i_IP2(:,:,1)
   i_IP2(:,:,k1)          = i_IP2(:,:,n(3))
   i_IP2(:,:,k1+1)        = i_IP2(:,:,n(3))
   i_IP2(:,:,k1+2)        = i_IP2(:,:,n(3))
   i_IP2(:,:,k1+3)        = i_IP2(:,:,n(3))
   i_IP2(:,:,k1+4)        = i_IP2(:,:,n(3))
   i_IP2(:,:,k1+5)        = i_IP2(:,:,n(3))

   j_IP1(:,:,0)           = j_IP1(:,:,1)
   j_IP1(:,:,-1)          = j_IP1(:,:,1)
   j_IP1(:,:,-2)          = j_IP1(:,:,1)
   j_IP1(:,:,-3)          = j_IP1(:,:,1)
   j_IP1(:,:,-4)          = j_IP1(:,:,1)
   j_IP1(:,:,-5)          = j_IP1(:,:,1)
   j_IP1(:,:,k1)          = j_IP1(:,:,n(3))
   j_IP1(:,:,k1+1)        = j_IP1(:,:,n(3))
   j_IP1(:,:,k1+2)        = j_IP1(:,:,n(3))
   j_IP1(:,:,k1+3)        = j_IP1(:,:,n(3))
   j_IP1(:,:,k1+4)        = j_IP1(:,:,n(3))
   j_IP1(:,:,k1+5)        = j_IP1(:,:,n(3))

   j_IP2(:,:,0)           = j_IP2(:,:,1)
   j_IP2(:,:,-1)          = j_IP2(:,:,1)
   j_IP2(:,:,-2)          = j_IP2(:,:,1)
   j_IP2(:,:,-3)          = j_IP2(:,:,1)
   j_IP2(:,:,-4)          = j_IP2(:,:,1)
   j_IP2(:,:,-5)          = j_IP2(:,:,1)
   j_IP2(:,:,k1)          = j_IP2(:,:,n(3))
   j_IP2(:,:,k1+1)        = j_IP2(:,:,n(3))
   j_IP2(:,:,k1+2)        = j_IP2(:,:,n(3))
   j_IP2(:,:,k1+3)        = j_IP2(:,:,n(3))
   j_IP2(:,:,k1+4)        = j_IP2(:,:,n(3))
   j_IP2(:,:,k1+5)        = j_IP2(:,:,n(3))


   k_IP1(:,:,0)           = k_IP1(:,:,1)
   k_IP1(:,:,-1)          = k_IP1(:,:,1)
   k_IP1(:,:,-2)          = k_IP1(:,:,1)
   k_IP1(:,:,-3)          = k_IP1(:,:,1)
   k_IP1(:,:,-4)          = k_IP1(:,:,1)
   k_IP1(:,:,-5)          = k_IP1(:,:,1)
   k_IP1(:,:,k1)          = k_IP1(:,:,n(3))
   k_IP1(:,:,k1+1)        = k_IP1(:,:,n(3))
   k_IP1(:,:,k1+2)        = k_IP1(:,:,n(3))
   k_IP1(:,:,k1+3)        = k_IP1(:,:,n(3))
   k_IP1(:,:,k1+4)        = k_IP1(:,:,n(3))
   k_IP1(:,:,k1+5)        = k_IP1(:,:,n(3))

   k_IP2(:,:,0)           = k_IP2(:,:,1)
   k_IP2(:,:,-1)          = k_IP2(:,:,1)
   k_IP2(:,:,-2)          = k_IP2(:,:,1)
   k_IP2(:,:,-3)          = k_IP2(:,:,1)
   k_IP2(:,:,-4)          = k_IP2(:,:,1)
   k_IP2(:,:,-5)          = k_IP2(:,:,1)
   k_IP2(:,:,k1)          = k_IP2(:,:,n(3))
   k_IP2(:,:,k1+1)        = k_IP2(:,:,n(3))
   k_IP2(:,:,k1+2)        = k_IP2(:,:,n(3))
   k_IP2(:,:,k1+3)        = k_IP2(:,:,n(3))
   k_IP2(:,:,k1+4)        = k_IP2(:,:,n(3))
   k_IP2(:,:,k1+5)        = k_IP2(:,:,n(3))

  if (myid.eq.0)  print*, '*** Existing IBM data loaded ***'

endif

  deallocate(i_mirror,j_mirror,k_mirror,x_mirror,y_mirror,z_mirror, &
             z_intersect,y_intersect,x_intersect,&
             x_IP1,x_IP2,y_IP1,y_IP2,z_IP1,z_IP2)
return
end subroutine initIBM
#endif
end module mod_initIBM