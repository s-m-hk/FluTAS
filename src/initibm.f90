module mod_initIBM
#if defined(_USE_IBM)
use mod_IBM
use mod_bound, only: boundp,bounduvw
use mod_load,  only: load,load_int,load_weight
use mod_param, only: l,ng,cbcpre,bcpre,cbcvel,bcvel,is_outflow
use mod_common_mpi, only: myid,ierr
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
                   nh_d,nh_b,halo,&
                   zc,zf,dzc,dzf,dl,dli)
implicit none
!
character(len=100) :: datadir,restart_dir
integer,  intent(in)                         :: nh_d,nh_b
integer,  dimension(3), intent(in)           :: dims,n,ng,halo
real(rp), dimension(1-nh_d:), intent(in)     :: zc,zf,dzc,dzf
real(rp), dimension(3), intent(in)           :: dl,dli
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(inout) :: Level_set,cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
! integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(inout) :: Level_set
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(inout) :: i_mirror,j_mirror,k_mirror
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(inout) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b)    , intent(inout) :: nx_surf,ny_surf,nz_surf,nabs_surf,deltan
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7), intent(inout) :: WP1,WP2
real(rp), dimension(:,:,:), allocatable :: x_intersect,y_intersect,z_intersect
real(rp), dimension(:,:,:), allocatable :: x_IP1,y_IP1,z_IP1
real(rp), dimension(:,:,:), allocatable :: x_IP2,y_IP2,z_IP2
real(rp), dimension(:,:,:), allocatable :: x_mirror,y_mirror,z_mirror
logical:: is_data
integer::i,j,k,q,n1,n2,n3
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
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
restart_dir = 'data/restart_dir/'
!
inquire(file=trim(restart_dir)//'stag.bin',exist=is_data)
!
if (.not.is_data) then
  allocate(x_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_mirror(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_IP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_IP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))
  allocate(x_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           y_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), &
           z_intersect(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b))

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
  !@cuf istat = cudaMemPrefetchAsync(x_intersect, size(x_intersect),   mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_intersect, size(y_intersect),   mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_intersect, size(z_intersect),   mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP1, size(x_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP1, size(y_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP1, size(z_IP1),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_IP2, size(x_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_IP2, size(y_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_IP2, size(z_IP2),               mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(x_mirror, size(x_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(y_mirror, size(y_mirror),         mydev, 0)
  !@cuf istat = cudaMemPrefetchAsync(z_mirror, size(z_mirror),         mydev, 0)

  call IBM_mask(dims,n,ng,nh_d,nh_b,cell_phi_tag,Level_set,zc,zf,dl,dzc)

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,Level_set)

  if (myid.eq.0)  print*, '*** Solid marker set ***'
  !---------------------------------------------------------------------

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_phi_tag)

   cell_u_tag(:,:,:)   = 1.0_rp
   cell_v_tag(:,:,:)   = 1.0_rp
   cell_w_tag(:,:,:)   = 1.0_rp

   do k=1,n3
    do j=1,n2
     do i=1,n1
      cell_u_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i-1,j,k)+cell_phi_tag(i,j,k))
      cell_v_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j-1,k)+cell_phi_tag(i,j,k))
      cell_w_tag(i,j,k) = 0.5_rp*(cell_phi_tag(i,j,k-1)+cell_phi_tag(i,j,k))
     enddo
    enddo
   enddo

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_u_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_v_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_w_tag)

  if (myid.eq.0)  print*, '*** Volume fractions calculated ***'
  !---------------------------------------------------------------------

  call normal_vectors(dims,n,ng,nh_d,nh_b,Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,dl,dli,zc,dzc)

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nabs_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nx_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,ny_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nz_surf)

  if (myid.eq.0)  print*, '*** Normal vectors calculated ***'
  !*********************************************************************

  call intersect(dims,n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dl,dzc,zc,zf)

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_intersect)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_intersect)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_intersect)
  
  if (myid.eq.0)  print*, '*** Intersection points calculated ***'
  !*********************************************************************

  call mirrorpoints(dims,n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf, &
                    x_intersect,y_intersect,z_intersect, &
                    x_mirror,y_mirror,z_mirror, &
                    deltan, &
                    x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                    dl,dzc,zc)

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,deltan)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_mirror)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_mirror)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_mirror)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_IP1)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_IP1)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_IP1)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,x_IP2)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,y_IP2)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,z_IP2)

  call mirrorpoints_ijk(dims,n,ng,nh_d,nh_b,nabs_surf,x_mirror,y_mirror,z_mirror, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        i_mirror,j_mirror,k_mirror, &
                        i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                        dl,zf)

  if (myid.eq.0)  print*, '*** Mirror points calculated ***'
  !------------------------------------------------------------
  
  call InterpolationWeights(dims,n,ng,nh_d,nh_b,nabs_surf,Level_set, &
                            x_IP1,y_IP1,z_IP1, &
                            x_IP2,y_IP2,z_IP2, &
                            i_IP1,j_IP1,k_IP1, &
                            i_IP2,j_IP2,k_IP2, &
                            WP1,WP2, &
                            dl,dzc,zc)

   do q = 1,7
     call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,WP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,q))
     call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,WP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,q))
   enddo

  if (myid.eq.0)  print*, '*** Interpolation Weights calculated ***'

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
  
  call load('w',trim(restart_dir)//'utag.bin',n,         cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'vtag.bin',n,         cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'wtag.bin',n,         cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'stag.bin',n,         cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'lvl.bin',n,      Level_set(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nx.bin',n,           nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ny.bin',n,           ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'nz.bin',n,           nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'ns.bin',n,           nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('w',trim(restart_dir)//'deltan.bin',n,       deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'imirror.bin',n,  i_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jmirror.bin',n,  j_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kmirror.bin',n,  k_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP1.bin',n,  i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP1.bin',n,  j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP1.bin',n,  k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'iIP2.bin',n,  i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'jIP2.bin',n,  j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('w',trim(restart_dir)//'kIP2.bin',n,  k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('w',trim(restart_dir)//'WP1.bin',n,    WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('w',trim(restart_dir)//'WP2.bin',n,    WP2(1:n(1),1:n(2),1:n(3),1:7))
  !
  deallocate(x_intersect,y_intersect,z_intersect,x_mirror,y_mirror,z_mirror,x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2)
  !
else
  !
  call load('r',trim(restart_dir)//'utag.bin',n,      cell_u_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'vtag.bin',n,      cell_v_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'wtag.bin',n,      cell_w_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'stag.bin',n,      cell_phi_tag(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'lvl.bin',n,   Level_set(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nx.bin',n,      nx_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ny.bin',n,      ny_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'nz.bin',n,      nz_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'ns.bin',n,      nabs_surf(1:n(1),1:n(2),1:n(3)))
  call load('r',trim(restart_dir)//'deltan.bin',n,      deltan(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'imirror.bin',n,      i_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jmirror.bin',n,      j_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kmirror.bin',n,      k_mirror(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP1.bin',n,  i_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP1.bin',n,  j_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP1.bin',n,  k_IP1(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'iIP2.bin',n,  i_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'jIP2.bin',n,  j_IP2(1:n(1),1:n(2),1:n(3)))
  call load_int('r',trim(restart_dir)//'kIP2.bin',n,  k_IP2(1:n(1),1:n(2),1:n(3)))
  call load_weight('r',trim(restart_dir)//'WP1.bin',n,    WP1(1:n(1),1:n(2),1:n(3),1:7))
  call load_weight('r',trim(restart_dir)//'WP2.bin',n,    WP2(1:n(1),1:n(2),1:n(3),1:7))

  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_u_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_v_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_w_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,cell_phi_tag)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,Level_set)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nabs_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nx_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,ny_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,nz_surf)
  call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,deltan)
  do q = 1,7
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,WP1(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,q))
   call boundp(cbcpre,n,bcpre,nh_d,nh_b,halo,dl,dzc,dzf,WP2(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,q))
  enddo
  !*********************************************************************

  if (myid.eq.0)  print*, '****** Existing IBM data loaded ******'

endif
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
end subroutine initIBM
#endif
end module mod_initIBM