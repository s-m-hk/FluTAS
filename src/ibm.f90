module mod_IBM
#if defined(_USE_IBM)
use mpi
use mod_param, only: l,cbcpre, &
                     surface_type, solid_height_ratio, &
                     Rotation_angle, &
                     xc_ibm, yc_ibm, zc_ibm, r_ibm, &
                     zmax_ibm, zmin_ibm, d1, d2, small
use mod_common_mpi, only: myid,ierr,ijk_start
use mod_types
!@cuf use cudafor
!@cuf use mod_common_mpi, only: mydev
implicit none
!
private
public IBM_Mask,normal_vectors,intersect,mirrorpoints, &
       Penalization_center,Penalization_face,interpolation_dphi, &
       mirrorpoints_ijk,interpolation_mirror,InterpolationWeights, &
       Wetting_radius
contains
!
subroutine IBM_Mask(n,ng,nh_d,nh_b,cell_phi_tag,Level_set,zc,zf,dl,dzc)
implicit none
integer, intent(in)                           :: nh_d,nh_b
integer, intent(in), dimension(3)             :: n,ng
real(rp), intent(in), dimension(1-nh_d:)                                               :: zc,zf,dzc
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: cell_phi_tag
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: Level_set
real(rp), intent(in), dimension(3)       :: dl
integer  :: i,j,k,ii,jj,ll,nn,m,number_of_divisions
real(rp) :: xxx,yyy,zzz,dx,dy,dz,dxx,dyy,dzz,lxl,lyl,lzl
real(rp) :: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp) :: RR,x_center,y_center,z_center,z_top,z_bottom
integer, dimension(3) :: iperiod
logical  :: inside,ghost
integer  :: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z,q,iperiod1,iperiod2,iperiod3
real(rp) :: sdist1,sdist2,sdistmin,eps
integer  :: counter
!
!@cuf attributes(managed) :: cell_phi_tag,Level_set,zc,zf,dzc
!
dx=dl(1)
dy=dl(2)
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
lxl = l(1)
lyl = l(2)
lzl = l(3)
!
RR =  r_ibm
x_center = xc_ibm
y_center = yc_ibm
z_center = zc_ibm
z_bottom = zmin_ibm
z_top    = zmax_ibm
!
do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo
iperiod1 = iperiod(1)
iperiod2 = iperiod(2)
iperiod3 = iperiod(3)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
! Wall Geometry
!
#if defined(_OPENACC)
number_of_divisions = 50
#else
number_of_divisions = 50
#endif
if (myid == 0) print*, '*** Calculating volume fractions ***'
!$acc enter data copyin(ijk_start_x,ijk_start_y,ijk_start_z,dx,dy,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
  do j=1,n2
    do i=1,n1
      xxx =  (i+ijk_start_x-0.5_rp)*dx
      yyy =  (j+ijk_start_y-0.5_rp)*dy
      zzz =  zc(k)
      dz  = dzc(k)
#if defined(_OPENACC)
      ! ghost = wall_mounted_cylinder_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                                          ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
      ghost = xGrv_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                         xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
#else
      call GuessGhostCells(xxx,yyy,zzz,dx,dy,dz,lxl,lyl,lzl,ghost)
#endif
      if (ghost) then
       Level_set(i,j,k) = 0
      endif
    enddo
  enddo
enddo
!$acc end kernels
!$acc wait
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!
!$acc enter data copyin(number_of_divisions,ijk_start_x,ijk_start_y,dx,dy,dz,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
  if (myid.eq.0) print*, '*** Calculating volume fractions at k = ', k
#endif
  do j=1,n2
    do i=1,n1
      dz = dzc(k)

      cell_start_x = (i+ijk_start_x-1.0_rp)*dx
      cell_end_x   = (i+ijk_start_x-0.0_rp)*dx

      cell_start_y = (j+ijk_start_y-1.0_rp)*dy
      cell_end_y   = (j+ijk_start_y-0.0_rp)*dy

      cell_start_z = zf(k-1)
      cell_end_z   = zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0
      !$acc loop seq
      do nn = 1,number_of_divisions
         zzz = cell_start_z + (nn-1)*dzz
       !$acc loop seq
       do m  = 1,number_of_divisions
          yyy = cell_start_y + (m-1)*dyy
        !$acc loop seq
        do ll = 1,number_of_divisions
           xxx = cell_start_x + (ll-1)*dxx
#if defined(_OPENACC)
            ! inside = wall_mounted_cylinder(RR,x_center,y_center,z_center,z_bottom,z_top, &
                                           ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
            inside = xGrv(RR,x_center,y_center,z_center,z_bottom,z_top, &
                          xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
#else
            call Solid_Surface(xxx,yyy,zzz,dx,dy,dz,lxl,lyl,lzl,inside)
#endif
            if (inside) counter = counter + 1
        enddo
        !$acc end loop
       enddo
       !$acc end loop
      enddo
      !$acc end loop
      cell_phi_tag(i,j,k) = 1.0_rp - (1._rp*counter)/(1._rp*number_of_divisions**3)
    enddo
  enddo
enddo
!$acc end kernels
!$acc wait
!$acc exit data copyout(number_of_divisions,ijk_start_x,ijk_start_y,dx,dy,dz,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!
end subroutine IBM_Mask

function wall_mounted_cylinder_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                                     xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(wall_mounted_cylinder_ghost) seq
logical :: wall_mounted_cylinder_ghost
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
real(rp):: sdist1,sdist2,sdistmin,eps
integer :: ii,jj

    wall_mounted_cylinder_ghost = .false.
    eps    = min(dx,dy)
    sdist1 = zzz - z_top

    do jj = -1,1
      do ii = -1,1
        sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                       (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
        if((sdist1.lt.-dz).and.(sdist2.lt.-eps)) then
          wall_mounted_cylinder_ghost = .true.
        elseif((sdist1.gt.dz).and.(sdist2.gt.eps)) then
          wall_mounted_cylinder_ghost = .false.
        endif
      enddo
    enddo
    if(zzz.lt.z_bottom) wall_mounted_cylinder_ghost = .true.

end function wall_mounted_cylinder_ghost
!
function wall_mounted_cylinder(RR,x_center,y_center,z_center,z_bottom,z_top, &
                               xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(wall_mounted_cylinder) seq
logical :: wall_mounted_cylinder
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
real(rp):: sdist1,sdist2,sdistmin ,eps
integer :: ii,jj

    wall_mounted_cylinder= .false.
    eps    = min(dx,dy)
    sdist1 = zzz - z_top

    do jj = -1,1
      do ii = -1,1
        sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                       (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
        if((sdist1.le.0.).and.(sdist2.le.0.)) then
          wall_mounted_cylinder = .true.
        elseif((sdist1.gt.0.).and.(sdist2.gt.0.)) then
          wall_mounted_cylinder = .false.
        endif
      enddo
    enddo
    if(zzz.le.z_bottom) wall_mounted_cylinder = .true.

end function wall_mounted_cylinder
!!
function flat_wall_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                         xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(flat_wall_ghost) seq
logical :: flat_wall_ghost
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
    flat_wall_ghost = .false.
    if (zzz.le.z_center) flat_wall_ghost = .true.
end function flat_wall_ghost
!
function flat_wall(RR,x_center,y_center,z_center,z_bottom,z_top, &
                   xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(flat_wall) seq
logical :: flat_wall
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
    flat_wall = .false.
    if (zzz.le.z_center) flat_wall = .true.
end function flat_wall
!!
function xGrv_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                    xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(xGrv_ghost) seq
logical :: xGrv_ghost
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
logical :: cond1,cond2,cond3,cond4
    cond1 = yyy.lt.RR
    cond2 = yyy.gt.(ly - RR)
    cond3 = zzz.lt.z_top
    cond4 = zzz.lt.z_bottom
    xGrv_ghost = .false.
    if ((cond1.and.cond3).or. &
        (cond2.and.cond3).or. &
        ((.not.cond1).and.cond4).or. &
        ((.not.cond2).and.cond4)) xGrv_ghost = .true.
end function xGrv_ghost
!
function xGrv(RR,x_center,y_center,z_center,z_bottom,z_top, &
              xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
!$acc routine(xGrv) seq
logical :: xGrv
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
logical  :: cond1,cond2,cond3,cond4
    cond1 = yyy.le.RR
    cond2 = yyy.ge.(ly - RR)
    cond3 = zzz.le.z_top
    cond4 = zzz.le.z_bottom
    xGrv = .false.
    if ((cond1.and.cond3).or. &
        (cond2.and.cond3).or. &
        ((.not.cond1).and.cond4).or. &
        ((.not.cond2).and.cond4)) xGrv = .true.
end function xGrv

Subroutine GuessGhostCells(xxx,yyy,zzz,dx,dy,dz,lx,ly,lz,ghost)
implicit none
logical, intent(out) :: ghost
real(rp),intent(in)  :: xxx,yyy,zzz,dx,dy,dz,lx,ly,lz
real(rp):: RR, x_center, y_center, z_center, z_top, z_bottom, eps, temp
real(rp):: sdist1,sdist2,sdistmin,l_y,l_x
logical:: cond1,cond2,cond3,cond4,cond5
real(rp):: LL, AA, BB
real(rp):: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
real(rp):: z_1,z_2,z_3,z_4,abscissa
real(rp):: y_begin, y_end,x_begin,x_end
integer :: ii,jj,q
real(rp):: A
real(rp):: y_12_mid,z_12_mid
real(rp):: Cos_angle,Sin_angle
real(rp):: edge_length_of_box
real(rp):: wall_inc = 0   !wall inclination
integer , dimension(3) :: iperiod
!
Cos_angle = cos(Rotation_angle)
Sin_angle = sin(Rotation_angle)
!
do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo

select case(surface_type)

case('WCyl')
     ghost=.false.
     RR =  r_ibm
     x_center = xc_ibm
     y_center = yc_ibm
     z_bottom = zmin_ibm
     z_top    = zmax_ibm
     eps      = min(dx,dy)
     sdist1 = zzz - z_top

     do jj = -1,1
       do ii = -1,1
         sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                        (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
         if((sdist1.lt.-dz).and.(sdist2.lt.-eps)) then
           ghost = .true.
         elseif((sdist1.gt.dz).and.(sdist2.gt.eps)) then
           ghost = .false.
         endif
       enddo
     enddo
     if(zzz.le.z_bottom) ghost = .true.

case('Sphe')
     ghost=.false.
     RR =  lz/4.0_rp
     x_center = lx/2.0_rp
     y_center = ly/2.0_rp
     z_center = lz/8.0_rp
     eps = 1.5*sqrt(dx**2 + dy**2 + dz**2)

!!!!!!!!!!!!!!!!!!! Added checks to avoid taking sqrt of negative numbers !!!!!!!!!!!!!!!!!!!!!!!!!!!!
     temp = RR**2 - (zzz-z_center)**2 - (yyy-y_center)**2
     if (sign(1.0_rp,temp) .eq. -1.0_rp) temp = 0
     x1 = x_center + sqrt(temp)
     x2 = x_center - sqrt(temp)

     temp = RR**2 - (zzz-z_center)**2 - (xxx-x_center)**2
     if (sign(1.0_rp,temp) .eq. -1.0_rp) temp = 0
     y1 = y_center + sqrt(temp)
     y2 = y_center - sqrt(temp)

     temp = RR**2 - (xxx-x_center)**2 - (yyy-y_center)**2
     if (sign(1.0_rp,temp) .eq. -1.0_rp) temp = 0
     z1 = z_center + sqrt(temp)
     z2 = z_center - sqrt(temp)

     cond1 = (abs(x1-xxx).lt.eps).or.(abs(x2-xxx).lt.eps)
     cond2 = (abs(y1-yyy).lt.eps).or.(abs(y2-yyy).lt.eps)
     cond3 = (abs(z1-zzz).lt.eps).or.(abs(z2-zzz).lt.eps)

     if (cond1.or.cond2.or.cond3) ghost=.true.
     
     LL = (zzz-z_center)**2 + (yyy-y_center)**2 + (xxx-x_center)**2
     if (LL.lt.RR**2) ghost=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

case('Plan')
     ghost=.false.
     if (zzz.lt.zc_ibm) ghost=.true.

case('xGrv')
     ghost=.false.
     cond1 = yyy.lt.r_ibm
     cond2 = yyy.gt.(ly - r_ibm)
     cond3 = zzz.lt.zmax_ibm
     cond4 = zzz.lt.zmin_ibm
     if ((cond1.and.cond3).or. &
         (cond2.and.cond3).or. &
         (.not.cond1.and.cond4).or. &
         (.not.cond2.and.cond4)) ghost=.true.

case('RotB')
     edge_length_of_box =0.5*ly
     ghost= .false.
     A = 0.5*(ly-edge_length_of_box*(Cos_angle+Sin_angle))
     y1 = A + edge_length_of_box*Sin_angle
     z1 = A
     y2 = ly-A
     z2 = A + edge_length_of_box*Sin_angle
     y3 = A + edge_length_of_box*Cos_angle
     z3= lz - A
     y4 = A
     z4 = lz - (A+edge_length_of_box*Sin_angle)
     z_1 = ( (z2-z1)/(y2-y1+small))*(yyy-y1)+z1
     z_2 = ( (z3-z2)/(y3-y2+small))*(yyy-y2)+z2
     z_3 = ( (z4-z3)/(y4-y3+small))*(yyy-y3)+z3
     z_4 = ( (z1-z4)/(y1-y4+small))*(yyy-y4)+z4


     y_12_mid = 0.5*(y1+y2)
     z_12_mid = 0.5*(z1+z2)

   if (Rotation_angle.gt.small) then
    cond1 = zzz.gt.z_1
    cond2 = zzz.lt.z_3
    cond3 = zzz.gt.z_4
    cond4 = zzz.lt.z_2
   else
    cond1 = zzz.gt.z1
    cond2 = yyy.gt.y1
    cond3 = zzz.lt.z3
    cond4 = yyy.lt.y2
   endif
   cond5 = cond1.and.cond2.and.cond3.and.cond4
   if (.not.cond5) ghost=.true.


case('IncW')
   ! ghost=.false.
   ! inside=.true.
   ! length = 0.8*lx
   ! height = (length - xxx)*tan(Rotation_angle)
   
   ! cond1 = (abs(height - zzz).lt.eps)
   ! if (cond1) ghost=.true.
   
   ! cond3 = (zzz.gt.height)
   ! cond4 = (xxx.gt.length)
   ! if (cond3.or.cond4) inside=.false.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ghost=.false.
   edge_length_of_box = sqrt((2.0*lx**2) - (2.0*(0.9375*lx)**2)) 
   A = 0.5*(lx-edge_length_of_box*(Cos_angle+Sin_angle))
   x1 = A + edge_length_of_box*Sin_angle
   z1 = A
   x2 = lx-A
   z2 = A + edge_length_of_box*Sin_angle
   x3 = A + edge_length_of_box*Cos_angle
   z3= lz - A
   x4 = A
   z4 = lz - (A+edge_length_of_box*Sin_angle)
   z_1 = ( (z2-z1)/(x2-x1+small))*(xxx-x1)+z1
   z_2 = ( (z3-z2)/(x3-x2+small))*(xxx-x2)+z2
   z_3 = ( (z4-z3)/(x4-x3+small))*(xxx-x3)+z3
   z_4 = ( (z1-z4)/(x1-x4+small))*(xxx-x4)+z4

   if (Rotation_angle.gt.small) then
    ! cond1 = zzz.gt.z_1
    ! cond2 = zzz.lt.z_3
    cond3 = zzz.gt.z_4
    ! cond4 = zzz.lt.z_2
   else
    ! cond1 = zzz.gt.z1
    ! cond2 = xxx.gt.x1
    ! cond3 = zzz.lt.z3
    ! cond4 = xxx.lt.x2
   endif
   if ((abs(z_4 - zzz)).lt.(dz)) ghost=.true.

end select

end Subroutine GuessGhostCells

Subroutine Solid_Surface(xxx,yyy,zzz,dx,dy,dz,lx,ly,lz,inside)
implicit none
logical, intent(out) :: inside
real(rp),intent(in)  :: xxx,yyy,zzz,dx,dy,dz,lx,ly,lz
real(rp):: LL, RR, x_begin, x_end, y_begin, y_end,x_center, y_center, z_center, z_top, z_bottom
real(rp):: sdist1,sdist2,sdistmin,eps
real(rp):: DD,L_periodic
integer ::j_start,j_end,number_of_blocks ,l
real(rp)::z_1,z_2,z_3,z_4,bb,abscissa,A
real(rp):: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
logical :: cond1,cond2,cond3,cond4,cond5
real(rp)::y_a,z_a,y_b,z_b,y_c,z_c,y_d,z_d,y_e,z_e,y_f,z_f,y_g,z_g,y_h,z_h
real(rp)::yc1,zc1,yc2,zc2,yc3,zc3,yc4,zc4
real(rp)::z_5,z_6,z_7,z_8,z_9,z_10,z_11,z_12,z_13,z_14,z_15,z_16,z_17,z_18,z_19,z_20
logical :: zone1,zone2,zone3,zone4
real(rp)::y_cent_1,y_cent_2,y_cent_3,z_cent_1,z_cent_2,z_cent_3
real(rp):: y_12_mid,z_12_mid
integer :: ii,jj,q
real(rp):: Cos_angle,Sin_angle
real(rp):: edge_length_of_box
real(rp):: wall_inc = 0   !wall inclination
integer , dimension(3) :: iperiod
!
Cos_angle = cos(Rotation_angle)
Sin_angle = sin(Rotation_angle)
!
do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo

select case(surface_type)

case('WCyl')

     inside=.false.
     RR =  r_ibm
     x_center = xc_ibm
     y_center = yc_ibm
     z_bottom = zmin_ibm
     z_top    = zmax_ibm
     eps      = min(dx,dy)
     sdist1 = zzz - z_top

     do jj = -1,1
       do ii = -1,1
         sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                        (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
         if((sdist1.lt.-dz).and.(sdist2.lt.-eps)) then
           inside = .true.
         elseif((sdist1.gt.dz).and.(sdist2.gt.eps)) then
           inside = .false.
         endif
       enddo
     enddo
     if(zzz.le.z_bottom) inside = .true.

case('Sphe')
     inside=.false.
     RR =  lz/4.0_rp
     x_center = lx/2.0_rp
     y_center = ly/2.0_rp
     z_center = lz/8.0_rp
     
     LL = (zzz-z_center)*(zzz-z_center)+(yyy-y_center)*(yyy-y_center)+(xxx-x_center)*(xxx-x_center)
     if (LL.le.RR**2) inside=.true.


case('Plan')
     inside=.false.
     if (zzz.le.zc_ibm) inside=.true.

case('xGrv')
     inside=.false.
     cond1 = yyy.le.r_ibm
     cond2 = yyy.ge.(ly - r_ibm)
     cond3 = zzz.le.zmax_ibm
     cond4 = zzz.le.zmin_ibm
     if ((cond1.and.cond3).or. &
         (cond2.and.cond3).or. &
         (.not.cond1.and.cond4).or. &
         (.not.cond2.and.cond4)) inside=.true.

case('RotB')
   inside= .false.
   edge_length_of_box =0.5*ly
   A = 0.5*(ly-edge_length_of_box*(Cos_angle+Sin_angle))
   y1 = A + edge_length_of_box*Sin_angle
   z1 = A
   y2 = ly-A
   z2 = A + edge_length_of_box*Sin_angle
   y3 = A + edge_length_of_box*Cos_angle
   z3= lz - A
   y4 = A
   z4 = lz - (A+edge_length_of_box*Sin_angle)
   z_1 = ( (z2-z1)/(y2-y1+small))*(yyy-y1)+z1
   z_2 = ( (z3-z2)/(y3-y2+small))*(yyy-y2)+z2
   z_3 = ( (z4-z3)/(y4-y3+small))*(yyy-y3)+z3
   z_4 = ( (z1-z4)/(y1-y4+small))*(yyy-y4)+z4


   y_12_mid = 0.5*(y1+y2)
   z_12_mid = 0.5*(z1+z2)

   if (Rotation_angle.gt.small) then
    cond1 = zzz.gt.z_1
    cond2 = zzz.lt.z_3
    cond3 = zzz.gt.z_4
    cond4 = zzz.lt.z_2
   else
    cond1 = zzz.gt.z1
    cond2 = yyy.gt.y1
    cond3 = zzz.lt.z3
    cond4 = yyy.lt.y2
   endif
   cond5 = cond1.and.cond2.and.cond3.and.cond4
   if (.not.cond5) inside=.true.


case('IncW')
   ! ghost=.false.
   ! inside=.true.
   ! length = 0.8*lx
   ! height = (length - xxx)*tan(Rotation_angle)
   
   ! cond1 = (abs(height - zzz).lt.treshhold)
   ! if (cond1) ghost=.true.
   
   ! cond3 = (zzz.gt.height)
   ! cond4 = (xxx.gt.length)
   ! if (cond3.or.cond4) inside=.false.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   inside= .false.
   edge_length_of_box = sqrt((2.0*lx**2) - (2.0*(0.9375*lx)**2)) 
   A = 0.5*(lx-edge_length_of_box*(Cos_angle+Sin_angle))
   x1 = A + edge_length_of_box*Sin_angle
   z1 = A
   x2 = lx-A
   z2 = A + edge_length_of_box*Sin_angle
   x3 = A + edge_length_of_box*Cos_angle
   z3= lz - A
   x4 = A
   z4 = lz - (A+edge_length_of_box*Sin_angle)
   z_1 = ( (z2-z1)/(x2-x1+small))*(xxx-x1)+z1
   z_2 = ( (z3-z2)/(x3-x2+small))*(xxx-x2)+z2
   z_3 = ( (z4-z3)/(x4-x3+small))*(xxx-x3)+z3
   z_4 = ( (z1-z4)/(x1-x4+small))*(xxx-x4)+z4

   if (Rotation_angle.gt.small) then
    ! cond1 = zzz.gt.z_1
    ! cond2 = zzz.lt.z_3
    cond3 = zzz.gt.z_4
    ! cond4 = zzz.lt.z_2
   else
    ! cond1 = zzz.gt.z1
    ! cond2 = xxx.gt.x1
    ! cond3 = zzz.lt.z3
    ! cond4 = xxx.lt.x2
   endif
   if (.not.cond3) inside=.true.

end select

end subroutine Solid_Surface

Subroutine normal_vectors(n,ng,nh_d,nh_b,Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,dl,dli,zc,dzc)
implicit none
integer,  intent(in)                          :: nh_d,nh_b
integer,  dimension(3), intent(in)            :: n,ng
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in)    :: Level_set
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in)    :: cell_phi_tag
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: nx_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: ny_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: nz_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout) :: nabs_surf
real(rp), dimension(1-nh_d:), intent(in) :: zc,dzc
real(rp), dimension(3), intent(in)  :: dl,dli
integer , dimension(:,:,:), allocatable :: ghost_cell_tag
real(rp), parameter :: eps = 10._rp**(-16)
integer, dimension(3) :: iperiod
integer  :: i,j,k,iperiod1,iperiod2,iperiod3
real(rp) :: nx, ny, nz, n_abs, n_abs_p, vol_frac
real(rp) :: m_av_x,m_av_y,m_av_z,normal_denum
real(rp) :: mx1,mx2,mx3,mx4,mx5,mx6,mx7,mx8
real(rp) :: my1,my2,my3,my4,my5,my6,my7,my8
real(rp) :: mz1,mz2,mz3,mz4,mz5,mz6,mz7,mz8
real(rp) :: dx, dy, dz, dli1, dli2, dli3, dl1, dl2, dl3, lxl, lyl, lzl
integer  :: ip,jp,kp
integer  :: im,jm,km
integer  :: iii,jjj,kkk
real(rp) :: xxx,yyy,zzz
real(rp) :: RR,x_center,y_center,z_center,z_top,z_bottom
integer  :: Level_set_all
logical  :: ghost,inside,ghost_cond
integer  :: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z,q,lvl_set
!
!@cuf attributes(managed) :: Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,zc,dzc,ghost_cell_tag
!@cuf integer :: istat
!
allocate(ghost_cell_tag(1:n(1),1:n(2),1:n(3)))
!@cuf istat = cudaMemAdvise(ghost_cell_tag, size(ghost_cell_tag), cudaMemAdviseSetReadMostly, 0)
ghost_cell_tag(:,:,:) = 0
!
dx=dl(1)
dy=dl(2)
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
lxl = l(1)
lyl = l(2)
lzl = l(3)
!
RR =  r_ibm
x_center = xc_ibm
y_center = yc_ibm
z_bottom = zmin_ibm
z_top    = zmax_ibm
!
do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo
iperiod1 = iperiod(1)
iperiod2 = iperiod(2)
iperiod3 = iperiod(3)
!
dli1=dli(1)
dli2=dli(2)
dli3=dli(3)
!
dl1=dl(1)
dl2=dl(2)
dl3=dl(3)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
! Preparing normal vectors to the surfaces
! Identifying the ghost cells
!
!$acc enter data copyin(ijk_start_x,ijk_start_y,ijk_start_z,dx,dy,dz,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!$acc kernels
do k=1,n3
 do j=1,n2
   do i=1,n1
       xxx   =  (ijk_start_x+i-0.5_rp)*dx
       yyy   =  (ijk_start_y+j-0.5_rp)*dy
       zzz   =  zc(k)
       dz    =  dzc(k)
#if defined(_OPENACC)
       ! ghost = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                     ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
       ghost = xGrv_ghost(RR,x_center,y_center,z_center,z_bottom,z_top, &
                          xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
#else
       call GuessGhostCells(xxx,yyy,zzz,dx,dy,dz,lxl,lyl,lzl,ghost)
#endif
       if (.not.ghost) cycle
       ip = i+1
       jp = j+1
       kp = k+1
       im = i-1
       jm = j-1
       km = k-1
       Level_set_all = Level_set(im,jm,kp) + Level_set(im,j,kp) + Level_set(im,jp,kp) + &
                       Level_set(i ,jm,kp) + Level_set(i ,j,kp) + Level_set(i ,jp,kp) + &
                       Level_set(ip,jm,kp) + Level_set(ip,j,kp) + Level_set(ip,jp,kp) + &
                       Level_set(im,jm,k ) + Level_set(im,j,k ) + Level_set(im,jp,k ) + &
                       Level_set(i ,jm,k ) + Level_set(i ,j,k ) + Level_set(i ,jp,k ) + &
                       Level_set(ip,jm,k ) + Level_set(ip,j,k ) + Level_set(ip,jp,k ) + &
                       Level_set(im,jm,kp) + Level_set(im,j,kp) + Level_set(im,jp,kp) + &
                       Level_set(i ,jm,km) + Level_set(i ,j,km) + Level_set(i ,jp,km) + &
                       Level_set(ip,jm,km) + Level_set(ip,j,km) + Level_set(ip,jp,km)

       if ((Level_set_all.gt.0).and.(Level_set(i,j,k).eq.0)) ghost_cell_tag(i,j,k) = 1

   enddo
 enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,ijk_start_z,dx,dy,dz,lxl,lyl,lzl,RR,x_center,y_center,z_center,z_bottom,z_top,iperiod1,iperiod2,iperiod3)
!
if (myid == 0) print*, '*** Calculating normal vectors ***'
!$acc parallel loop collapse(3) private(mx1,mx2,mx3,mx4,mx5,mx6,mx7,my1,my2,my3,my4,my5,my6,my7,mz1,mz2,mz3,mz4,mz5,mz6,mz7,m_av_x,m_av_y,m_av_z,normal_denum)
do k=1,n3
#if !defined(_OPENACC)
 if (myid == 0) print*, '*** Calculating normal vectors at k = ', k
#endif
 do j=1,n2
  do i=1,n1
   if (ghost_cell_tag(i,j,k).eq.1) then

          !i+1/2 j+1/2 k+1/2
          mx1 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mx2 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i+1,j-1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)))*dli1*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mx3 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i+1,j+1,k-1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)))*dli1*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mx4 = ((cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i+1,j-1,k-1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mx5 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)) - &
                 (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i-1,j+1,k+1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mx6 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)) - &
                 (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i-1,j-1,k+1)))*dli1*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mx7 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)) - &
                 (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i-1,j+1,k-1)))*dli1*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mx8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)) - &
                 (cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli1*0.25_rp
          !
          m_av_x = 0.125_rp*(mx1+mx2+mx3+mx4+mx5+mx6+mx7+mx8)
          ! 
          !i+1/2 j+1/2 k+1/2
          my1 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k+1/2
          my2 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)) - &
                 (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i+1,j-1,k+1)))*dli2*0.25_rp
          !i+1/2 j+1/2 k-1/2
          my3 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i+1,j+1,k-1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)))*dli2*0.25_rp
          !i+1/2 j-1/2 k-1/2
          my4 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)) - &
                 (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i+1,j-1,k-1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k+1/2
          my5 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i-1,j+1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k+1/2
          my6 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)) - &
                 (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i-1,j-1,k+1)))*dli2*0.25_rp
          !i-1/2 j+1/2 k-1/2
          my7 = ((cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i-1,j+1,k-1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)))*dli2*0.25_rp
          !i-1/2 j-1/2 k-1/2
          my8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)) - &
                 (cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli2*0.25_rp
          !
          m_av_y = 0.125_rp*(my1+my2+my3+my4+my5+my6+my7+my8)
          !
          !i+1/2 j+1/2 k+1/2
          mz1 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i+1,j+1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )))*dli3*0.25_rp
          !i+1/2 j-1/2 k+1/2
          mz2 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i+1,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i+1,j-1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )))*dli3*0.25_rp
          !i+1/2 j+1/2 k-1/2
          mz3 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i+1,j+1,k  )) - &
                 (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i+1,j+1,k-1)))*dli3*0.25_rp
          !i+1/2 j-1/2 k-1/2
          mz4 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i+1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i+1,j-1,k  )) - &
                 (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i+1,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i+1,j-1,k-1)))*dli3*0.25_rp
          !i-1/2 j+1/2 k+1/2
          mz5 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i  ,j+1,k+1)+cell_phi_tag(i-1,j+1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )))*dli3*0.25_rp
          !i-1/2 j-1/2 k+1/2
          mz6 = ((cell_phi_tag(i  ,j  ,k+1)+cell_phi_tag(i-1,j  ,k+1)+cell_phi_tag(i  ,j-1,k+1)+cell_phi_tag(i-1,j-1,k+1)) - &
                 (cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )))*dli3*0.25_rp
          !i-1/2 j+1/2 k-1/2
          mz7 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j+1,k  )+cell_phi_tag(i-1,j+1,k  )) - &
                 (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i  ,j+1,k-1)+cell_phi_tag(i-1,j+1,k-1)))*dli3*0.25_rp
          !i-1/2 j-1/2 k-1/2
          mz8 = ((cell_phi_tag(i  ,j  ,k  )+cell_phi_tag(i-1,j  ,k  )+cell_phi_tag(i  ,j-1,k  )+cell_phi_tag(i-1,j-1,k  )) - &
                 (cell_phi_tag(i  ,j  ,k-1)+cell_phi_tag(i-1,j  ,k-1)+cell_phi_tag(i  ,j-1,k-1)+cell_phi_tag(i-1,j-1,k-1)))*dli3*0.25_rp
          !
          m_av_z = 0.125_rp*(mz1+mz2+mz3+mz4+mz5+mz6+mz7+mz8)

      normal_denum = sqrt(m_av_x**2 + m_av_y**2 + m_av_z**2 + eps)

      nx_surf(i,j,k) = m_av_x/normal_denum
      ny_surf(i,j,k) = m_av_y/normal_denum
      nz_surf(i,j,k) = m_av_z/normal_denum

      nabs_surf(i,j,k) = sqrt(nx_surf(i,j,k)**2 + ny_surf(i,j,k)**2 + nz_surf(i,j,k)**2)

     endif

    enddo
  enddo
enddo
!$acc end parallel loop
!
deallocate(ghost_cell_tag)
!
end subroutine normal_vectors
!
Subroutine intersect(n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dl,dzc,zc,zf)
implicit none
integer,  intent(in)                                           :: nh_d,nh_b
integer,  dimension(3), intent(in)                             :: n,ng
real(rp), dimension(3),       intent(in) :: dl
real(rp), dimension(1-nh_d:), intent(in) :: dzc,zc,zf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nx_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: ny_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nz_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nabs_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: x_intersect
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: y_intersect
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: z_intersect
real(rp), parameter :: eps = 10._rp**(-16)
integer, dimension(3) :: iperiod
integer :: i,j,k,ll,iperiod1,iperiod2,iperiod3
real(rp):: nx, ny, nz, n_abs
real(rp):: step
real(rp):: xxx,yyy,zzz,dx,dy,dz,dxx,dyy,dzz,lxl,lyl,lzl
real(rp):: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
logical :: inside
logical :: confirmation
integer :: lmax
real(rp):: distance_ghost_intersect
real(rp):: x_ghost,x_ghost_end,y_ghost,y_ghost_end,z_ghost
real(rp):: RR,x_center,y_center,z_center,z_top,z_bottom
integer :: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z,q
!
!@cuf attributes(managed) :: nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dzc,zc,zf
!
dx=dl(1)
dy=dl(2)
!
n1=n(1)
n2=n(2)
n3=n(3)
!
lxl = l(1)
lyl = l(2)
lzl = l(3)

do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo
iperiod1 = iperiod(1)
iperiod2 = iperiod(2)
iperiod3 = iperiod(3)
!
RR =  r_ibm
x_center = xc_ibm
y_center = yc_ibm
z_center = zc_ibm
z_bottom = zmin_ibm
z_top    = zmax_ibm
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
distance_ghost_intersect= 1000.0_rp
!
!***************************************************************************************
#if defined(_OPENACC)
 if (myid == 0) print*, '*** Calculating intersection points ***'
#endif
!$acc enter data copyin(ijk_start_x,ijk_start_y,dx,dy,dz,lxl,lyl,lzl,step,lmax,confirmation,nx,ny,nz,n_abs,x_ghost,y_ghost,z_ghost,xxx,yyy,zzz,inside,RR,x_center,y_center,z_center,z_bottom,z_top)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
 if (myid == 0) print*, '*** Calculating intersection points at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       dz = dzc(k)
       step = 1.e-5_rp*dz
       lmax = 100*int(2.0_rp*sqrt(3.0_rp)*dz/step)
       confirmation = .false.
       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       x_ghost = (ijk_start_x+i-0.5_rp)*dx
       y_ghost = (ijk_start_y+j-0.5_rp)*dy
       z_ghost = zc(k)
       !$acc loop seq
       do ll=0,lmax
         xxx = x_ghost+ll*(nx/n_abs)*step
         yyy = y_ghost+ll*(ny/n_abs)*step
         zzz = z_ghost+ll*(nz/n_abs)*step
#if defined(_OPENACC)
         inside = xGrv(RR,x_center,y_center,z_center,z_bottom,z_top, &
                       xxx,yyy,zzz,dx,dy,lxl,lyl,dz,(/iperiod1,iperiod2,iperiod3/))
#else
         call Solid_Surface(xxx,yyy,zzz,dl(1),dl(2),dl(3),l(1),l(2),l(3),inside)
#endif
         if (.not.inside) then
           x_intersect(i,j,k) = xxx
           y_intersect(i,j,k) = yyy
           z_intersect(i,j,k) = zzz
           confirmation = .true.
           exit
         endif
      enddo
      !$acc end loop
#if !defined(_OPENACC)
      if (.not.confirmation) then
       print*, '--------------------------------------------------------------------------'
       print*,'Error in detecting intersect point at  i, j, k =' &
              ,i,j,k,'at processor ',myid, 'where the normal vector components are ', &
               nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
      endif

       distance_ghost_intersect =  sqrt((x_intersect(i,j,k)-x_ghost)**2 + &
                                        (y_intersect(i,j,k)-y_ghost)**2 + &
                                        (z_intersect(i,j,k)-z_ghost)**2)

      if (distance_ghost_intersect.gt.sqrt(dx**2 + dy**2 + dz**2)) then
        print*, '--------------------------------------------------------------------------'
        print*, ' Error in detecting intersect point  processor   ', &
        myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
       'where the normal vector components are ', nx,ny,nz,n_abs
        print*, '--------------------------------------------------------------------------'
      endif
#endif
    endif
   enddo
  enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,dz,lxl,lyl,lzl,step,lmax,confirmation,nx,ny,nz,n_abs,x_ghost,y_ghost,z_ghost,xxx,yyy,zzz,inside,RR,x_center,y_center,z_center,z_bottom,z_top)
!
end subroutine intersect

Subroutine mirrorpoints(n,ng,nh_d,nh_b,nx_surf,ny_surf,nz_surf,nabs_surf, &
                        x_intersect,y_intersect,z_intersect, &
                        x_mirror,y_mirror,z_mirror, &
                        deltan, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        dl,dzc,zc)
implicit none
integer,  intent(in)                                           :: nh_d,nh_b
integer,  dimension(3), intent(in)                             :: n,ng
real(rp), dimension(3), intent(in)                             :: dl
real(rp), dimension(1-nh_d:), intent(in)                       :: dzc,zc
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nx_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: ny_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nz_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nabs_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: x_intersect
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: y_intersect
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: z_intersect
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: x_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: y_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: z_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: deltan
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: x_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: y_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: z_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: x_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: y_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: z_IP2
integer :: i,j,k
real(rp):: nx, ny, nz, n_abs
real(rp):: step,step_2
real(rp):: xxx,yyy,zzz
real(rp):: dx,dy,dz,d_1,d_2
real(rp):: cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp):: distance_ghost_intersect
integer :: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z
!@cuf attributes(managed) :: nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect
!@cuf attributes(managed) :: x_mirror,y_mirror,z_mirror,x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2,deltan
!@cuf attributes(managed) :: dzc,zc
!
dx=dl(1)
dy=dl(2)
!
d_1 = d1
d_2 = d2
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
!***************************************************************************************
#if defined(_OPENACC)
if (myid.eq.0) print*, '*** Setting mirror points ***'
#endif
!$acc enter data copyin(ijk_start_x,ijk_start_y,dx,dy,dz,step,step_2,xxx,yyy,zzz,nx,ny,nz,n_abs)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Setting mirror points at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       dz  = dzc(k)
       step   = d_1*sqrt(dx**2 + dy**2 + dz**2)
       step_2 = d_2*sqrt(dx**2 + dy**2 + dz**2)
       !
       xxx = (ijk_start_x+i-0.5_rp)*dx
       yyy = (ijk_start_y+j-0.5_rp)*dy
       zzz = zc(k)
       !
       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       !
       deltan(i,j,k) = sqrt( (x_intersect(i,j,k)-xxx)**2 + &
                             (y_intersect(i,j,k)-yyy)**2 + &
                             (z_intersect(i,j,k)-zzz)**2    )
       !
       if  (deltan(i,j,k).gt.sqrt(dx**2 + dy**2 + dz**2)) then
       print*, '--------------------------------------------------------------------------'
           print*, ' Error: in mirror point detection at cell-center at processor ', &
            myid,' : check IBM.f90 - distance_ghost_intesect is ',deltan(i,j,k), &
            'where the normal vector components are ', nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
       endif
       !
       x_mirror(i,j,k) = x_intersect(i,j,k)+(nx/n_abs)*step
       y_mirror(i,j,k) = y_intersect(i,j,k)+(ny/n_abs)*step
       z_mirror(i,j,k) = z_intersect(i,j,k)+(nz/n_abs)*step
       !
       x_IP1(i,j,k)    = x_mirror(i,j,k)+(nx/n_abs)*step_2
       y_IP1(i,j,k)    = y_mirror(i,j,k)+(ny/n_abs)*step_2
       z_IP1(i,j,k)    = z_mirror(i,j,k)+(nz/n_abs)*step_2
       !
       x_IP2(i,j,k)    = x_IP1(i,j,k)+(nx/n_abs)*step_2
       y_IP2(i,j,k)    = y_IP1(i,j,k)+(ny/n_abs)*step_2
       z_IP2(i,j,k)    = z_IP1(i,j,k)+(nz/n_abs)*step_2
    endif
   enddo
  enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,dz,step,step_2,xxx,yyy,zzz,nx,ny,nz,n_abs)
!
end subroutine mirrorpoints

Subroutine mirrorpoints_ijk(n,ng,nh_d,nh_b,nabs_surf,x_mirror,y_mirror,z_mirror, &
                            x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                            i_mirror,j_mirror,k_mirror, &
                            i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                            dl,zf)

implicit none
integer,  intent(in)                                           :: nh_d,nh_b
integer,  dimension(3), intent(in)                             :: n,ng
real(rp), dimension(1-nh_d:), intent(in)                       :: zf
real(rp), dimension(3)      , intent(in)                       :: dl
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: nabs_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: x_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: y_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: z_mirror
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: x_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: y_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: z_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: x_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: y_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(in) :: z_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: i_mirror
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: j_mirror
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: k_mirror
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: i_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: j_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: k_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: i_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: j_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b), intent(inout):: k_IP2
integer :: i,j,k,ll,m,nn
real(rp):: nx, ny, nz, n_abs
real(rp):: xxx,yyy,zzz
real(rp):: dx,dy
real(rp):: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp):: distance_ghost_intesect
integer :: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z
!@cuf attributes(managed) :: nabs_surf,x_mirror,y_mirror,z_mirror
!@cuf attributes(managed) :: x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2
!@cuf attributes(managed) :: i_mirror,j_mirror,k_mirror,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
!@cuf attributes(managed) :: zf
!
dx=dl(1)
dy=dl(2)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
n1=n(1)
n2=n(2)
n3=n(3)
!
!***************************************************************************************
#if defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining mirror point indices ***'
#endif
!$acc enter data copyin(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining mirror point indices at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       xxx = x_mirror(i,j,k)
       yyy = y_mirror(i,j,k)
       zzz = z_mirror(i,j,k)
       !$acc loop seq
       do ll = 1,n3
         !$acc loop seq
         do m = 1-nh_b,n2+nh_b
           !$acc loop seq
           do nn = 1-nh_b,n1+nh_b
              cell_start_x = (ijk_start_x+nn-1.0_rp)*dx
              cell_end_x   = (ijk_start_x+nn-0.0_rp)*dx
              cell_start_y = (ijk_start_y+m -1.0_rp)*dy
              cell_end_y   = (ijk_start_y+m -0.0_rp)*dy
              cell_start_z = zf(ll-1)
              cell_end_z   = zf(ll)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_mirror(i,j,k) = nn
                  j_mirror(i,j,k) = m
                  k_mirror(i,j,k) = ll
              exit
             endif
          enddo
          !$acc end loop
        enddo
        !$acc end loop
       enddo
       !$acc end loop
       !
#if !defined(_OPENACC)
       if ((i_mirror(i,j,k).eq.i).and.(j_mirror(i,j,k).eq.j).and.(k_mirror(i,j,k).eq.k)) then
             print*, 'Error: Ghost and mirror point are the same'
       endif
       !
       if ((i_mirror(i,j,k).eq.-1000).or.(j_mirror(i,j,k).eq.-1000).or.(k_mirror(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for mirror point at center i = ', &
                 i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif
#endif
    endif
  enddo
 enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!
#if defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining interpolation point 1 indices ***'
#endif
!$acc enter data copyin(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining interpolation point 1 indices at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       xxx = x_IP1(i,j,k)
       yyy = y_IP1(i,j,k)
       zzz = z_IP1(i,j,k)
       !$acc loop seq
       do ll = 1,n3
         !$acc loop seq
         do m = 1-nh_b,n2+nh_b
           !$acc loop seq
           do nn = 1-nh_b,n1+nh_b
              cell_start_x = (ijk_start_x+nn-1.0_rp)*dx
              cell_end_x   = (ijk_start_x+nn-0.0_rp)*dx
              cell_start_y = (ijk_start_y+m -1.0_rp)*dy
              cell_end_y   = (ijk_start_y+m -0.0_rp)*dy
              cell_start_z = zf(ll-1)
              cell_end_z   = zf(ll)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP1(i,j,k) = nn
                  j_IP1(i,j,k) = m
                  k_IP1(i,j,k) = ll
              exit
             endif
          enddo
          !$acc end loop
        enddo
        !$acc end loop
       enddo
       !$acc end loop
       !
#if !defined(_OPENACC)
       if ((i_IP1(i,j,k).eq.i).and.(j_IP1(i,j,k).eq.j).and.(k_IP1(i,j,k).eq.k)) then
             print*, 'Error: Ghost and IP1 are the same'
       endif
       !
       if ((i_IP1(i,j,k).eq.-1000).or.(j_IP1(i,j,k).eq.-1000).or.(k_IP1(i,j,k).eq.-1000)) then
        print*, '--------------------------------------------------------------------------'
          print*,'Error: no grid point detected for IP1 at center i = ', &
                 i, ' j= ',j,' k= ',k, 'at processor ',myid
        print*, '--------------------------------------------------------------------------'
      endif
#endif
    endif
  enddo
 enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!
#if defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining interpolation point 2 indices ***'
#endif
!$acc enter data copyin(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!$acc kernels
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Determining interpolation point 2 indices at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       xxx = x_IP2(i,j,k)
       yyy = y_IP2(i,j,k)
       zzz = z_IP2(i,j,k)
       !$acc loop seq
       do ll = 1,n3
         !$acc loop seq
         do m = 1-nh_b,n2+nh_b
           !$acc loop seq
           do nn = 1-nh_b,n1+nh_b
              cell_start_x = (ijk_start_x+nn-1.0_rp)*dx
              cell_end_x   = (ijk_start_x+nn-0.0_rp)*dx
              cell_start_y = (ijk_start_y+m -1.0_rp)*dy
              cell_end_y   = (ijk_start_y+m -0.0_rp)*dy
              cell_start_z = zf(ll-1)
              cell_end_z   = zf(ll)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP2(i,j,k) = nn
                  j_IP2(i,j,k) = m
                  k_IP2(i,j,k) = ll
              exit
             endif
          enddo
          !$acc end loop
        enddo
        !$acc end loop
       enddo
       !$acc end loop
       !
#if !defined(_OPENACC)
       if ((i_IP2(i,j,k).eq.i).and.(j_IP2(i,j,k).eq.j).and.(k_IP2(i,j,k).eq.k)) then
             print*, 'Error: Ghost and IP2 point are the same'
       endif
       !
       if ((i_IP2(i,j,k).eq.-1000).or.(j_IP2(i,j,k).eq.-1000).or.(k_IP2(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for IP2 at center i = ', &
                 i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif
#endif
    endif
  enddo
 enddo
enddo
!$acc end kernels
!$acc exit data copyout(ijk_start_x,ijk_start_y,dx,dy,xxx,yyy,zzz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z)
!
end Subroutine mirrorpoints_ijk
!
Subroutine InterpolationWeights(n,ng,nh_d,nh_b,nabs_surf,Level_set, &
                                x_IP1,y_IP1,z_IP1, &
                                x_IP2,y_IP2,z_IP2, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2, &
                                dl,dzc,zc)
implicit none
integer,  intent(in)                                          :: nh_d,nh_b
integer,  dimension(3), intent(in)                            :: n,ng
real(rp), dimension(3), intent(in)                            :: dl
real(rp), dimension(1-nh_d:), intent(in)                      :: zc,dzc
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: Level_set
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: nabs_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: x_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: y_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: z_IP1
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: x_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: y_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: z_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: i_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: j_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: k_IP1
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: i_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: j_IP2
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in   ) :: k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7), intent(inout) :: WP1,WP2
real(rp):: x1(7),y1(7),z1(7),x2(7),y2(7),z2(7)
real(rp):: h1(7),h2(7)
real(rp):: xx1,xx2,yy1,yy2,zz1,zz2,dx,dy
real(rp):: contribution1(7),contribution2(7)
integer :: i_p_1(7),j_p_1(7),k_p_1(7),i_p_2(7),j_p_2(7),k_p_2(7)
integer :: ii,i,j,k,ii1,ii2,jj1,jj2,kk1,kk2
integer :: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z
logical :: cond11,cond12,cond21,cond22
!@cuf attributes(managed) :: nabs_surf,Level_set
!@cuf attributes(managed) :: x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2
!@cuf attributes(managed) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
!@cuf attributes(managed) :: WP1,WP2,dzc,zc
!
dx=dl(1)
dy=dl(2)
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
#if defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating Interpolation Weights ***'
#endif
!$acc enter data copyin(solid_height_ratio,n1,n2,n3,ijk_start_x,ijk_start_y,dx,dy)
!$acc enter data create(xx1,xx2,yy1,yy2,zz1,zz2,cond11,cond12,cond21,cond22,ii1,ii2,jj1,jj2,kk1,kk2)
!$acc enter data create(x1(:),y1(:),z1(:),x2(:),y2(:),z2(:),h1(:),h2(:),contribution1(:),contribution2(:),i_p_1(:),j_p_1(:),k_p_1(:),i_p_2(:),j_p_2(:),k_p_2(:))
!$acc kernels
do k = 1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating Interpolation Weights at k = ', k
#endif
 do j = 1,n2
   do i = 1,n1
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
     !Initialization
     xx1               = 0.0_rp
     xx2               = 0.0_rp
     yy1               = 0.0_rp
     yy2               = 0.0_rp
     zz1               = 0.0_rp
     zz2               = 0.0_rp
     cond11            = .false.
     cond21            = .false.
     cond12            = .false.
     cond22            = .false.
     contribution1(:)  =  0.0_rp
     contribution2(:)  =  0.0_rp
     x1(:)             =  0.0_rp
     y1(:)             =  0.0_rp
     z1(:)             =  0.0_rp
     x2(:)             =  0.0_rp
     y2(:)             =  0.0_rp
     z2(:)             =  0.0_rp
     h1(:)             =  0.0_rp
     h2(:)             =  0.0_rp
     i_p_1(:)          =  0
     j_p_1(:)          =  0
     k_p_1(:)          =  0
     i_p_2(:)          =  0
     j_p_2(:)          =  0
     k_p_2(:)          =  0
     !
     if (nabs_surf(i,j,k).gt.small) then
         xx1  =  x_IP1(i,j,k)
         yy1  =  y_IP1(i,j,k)
         zz1  =  z_IP1(i,j,k)
         ii1  =  i_IP1(i,j,k)
         jj1  =  j_IP1(i,j,k)
         kk1  =  k_IP1(i,j,k)
         
         xx2  =  x_IP2(i,j,k)
         yy2  =  y_IP2(i,j,k)
         zz2  =  z_IP2(i,j,k)
         ii2  =  i_IP2(i,j,k)
         jj2  =  j_IP2(i,j,k)
         kk2  =  k_IP2(i,j,k)
         
         i_p_1(1) = ii1
         j_p_1(1) = jj1+1
         k_p_1(1) = kk1
         
         i_p_2(1) = ii2
         j_p_2(1) = jj2+1
         k_p_2(1) = kk2
         
         i_p_1(2) = ii1
         j_p_1(2) = jj1
         k_p_1(2) = kk1+1
         
         i_p_2(2) = ii2
         j_p_2(2) = jj2
         k_p_2(2) = kk2+1
         
         i_p_1(3) = ii1
         j_p_1(3) = jj1-1
         k_p_1(3) = kk1
         
         i_p_2(3) = ii2
         j_p_2(3) = jj2-1
         k_p_2(3) = kk2
         
         i_p_1(4) = ii1
         j_p_1(4) = jj1
         k_p_1(4) = kk1-1
         
         i_p_2(4) = ii2
         j_p_2(4) = jj2
         k_p_2(4) = kk2-1
         
         i_p_1(5) = ii1
         j_p_1(5) = jj1
         k_p_1(5) = kk1
         
         i_p_2(5) = ii2
         j_p_2(5) = jj2
         k_p_2(5) = kk2
         
         i_p_1(6) = ii1-1
         j_p_1(6) = jj1
         k_p_1(6) = kk1
         
         i_p_2(6) = ii2-1
         j_p_2(6) = jj2
         k_p_2(6) = kk2
         
         i_p_1(7) = ii1+1
         j_p_1(7) = jj1
         k_p_1(7) = kk1
         
         i_p_2(7) = ii2+1
         j_p_2(7) = jj2
         k_p_2(7) = kk2

        do ii = 1,7
           cond11 = (Level_set(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.1)
           cond12 = (nabs_surf(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.0.0_rp)
           cond21 = (Level_set(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.1)
           cond22 = (nabs_surf(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.0.0_rp)

           if (cond11.and.cond12) contribution1(ii) = 1.0_rp
           if (cond21.and.cond22) contribution1(ii) = 1.0_rp
           ! if (cond11) contribution1(ii) = 1.0_rp
           ! if (cond21) contribution2(ii) = 1.0_rp

          x1(ii) = (ijk_start_x+i_p_1(ii)-0.5_rp)*dx
          y1(ii) = (ijk_start_y+j_p_1(ii)-0.5_rp)*dy
          z1(ii) =  zc(k_p_1(ii))
          x2(ii) = (ijk_start_x+i_p_2(ii)-0.5_rp)*dx
          y2(ii) = (ijk_start_y+j_p_2(ii)-0.5_rp)*dy
          z2(ii) =  zc(k_p_2(ii))
        enddo

        do ii = 1,7
         h1(ii) = sqrt( (x1(ii)-xx1)**2 + (y1(ii)-yy1)**2 + (z1(ii)-zz1)**2)
         h2(ii) = sqrt( (x2(ii)-xx2)**2 + (y2(ii)-yy2)**2 + (z2(ii)-zz2)**2)
        enddo

        do ii = 1,7
         WP1(i,j,k,ii)  = (1.0_rp/(h1(ii)*h1(ii)))*contribution1(ii)
         WP2(i,j,k,ii)  = (1.0_rp/(h2(ii)*h2(ii)))*contribution2(ii)
        enddo

        !-------- Exceptional cases ---------------------
        do ii = 1,7
         if ((h1(ii).lt.(1.0e-8_rp)).and.(contribution1(ii).eq.1.0_rp)) then
             WP1(i,j,k,1)  = 0.0_rp
             WP1(i,j,k,2)  = 0.0_rp
             WP1(i,j,k,3)  = 0.0_rp
             WP1(i,j,k,4)  = 0.0_rp
             WP1(i,j,k,5)  = 0.0_rp
             WP1(i,j,k,6)  = 0.0_rp
             WP1(i,j,k,7)  = 0.0_rp
             WP1(i,j,k,ii) = 1.0_rp
         endif
         if ((h2(ii).lt.(1.0e-8_rp)).and.(contribution2(ii).eq.1.0_rp)) then
             WP2(i,j,k,1)  = 0.0_rp
             WP2(i,j,k,2)  = 0.0_rp
             WP2(i,j,k,3)  = 0.0_rp
             WP2(i,j,k,4)  = 0.0_rp
             WP2(i,j,k,5)  = 0.0_rp
             WP2(i,j,k,6)  = 0.0_rp
             WP2(i,j,k,7)  = 0.0_rp
             WP2(i,j,k,ii) = 1.0_rp
         endif
       enddo
       !-------------------------------------------
     endif

   enddo
  enddo
enddo
!$acc end kernels
!$acc exit data copyout(solid_height_ratio,n1,n2,n3,ijk_start_x,ijk_start_y,dx,dy)
!$acc exit data copyout(xx1,xx2,yy1,yy2,zz1,zz2,cond11,cond12,cond21,cond22,ii1,ii2,jj1,jj2,kk1,kk2)
!$acc exit data copyout(x1(:),y1(:),z1(:),x2(:),y2(:),z2(:),h1(:),h2(:),contribution1(:),contribution2(:),i_p_1(:),j_p_1(:),k_p_1(:),i_p_2(:),j_p_2(:),k_p_2(:))
!
end subroutine InterpolationWeights

subroutine interpolation_mirror(n,nh_v,nh_b,vof,iii,jjj,kkk, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2,BB)
implicit none
!$acc routine(interpolation_mirror) seq
integer,  dimension(3), intent(in)                                 :: n
integer,  intent(in)                                               :: nh_v,nh_b
integer,  intent(in)                                               :: iii,jjj,kkk
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),    intent(in) :: vof
integer,  dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),    intent(in) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7),intent(in) :: WP1,WP2
real(rp), intent(out) :: BB
real(rp), parameter   :: eps = 10._rp**(-16)
real(rp) :: q_1,WW1_1,WW2_1,WW3_1,WW4_1,WW5_1,WW6_1,WW7_1,B1_1,B2_1,B3_1,B4_1,B5_1,B6_1,B7_1
real(rp) :: q_2,WW1_2,WW2_2,WW3_2,WW4_2,WW5_2,WW6_2,WW7_2,B1_2,B2_2,B3_2,B4_2,B5_2,B6_2,B7_2
real(rp) :: BB1,BB2
integer  :: ii1,jj1,kk1,ii2,jj2,kk2
!@cuf attributes(managed) :: vof,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2
!
ii1   = i_IP1(iii,jjj,kkk) ! Coordinates for cell of interpolation point 1 (IP1)
jj1   = j_IP1(iii,jjj,kkk)
kk1   = k_IP1(iii,jjj,kkk)
ii2   = i_IP2(iii,jjj,kkk) ! Coordinates for cell of interpolation point 2  (IP2)
jj2   = j_IP2(iii,jjj,kkk)
kk2   = k_IP2(iii,jjj,kkk)
!
WW1_1 = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
WW2_1 = WP1(iii,jjj,kkk,2)
WW3_1 = WP1(iii,jjj,kkk,3)
WW4_1 = WP1(iii,jjj,kkk,4)
WW5_1 = WP1(iii,jjj,kkk,5)
WW6_1 = WP1(iii,jjj,kkk,6)
WW7_1 = WP1(iii,jjj,kkk,7)
WW1_2 = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP1
WW2_2 = WP2(iii,jjj,kkk,2)
WW3_2 = WP2(iii,jjj,kkk,3)
WW4_2 = WP2(iii,jjj,kkk,4)
WW5_2 = WP2(iii,jjj,kkk,5)
WW6_2 = WP2(iii,jjj,kkk,6)
WW7_2 = WP2(iii,jjj,kkk,7)
!
B1_1  = vof(ii1,jj1+1,kk1)  ! Values for interpolating to IP1
B2_1  = vof(ii1,jj1,kk1+1)
B3_1  = vof(ii1,jj1-1,kk1)
B4_1  = vof(ii1,jj1,kk1-1)
B5_1  = vof(ii1,jj1,kk1)
B6_1  = vof(ii1-1,jj1,kk1)
B7_1  = vof(ii1+1,jj1,kk1)
B1_2  = vof(ii2,jj2+1,kk2)  ! Values for interpolating to IP1
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

BB = 2.0_rp*BB1-BB2 ! Value at mirror point

end subroutine interpolation_mirror

Subroutine interpolation_dphi(n,nh_v,nh_b,vof,iii,jjj,kkk, &
                                        i_IP1,j_IP1,k_IP1, &
                                        i_IP2,j_IP2,k_IP2, &
                                                  WP1,WP2, &
                                     dphidx,dphidy,dphidz, &
                                     dx,dy,dz,dxi,dyi,dzi   )
implicit none
!$acc routine(interpolation_dphi) seq
integer,  dimension(3), intent(in)                                 :: n
integer,  intent(in)                                               :: nh_v,nh_b
integer , intent(in)                                               :: iii,jjj,kkk
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),     intent(in)  :: vof
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in)  :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7), intent(in)  :: WP1,WP2
real(rp),                                                                    intent(in)  :: dx,dy,dz,dxi,dyi,dzi
real(rp)                                                                   , intent(out) :: dphidx,dphidy,dphidz
real(rp), parameter :: eps = 10._rp**(-16)
integer  :: ii1,jj1,kk1,ii2,jj2,kk2
real(rp) :: q_1,WW1_1,WW2_1,WW3_1,WW4_1,WW5_1,WW6_1,WW7_1
real(rp) :: q_2,WW1_2,WW2_2,WW3_2,WW4_2,WW5_2,WW6_2,WW7_2
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
!@cuf attributes(managed) :: vof,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2
!
ii1 = i_IP1(iii,jjj,kkk)
jj1 = j_IP1(iii,jjj,kkk)
kk1 = k_IP1(iii,jjj,kkk)
ii2 = i_IP2(iii,jjj,kkk)
jj2 = j_IP2(iii,jjj,kkk)
kk2 = k_IP2(iii,jjj,kkk)

WW1_1 = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
WW2_1 = WP1(iii,jjj,kkk,2)
WW3_1 = WP1(iii,jjj,kkk,3)
WW4_1 = WP1(iii,jjj,kkk,4)
WW5_1 = WP1(iii,jjj,kkk,5)
WW6_1 = WP1(iii,jjj,kkk,6)
WW7_1 = WP1(iii,jjj,kkk,7)
WW1_2 = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP2
WW2_2 = WP2(iii,jjj,kkk,2)
WW3_2 = WP2(iii,jjj,kkk,3)
WW4_2 = WP2(iii,jjj,kkk,4)
WW5_2 = WP2(iii,jjj,kkk,5)
WW6_2 = WP2(iii,jjj,kkk,6)
WW7_2 = WP2(iii,jjj,kkk,7)

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

q_1 = WW1_1 + WW2_1 + WW3_1 + WW4_1 + WW5_1 + WW6_1 + WW7_1
q_2 = WW1_2 + WW2_2 + WW3_2 + WW4_2 + WW5_2 + WW6_2 + WW7_2

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

end subroutine interpolation_dphi

subroutine Penalization_center(n,ng,nh_d,nh_b,uu,vv,ww,cell_phi_tag,nabs_surf)
implicit none
integer,  intent(in)                                             :: nh_d,nh_b
integer,  dimension(3), intent(in)                               :: n,ng
real(rp),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),intent(in)    :: cell_phi_tag,nabs_surf
real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1),intent(inout) :: uu,vv,ww
integer::i,j,k,n1,n2,n3
!@cuf attributes(managed) :: uu,vv,ww,cell_phi_tag,nabs_surf
!
n1=n(1)
n2=n(2)
n3=n(3)
!
!$acc parallel loop collapse(3)
do k=1,n3
 do j=1,n2
  do i=1,n1
    ! if (slip_length.gt.small) then
     ! if (nabs_surf(i,j,k).lt.small) then
      ! uu(i,j,k)=cell_phi_tag(i,j,k)*uu(i,j,k)
      ! vv(i,j,k)=cell_phi_tag(i,j,k)*(vv(i,j,k))
     ! endif
     ! if (nabs_surf(i,j,k).lt.small) then
      ! ww(i,j,k)=cell_phi_tag(i,j,k)*ww(i,j,k)
     ! endif
    ! else
      uu(i,j,k)= cell_phi_tag(i,j,k)*uu(i,j,k)
      vv(i,j,k)= cell_phi_tag(i,j,k)*vv(i,j,k)
      ww(i,j,k)= cell_phi_tag(i,j,k)*ww(i,j,k)
    ! endif
   enddo
  enddo
enddo
!$acc end parallel loop
end subroutine Penalization_center


subroutine Penalization_face(n,ng,nh_d,nh_b,uu,vv,ww,cell_u_tag,cell_v_tag,cell_w_tag)
implicit none
integer,  intent(in)                                             :: nh_d,nh_b
integer,  dimension(3), intent(in)                               :: n,ng
real(rp),dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),intent(in)::cell_u_tag,cell_v_tag,cell_w_tag
real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1),intent(inout)::uu,vv,ww
integer::i,j,k,n1,n2,n3
!@cuf attributes(managed) :: uu,vv,ww,cell_u_tag,cell_v_tag,cell_w_tag
!
n1=n(1)
n2=n(2)
n3=n(3)
!
!$acc parallel loop collapse(3)
do k=0,n3
 do j=0,n2
  do i=0,n1
    ! if (slip_length.gt.small) then
       ! if (.not. ((cell_u_tag(i,j,k).lt.1).and.(cell_u_tag(i,j,k).gt.small))) then
          ! uu(i,j,k)=cell_u_tag(i,j,k)*uu(i,j,k)
       ! endif
       ! if (.not. ((cell_v_tag(i,j,k).lt.1).and.(cell_v_tag(i,j,k).gt.small))) then
          ! vv(i,j,k)=cell_v_tag(i,j,k)*(vv(i,j,k))
       ! endif
       ! if (.not. ((cell_w_tag(i,j,k).lt.1).and.(cell_w_tag(i,j,k).gt.small))) then
           ! ww(i,j,k)=cell_w_tag(i,j,k)*ww(i,j,k)
       ! endif
    ! else
           uu(i,j,k)= cell_u_tag(i,j,k)*uu(i,j,k)
           vv(i,j,k)= cell_v_tag(i,j,k)*vv(i,j,k)
           ww(i,j,k)= cell_w_tag(i,j,k)*ww(i,j,k)

    ! endif
   enddo
  enddo
enddo
!$acc end parallel loop
end subroutine Penalization_face
!
subroutine Wetting_radius(n,ng,nh_v,nh_b, &
                          time,nabs_surf,vof,deltan, &
                          i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                          WP1,WP2, &
                          dzc,dl,dims)

implicit none
integer , dimension(3), intent(in)                                                      :: n,ng
integer , intent(in)                                                                    :: nh_v,nh_b
real(rp), intent(in)                                                                    :: time
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in) :: nabs_surf
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in) :: deltan
real(rp), dimension(1-nh_v:n(1)+nh_v,1-nh_v:n(2)+nh_v,1-nh_v:n(3)+nh_v),     intent(in) :: vof
integer , dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b),     intent(in) :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
real(rp), dimension(1-nh_b:n(1)+nh_b,1-nh_b:n(2)+nh_b,1-nh_b:n(3)+nh_b,1:7), intent(in) :: WP1,WP2
real(rp), dimension(0:), intent(in)                                                     :: dzc
real(rp), dimension(3), intent(in)                                                      :: dl
integer , dimension(3), intent(in)                                                      :: dims
real(rp):: starting_point_j,end_point_j,starting_point_k,end_point_k,phi_wall
real(rp):: time_stari_visc,time_star_inert
real(rp):: Wetting_rad, phi_m,delta_y,delta_z
real(rp), dimension(1:dims(2)) :: j_min,j_max,k_min,k_max
real(rp):: j_min_all,j_max_all,dn_1,dn_2,dx,dy,alpha,eps
integer :: n1,n2,n3,nv,nb,k,j,i,ijk_start_x,ijk_start_y,ijk_start_z
integer :: dims1, dims2, min_processor,max_processor
integer :: error
logical :: confirmation
!@cuf attributes(managed) :: nabs_surf,vof,deltan,vof,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2,dzc
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
nv = nh_v
nb = nh_b
!
dims1 = dims(1)
dims2 = dims(2)
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
dx=dl(1)
dy=dl(2)
!
eps = small
!
j_min(:) =  10000._rp
j_max(:) = -10000._rp
k_min(:) =  10000._rp
k_max(:) = -10000._rp
Wetting_rad = 0._rp
j_min_all =  10000._rp
j_max_all = -10000._rp
min_processor = 0
max_processor = 0
confirmation = .false.
delta_y = 0._rp
delta_z = 0._rp

  i=int(n1/2)
  do j=1,n2
    do k=1,n3
       if (nabs_surf(i,j,k).gt.eps) then
          call interpolation_mirror((/n1,n2,n3/),nv,nb,vof,i,j,k, &
                                    i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                                    WP1,WP2,phi_m)



          if (phi_m.gt.1.0e-8) then
           j_min(ijk_start_y+1) = real((j+ijk_start_y),rp)
           k_min(ijk_start_y+1) = real(k,rp)
           confirmation = .true.
           exit
          endif
       endif
    enddo
    if (confirmation) exit
  enddo

  confirmation = .false.
  do j=n(2),1,-1
   do k=n(3),1,-1
      if  (nabs_surf(i,j,k).gt.small) then
          call interpolation_mirror((/n1,n2,n3/),nv,nb,vof,i,j,k, &
                                    i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                                    WP1,WP2,phi_m)
          if ((phi_m).gt.1.0e-8) then
           j_max(ijk_start_y+1) = real((j+ijk_start_y),rp)
           k_max(ijk_start_y+1) = real(k,rp)
           confirmation = .true.
           exit
          endif
      endif
   enddo
   if (confirmation) exit
  enddo

  call mpi_allreduce(MPI_IN_PLACE,j_min(1),dims(2),MPI_REAL_RP,mpi_min,MPI_COMM_WORLD,error)
  call mpi_allreduce(MPI_IN_PLACE,j_max(1),dims(2),MPI_REAL_RP,mpi_max,MPI_COMM_WORLD,error)
  call mpi_allreduce(MPI_IN_PLACE,k_min(1),dims(2),MPI_REAL_RP,mpi_min,MPI_COMM_WORLD,error)
  call mpi_allreduce(MPI_IN_PLACE,k_max(1),dims(2),MPI_REAL_RP,mpi_max,MPI_COMM_WORLD,error)


j_min_all =  10000.0_rp
j_max_all = -10000.0_rp
  do i = 1, dims(2)
    if (j_min(i).lt.j_min_all) then
        j_min_all = j_min(i)
        min_processor = i
    endif
    if (j_max(i).gt.j_max_all) then
        j_max_all = j_max(i)
        max_processor = i
    endif
enddo

  delta_y = (j_max_all-j_min_all)*dy
  delta_z = (k_max(max_processor)-k_min(min_processor))*dzc(k)
  Wetting_rad = 0.5_rp*sqrt(delta_y**2 + delta_z**2)
  if (myid.eq.0) then
    open(1329,file='data/spreading_radius.txt',position='append')
    write(1329,'(2E16.8)' ) time, Wetting_rad
    close(1329)
  endif

end subroutine Wetting_radius
!
#endif
end module mod_IBM