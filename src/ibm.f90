module mod_IBM
#if defined(_USE_IBM)
use mpi
use mod_param
use mod_common_mpi, only: myid,ierr,ijk_start
use mod_types
!@cuf use cudafor
!@cuf use mod_common_mpi, only: mydev
implicit none
private
public IBM_Mask,normal_vectors,intersect,mirrorpoints, &
       Penalization_center,Penalization_face,interpolation_2D_velocity,interpolation_dphi, &
       mirrorpoints_ijk, interpolation_mirror,InterpolationWeights,Wetting_radius
contains
!

subroutine IBM_Mask(cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag,Level_set,zc,zf,dl,dzc)
implicit none
real(rp), intent(in), dimension(0:)        :: zc,zf,dzc
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: cell_u_tag,cell_v_tag,cell_w_tag,cell_phi_tag
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: Level_set
real(rp), intent(in), dimension(3)       :: dl
integer i,j,k,l,nn,m,number_of_divisions
real(rp):: xxx,yyy,zzz,dx,dy,dz,dxx,dyy,dzz,lxl,lyl
real(rp):: cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp):: counter
real(rp):: RR,x_center,y_center,z_center,z_top,z_bottom
integer, dimension(3) :: iperiod
logical:: inside,ghost
integer:: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z,q
#if defined(_OPENACC)
!@cuf integer :: istat
!@cuf attributes(managed) :: cell_phi_tag
#endif
!
dx=dl(1)
dy=dl(2)
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
lxl = lx
lyl = ly
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
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
! Wall Geometry
#if defined(_OPENACC)
number_of_divisions = 100
#else
number_of_divisions = 50
#endif
!$acc enter data create(xxx,yyy,zzz,dxx,dyy,dzz,dz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
!$acc enter data copyin(iperiod) async(1)
#if defined(_OPENACC)
if (myid == 0) print*, '*** Calculating volume fractions ***'
#endif
!$acc parallel loop gang collapse(3) default(present) private(xxx,yyy,zzz,dxx,dyy,dzz,dz,cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z,ghost,inside) async(1)
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
  if (myid.eq.0) print*, '*** Calculating volume fractions at k = ', k
#endif
  do j=1,n2
    do i=1,n1
      xxx =  (i+ijk_start_x)*dx-0.5_rp*dx
      yyy =  (j+ijk_start_y)*dy-0.5_rp*dy
      zzz =  zc(k)
      dz  = dzc(k)
#if defined(_OPENACC)
      ghost = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                    xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
#else
      call GuessGhostCells(xxx,yyy,zzz,ghost)
#endif
      if (ghost) then
         cell_phi_tag(i,j,k) = 1.0_rp
         cell_u_tag(i,j,k)   = 1.0_rp
         cell_v_tag(i,j,k)   = 1.0_rp
         cell_w_tag(i,j,k)   = 1.0_rp
         Level_set(i,j,k)    = 1
      endif

! Cell Center
      cell_start_x = (i+ijk_start_x-1)*dx
      cell_end_x   = (i+ijk_start_x  )*dx

      cell_start_y = (j+ijk_start_y-1)*dy
      cell_end_y   = (j+ijk_start_y  )*dy

      cell_start_z = zf(k-1)
      cell_end_z   = zf(k)

      dxx = (cell_end_x-cell_start_x)/number_of_divisions
      dyy = (cell_end_y-cell_start_y)/number_of_divisions
      dzz = (cell_end_z-cell_start_z)/number_of_divisions

      counter = 0._rp
     !$acc loop seq
      do nn= 1,number_of_divisions
          zzz = cell_start_z+(nn-1)*dzz
       !$acc loop seq
        do m = 1,number_of_divisions
            yyy = cell_start_y + (m-1)*dyy
           !$acc loop seq
            do l = 1,number_of_divisions
              xxx = cell_start_x + (l-1)*dxx
#if defined(_OPENACC)
              inside = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                             xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
#else
              call Solid_Surface(xxx,yyy,zzz,inside)
#endif
              if (inside) counter = counter +1
            enddo
        enddo
      enddo
      cell_phi_tag(i,j,k) = counter/(1.0_rp*number_of_divisions**3)

! u cells
      ! cell_start_x = (i+ijk_start_x)*dx-0.5_rp*dx
      ! cell_end_x   = (i+ijk_start_x)*dx-0.5_rp*dx

      ! cell_start_y = (j+ijk_start_y-1)*dy
      ! cell_end_y   = (j+ijk_start_y  )*dy

      ! cell_start_z = zf(k-1)
      ! cell_end_z   = zf(k)

      ! dxx = (cell_end_x-cell_start_x)/number_of_divisions
      ! dyy = (cell_end_y-cell_start_y)/number_of_divisions
      ! dzz = (cell_end_z-cell_start_z)/number_of_divisions

      ! counter = 0._rp
!      !$acc loop seq
      ! do nn= 1,number_of_divisions
          ! xxx = 0._rp
          ! zzz = cell_start_z+(nn-1)*dzz
!        !$acc loop seq
        ! do m = 1,number_of_divisions
            ! yyy = cell_start_y + (m-1)*dyy
!            !$acc loop seq
            ! do l = 1,number_of_divisions
              ! xxx = cell_start_x + (l-1)*dxx
              ! inside = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                             ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
              ! counter = counter + real(inside,rp)
              ! call Solid_Surface(xxx,yyy,zzz,inside)
              ! if (inside) counter = counter +1
            ! enddo
        ! enddo
      ! enddo

      ! cell_u_tag(i,j,k) = counter/(1.0_rp*number_of_divisions**3)

! v cells
      ! cell_start_x = (i+ijk_start_x-1)*dx
      ! cell_end_x   = (i+ijk_start_x  )*dx
      ! cell_start_y = (j+ijk_start_y)*dy-0.5_rp*dy
      ! cell_end_y   = (j+ijk_start_y)*dy-0.5_rp*dy
      ! cell_start_z = zf(k-1)
      ! cell_end_z   = zf(k)
      ! dxx = (cell_end_x-cell_start_x)/number_of_divisions
      ! dyy = (cell_end_y-cell_start_y)/number_of_divisions
      ! dzz = (cell_end_z-cell_start_z)/number_of_divisions
      ! counter = 0
!      !$acc loop seq
      ! do nn= 1,number_of_divisions
          ! xxx = 0._rp
          ! zzz = cell_start_z+(nn-1)*dzz
!        !$acc loop seq
        ! do m = 1,number_of_divisions
            ! yyy = cell_start_y + (m-1)*dyy
!            !$acc loop seq
            ! do l = 1,number_of_divisions
              ! xxx = cell_start_x + (l-1)*dxx
              ! inside = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                             ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
              ! counter = counter + real(inside,rp)
              ! call Solid_Surface(xxx,yyy,zzz,inside)
              ! if (inside) counter = counter +1
            ! enddo
        ! enddo
      ! enddo

    ! cell_v_tag(i,j,k) = counter/(1.0_rp*number_of_divisions**3)

! w cells
      ! cell_start_x = (i+ijk_start_x-1)*dx
      ! cell_end_x   = (i+ijk_start_x  )*dx
      ! cell_start_y = (j+ijk_start_y-1)*dy
      ! cell_end_y   = (j+ijk_start_y  )*dy
      ! cell_start_z = zc(k-1)
      ! cell_end_z   = zc(k)
      ! dxx = (cell_end_x-cell_start_x)/number_of_divisions
      ! dyy = (cell_end_y-cell_start_y)/number_of_divisions
      ! dzz = (cell_end_z-cell_start_z)/number_of_divisions
      ! counter = 0._rp
!      !$acc loop seq
      ! do nn= 1,number_of_divisions
          ! xxx = 0._rp
          ! zzz = cell_start_z+(nn-1)*dzz
!        !$acc loop seq
        ! do m = 1,number_of_divisions
            ! yyy = cell_start_y + (m-1)*dyy
!            !$acc loop seq
            ! do l = 1,number_of_divisions
              ! xxx = cell_start_x + (l-1)*dxx
              ! inside = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                             ! xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
              ! counter = counter + real(inside,rp)
              ! call Solid_Surface(xxx,yyy,zzz,inside)
              ! if (inside) counter = counter +1
           ! enddo
        ! enddo
      ! enddo
      ! cell_w_tag(i,j,k) = counter/(1.0_rp*number_of_divisions**3)

    enddo
  enddo
enddo

end subroutine IBM_Mask

function wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                               xxx,yyy,zzz,dx,dy,lx,ly,dz,iperiod)
implicit none
logical :: wall_mounted_cylinder
real(rp), intent(in):: xxx,yyy,zzz,dx,dy,dz,RR,x_center,y_center,z_top,z_bottom,lx,ly
integer, dimension(3), intent(in) :: iperiod
real(rp):: sdist1,sdist2
integer :: ii,jj

wall_mounted_cylinder = .false.
sdist1 = (zzz - z_top)

     do jj = -1,1
       do ii = -1,1
         sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                        (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
        if( sdist1.lt.-dz .and. sdist2.lt.-min(dx,dy) ) then
          wall_mounted_cylinder = .true.
        elseif( sdist1.gt.dz  .and. sdist2.gt.min(dx,dy)  ) then
          wall_mounted_cylinder = .false.
        endif
        if ((zzz - z_bottom).le.dz) wall_mounted_cylinder = .true.
       enddo
     enddo

end function wall_mounted_cylinder

Subroutine GuessGhostCells(xxx,yyy,zzz,ghost)
implicit none
logical, intent(out) :: ghost
real(rp),intent(in)  :: xxx,yyy,zzz
real(rp):: RR, x_center, y_center, z_center, z_top, z_bottom, treshhold, temp
real(rp):: sdist1,sdist2,l_y,l_x
logical:: cond1,cond2,cond3,cond4,cond5
real(rp):: LL, AA, BB
real(rp):: x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4
real(rp):: z_1,z_2,z_3,z_4,abscissa
real(rp):: y_begin, y_end,x_begin,x_end
integer :: i1,j1,k1,ii,jj,q
real(rp):: A
real(rp):: y_12_mid,z_12_mid
real(rp):: Cos_angle,Sin_angle
real(rp):: edge_length_of_box
real(rp):: wall_inc = 0   !wall inclination
integer , dimension(3) :: iperiod

Cos_angle = cos(Rotation_angle)
Sin_angle = sin(Rotation_angle)

treshhold= 1.5*sqrt(dx**2 + dy**2 + dz**2)

i1=n(1)+1
j1=n(2)+1
k1=n(3)+1

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

     sdist1 = (zzz - z_top)
     do jj = -1,1
       do ii = -1,1
         sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                        (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
        if( sdist1.lt.-dz .and. sdist2.lt.-min(dx,dy) ) then
          ghost = .true.
        elseif( sdist1.gt.dz  .and. sdist2.gt.min(dx,dy)  ) then
          ghost = .false.
        endif
        if ((zzz - z_bottom).le.dz) ghost = .true.
       enddo
     enddo

case('Sphe')
     ghost=.false.
     RR =  lz/4.0_rp
     x_center = lx/2.0_rp
     y_center = ly/2.0_rp
     z_center = lz/8.0_rp

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

     cond1 = (abs(x1-xxx).lt.treshhold).or.(abs(x2-xxx).lt.treshhold)
     cond2 = (abs(y1-yyy).lt.treshhold).or.(abs(y2-yyy).lt.treshhold)
     cond3 = (abs(z1-zzz).lt.treshhold).or.(abs(z2-zzz).lt.treshhold)

     if (cond1.or.cond2.or.cond3) ghost=.true.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     LL = (zzz-z_center)**2 + (yyy-y_center)**2 + (xxx-x_center)**2
     if (LL.lt.RR**2) ghost=.true.

case('Plan')
     ghost=.false.
     x_center = xc_ibm
     y_center = yc_ibm
     z_center = zc_ibm

     cond3 = (abs(z_center-zzz).lt.treshhold)
     if (cond3) ghost=.true.

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
   
   ! cond1 = (abs(height - zzz).lt.treshhold)
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

Subroutine Solid_Surface(xxx,yyy,zzz,inside)
implicit none
logical, intent(out) :: inside
real(rp),intent(in)  :: xxx,yyy,zzz
real(rp):: LL, RR,x_begin,x_end, y_begin, y_end,x_center, y_center, z_center, z_top, z_bottom
real(rp):: sdist1,sdist2
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

Cos_angle = cos(Rotation_angle)
Sin_angle = sin(Rotation_angle)

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

     sdist1 = (zzz - z_top)
     do jj = -1,1
       do ii = -1,1
         sdist2 = sqrt( (xxx+ii*iperiod(1)*lx-x_center)**2 + &
                        (yyy+jj*iperiod(2)*ly-y_center)**2) - RR
        if( sdist1.lt.-dz .and. sdist2.lt.-min(dx,dy) ) then
          inside = .true.
        elseif( sdist1.gt.dz  .and. sdist2.gt.min(dx,dy)  ) then
          inside = .false.
        endif
        if ((zzz - z_bottom).le.dz) inside = .true.
       enddo
     enddo


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
     x_center = xc_ibm
     y_center = yc_ibm
     z_center = zc_ibm
     if (zzz.le.z_center) inside=.true.


case('RotB')
   inside= 0
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
   inside= 0
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
   if (.not.cond3) inside=1

end select

return

end subroutine Solid_Surface

Subroutine normal_vectors(Level_set,cell_phi_tag,nx_surf,ny_surf,nz_surf,nabs_surf,dl,dli,zc,dzc)
implicit none
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in)  :: cell_phi_tag
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in)  :: Level_set
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out) :: nx_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out) :: ny_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out) :: nz_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out) :: nabs_surf
real(rp), dimension(0:), intent(in) :: zc,dzc
real(rp), dimension(3), intent(in)  :: dl,dli
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: ghost_cell_tag
integer, dimension(3) :: iperiod
integer::i,j,k,l
real(rp):: nx, ny, nz, n_abs,n_abs_p
real(rp):: m_av_x,m_av_y,m_av_z,normal_denum
real(rp):: mx1,mx2,mx3,mx4,mx5,mx6,mx7,mx8
real(rp):: my1,my2,my3,my4,my5,my6,my7,my8
real(rp):: mz1,mz2,mz3,mz4,mz5,mz6,mz7,mz8
real(rp):: cell_tag_wall = 0.5_rp
real(rp) :: dx, dy, dz, dli1, dli2, dli3, dl1, dl2, dl3, lxl, lyl
integer:: ip,jp,kp
integer:: im,jm,km
integer:: iii,jjj,kkk
real(rp)::xxx,yyy,zzz
real(rp):: RR,x_center,y_center,z_center,z_top,z_bottom
integer:: Level_set_all
logical:: ghost,inside,ghost_cond
integer:: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z,q
dx=dl(1)
dy=dl(2)
!
n1 = n(1)
n2 = n(2)
n3 = n(3)
!
lxl = lx
lyl = ly
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
!! Preparing normal vectors to the surfaces
! Identifying the ghost cells
#if defined(_OPENACC)
if (myid == 0) print*, '*** Calculating Normal Vectors ***'
#endif
ghost_cell_tag(:,:,:) = 0._rp
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid == 0) print*, '*** Calculating Normal Vectors at k = ', k
#endif
 do j=1,n2
   do i=1,n1
       ghost = .false.
       xxx   =  (i+ijk_start_x)*dx-0.5_rp*dx
       yyy   =  (j+ijk_start_y)*dy-0.5_rp*dy
       zzz   =  zc(k)
       ! dz    =  dzc(k)
#if defined(_OPENACC)
       ghost = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                     xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
#else
       call GuessGhostCells(xxx,yyy,zzz,ghost)
#endif
       if (.not.ghost)  cycle
       ip = i+1
       jp = j+1
       kp = k+1
       im = i-1
       jm = j-1
       km = k-1
       ghost_cond = .false.
       Level_set_all = Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)  + &
                       Level_set(i,jm,kp) + Level_set(i,j,kp)+Level_set(i,jp,kp) + &
                       Level_set(ip,jm,kp) + Level_set(ip,j,kp)+Level_set(ip,jp,kp)+ &
                       Level_set(im,jm,k) + Level_set(im,j,k)+Level_set(im,jp,k)  + &
                       Level_set(i,jm,k) + Level_set(i,j,k)+Level_set(i,jp,k) + &
                       Level_set(ip,jm,k) + Level_set(ip,j,k)+Level_set(ip,jp,k)+ &
                       Level_set(im,jm,kp) + Level_set(im,j,kp)+Level_set(im,jp,kp)  + &
                       Level_set(i,jm,km) + Level_set(i,j,km)+Level_set(i,jp,km) + &
                       Level_set(ip,jm,km) + Level_set(ip,j,km)+Level_set(ip,jp,km)

       if ((Level_set_all.gt.1).and.(Level_set(i,j,k)).eq.1) ghost_cond =.true.
       ghost_cell_tag(i,j,k) = ghost_cond
   enddo
 enddo
enddo
do k=1,n3
 do j=1,n2
  do i=1,n1
   if (ghost_cell_tag(i,j,k) == 1) then

      nx     = 0.0_rp
      ny     = 0.0_rp
      nz     = 0.0_rp
      n_abs  = 0.0_rp
      m_av_x = 0.0_rp
      m_av_y = 0.0_rp
      m_av_z = 0.0_rp
      normal_denum = 0.0_rp

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
      m_av_x= 0.125_rp*(mx1+mx2+mx3+mx4+mx5+mx6+mx7+mx8)
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
      m_av_y= 0.125_rp*(my1+my2+my3+my4+my5+my6+my7+my8)
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
      m_av_z= 0.125_rp*(mz1+mz2+mz3+mz4+mz5+mz6+mz7+mz8)

      normal_denum = sqrt(m_av_x**2 + m_av_y**2 + m_av_z**2)

      nx_surf(i,j,k) = m_av_x/max(normal_denum, small)
      ny_surf(i,j,k) = m_av_y/max(normal_denum, small)
      nz_surf(i,j,k) = m_av_z/max(normal_denum, small)

      nabs_surf(i,j,k) = sqrt(nx_surf(i,j,k)**2 + ny_surf(i,j,k)**2 + nz_surf(i,j,k)**2)

     endif

    enddo
  enddo
enddo
return

end subroutine normal_vectors
!
Subroutine intersect(nx_surf,ny_surf,nz_surf,nabs_surf,x_intersect,y_intersect,z_intersect,dl,dzc,zc)
implicit none
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nx_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: ny_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nz_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nabs_surf
real(rp), dimension(0:), intent(in) :: dzc,zc
real(rp), dimension(3) , intent(in) :: dl
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: x_intersect
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: y_intersect
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: z_intersect
integer, dimension(3) :: iperiod
integer::i,j,k,l
real(rp):: nx, ny, nz, n_abs
real(rp):: step
real(rp):: xxx,yyy,zzz,dx,dy,dz,lxl,lyl
logical :: inside
logical :: confirmation
integer:: lmax
real(rp)::distance_ghost_intersect
real(rp):: x_ghost,y_ghost,z_ghost
real(rp):: RR,x_center,y_center,z_center,z_top,z_bottom
integer:: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z,q
!
dx=dl(1)
dy=dl(2)
!
n1=n(1)
n2=n(2)
n3=n(3)
!
lxl = lx
lyl = ly
!
do q=1,3
  if(cbcpre(0,q)//cbcpre(1,q).eq.'PP') iperiod(q) = 1
enddo
!
RR =  r_ibm
x_center = xc_ibm
y_center = yc_ibm
z_bottom = zmin_ibm
z_top    = zmax_ibm
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
distance_ghost_intersect= 1000.0_rp
!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
#if defined(_OPENACC)
if (myid == 0) print*, '*** Calculating Intersect Points ***'
#endif
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
 if (myid == 0) print*, '*** Calculating Intersect Points at k = ', k
#endif
 step = 1.0e-5_rp*dzc(k)
 lmax=100*int(2.0_rp*sqrt(3.0_rp)*dzc(k)/step)
 dz = dzc(k)
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then
       confirmation = .false.
       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       x_ghost = (i+ijk_start_x)*dx-0.5_rp*dx
       y_ghost = (j+ijk_start_y)*dy-0.5_rp*dy
       z_ghost = zc(k)
       do l=0,lmax
         yyy = y_ghost+l*(ny/(n_abs+small))*step
         zzz = z_ghost+l*(nz/(n_abs+small))*step
         xxx = x_ghost+l*(nx/(n_abs+small))*step
#if defined(_OPENACC)
         inside = wall_mounted_cylinder(RR,x_center,y_center,z_bottom,z_top, &
                                        xxx,yyy,zzz,dx,dy,lxl,lyl,dz,iperiod)
#else
         call Solid_Surface(xxx,yyy,zzz,inside)
#endif
         if (.not.inside) then
           x_intersect(i,j,k) = xxx
           y_intersect(i,j,k) = yyy
           z_intersect(i,j,k) = zzz
           confirmation = .true.
          exit
         endif
      enddo
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

    endif
   enddo
  enddo
enddo

return
end subroutine intersect

Subroutine mirrorpoints(nx_surf,ny_surf,nz_surf,nabs_surf, &
                        x_intersect,y_intersect,z_intersect, &
                        x_mirror,y_mirror,z_mirror, &
                        deltan, &
                        x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                        dl,dzc,zc)
implicit none
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nx_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: ny_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nz_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: nabs_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: x_intersect
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: y_intersect
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in) :: z_intersect
real(rp), dimension(0:), intent(in) :: dzc,zc
real(rp), dimension(3) , intent(in) :: dl
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: x_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: y_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: z_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: deltan
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: x_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: y_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: z_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: x_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: y_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: z_IP2
integer::i,j,k,l,m,nn
real(rp):: nx, ny, nz, n_abs
real(rp):: step,step_2
real(rp):: xxx,yyy,zzz
real(rp):: cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp)::distance_ghost_intersect
integer:: n1,n2,n3,i1,j1,k1,ijk_start_x,ijk_start_y,ijk_start_z
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
!***************************************************************************************
! Part one: cell-centered
!***************************************************************************************
#if defined(_OPENACC)
if (myid == 0) print*, '*** Calculating Mirror Points ***'
#endif
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating Mirror Points at k = ', k
#endif
step   = sqrt(dx**2 + dy**2 + dzc(k)**2)
step_2 = 0.5_rp*sqrt(dx**2 + dy**2 + dzc(k)**2)
 do j= 1,n2
   do i= 1,n1
    if (nabs_surf(i,j,k).gt.small) then

       xxx = (i+ijk_start_x)*dx-0.5_rp*dx
       yyy = (j+ijk_start_y)*dy-0.5_rp*dy
       zzz = zc(k)

       nx =   nx_surf(i,j,k)
       ny =   ny_surf(i,j,k)
       nz =   nz_surf(i,j,k)
       n_abs = nabs_surf(i,j,k)
       distance_ghost_intersect = sqrt((x_intersect(i,j,k)-xxx)**2 + &
                                       (y_intersect(i,j,k)-yyy)**2 + &
                                       (z_intersect(i,j,k)-zzz)**2)


       if  (distance_ghost_intersect.gt.sqrt(dx**2 + dy**2 + dzc(k)**2)) then
       print*, '--------------------------------------------------------------------------'
           print*, ' Error: in mirror point detection at cell-center at processor ', &
            myid,' : check IBM.f90 - distance_ghost_intesect is ',distance_ghost_intersect, &
            'where the normal vector components are ', nx,ny,nz,n_abs
       print*, '--------------------------------------------------------------------------'
       endif

       x_mirror(i,j,k) = x_intersect(i,j,k)+(nx/(n_abs+small))*step
       y_mirror(i,j,k) = y_intersect(i,j,k)+(ny/(n_abs+small))*step !change for 3D
       z_mirror(i,j,k) = z_intersect(i,j,k)+(nz/(n_abs+small))*step !change for 3D

       deltan(i,j,k)   = distance_ghost_intersect

       x_IP1(i,j,k)    = x_mirror(i,j,k)+(nx/(n_abs+small))*step_2
       y_IP1(i,j,k)    = y_mirror(i,j,k)+(ny/(n_abs+small))*step_2
       z_IP1(i,j,k)    = z_mirror(i,j,k)+(nz/(n_abs+small))*step_2

       x_IP2(i,j,k)    = x_IP1(i,j,k)+(nx/(n_abs+small))*step_2
       y_IP2(i,j,k)    = y_IP1(i,j,k)+(ny/(n_abs+small))*step_2
       z_IP2(i,j,k)    = z_IP1(i,j,k)+(nz/(n_abs+small))*step_2

    endif
   enddo
  enddo
enddo

return
end subroutine mirrorpoints

Subroutine mirrorpoints_ijk(nabs_surf,x_mirror,y_mirror,z_mirror, &
                            x_IP1,y_IP1,z_IP1,x_IP2,y_IP2,z_IP2, &
                            i_mirror,j_mirror,k_mirror, &
                            i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                            dl,zf)

implicit none
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: nabs_surf
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: x_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: y_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: z_mirror
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: x_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: y_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: z_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: x_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: y_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: z_IP2
real(rp), dimension(0:), intent(in) :: zf
real(rp), dimension(3) , intent(in) :: dl
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: i_mirror
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: j_mirror
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: k_mirror
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: i_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: j_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: k_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: i_IP2
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: j_IP2
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(out):: k_IP2
integer::i,j,k,l,m,nn
real(rp):: nx, ny, nz, n_abs
real(rp):: xxx,yyy,zzz
real(rp)::cell_start_x,cell_end_x,cell_start_y,cell_end_y,cell_start_z,cell_end_z
real(rp)::distance_ghost_intesect
integer:: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z
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
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating mirror point indices ***'
#endif
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating mirror points indices at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
   !***************************************************************************************
   ! Part one: cell-centered
   !***************************************************************************************
    if (nabs_surf(i,j,k).gt.small) then
    ! Mirror points
       xxx = x_mirror(i,j,k)
       yyy = y_mirror(i,j,k)
       zzz = z_mirror(i,j,k)
       do l = 1,n3
         do m = -5,n2+6
           do nn = -5,n1+6
             cell_start_x = (nn+ijk_start_x-1)*dx
             cell_end_x   = (nn+ijk_start_x  )*dx
             cell_start_y = (m +ijk_start_y-1)*dy
             cell_end_y   = (m +ijk_start_y  )*dy
             cell_start_z = zf(k-1)
             cell_end_z   = zf(k)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_mirror(i,j,k) = nn
                  j_mirror(i,j,k) = m
                  k_mirror(i,j,k) = l
              exit
             endif
            enddo
          enddo
       enddo
       if ((i_mirror(i,j,k).eq.i).and.(j_mirror(i,j,k).eq.j).and.(k_mirror(i,j,k).eq.k)) then
             print*, 'Error: Ghost and mirror point are the same'
       endif
       if ((i_mirror(i,j,k).eq.-1000).or.(j_mirror(i,j,k).eq.-1000).or.(k_mirror(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for mirror point at center i = ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! First interpolation points
       xxx = x_IP1(i,j,k)
       yyy = y_IP1(i,j,k)
       zzz = z_IP1(i,j,k)
       do l= 1,n3
         do m= -5,n2+6
           do nn= -5,n1+6
             cell_start_x = (nn+ijk_start_x-1)*dx
             cell_end_x   = (nn+ijk_start_x  )*dx
             cell_start_y = (m +ijk_start_y-1)*dy
             cell_end_y   = (m +ijk_start_y  )*dy
             cell_start_z = zf(k-1)
             cell_end_z   = zf(k)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP1(i,j,k) = nn
                  j_IP1(i,j,k) = m
                  k_IP1(i,j,k) = l
              exit
             endif
            enddo
          enddo
       enddo


       if ((i_IP1(i,j,k).eq.-1000).or.(j_IP1(i,j,k).eq.-1000).or.(k_IP1(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the first interpolation point at center i = '&
         ,i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    ! Second interpolation points
       xxx = x_IP2(i,j,k)
       yyy = y_IP2(i,j,k)
       zzz = z_IP2(i,j,k)
       do l= 1,n3
         do m= -5,n2+6
           do nn= -5,n1+6
             cell_start_x = (nn+ijk_start_x-1)*dx
             cell_end_x   = (nn+ijk_start_x  )*dx
             cell_start_y = (m +ijk_start_y-1)*dy
             cell_end_y   = (m +ijk_start_y  )*dy
             cell_start_z = zf(k-1)
             cell_end_z   = zf(k)
             if ((yyy.ge.cell_start_y).and.(yyy.lt.cell_end_y).and.&
                 (zzz.ge.cell_start_z).and.(zzz.lt.cell_end_z).and.&
                 (xxx.ge.cell_start_x).and.(xxx.lt.cell_end_x)) then
                  i_IP2(i,j,k) = nn
                  j_IP2(i,j,k) = m
                  k_IP2(i,j,k) = l
              exit
             endif

            enddo
          enddo
       enddo

       if ((i_IP2(i,j,k).eq.-1000).or.(j_IP2(i,j,k).eq.-1000).or.(k_IP2(i,j,k).eq.-1000)) then
       print*, '--------------------------------------------------------------------------'
         print*,'Error: no grid point detected for the second interpolation  point at center i = ', &
         i, ' j= ',j,' k= ',k, 'at processor ',myid
       print*, '--------------------------------------------------------------------------'
      endif

    endif

  enddo
 enddo
enddo
return
end Subroutine mirrorpoints_ijk
!
Subroutine InterpolationWeights(nabs_surf,Level_set, &
                                x_IP1,y_IP1,z_IP1, &
                                x_IP2,y_IP2,z_IP2, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2, &
                                dl,dzc)
implicit none
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: nabs_surf
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: Level_set
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: x_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: y_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: z_IP1
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: x_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: y_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: z_IP2
real(rp), dimension(0:), intent(in):: dzc
real(rp), dimension(3) , intent(in):: dl
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: i_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: j_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: k_IP1
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: i_IP2
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: j_IP2
integer , dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: k_IP2
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7), intent(out) :: WP1,WP2
real(rp):: x1(7),y1(7),z1(7),x2(7),y2(7),z2(7)
real(rp):: h1(7),h2(7)
integer :: i_p_1(7),j_p_1(7),k_p_1(7),i_p_2(7),j_p_2(7),k_p_2(7)
integer :: ii,i,j,k,ii1,ii2,jj1,jj2,kk1,kk2
real(rp):: xx1,xx2,yy1,yy2,zz1,zz2
logical :: cond11,cond21,cond12,cond22
real(rp):: contribution1(7),contribution2(7)
integer :: n1,n2,n3,ijk_start_x,ijk_start_y,ijk_start_z
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
do k=1,int(solid_height_ratio*n3)
#if !defined(_OPENACC)
if (myid.eq.0) print*, '*** Calculating Interpolation Weights at k = ', k
#endif
 do j= 1,n2
   do i= 1,n1
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
           cond21 = (nabs_surf(i_p_1(ii),j_p_1(ii),k_p_1(ii)).eq.0.0_rp)
           cond12 = (Level_set(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.1)
           cond22 = (nabs_surf(i_p_2(ii),j_p_2(ii),k_p_2(ii)).eq.0.0_rp)

           if (cond11) contribution1(ii) = 1.0_rp
           if (cond12) contribution2(ii) = 1.0_rp

           x1(ii) = (i_p_1(ii)+ijk_start_x-0.5_rp)*dx
           y1(ii) = (j_p_1(ii)+ijk_start_y-0.5_rp)*dy
           z1(ii) = zc(k_p_1(ii))
           x2(ii) = (i_p_2(ii)+ijk_start_x-0.5_rp)*dx
           y2(ii) = (j_p_2(ii)+ijk_start_y-0.5_rp)*dy
           z2(ii) = zc(k_p_2(ii))
        enddo

        do ii = 1,7
         h1(ii) = sqrt( (x1(ii)-xx1)**2 + (y1(ii)-yy1)**2 + (z1(ii)-zz1)**2 )
         h2(ii) = sqrt( (x2(ii)-xx2)**2 + (y2(ii)-yy2)**2 + (z2(ii)-zz2)**2 )
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
return
end subroutine InterpolationWeights

subroutine interpolation_mirror(AA,iii,jjj,kkk, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2,BB)
implicit none
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7),intent(in):: WP1,WP2
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: AA
integer, dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)    :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
integer, intent(in)  :: iii,jjj,kkk
real(rp),intent(out) :: BB
real(rp),dimension(7) ::WW1,WW2,B1,B2
real(rp) :: q1,q2,B_1,B_2
integer  :: ii1,jj1,kk1,ii2,jj2,kk2,l

!initialization
q1       = 0.0_rp
q2       = 0.0_rp
WW1(:)   = 0.0_rp
WW2(:)   = 0.0_rp
B1(:)    = 0.0_rp
B2(:)    = 0.0_rp
B_1      = 0.0_rp
B_2      = 0.0_rp
BB       = 0.0_rp

ii1    = i_IP1(iii,jjj,kkk) ! Coordinates for cell of interpolation point 1 (IP1)
jj1    = j_IP1(iii,jjj,kkk)
kk1    = k_IP1(iii,jjj,kkk)
ii2    = i_IP2(iii,jjj,kkk) ! Coordinates for cell of interpolation point 2  (IP2)
jj2    = j_IP2(iii,jjj,kkk)
kk2    = k_IP2(iii,jjj,kkk)


WW1(1) = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
WW1(2) = WP1(iii,jjj,kkk,2)
WW1(3) = WP1(iii,jjj,kkk,3)
WW1(4) = WP1(iii,jjj,kkk,4)
WW1(5) = WP1(iii,jjj,kkk,5)
WW1(6) = WP1(iii,jjj,kkk,6)
WW1(7) = WP1(iii,jjj,kkk,7)
WW2(1) = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP2
WW2(2) = WP2(iii,jjj,kkk,2)
WW2(3) = WP2(iii,jjj,kkk,3)
WW2(4) = WP2(iii,jjj,kkk,4)
WW2(5) = WP2(iii,jjj,kkk,5)
WW2(6) = WP2(iii,jjj,kkk,6)
WW2(7) = WP2(iii,jjj,kkk,7)
B1(1)   = AA(ii1,jj1+1,kk1)  ! Values for interpolating to IP1
B1(2)   = AA(ii1,jj1,kk1+1)
B1(3)   = AA(ii1,jj1-1,kk1)
B1(4)   = AA(ii1,jj1,kk1-1)
B1(5)   = AA(ii1,jj1,kk1)
B1(6)   = AA(ii1-1,jj1,kk1)
B1(7)   = AA(ii1+1,jj1,kk1)
B2(1)   = AA(ii2,jj2+1,kk2)  ! Values for interpolating to IP2
B2(2)   = AA(ii2,jj2,kk2+1)
B2(3)   = AA(ii2,jj2-1,kk2)
B2(4)   = AA(ii2,jj2,kk2-1)
B2(5)   = AA(ii2,jj2,kk2)
B2(6)   = AA(ii2-1,jj2,kk2)
B2(7)   = AA(ii2+1,jj2,kk2)

do l = 1,7 ! Sum weights
   q1 =  max(WW1(l) + q1, small)
   q2 =  max(WW2(l) + q2, small)
enddo

do l = 1,7 ! Compute value at IP1 and IP2
   B_1 = (1.0_rp/q1)*(WW1(l)*B1(l))+B_1
   B_2 = (1.0_rp/q2)*(WW2(l)*B2(l))+B_2
enddo

BB = 2.0_rp*B_1-B_2 ! Value at mirror point

end subroutine interpolation_mirror
!
Subroutine interpolation_2D_velocity(UU,VV,WW,iii,jjj,kkk, &
                                     i_IP1,j_IP1,k_IP1, &
                                     i_IP2,j_IP2,k_IP2, &
                                     WP1,WP2,U_m,V_m,W_m)

implicit none
integer,intent(in) :: iii,jjj,kkk
integer,intent(in)::i_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::i_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
real(rp),intent(in)::     WP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
real(rp),intent(in)::     WP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
!
real(rp),intent(out) ::U_m, V_m, W_m
real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), intent(in):: UU ,VV ,WW
real(rp) :: q1,q2,WW1(7),WW2(7),U1(7),U2(7),V1(7),V2(7),W1(7),W2(7)
real(rp)::  U_1,U_2,V_1,V_2,W_1,W_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2,l
integer:: i1,j1,k1
i1=n(1)+1
j1=n(2)+1
k1=n(3)+1
!
!initialization
q1      = 0.0_rp
q2      = 0.0_rp
WW1(:)  = 0.0_rp
WW2(:)  = 0.0_rp
U1(:)   = 0.0_rp
U2(:)   = 0.0_rp
V1(:)   = 0.0_rp
V2(:)   = 0.0_rp
W1(:)   = 0.0_rp
W2(:)   = 0.0_rp
U_1     = 0.0_rp
U_2     = 0.0_rp
V_1     = 0.0_rp
V_2     = 0.0_rp
W_1     = 0.0_rp
W_2     = 0.0_rp

ii1     = i_IP1(iii,jjj,kkk)
jj1     = j_IP1(iii,jjj,kkk)
kk1     = k_IP1(iii,jjj,kkk)
ii2     = i_IP2(iii,jjj,kkk)
jj2     = j_IP2(iii,jjj,kkk)
kk2     = k_IP2(iii,jjj,kkk)

WW1(1) = WP1(iii,jjj,kkk,1)
WW1(2) = WP1(iii,jjj,kkk,2)
WW1(3) = WP1(iii,jjj,kkk,3)
WW1(4) = WP1(iii,jjj,kkk,4)
WW1(5) = WP1(iii,jjj,kkk,5)
WW1(6) = WP1(iii,jjj,kkk,6)
WW1(7) = WP1(iii,jjj,kkk,7)
WW2(1) = WP2(iii,jjj,kkk,1)
WW2(2) = WP2(iii,jjj,kkk,2)
WW2(3) = WP2(iii,jjj,kkk,3)
WW2(4) = WP2(iii,jjj,kkk,4)
WW2(5) = WP2(iii,jjj,kkk,5)
WW2(6) = WP2(iii,jjj,kkk,6)
WW2(7) = WP2(iii,jjj,kkk,7)

U1(1)    = UU(ii1,jj1+1,kk1)
U1(2)    = UU(ii1,jj1,  kk1+1)
U1(3)    = UU(ii1,jj1-1,kk1)
U1(4)    = UU(ii1,jj1,  kk1-1)
U1(5)    = UU(ii1,jj1,  kk1)
U1(6)    = UU(ii1-1,jj1,  kk1)
U1(7)    = UU(ii1+1,jj1,  kk1)
U2(1)    = UU(ii2,jj2+1,kk2)
U2(2)    = UU(ii2,jj2,  kk2+1)
U2(3)    = UU(ii2,jj2-1,kk2)
U2(4)    = UU(ii2,jj2,  kk2-1)
U2(5)    = UU(ii2,jj2,  kk2)
U2(6)    = UU(ii2-1,jj2,  kk2)
U2(7)    = UU(ii2+1,jj2,  kk2)

V1(1)    = VV(ii1,jj1+1,kk1)
V1(2)    = VV(ii1,jj1,  kk1+1)
V1(3)    = VV(ii1,jj1-1,kk1)
V1(4)    = VV(ii1,jj1,  kk1-1)
V1(5)    = VV(ii1,jj1,  kk1)
V1(6)    = VV(ii1-1,jj1,  kk1)
V1(7)    = VV(ii1+1,jj1,  kk1)
V2(1)    = VV(ii2,jj2+1,kk2)
V2(2)    = VV(ii2,jj2,  kk2+1)
V2(3)    = VV(ii2,jj2-1,kk2)
V2(4)    = VV(ii2,jj2,  kk2-1)
V2(5)    = VV(ii2,jj2,  kk2)
V2(6)    = VV(ii2-1,jj2,  kk2)
V2(7)    = VV(ii2+1,jj2,  kk2)

W1(1)    = WW(ii1,jj1+1,kk1)
W1(2)    = WW(ii1,jj1,  kk1+1)
W1(3)    = WW(ii1,jj1-1,kk1)
W1(4)    = WW(ii1,jj1,  kk1-1)
W1(5)    = WW(ii1,jj1,  kk1)
W1(6)    = WW(ii1-1,jj1,  kk1)
W1(7)    = WW(ii1+1,jj1,  kk1)
W2(1)    = WW(ii2,jj2+1,kk2)
W2(2)    = WW(ii2,jj2,  kk2+1)
W2(3)    = WW(ii2,jj2-1,kk2)
W2(4)    = WW(ii2,jj2,  kk2-1)
W2(5)    = WW(ii2,jj2,  kk2)
W2(6)    = WW(ii2-1,jj2,  kk2)
W2(7)    = WW(ii2+1,jj2,  kk2)

do l = 1,7
   q1 =  max(WW1(l) + q1, small)
   q2 =  max(WW2(l) + q2, small)
enddo

do l = 1,7
 U_1 = (1.0_rp/q1)*(WW1(l)*U1(l))+U_1
 U_2 = (1.0_rp/q2)*(WW2(l)*U2(l))+U_2
 V_1 = (1.0_rp/q1)*(WW1(l)*V1(l))+V_1
 V_2 = (1.0_rp/q2)*(WW2(l)*V2(l))+V_2
 W_1 = (1.0_rp/q1)*(WW1(l)*W1(l))+W_1
 W_2 = (1.0_rp/q2)*(WW2(l)*W2(l))+W_2
enddo

u_m = 2.0_rp*U_1-U_2
V_m = 2.0_rp*V_1-V_2
W_m = 2.0_rp*W_1-W_2

return
end subroutine interpolation_2D_velocity


Subroutine interpolation_dphi(PFM_phi,iii,jjj,kkk, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                          WP1,WP2, &
                              dphidx,dphidy,dphidz)
implicit none
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6):: PFM_phi
integer,intent(in) :: iii,jjj,kkk
integer,intent(in)::i_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::i_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
real(rp),intent(in)::     WP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
real(rp),intent(in)::     WP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
!
real(rp),intent(out):: dphidx,dphidy,dphidz
integer::ii,l
real(rp) :: q1,q2,WW1(7),WW2(7)!,B1(5),B2(5),B_1,B_2
integer:: ii1,jj1,kk1,ii2,jj2,kk2
real(rp):: phi_ip1(7),phi_im1(7), phi_jp1(7),phi_jm1(7),phi_kp1(7),phi_km1(7)
real(rp):: phi_ip2(7),phi_im2(7),phi_jp2(7),phi_jm2(7),phi_kp2(7),phi_km2(7)
real(rp):: phi_x1(7),phi_y1(7),phi_z1(7),phi_x2(7),phi_y2(7),phi_z2(7)
real(rp):: dphidx1,dphidy1,dphidz1,dphidx2,dphidy2,dphidz2
integer:: i1,j1,k1
i1=n(1)+1
j1=n(2)+1
k1=n(3)+1
!
!initialization
q1          = 0.0_rp
q2          = 0.0_rp
WW1(:)      = 0.0_rp
WW2(:)      = 0.0_rp
phi_ip1(:)  = 0.0_rp
phi_im1(:)  = 0.0_rp
phi_jp1(:)  = 0.0_rp
phi_jm1(:)  = 0.0_rp
phi_kp1(:)  = 0.0_rp
phi_km1(:)  = 0.0_rp
phi_ip2(:)  = 0.0_rp
phi_im2(:)  = 0.0_rp
phi_jp2(:)  = 0.0_rp
phi_jm2(:)  = 0.0_rp
phi_kp2(:)  = 0.0_rp
phi_km2(:)  = 0.0_rp
phi_x1(:)   = 0.0_rp
phi_y1(:)   = 0.0_rp
phi_z1(:)   = 0.0_rp
phi_x2(:)   = 0.0_rp
phi_y2(:)   = 0.0_rp
phi_z2(:)   = 0.0_rp
dphidx1     = 0.0_rp
dphidy1     = 0.0_rp
dphidz1     = 0.0_rp
dphidx2     = 0.0_rp
dphidy2     = 0.0_rp
dphidz2     = 0.0_rp

  ii1     = i_IP1(iii,jjj,kkk)
  jj1     = j_IP1(iii,jjj,kkk)
  kk1     = k_IP1(iii,jjj,kkk)
  ii2     = i_IP2(iii,jjj,kkk)
  jj2     = j_IP2(iii,jjj,kkk)
  kk2     = k_IP2(iii,jjj,kkk)


  WW1(1) = WP1(iii,jjj,kkk,1) ! Weights for interpolating to IP1
  WW1(2) = WP1(iii,jjj,kkk,2)
  WW1(3) = WP1(iii,jjj,kkk,3)
  WW1(4) = WP1(iii,jjj,kkk,4)
  WW1(5) = WP1(iii,jjj,kkk,5)
  WW1(6) = WP1(iii,jjj,kkk,6)
  WW1(7) = WP1(iii,jjj,kkk,7)
  WW2(1) = WP2(iii,jjj,kkk,1) ! Weights for interpolating to IP2
  WW2(2) = WP2(iii,jjj,kkk,2)
  WW2(3) = WP2(iii,jjj,kkk,3)
  WW2(4) = WP2(iii,jjj,kkk,4)
  WW2(5) = WP2(iii,jjj,kkk,5)
  WW2(6) = WP2(iii,jjj,kkk,6)
  WW2(7) = WP2(iii,jjj,kkk,7)

  ! Cell edge values for AP for IP1
  phi_ip1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_im1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1-1,jj1+1,kk1))
  phi_jp1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+2,kk1))
  phi_jm1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_kp1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_km1(1) = 0.5_rp*(PFM_phi(ii1,jj1+1,kk1)+PFM_phi(ii1,jj1+1,kk1-1))


  phi_ip1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_im1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1-1,jj1-1,kk1+1))
  phi_jp1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1+1,kk1+1))
  phi_jm1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_kp1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1+2))
  phi_km1(2) = 0.5_rp*(PFM_phi(ii1,jj1,kk1+1)+PFM_phi(ii1,jj1,kk1))


  phi_ip1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_im1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_jp1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jm1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-2,kk1))
  phi_kp1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1+1))
  phi_km1(3) = 0.5_rp*(PFM_phi(ii1,jj1-1,kk1)+PFM_phi(ii1,jj1-1,kk1-1))

  phi_ip1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1+1,jj1,kk1-1))
  phi_im1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1-1,jj1-1,kk1-1))
  phi_jp1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1+1,kk1-1))
  phi_jm1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1-1,kk1-1))
  phi_kp1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1))
  phi_km1(4) = 0.5_rp*(PFM_phi(ii1,jj1,kk1-1)+PFM_phi(ii1,jj1,kk1-2))


  phi_ip1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1))
  phi_im1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1))
  phi_jp1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1+1,kk1))
  phi_jm1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1-1,kk1))
  phi_kp1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1+1))
  phi_km1(5) = 0.5_rp*(PFM_phi(ii1,jj1,kk1)+PFM_phi(ii1,jj1,kk1-1))


  phi_ip1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_im1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-2,jj1,kk1))
  phi_jp1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1+1,kk1))
  phi_jm1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1-1,kk1))
  phi_kp1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1+1))
  phi_km1(6) = 0.5_rp*(PFM_phi(ii1-1,jj1,kk1)+PFM_phi(ii1-1,jj1,kk1-1))


  phi_ip1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+2,jj1,kk1))
  phi_im1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1,jj1,kk1))
  phi_jp1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1+1,kk1))
  phi_jm1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1-1,kk1))
  phi_kp1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1+1))
  phi_km1(7) = 0.5_rp*(PFM_phi(ii1+1,jj1,kk1)+PFM_phi(ii1+1,jj1,kk1-1))


  ! Cell edge values for AP for IP2
  phi_ip2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_im2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jp2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+2,kk2))
  phi_jm2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_kp2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_km2(1) = 0.5_rp*(PFM_phi(ii2,jj2+1,kk2)+PFM_phi(ii2,jj2+1,kk2-1))


  phi_ip2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_im2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2-1,jj2-1,kk2+1))
  phi_jp2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2+1,kk2+1))
  phi_jm2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_kp2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2+2))
  phi_km2(2) = 0.5_rp*(PFM_phi(ii2,jj2,kk2+1)+PFM_phi(ii2,jj2,kk2))


  phi_ip2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_im2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_jp2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jm2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-2,kk2))
  phi_kp2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2+1))
  phi_km2(3) = 0.5_rp*(PFM_phi(ii2,jj2-1,kk2)+PFM_phi(ii2,jj2-1,kk2-1))

  phi_ip2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2+1,jj2,kk2-1))
  phi_im2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2-1,jj2-1,kk2-1))
  phi_jp2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2+1,kk2-1))
  phi_jm2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2-1,kk2-1))
  phi_kp2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2))
  phi_km2(4) = 0.5_rp*(PFM_phi(ii2,jj2,kk2-1)+PFM_phi(ii2,jj2,kk2-2))


  phi_ip2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2))
  phi_im2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2))
  phi_jp2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2+1,kk2))
  phi_jm2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2-1,kk2))
  phi_kp2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2+1))
  phi_km2(5) = 0.5_rp*(PFM_phi(ii2,jj2,kk2)+PFM_phi(ii2,jj2,kk2-1))


  phi_ip2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_im2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-2,jj2,kk2))
  phi_jp2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2+1,kk2))
  phi_jm2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2-1,kk2))
  phi_kp2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2+1))
  phi_km2(6) = 0.5_rp*(PFM_phi(ii2-1,jj2,kk2)+PFM_phi(ii2-1,jj2,kk2-1))


  phi_ip2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+2,jj2,kk2))
  phi_im2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2,jj2,kk2))
  phi_jp2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2+1,kk2))
  phi_jm2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2-1,kk2))
  phi_kp2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2+1))
  phi_km2(7) = 0.5_rp*(PFM_phi(ii2+1,jj2,kk2)+PFM_phi(ii2+1,jj2,kk2-1))


  ! Derivatives at AP
  do ii=1,7
    phi_x1(ii) = (phi_ip1(ii)-phi_im1(ii))*dxi
    phi_y1(ii) = (phi_jp1(ii)-phi_jm1(ii))*dyi
    phi_z1(ii) = (phi_kp1(ii)-phi_km1(ii))*dzi
    phi_x2(ii) = (phi_ip2(ii)-phi_im2(ii))*dxi
    phi_y2(ii) = (phi_jp2(ii)-phi_jm2(ii))*dyi
    phi_z2(ii) = (phi_kp2(ii)-phi_km2(ii))*dzi
    q1 =  max(WW1(ii) + q1, small)
    q2 =  max(WW2(ii) + q2, small)
  enddo

! Derivatives at IP1 and IP2
  do l = 1,7
   dphidx1 = (1.0_rp/q1)*(WW1(l)*phi_x1(l))+dphidx1
   dphidy1 = (1.0_rp/q1)*(WW1(l)*phi_y1(l))+dphidy1
   dphidz1 = (1.0_rp/q1)*(WW1(l)*phi_z1(l))+dphidz1
   dphidx2 = (1.0_rp/q2)*(WW2(l)*phi_x2(l))+dphidx2
   dphidy2 = (1.0_rp/q2)*(WW2(l)*phi_y2(l))+dphidy2
   dphidz2 = (1.0_rp/q2)*(WW2(l)*phi_z2(l))+dphidz2
  enddo

  ! Extrapolate to mirror point
  dphidx = 2.0_rp*dphidx1-dphidx2
  dphidy = 2.0_rp*dphidy1-dphidy2
  dphidz = 2.0_rp*dphidz1-dphidz2

return
end subroutine interpolation_dphi

subroutine Penalization_center(uu,vv,ww,cell_phi_tag,nabs_surf)
implicit none
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)::cell_phi_tag,nabs_surf
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(inout)::uu,vv,ww
integer::i,j,k,n1,n2,n3
!@cuf attributes(managed) :: uu,vv,ww,cell_phi_tag,nabs_surf
!
n1=n(1)
n2=n(2)
n3=n(3)

   !$acc kernels
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
         uu(i,j,k)=(1.0_rp - cell_phi_tag(i,j,k))*uu(i,j,k)
         vv(i,j,k)=(1.0_rp - cell_phi_tag(i,j,k))*vv(i,j,k)
         ww(i,j,k)=(1.0_rp - cell_phi_tag(i,j,k))*ww(i,j,k)
       ! endif
      enddo
     enddo
   enddo
   !$acc end kernels
end subroutine Penalization_center


subroutine Penalization_face(uu,vv,ww,cell_u_tag,cell_v_tag,cell_w_tag)
implicit none
real(rp),dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),intent(in)::cell_u_tag,cell_v_tag,cell_w_tag
real(rp),dimension(0:n(1)+1,0:n(2)+1,0:n(3)+1),intent(inout)::uu,vv,ww
integer::i,j,k,n1,n2,n3
!@cuf attributes(managed) :: uu,vv,ww,cell_u_tag,cell_v_tag,cell_w_tag
!
n1=n(1)
n2=n(2)
n3=n(3)
   !$acc kernels
   do k=1,n3
    do j=1,n2
     do i=1,n1
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
              uu(i,j,k)=(1.0_rp - cell_u_tag(i,j,k))*uu(i,j,k)
              vv(i,j,k)=(1.0_rp - cell_v_tag(i,j,k))*vv(i,j,k)
              ww(i,j,k)=(1.0_rp - cell_w_tag(i,j,k))*ww(i,j,k)

       ! endif
      enddo
     enddo
   enddo
   !$acc end kernels
end subroutine Penalization_face
!
subroutine Wetting_radius(time,nabs_surf,PFM_phi,deltan, &
                          i_IP1,j_IP1,k_IP1, &
                          i_IP2,j_IP2,k_IP2,WP1,WP2, &
                          dzc,dims)

implicit none
real(rp), intent(inout), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6) :: PFM_phi
real(rp),intent(in)::    deltan(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
real(rp),intent(in)::     nabs_surf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::i_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::i_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::j_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
integer,intent(in)::k_IP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6)
real(rp),intent(in)::     WP1(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
real(rp),intent(in)::     WP2(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7)
real(rp),intent(in)::     dzc(0:)
integer ,intent(in), dimension(3)     :: dims
integer::k,j,i,ijk_start_x,ijk_start_y,ijk_start_z
real(rp):: starting_point_j,end_point_j,starting_point_k,end_point_k,phi_wall
real(rp),intent(in)::time
real(rp):: time_stari_visc,time_star_inert
real(rp)::Wetting_rad, phi_m,delta_y,delta_z
logical:: confirmation
real(rp),dimension(1:dims(2)):: j_min,j_max,k_min,k_max
integer::min_processor,max_processor
real(rp)::j_min_all,j_max_all,dn_1,dn_2,alpha
integer :: error
!
ijk_start_x = ijk_start(1)
ijk_start_y = ijk_start(2)
ijk_start_z = ijk_start(3)
!
j_min(:) =  10000.0_rp
j_max(:) = -10000.0_rp
k_min(:) =  10000.0_rp
k_max(:) = -10000.0_rp
Wetting_rad = 0.0_rp
j_min_all =  10000.0_rp
j_max_all = -10000.0_rp
min_processor = 0
max_processor = 0
confirmation = .false.
delta_y = 0.0_rp
delta_z = 0.0_rp
j_min(:) =  10000.0_rp
j_max(:) = -10000.0_rp
k_min(:) =  10000.0_rp
k_max(:) = -10000.0_rp
Wetting_rad = 0.0_rp
j_min_all = 0.0_rp
j_max_all = 0.0_rp
min_processor = 0
max_processor = 0

  i=int(n(1)/2)
  do j=1,n(2)
    do k=1,n(3)
      if  (nabs_surf(i,j,k).gt.small) then
          call interpolation_mirror(PFM_phi,i,j,k, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2, phi_m)



          if ((phi_m).gt.0.5) then
           j_min(ijk_start_y+1) = (j+ijk_start_y)*1.0_rp
           k_min(ijk_start_y+1) = k*1.0_rp
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
          call interpolation_mirror(PFM_phi,i,j,k, &
                                i_IP1,j_IP1,k_IP1, &
                                i_IP2,j_IP2,k_IP2, &
                                WP1,WP2, phi_m)
          if ((phi_m).gt.0.5) then
           j_max(ijk_start_y+1) = (j+ijk_start_y)*1.0_rp
           k_max(ijk_start_y+1) = k*1.0_rp
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
  Wetting_rad = 0.5*sqrt(delta_y**2 + delta_z**2)
  if (myid.eq.0) then
    open(1329,file='data/spreading_radius.txt',position='append')
    write(1329,'(2E16.8)' ) time, Wetting_rad
    close(1329)
  endif

return
end subroutine Wetting_radius
#endif
end module mod_IBM