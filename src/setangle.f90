module mod_setangle
#if defined(_USE_IBM)
  use mpi
  use mod_IBM, only : interpolation_mirror,interpolation_dphi
  use mod_contactangle, only: getTheta
  use mod_param, only: small,dynamic_angle,theta1,theta2,zmin_ibm
  use mod_types
  !@cuf use cudafor
  !@cuf use mod_common_mpi, only: mydev
  !
  implicit none
  real(rp), parameter :: pi = acos(-1._rp)
  !
  private
  public  :: setangle
  !
  contains
  !
  subroutine setangle(n,dl,vof,nabs_surf,nx_surf,ny_surf,nz_surf,deltan, &
                      i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2, &
                      WP1,WP2,dir,vel_x,vel_y,vel_z,zc,dzc)
!
! Routine for imposing static/dynamic contact angle using IBM's ghost cell method for imposing BCs
! Contact line velocity calculation adapted from Patel et al. 2017 Chem. Eng. Sci. 
!
   implicit none
   integer, dimension(3), intent(in) :: n
   real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),     intent(inout) :: vof
   real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),     intent(in)    :: nabs_surf,nx_surf,ny_surf,nz_surf,deltan
   integer,  dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),     intent(in)    :: i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2
   real(rp), dimension(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6,1:7), intent(in)    :: WP1,WP2
   real(rp), dimension(0:,0:,0:), intent(in) :: vel_x, vel_y, vel_z
   real(rp), dimension(3) , intent(in) :: dl
   real(rp), dimension(0:), intent(in) :: zc,dzc
   real(rp), dimension(3) :: tang, vel_cl_xyz
   integer, intent(in) :: dir
   integer :: i,j,k,n1,n2,n3
   real(rp):: temp, chem_m, dphidx, dphidy, dphidz, dPhidN, dPhidTAbs2, cosTheta, d, thetad, vel_cl, intdir, z, dz, z_bottom
   !@cuf integer :: istat
   !@cuf attributes(managed) :: vof,vel_x,vel_y,vel_z,zc
   !@cuf attributes(managed) :: nabs_surf,nx_surf,ny_surf,nz_surf,deltan,i_IP1,j_IP1,k_IP1,i_IP2,j_IP2,k_IP2,WP1,WP2
   !
   !@cuf istat = cudaMemPrefetchAsync(vof, size(vof), cudaCpuDeviceId, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_x, size(vel_x), cudaCpuDeviceId, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_y, size(vel_y), cudaCpuDeviceId, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_z, size(vel_z), cudaCpuDeviceId, 0)
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
   !
   n1=n(1)
   n2=n(2)
   n3=n(3)
   !
   d = sqrt(dl(1)**2 + dl(2)**2 + dl(3)**2)
   z_bottom = zmin_ibm
   !
    do k=1,n3
      z = zc(k)
      dz = dzc(k)
       do j=1,n2
         do i=1,n1
  
           chem_m = 0._rp
           vel_cl = 0._rp
           vel_cl_xyz(:) = 0._rp
  
           if (nabs_surf(i,j,k).gt.small) then ! Interface cell (ghost points)
  
             call interpolation_mirror(vof,i,j,k, &
                                       i_IP1,j_IP1,k_IP1, &
                                       i_IP2,j_IP2,k_IP2, &
                                       WP1,WP2,chem_m)   ! Compute VoF marker at mirror point associated with ghost point at i,j,k
  
             call interpolation_dphi(vof,i,j,k, &
                                     i_IP1,j_IP1,k_IP1, &
                                     i_IP2,j_IP2,k_IP2, &
                                     WP1,WP2, &
                                     dphidx,dphidy,dphidz) ! Compute VoF marker gradient at mirror point associated with ghost point at i,j,k
  
             dPhidN = dphidx*nx_surf(i,j,k) + dphidy*ny_surf(i,j,k) + dphidz*nz_surf(i,j,k)  ! Project gradient onto surface normal
  
             dPhidTAbs2 = (dphidx - dPhidN*nx_surf(i,j,k))**2 + &
                          (dphidy - dPhidN*ny_surf(i,j,k))**2 + &
                          (dphidz - dPhidN*nz_surf(i,j,k))**2   ! Find absolute value squared of surface tangent
  
             if (dynamic_angle) then
              tang(1) = (dphidx - dPhidN*nx_surf(i,j,k))/sqrt(dPhidTAbs2) ! Components of unit surface tangent
              tang(2) = (dphidy - dPhidN*ny_surf(i,j,k))/sqrt(dPhidTAbs2)
              tang(3) = (dphidz - dPhidN*nz_surf(i,j,k))/sqrt(dPhidTAbs2)
  
              temp = (0.5_rp*(vel_x(i,j,k) + vel_x(i-1,j,k)))*tang(1) + &
                     (0.5_rp*(vel_y(i,j,k) + vel_y(i,j-1,k)))*tang(2) + &
                     (0.5_rp*(vel_z(i,j,k) + vel_z(i,j,k-1)))*tang(3)   ! Projected velocity along surface tangent
  
              vel_cl_xyz(1) = temp*tang(1)/sqrt(dPhidTAbs2)
              vel_cl_xyz(2) = temp*tang(2)/sqrt(dPhidTAbs2)
              vel_cl_xyz(3) = temp*tang(3)/sqrt(dPhidTAbs2)
  
              vel_cl = sqrt(vel_cl_xyz(1)**2 + vel_cl_xyz(2)**2 + vel_cl_xyz(3)**2)
  
              intdir = sign(1._rp,( vel_cl_xyz(1)*dphidx/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) + &
                                     vel_cl_xyz(2)*dphidy/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) + &
                                     vel_cl_xyz(3)*dphidz/(sqrt(dphidx**2 + dphidy**2 + dphidz**2)) ) ) !Interface movement direction
  
              vel_cl = intdir*vel_cl ! Contact line velocity
  
              if( (z-z_bottom).le.dz) call getTheta(theta1, dl(dir), thetad, vel_cl)  ! Get dynamic contact angle
              if( (z-z_bottom).gt.dz) call getTheta(theta2, dl(dir), thetad, vel_cl)
             else
              if( (z-z_bottom).le.dz) thetad = theta1
              if( (z-z_bottom).gt.dz) thetad = theta2
             endif
  
             dPhidN = cos(thetad*pi/180._rp)*sqrt((1._rp/(1._rp - cos(thetad*pi/180._rp)**2))*dPhidTAbs2) ! Adjust the gradient
  
             chem_m = chem_m + dPhidN*(d + deltan(i,j,k))    ! Updated value for VoF marker at ghost point
  
             if (chem_m == chem_m) vof(i,j,k) = chem_m ! NaN Check
  
           endif
  
         enddo
       enddo
    enddo
   !
   !@cuf istat = cudaMemPrefetchAsync(vof, size(vof), mydev, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_x, size(vel_x), mydev, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_y, size(vel_y), mydev, 0)
   !@cuf istat = cudaMemPrefetchAsync(vel_z, size(vel_z), mydev, 0)
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
   !
   return
  end subroutine setangle
#endif
end module mod_setangle