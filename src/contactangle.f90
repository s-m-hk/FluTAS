module mod_contactangle
#if defined(_USE_IBM) && defined(_USE_CONTACTANGLE_DYNAMIC)
use mod_param, only: mu1, mu2, sigma, pi
use mod_types
!
implicit none
real(rp), parameter :: gridCorrectionK = 1.e-9_rp
private
public getTheta
!
contains
#if !defined(_OPENACC)
!
! Cox hydrodynamic model for dynamic contact angle, see eq. 11 of Legendre and Maglio 2015 Computers & Fluids
!
subroutine getTheta(thetaIn, Delta, thetaOut, velInterf)
  implicit none
  !$acc routine(getTheta) seq
  real(rp), intent(in) :: thetaIn, Delta, velInterf
  real(rp), intent(out) :: thetaOut
  real(rp) :: Ca, NSconst
  
  Ca = sqrt(mu1*mu2)*velInterf/sigma

  NSconst = gStar(thetaIn) + Ca*log((Delta/2.0_rp)/gridCorrectionK)
  if (NSConst .le. 0._rp) then ! No solution exist in this case
    ! thetaOut = 1.0d-8; ! print*, "No solution for theta exists: theta set to 0"
    thetaOut = 30._rp*pi/180._rp;
  else if (NSConst .ge. gStar(pi)) then
    ! thetaOut = pi - 1.0d-8; ! print*, "No solution for theta exists: theta set to pi"
    thetaOut = 150._rp*pi/180._rp;
  else
    ! call NewtonSolver(gStar, fStarInv, thetaIn, NSconst, theta)
    call NewtonSolver(gStar, fStarInv, 0.5_rp*pi, NSconst, thetaOut)
            ! 0.5*pi is a good initial condition if K differs a lot from dx
    ! thetaOut = acos(cos(theta_stat*angleToRadian) + 5.63*Ca*log(0.02*2.0/dz(ks)))

    ! Set boundaries
    if (thetaOut < 30.0_rp*pi/180.0_rp) then
        thetaOut =  30.0_rp*pi/180.0_rp
    elseif (thetaOut > 150.0_rp*pi/180.0_rp) then
        thetaOut = 150.0_rp*pi/180.0_rp
    endif
  endif
end subroutine getTheta


subroutine NewtonSolver(F, Fp, x0, Const, xVar) !Numerical solver using Newton's method

! Estimate the zero of f(x) using Newton's method.
! Input:
!   F:  the function to find root of F = Const
!   Fp: function returning the derivative f'
!   x0: the initial guess
! Returns:
!   the estimate x satisfying F(x) = Const (assumes Newton converged!)

implicit none
real(rp), intent(in) :: x0, Const
real(rp), external :: F, Fp
real(rp), intent(out) :: xVar

real(rp) :: deltax, fx, fxprime, iter, tol = 1.e-6_rp
integer :: k, maxiter = 40

! Initial guess
xVar = x0

! Newton iteration to find a zero of f(x)
!$acc seq
do k = 1,maxiter
    ! Calculate function and its derivative:
    fx = F(xVar) - Const
    fxprime = Fp(xVar)

    if (abs(fx) < tol) then
        ! print*,"k =", k
        exit  ! jump out of do loop
    endif

    ! Compute Newton increment x:
    deltax = fx/fxprime

    ! Update x:
    xVar = xVar - deltax
enddo

if (k .ge. maxiter) then
    ! Might not have converged
    fx = F(xVar) - Const
    if (abs(fx) > tol) then
        print *, '****** Warning: Newton solver has not yet converged ******', abs(fx), Const
        xVar = x0 ! return initial value
    endif
endif

end subroutine NewtonSolver

real(rp) function fStarInv(phi)
  implicit none
  real(rp), intent(in) :: phi

  real(rp) :: sqrtq, sin2phi,q
  if (phi == 0._rp) then
    fStarInv = 0._rp
    return
  endif

  ! The angle is measured in fluid 2
  q = mu1/mu2; sqrtq = sqrt(q); sin2phi = sin(phi)**2

  fStarInv = (sqrtq* ( phi**2          - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
             1/sqrtq*((pi - phi)**2 - sin2phi)*( phi          - sin(phi)*cos(phi)))/ &
             (2._rp*sin(phi)*((2._rp - q - 1._rp/q)*sin2phi + q*phi**2 + 1._rp/q*(pi - phi)**2 + 2._rp*phi*(pi - phi) ))
end function fStarInv

real(rp) function gStar(theta)
  ! Approximation of the definite integral of fStarInv with the
  ! composite trapezoidal rule, using n subintervals
  ! See http://en.wikipedia.org/wiki/Trapezoid_rule

  implicit none
  real(rp), intent(in) :: theta

  ! [a, b] : the interval of integration
  real(rp) :: a = 0.0_rp, b
  real(rp) :: summation, invN
  ! n : number of subintervals used
  integer :: n = 500, k

  invN = 1._rp/n
  b = theta
  summation = 0._rp
  !$acc seq
  do k = 1, n-1
      summation = summation + fStarInv(a + (b-a)*k*invN)
  end do

  gStar = (b-a) * (0.5_rp*fStarInv(a) + 0.5_rp*fStarInv(b) + summation) * invN
end function gStar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GPU accelerated version!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getTheta(thetaIn, Delta, thetaOut, velInterf)
  implicit none
  real(rp), intent(in) :: thetaIn, Delta, velInterf
  real(rp), intent(out) :: thetaOut
  real(rp) :: Ca, NSconst
  real(rp) :: gStar, gStar_pi, fStarInv
  real(rp) :: summation, phi, q, sin2phi,a = 0._rp
  integer :: k, n = 500
  
  Ca = sqrt(mu1*mu2)*velInterf/sigma

  !
  q = mu1/mu2
  summation = 0._rp
  fStarInv = 0._rp
  !
  !$acc kernels
  do k = 1, n-1
   phi = a + (thetaIn-a)*k*1._rp/n
   sin2phi = sin(phi)**2
   if (phi == 0._rp) exit
   !
   fStarInv = (sqrt(q)*(phi**2 - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
              (1._rp/sqrt(q))*((pi - phi)**2 - sin2phi)*(phi - sin(phi)*cos(phi)))/ &
              (2._rp*sin(phi)*((2._rp - q - 1._rp/q)*sin2phi + q*phi**2 + 1._rp/q*(pi - phi)**2 + 2._rp*phi*(pi - phi) ))
   !
   summation = summation + fStarInv
  end do
  !$acc end kernels
  gStar = (thetaIn-a) * (0.5_rp*((sqrt(q)*(a**2 - sin(a)**2)*((pi - a) + sin(a)*cos(a)) + &
                                 (1._rp/sqrt(q))*((pi - a)**2 - sin(a)**2)*(a - sin(a)*cos(a)))/ &
                                 (2._rp*sin(a)*((2._rp - q - 1._rp/q)*sin(a)**2 + q*a**2 + 1._rp/q*(pi - a)**2 + 2._rp*a*(pi - a) ))) &
                         + 0.5_rp*((sqrt(q)*(thetaIn**2 - sin(thetaIn)**2)*((pi - thetaIn) + sin(thetaIn)*cos(thetaIn)) + &
                                   (1._rp/sqrt(q))*((pi - thetaIn)**2 - sin(thetaIn)**2)*(thetaIn - sin(thetaIn)*cos(thetaIn)))/ &
                                   (2._rp*sin(thetaIn)*((2._rp - q - 1._rp/q)*sin(thetaIn)**2 + q*thetaIn**2 + 1._rp/q*(pi - thetaIn)**2 + 2._rp*thetaIn*(pi - thetaIn) ))) &
                         + summation) * 1._rp/n
  !
  !$acc kernels
  do k = 1, n-1
   phi = a + (pi-a)*k*1._rp/n
   sin2phi = sin(phi)**2
   if (phi == 0._rp) exit
   !
   fStarInv = (sqrt(q)*(phi**2 - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
              (1._rp/sqrt(q))*((pi - phi)**2 - sin2phi)*(phi - sin(phi)*cos(phi)))/ &
              (2._rp*sin(phi)*((2._rp - q - 1._rp/q)*sin2phi + q*phi**2 + 1._rp/q*(pi - phi)**2 + 2._rp*phi*(pi - phi) ))
   !
   summation = summation + fStarInv
  end do
  !$acc end kernels
  gStar_pi = (pi-a) * (0.5_rp*((sqrt(q)*(a**2 - sin(a)**2)*((pi - a) + sin(a)*cos(a)) + &
                                 (1._rp/sqrt(q))*((pi - a)**2 - sin(a)**2)*(a - sin(a)*cos(a)))/ &
                                 (2._rp*sin(a)*((2._rp - q - 1._rp/q)*sin(a)**2 + q*a**2 + 1._rp/q*(pi - a)**2 + 2._rp*a*(pi - a) ))) &
                       + 0.5_rp*((sqrt(q)*(pi**2 - sin(pi)**2)*((pi - pi) + sin(pi)*cos(pi)) + &
                                   (1._rp/sqrt(q))*((pi - pi)**2 - sin(pi)**2)*(pi - sin(pi)*cos(pi)))/ &
                                   (2._rp*sin(pi)*((2._rp - q - 1._rp/q)*sin(pi)**2 + q*pi**2 + 1._rp/q*(pi - pi)**2 + 2._rp*pi*(pi - pi) ))) &
                       + summation) * 1._rp/n
  !

  NSconst = gStar + Ca*log((Delta/2.0_rp)/gridCorrectionK)
  if (NSConst .le. 0._rp) then ! No solution exist in this case
    ! thetaOut = 1.0d-8; ! print*, "No solution for theta exists: theta set to 0"
    thetaOut = 30._rp*pi/180._rp;
  else if (NSConst .ge. gStar_pi) then
    ! thetaOut = pi - 1.0d-8; ! print*, "No solution for theta exists: theta set to pi"
    thetaOut = 150._rp*pi/180._rp;
  else
    ! call NewtonSolver(gStar, fStarInv, thetaIn, NSconst, theta)
    call NewtonSolver(0.5_rp*pi, NSconst, thetaOut)
            ! 0.5*pi is a good initial condition if K differs a lot from dx
    ! thetaOut = acos(cos(theta_stat*angleToRadian) + 5.63*Ca*log(0.02*2.0/dz(ks)))

    ! Set boundaries
    if (thetaOut < 30.0_rp*pi/180.0_rp) then
        thetaOut =  30.0_rp*pi/180.0_rp
    elseif (thetaOut > 150.0_rp*pi/180.0_rp) then
        thetaOut = 150.0_rp*pi/180.0_rp
    endif
  endif
end subroutine getTheta
!
subroutine NewtonSolver(x0, Const, xVar) !Numerical solver using Newton's method

! Estimate the zero of f(x) using Newton's method.
! Input:
!   F:  the function to find root of F = Const
!   Fp: function returning the derivative f'
!   x0: the initial guess
! Returns:
!   the estimate x satisfying F(x) = Const (assumes Newton converged!)

implicit none
real(rp), intent(in) :: x0, Const
real(rp), intent(out) :: xVar
real(rp) :: gStar, fStarInv
real(rp) :: phi, sin2phi, deltax, fx, fxprime, iter, summation, q, a = 0._rp, tol = 1.e-6_rp
integer :: k, kk, maxiter = 40, n = 500

q = mu1/mu2

! Initial guess
xVar = x0

! Newton iteration to find a zero of f(x)
!$acc kernels
do k = 1,maxiter
  ! Calculate function and its derivative:
  summation = 0._rp
  fStarInv = 0._rp
  !
  do kk = 1, n-1
   phi = a + (xVar-a)*kk*1._rp/n
   sin2phi = sin(phi)**2
   if (phi == 0._rp) exit
   !
   fStarInv = (sqrt(q)*(phi**2 - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
              (1._rp/sqrt(q))*((pi - phi)**2 - sin2phi)*(phi - sin(phi)*cos(phi)))/ &
              (2._rp*sin(phi)*((2._rp - q - 1._rp/q)*sin2phi + q*phi**2 + 1._rp/q*(pi - phi)**2 + 2._rp*phi*(pi - phi) ))
   !
   summation = summation + fStarInv
  end do
  gStar = (xVar-a) * (0.5_rp*((sqrt(q)*(a**2 - sin(a)**2)*((pi - a) + sin(a)*cos(a)) + &
                                 (1._rp/sqrt(q))*((pi - a)**2 - sin(a)**2)*(a - sin(a)*cos(a)))/ &
                                 (2._rp*sin(a)*((2._rp - q - 1._rp/q)*sin(a)**2 + q*a**2 + 1._rp/q*(pi - a)**2 + 2._rp*a*(pi - a) ))) &
                      + 0.5_rp*((sqrt(q)*(xVar**2 - sin(xVar)**2)*((pi - xVar) + sin(xVar)*cos(xVar)) + &
                                   (1._rp/sqrt(q))*((pi - xVar)**2 - sin(xVar)**2)*(xVar - sin(xVar)*cos(xVar)))/ &
                                   (2._rp*sin(xVar)*((2._rp - q - 1._rp/q)*sin(xVar)**2 + q*xVar**2 + 1._rp/q*(pi - xVar)**2 + 2._rp*xVar*(pi - xVar) ))) &
                      + summation) * 1._rp/n
  !
  fx = gStar - Const
  fxprime = (sqrt(q)*(xVar**2 - sin2phi)*((pi - xVar) + sin(xVar)*cos(xVar)) + &
              (1._rp/sqrt(q))*((pi - xVar)**2 - sin2phi)*(xVar - sin(xVar)*cos(xVar)))/ &
              (2._rp*sin(xVar)*((2._rp - q - 1._rp/q)*sin2phi + q*xVar**2 + 1._rp/q*(pi - xVar)**2 + 2._rp*xVar*(pi - xVar) ))

  if (abs(fx) < tol) then
      ! print*,"kk =", kk
      exit  ! jump out of do loop
  endif

  ! Compute Newton increment x:
  deltax = fx/fxprime

  ! Update x:
  xVar = xVar - deltax
enddo
!$acc end kernels

if (k .ge. maxiter) then
    ! Might not have converged
    summation = 0._rp
    fStarInv = 0._rp
    !
    do kk = 1, n-1
     phi = a + (xVar-a)*kk*1._rp/n
     sin2phi = sin(phi)**2
     if (phi == 0._rp) exit
     !
     fStarInv = (sqrt(q)*(phi**2 - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
                (1._rp/sqrt(q))*((pi - phi)**2 - sin2phi)*(phi - sin(phi)*cos(phi)))/ &
                (2._rp*sin(phi)*((2._rp - q - 1._rp/q)*sin2phi + q*phi**2 + 1._rp/q*(pi - phi)**2 + 2._rp*phi*(pi - phi) ))
     !
     summation = summation + fStarInv
    end do
    gStar = (xVar-a) * (0.5_rp*((sqrt(q)*(a**2 - sin(a)**2)*((pi - a) + sin(a)*cos(a)) + &
                                 (1._rp/sqrt(q))*((pi - a)**2 - sin(a)**2)*(a - sin(a)*cos(a)))/ &
                                 (2._rp*sin(a)*((2._rp - q - 1._rp/q)*sin(a)**2 + q*a**2 + 1._rp/q*(pi - a)**2 + 2._rp*a*(pi - a) ))) &
                        + 0.5_rp*((sqrt(q)*(xVar**2 - sin(xVar)**2)*((pi - xVar) + sin(xVar)*cos(xVar)) + &
                                   (1._rp/sqrt(q))*((pi - xVar)**2 - sin(xVar)**2)*(xVar - sin(xVar)*cos(xVar)))/ &
                                   (2._rp*sin(xVar)*((2._rp - q - 1._rp/q)*sin(xVar)**2 + q*xVar**2 + 1._rp/q*(pi - xVar)**2 + 2._rp*xVar*(pi - xVar) ))) &
                        + summation) * 1._rp/n
    !
    fx = gStar - Const
    if (abs(fx) > tol) then
        print *, '****** Warning: Newton solver has not yet converged ******', abs(fx), Const
        xVar = x0 ! return initial value
    endif
endif

end subroutine NewtonSolver
#endif
#endif
end module mod_contactangle