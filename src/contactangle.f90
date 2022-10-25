module mod_contactangle
#if defined(_USE_IBM)
use mod_param, only: mu1, mu2, sigma
use mod_types

implicit none
real(rp), parameter :: pi = acos(-1._rp)
real(rp), parameter :: gridCorrectionK = 1.0e-9_rp
private
public getTheta

contains
!
! Cox hydrodynamic model for dynamic contact angle, see eq. 11 of Legendre and Maglio 2015 Computers & Fluids
!
subroutine getTheta(thetaIn, Delta, thetaOut, velInterf)
  implicit none
  real(rp), intent(in) :: thetaIn, Delta, velInterf
  real(rp), intent(out) :: thetaOut
  real(rp) :: Ca, NSconst
  
  Ca = sqrt(mu1*mu2)*velInterf/sigma

  NSconst = gStar(thetaIn) + Ca*log((Delta/2.0_rp)/gridCorrectionK)
  if (NSConst .le. 0.0_rp) then ! No solution exist in this case
    ! thetaOut = 1.0d-8; ! print*, "No solution for theta exists: theta set to 0"
    thetaOut = 30.0_rp*pi/180.0_rp;
  else if (NSConst .ge. gStar(pi)) then
    ! thetaOut = pi - 1.0d-8; ! print*, "No solution for theta exists: theta set to pi"
    thetaOut = 150.0_rp*pi/180.0_rp;
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

real(rp) :: deltax, fx, fxprime, iter, tol = 1.0e-5_rp
integer :: k, maxiter = 20


! Initial guess
xVar = x0

! Newton iteration to find a zero of f(x)
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
  if (phi == 0.0_rp) then
    fStarInv = 0.0_rp
    return
  endif

  ! The angle is measured in fluid 2
  q = mu1/mu2; sqrtq = sqrt(q); sin2phi = sin(phi)**2

  fStarInv = (sqrtq* ( phi**2          - sin2phi)*((pi - phi) + sin(phi)*cos(phi)) + &
             1/sqrtq*((pi - phi)**2 - sin2phi)*( phi          - sin(phi)*cos(phi)))/ &
             (2*sin(phi)*((2 - q - 1/q)*sin2phi + q*phi**2 + 1/q*(pi - phi)**2 + 2*phi*(pi - phi) ))
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
  integer :: n = 400, k

  invN = 1.0_rp/n
  b = theta
  summation = 0.0_rp
  do k = 1, n-1
      summation = summation + fStarInv(a + (b-a)*k*invN)
  end do
  gStar = (b-a) * (0.5_rp*fStarInv(a) + 0.5_rp*fStarInv(b) + summation) * invN
end function gStar
#endif
end module mod_contactangle