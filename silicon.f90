module silicon

real(8), public:: freq
real(8), dimension(3), parameter::g=(/9.0e31_8, 2.8e32_8, 8.2e31_8/)
real(8), dimension(3), parameter::gammam=(/1.6e15_8, 1.05e15_8, 5.4e14_8/)
real(8), dimension(3), parameter::wm=(/8.32e15_8, 6.37e15_8, 5.33e15_8/)

contains

function oscillators(x)

real(4):: x, oscillators

	do i=1,size(g)
	oscillators=oscillators+g(i)*gammam(i)*x/((wm(i)**2-x**2)**2+gammam(i)**2*x**2)
	end do

end function

!-----------------------------------------------------------------------

function KramersKr(x)

real(4):: x, KramersKr

KramersKr=x/(x**2+freq**2)*oscillators(x)

end function


!---------------------------------------------------------------------

subroutine monte_carlo ( func, a, b, n, result )

!*****************************************************************************80
!
!! MONTE_CARLO estimates the integral of a function by Monte Carlo.
!
!  Modified:
!
!    18 December 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, name of external function to be
!    integrated.  This name must be in an external statement in the calling
!    program.  FUNC must be a function of one real argument.  The value
!    of the argument to FUNC is the variable of integration
!    which ranges from A to B.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.
!
!    Input, integer ( kind = 4 ) N, the number of points to use.
!
!    Output, real ( kind = 8 ) RESULT, the computed value of the integral.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ), external :: func
  integer ( kind = 4 ) i
  real ( kind = 8 ) result
  real ( kind = 8 ) x(n)

  call random_number ( harvest = x(1:n) )

  x(1:n) = a + ( b - a ) * x(1:n)

  result = 0.0D+00
  do i = 1, n
    result = result + func ( x(i) )
  end do

  result = result / real ( n, kind = 8 )

  return
end subroutine


end module
