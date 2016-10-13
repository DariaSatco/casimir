module dopcasimir

    real(8), public :: w,dist
	! w - Matsubara frequency
	! dist - the distance between interacting plates
    real(8), public :: eps(3)
	! eps - dielectric permettivity vector
        ! eps(1),eps(2) describes materials
	! eps(3) describes gap permettivity
    real(8), parameter :: pi=3.1415926
    real(8), parameter :: eV=1.519e15  ! rad/s
    real(8), parameter :: kb=8.617e-5  ! eV/K kb=1.38e-23 J/K  Boltzman constant
    real(8), parameter :: h=6.585e-16  ! eV*s h=1.054e-34J*s  Plank constant
    real(8), parameter :: c=299792458. ! c - velocity of light


contains


	function eq(a,x,y)
	! function needed to integrate in formula for Casimir free energy 
	! x=dist*r, y=w*dist/c change the variables to dimensioness
        ! the constants in front of the integral: rdr --> 1/dist**2 xdx;  dw --> c/dist

	real(8):: a,x,q,y,rTMn(2),rTEn(2),eq

		rTMn=rTM(x/a,y*c/a)
!calculate the TM Frenel coefficient for the concrete dist as a function of x,y
		rTEn=rTE(x/a,y*c/a)
!calculate the TE Frenel coefficient for the concrete dist as a function of x,y
			
		q=sqrt(x**2+y**2)
		eq=x*(log(1._8-rTMn(1)*rTMn(1)*exp(-2.*q))+log(1._8-rTEn(1)*rTEn(1)*exp(-2.*q)))

					
	end function eq


!-----------------------------------------------------------------------------------------------
	function zerotempint(arg1,arg2)
	! function - argument of 2-dimendional integral by frequencies and k-frequencies at T=0 
        ! limits [1, Inf) --> (0,1]
        ! change variables xdx -->  1/(arg1**3)darg1;  dy --> 1/(arg2**2)darg2

	real(8):: arg1,arg2,zerotempint

	if (arg1.eq.(0._8)) then
		
		zerotempint=0._8

		else

		zerotempint=(1._8/arg1**2)*(1._8/arg2**2)*eq(dist, 1._8/arg1, 1._8/arg2)

	endif

	end function

!-----------------------------------------------------------------------------------------------

	
	function to_int(arg)
	!integrable function for nonzero temperatures - function of one argument - x <--> k-frequency
        !interval [0,1]

	real(8):: arg, to_int

	to_int = eq(dist, arg, w*dist/c)

	end function

!------------------------------------------------------------------------------------------------

	function to_int1(arg)
	!integrable function for nonzero temperatures
        !make change of interval for integration [1, Inf) --> (0,1]

	real(8):: arg, to_int1
	
	if (arg.eq.(0._8)) then
		to_int1=0._8
	else

		to_int1 = (1._8/arg**2)*eq(dist, 1._8/arg, w*dist/c)

	endif

	end function

!--------------------------------------------------------------------------------------------

	function rTM(kf,freq)
	!the Frenel reflection coeffitient for TM-field

	  real(8):: kf,freq, rTM(2),q,k(2)
	  integer:: i


		q=sqrt(kf**2+(freq/c)**2)
		do i=1,2
			k(i)=sqrt(kf**2+eps(i)*(freq/c)**2)
			rTM(i)=(eps(i)*q-k(i))/(eps(i)*q+k(i))
		end do
	  
	end function rTM



!--------------------------------------------------------------------------------

	function rTE(kf,freq)
	! the Frenel reflection coeffitien for TE-field

	  real(8):: kf,freq,rTE(2),q,k(2)
	  integer:: i

	  
		q=sqrt(kf**2+(freq/c)**2)
		do i=1,2
			k(i)=sqrt(kf**2+eps(i)*(freq/c)**2)
			rTE(i)=(q-k(i))/(q+k(i))
		end do
	  
	end function rTE

!-------------------------------------------------------------------------------------------------	 
	 
	subroutine trapzd(func,lowl,upl,res)
	! approximates the integral with the trapezoidal method
        
        real(8):: func
	real(8):: res
	real(8):: lowl,upl
	real(8):: s,it
	integer:: i,j,n
	
	s=0._8
	res=0.5*abs(upl-lowl)*(func(lowl)+func(upl))
	i=0

		do while (abs(res-s).gt.(1e-3*abs(res)))
			s=res
			res=0._8
			i=i+1
			n=2**i
			it=abs(upl-lowl)/n
				
				do j=0,n-1
					res=res+0.5*it*(func(lowl+j*it)+func(lowl+(j+1)*it))
				enddo

		enddo

	!print *, 'step number=', 2**i
	return

	end subroutine

!------------------------------------------------------------------------------------------------

	subroutine trapzd2d(func,lowl1,upl1,lowl2,upl2,res,stepn)
	! approximates the integral with the trapezoidal method

	real(8):: func
	real(8):: res
	real(8), allocatable:: q(:,:)
	real(8):: lowl1,upl1,lowl2,upl2
	real(8):: s,it1, it2
	integer:: i,j,n,stepn
	
	s=0.
	res=0.25*abs(upl1-lowl1)*abs(upl2-lowl2)*(func(lowl1,lowl2)+func(upl1,upl2))
	i=0

		do while (abs(res-s).gt.(1e-2*abs(res)))
			s=res
			res=0.
			i=i+1
			n=2**i
			it1=abs(upl1-lowl1)/n	
			it2=abs(upl2-lowl2)/n
			allocate(q(0:n-1,0:n-1))
				
				do j=0,n-1

				do k=0,n-1

				if ((j==0).or.(j==n-1)) then
					if ((k==0).or.(k==n-1)) then
					q(j,k)=0.25
						else
					q(i,j)=0.5
					endif
				else
					if ((k==0).or.(k==n-1)) then
					q(j,k)=0.5
						else
					q(j,k)=1.
					endif
				endif
				
					res=res+q(j,k)*(func(lowl1+j*it1,lowl2+k*it2))
				enddo

				enddo
			deallocate(q)

		res=res*it1*it2

		enddo
	
	stepn=n

	return

	end subroutine
	
!-----------------------------------------------------------------------------------------


subroutine cspint ( ntab, xtab, ftab, a, b, y, e, work, result )

!*****************************************************************************
!
!! CSPINT estimates the integral of a tabulated function.
!
!  Discussion:
!
!    The routine is given the value of a function F(X) at a set of 
!    nodes XTAB, and estimates
!
!      Integral ( A <= X <= B ) F(X) DX
!
!    by computing the cubic natural spline S(X) that interpolates
!    F(X) at the nodes, and then computing
!
!      Integral ( A <= X <= B ) S(X) DX
!
!    exactly.
!
!    Other output from the program includes the definite integral
!    from X(1) to X(I) of S(X), and the coefficients necessary for
!    the user to evaluate the spline S(X) at any point.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Carl DeBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was evaluated.  The XTAB's must be distinct and
!    in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated values of
!    the function, FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) A, lower limit of integration.
!
!    Input, real ( kind = 8 ) B, upper limit of integration.
!
!    Output, real ( kind = 8 ) Y(3,NTAB), will contain the coefficients
!    of the interpolating natural spline over each subinterval.
!    For XTAB(I) <= X <= XTAB(I+1),
!      S(X) = FTAB(I) + Y(1,I)*(X-XTAB(I))
!                   + Y(2,I)*(X-XTAB(I))**2
!                   + Y(3,I)*(X-XTAB(I))**3
!
!    Output, real ( kind = 8 ) E(NTAB), E(I) = the definite integral from
!    XTAB(1) to XTAB(I) of S(X).
!
!    Workspace, real ( kind = 8 ) WORK(NTAB).
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) e(ntab)
  real ( kind = 8 ) ftab(ntab)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) r
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) term
  real ( kind = 8 ) u
  real ( kind = 8 ) work(ntab)
  real ( kind = 8 ) xtab(ntab)
  real ( kind = 8 ) y(3,ntab)

  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CSPINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB must be at least 3, but input NTAB = ', ntab
    stop 1
  end if
 
  do i = 1, ntab-1
 
    if ( xtab(i+1) <= xtab(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CSPINT - Fatal error!'
      write ( *, '(a)' ) '  Nodes not in strict increasing order.'
      write ( *, '(a,i8)' ) '  XTAB(I) <= XTAB(I-1) for I=',i
      write ( *, '(a,g14.6)' ) '  XTAB(I) = ',xtab(i)
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ',xtab(i-1)
      stop 1
    end if
 
  end do
 
  s = 0.0D+00
  do i = 1, ntab-1
    r = ( ftab(i+1) - ftab(i) ) / ( xtab(i+1) - xtab(i) )
    y(2,i) = r - s
    s = r
  end do
 
  result = 0.0D+00
  s = 0.0D+00
  r = 0.0D+00
  y(2,1) = 0.0D+00
  y(2,ntab) = 0.0D+00
 
  do i = 2, ntab-1
    y(2,i) = y(2,i) + r * y(2,i-1)
    work(i) = 2.0D+00 * ( xtab(i-1) - xtab(i+1) ) - r * s
    s = xtab(i+1) - xtab(i)
    r = s / work(i)
  end do
 
  do j = 2, ntab-1
    i = ntab+1-j
    y(2,i) = ( ( xtab(i+1) - xtab(i) ) * y(2,i+1) - y(2,i) ) / work(i)
  end do
 
  do i = 1, ntab-1
    s = xtab(i+1) - xtab(i)
    r = y(2,i+1) - y(2,i)
    y(3,i) = r / s
    y(2,i) = 3.0D+00 * y(2,i)
    y(1,i) = ( ftab(i+1) - ftab(i) ) / s - ( y(2,i) + r ) * s
  end do
 
  e(1) = 0.0D+00
  do i = 1, ntab-1
    s = xtab(i+1)-xtab(i)
    term = ((( y(3,i) * 0.25D+00 * s + y(2,i) / 3.0D+00 ) * s &
      + y(1,i) * 0.5D+00 ) * s + ftab(i) ) * s
    e(i+1) = e(i) + term
  end do
!
!  Determine where the endpoints A and B lie in the mesh of XTAB's.
!
  r = a
  u = 1.0D+00
 
  do j = 1, 2
!
!  The endpoint is less than or equal to XTAB(1).
!
    if ( r <= xtab(1) ) then
      result = result - u * ( ( r - xtab(1) ) * y(1,1) * 0.5D+00 &
        + ftab(1) ) * ( r - xtab(1) )
!
!  The endpoint is greater than or equal to XTAB(NTAB).
!
    else if ( xtab(ntab) <= r ) then

      result = result -u * ( e(ntab) + ( r - xtab(ntab) ) &
        * ( ftab(ntab) + 0.5D+00 * ( ftab(ntab-1) &
        + ( xtab(ntab) - xtab(ntab-1) ) * y(1,ntab-1) ) &
        * ( r - xtab(ntab) )))
!
!  The endpoint is strictly between XTAB(1) and XTAB(NTAB).
!
    else

      do i = 1, ntab-1
 
        if ( r <= xtab(i+1) ) then
          r = r - xtab(i)
          result = result - u * ( e(i) + ( ( ( &
              y(3,i) * 0.25D+00  * r &
            + y(2,i) / 3.0D+00 ) * r &
            + y(1,i) * 0.5D+00 ) * r + ftab(i) ) * r )
          go to 120
        end if
 
      end do
 
    end if
 
  120   continue
 
    u = -1.0D+00
    r = b
 
  end do
 
  return
end subroutine

!------------------------------------------------------------------------------------

subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

!*****************************************************************************
!
!! LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!
!    Input, real ( kind = 8 ) XD(ND), the interpolation nodes.
!
!    Input, integer ( kind = 4 ) NI, the number of evaluation points.
!
!    Input, real ( kind = 8 ) XI(NI), the evaluation points.
!
!    Output, real ( kind = 8 ) LB(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      lb(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end subroutine

!--------------------------------------------------------------------------

subroutine lagrange_value_1d ( nd, xd, yd, ni, xi, yi )

!*****************************************************************************
!
!! LAGRANGE_VALUE_1D evaluates the Lagrange interpolant.
!
!  Discussion:
!
!    The Lagrange interpolant L(ND,XD,YD)(X) is the unique polynomial of
!    degree ND-1 which interpolates the points (XD(I),YD(I)) for I = 1
!    to ND.
!
!    The Lagrange interpolant can be constructed from the Lagrange basis
!    polynomials.  Given ND distinct abscissas, XD(1:ND), the I-th Lagrange 
!    basis polynomial LB(ND,XD,I)(X) is defined as the polynomial of degree 
!    ND - 1 which is 1 at  XD(I) and 0 at the ND - 1 other abscissas.
!
!    Given data values YD at each of the abscissas, the value of the
!    Lagrange interpolant may be written as
!
!      L(ND,XD,YD)(X) = sum ( 1 <= I <= ND ) LB(ND,XD,I)(X) * YD(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ND, the number of data points.
!    ND must be at least 1.
!
!    Input, real ( kind = 8 ) XD(ND), the data points.
!
!    Input, real ( kind = 8 ) YD(ND), the data values.
!
!    Input, integer ( kind = 4 ) NI, the number of interpolation points.
!
!    Input, real ( kind = 8 ) XI(NI), the interpolation points.
!
!    Output, real ( kind = 8 ) YI(NI), the interpolated values.
!
  implicit none

  integer ( kind = 4 ) nd
  integer ( kind = 4 ) ni

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) lb(ni,nd)
  real ( kind = 8 ) xd(nd)
  real ( kind = 8 ) yd(nd)
  real ( kind = 8 ) xi(ni)
  real ( kind = 8 ) yi(ni)

  call lagrange_basis_1d ( nd, xd, ni, xi, lb )

  yi = matmul ( lb, yd )

  return
end subroutine


end module dopcasimir
