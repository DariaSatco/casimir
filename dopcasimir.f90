module dopcasimir

!contains subroutines for integration and some functions for Lifshitz formula

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
    real(8), parameter :: c=299792458._8 ! c - velocity of light

    integer, public :: model(2)

contains

!------------------------------------------------------------------------------------------------

	function eqf(a,x,y)
	! function needed to integrate in formula for Casimir free energy 
	! x=dist*r, y=w*dist/c change the variables to dimensioness
        ! the constants in front of the integral: rdr --> 1/dist**2 xdx;  dw --> c/dist

	real(8):: a,x,q,y,rTMn(2),rTEn(2),eqf

		rTMn=rTM(x/a,y)
!calculate the TM Frenel coefficient for the concrete dist as a function of x,y
		rTEn=rTE(x/a,y)
!calculate the TE Frenel coefficient for the concrete dist as a function of x,y
			
		q=sqrt(x**2+(a*y/c)**2)
		eqf=x*(log(1._8-rTMn(1)*rTMn(2)*exp(-2.*q))+log(1._8-rTEn(1)*rTEn(2)*exp(-2.*q)))

					
	end function eqf

!-----------------------------------------------------------------------------------------------

    function epsmodel_func(x)
    use gold
    use silicon
    !epsmodel_func choose from the list of models desired one

    real(8):: epsmodel_func(2)
    real(8):: x

    do i=1,2

    select case(model(i))
        case (1)
        epsmodel_func(i)=epsMar(x)
        case (2)
        epsmodel_func(i)=Drude(x, parMost)
        case (3)
        epsmodel_func(i)=epsSil_art(x)
    end select

    end do

    end function


!-----------------------------------------------------------------------------------------------

    function zerotempint(arg)
	! function - argument of 2-dimendional integral by frequencies and k-frequencies at T=0 
        ! limits [1, Inf) --> (0,1]
        ! change variables xdx -->  1/(arg1**3)darg1;  dy --> 1/(arg2**2)darg2

	real(8):: arg(2)
	real(8):: q, k(2), a
	real(8):: zerotempint
	real(8):: rTMn(2), rTEn(2), epsmodel(2)

    a=dist
	q=sqrt(arg(1)**2+arg(2)**2)

    epsmodel=epsmodel_func(arg(2)*c/a)

	!need to write rTE(arg1,arg2), rTM(arg1,arg2)
	do i=1,2
	   k(i)=sqrt((arg(1)/a)**2+epsmodel(i)*(arg(2)/a)**2)
       rTEn(i)=(q-k(i))/(q+k(i))
       rTMn(i)=(epsmodel(i)*q-k(i))/(epsmodel(i)*q+k(i))
    end do

    if (arg(1)*arg(2) .eq. 0.0_8) then

        zerotempint=0.0_8
        else
	    zerotempint=arg(1)*(log(1._8-rTMn(1)*rTMn(2)*exp(-2.*q))+log(1._8-rTEn(1)*rTEn(2)*exp(-2.*q)))
    end if

	end function

!-----------------------------------------------------------------------------------------------
    function zerotempint1(x)

    real(8):: x(2), zerotempint1

    if (x(1)*x(2) .eq. 0._8) then

        zerotempint1=0.0_8
        else
        zerotempint1=(1._8/x(1))**2*(1._8/x(2))**2*zerotempint(x**(-1))

    end if

    end function
!-----------------------------------------------------------------------------------------------

	
	function to_int(arg)
	!integrable function for nonzero temperatures - function of one argument - x <--> k-frequency
        !interval [0,1]

	real(8):: arg, to_int

	to_int = eqf(dist, arg, w)

	end function

!------------------------------------------------------------------------------------------------

	function to_int1(arg)
	!integrable function for nonzero temperatures
        !make change of interval for integration [1, Inf) --> (0,1]

	real(8):: arg, to_int1
	
	if (arg.eq.(0._8)) then
		to_int1=0._8
	else

		to_int1 = (1._8/arg**2)*eqf(dist, 1._8/arg, w)

	endif

	end function

!----------------------------------------------------------------------------------

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

!---------------------------------------------------------------------------------

!*****************************
!INTEGRATION SUBROUTINES
!*****************************

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

	subroutine trapzd2d(func,a,b,res,stepn)
	! approximates the integral with the trapezoidal method

	real(8):: func
	real(8):: res
	real(8), allocatable:: q(:,:)
	real(8):: a(2), b(2), x(2)
	real(8):: s,it(2)
	integer:: i,j,n,stepn
	
	s=0.
	res=0.25*abs(b(1)-a(1))*abs(b(2)-a(2))*(func(a)+func(b))
	i=0

		do while (abs(res-s).gt.(1e-2*abs(res)))
			s=res
			res=0._8
			i=i+1
			n=2**i
			it=abs(b-a)/n

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
				
				    x(1)=a(1)+j*it(1)
				    x(2)=a(2)+k*it(2)
					res=res+q(j,k)*(func(x))
				enddo

				enddo
			deallocate(q)

		res=res*it(1)*it(2)

		enddo
	
	stepn=n

	return

	end subroutine
	
!-----------------------------------------------------------------------------------------


end module dopcasimir
