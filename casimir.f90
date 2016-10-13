program casimir

use dopcasimir
use spline
use silicon
use gold
use quadpack

implicit none

!counters
integer:: i,j,k,l,s
!---------------------------------------
!casimir energy variables

integer:: N,stepn
real(8):: res, rest0, ssum0, ssum01
real(8), allocatable:: ssum(:), fenergy(:), ssum1(:), fenergy1(:)
real(8):: T
real(8):: x3(30), y3(30),z3(30,30)

!---------------------------------------
!silicon permittivity variables

integer, parameter:: nsil=483 !number of rows in the file
real(8):: matrixSil(nsil,3), epsteorSil(100)
real(8):: x(100) !frequency vector
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funval(nsil),Spcoef(3,nsil), exint(nsil), work(nsil)
!variables for Kr-K integral in qags subroutine
real(4):: abserr,abserr2
integer(4)::neval,ier,neval2,ier2
!integral-results in Kr-Kr
real(8):: integral1,integral2,integral3
!Silicon permittivity
real(8), allocatable:: epsSil(:)
!parameters for article formula
real(8):: eps0, epsinf, w0

!--------------------------------------------------------------------
!gold permettivity variables
integer,parameter:: n1=310 !rows number in file
!Palik data table
real(8):: matrixAu(n1,3), epsAurDr(100), epsAurMar(100), x1(100)
!permettivity vector
real(8), allocatable:: epsAur(:)
!integration results
real(8):: integralA1,integralA2,integralA3
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funvalAu(n1),SpcoefA(3,n1), exintA(n1), workA(n1)
!variables for Kr-K integral in qags subroutine
real(4):: abserrA,abserrA2
integer(4)::nevalA,ierA,nevalA2,ierA2



!---------------------------------------------------------------------
!code for calculation of silicon permettivity

open(unit=10, file='resultSi.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,nsil
read(10,*), (matrixSil(nsil+1-i,j),j=1,3)
end do


!Kramers-Kronig for 'freq'-vector

allocate(epsSil(100))
freq=0._8

	do i=1,100
	
	!integrates the oscillator-function
	call qags(KramersKr, 0.0, sngl(matrixSil(1,1)), 1.0e-5, 1.0e-3, sngl(integral1), abserr, neval,ier)
	
	!builds vector to integrate from the Palik data
		do k=1,nsil
		funval(k)=matrixSil(k,3)*matrixSil(k,1)/(matrixSil(k,1)**2+freq**2)
		end do
	!integrates the vector of data
 	call cspint(nsil, matrixSil(1:483,1), funval, matrixSil(1,1), matrixSil(nsil,1), Spcoef, exint, work, integral2)

	!integrates the oscillator-function
	call qags(KramersKr, sngl(matrixSil(nsil,1)), sngl(matrixSil(nsil,1)*1e4_8), 1.0e-5, &
 1.0e-3,sngl(integral3), abserr2, neval2, ier2)

	
	!Kr-Kr formula for permittivity 
	epsSil(i)=(integral1+integral2+integral3)*2/pi+1

	x(i)=freq
	freq=(freq+20._8)*10**(0.2_8)

	end do

!print *, epsSil

!theory - Drude-Loretz model
!parameters from the article
eps0=11.87_8
epsinf=1.035_8
w0=6.6e+15_8

!article formula

do i=1,100
    epsteorSil(i)=epsinf+(eps0-epsinf)*w0**2/(w0**2+x(i)**2);
end do

!call spline_b_val(100, x, epsSil, x(30)+1.0e4_8, splineint)
!print *, 'splineint=', splineint, epsSil(30)

open(unit=17, file='silicondata.txt', status='replace')

!write data in file
do i=1,100
write(17,*) x(i), epsSil(i), epsteorSil(i)
end do



!---------------------------------------------------------------------
!code for calculation of gold permettivity

open(unit=15, file='resultAu.txt', status='old')
!read matrix from file
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,n1
read(15,*), (matrixAu(n1+1-i,j),j=1,3)
end do

200 format(' ', e15.3, f15.7, f15.7)

!do i=1,n1
!write(*,200), (matrixAu(i,j),j=1,3)
!end do

!Kramers-Kronig for 'freq'-vector

allocate(epsAur(100))
freqA=5._8

	do i=1,100
	
	!integrates the Drude-model approximation
	call qnc79(drude_to_int, 0.0_8, matrixAu(1,1), 1.0e-3_8, integralA1, ierA, nevalA)
	!PRINT*, 'INTEGRAL1=', integralA1
	!print*, 'ier=', ierA

	
	!builds vector to integrate from the Palik data
		do k=1,n1
		funvalAu(k)=matrixAu(k,3)*matrixAu(k,1)/(matrixAu(k,1)**2+freqA**2)
		end do
	!integrates the vector of data
 	call cspint(n1, matrixAu(1:n1,1), funvalAu, matrixAu(1,1), matrixAu(n1,1), SpcoefA, exintA, workA, integralA2)
	!print*, 'integral2=', integralA2

	!integrates the oscillator-function
	call monte_carlo(oscil_to_int, matrixAu(n1,1), matrixAu(n1,1)*1.0e4_8, 1000000, integralA3)
	!print*, 'integral3=', integralA3
	!print*, 'ier3=', ierA2

	
	!Kr-Kr formula for permittivity 
	epsAur(i)=(integralA1+integralA2+integralA3)*2/pi+1

	x1(i)=freqA
	freqA=(freqA+20.)*10**(0.2_8)

	end do

open(unit=18, file='golddata.txt', status='replace')

do i=1,100
epsAurDr(i)=Drude(x1(i),parAproxIm)
epsAurMar(i)=epsMar(x1(i))
end do

do i=1,100
write(18,*) x1(i), epsAur(i), epsAurDr(i), epsAurMar(i)
end do


!---------------------------------------------------------------------
!code for calulation of casimir free energy

T=300._8                            
!T - Temperature of the material

eps(3)=1._8
!fix the dielectric permettivity of the gap (vacuum)

print *, 'type number of points '    
read *, N
! N - number of points on the plot

!open(unit=11, file='casimirsilicon.txt', status='replace')
open(unit=11, file='casimirgold.txt', status='replace')
!casimir.txt - file with computational results

allocate(ssum(N))
allocate(fenergy(N))
allocate(ssum1(N))
allocate(fenergy1(N))

	do j=1,N     
	!do-cycle for distance
	dist=1e-7*j	
	ssum(j)=0._8

	call trapzd2d(zerotempint,0._8,1._8,1e-8_8,1._8,rest0,stepn) 
x3(1)=0._8
y3(1)=0._8
do i=1,29
x3(i+1)=(1._8)/30*i
y3(i+1)=(1._8)/30*i
end do
do i=1,30
do s=1,30
z3(i,s)=zerotempint(x3(i),y3(s))
end do
end do
	!call surf(x3,y3,z3)
	!integrate for T=0 to find numberof sum members
	print *, 'n=', stepn                                        
	!stepn - number of steps in trapezoidal integration method

	ssum(j)=0.
	k=int((stepn*c/dist)/(2*pi*kb*T/h))      
	!number of sum members - Matsubara frequencies
	print *, 'k=', k
		
		!calculate casimir energy
		!Kramers-Kronig model for silicon
		do i=1,k
		!do-cycle for frequencies

			w=2*i*pi*kb*T/h
			!w - Matsubara frequency
			!print *, 'w=', w
			!call spline_b_val(100, x, epsSil, w, eps(1))
			call spline_b_val(100, x1, epsAur, w, eps(1))
			!print*, 'eps=', eps(1)
			!call spline_b_val(100, x, epsSil, w, eps(2))
			call spline_b_val(100, x1, epsAur, w, eps(1))

			call trapzd(to_int,0._8,1._8,res) 
			!print*, 'to_int=', to_int(0.5_8)
			!print*, 'res=', res
			ssum(j)=ssum(j)+res
			call trapzd(to_int1,0._8,1._8,res)
			!print*,'res=', res
			ssum(j)=ssum(j)+res
			
	
		end do
	
	ssum0=0._8
	!calculate the zeroth summand
	!eps(1)=epsSil(1)
	!eps(2)=epsSil(1)
	eps(1)=epsAur(1)
	eps(2)=epsAur(1)
	call trapzd(to_int,0._8,1._8,res) 
	ssum0=ssum0+res
	call trapzd(to_int1,0._8,1._8,res)
	ssum0=ssum0+res
	

	fenergy(j)=kb*1.602e-19*T/(2*pi*dist**2)*(0.5*ssum0+ssum(j))
	!free casimir energy

		!calculate casimir free energy
		!article model fo silicon
		do i=1,k
		!do-cycle for frequencies

			w=2*i*pi*kb*T/h
			!w - Matsubara frequency
			!print *, 'w=', w
			!call spline_b_val(100, x, epsteorSil, w, eps(1))
			call spline_b_val(100, x1, epsAurMar, w, eps(1))			
			!print*, 'eps=', eps(1)
			!call spline_b_val(100, x, epsteorSil, w, eps(2))
			call spline_b_val(100, x1, epsAurMar, w, eps(2))

			call trapzd(to_int,0._8,1._8,res) 
			!print*, 'to_int=', to_int(0.5_8)
			!print*, 'res=', res
			ssum1(j)=ssum1(j)+res
			call trapzd(to_int1,0._8,1._8,res)
			!print*,'res=', res
			ssum1(j)=ssum1(j)+res
			
	
		end do
	
	ssum01=0._8
	!calculate the zeroth summand
	!eps(1)=epsteorSil(1)
	!eps(2)=epsteorSil(1)
	eps(1)=epsAurMar(1)
	eps(2)=epsAurMar(1)
	call trapzd(to_int,0._8,1._8,res) 
	ssum01=ssum01+res
	call trapzd(to_int1,0._8,1._8,res)
	ssum01=ssum01+res
	

	fenergy1(j)=kb*1.602e-19*T/(2*pi*dist**2)*(0.5*ssum01+ssum1(j))



	write(11,*) dist*1e7, abs((fenergy(j))*1e9), abs((fenergy1(j))*1e9)
	
	end do

deallocate(ssum)
deallocate(fenergy)
deallocate(ssum1)
deallocate(fenergy1)

close(11)
close(10)
close(17)
deallocate(epsSil)

deallocate(epsAur)
close(15)
close(18)

end program


	
