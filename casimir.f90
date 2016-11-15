program casimir

use dopcasimir
use spline
use silicon
use gold
!use quadpack

implicit none

!counters
integer:: i,j,k,l,s
!---------------------------------------
!casimir energy variables

integer:: N,stepn
real(8):: res, rest0, ssum0, ssum01, res1
real(8), allocatable:: ssum(:), fenergy(:), ssum1(:), fenergy1(:), ecasimir(:)
real(8):: T, wmax

!---------------------------------------
!silicon permittivity variables

integer, parameter:: nsil=483 !number of rows in the file
real(8):: matrixSil(nsil,3)
real(8):: x(100) !frequency vector
!variables needed fo Kr-Kr integral in cspint subroutine
real(8)::funval(nsil),Spcoef(3,nsil), exint(nsil), work(nsil)
!variables for Kr-K integral in qags subroutine
real(4):: abserr,abserr2
integer(4)::neval,ier,neval2,ier2, ierr
!integral-results in Kr-Kr
real(8):: integral1,integral2,integral3
!Silicon permittivity
real(8), allocatable:: epsSil(:),  epsteorSil(:)
real(8):: a(2), b(2)

!--------------------------------------------------------------------
!gold permettivity variables
!integer,parameter:: n1=310 !rows number in file
!!Palik data table
!real(8):: matrixAu(n1,3), epsAurDr(100), epsAurMar(100), x1(100)
!!permettivity vector
!real(8), allocatable:: epsAur(:)
!!integration results
!real(8):: integralA1,integralA2,integralA3
!!variables needed fo Kr-Kr integral in cspint subroutine
!real(8)::funvalAu(n1),SpcoefA(3,n1), exintA(n1), workA(n1)
!!variables for Kr-K integral in qags subroutine
!real(4):: abserrA,abserrA2
!integer(4)::nevalA,ierA,nevalA2,ierA2


!-----------------------------------------------------------------------
open(unit=11, file='casimirsilicon.txt', status='replace')
!casimirsilicon.txt, casimirgold.txt - files with computational results

!-----------------------------------------------------------------------
!code for calulation of casimir free energy

open(unit=10, file='resultSi.txt', status='old')
!read matrix from file (silicon)
!1st column - frequencies
!2nd column - real part of permittivity
!3rd column - imaginary part of permittivity

do i=1,nsil
read(10,*), (matrixSil(nsil+1-i,j),j=1,3)
end do

T=300._8
!T - Temperature of the material

eps(3)=1._8
!fix the dielectric permettivity of the gap (vacuum)

print *, 'type number of points '
read *, N
! N - number of points on the plot

allocate(ssum(N))
allocate(fenergy(N))
allocate(ssum1(N))
allocate(fenergy1(N))
allocate(ecasimir(N))

model(1)=3
model(2)=3

    do j=1,N
    !do-cycle for distance
    dist=1.0e-7*j
    ssum(j)=0._8
    ssum1(j)=0._8

    ssum0=0._8
    ssum01=0._8

    !first calculate zero-temperature casimir energy
    a=(/0._8, 0._8/)
    b=(/1._8, 1._8/)
    call monte_carlo_nd (zerotempint, 2, a,b, 1000500, res)
    call trapzd2d(zerotempint,a,b,res1,stepn)
    rest0=res
    call monte_carlo_nd (zerotempint1, 2,  a,b, 1000500, res)
    call trapzd2d(zerotempint1,a,b,res1,stepn)
    rest0=(rest0+res)*c/dist**3

    k=1
    wmax=stepn*c/dist
    w=2*pi*kb*T/h
    do while (w.le.wmax)
    k=k+1
    w=2*k*pi*kb*T/h
    end do

print*, 'k=', k

allocate(epsSil(k))
allocate(epsteorSil(k))

        !calculate casimir energy
        do i=1,k
        !do-cycle for frequencies

            w=2*i*pi*kb*T/h
            !w - Matsubara frequency
            freq=w

!-----------------------------------------------------------------------
!silicon dielectric permettivity function in point = freq = w

    !integrates the oscillator-function
    call qags(KramersKr, 0.0, sngl(matrixSil(1,1)), 1.0e-5, 1.0e-3, &
    sngl(integral1), abserr, neval,  ierr)

    !builds vector to integrate from the Palik data
        do k=1,nsil
        funval(k)=matrixSil(k,3)*matrixSil(k,1)/(matrixSil(k,1)**2+freq**2)
        end do
    !integrates the vector of data
    call cspint(nsil, matrixSil(1:483,1), funval, matrixSil(1,1), matrixSil(nsil,1),&
    Spcoef, exint, work, integral2)

    !integrates the oscillator-function
    call qags(KramersKr, sngl(matrixSil(nsil,1)), sngl(matrixSil(nsil,1)*1e4_8), 1.0e-5, &
 1.0e-3,sngl(integral3), abserr2, neval2, ier2)


    !Kr-Kr formula for permittivity
    epsSil(i)=(integral1+integral2+integral3)*2/pi+1.0_8

    !theory (article)
    epsteorSil(i)=epsSil_art(w)

!-----------------------------------------------------------------------
    !calculate casimir (my calculation)
    eps(1)=epsSil(i)
    eps(2)=epsSil(i)
    call qnc79(to_int,0._8,1._8, 1.0e-3_8, res, ierr, s)

            ssum(j)=ssum(j)+res

    call qnc79(to_int1,0._8,1._8, 1.0e-3_8, res, ierr, s)

            ssum(j)=ssum(j)+res


!-----------------------------------------------------------------------

    !calculate casimir (article)
    eps(1)=epsteorSil(i)
    eps(2)=epsteorSil(i)
    call trapzd(to_int,0._8,1._8,res)

            ssum1(j)=ssum1(j)+res

    call trapzd(to_int1,0._8,1._8,res)

            ssum1(j)=ssum1(j)+res

     end do

!*********
!ssum0
    call qnc79(zero_sum_int_calc,0._8,1._8, 1.0e-3_8, res, ierr, s)

            ssum0=ssum0+res

    call qnc79(zero_sum_int1_calc,0._8,1._8, 1.0e-3_8, res, ierr, s)

            ssum0=ssum0+res

 !ssum01
   call trapzd(zero_sum_int_art,0._8,1._8,res)

            ssum01=ssum01+res

   call trapzd(zero_sum_int1_art,0._8,1._8,res)

            ssum01=ssum01+res

!****************

fenergy(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum0+ssum(j))

fenergy1(j)=1.602e-19*kb*T/(2*pi*dist**2)*(0.5*ssum01+ssum1(j))

rest0=1.602e-19*h*rest0/(2*pi)**2

ecasimir(j)=1.602e-19*(pi)**2*h*c/(720*dist**3)

write(11,*) dist*1.0e7, abs((fenergy(j))*1.0e9), abs((fenergy1(j))*1.0e9), &
abs(rest0)*1.0e9, ecasimir(j)*1.0e9
!*********

deallocate(epsSil)
deallocate(epsteorSil)

    end do


!**********************************************************************************************



deallocate(ssum)
deallocate(fenergy)
deallocate(ssum1)
deallocate(fenergy1)

close(11)
close(10)

end program


	
