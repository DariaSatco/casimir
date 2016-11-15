module silicon

real(8), public:: freq
real(8), dimension(3), parameter::g= (/9.0e31_8, 2.8e32_8, 8.2e31_8/) !(/5.92e16_8, 1.84e17_8, 5.4e16_8/) eV/s
real(8), dimension(3), parameter::gammam=(/1.6e15_8, 1.05e15_8, 5.4e14_8/)  !(/1.05_8, 0.69_8, 0.35_8/)    eV/s
real(8), dimension(3), parameter::wm= (/8.32e15_8, 6.37e15_8, 5.33e15_8/)   !(/5.48_8, 4.19_8, 3.51_8/)   eV/s


contains

!--------------------------------------------------------------------------

function oscillators(x)

real(4):: x, oscillators

oscillators=0.0
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

function epsSil_art(x)

real(8):: epsSil_art
real(8):: x
real(8):: eps0, epsinf, w0

!silicon parameters
!theory - Drude-Loretz model
!parameters from the article
eps0=11.87_8
epsinf=1.035_8
w0=6.6e15_8

epsSil_art=epsinf+(eps0-epsinf)*w0**2/(w0**2+x**2)

end function

!---------------------------------------------------------------------

!Want to calculate zero summand

function zero_sum1(a,x)

real(8):: zero_sum1, x
real(8):: rTm, rTe, a
real(8):: eps0, epsinf, w0

!silicon parameters
!theory - Drude-Loretz model
!parameters from the article
eps0=11.87_8

rTe=0._8
rTm=(eps0-1._8)/(eps0+1._8)

zero_sum1=x*(log(1._8-rTm**2*exp(-2.*x))+log(1._8-rTe**2*exp(-2.*x)))

end function

!----------------------------------------------------------------------

function zero_sum2(a,x)

real(8):: zero_sum2, x
real(8):: rTm, rTe, a
real(8):: eps0


rTe=0._8
eps0=1._8

do i=1,3
eps0=eps0+g(i)/wm(i)**2
end do

rTm=(eps0-1)/(eps0+1)

zero_sum2=x*(log(1._8-rTm**2*exp(-2.*x))+log(1._8-rTe**2*exp(-2.*x)))

end function
!---------------------------------------------------------------------

function zero_sum_int_art(x)
use dopcasimir

real(8):: zero_sum_int_art, x

if (x.eq.(0._8)) then
    zero_sum_int=0._8
    else
zero_sum_int_art=zero_sum1(dist,x)
endif

end function
!---------------------------------------------------------------------

function zero_sum_int1_art(x)
use dopcasimir

real(8):: zero_sum_int1_art, x

if (x.eq.(0._8)) then

    zero_sum_int1=0._8
    else

    zero_sum_int1_art=(1._8/x)**2*zero_sum1(dist,1._8/x)
endif

end function
!---------------------------------------------------------------------


function zero_sum_int_calc(x)
use dopcasimir

real(8):: zero_sum_int_calc, x

if (x.eq.(0._8)) then
    zero_sum_int=0._8
    else
zero_sum_int_calc=zero_sum2(dist,x)
endif

end function
!---------------------------------------------------------------------

function zero_sum_int1_calc(x)
use dopcasimir

real(8):: zero_sum_int1_calc, x

if (x.eq.(0._8)) then

    zero_sum_int1=0._8
    else

    zero_sum_int1_calc=(1._8/x)**2*zero_sum2(dist,1._8/x)
endif

end function
!---------------------------------------------------------------------

end module
