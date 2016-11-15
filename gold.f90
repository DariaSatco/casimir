module gold

real(8), public:: freqA
!parameters (Drude model)
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter, public:: parMost=(/13.671e15, 0.5317e14/)
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter, public:: parAproxIm=(/10.06e15, 1.36e14/)
!approximation data of the real part (Palik)
real(8), dimension(2), parameter, public:: parAproxRe=(/12.06e15, 1.219e14/)

!oscillator model, Mostepanenko parameters
real(8), dimension(6), parameter:: gAu=(/7.091, 41.46, 2.7, 154.7, 44.55, 309.6/)*(1.519e15)**2
real(8), dimension(6), parameter:: gammaAu=(/0.75, 1.85, 1.0, 7.0, 6.0, 9.0/)*1.519e15
real(8), dimension(6), parameter:: wAu=(/3.05, 4.15, 5.4, 8.5, 13.5, 21.5/)*1.519e15


contains

function Drude(x, param)
!Drude-model permettivity function eps(iw)
!param - vector of parameters (wp, damp)
real(8):: x, Drude, param(2)

	Drude=1+param(1)**2/(x**2+param(2)*x);

end function

!----------------------------------------------------

function epsMar(x)
!Marachevsky formula fo eps(iw)
real(8):: x,epsMar,num,denum
real(8), dimension(2,4):: paramMar

!Marachevsky model parameters (wl1,wl2//gl1,gl2//gt1,gt2//wt2,0)
paramMar(1,1)=exp(-0.96)*1e16
paramMar(2,1)=exp(0.2866)*1e16
paramMar(1,2)=exp(-2.536)*1e16
paramMar(2,2)=exp(1.255)*1e16
paramMar(1,3)=exp(-4.7922)*1e16
paramMar(2,3)=exp(-0.957)*1e16
paramMar(1,4)=exp(-0.8359)*1e16
paramMar(2,4)=0.0_8

num=(paramMar(1,1)**2+x**2+paramMar(1,2)*x)*(paramMar(2,1)**2+x**2+paramMar(2,2)*x)
denum=(x**2+paramMar(1,3)*x)*(paramMar(1,4)**2+x**2+paramMar(2,3)*x)

epsMar=num/denum

end function

!---------------------------------------------------

function oscil_to_int(x)

real(8):: oscil,x,oscil_to_int

oscil=0.0_8

do i=1,6
oscil=oscil+gAu(i)*gammaAu(i)*x/((wAu(i)**2-x**2)**2+(gammaAu(i)*x)**2)
end do

oscil_to_int=oscil+parAproxRe(1)**2*parAproxRe(2)/((x**2+parAproxRe(2)**2)*(x**2+freqA**2))

end function

!--------------------------------------------------

function drude_to_int(x)

real(8):: x, drude_to_int

drude_to_int=parAproxRe(1)**2*parAproxRe(2)/((x**2+parAproxRe(2)**2)*(x**2+freqA**2))

end function

!-------------------------------------------------


end module



