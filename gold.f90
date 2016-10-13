module gold

real(8), public:: freqA
!parameters (Drude model)
!vector format: (wp, damp)
!Mostepanenko
real(8), dimension(2), parameter:: parMost=(/13.671e15, 0.5317e14/)
!approximation data of the imaginary part (Palik)
real(8), dimension(2), parameter:: parAproxIm=(/10.06e15, 1.36e14/)
!approximation data of the real part (Palik)
real(8), dimension(2), parameter:: parAproxRe=(/12.06e15, 1.219e14/)

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

subroutine qnc79 ( func, a, b, err, result, ierr, k )

!*****************************************************************************80
!
!! QNC79 approximates the integral of F(X) using Newton-Cotes quadrature.
!
!  Discussion:
!
!    QNC79 is a general purpose program for evaluation of one
!    dimensional integrals  of user defined functions.  QNC79 will
!    pick its own points for evaluation of the integrand and these
!    will vary from problem to problem.
!
!    Thus QNC79 is not designed to integrate over data sets.
!
!    Moderately smooth integrands will be integrated efficiently
!    and reliably.  For problems with strong singularities,
!    oscillations etc., the user may wish to use more sophisticated
!    routines such as those in QUADPACK.
!
!    One measure of the reliability of QNC79 is the output parameter
!    K, giving the number of integrand evaluations that were needed.
!
!  Modified:
!
!    30 October 2000
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
!  Parameters:
!
!    Input, real ( kind = 8 ), external FUNC, the name of the function to be
!    integrated.  The user must declare the name an external
!    parameter in the calling program, pass the name of the
!    function in FUNC, and write a function of the form
!      FUNCTION FUNC ( X )
!    which evaluates the function at the point X.
!
!    Input, real ( kind = 8 ) A, lower limit of integral.
!
!    Input, real ( kind = 8 ) B, upper limit of integral.
!
!    Input, real ( kind = 8 ) ERR, is a requested error tolerance.
!    Normally pick a value, 0 < ERR < 1.0D-03.
!
!    Output, real ( kind = 8 ) RESULT, computed value of the integral.
!    Hopefully RESULT is accurate to within ERR times the
!    integral of ABS ( FUNC ( X ) ).
!
!    Output, integer ( kind = 4 ) IERR, a status code
!     1 RESULT most likely meets requested error tolerance.
!    -1 A and B are too nearly equal to allow normal integration.
!     2 RESULT probably does not meet requested error tolerance.
!
!    Output, integer ( kind = 4 ) K, the number of function evaluations
!    actually used to do the integration.
!    A value of K .GT. 1000 indicates a difficult problem.
!    Other programs may be more efficient.
!    QNC79 will gracefully give up if K exceeds 2000.
!
  implicit none

  integer ( kind = 4 ), parameter :: kmx = 2000

  real ( kind = 8 ) a
  real ( kind = 8 ) aa(40)
  real ( kind = 8 ) ae
  real ( kind = 8 ) area
  real ( kind = 8 ) b
  real ( kind = 8 ) bank
  real ( kind = 8 ) blocal
  real ( kind = 8 ) c
  real ( kind = 8 ) ce
  real ( kind = 8 ) ee
  real ( kind = 8 ) ef
  real ( kind = 8 ) eps
  real ( kind = 8 ) err
  real ( kind = 8 ) f(13)
  real ( kind = 8 ) f1(40)
  real ( kind = 8 ) f2(40)
  real ( kind = 8 ) f3(40)
  real ( kind = 8 ) f4(40)
  real ( kind = 8 ) f5(40)
  real ( kind = 8 ) f6(40)
  real ( kind = 8 ) f7(40)
  real ( kind = 8 ), external :: func
  real ( kind = 8 ) hh(40)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save :: icall = 0
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ), save :: kml = 7
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lmn
  integer ( kind = 4 ) lmx
  integer ( kind = 4 ) lr(40)
  integer ( kind = 4 ) nbits
  integer ( kind = 4 ) nib
  integer ( kind = 4 ), save :: nlmn = 2
  integer ( kind = 4 ) nlmx
  real ( kind = 8 ) q13
  real ( kind = 8 ) q7
  real ( kind = 8 ) q7l
  real ( kind = 8 ) q7r(40)
  real ( kind = 8 ) result
  real ( kind = 8 ) test
  real ( kind = 8 ) tol
  real ( kind = 8 ) vl(40)
  real ( kind = 8 ) vr
  real ( kind = 8 ) w1
  real ( kind = 8 ) w2
  real ( kind = 8 ) w3
  real ( kind = 8 ) w4

  if ( a == b ) then
    result = 0.0D+00
    return
  end if
 
  if ( icall /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  QNC79 was called recursively!'
    stop 1
  end if
 
  icall = 1
  w1 = 41.0D+00 / 140.0D+00
  w2  = 216.0D+00 / 140.0D+00
  w3 = 27.0D+00 / 140.0D+00
  w4  = 272.0D+00 / 140.0D+00
!
!  DIGITS ( X ) = number of base 2 digits in representation of X.
!
  nbits = int ( log10 ( 2.0D+00 ) &
    * real ( digits ( result ), kind = 8 ) / 0.30102000 )

  nlmx = min ( 40, ( nbits * 4 ) / 5 )
  result = 0.0D+00
  ierr = 1
  ce = 0.0D+00
  lmx = nlmx
  lmn = nlmn

  if ( b == 0.0D+00 ) then
    go to 3
  end if

  if ( sign ( 1.0D+00, b ) * a <= 0.0D+00 ) then
    go to 3
  end if

  c = abs ( 1.0D+00 - a / b )

  if ( 0.1D+00 < c ) then
    go to 3
  end if
 
  if ( c <= 0.0D+00 ) then
    ierr  = -1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Fatal error!'
    write ( *, '(a)' ) '  A and B are too close.'
    stop 1
  end if
 
  nib = int ( 0.5D+00 - log ( c ) / log(2.0D+00) )
  lmx = min ( nlmx, nbits-nib-4 )

  if ( lmx < 2 ) then
    go to 32
  end if

  lmn = min ( lmn, lmx )
 
3 continue
 
  tol = max ( abs ( err ), 2.0D+00**(5-nbits) )
  if ( err == 0.0D+00 ) then
    tol = sqrt ( epsilon ( tol ) )
  end if

  eps = tol
  hh(1) = ( b - a ) / 12.0D+00
  aa(1) = a
  lr(1) = 1
 
  do i = 1, 11, 2
    f(i) = func ( a + real ( i - 1, kind = 8 ) * hh(1) )
  end do
 
  blocal = b
  f(13) = func ( blocal )
  k = 7
  l = 1
  area = 0.0D+00
  q7 = 0.0D+00
  ef = 256.0D+00 / 255.0D+00
  bank = 0.0D+00
!
!  Compute refined estimates, estimate the error, etc.
!
5 continue
 
  do i = 2, 12, 2
    f(i) = func ( aa(l) + real ( i - 1, kind = 8 ) * hh(l) )
  end do
 
  k = k + 6
!
!  Compute left and right half estimates.
!
  q7l = hh(l) * ( ( w1 * ( f(1) + f(7) )   &
                  + w2 * ( f(2) + f(6) ) ) &
                + ( w3 * ( f(3) + f(5) ) + w4 * f(4) ) )

  q7r(l) = hh(l) * ( ( w1 * ( f(7) + f(13) ) + w2 * ( f(8) + f(12) ) ) &
                + ( w3 * ( f(9) + f(11) ) + w4 * f(10) ) )
!
!  Update estimate of integral of absolute value.
!
  area = area + ( abs ( q7l ) + abs ( q7r(l) ) - abs ( q7 ) )
!
!  Do not bother to test convergence before minimum refinement level.
!
  if ( l < lmn ) then
    go to 11
  end if
!
!  Estimate the error in new value for whole interval, Q13.
!
  q13 = q7l + q7r(l)
  ee = abs ( q7 - q13 ) * ef
!
!  Compute nominal allowed error.
!
  ae = eps * area
!
!  Borrow from bank account, but not too much.
!
  test = min ( ae + 0.8D+00 * bank, 10.0D+00 * ae )
!
!  Don't ask for excessive accuracy.
!
  test = max ( test, tol * abs ( q13 ), 0.00003D+00 * tol * area )
!
!  Now, did this interval pass or not?
!
  if ( ee <= test ) then
    go to 8
  end if

  go to 10
!
!  Have hit max refinement level - penalize the cumulative error.
!
7 continue
 
  ce = ce + ( q7 - q13 )
  go to 9
!
!  On good intervals accumulate the theoretical estimate.
!
8 continue

  ce = ce + ( q7 - q13 ) / 255.0D+00
!
!  Update the bank account.  Don't go into debt.
!
9 continue

  bank = bank + ( ae - ee )
  if ( bank < 0.0D+00 ) then
    bank = 0.0D+00
  end if
!
!  Did we just finish a left half or a right half?
!
  if ( lr(l) <= 0 ) then
    go to 15
  end if

  go to 20
!
!  Consider the left half of the next deeper level.
!
10 continue
 
  if ( kmx < k ) then
    lmx = min ( kml, lmx )
  end if

  if ( lmx <= l ) then
    go to 7
  end if

11 continue

  l = l + 1
  eps = eps * 0.5D+00

  if ( l <= 17 ) then
    ef = ef / sqrt ( 2.0D+00 )
  end if

  hh(l) = hh(l-1) * 0.5D+00
  lr(l) = -1
  aa(l) = aa(l-1)
  q7 = q7l
  f1(l) = f(7)
  f2(l) = f(8)
  f3(l) = f(9)
  f4(l) = f(10)
  f5(l) = f(11)
  f6(l) = f(12)
  f7(l) = f(13)
  f(13) = f(7)
  f(11) = f(6)
  f(9)  = f(5)
  f(7)  = f(4)
  f(5)  = f(3)
  f(3)  = f(2)
  go to 5
!
!  Proceed to right half at this level
!
15 continue

  vl(l) = q13

16 continue

  q7 = q7r(l-1)
  lr(l) = 1
  aa(l) = aa(l) + 12.0D+00 * hh(l)
  f(1)  = f1(l)
  f(3)  = f2(l)
  f(5)  = f3(l)
  f(7)  = f4(l)
  f(9)  = f5(l)
  f(11) = f6(l)
  f(13) = f7(l)
  go to 5
!
!  Left and right halves are done, so go back up a level
!
20 continue

  vr = q13

22 continue

  if ( l <= 1 ) then
    go to 30
  end if

  if ( l <= 17 ) then
    ef = ef * sqrt ( 2.0D+00 )
  end if

  eps = eps * 2.0D+00
  l = l-1
 
  if ( lr(l) <= 0 ) then
    vl(l) = vl(l+1) + vr
    go to 16
  else
    vr = vl(l+1) + vr
    go to 22
  end if
 
   30 continue
 
  if ( 2.0D+00 * tol * area < abs ( ce ) ) then
    ierr = 2
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'QNC79 - Warning!'
    write ( *, '(a)' ) '  RESULT is probably insufficiently accurate.'
  end if
 
32    continue
 
  result = vr
 
  icall = 0
 
  return
end subroutine

end module



