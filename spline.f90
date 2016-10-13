module spline

contains 

subroutine basis_function_b_val ( tdata, tval, yval )

!*****************************************************************************80
!
!! BASIS_FUNCTION_B_VAL evaluates the B spline basis function.
!
!  Discussion:
!
!    The B spline basis function is a piecewise cubic which
!    has the properties that:
!
!    * it equals 2/3 at TDATA(3), 1/6 at TDATA(2) and TDATA(4);
!    * it is 0 for TVAL <= TDATA(1) or TDATA(5) <= TVAL;
!    * it is strictly increasing from TDATA(1) to TDATA(3),
!      and strictly decreasing from TDATA(3) to TDATA(5);
!    * the function and its first two derivatives are continuous
!      at each node TDATA(I).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Alan Davies, Philip Samuels,
!    An Introduction to Computational Geometry for Curves and Surfaces,
!    Clarendon Press, 1996,
!    ISBN: 0-19-851478-6,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) TDATA(5), the nodes associated with the
!    basis function.  The entries of TDATA are assumed to be distinct
!    and increasing.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the B spline basis
!    function is to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ), parameter :: ndata = 5

  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) tdata(ndata)
  real ( kind = 8 ) tval
  real ( kind = 8 ) u
  real ( kind = 8 ) yval

  if ( tval <= tdata(1) .or. tdata(ndata) <= tval ) then
    yval = 0.0D+00
    return
  end if
!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] containing TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  U is the normalized coordinate of TVAL in this interval.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
!
!  Now evaluate the function.
!
  if ( tval < tdata(2) ) then
    yval = u**3 / 6.0D+00
  else if ( tval < tdata(3) ) then
    yval = ( ( (    - 3.0D+00   &
                * u + 3.0D+00 ) &
                * u + 3.0D+00 ) &
                * u + 1.0D+00 ) / 6.0D+00
  else if ( tval < tdata(4) ) then
    yval = ( ( (    + 3.0D+00   &
                * u - 6.0D+00 ) &
                * u + 0.0D+00 ) &
                * u + 4.0D+00 ) / 6.0D+00
  else if ( tval < tdata(5) ) then
    yval = ( 1.0D+00 - u )**3 / 6.0D+00
  end if

  return

end subroutine

!-------------------------------------------------------------------

subroutine spline_b_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_B_VAL evaluates a cubic B spline approximant.
!
!  Discussion:
!
!    The cubic B spline will approximate the data, but is not
!    designed to interpolate it.
!
!    In effect, two "phantom" data values are appended to the data,
!    so that the spline will interpolate the first and last data values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data values.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, a point at which the spline is
!    to be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the function at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real ( kind = 8 ) bval
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) tdata(ndata)
  real ( kind = 8 ) tval
  real ( kind = 8 ) u
  real ( kind = 8 ) ydata(ndata)
  real ( kind = 8 ) yval
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the 5 nonzero B spline basis functions in the interval,
!  weighted by their corresponding data values.
!
  u = ( tval - tdata(left) ) / ( tdata(right) - tdata(left) )
  yval = 0.0D+00
!
!  B function associated with node LEFT - 1, (or "phantom node"),
!  evaluated in its 4th interval.
!
  bval = ( ( (     - 1.0D+00   &
               * u + 3.0D+00 ) &
               * u - 3.0D+00 ) &
               * u + 1.0D+00 ) / 6.0D+00

  if ( 0 < left-1 ) then
    yval = yval + ydata(left-1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(1) - ydata(2) ) * bval
  end if
!
!  B function associated with node LEFT,
!  evaluated in its third interval.
!
  bval = ( ( (       3.0D+00   &
               * u - 6.0D+00 ) &
               * u + 0.0D+00 ) &
               * u + 4.0D+00 ) / 6.0D+00

  yval = yval + ydata(left) * bval
!
!  B function associated with node RIGHT,
!  evaluated in its second interval.
!
  bval = ( ( (     - 3.0D+00   &
               * u + 3.0D+00 ) &
               * u + 3.0D+00 ) &
               * u + 1.0D+00 ) / 6.0D+00

  yval = yval + ydata(right) * bval
!
!  B function associated with node RIGHT+1, (or "phantom node"),
!  evaluated in its first interval.
!
  bval = u**3 / 6.0D+00

  if ( right+1 <= ndata ) then
    yval = yval + ydata(right+1) * bval
  else
    yval = yval + ( 2.0D+00 * ydata(ndata) - ydata(ndata-1) ) * bval
  end if

  return

end subroutine

!-------------------------------------------------------------------------------

subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return

end subroutine

!------------------------------------------------------------------------------

subroutine spline_overhauser_nonuni_val ( ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! SPLINE_OVERHAUSER_NONUNI_VAL evaluates the nonuniform Overhauser spline.
!
!  Discussion:
!
!    The nonuniformity refers to the fact that the abscissa values
!    need not be uniformly spaced.
!
!    Thanks to Doug Fortune for pointing out that the point distances
!    used to define ALPHA and BETA should be the Euclidean distances
!    between the points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NDATA, the number of data points.
!    3 <= NDATA is required.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissas of the data points.
!    The values of TDATA are assumed to be distinct and increasing.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values.
!
!    Input, real ( kind = 8 ) TVAL, the value where the spline is to
!    be evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ) ndata

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) d21
  real ( kind = 8 ) d32
  real ( kind = 8 ) d43
  integer ( kind = 4 ) left
  real ( kind = 8 ) mbasis(4,4)
  real ( kind = 8 ) mbasis_l(3,3)
  real ( kind = 8 ) mbasis_r(3,3)
  integer ( kind = 4 ) right
  real ( kind = 8 ) tdata(ndata)
  real ( kind = 8 ) tval
  real ( kind = 8 ) ydata(ndata)
  real ( kind = 8 ) yval
!
!  Check NDATA.
!
  if ( ndata < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
    write ( *, '(a)' ) '  NDATA < 3.'
    stop 1
  end if
!
!  Find the nearest interval [ TDATA(LEFT), TDATA(RIGHT) ] to TVAL.
!
  call r8vec_bracket ( ndata, tdata, tval, left, right )
!
!  Evaluate the spline in the given interval.
!
  if ( left == 1 ) then

    d21 = sqrt ( ( tdata(2) - tdata(1) )**2 &
               + ( ydata(2) - ydata(1) )**2 )

    d32 = sqrt ( ( tdata(3) - tdata(2) )**2 &
               + ( ydata(3) - ydata(2) )**2 )

    alpha = d21 / ( d32 + d21 )

    call basis_matrix_overhauser_nul ( alpha, mbasis_l )

    call basis_matrix_tmp ( left, 3, mbasis_l, ndata, tdata, ydata, tval, yval )

  else if ( left < ndata - 1 ) then

    d21 = sqrt ( ( tdata(left) - tdata(left-1) )**2 &
               + ( ydata(left) - ydata(left-1) )**2 )

    d32 = sqrt ( ( tdata(left+1) - tdata(left) )**2 &
               + ( ydata(left+1) - ydata(left) )**2 )

    d43 = sqrt ( ( tdata(left+2) - tdata(left+1) )**2 &
               + ( ydata(left+2) - ydata(left+1) )**2 )

    alpha = d21 / ( d32 + d21 )
    beta  = d32 / ( d43 + d32 )

    call basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

    call basis_matrix_tmp ( left, 4, mbasis, ndata, tdata, ydata, tval, yval )

  else if ( left == ndata - 1 ) then

    d32 = sqrt ( ( tdata(ndata-1) - tdata(ndata-2) )**2 &
               + ( ydata(ndata-1) - ydata(ndata-2) )**2 )

    d43 = sqrt ( ( tdata(ndata) - tdata(ndata-1) )**2 &
               + ( ydata(ndata) - ydata(ndata-1) )**2 )

    beta  = d32 / ( d43 + d32 )

    call basis_matrix_overhauser_nur ( beta, mbasis_r )

    call basis_matrix_tmp ( left, 3, mbasis_r, ndata, tdata, ydata, tval, yval )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_OVERHAUSER_NONUNI_VAL - Fatal error!'
    write ( *, '(a,i8)' ) '  Nonsensical value of LEFT = ', left
    write ( *, '(a,i8)' ) '  but 0 < LEFT < NDATA = ', ndata
    write ( *, '(a)' ) '  is required.'
    stop 1

  end if

  return
end subroutine

!-------------------------------------------------------------------------------

subroutine basis_matrix_overhauser_nul ( alpha, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NUL sets the nonuniform left Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, and
!    P3 are not uniformly spaced in T, and that P1 corresponds to T = 0,
!    and P2 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA.
!    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   1.0D+00 / alpha
  mbasis(1,2) = - 1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
  mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - alpha )

  mbasis(2,1) = - ( 1.0D+00 + alpha ) / alpha
  mbasis(2,2) =   1.0D+00 / ( alpha * ( 1.0D+00 - alpha ) )
  mbasis(2,3) = - alpha / ( 1.0D+00 - alpha )

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end subroutine

!-------------------------------------------------------------------------------

subroutine basis_matrix_overhauser_nur ( beta, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NUR: the nonuniform right Overhauser basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points PN-2, PN-1, and
!    PN are not uniformly spaced in T, and that PN-1 corresponds to T = 0,
!    and PN to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 August 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) BETA.
!    BETA = || P(N) - P(N-1) || 
!         / ( || P(N) - P(N-1) || + || P(N-1) - P(N-2) || )
!
!    Output, real ( kind = 8 ) MBASIS(3,3), the basis matrix.
!
  implicit none

  real ( kind = 8 ) beta
  real ( kind = 8 ) mbasis(3,3)

  mbasis(1,1) =   1.0D+00 / beta
  mbasis(1,2) = - 1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
  mbasis(1,3) =   1.0D+00 / ( 1.0D+00 - beta )

  mbasis(2,1) = - ( 1.0D+00 + beta ) / beta
  mbasis(2,2) =   1.0D+00 / ( beta * ( 1.0D+00 - beta ) )
  mbasis(2,3) = - beta / ( 1.0D+00 - beta )

  mbasis(3,1) =   1.0D+00
  mbasis(3,2) =   0.0D+00
  mbasis(3,3) =   0.0D+00

  return
end subroutine

!-----------------------------------------------------------------------------------

subroutine basis_matrix_tmp ( left, n, mbasis, ndata, tdata, ydata, tval, yval )

!*****************************************************************************80
!
!! BASIS_MATRIX_TMP computes Q = T * MBASIS * P
!
!  Discussion:
!
!    YDATA is a vector of data values, most frequently the values of some
!    function sampled at uniformly spaced points.  MBASIS is the basis
!    matrix for a particular kind of spline.  T is a vector of the
!    powers of the normalized difference between TVAL and the left
!    endpoint of the interval.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) LEFT, indicats that TVAL is in the interval
!    [ TDATA(LEFT), TDATA(LEFT+1) ], or that this is the "nearest"
!    interval to TVAL.
!    For TVAL < TDATA(1), use LEFT = 1.
!    For TDATA(NDATA) < TVAL, use LEFT = NDATA - 1.
!
!    Input, integer ( kind = 4 ) N, the order of the basis matrix.
!
!    Input, real ( kind = 8 ) MBASIS(N,N), the basis matrix.
!
!    Input, integer ( kind = 4 ) NDATA, the dimension of the vectors TDATA
!    and YDATA.
!
!    Input, real ( kind = 8 ) TDATA(NDATA), the abscissa values.  This routine
!    assumes that the TDATA values are uniformly spaced, with an
!    increment of 1.0.
!
!    Input, real ( kind = 8 ) YDATA(NDATA), the data values to be
!    interpolated or approximated.
!
!    Input, real ( kind = 8 ) TVAL, the value of T at which the spline is to be
!    evaluated.
!
!    Output, real ( kind = 8 ) YVAL, the value of the spline at TVAL.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxn = 4
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ndata

  real ( kind = 8 ) arg
  integer ( kind = 4 ) first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) left
  real ( kind = 8 ) mbasis(n,n)
  real ( kind = 8 ) tdata(ndata)
  real ( kind = 8 ) tval
  real ( kind = 8 ) tvec(maxn)
  real ( kind = 8 ) ydata(ndata)
  real ( kind = 8 ) yval

  if ( left == 1 ) then
    arg = 0.5D+00 * ( tval - tdata(left) )
    first = left
  else if ( left < ndata - 1 ) then
    arg = tval - tdata(left)
    first = left - 1
  else if ( left == ndata - 1 ) then
    arg = 0.5D+00 * ( 1.0D+00 + tval - tdata(left) )
    first = left - 1
  end if
!
!  TVEC(I) = ARG^(N-I).
!
  tvec(n) = 1.0D+00
  do i = n-1, 1, -1
    tvec(i) = arg * tvec(i+1)
  end do

  yval = 0.0D+00
  do j = 1, n
    yval = yval + dot_product ( tvec(1:n), mbasis(1:n,j) ) &
      * ydata(first - 1 + j)
  end do

  return
end subroutine

!---------------------------------------------------------------------------------
subroutine basis_matrix_overhauser_nonuni ( alpha, beta, mbasis )

!*****************************************************************************80
!
!! BASIS_MATRIX_OVERHAUSER_NONUNI: nonuniform Overhauser spline basis matrix.
!
!  Discussion:
!
!    This basis matrix assumes that the data points P1, P2, P3 and
!    P4 are not uniformly spaced in T, and that P2 corresponds to T = 0,
!    and P3 to T = 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) ALPHA, BETA.
!    ALPHA = || P2 - P1 || / ( || P3 - P2 || + || P2 - P1 || )
!    BETA  = || P3 - P2 || / ( || P4 - P3 || + || P3 - P2 || ).
!
!    Output, real ( kind = 8 ) MBASIS(4,4), the basis matrix.
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) mbasis(4,4)

  mbasis(1,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(1,2) =   beta + ( 1.0D+00 - alpha ) / alpha
  mbasis(1,3) =   alpha - 1.0D+00 / ( 1.0D+00 - beta )
  mbasis(1,4) =   beta * beta / ( 1.0D+00 - beta )

  mbasis(2,1) =   2.0D+00 * ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(2,2) = ( - 2.0D+00 * ( 1.0D+00 - alpha ) - alpha * beta ) / alpha
  mbasis(2,3) = ( 2.0D+00 * ( 1.0D+00 - alpha ) &
    - beta * ( 1.0D+00 - 2.0D+00 * alpha ) ) / ( 1.0D+00 - beta )
  mbasis(2,4) = - beta * beta / ( 1.0D+00 - beta )

  mbasis(3,1) = - ( 1.0D+00 - alpha ) * ( 1.0D+00 - alpha ) / alpha
  mbasis(3,2) =   ( 1.0D+00 - 2.0D+00 * alpha ) / alpha
  mbasis(3,3) =   alpha
  mbasis(3,4) =   0.0D+00

  mbasis(4,1) =   0.0D+00
  mbasis(4,2) =   1.0D+00
  mbasis(4,3) =   0.0D+00
  mbasis(4,4) =   0.0D+00

  return
end subroutine

end module spline
