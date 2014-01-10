////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Un_series.c                                               //
// Routine(s):                                                                //
//    xChebyshev_Un_Series                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Un_Series(long double x, long double a[],int degree)//
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Chebyshev         //
//     polynomials of the second kind, Un, at a point x,                      //
//        p(x) = a[0] + a[1]*U[1](x) + a[2]*U[2](x) + ... + a[deg]*U[deg](x). //
//                                                                            //
//     Clenshaw's recursion formula applied to Chebyshev polynomials of the   //
//     second kind is:                                                        //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = 2 * x * y[k+1] - y[k+2] + a[k].  Then p(x) = y[0].          //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial p(x).                 //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Chebyshev polynomials,//
//        i.e. a[k] is the coefficient of U[k](x).  Note that in the calling  //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial at p(x).                                   //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a[N], p;                                                //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = xChebyshev_Un_Series(x, a, deg);                                   //
////////////////////////////////////////////////////////////////////////////////

long double xChebyshev_Un_Series(long double x, long double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double two_x = x + x;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

                  // Apply Clenshaw's recursion algorithm. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) 
      y = two_x * yp1 - yp2 + a[k];

   return y; 
}
