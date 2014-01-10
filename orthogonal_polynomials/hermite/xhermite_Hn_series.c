////////////////////////////////////////////////////////////////////////////////
// File: xhermite_Hn_series.c                                                 //
// Routine(s):                                                                //
//    xHermite_Hn_Series                                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xHermite_Hn_Series(long double x, long double a[], int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Hermite           //
//     polynomials, Hn, at a point x,                                         //
//        p(x) = a[0] + a[1]*H[1](x) + a[2]*H[2](x) + ... + a[deg]*H[deg](x). //
//                                                                            //
//     Clenshaw's recursion formula applied to Hermite polynomials Hn is:     //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = 2 * x * y[k+1] - 2 * (k+1) * y[k+2] + a[k].  Then           //
//     p(x) = y[0].                                                           //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial p(x).                 //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Hermite polynomials,  //
//        i.e. a[k] is the coefficient of H[k](x).  Note that in the calling  //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a[N], p;                                                //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = xHermite_Hn_Series(x, a, deg);                                     //
////////////////////////////////////////////////////////////////////////////////
long double xHermite_Hn_Series(long double x, long double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double two_x = x + x;
   int k;
   int two_k = degree + degree + 2;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y, two_k -= 2)
      y = two_x * yp1 - two_k * yp2 + a[k];

   return y; 
}
