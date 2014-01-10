////////////////////////////////////////////////////////////////////////////////
// File: xgegenbauer_Cn_series.c                                              //
// Routine(s):                                                                //
//    xGegenbauer_Cn_Series                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xGegenbauer_Cn_Series(long double x, long double alpha,        //
//                                               long double a[], int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Gegenbauer        //
//     polynomials, Cn, with common parameter alpha evaluated at x,           //
//        p(x) = a[0] + a[1]*C[1](x) + a[2]*C[2](x) + ... + a[deg]*C[deg](x)  //
//     where C[n] is the n-th Gegenbauer polynomial with parameter alpha.     //
//                                                                            //
//     Clenshaw's recursion formula applied to Gegenbauer polynomials is:     //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = 2(k+alpha)*x*y[k+1]/(k + 1) - (k+2alpha) y[k+2])/(k+2)+a[k].//
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial p(x).                 //
//     long double alpha                                                      //
//        The common parameter of the Gegenbauer polynomials.                 //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Gegenbauer            //
//        polynomials, i.e. a[k] is the coefficient of C[k](x).  Note that in //
//        the calling routine a must be defined long double a[N] where        //
//        N >= degree + 1.                                                    //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, alpha, a[N], p;                                         //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, alpha, and a[i] i = 0, ... , a[deg] )          //
//                                                                            //
//     p = xGegenbauer_Cn_Series(x, alpha, a, deg);                           //
////////////////////////////////////////////////////////////////////////////////
long double xGegenbauer_Cn_Series(long double x, long double alpha,
                                                   long double a[], int degree)
{
   long double two_alpha = 2.0L * alpha;
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double kp1;
   long double kp2;
   long double k_two_alpha;
   long double beta;
   long double gamma;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      kp2 = (long double) (k + 2);
      k_two_alpha = (long double)k + two_alpha;
      beta = ( ( (long double)k + k_two_alpha ) * x ) / kp1;
      gamma = k_two_alpha / kp2;
      y = beta * yp1 - gamma * yp2 + a[k];
   }

   return y; 
}
