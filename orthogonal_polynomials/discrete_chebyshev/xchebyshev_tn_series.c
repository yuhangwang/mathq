////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_tn_series.c                                                //
// Routine(s):                                                                //
//    xChebyshev_tn_Series                                                    //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Chebyshev_Step(a, tn1, g, tn2, kp1, c) \
                     ((( a * tn1 ) / kp1  - ( g * tn2 ) / ( kp1 + 1.0L ) )+ c )
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_tn_Series(long double x, int N, long  double c[],   //
//                                                                int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion formula to evaluate a given     //
//     polynomial p(x) expressed as a linear combination of discrete          //
//     Chebyshev polynomials, t[n](x), at a point x where n < N, N being      //
//     the support points {0,1,...,N-1} of the weight function.               //
//     The polynomial p(x) is                                                 //
//        p(x) = c[0] + c[1]*t[1](x) + c[2]*t[2](x) + ... + c[deg]*t[deg](x)  //
//     where t[n] is the n-th discrete Chebyshev polynomial with n < N.       //
//                                                                            //
//     Clenshaw's recursion formula applied to the discrete Chebyshev         //
//     polynomials is:                                                        //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = [(2k+1)(2x - (N-1)) y[k+1] - (k+1)(N^2 - (k+1)^2 ) y[k+2]]  //
//                                                            / (k+1) + c[k]. //
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the discrete Chebyshev polynomials with weight      //
//        function support points {0,1,2,...,N-1}.                            //
//     int    N                                                               //
//        The number of equally spaced support points of the weight function. //
//     long double c[]                                                        //
//        The coefficients of the expansion in terms of discrete Chebyshev    //
//        polynomials i.e. c[k] is the coefficient of t[k](x).  Note that in  //
//        the calling routine c must be defined double c[K] where             //
//        K >= degree + 1.                                                    //
//     int    degree                                                          //
//        The degree of the polynomial p(x).  Note that degree <= N-1.        //
//        If degree > N - 1, degree is reset internally to N-1.               //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a, c[K], p;                                             //
//     int    deg = K - 1;                                                    //
//     int    N;                                                              //
//                                                                            //
//     ( code to initialize x, N, and c[i] i = 0, ... , deg )                 //
//     ( where K < N. )                                                       //
//                                                                            //
//     p = xChebyshev_tn_Series(x, N, c, deg);                                //
////////////////////////////////////////////////////////////////////////////////
long double xChebyshev_tn_Series(long double x, int N, long double c[],
                                                                    int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double a = 2.0L * x - (long double)(N-1);
   long double kp1;
   long double N2 = (long double)(N * N);
   int k;

          // Check that degree < N.  If not, then set degree = N-1. //

   if ( degree >= N ) degree = N - 1;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

                       // Apply Clenshaw's recursion. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      y = Clenshaw_Chebyshev_Step(( a * (long double)(k + k + 1) ), yp1,
                                    (kp1 * (N2 - kp1 * kp1)), yp2, kp1, c[k]);
   }

   return y; 
}
