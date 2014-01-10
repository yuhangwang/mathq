////////////////////////////////////////////////////////////////////////////////
// File: xkrawtchouk_Kn_series.c                                              //
// Routine(s):                                                                //
//    xKrawtchouk_Kn_Series_Evaluation                                        //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Krawtchouk_Step(a, kn1, g, kn2, kp1, c) \
                     ((( a * kn1 ) / kp1  - ( g * kn2 ) / ( kp1 + 1.0L ) )+ c )
////////////////////////////////////////////////////////////////////////////////
// long double Krawtchouk_Kn_Series_Evaluation(long double x, long double p,  //
//                                       int N, long double c[], int degree)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion formula to evaluate a given     //
//     polynomial p(x) expressed as a linear combination of Krawtchouk        //
//     polynomials, K[n](x), at a point x where n <= N, N being the maximum   //
//     of the support points {0,1,...,N} of the weight function w(x) being    //
//     the binomial distribution Bin(N,p).                                    //
//     The polynomial p(x) is                                                 //
//        p(x) = c[0] + c[1]*K[1](x) + c[2]*K[2](x) + ... + c[deg]*K[deg](x)  //
//     where K[n] is the n-th Krawtchouk polynomial with deg <= N.            //
//                                                                            //
//     Clenshaw's recursion formula applied to the Krawtchouk polynomials     //
//     is:                                                                    //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = {[x - (k + p*(N-2k))] / (k+1)} y[k+1]                       //
//                                - [(N - k)*p*(1-p) / (k+1)] y[k+2]] + c[k]. //
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Krawtchouk polynomials with weight function     //
//        function support points {0,1,2,...,N} and values Bin(N,p,x).        //
//     long double p                                                          //
//        A parameter of the Krawtchouk polynomial, 0 < p < 1.                //
//     int    N                                                               //
//        N + 1 is the number of equally spaced support points of the weight  //
//     long double c[]                                                        //
//        The coefficients of the expansion in terms of Krawtchouk polynomials//
//        i.e. c[k] is the coefficient of K[k](x).  Note that in the calling  //
//        routine c must be defined double c[L] where L >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).  Note that degree <= N.          //
//        If degree > N, degree is reset internally to N.                     //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a, c[K], p, Kn;                                         //
//     int    deg = K - 1;                                                    //
//     int    N;                                                              //
//                                                                            //
//     ( code to initialize x, N, p, and c[i] i = 0, ... , deg )              //
//                                                                            //
//     Kn = xKrawtchouk_Kn_Series(x, p, N, c, deg);                           //
////////////////////////////////////////////////////////////////////////////////
long double xKrawtchouk_Kn_Series(long double x, long double p, int N,
                                                   long double c[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double a;
   long double g;
   long double pq = p * ( 1.0L - p);
   long double kp1;
   int k;

          // Check that degree <= N.  If not, then set degree = N. //

   if ( degree > N ) degree = N;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

                       // Apply Clenshaw's recursion. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      a = ( x - (long double) k - p * (long double)( N - k - k) );
      g = (long double)(N - k) * pq;
      y = Clenshaw_Krawtchouk_Step( a, yp1, g, yp2, kp1, c[k]);
   }

   return y; 
}
