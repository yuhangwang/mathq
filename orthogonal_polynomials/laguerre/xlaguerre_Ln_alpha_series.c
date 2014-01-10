////////////////////////////////////////////////////////////////////////////////
// File: xlaguerre_Ln_alpha_series.c                                          //
// Routine(s):                                                                //
//    xLaguerre_Ln_alpha_Series                                               //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Laguerre_Step(beta, ln1, kp1, akp1, kp2, ln2, a) \
                              ( (beta * ln1) / kp1 - (akp1 * ln2) / kp2 + a )
////////////////////////////////////////////////////////////////////////////////
// long double xLaguerre_Ln_alpha_Series(long double x, long double alpha,    //
//                                               long double a[], int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of generalized       //
//     Laguerre polynomials, L^(alpha)[n](x), at a point x.                   //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     L^(alpha)[n](x) will be denoted by L[n](x), i.e. the alpha is to be    //
//     understood.                                                            //
//     The polynomial p(x) is                                                 //
//        p(x) = a[0] + a[1]*L[1](x) + a[2]*L[2](x) + ... + a[deg]*L[deg](x)  //
//     where L[n] is the n-th generalized Laguerre polynomial with parameter  //
//     alpha.                                                                 //
//                                                                            //
//     Clenshaw's recursion formula applied to generalized Laguerre           //
//     polynomials with parameter alpha is:                                   //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = ((2k+alpha+1-x))/(k+1) y[k+1] - ((k+alpha+1)/(k+2)) y[k+2]  //
//                                                                    + a[k]. //
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the generalized Laguerre polynomials with parameter //
//        alpha.                                                              //
//     double alpha                                                           //
//        The parameter of the generalized Laguerre polynomials.              //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of generalized Laguerre  //
//        polynomials, i.e. a[k] is the coefficient of L[k](x).  Note that in //
//        the calling routine a must be defined long double a[N] where        //
//        N >= degree + 1.                                                    //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p at x.                                    //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, alpha, a[N], p;                                         //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, alpha, and a[i] i = 0, ... , a[deg] )          //
//                                                                            //
//     p = xLaguerre_Ln_alpha_Series(x, alpha, a, deg);                       //
////////////////////////////////////////////////////////////////////////////////

long double xLaguerre_Ln_alpha_Series(long double x, long double alpha,
                                                   long double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double beta;
   long double kp1;
   long double akp1;
   long double kp2;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

                       // Apply Clenshaw's recursion. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      akp1 = alpha + kp1;
      kp2 = (long double) (k + 2);
      beta = ((long double)(k + k + 1) + alpha - x);
      y = Clenshaw_Laguerre_Step(beta, yp1, kp1, akp1, kp2, yp2, a[k]);
   }
   return y; 
}
