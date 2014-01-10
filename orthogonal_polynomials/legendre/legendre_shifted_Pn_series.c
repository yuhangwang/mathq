////////////////////////////////////////////////////////////////////////////////
// File: legendre_shifted_Pn_series.c                                         //
// Routine(s):                                                                //
//    Legendre_Shifted_Pn_Series                                              //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Legendre_Step(alpha, pn1, gamma, pn2, a) \
                                               (alpha * pn1 - gamma * pn2 + a)
////////////////////////////////////////////////////////////////////////////////
// double Legendre_Shifted_Pn_Series(double x, double a[], int degree)        //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of shifted Legendre  //
//     polynomials, P*n, at a point x,                                        //
//           p(x) = a[0] + a[1]*P*[1](x) + ... + a[deg]*P*[deg](x).           //
//                                                                            //
//     Clenshaw's recursion formula applied to shifted Legendre polynomials   //
//     is:                                                                    //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = (2k+1)*(2x-1)*y[k+1]/(k+1) - (k+1)*y[k+2]/(k+2) + a[k].     //
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The point at which to evaluate the polynomial p(x).                 //
//     double a[]                                                             //
//        The coefficients of the expansion in terms of Legendre shifted      //
//        polynomials i.e. a[k] is the coefficient of P*[k](x).  Note that in //
//        the calling routine a must be defined double a[N]                   //
//        where N >= degree + 1.                                              //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p at x.                                    //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     double x, a[N], p;                                                     //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = Legendre_Shifted_Pn_Series(x, a, deg);                             //
////////////////////////////////////////////////////////////////////////////////

double Legendre_Shifted_Pn_Series(double x, double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double two_x_m1 = (long double)x + (long double)x - 1.0L;
   long double kp1;
   long double kp2;
   long double alpha;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      kp2 = (long double) (k + 2);
      alpha = (kp1 + (long double) k) * two_x_m1 / kp1;
      y = Clenshaw_Legendre_Step(alpha, yp1, kp1/kp2, yp2, (long double) a[k]);
   }

   return (double) y; 
}
