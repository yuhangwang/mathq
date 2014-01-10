////////////////////////////////////////////////////////////////////////////////
// File: xlegendre_Pn_series.c                                                //
// Routine(s):                                                                //
//    xLegendre_Pn_Series                                                     //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Legendre_Step(alpha, pn1, gamma, pn2, a) \
                                               (alpha * pn1 - gamma * pn2 + a)
////////////////////////////////////////////////////////////////////////////////
// long double xLegendre_Pn_Series(long double x, long double a[],int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Legendre          //
//     polynomials, Pn, at a point x,                                         //
//        p(x) = a[0] + a[1]*P[1](x) + a[2]*P[2](x) + ... + a[deg]*P[deg](x). //
//                                                                            //
//     Clenshaw's recursion formula applied to Legendre polynomials is:       //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = (2k+1)*x*y[k+1]/(k+1) - (k+1)*y[k+2]/(k+2) + a[k].  Then    //
//     p(x) = y[0].                                                           //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial p(x).                 //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Legendre polynomials,  /
//        i.e. a[k] is the coefficient of P[k](x).  Note that in the calling  //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p at x.                                    //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a[N], p;                                                //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = xLegendre_Pn_Series(x, a, deg);                                    //
////////////////////////////////////////////////////////////////////////////////
long double xLegendre_Pn_Series(long double x, long double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double kp1, kp2;
   int k;
             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      kp1 = (long double) (k + 1);
      kp2 = (long double) (k + 2);
      y = Clenshaw_Legendre_Step( (kp1 + (long double)k) * x / kp1, yp1,
                                                       kp1 / kp2, yp2, a[k]);
   }
   return y; 
}
