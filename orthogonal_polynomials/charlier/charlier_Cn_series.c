////////////////////////////////////////////////////////////////////////////////
// File: charlier_Cn_series.c                                                 //
// Routine(s):                                                                //
//    Charlier_Cn_Series                                                      //
////////////////////////////////////////////////////////////////////////////////
#define Clenshaw_Charlier_Step(amx, cn1, k, cn2, a, c) \
                       ( ( (k + amx) * cn1  - (k + 1.0L) * cn2 ) / a + c )
////////////////////////////////////////////////////////////////////////////////
// double Charlier_Cn_Series(double x, double a, double c[], int degree)      //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion formula to evaluate a given     //
//     polynomial p(x) expressed as a linear combination of Charlier          //
//     polynomials, C[n](x), at a point x.                                    //
//     The polynomial p(x) is                                                 //
//        p(x) = c[0] + c[1]*C[1](x) + c[2]*C[2](x) + ... + c[deg]*C[deg](x)  //
//     where C[n] is the n-th Charlier polynomial with parameter a.           //
//                                                                            //
//     Clenshaw's recursion formula applied to the Charlier polynomials with  //
//     parameter a is:                                                        //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = ( (k+a-x)/a ) y[k+1] - ( (k+1)/a ) y[k+2] + c[k].           //
//     Then p(x) = y[0].                                                      //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Charlier polynomials with parameter a.          //
//     double a                                                               //
//        The parameter of the Charlier polynomials.                          //
//     double c[]                                                             //
//        The coefficients of the expansion in terms of Charlier polynomials. //
//        i.e. c[k] is the coefficient of C[k](x).  Note that in the calling  //
//        routine c must be defined double c[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     double x, a, c[N], p;                                                  //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, a, and c[i] i = 0, ... , deg )                 //
//                                                                            //
//     p = Charlier_Cn_Series(x, a, c, deg);                                  //
////////////////////////////////////////////////////////////////////////////////

double Charlier_Cn_Series(double x, double a, double c[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double amx = (long double)a - (long double)x;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

                       // Apply Clenshaw's recursion. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
      y = Clenshaw_Charlier_Step(amx, yp1, (long double) k, yp2,
                                          (long double) a, (long double) c[k]);
   }
   return (double) y; 
}
