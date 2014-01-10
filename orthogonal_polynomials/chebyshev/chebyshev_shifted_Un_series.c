////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Un_series_evaluation.c                             //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Un_Series                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Shifted_Un_Series(double x, double a[], int degree)       //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of shifted Chebyshev //
//     polynomials of the second kind, U*n, at a point x,                     //
//     p(x) = a[0] + a[1]*U*[1](x) + a[2]*U*[2](x) + ... + a[deg]*U*[deg](x). //
//                                                                            //
//     Clenshaw's recursion formula applied to Chebyshev shifted polynomials  //
//     of the second kind is:                                                 //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = (4*x - 2) * y[k+1] - y[k+2] + a[k].  Then p(x) = y[0].      //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The point at which to evaluate the polynomial p(x).                 //
//     double a[]                                                             //
//        The coefficients of the expansion in terms of Chebyshev polynomials,//
//        i.e. a[k] is the coefficient of U*[k](x).  Note that in the calling //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     double x, a[N], p;                                                     //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = Chebyshev_Shifted_Un_Series(x, a, deg);                            //
////////////////////////////////////////////////////////////////////////////////

double Chebyshev_Shifted_Un_Series(double x, double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double two_x_m1 = (long double)x + (long double)x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) 
      y = four_x_m2 * yp1 - yp2 + (long double) a[k];

   return (double) y; 
}
