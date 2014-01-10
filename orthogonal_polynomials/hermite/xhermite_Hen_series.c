////////////////////////////////////////////////////////////////////////////////
// File: xhermite_Hen_series.c                                                //
// Routine(s):                                                                //
//    xHermite_Hen_Series                                                     //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xHermite_Hen_Series(long double x, long double a[], int degree)//
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Hermite           //
//     polynomials, He_n, at a point x,                                       //
//      p(x) = a[0] + a[1]*He[1](x) + a[2]*He[2](x) + ... + a[deg]*He[deg](x).//
//                                                                            //
//     Clenshaw's recursion formula applied to Hermite polynomials He_n is:   //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = x * y[k+1] - (k + 1) * y[k+2] + a[k].  Then p(x) = y[0].    //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial p(x).                 //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Hermite polynomials,  //
//        i.e. a[k] is the coefficient of He[k](x).  Note that in the calling //
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
//     p = xHermite_Hen_Series(x, a, deg);                                    //
////////////////////////////////////////////////////////////////////////////////
long double xHermite_Hen_Series(long double x, long double a[], int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0L;

           // Apply Clenshaw's recursion save the last iteration. //
 
   for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) 
      y = x * yp1 - (long double)(k+1) * yp2 + a[k];

   return y; 
}
