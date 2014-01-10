////////////////////////////////////////////////////////////////////////////////
// File: legendre_shifted_Pn.c                                                //
// Routine(s):                                                                //
//    Legendre_Shifted_Pn                                                     //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xLegendre_Shifted_Pn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Legendre_Shifted_Pn(double x, int n)                                //
//                                                                            //
//  Description:                                                              //
//     The shifted Legendre polynomials, Pn*(x) = Pn(2x-1), are orthogonal on //
//     the interval [0,1] with weight function w(x) = 1 for 0 <= x <= 1 and 0 //
//     elsewhere.  They are normalized so that Pn*(1) = 1.  The inner products//
//     are:                                                                   //
//             <Pn*,Pm*> = 0        if n != m,                                //
//             <Pn*,Pn*> = 1/(2n+1) if n >= 0.                                //
//     This routine calculates Pn*(x) using the function xLegendre_Shifted_Pn //
//     in the file xlegendre_shifted_Pn which in turn uses the following      //
//     recursion:                                                             //
//     (k+1) P[k+1]*(x) = (2k+1)(2x-1) P[k]*(x) - k P[k-1]*(x), k = 1,...,n-1 //
//              P[0]*(x) = 1, P[1]*(x) = 2x-1.                                //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the shifted Legendre polynomial Pn*.                //
//     int    n                                                               //
//        The degree of the shifted Legendre polynomial Pn*.                  //
//                                                                            //
//  Return Value:                                                             //
//     Pn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//     If Pn*(x) > DBL_MAX, then DBL_MAX is returned and if Pn*(x) < -DBL_MAX //
//     then -DBL_MAX is returned (this applies only if x < 0 or x > 1.)       //
//                                                                            //
//  Example:                                                                  //
//     double Pn;                                                             //
//     double x;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Pn = Legendre_Shifted_Pn(x, n);                                        //
////////////////////////////////////////////////////////////////////////////////
double Legendre_Shifted_Pn(double x, int n)
{
   long double Pn;

   if (n < 0) return 0.0;
   Pn = xLegendre_Shifted_Pn((long double)x, n);
   if (fabsl(Pn) < DBL_MAX) return (double) Pn;
   return (Pn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
