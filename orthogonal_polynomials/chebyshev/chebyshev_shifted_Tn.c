////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Tn.c                                               //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Tn                                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //
long double xChebyshev_Shifted_Tn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Shifted_Tn(double x, int n)                               //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the first kind, Tn*(x) = Tn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = 1 / sqrt(x-x^2) on [0,1] and 0 elsewhere.  Tn* is normalized    //
//     so that T*n(1) = 1.                                                    //
//             <Tn*,Tm*> = 0    if n != m,                                    //
//             <Tn*,Tn*> = pi/2 if n > 0,                                     //
//             <T0*,T0*> = pi.                                                //
//     This routine uses xChebyhev_Shift_Tn to calculate Tn*(x) which in turn //
//     uses the following techniques:                                         //
//     If n < N or if |2x-1| > 1, where N is defined in xchebyshev_shift_Tn.c,//
//     then use the recursion                                                 //
//              T*[k+1](x) = 2(2x - 1) T*[k](x) - T*[k-1](x), k = 1,...,n-1   //
//              T*[0](x) = 1, T*[1](x) = 2x - 1.                              //
//     otherwise use the explicit formula T*[n](x) = cos(n*acos(2x-1))        //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomial Tn*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Tn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//     If Tn*(x) > DBL_MAX, then DBL_MAX is returned and if Tn*(x) < -DBL_MAX //
//     then -DBL_MAX is returned (this applies only if x < 0 or x > 1.)       //
//                                                                            //
//  Example:                                                                  //
//     double x, Tn;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Tn = Chebyshev_Shifted_Tn(x, n);                                       //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Shifted_Tn(double x, int n)
{
   long double Tn;

   if (n < 0) return 0.0;
   Tn = xChebyshev_Shifted_Tn((long double)x, n);
   if (fabsl(Tn) < DBL_MAX) return (double) Tn;
   return (Tn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
