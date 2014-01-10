////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Vn.c                                               //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Vn                                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>         // required for fabsl()                              
#include <float.h>        // required for DBL_MAX

//                        Externally Defined Routines                         //
extern long double xChebyshev_Shifted_Vn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Shifted_Vn(double x, int n)                               //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the third kind, Vn*(x) = Vn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt(x/(1-x)) on [0,1] and 0 elsewhere.  Further Vn* is         //
//     normalized so that Vn*(1) = 1.                                         //
//             <Vn*,Vm*> = 0    if n != m,                                    //
//             <Vn*,Vn*> = pi/2 if n >= 0.                                    //
//     This routine uses xChebyhev_Shift_Vn to calculate Vn*(x) which in turn //
//     uses the following techniques:                                         //
//     If n < N or if |2x-1| > 1, where N is defined in xchebyshev_shift_Vn.c,//
//     then use the recursion                                                 //
//              V*[k+1](x) = 2(2x - 1) V*[k](x) - V*[k-1](x), k = 1,...,n-1   //
//              V*[0](x) = 1, V*[1](x) = 4x - 3.                              //
//     otherwise use V*[n](x) = cos((n+1/2)*theta)/cos(theta/2),              //
//     where theta = acos(2x-1).                                              //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomials Vn*.              //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Vn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//     If Vn*(x) > DBL_MAX, then DBL_MAX is returned and if Vn*(x) < -DBL_MAX //
//     then -DBL_MAX is returned (this applies only if x < 0 or x > 1.)       //
//                                                                            //
//  Example:                                                                  //
//     double x, Vn;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Vn = Chebyshev_Shifted_Vn(x, n);                                       //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Shifted_Vn(double x, int n)
{
   long double Vn;

   if (n < 0) return 0.0;
   Vn = xChebyshev_Shifted_Vn((long double)x, n);
   if (fabsl(Vn) < DBL_MAX) return (double) Vn;
   return (Vn > 0.0L) ? DBL_MAX : -DBL_MAX;

}
