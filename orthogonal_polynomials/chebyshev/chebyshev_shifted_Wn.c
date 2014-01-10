////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Wn.c                                               //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Wn                                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>        // required for fabsl()                             
#include <float.h>       // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xChebyshev_Shifted_Wn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Shifted_Wn(double x, int n)                               //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the fourth kind, Wn*(x)=Wn(2x-1), //
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = sqrt((1-x)/x) on [0,1] and 0 elsewhere.  Further Wn*(x) is      //
//     normalized so that W*n(1) = 2n+1.                                      //
//             <Wn*,Wm*> = 0    if n != m,                                    //
//             <Wn*,Wn*> = pi/2 if n >= 0,                                    //
//     This routine uses xChebyhev_Shift_Wn to calculate Wn*(x) which in turn //
//     uses the following techniques:                                         //
//     If n < N or if |2x-1| > 1, where N is defined in xchebyshev_shift_Wn.c,//
//     then use the recursion                                                 //
//              W*[k+1](x) = 2(2x - 1) W*[k](x) - W*[k-1](x), k = 1,...,n-1   //
//              W*[0](x) = 1, W*[1](x) = 4x - 1.                              //
//     otherwise use the explicit formula                                     //
//      W*[n](x) = sin((n+1/2)*theta)/sin(theta/2), where theta = acos(2x-1). //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomial Wn*.               //
//     int    n                                                               //
//        The degree of the shifted Chebyshev polynomial.                     //
//                                                                            //
//  Return Value:                                                             //
//     Wn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//     If Wn*(x) > DBL_MAX, then DBL_MAX is returned and if Wn*(x) < -DBL_MAX //
//     then -DBL_MAX is returned (this applies only if x < 0 or x > 1.)       //
//                                                                            //
//  Example:                                                                  //
//     double x, Wn;                                                          //
//     double x;                                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Wn = Chebyshev_Shifted_Wn(x, n);                                       //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Shifted_Wn(double x, int n)
{
   long double Wn;

   if (n < 0) return 0.0;
   Wn = xChebyshev_Shifted_Wn((long double)x, n);
   if (fabsl(Wn) < DBL_MAX) return (double) Wn;
   return (Wn > 0.0L) ? DBL_MAX : -DBL_MAX;

}
