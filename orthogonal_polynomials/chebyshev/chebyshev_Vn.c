////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Vn.c                                                       //
// Routine(s):                                                                //
//    Chebyshev_Vn                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>         // required for fabsl() 
#include <float.h>        // required for DBL_MAX

//                        Externally Defined Routines                         //

long double xChebyshev_Vn(long double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double Chebyshev_Vn(double x, int n)                                       //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the third kind, Vn(x) are orthogonal on   //
//     the interval [-1,1] with weight function w(x) = sqrt((1+x)/(1-x)) on   //
//     interval [-1,1] and 0 elsewhere and are normalized so that Vn(1) = 1.  //
//             <Vn,Vm> = 0    if n != m,                                      //
//             <Vn,Vn> = pi   if n >= 0.                                      //
//     This routine uses xChebyhev_Vn to calculate Vn(x) which in turn uses   //
//     the following techniques:                                              //
//     If n < N or if |x| > 1, where N is defined in xchebyshev_Vn.c, then    //
//     use the recursion                                                      //
//              V[k+1](x) = 2x V[k](x) - V[k-1](x), k = 1,2,...,n-1           //
//              V[0](x) = 1, V[1](x) = 2x-1.                                  //
//     otherwise use V[n](x) = cos((n+1/2)*theta)/cos(theta/2), where         //
//     theta = acos(x).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the Chebyshev polynomial Vn.                        //
//     int    n                                                               //
//        The degree of the Chebyshev polynomial.                             //
//                                                                            //
//  Return Value:                                                             //
//     Vn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Vn(x) > DBL_MAX, then DBL_MAX is returned and if Vn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned (this applies only if |x| > 1).              //
//                                                                            //
//  Example:                                                                  //
//     double x, Vn;                                                          //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Vn = Chebyshev_Vn(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
double Chebyshev_Vn(double x, int n)
{
   long double Vn;

   if (n < 0) return 0.0;
   Vn = xChebyshev_Vn((long double)x, n);
   if (fabsl(Vn) < DBL_MAX) return (double) Vn;
   return (Vn > 0.0L) ? DBL_MAX : -DBL_MAX;

}
