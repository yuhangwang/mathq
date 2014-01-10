////////////////////////////////////////////////////////////////////////////////
// File: xlaguerre_Ln.c                                                       //
// Routine(s):                                                                //
//    xLaguerre_Ln                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h> 

////////////////////////////////////////////////////////////////////////////////
// long double xLaguerre_Ln(long double x, int n)                             //
//                                                                            //
//  Description:                                                              //
//     The Laguerre polynomials, Ln,  are orthogonal on the interval [0,inf)  //
//     with weight function w(x) = exp(-x) on [0, inf) and 0 elsewhere.  The  //
//     Laguerre polynomials are normalized so that Ln(0) = 1.                 //
//             <Ln,Lm> = 0    if n != m,                                      //
//             <Ln,Ln> = 1    if n >= 0.                                      //
//     This routine calculates Ln(x) using the recursion formula:             //
//       (k+1) L[k+1](x) = (2k+1-x) L[k](x) - k L[k-1](x), k = 1,...,n-1      //
//              L[0](x) = 1, L[1](x) = 1-x.                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Laguerre polynomial Ln.                         //
//     int    n                                                               //
//        The degree of the Laguerre polynomial Ln.                           //
//                                                                            //
//  Return Value:                                                             //
//     Ln(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Ln;                                                        //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Ln = xLaguerre_Ln(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
long double xLaguerre_Ln(long double x, int n)
{
   long double one_mx = 1.0L - x;
   long double L0, L1, Ln;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 0.0L ) return 1.0L;
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return one_mx;
   L0 = 1.0L;
   L1 = one_mx;

                             // Calculate Ln(x) //

   for (k = 1; k < n; k++) {
      Ln = (((long double)(k+k) + one_mx) * L1 - (long double)k * L0)
                                                      / (long double)(k + 1);
      L0 = L1;
      L1 = Ln;
   }

   return Ln;
}
