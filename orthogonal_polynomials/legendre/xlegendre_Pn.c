////////////////////////////////////////////////////////////////////////////////
// File: xlegendre_Pn.c                                                       //
// Routine(s):                                                                //
//    xLegendre_Pn                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h> 

////////////////////////////////////////////////////////////////////////////////
// long double xLegendre_Pn(long double x, int n)                             //
//                                                                            //
//  Description:                                                              //
//     The Legendre polynomials, Pn(x), are orthogonal on the interval [-1,1] //
//     with weight function w(x) = 1 for -1 <= x <= 1 and 0 elsewhere.  They  //
//     are normalized so that Pn(1) = 1.  The inner products are:             //
//             <Pn,Pm> = 0        if n != m,                                  //
//             <Pn,Pn> = 2/(2n+1) if n >= 0.                                  //
//     This routine calculates Pn(x) using the following recursion:           //
//        (k+1) P[k+1](x) = (2k+1)x P[k](x) - k P[k-1](x), k = 1,...,n-1      //
//              P[0](x) = 1, P[1](x) = x.                                     //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Legendre polynomial Pn.                         //
//     int    n                                                               //
//        The degree of the Legendre polynomial Pn.                           //
//                                                                            //
//  Return Value:                                                             //
//     Pn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Pn;                                                        //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Pn = xLegendre_Pn(x, n);                                               //
////////////////////////////////////////////////////////////////////////////////
long double xLegendre_Pn(long double x, int n)
{
   long double P0, P1, Pn;
   int k;

   if (n < 0) return 0.0L;
   if ( fabsl(x) == 1.0L ) { 
      if (x > 0.0L) return 1.0L;
      else if (n % 2 == 0) return 1.0L;
      else return -1.0L;
   }
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return x;
   P0 = 1.0L;
   P1 = x;

                             // Calculate Pn(x) //

   for (k = 1; k < n; k++, P0 = P1, P1 = Pn) 
      Pn = ( (long double)( k + k + 1) * x * P1 - (long double)k * P0)
                                                       / (long double)(k + 1);

   return Pn;
}
