////////////////////////////////////////////////////////////////////////////////
// File: xcharlier_Cn.c                                                       //
// Routine(s):                                                                //
//    xCharlier_Cn                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xCharlier_Cn(long double x, long double a, int n)              //
//                                                                            //
//  Description:                                                              //
//     The Charlier polynomials are orthogonal on the interval [0,inf)        //
//     with Poisson weight function w(x) = Sum [( (exp(-a) a^j / j! ) D(x-j)],//
//     where a > 0 is the mean of the Poisson distribution, D() is the Dirac  //
//     delta function, and the sum is over j = 0, 1, ... .                    //
//     This routine calculates, C[n](x), the n-th Charlier polynomial         //
//     evaluated at x using the following recursion formula:                  //
//               a C[k+1](x) = (k + a - x) C[k](x) - k C[k-1]                 //
//               C[0](x) = 1, C[1](x) = 1 - x/a.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Charlier polynomial with parameter a.           //
//     long double a                                                          //
//        The parameter of the Charlier polynomial, a > 0.                    //
//     int    n                                                               //
//        The degree of the Charlier polynomial with parameter a.             //
//                                                                            //
//  Return Value:                                                             //
//     Cn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Cn;                                                        //
//     long double x;                                                         //
//     long double a;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, a, and n)                                         //
//                                                                            //
//     Cn = xCharlier_Cn( x, a, n);                                           //
////////////////////////////////////////////////////////////////////////////////
long double xCharlier_Cn(long double x, long double a, int n)
{
   long double amx = a - x;
   long double C0, C1, Cn;
   int k;
                    // Initialize the recursion process. //
 
   if (n < 0) return 0.0L;
   C0 = 1.0L;
   if (n == 0) return C0;
   C1 = amx / a;
   if (n == 1) return C1;

                             // Calculate Cn(x) //

   for (k = 1; k < n; k++) {
      Cn = ( ( (long double)k + amx ) * C1 - (long double)k * C0 ) / a;
      C0 = C1;
      C1 = Cn;
   }

   return Cn;
}
