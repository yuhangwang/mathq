////////////////////////////////////////////////////////////////////////////////
// File: xcharlier_Cn_sequence.c                                              //
// Routine(s):                                                                //
//    xCharlier_Cn_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xCharlier_Cn_Sequence(long double C[], long double x, long double a,  //
//                                                                 int max_n) //
//  Description:                                                              //
//     The Charlier polynomials are orthogonal on the interval [0,inf)        //
//     with Poisson weight function w(x) = Sum [( (exp(-a) a^j / j! ) D(x-j)],//
//     where a > 0 is the mean of the Poisson distribution, D() is the Dirac  //
//     delta function, and the sum is over j = 0, 1, ... .                    //
//     This routine calculates, C[n](x), the n-th Charlier polynomial         //
//     evaluated at x for n = 0, ,,,, max_n using the following recursion     //
//     formula:                                                               //
//               a C[k+1](x) = (k + a - x) C[k](x) - k C[k-1]                 //
//               C[0](x) = 1, C[1](x) = 1 - x/a.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double C[]                                                        //
//        On output, C[n] is the value of the Charlier polynomial with        //
//        parameter a, Cn(), evaluated at x where 0 <= n <= max_n.            //
//        The calling routine must have defined C as long double C[N] where   //
//        N >= max_n + 1.                                                     //
//     long double x                                                          //
//        The argument of the Charlier polynomials with parameter a, Cn       //
//        n = 0,...,max_n.                                                    //
//     long double a                                                          //
//        The parameter of the Charlier polynomials, Cn, n = 0,...,max_n.     //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Cn[N];                                                     //
//     long double x;                                                         //
//     long double a;                                                         //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x and a)                                             //
//                                                                            //
//     xCharlier_Cn_Sequence(Cn, x, a, max_deg);                              //
////////////////////////////////////////////////////////////////////////////////
void xCharlier_Cn_Sequence(long double C[], long double x, long double a,
                                                                     int max_n)
{
   long double amx = a - x;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   C[0] = 1.0L;
   if (max_n == 0) return;
   C[1] = amx / a;
   if (max_n == 1) return;

                   // Calculate Cn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++)
      C[k+1] = ((amx + (long double)k) * C[k] - (long double) k * C[k-1]) / a;
}
