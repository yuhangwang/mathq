////////////////////////////////////////////////////////////////////////////////
// File: xhermite_Hen_sequence.c                                              //
// Routine(s):                                                                //
//    xHermite_Hen_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void xHermite_Hen_Sequence(long double He[], long double x, int max_n)     //
//                                                                            //
//  Description:                                                              //
//     There are several systems of orthogonal polynomials called Hermite     //
//     polynomials.  The version programmed here is the system with weight    //
//     function w(x) = exp(-x^2 / 2), denoted by He_n, and normalized so that //
//     the coefficient of the leading term of He_n is 1.  The inner product is//
//     given by:                                                              //
//             <He_n,He_m> = 0             if n != m,                         //
//             <He_n,He_n> = sqrt(2pi) n!  if n >= 0.                         //
//     This routine calculates He[n](x) using the recursion formula:          //
//            He[k+1](x) = x He[k](x) - k He[k-1](x), k = 1,...,n-1           //
//            He[0](x) = 1, He[1](x) = x                                      //
//     to evaluate He_n at x for n = 0, ..., max_n.                           //
//                                                                            //
//  Arguments:                                                                //
//     long double He[]                                                       //
//        On output, He[n] is the value of the Hermite polynomial, He_n(),    //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have be defined Hen as long double He[N] where N >= max_n + 1.      //
//     long double x                                                          //
//        The argument of the Hermite polynomials He_n, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Hen[N];                                                    //
//     long double x;                                                         //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xHermite_Hen_Sequence(Hen, x, max_deg);                                //
////////////////////////////////////////////////////////////////////////////////
void xHermite_Hen_Sequence(long double He[], long double x, int max_n)
{
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return; 
   He[0] = 1.0L;
   if (max_n == 0) return;
   He[1] = x;
   if (max_n == 1) return;

                  // Calculate He_n(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++) 
      He[k] = x * He[k-1] - (long double)(k-1) * He[k-2];
}
