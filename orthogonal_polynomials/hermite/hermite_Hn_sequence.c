////////////////////////////////////////////////////////////////////////////////
// File: hermite_Hn_sequence.c                                                //
// Routine(s):                                                                //
//    Hermite_Hn_Sequence                                                     //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Hermite_Hn_Sequence(double H[], double x, int max_n)                  //
//                                                                            //
//  Description:                                                              //
//     There are several systems of orthogonal polynomials called Hermite     //
//     polynomials.  The version programmed here is the system with weight    //
//     function w(x) = exp(-x^2).  Normalized so that the coefficient of the  //
//     leading term of Hn is 2^n.  The inner product is given by:             //
//             <Hn,Hm> = 0    if n != m,                                      //
//             <Hn,Hn> = n! 2^n sqrt(pi) if n >= 0.                           //
//     This routine calculates H[n](x) using the recursion formula:           //
//            H[k+1](x) = 2x H[k](x) - 2k H[k-1](x), k = 1,...,n-1            //
//            H[0](x) = 1, H[1](x) = 2x                                       //
//     to evaluate Hn at x for n = 0, ..., max_n.                             //
//                                                                            //
//  Arguments:                                                                //
//     double H[]                                                             //
//        On output, H[n] is the value of the Hermite polynomial, Hn(),       //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have be defined H as double H[N] where N >= max_n + 1.              //
//     double x                                                               //
//        The argument of the Hermite polynomials Hn, n = 0,...,max_n.        //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double Hn[N];                                                          //
//     double x;                                                              //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Hermite_Hn_Sequence(Hn, x, max_deg);                                   //
////////////////////////////////////////////////////////////////////////////////
void Hermite_Hn_Sequence(double H[], double x, int max_n)
{
   long double hn, hn1, hn2;
   long double two_x = (long double) x + (long double) x;
   int k;
   int two_k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   hn2 = 1.0L;
   hn1 = two_x;
   H[0] = (double) hn2;
   if (max_n == 0) return;
   H[1] = (double) hn1;
   if (max_n == 1) return;

                   // Calculate Hn(x) for n = 2,...,max_n //

   for (k = 2, two_k = 2; k <= max_n; k++, two_k += 2, hn2 = hn1, hn1 = hn) {
      hn = two_x * hn1 - (long double)two_k * hn2;
      H[k] = (double)hn;
   }
}
