////////////////////////////////////////////////////////////////////////////////
// File: laguerre_Ln_sequence.c                                               //
// Routine(s):                                                                //
//    Laguerre_Ln_Sequence                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Laguerre_Ln_Sequence(double L[], double x, int max_n)                 //
//                                                                            //
//  Description:                                                              //
//     The Laguerre polynomials, Ln,  are orthogonal on the interval [0,inf)  //
//     with weight function w(x) = exp(-x) on [0, inf) and 0 elsewhere.  The  //
//     Laguerre polynomials are normalized so that Ln(0) = 1.                 //
//             <Ln,Lm> = 0    if n != m,                                      //
//             <Ln,Ln> = 1    if n >= 0.                                      //
//     This routine calculates, L[n](x), Laguerre polynomials evaluated at    //
//     x for n = 0, ,,,, max_n using the recursion formula:                   //
//       (k+1) L[k+1](x) = (2k+1-x) L[k](x) - k L[k-1](x), k = 1,...,max_n-1  //
//              L[0](x) = 1, L[1](x) = 1-x.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double L[]                                                             //
//        On output, L[n] is the value of the Laguerre polynomial, Ln(),      //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        must have defined L as double L[N] where N >= max_n + 1.            //
//     double x                                                               //
//        The argument of the Laguerre polynomials Ln, n = 0,...,max_n.       //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double x,Ln[N];                                                        //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Laguerre_Ln_Sequence(Ln, x, max_deg);                                  //
////////////////////////////////////////////////////////////////////////////////
void Laguerre_Ln_Sequence(double L[], double x, int max_n)
{
   long double ln, ln1, ln2;
   long double one_mx = 1.0L - (long double) x;
   long double alpha, gamma;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return; 
   ln2 = 1.0L;
   ln1 = one_mx;
   L[0] = (double) ln2;
   if (max_n == 0) return;
   L[1] = (double) ln1;
   if (max_n == 1) return;

                   // Calculate Ln(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++, ln2 = ln1, ln1 = ln) {
      alpha = ( (long double)(k + k) + one_mx );
      gamma = (long double) k;
      ln = (alpha * ln1 -  gamma * ln2) / (long double)(k+1);
      L[k+1] = (double)ln;
   }
}
