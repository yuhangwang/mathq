////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_shifted_Tn_sequence.c                                      //
// Routine(s):                                                                //
//    Chebyshev_Shifted_Tn_Sequence                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_Shifted_Tn_Sequence(double T[], double x, int max_n)        //
//                                                                            //
//  Description:                                                              //
//     The shifted Chebyshev polynomials of the first kind, Tn*(x) = Tn(2x-1),//
//     are orthogonal on the interval [0,1] with weight function              //
//     w(x) = 1 / sqrt(x-x^2) on [0,1] and 0 elsewhere.  Tn* is normalized    //
//     so that T*n(1) = 1.                                                    //
//             <Tn*,Tm*> = 0    if n != m,                                    //
//             <Tn*,Tn*> = pi/2 if n > 0,                                     //
//             <T0*,T0*> = pi.                                                //
//     This routine calculates T*[n](x), the shifted Chebyshev polynomials    //
//     of the first kind evaluated at x for n = 0, ..., max_n.                //
//     This procedure uses the recursion:                                     //
//         T*[n+1](x) = 2(2x - 1) T*[n](x) - T*[n-1](x), n = 1,...,max_n - 1  //
//              T*[0](x) = 1, T*[1](x) = 2x - 1.                              //
//                                                                            //
//  Arguments:                                                                //
//     double T[]                                                             //
//        On output, T[n] is the value of the shifted Chebyshev polynomial,   //
//        Tn*(), evaluated at x, where 0 <= n <= max_n.  The calling routine  //
//        must have defined T as double T[N] where N >= max_n + 1.            //
//     double x                                                               //
//        The argument of the shifted Chebyshev polynomials Tn*, n = 0,...,   //
//        max_n.                                                              //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double x, Tn[N];                                                       //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Chebyshev_Shifted_Tn_Sequence(Tn, x, max_deg);                         //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_Shifted_Tn_Sequence(double T[], double x, int max_n)
{
   long double tn, tn1, tn2;
   long double two_x_m1 = (long double) x + (long double) x - 1.0L;
   long double four_x_m2 = two_x_m1 + two_x_m1;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   tn2 = 1.0L;
   T[0] = (double) tn2;
   if (max_n == 0) return;
   tn1 = (long double) two_x_m1;
   T[1] = (double) tn1;
   if (max_n == 1) return;

                   // Calculate Tn*(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, tn2 = tn1, tn1 = tn) {
      tn = four_x_m2 * tn1 - tn2;
      T[k] = (double)tn;
   }
}
