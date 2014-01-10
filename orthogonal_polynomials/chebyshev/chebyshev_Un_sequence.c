////////////////////////////////////////////////////////////////////////////////
// File: chebyshev_Un_sequence.c                                              //
// Routine(s):                                                                //
//    Chebyshev_Un_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Chebyshev_Un_Sequence(double U[], double x, int max_n)                //
//                                                                            //
//  Description:                                                              //
//     The Chebyshev polynomials of the second kind, Un(x), are orthogonal on //
//     the interval [-1,1] with weight function w(x) = sqrt(1-x^2) on the     //
//     interval [-1,1] and 0 elsewhere and are normalized so that Un(1) = n+1.//
//             <Un,Um> = 0    if n != m,                                      //
//             <Un,Un> = pi/2 if n >= 0.                                      //
//     This routine calculates, U[n](x), the Chebyshev polynomials of the     //
//     second kind evaluated at x for n = 0, ..., max_n using the recursion   //
//              U[n+1](x) = 2x U[n](x) - U[n-1](x), n = 1,...,max_n - 1       //
//              U[0](x) = 1, U[1](x) = 2x.                                    //
//                                                                            //
//  Arguments:                                                                //
//     double U[]                                                             //
//        On output, U[n] is the value of the Chebyshev polynomial, Un(),     //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must have//
//        defined U as double U[N] where N >= max_n + 1.                      //
//     double x                                                               //
//        The argument of the Chebyshev polynomials Un, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double x, Un[N];                                                       //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Chebyshev_Un_Sequence(Un, x, max_deg);                                 //
////////////////////////////////////////////////////////////////////////////////
void Chebyshev_Un_Sequence(double U[], double x, int max_n)
{
   long double un, un1, un2;
   long double two_x = (long double) x + (long double) x;
   int k;
                    // Initialize the recursion process. //
   
   if (max_n < 0) return; 
   un2 = 1.0L;
   U[0] = (double) un2;
   if (max_n == 0) return;
   un1 = two_x;
   U[1] = (double) un1;
   if (max_n == 1) return;

                   // Calculate Un(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, un2 = un1, un1 = un) {
      un = two_x * un1 - un2;
      U[k] = (double)un;
   }
}
