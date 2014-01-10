////////////////////////////////////////////////////////////////////////////////
// File: hermite_Hen_sequence.c                                               //
// Routine(s):                                                                //
//    Hermite_Hen_Sequence                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Hermite_Hen_Sequence(double He[], double x, int max_n)                //
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
//     double He[]                                                            //
//        On output, He[n] is the value of the Hermite polynomial, He_n(),    //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have be defined Hen as double He[N] where N >= max_n + 1.           //
//     double x                                                               //
//        The argument of the Hermite polynomials He_n, n = 0,...,max_n.      //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double Hen[N];                                                         //
//     double x;                                                              //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Hermite_Hen_Sequence(Hen, x, max_deg);                                 //
////////////////////////////////////////////////////////////////////////////////
void Hermite_Hen_Sequence(double He[], double x, int max_n)
{
   long double hen, hen1, hen2;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   hen2 = 1.0L;
   hen1 = (long double) x;
   He[0] = (double) hen2;
   if (max_n == 0) return;
   He[1] = x;
   if (max_n == 1) return;

                   // Calculate He_n(x) for n = 2,...,max_n //

   for (k = 2; k <= max_n; k++, hen2 = hen1, hen1 = hen) {
      hen = (long double) x * hen1 - (long double)(k-1) * hen2;
      He[k] = (double)hen;
   }
}
