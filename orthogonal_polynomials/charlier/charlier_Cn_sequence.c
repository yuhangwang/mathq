////////////////////////////////////////////////////////////////////////////////
// File: charlier_Cn_sequence.c                                               //
// Routine(s):                                                                //
//    Charlier_Cn_Sequence                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// void Charlier_Cn_Sequence(double C[], double x, double a, int max_n)       //
//                                                                            //
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
//     double C[]                                                             //
//        On output, C[n] is the value of the Charlier polynomial of degree n //
//        with parameter a, Cn(), evaluated at x where 0 <= n <= max_n.       //
//        The calling routine must have defined C as double C[N] where        //
//        N >= max_n + 1.                                                     //
//     double x                                                               //
//        The argument of the Charlier polynomials with parameter a, Cn       //
//        n = 0,...,max_n.                                                    //
//     double a                                                               //
//        The parameter of the Charlier polynomials, Cn, n = 0,...,max_n.     //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     double Cn[N];                                                          //
//     double x;                                                              //
//     double a;                                                              //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x and a)                                             //
//                                                                            //
//     Charlier_Cn_Sequence(Cn, x, a, max_deg);                               //
////////////////////////////////////////////////////////////////////////////////
void Charlier_Cn_Sequence(double C[], double x, double a, int max_n)
{
   long double cn, cn1, cn2;
   long double amx = (long double)a - (long double)x;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   cn2 = 1.0L;
   C[0] = (double) cn2;
   if (max_n == 0) return;
   cn1 = amx / a;
   C[1] = (double) cn1;
   if (max_n == 1) return;

                   // Calculate Cn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++, cn2 = cn1, cn1 = cn) {
      cn = ((amx + (long double)k) * cn1 - (long double) k * cn2)
                                                           / (long double) a;
      C[k+1] = (double)cn;
   }
}
