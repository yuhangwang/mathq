////////////////////////////////////////////////////////////////////////////////
// File: xlegendre_Pn_sequence.c                                              //
// Routine(s):                                                                //
//    xLegendre_Pn_Sequence                                                   //
////////////////////////////////////////////////////////////////////////////////
#define Legendre_Step(kp1, two_kp1_x, pn1, k, pn2)  \
                                    ((two_kp1_x * pn1 - k * pn2) / kp1)
////////////////////////////////////////////////////////////////////////////////
// void xLegendre_Pn_Sequence(long double P[], long double x, int max_n)      //
//                                                                            //
//  Description:                                                              //
//     The Legendre polynomials, Pn(x), are orthogonal on the interval [-1,1] //
//     with weight function w(x) = 1 for -1 <= x <= 1 and 0 elsewhere.  They  //
//     are normalized so that Pn(1) = 1.  The inner products are:             //
//             <Pn,Pm> = 0        if n != m,                                  //
//             <Pn,Pn> = 2/(2n+1) if n >= 0.                                  //
//     This routine calculates Pn(x) using the following recursion:           //
//        (k+1) P[k+1](x) = (2k+1)x P[k](x) - k P[k-1](x), k = 1,...,max_n-1  //
//              P[0](x) = 1, P[1](x) = x.                                     //
//                                                                            //
//  Arguments:                                                                //
//     long double P[]                                                        //
//        On output, P[n] is the value of the Legendre polynomial, Pn(),      //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have defined P as double P[N] where N >= max_n + 1.                 //
//     long double x                                                          //
//        The argument of the Legendre polynomials Pn for n = 0,...,max_n.    //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Pn[N];                                                     //
//     long double x;                                                         //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     xLegendre_Pn_Sequence(Pn, x, max_deg);                                 //
////////////////////////////////////////////////////////////////////////////////
void xLegendre_Pn_Sequence(long double P[], long double x, int max_n)
{
   int k;
                    // Initialize the recursion process. //

   if (max_n < 0) return; 
   P[0] = 1.0L;
   if (max_n == 0) return;
   P[1] = x;
   if (max_n == 1) return;

                   // Calculate Pn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) 
      P[k+1] = Legendre_Step((long double)(k+1), (long double)(k+k+1)*x,
                                               P[k], (long double)k, P[k-1]);
}
