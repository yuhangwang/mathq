////////////////////////////////////////////////////////////////////////////////
// File: xlegendre_shifted_Pn_sequence.c                                      //
// Routine(s):                                                                //
//    xLegendre_Shifted_Pn_Sequence                                           //
////////////////////////////////////////////////////////////////////////////////
#define Legendre_Step(divisor, alpha, pn1, gamma, pn2) \
                               ( (alpha * pn1 - gamma * pn2) / divisor )
////////////////////////////////////////////////////////////////////////////////
// void xLegendre_Shifted_Pn_Sequence(long double P[], long double x,         //
//                                                                int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The shifted Legendre polynomials, Pn*(x) = Pn(2x-1), are orthogonal on //
//     the interval [0,1] with weight function w(x) = 1 for 0 <= x <= 1 and 0 //
//     elsewhere.  They are normalized so that Pn*(1) = 1.  The inner products//
//     are:                                                                   //
//             <Pn*,Pm*> = 0        if n != m,                                //
//             <Pn*,Pn*> = 1/(2n+1) if n >= 0.                                //
//     This routine calculates Pn*(x) using the following recursion:          //
//     (k+1) P[k+1]*(x) = (2k+1)(2x-1) P[k]*(x) - k P[k-1]*(x),               //
//              P[0]*(x) = 1, P[1]*(x) = 2x-1, k = 1,...,max_n - 1.           //
//                                                                            //
//  Arguments:                                                                //
//     long double P[]                                                        //
//        On output, P[n] is the value of the shifted Legendre polynomial,    //
//        P*n(), evaluated at x where 0 <= n <= max_n.  The calling routine   //
//        must have be defined P as long double P[N] where N >= max_n + 1.    //
//     double x                                                               //
//        The argument of the shifted Legendre polynomials P*n, n = 0, ...,   //
//        max_n.                                                              //
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
//     xLegendre_Shifted_Pn_Sequence(Pn, x, max_deg);                         //
////////////////////////////////////////////////////////////////////////////////
void xLegendre_Shifted_Pn_Sequence(long double P[], long double x, int max_n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double alpha;
   long double gamma;
   long double divisor;
   int k;
   int k1;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   P[0] = 1.0L;
   if (max_n == 0) return;
   P[1] = two_x_m1;
   if (max_n == 1) return;

                   // Calculate P*n(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      k1 = k + 1;
      divisor = (long double) k1;
      alpha = (long double) (k + k1) * two_x_m1;
      gamma = (long double) k;
      P[k1] = Legendre_Step(divisor, alpha, P[k], gamma, P[k-1]);
   }
}
