////////////////////////////////////////////////////////////////////////////////
// File: legendre_shifted_Pn_sequence.c                                       //
// Routine(s):                                                                //
//    Legendre_Shifted_Pn_Sequence                                            //
////////////////////////////////////////////////////////////////////////////////
#define Legendre_Step(divisor, alpha, pn1, gamma, pn2) \
                               ( (alpha * pn1 - gamma * pn2) / divisor )
////////////////////////////////////////////////////////////////////////////////
// void Legendre_Shifted_Pn_Sequence(double P[], double x, int max_n)         //
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
//     double P[]                                                             //
//        On output, P[n] is the value of the shifted Legendre polynomial,    //
//        P*n(), evaluated at x where 0 <= n <= max_n.  The calling routine   //
//        must have be defined P as double P[N] where N >= max_n + 1.         //
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
//     double Pn[N];                                                          //
//     double x;                                                              //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x)                                                   //
//                                                                            //
//     Legendre_Shifted_Pn_Sequence(Pn, x, max_deg);                          //
////////////////////////////////////////////////////////////////////////////////
void Legendre_Shifted_Pn_Sequence(double P[], double x, int max_n)
{
   long double pn, pn1, pn2;
   long double two_x_m1 = (long double) x + (long double) x - 1.0L;
   long double alpha;
   long double gamma;
   long double divisor;
   int k;
   int k1;

                    // Initialize the recursion process. //
 
   pn2 = 1.0L;
   P[0] = (double) pn2;
   if (max_n == 0) return;
   pn1 = (long double) two_x_m1;
   P[1] = (double) pn1;
   if (max_n == 1) return;

                   // Calculate P*n(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      k1 = k + 1;
      divisor = (long double) k1;
      alpha = (long double) (k + k1) * two_x_m1;
      gamma = (long double) k;
      pn = Legendre_Step(divisor, alpha, pn1, gamma, pn2);
      P[k1] = (double)pn;
      pn2 = pn1;
      pn1 = pn;
   }
}
