////////////////////////////////////////////////////////////////////////////////
// File: xjacobi_Pn_sequence.c                                                //
// Routine(s):                                                                //
//    xJacobi_Pn_Sequence                                                     //
////////////////////////////////////////////////////////////////////////////////

#define Jacobi_Step_First_Term(two_k_gam, two_k_gam_2, x, a2mb2, pn1) \
           ((two_k_gam + 1.0L) * (two_k_gam_2 * two_k_gam * x + a2mb2) * pn1)
#define Jacobi_Step_Second_Term(k, alpha, beta, two_k_gam_2, pn2) \
                        (2.0L * (k + alpha) * (k + beta) * two_k_gam_2 * pn2)
#define Jacobi_Step_Divisor(k, gam, two_k_gam) \
                           (2.0L * (k + 1.0L) * (k + gam + 1.0L) * two_k_gam)

////////////////////////////////////////////////////////////////////////////////
// void xJacobi_Pn_Sequence(long double P[], long double x,                   //
//                           long double alpha, long double beta, int max_n)  //
//                                                                            //
//  Description:                                                              //
//     The Jacobi polynomials with parameters alpha > -1 and beta > -1,       //
//     P^(alpha,beta)n(x), are orthogonal on the interval [-1,1] with         //
//     weight function w(x) = (1-x)^alpha (1+x)^beta and 0 elsewhere.         //
//     For convenience in typing, the superscript (alpha,beta) is dropped so  //
//     that P^(alpha,beta)n(x) will be denoted by P[n](x), i.e. the parameters//
//     alpha and beta are to be understood.  The Jacobi polynomials are norm- //
//     alized so that Pn(1) = gamma(n+alpha+1)/(n! gamma(alpha+1)).           //
//      <Pn,Pm> = 0                                              if n != m,   //
//      <Pn,Pn> = 2^(alpha+beta+1) gamma(n+alpha+1) gamma(n+beta+1)           //
//                / (n! (2n+alpha+beta+1) gamma(n+alpha+beta+1)) if n >= 0.   //
//     These routines calculate, Pn(x), the Jacobi polynomial of degree n     //
//     with parameters alpha and beta of degree n, n = 0,...,max_n, evaluated //
//     at x using the recursion formula:                                      //
//     Set a = alpha, b = beta, and c = alpha + beta,                         //
//     2 (k+1) (k + c + 1) (2k + c) P[k+1](x)                                 //
//             = (2k + c + 1)[(2k + c + 2) (2k + c) x + a^2 - b^2] P[k](x)    //
//               - (k + a) (k + b) (2k + c + 2) P[k-1](x), k = 1,...,n-1.     //
//     P[0](x) = 1, P[1](x) = [ (c + 2) * x + (a - b) ] / 2.                  //
//                                                                            //
//  Arguments:                                                                //
//     long double P[]                                                        //
//        On output, P[n] is the value of the Jacobi polynomial, Pn(),        //
//        evaluated at x where 0 <= n <= max_n.  The calling routine must     //
//        have defined P as long double P[N] where N >= max_n + 1.            //
//     long double x                                                          //
//        The argument of the Jacobi polynomials with parameters alpha and    //
//        beta of degree n, 0 <= n <= max_n.                                  //
//     long double alpha                                                      //
//        The first parameter of the Jacobi polynomial, the exponent of (1-x) //
//        in the weight function.  Note that alpha > -1.                      //
//     long double beta                                                       //
//        The second parameter of the Jacobi polynomial, the exponent of (1+x)//
//        in the weight function.  Note that beta > -1.                       //
//     int    max_n                                                           //
//        The maximum order of the sequence.                                  //
//                                                                            //
//  Return Value:                                                             //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     long double Pn[N];                                                     //
//     long double x;                                                         //
//     long double alpha, beta;                                               //
//     int    max_deg = N - 1;                                                //
//                                                                            //
//     (user code to set x, alpha, and beta)                                  //
//                                                                            //
//     xJacobi_Pn_Sequence(Pn, x, alpha, beta, max_deg);                      //
////////////////////////////////////////////////////////////////////////////////
void xJacobi_Pn_Sequence(long double P[], long double x, long double alpha,
                                                   long double beta, int max_n)
{
   long double pn;
   long double gam = alpha + beta;
   long double two_k_gam;
   long double two_k_gam_2;
   long double a2mb2 = alpha*alpha - beta*beta;
   int k;
                    // Initialize the recursion process. //
 
   if (max_n < 0) return;
   P[0] = 1.0L;
   if (max_n == 0) return;
   P[1] = ((gam + 2.0L) * x + (alpha - beta)) / 2.0L;

                   // Calculate Pn(x) for n = 2,...,max_n //

   for (k = 1; k < max_n; k++) {
      two_k_gam = (long double)(k+k) + gam;
      two_k_gam_2 = two_k_gam + 2.0L;
      pn = Jacobi_Step_First_Term(two_k_gam, two_k_gam_2, x, a2mb2,P[k]);
      pn -= Jacobi_Step_Second_Term(k, alpha, beta, two_k_gam_2, P[k-1]);
      pn /= Jacobi_Step_Divisor(k, gam, two_k_gam);
      P[k+1] = pn;
   }
}
