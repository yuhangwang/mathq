////////////////////////////////////////////////////////////////////////////////
// File: xjacobi_Pn.c                                                         //
// Routine(s):                                                                //
//    xJacobi_Pn                                                              //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xJacobi_Pn(long double x, long double alpha, long double beta, //
//                                                                     int n) //
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
//     This routine calculates, Pn(x), the Jacobi polynomial of degree n      //
//     with parameters alpha and beta of degree n evaluated at x using the    //
//     recursion formula:                                                     //
//     Set a = alpha, b = beta, and c = alpha + beta,                         //
//     2 (k+1) (k + c + 1) (2k + c) P[k+1](x)                                 //
//             = (2k + c + 1)[(2k + c + 2) (2k + c) x + a^2 - b^2] P[k](x)    //
//               - (k + a) (k + b) (2k + c + 2) P[k-1](x), k = 1,...,n-1.     //
//     P[0](x) = 1, P[1](x) = [ (c + 2) * x + (a - b) ] / 2.                  //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the nth-degree Jacobi polynomial with parameters    //
//        alpha and beta.                                                     //
//     long double alpha                                                      //
//        The first parameter of the Jacobi polynomial, the exponent of (1-x) //
//        in the weight function.  Note that alpha > -1.                      //
//     long double beta                                                       //
//        The second parameter of the Jacobi polynomial, the exponent of (1+x)//
//        in the weight function.  Note that beta > -1.                       //
//     int    n                                                               //
//        The degree of the Jacobi polynomial.                                //
//                                                                            //
//  Return Value:                                                             //
//     Pn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Pn;                                                        //
//     long double x;                                                         //
//     long double alpha, beta;                                               //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, beta, and n)                               //
//                                                                            //
//     Pn = xJacobi_Pn(x, alpha, beta, n);                                    //
////////////////////////////////////////////////////////////////////////////////
long double xJacobi_Pn(long double x, long double alpha, long double beta,  
                                                                         int n)
{
   long double P0, P1, Pn;
   long double gam = alpha + beta;
   long double two_k_gam;
   long double two_k_gam_2;
   long double a2mb2 = alpha*alpha - beta*beta;
   int k;

   if (n < 0) return 0.0L;
      
                    // Initialize the recursion process. //
 
   P0 = 1.0L;
   if (n == 0) return P0;
   P1 = ((gam + 2.0L) * x + (alpha - beta)) / 2.0L;
   if (n == 1) return P1;

                             // Calculate Pn(x) //

   for (k = 1; k < n; k++) {
      two_k_gam = (long double)(k+k) + gam;
      two_k_gam_2 = two_k_gam + 2.0L;
      Pn = (two_k_gam + 1.0L) * (two_k_gam_2 * two_k_gam * x + a2mb2) * P1;

      Pn -= 2.0L*((long double)k+alpha)*((long double)k+beta)*(two_k_gam_2)*P0;
      Pn /= (2.0L *(long double)(k+1)*((long double)(k+1) + gam) * two_k_gam);
      P0 = P1;
      P1 = Pn;
   }

   return Pn;
}
