////////////////////////////////////////////////////////////////////////////////
// File: jacobi_Pn.c                                                          //
// Routine(s):                                                                //
//    Jacobi_Pn                                                               //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>             // required for fabsl()
#include <float.h>            // required for DBL_MAX

//                        Externally Defined Routines                         //

extern long double xJacobi_Pn(long double x, long double alpha,
                                                     long double beta, int n);

////////////////////////////////////////////////////////////////////////////////
// double Jacobi_Pn(double x, double alpha, double beta, int n)               //
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
//     function xJacobi_Pn in the file xjacobi_Pn.c which in turn uses the    //
//     recursion formula:                                                     //
//     Set a = alpha, b = beta, and c = alpha + beta,                         //
//     2 (k+1) (k + c + 1) (2k + c) P[k+1](x)                                 //
//             = (2k + c + 1)[(2k + c + 2) (2k + c) x + a^2 - b^2] P[k](x)    //
//               - (k + a) (k + b) (2k + c + 2) P[k-1](x), k = 1,...,n-1.     //
//     P[0](x) = 1, P[1](x) = [ (c + 2) * x + (a - b) ] / 2.                  //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The argument of the nth-degree Jacobi polynomial with parameters    //
//        alpha and beta.                                                     //
//     double alpha                                                           //
//        The first parameter of the Jacobi polynomial, the exponent of (1-x) //
//        in the weight function.  Note that alpha > -1.                      //
//     double beta                                                            //
//        The second parameter of the Jacobi polynomial, the exponent of (1+x)//
//        in the weight function.  Note that beta > -1.                       //
//     int    n                                                               //
//        The degree of the Jacobi polynomial.                                //
//                                                                            //
//  Return Value:                                                             //
//     Pn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//     If Pn(x) > DBL_MAX, then DBL_MAX is returned and if Pn(x) < -DBL_MAX   //
//     then -DBL_MAX is returned (this applies only if |x| > 1).              //
//                                                                            //
//  Example:                                                                  //
//     double Pn;                                                             //
//     double x;                                                              //
//     double alpha, beta;                                                    //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, beta, and n)                               //
//                                                                            //
//     Pn = Jacobi_Pn(x, alpha, beta, n);                                     //
////////////////////////////////////////////////////////////////////////////////
double Jacobi_Pn(double x, double alpha, double beta, int n)
{
   long double Pn;

   if (n < 0) return 0.0;
   Pn = xJacobi_Pn((long double)x, (long double) alpha, (long double) beta, n);
   if (fabsl(Pn) < DBL_MAX) return (double) Pn;
   return (Pn > 0.0L) ? DBL_MAX : -DBL_MAX;
}
