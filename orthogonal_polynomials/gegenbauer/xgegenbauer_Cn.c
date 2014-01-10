////////////////////////////////////////////////////////////////////////////////
// File: xgegenbauer_Cn.c                                                     //
// Routine(s):                                                                //
//    xGegenbauer_Cn                                                          //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xGegenbauer_Cn(long double x, long alpha, int n)               //
//                                                                            //
//  Description:                                                              //
//     The Gegenbauer polynomials with parameter alpha > -1/2, C^(alpha)n(x)  //
//     are orthogonal on the interval [-1,1] with weight function             //
//     w(x) = (1-x^2)^(alpha-1/2) on [-1,1] and 0 elsewhere.  Further they are//
//     normalized so that C^(alpha)n(1) = gamma(n+2alpha)/(n! gamma(2*alpha)).//
//     Gegenbauer polynomials are also called the ultraspherical polynomials. //
//     For convenience in typing, the superscript (alpha) is dropped so that  //
//     C^(alpha)n(x) will be denoted by Cn(x), i.e. the parameter alpha is    //
//     to be understood.                                                      //
//      <Cn,Cm> = 0                                      if n != m,           //
//      <Cn,Cn> = 2^(1-2alpha) gamma(n+2alpha) pi                             //
//                / (n! (n+alpha)gamma^2(alpha) )        if n >= 0.           //
//     This routine calculates, Cn(x), the Gegenbauer polynomial of degree n  //
//     with parameter alpha evaluated at x  using recursion formula:          //
//     (k+1) C[k+1](x) = 2(k + alpha) x C[k](x) - (k + 2 alpha - 1) C[k-1](x),//
//                                                         k = 1,...,n-1.     //
//     C[0](x) = 1, C[1](x) = 2*alpha*x.                                      //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Gegenbauer polynomial with parameter alpha of   //
//        degree n.                                                           //
//     long double alpha                                                      //
//        The parameter of the Gegenbauer polynomial.                         //
//     int    n                                                               //
//        The degree of the Gegenbauer polynomial.                            //
//                                                                            //
//  Return Value:                                                             //
//     Cn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double x, alpha, Cn;                                              //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x, alpha, and n)                                     //
//                                                                            //
//     Cn = xGegenbauer_Cn(x, alpha, n);                                      //
////////////////////////////////////////////////////////////////////////////////
long double xGegenbauer_Cn(long double x, long double alpha, int n)
{
   long double C0, C1, Cn;
   int k;

   if (n < 0) return 0.0L;
      
                    // Initialize the recursion process. //
 
   C0 = 1.0L;
   if (n == 0) return C0;
   C1 = 2.0L * alpha * x;
   if (n == 1) return C1;

                             // Calculate Cn(x) //

   for (k = 1; k < n; k++) {
      Cn = (2.0L*((long double) k + alpha) * x * C1
                       - ((long double)k + 2.0L*alpha - 1.0L) * C0)
                                                       / (long double)(k + 1);
      C0 = C1;
      C1 = Cn;
   }

   return Cn;
}
