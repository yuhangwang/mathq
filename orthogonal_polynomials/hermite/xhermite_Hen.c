////////////////////////////////////////////////////////////////////////////////
// File: xhermite_Hen.c                                                       //
// Routine(s):                                                                //
//    xHermite_Hen                                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xHermite_Hen(long double x, int n)                             //
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
//     to evaluate He_n at x.                                                 //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Hermite polynomial, He_n, of degree n.          //
//     int    n                                                               //
//        The degree of the Hermite polynomial.                               //
//                                                                            //
//  Return Value:                                                             //
//     He_n(x) if n is a nonnegative integer.  If n is negative, then         //
//     0 is returned.                                                         //
//                                                                            //
//  Example:                                                                  //
//     long double Hen;                                                       //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Hen = xHermite_Hen(x, n);                                              //
////////////////////////////////////////////////////////////////////////////////
long double xHermite_Hen(long double x, int n)
{
   long double He0, He1, Hen;
   int k;

   if (n < 0) return 0.0L;
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return x;
   He0 = 1.0L;
   He1 = x;

                            // Calculate He_n(x) //

   for (k = 1; k < n; k++) {
      Hen = x * He1 - k * He0;
      He0 = He1;
      He1 = Hen;
   }

   return Hen;
}
