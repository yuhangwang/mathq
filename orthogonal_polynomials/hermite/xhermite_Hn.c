////////////////////////////////////////////////////////////////////////////////
// File: xhermite_Hn.c                                                        //
// Routine(s):                                                                //
//    xHermite_Hn                                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xHermite_Hn(long double x, int n)                              //
//                                                                            //
//  Description:                                                              //
//     There are several systems of orthogonal polynomials called Hermite     //
//     polynomials.  The version programmed here is the system with weight    //
//     function w(x) = exp(-x^2).  Normalized so that the coefficient of the  //
//     leading term of Hn is 2^n.  The inner product is given by:             //
//             <Hn,Hm> = 0    if n != m,                                      //
//             <Hn,Hn> = n! 2^n sqrt(pi) if n >= 0.                           //
//     This routine calculates H[n](x) using the recursion formula:           //
//            H[k+1](x) = 2x H[k](x) - 2k H[k-1](x), k = 1,...,n-1            //
//            H[0](x) = 1, H[1](x) = 2x                                       //
//     to evaluate Hn at x.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the Hermite polynomial, Hn, of degree n.            //
//     int    n                                                               //
//        The degree of the Hermite polynomial.                               //
//                                                                            //
//  Return Value:                                                             //
//     Hn(x) if n is a nonnegative integer.  If n is negative, 0 is returned. //
//                                                                            //
//  Example:                                                                  //
//     long double Hn;                                                        //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Hn = xHermite_Hn(x, n);                                                //
////////////////////////////////////////////////////////////////////////////////
long double xHermite_Hn(long double x, int n)
{
   long double two_x = x + x;
   long double H0, H1, Hn;
   int two_n = n + n;
   int k;

   if (n < 0) return 0.0L;
      
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x;
   H0 = 1.0L;
   H1 = two_x;

                             // Calculate Hn(x) //

   for (k = 2; k < two_n; k+=2) {
      Hn = two_x * H1 - k * H0;
      H0 = H1;
      H1 = Hn;
   }

   return Hn;
}
