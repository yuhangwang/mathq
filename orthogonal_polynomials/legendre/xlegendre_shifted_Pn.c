////////////////////////////////////////////////////////////////////////////////
// File: xlegendre_shifted_Pn.c                                               //
// Routine(s):                                                                //
//    xLegendre_Shifted_Pn                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xLegendre_Shifted_Pn(long double x, int n)                     //
//                                                                            //
//  Description:                                                              //
//     The shifted Legendre polynomials, Pn*(x) = Pn(2x-1), are orthogonal on //
//     the interval [0,1] with weight function w(x) = 1 for 0 <= x <= 1 and 0 //
//     elsewhere.  They are normalized so that Pn*(1) = 1.  The inner products//
//     are:                                                                   //
//             <Pn*,Pm*> = 0        if n != m,                                //
//             <Pn*,Pn*> = 1/(2n+1) if n >= 0.                                //
//     This routine calculates Pn*(x) using the following recursion:          //
//     (k+1) P[k+1]*(x) = (2k+1)(2x-1) P[k]*(x) - k P[k-1]*(x), k = 1,...,n-1 //
//              P[0]*(x) = 1, P[1]*(x) = 2x-1.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The argument of the shifted Legendre polynomial Pn*.                //
//     int    n                                                               //
//        The degree of the shifted Legendre polynomial Pn*.                  //
//                                                                            //
//  Return Value:                                                             //
//     Pn*(x) if n is a nonnegative integer.  If n is negative, 0 is returned.//
//                                                                            //
//  Example:                                                                  //
//     long double Pn;                                                        //
//     long double x;                                                         //
//     int    n;                                                              //
//                                                                            //
//     (user code to set x and n)                                             //
//                                                                            //
//     Pn = xLegendre_Shifted_Pn(x, n);                                       //
////////////////////////////////////////////////////////////////////////////////
long double xLegendre_Shifted_Pn(long double x, int n)
{
   long double two_x_m1 = x + x - 1.0L;
   long double P0, P1, Pn;
   int k;

   if (n < 0) return 0.0L;
   if ( x == 1.0L ) return 1.0L;
   if ( x == 0.0L ) return (n % 2 == 0) ? 1.0L : -1.0L;
    
                    // Initialize the recursion process. //
 
   if (n == 0) return 1.0L;
   if (n == 1) return two_x_m1;
   P0 = 1.0L;
   P1 = two_x_m1;

                             // Calculate P*n(x) //

   for (k = 1; k < n; k++) {
      Pn = ((long double) (k + k + 1) * two_x_m1 * P1 - (long double) k * P0)
                                                       / (long double)(k + 1);
      P0 = P1;
      P1 = Pn;
   }
   return Pn;
}
