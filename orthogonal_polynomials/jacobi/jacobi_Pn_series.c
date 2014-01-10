////////////////////////////////////////////////////////////////////////////////
// File: jacobi_Pn_series.c                                                   //
// Routine(s):                                                                //
//    Jacobi_Pn_Series                                                        //
////////////////////////////////////////////////////////////////////////////////

#define Clenshaw_Jacobi_Step(b, pn1, c, pn2, a) (b * pn1 - c * pn2 + a)

////////////////////////////////////////////////////////////////////////////////
// double Jacobi_Pn_Series(double x, double alpha, double beta, double a[],   //
//                                                                int degree) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Jacobi            //
//     polynomials, Pn, with common parameters alpha > -1 and beta > -1       //
//     evaluated at x,                                                        //
//        p(x) = a[0] + a[1]*P[1](x) + a[2]*P[2](x) + ... + a[deg]*P[deg](x). //
//                                                                            //
//     Clenshaw's recursion formula applied to Jacobi polynomials is:         //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 0  //
//     set y[k] = b[k] * (x - c[k]) * y[k+1] - d[k] * y[k+2]) + a[k], where   //
//     after setting g = alpha + beta,                                        //
//        b[k] = (2k + g + 1) * (2k + g + 2) / [2 * (k + 1) * (k + g + 1)],   //
//        c[k] = g * (beta - alpha) / [(2k + g) * (2k + g + 2)],              //
//        d[k] = (k + alpha + 1) * (k + beta + 1) * (2k + g + 4) /            //
//                           [(k + 2) * (k + g + 2) * (2k + g + 2)].          //
//      Then p(x) = y[0].                                                     //
//                                                                            //
//  Arguments:                                                                //
//     double x                                                               //
//        The point at which to evaluate the polynomial p(x).                 //
//     double alpha                                                           //
//        The first parameter of the Jacobi polynomial, the exponent of (1-x) //
//        in the weight function.  Note that alpha > -1.                      //
//     double beta                                                            //
//        The second parameter of the Jacobi polynomial, the exponent of (1+x)//
//        in the weight function.  Note that beta > -1.                       //
//     double a[]                                                             //
//        The coefficients of the expansion in terms of Jacobi polynomials,   //
//        i.e. a[k] is the coefficient of P[k](x).  Note that in the calling  //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial p(x).                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     double x, alpha, beta, a[N], p;                                        //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, alpha, beta, and a[i] i = 0, ... , a[deg] )    //
//                                                                            //
//     p = Jacobi_Pn_Series(x, alpha, beta, a, deg);                          //
////////////////////////////////////////////////////////////////////////////////
double Jacobi_Pn_Series(double x, double alpha, double beta, double a[],
                                                                    int degree)
{
   long double yp2 = 0.0L;
   long double yp1 = 0.0L;
   long double y = 0.0L;
   long double gamma = (long double) alpha + (long double) beta;
   long double xx = (long double) x;
   long double b;
   long double c;
   long double kp1;
   long double kp2;
   long double a2mb2;
   long double two_k;
   long double two_k_p1;
   long double two_k_p2;
   long double two_k_p4;
   int k;

             // Check that degree >= 0.  If not, then return 0. //

   if ( degree < 0 ) return 0.0;

                       // Apply Clenshaw's recursion. //

   if (gamma == 0.0L) {
      for (k = degree; k >= 1; k--, yp2 = yp1, yp1 = y) {
         kp1 = (long double)(k + 1);
         two_k = (long double)(k + k);
         b = ( (two_k + 1.0L) * xx ) / kp1;
         c = ((kp1 + alpha) * (kp1 + beta)) / (kp1 * (kp1 + 1.0L));
         y = Clenshaw_Jacobi_Step(b, yp1, c, yp2, (long double) a[k]);
      }
      b = x + ( (long double)alpha - (long double) beta) / 2.0L;
      c = (alpha + 1.0L) * (beta + 1.0L) / 2.0L;
      y = Clenshaw_Jacobi_Step(b, yp1, c, yp2, a[k]);
   }
   else if (gamma == -1.0L) {
      for (k = degree; k >= 1; k--, yp2 = yp1, yp1 = y) {
         kp1 = (long double)(k + 1);
         two_k = (long double)(k + k);
         b = (two_k + 1.0L) * xx;
         b -= ((long double)alpha - (long double) beta) / (two_k - 1.0L);
         b /= kp1;
         c = (kp1 + alpha) * (kp1 + beta) * ( two_k + 3.0L);
         c /= (kp1 * (kp1 + 1.0L) * ( two_k + 1.0L ) );
         y = Clenshaw_Jacobi_Step(b, yp1, c, yp2, (long double) a[k]);
      }
      b = (x + (long double)alpha - (long double) beta) / 2.0L;
      c = 3.0L * (1.0L + (long double) alpha) * (1.0L + (long double) beta)
                                                                     / 2.0L;
      y = Clenshaw_Jacobi_Step(b, yp1, c, yp2, (long double) a[k]);
   }
 else {
      a2mb2 = alpha * alpha - beta * beta;
      for (k = degree; k >= 0; k--, yp2 = yp1, yp1 = y) {
         kp1 = (long double) (k + 1);
         kp2 = (long double) (k + 2);
         two_k = (long double) (k + k);
         two_k_p1 = two_k + 1.0L;
         two_k_p2 = two_k + 2.0L;
         two_k_p4 = two_k + 4.0L;
         b = (two_k_p1 + gamma)*((two_k_p2 + gamma)*xx + a2mb2/(two_k + gamma));
         b /=  ( 2.0L * kp1 * (kp1 + gamma) );
         c = (kp1 + (long double) alpha) * (kp1 + (long double) beta)
                                                          * (two_k_p4 + gamma);
         c /= ( kp2 * (kp2 + gamma) * (two_k_p2 + gamma) );
         y = Clenshaw_Jacobi_Step(b, yp1, c, yp2, (long double) a[k]);
      }
   }
   return y; 
}
