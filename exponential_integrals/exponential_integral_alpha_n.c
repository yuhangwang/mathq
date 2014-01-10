////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_alpha_n.c                                       //
// Routine(s):                                                                //
//    Exponential_Integral_Alpha_n                                            //
//    Exponential_Integral_Alpha_n_Sequence                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for expl()
#include <float.h>                            // required for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Alpha_n( double x, int n )                     //
//                                                                            //
//  Description:                                                              //
//     The alpha exponential integral, alpha_n(x), is the integral of         //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from 1 to inf.  The parameter n is called   //
//     the order of the integral.                                             //
//     For x <= 0, alpha_n = inf for all n >= 0.                              //
//     For x > 0, integration by parts yields the following recursive formula //
//            alpha_n(x) = ( n / x ) alpha_(n-1)(x) + alpha_0(x)              //
//     where alpha_0(x) is readily integrated to yield                        //
//                          alpha_0 = exp(-x) / x.                            //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of alpha exponential integral alpha_n().       //
//     int     n  The order of the alpha exponential integral.                //
//                                                                            //
//  Return Value:                                                             //
//     The value of alpha exponential integral of order n evaluated at x.     //
//     If x <= 0, then the integrals diverge and DBL_MAX is returned.         //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     ( code to initialize x and n )                                         //
//                                                                            //
//     y = Exponential_Integral_Alpha_n( x, n);                               //
////////////////////////////////////////////////////////////////////////////////

double Exponential_Integral_Alpha_n( double x, int n )
{
   long double xx = (long double) x;
   long double a0, an;
   long double xm;
   int m;
   int nx;

      // If x <= 0.0, then the integral diverges so we return DBL_MAX. //

   if (x <= 0.0) return DBL_MAX;

                          // Calculate alpha_0(x). //
                         
   a0 = expl(-xx) / xx;
   if (n == 0) return (double) a0;

                          // Calculate alpha_1(x). //
                         
   if ( x < 1.0 ) {
      an = a0 + xx * a0;
      an = 1.0L / an;
      an = xx * an;
      an = 1.0L / an;
   } else  an =  a0 / xx + a0;         
   if ( n == 1) return (double) an;

                    // Calculate alpha_n(x) for n <= x. //

   nx = (int) x;
   if ( n < nx ) nx = n;

   for (m = 2; m <= nx; m++)
      an = ( (long double) m /  xx) * an + a0;

                    // Calculate alpha_n(x) for n > x. //

   for (; m <= n; m++) {
      xm = xx / (long double) m;
      an += xm * a0;
      an = 1.0L / an;
      an = xm * an;
      an = 1.0L / an;
   }

   return (double) an;
}


////////////////////////////////////////////////////////////////////////////////
// void Exponential_Integral_Alpha_n_En_Sequence(double* an, double x, int N )//
//                                                                            //
//  Description:                                                              //
//     This routine returns the values of the alpha exponential integrals     //
//     alpha_n(x) for n = 0,..., N in the array an.                           //
//                                                                            //
//  Arguments:                                                                //
//     double* an On return, an contains the values of alpha exponential      //
//                integrals alpha_n(x), for n = 0, 1, ..., N.                 //
//                The array en must be dimensioned at least N+1 in the        //
//                calling routine.                                            //
//     double  x  The argument of alpha exponential integrals alpha_n().      //
//     int     N  The maximum order of alpha exponential integrals to         //
//                evaluate.                                                   //
//                                                                            //
//  Return Value:                                                             //
//     The function is declared void.  The values of the alpha exponential    //
//     integrals are returned in the array an.                                //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double x;                                                              //
//     double an[N+1];                                                        //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     Exponential_Integral_Alpha_n_Sequence( an, x, N);                      //
////////////////////////////////////////////////////////////////////////////////

void Exponential_Integral_Alpha_n_Sequence( double *an, double x, int N )
{
   long double xx = (long double) x;
   long double a0, am;
   long double xm;
   int m;
   int nx;

      // If x <= 0.0, then the integrals diverge so we return DBL_MAX. //

   if (x <= 0.0) {
      for (m = 0; m <= N; m++) an[m] = DBL_MAX;
      return;
   }

                          // Calculate alpha_0(x). //

   a0 = expl(-xx) / xx;
   an[0] = (double) a0;
   if (N == 0) return;

                          // Calculate alpha_1(x). //
                         
   if ( x < 1.0 ) {
      am = a0 + xx * a0;
      am = 1.0L / am;
      am = xx * am;
      am = 1.0L / am;
   } else  am =  a0 / xx + a0;         
   an[1] = (double) am;
   if ( N == 1) return;

                    // Calculate alpha_n(x) for n <= x. //

   nx = (int) x;
   if ( N < nx ) nx = N;

   for (m = 2; m <= nx; m++) {
      am = ( (long double) m /  xx) * am + a0;
      an[m] = (double) am;
   }

                    // Calculate alpha_n(x) for n > x. //

   for (; m <= N; m++) {
      xm = xx / (long double) m;
      am += xm * a0;
      am = 1.0L / am;
      am = xm * am;
      am = 1.0L / am;
      an[m] = (double) am;
   }

   return;
}
