////////////////////////////////////////////////////////////////////////////////
// File: exponential_integral_beta_n.c                                        //
// Routine(s):                                                                //
//    Exponential_Integral_Beta_n                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                 // required for fabs(), sinhl() and coshl()
#include <float.h>                // required for LDBL_EPSILON 

static long double Beta_n_Recursion( double x, int n );
static long double Beta_n_Power_Series( double x, int n);
static long double Beta_n_Power_Series_Even_n( double x, int n );
static long double Beta_n_Power_Series_Odd_n( double x, int n );

////////////////////////////////////////////////////////////////////////////////
// double Exponential_Integral_Beta_n( double x, int n )                      //
//                                                                            //
//  Description:                                                              //
//     The beta exponential integral, beta_n(x), is the integral of           //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from -1 to 1.  The parameter n is called    //
//     the order of the integral.                                             //
//     Integration by parts yields the following recursive formula            //
//        beta_n(x) = ( n / x ) beta_(n-1)(x) - (2/x) cosh(x), if n is odd    //
//        beta_n(x) = ( n / x ) beta_(n-1)(x) + (2/x) sinh(x), if n is even.  //
//     where beta_0(x) is readily integrated                                  //
//                         beta_0 = 2 sinh(x) / x.                            //
//     The i-th deriviative of beta_n(x) is (-1)^i beta_(n+i)(x).  The        //
//     power series expansion of beta_n(x) about x = 0.0 is:                  //
//     For n even,                                                            //
//         beta_n(x) = 2/(n+1) + 2 x^2/[(2)!(n+3)] + ... +                    //
//                                           2 x^(2i)/[(2i)!(n+2i+1)] + ...   //
//     and for n odd,                                                         //
//         beta_n(x) = -2 x/(n+2) + 2 x^3/[(3)!(n+4)] + ... +                 //
//                                       2 x^(2i+1)/[(2i+1)!(n+2i+2)] + ...   //
//                                                                            //
//     For 0 < n <= 4, the recursive formula is used for |x| > .1 and the     //
//                     power series representation for |x| <= 0.1.            //
//     For n = 5,      the recursive formula is used for |x| > .4 and the     //
//                     power series representation for |x| <= 0.4.            //
//     For n = 6,      the recursive formula is used for |x| > .6 and the     //
//                     power series representation for |x| <= 0.6.            //
//     For n = 7,      the recursive formula is used for |x| > 1.0 and the    //
//                     power series representation for |x| <= 1.0.            //
//     For n = 8,      the recursive formula is used for |x| > 1.4 and the    //
//                     power series representation for |x| <= 1.4.            //
//     For n = 9,      the recursive formula is used for |x| > 1.7 and the    //
//                     power series representation for |x| <= 1.7.            //
//     For n = 10,     the recursive formula is used for |x| > 2.4 and the    //
//                     power series representation for |x| <= 2.4.            //
//     For n >= 11,    the power series representation is used for all x.     //
//     Convergence can be slow for large |x|.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of beta exponential integral beta_n().         //
//     int     n  The order of the beta exponential integral.                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of beta exponential integral of order n evaluated at x.      //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     ( code to initialize x and n )                                         //
//                                                                            //
//     y = Exponential_Integral_Beta_n( x, n);                                //
////////////////////////////////////////////////////////////////////////////////

double Exponential_Integral_Beta_n( double x, int n )
{
   long double xx = (long double) x;
   long double bn;
   long double c;
   long double s;
   int m = 2;

   if (x == 0.0) {
      if ( n % 2 == 0 ) return 2.0 / (double)(n+1);
      else return 0.0;
   }

   s = 2.0L * sinhl(x) / xx;
   c = -2.0L * coshl(x) / xx;
 
   if (n == 0) return (double) s;

   if (n <= 4)
      if ( fabs(x) <= 0.1 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 5)
      if ( fabs(x) <= 0.4 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 6)
      if ( fabs(x) <= 0.6 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 7)
      if ( fabs(x) <= 1.0 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 8)
      if ( fabs(x) <= 1.4 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 9)
      if ( fabs(x) <= 1.7 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else if (n <= 10)
      if ( fabs(x) <= 2.4 ) bn = Beta_n_Power_Series( x, n );
      else bn = Beta_n_Recursion( x, n );
   else bn = Beta_n_Power_Series( x, n );

   return (double) bn;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Beta_n_Recursion( double x, int n )                     //
//                                                                            //
//  Description:                                                              //
//     The beta exponential integral, beta_n(x), is the integral of           //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from -1 to 1.  The parameter n is called    //
//     the order of the integral.                                             //
//     Integration by parts yields the following recursive formula            //
//        beta_n(x) = ( n / x ) beta_(n-1)(x) - (2/x) cosh(x), if n is odd    //
//        beta_n(x) = ( n / x ) beta_(n-1)(x) + (2/x) sinh(x), if n is even.  //
//     where beta_0(x) is readily integrated                                  //
//                         beta_0 = 2 sinh(x) / x.                            //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of beta exponential integral beta_n().         //
//     int     n  The order of the beta exponential integral.                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of beta exponential integral of order n evaluated at x.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static long double Beta_n_Recursion( double x, int n )
{
   long double xx = (long double) x;
   long double bn;
   long double c;
   long double s;
   int m = 2;

   if (x == 0.0) {
      if ( n % 2 == 0 ) return 2.0L / (long double)(n+1);
      else return 0.0L;
   }

   s = 2.0L * sinhl(x) / xx;
   c = -2.0L * coshl(x) / xx;
 
   if (n == 0) return (double) s;

   bn = s;
   for (m = 1; m <= n; m++) {
      if ( x > m )  bn = ( (long double) m /  xx) * bn + c;
      else {
         bn += ( xx / (long double) m ) * c;
         bn = 1.0L / bn;
         bn *= (xx / (long double) m );
         bn = 1.0L / bn;
      }
      m++;
      if (m > n) break;
      if ( x > m )  bn = ( (long double) m /  xx) * bn + s;
      else {
         bn += ( xx / (long double) m ) * s;
         bn = 1.0L / bn;
         bn *= (xx / (long double) m );
         bn = 1.0L / bn;
      }
   }
      
   return bn;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Beta_n_Power_Series( double x, int n )                  //
//                                                                            //
//  Description:                                                              //
//     The beta exponential integral, beta_n(x), is the integral of           //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from -1 to 1.  The parameter n is called    //
//     the order of the integral.                                             //
//     The i-th deriviative of beta_n(x) is (-1)^i beta_(n+i)(x).  The        //
//     power series expansion of beta_n(x) about x = 0 is:                    //
//     For n even,                                                            //
//         beta_n(x) = 2/(n+1) + 2 x^2/[(2)!(n+3)] + ... +                    //
//                                           2 x^(2i)/[(2i)!(n+2i+1)] + ...   //
//     and for n odd,                                                         //
//         beta_n(x) = -2 x/(n+2) + 2 x^3/[(3)!(n+4)] + ... +                 //
//                                       2 x^(2i+1)/[(2i+1)!(n+2i+2)] + ...   //
//                                                                            //
//     Convergence can be slow for large |x|.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of beta exponential integral beta_n().         //
//     int     n  The order of the beta exponential integral.                 //
//                                                                            //
//  Return Value:                                                             //
//     The value of beta exponential integral of order n evaluated at x.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static long double Beta_n_Power_Series( double x, int n)
{
   if ( n % 2 == 0) return Beta_n_Power_Series_Even_n( x, n );
   return Beta_n_Power_Series_Odd_n( x, n );
}

  
////////////////////////////////////////////////////////////////////////////////
// static long double Beta_n_Power_Series_Even_n( double x, int n )           //
//                                                                            //
//  Description:                                                              //
//     The beta exponential integral, beta_n(x), is the integral of           //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from -1 to 1.  The parameter n is called    //
//     the order of the integral.                                             //
//     The i-th deriviative of beta_n(x) is (-1)^i beta_(n+i)(x).  The        //
//     power series expansion of beta_n(x) about x = 0 is:                    //
//     For n even,                                                            //
//         beta_n(x) = 2/(n+1) + 2 x^2/[(2)!(n+3)] + ... +                    //
//                                           2 x^(2i)/[(2i)!(n+2i+1)] + ...   //
//     Convergence can be slow for large |x|.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of beta exponential integral beta_n().         //
//     int     n  The order of the beta exponential integral n must be even.  //
//                                                                            //
//  Return Value:                                                             //
//     The value of beta exponential integral of order n evaluated at x.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static long double Beta_n_Power_Series_Even_n( double x, int n )
{ 
   long double xx = (long double) x;
   long double x2 = xx * xx;
   long double xn = 1.0L;
   long double np1 = (long double) (n + 1);
   long double Sn = xn / np1;
   long double Sm1 = 0.0L;
   long double factorial = 1.0L;
   long double y = 0.0L;
   
   while ( fabsl(Sn - Sm1) > LDBL_EPSILON * fabsl(Sm1) ) {
      Sm1 = Sn;
      y += 1.0L;
      factorial *= y;
      y += 1.0L;
      factorial *= y;
      xn *= x2;
      Sn += ( xn / (factorial * (np1 + y)) );
   }
   return Sn + Sn;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Beta_n_Power_Series_Odd_n( double x, int n )            //
//                                                                            //
//  Description:                                                              //
//     The beta exponential integral, beta_n(x), is the integral of           //
//                             t^n exp(-xt) dt                                //
//     where the integral extends from -1 to 1.  The parameter n is called    //
//     the order of the integral.                                             //
//     The i-th deriviative of beta_n(x) is (-1)^i beta_(n+i)(x).  The        //
//     power series expansion of beta_n(x) about x = 0 is:                    //
//     For n odd,                                                             //
//         beta_n(x) = -2 x/(n+2) + 2 x^3/[(3)!(n+4)] + ... +                 //
//                                       2 x^(2i+1)/[(2i+1)!(n+2i+2)] + ...   //
//     Convergence can be slow for large |x|.                                 //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of beta exponential integral beta_n().         //
//     int     n  The order of the beta exponential integral n must be odd.   //
//                                                                            //
//  Return Value:                                                             //
//     The value of beta exponential integral of order n evaluated at x.      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static long double Beta_n_Power_Series_Odd_n( double x, int n )
{
   long double xx = (long double) x;
   long double x2 = xx * xx;
   long double xn = xx;
   long double np1 = (long double) (n + 1);
   long double np2 = (long double) (n + 2);
   long double Sn = xn / np2;
   long double Sm1 = 0.0L;
   long double factorial = 1.0L;
   long double y = 1.0L;
  
   if (x == 0.0) Sn = 0.0L; 
   while ( fabsl(Sn - Sm1) > LDBL_EPSILON * fabsl(Sm1) ) {
      Sm1 = Sn;
      y += 1.0L;
      factorial *= y;
      y += 1.0L;
      factorial *= y;
      xn *= x2;
      Sn += ( xn / (factorial * (np1 + y)) );
   }
   return -(Sn + Sn);
}
