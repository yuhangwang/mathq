////////////////////////////////////////////////////////////////////////////////
// File: rising_factorial.c                                                   //
// Routine(s):                                                                //
//    Rising_Factorial                                                        //
//    xRising_Factorial                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The rising factorial, (n)_sub_m, is defined as                         //
//                   (n)_sub_m = n (n + 1) ... (n + m - 1)                    //
//     i.e. the rising factorial is defined as Pochhammer's function          //
//     restricted to integer arguments n > 0 and m >= 0, and where by         //
//     convention, (n)_sub_0 = 1.                                             //
//     Using this convention,                                                 //
//          (n)_sub_m = Gamma(n+m) / Gamma(n) = (n+m-1)! / (n-1)!.            //
//                                                                            //
//     The functions Rising_Factorial(n,m) and xRising_Factorial(n,m) return  //
//     (n)_sub_m provided (n + m - 1)! / (n - 1)! < DBL_MAX.  If not, then    //
//     DBL_MAX is returned.                                                   //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                          // required for expl()
#include <float.h>                         // required for DBL_MAX
#include <limits.h>                        // required for ULONG_MAX

//                        Externally Defined Routines                         //

extern long double xFactorial( int n );
extern int    Factorial_Max_Arg( void );
extern long double xLn_Factorial( int n );

//                        Internally Defined Routines                         //

double Rising_Factorial(int n, int m);
long double xRising_Factorial(int n, int m);

//                       Internally Defined Constant(s)                       //

static const long double ln_DBL_MAX = 7.097827128933839967321e2L;

////////////////////////////////////////////////////////////////////////////////
// double Rising_Factorial( int n, int m )                                    //
//                                                                            //
//  Description:                                                              //
//     This function computes (n)_sub_m = n (n + 1) ... (n + m - 1).          //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Base argument of the rising factorial function. n must be   //
//                positive.                                                   //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n is nonpositive, then 0 is returned.  If (n+m-1)! / (n-1)! >       //
//     DBL_MAX then DBL_MAX is returned.  Otherwise (n+m-1)! / (n-1)! is      //
//     returned.                                                              //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     double x;                                                              //
//                                                                            //
//     x = Rising_Factorial( n, m );                                          //
////////////////////////////////////////////////////////////////////////////////
double Rising_Factorial(int n, int m) {

   return (double) xRising_Factorial(n,m);

}


////////////////////////////////////////////////////////////////////////////////
// long double xRising_Factorial( int n, int m )                              //
//                                                                            //
//  Description:                                                              //
//     This function computes (n)_sub_m = n (n + 1) ... (n + m - 1).          //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Base argument of the rising factorial function. n must be   //
//                positive.                                                   //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n is nonpositive, then 0 is returned.  If (n+m-1)! / (n-1)! >       //
//     DBL_MAX then DBL_MAX is returned.  Otherwise (n+m-1)! / (n-1)! is      //
//     returned.                                                              //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     long double x;                                                         //
//                                                                            //
//     x = xRising_Factorial( n, m );                                         //
////////////////////////////////////////////////////////////////////////////////
long double xRising_Factorial(int n, int m) {

   long double ln_poch;
   long double poch;
   unsigned long c;

              // For a nonpositive base argument (n <= 0) or a //
              // negative length product (m < 0) return 0.0    //

   if ( n <= 0 || m < 0 ) return 0.0L;

              // For a zero length product (m = 0), return 1.0 //

   if ( m == 0 ) return 1.0L;

            // If n + m - 1 <= Factorial_Max_Arg(), then return //
            // the quotient (n + m - 1)! / (n - 1)!.            //

   if ( (n + m - 1) <= Factorial_Max_Arg() )
      return xFactorial(n + m - 1) / xFactorial(n - 1);

             // If ln (n+m-1)! - ln (n-1)! < ln DBL_MAX return  //
             // (n + m - 1)! / (n-1)! otherwise return DBL_MAX. //

   ln_poch = xLn_Factorial(n+m-1) - xLn_Factorial(n-1);
   if ( ln_poch < ln_DBL_MAX ) {
      poch = expl(ln_poch);
      if (poch < (long double) ULONG_MAX) {
         c = (unsigned long) (poch + 0.5L);
         poch = (long double) c;
      }
      return poch;
   }

   return (long double) DBL_MAX;
}
