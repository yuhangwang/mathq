////////////////////////////////////////////////////////////////////////////////
// File: multinomial_coefficient.c                                            //
// Routine(s):                                                                //
//    Multinomial_Coefficient                                                 //
//    xMultinomial_Coefficient                                                //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The multinomial coefficient is number of ways of partitioning n >= 0   //
//     objects into m sets of size x[0], x[1],...,x[m-2], n-x[0]-...-x[m-2]   //
//     respectively where 0 <= x[i] <= n and x[0] + ... x[m-2] <= n,          //
//     C(n;x) = C(n;x[0],x[1],...,x[m-2])                                     //
//                                    n!                                      //
//            = ------------------------------------------------              //
//               x[0]! x[1]! ... x[m-2]! ( n-x[0]-...-x[m-2] )!               //
//                                                                            //
//     The functions Multinomial_Coefficient(n,x,m) and                       //
//     xMultinomial_Coefficient(n,x,m) return C(n;x).                         //
//     If C(n;x) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     if x[i] < 0 for some i, 0 <= i <= m-2 or if n < x[0] + ... + x[m-2],   //
//     then 0 is returned.                                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                          // required for expl()
#include <float.h>                         // required for DBL_MAX
#include <limits.h>                        // required for ULONG_MAX

//                        Externally Defined Routines                         //

extern long double xFactorial( int n );
extern int    Factorial_Max_Arg( void );
extern long double xLn_Factorial( int n );

//                        Internally Defined Routines                         //

double Multinomial_Coefficient(int n, int x[], int m);
long double xMultinomial_Coefficient(int n, int x[], int m);

//                       Internally Defined Constant(s)                       //

static const long double ln_DBL_MAX = 7.097827128933839967321e2L;

////////////////////////////////////////////////////////////////////////////////
// double Multinomial_Coefficient( int n, int x[], int m )                    //
//                                                                            //
//  Description:                                                              //
//     If C(n;x) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     if x[i] < 0 for some i, 0 <= i <= m-2 or if n < x[0] + ... + x[m-2],   //
//     then 0 is returned.  Otherwise C(n;x) is returned.                     //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    x[] The ith component,x[i], is the number in the ith set,       //
//                0 <= i < m-2. The number in the (m-1)st set is              //
//                n - x[0]-...-x[m-2].                                        //
//     int    m   The number of factorial terms in the quotient,              //
//                m = dim x + 1.                                              //
//                                                                            //
//  Return Values:                                                            //
//     If n or x[i] for some i are negative, then 0 is returned.  If          //
//     n < x[0] + ... + x[m-2] then 0 is returned.  If C(n;x) >= DBL_MAX,     //
//     then DBL_MAX is returned otherwise, C(n;x) is returned.                //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     int n;                                                                 //
//     int int x[M]                                                           //
//     double c;                                                              //
//                                                                            //
//     c = Multinomial_Coefficient( n, x, M );                                //
////////////////////////////////////////////////////////////////////////////////
double Multinomial_Coefficient(int n, int x[], int m) {
   return (double)xMultinomial_Coefficient(n,x,m);
}


////////////////////////////////////////////////////////////////////////////////
// long double xMultinomial_Coefficient( int n, int x[], int m )              //
//                                                                            //
//  Description:                                                              //
//     If C(n;x) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     if x[i] < 0 for some i, 0 <= i <= m-2 or if n < x[0] + ... + x[m-2],   //
//     then 0 is returned.  Otherwise C(n;x) is returned.                     //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    x[] The ith component,x[i], is the number in the ith set,       //
//                0 <= i < m-2. The number in the (m-1)st set is              //
//                n - x[0]-...-x[m-2].                                        //
//     int    m   The number of factorial terms in the quotient,              //
//                m = dim x + 1.                                              //
//                                                                            //
//  Return Values:                                                            //
//     If n or x[i] for some i are negative, then 0 is returned.  If          //
//     n < x[0] + ... + x[m-2] then 0 is returned.  If C(n;x) >= DBL_MAX,     //
//     then DBL_MAX is returned otherwise, C(n;x) is returned.                //
//                                                                            //
//  Example:                                                                  //
//     #define M                                                              //
//     int n;                                                                 //
//     int int x[M]                                                           //
//     long double c;                                                         //
//                                                                            //
//     c = xMultinomial_Coefficient( n, x, M );                               //
////////////////////////////////////////////////////////////////////////////////
long double xMultinomial_Coefficient(int n, int x[], int m) {

   long double ln_combination;
   long double combination;
   unsigned long c;
   int    nx = 0;
   int    i;

                          // For n < 0 return 0.0. //

   if (n < 0)  return 0.0L;

      // If n <= Factorial_Max_Arg(), then if all x[i] >= 0 and if     //
      // n >= x[0] + x[1] + ... + x[m-2] then return the quotient      //
      // n! / x[0]! / x[1]! / ... x[m-2]! / (n - x[0] - ... - x[m-2])! //
      // otherwise return 0.                                           //

   if ( n <= Factorial_Max_Arg() ) {
      combination = xFactorial(n);
      for (i = 0; i < (m-1); i++) {
         if (x[i] < 0) return 0.0L;
         nx += x[i];
         combination /= xFactorial(x[i]);
      }
      if (nx > n) return 0.0L;
      return combination / xFactorial(n-nx);
   }

      // If n > Factorial_Max_Arg(), then if all x[i] >= 0 and if       //
      // n >= x[0] + x[1] + ... + x[m-2] then if the difference         //
      // ln(n!) - ln(x[0]!) - ... - ln(x[m-2]!)                         //
      //                 - ln((n - x[0] - ... - x[m-2])! > ln(DBL_MAX)  //
      // then return DBL_MAX, otherwise return C(n;x).                  //
      // Further since C(n;x) is an integer, if possible correct        //
      // roundoff errors.                                               //


   ln_combination = xLn_Factorial( n );
   for (i = 0; i < (m-1); i++) {
      if (x[i] < 0) return 0.0L;
      nx += x[i];
      ln_combination -= xLn_Factorial( x[i] );
   }
   if (nx > n) return 0.0L;
   ln_combination -= xLn_Factorial( n - nx );
   if ( ln_combination < ln_DBL_MAX ) {
      combination = expl(ln_combination);
      if (combination < (double) ULONG_MAX) {
         c = (unsigned long)(combination + 0.5L);
         combination = (long double) c;
      }
      return combination;
   }

   return (long double) DBL_MAX;
}
