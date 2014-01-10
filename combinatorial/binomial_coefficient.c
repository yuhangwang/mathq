////////////////////////////////////////////////////////////////////////////////
// File: binomial_coefficient.c                                               //
// Routine(s):                                                                //
//    Binomial_Coefficient                                                    //
//    xBinomial_Coefficient                                                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The binomial coefficients are combination of n >= 0 things taken m,    //
//     0 <= m <= n, at a time, i.e. C(n,m) =  n! / [m! (n-m)!].               //
//                                                                            //
//     The functions Binomial_Coefficient(n,m) and xBinomial_Coefficient(n,m) //
//     return C(n,m) if n >= 0 and 0 <= m <= n and C(n,m) < DBL_MAX.          //
//     If C(n,m) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     m < 0 or m > n, then 0 is returned.                                    //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                          // required for expl()
#include <float.h>                         // required for DBL_MAX
#include <limits.h>                        // required for ULONG_MAX

//                        Externally Defined Routines                         //

extern long double xFactorial( int n );
extern int    Factorial_Max_Arg( void );
extern long double xLn_Factorial( int n );

//                        Internally Defined Routines                         //

double Binomial_Coefficient(int n, int m);
long double xBinomial_Coefficient(int n, int m);

//                        Internally Defined Constants                        //

static const long double ln_DBL_MAX = 7.097827128933839967321e2L;

////////////////////////////////////////////////////////////////////////////////
// double Binomial_Coefficient( int n, int m )                                //
//                                                                            //
//  Description:                                                              //
//     This function returns C(n,m) for n >= 0, 0 <= m <= n and               //
//     C(n,m) < DBL_MAX.                                                      //
//                                                                            //
//     If C(n,m) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     m < 0 or m > n, then 0 is returned.                                    //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n or m or both are negative,  then 0 is returned.  If m > n, then   //
//     0 is returned.  C(n,m) >= DBL_MAX, then DBL_MAX is returned.           //
//     Otherwise C(n,m) is returned.                                          //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     double x;                                                              //
//                                                                            //
//     x = Binomial_Coefficient( n, m );                                      //
////////////////////////////////////////////////////////////////////////////////
double Binomial_Coefficient(int n, int m) {
   return (double) xBinomial_Coefficient(n,m);
}


////////////////////////////////////////////////////////////////////////////////
// long double xBinomial_Coefficient( int n, int m )                          //
//                                                                            //
//  Description:                                                              //
//     This function returns C(n,m) for n >= 0, 0 <= m <= n and               //
//     C(n,m) < DBL_MAX.                                                      //
//                                                                            //
//     If C(n,m) >= DBL_MAX, then DBL_MAX is returned.  If n < 0 or           //
//     m < 0 or m > n, then 0 is returned.                                    //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of objects.                                      //
//     int    m   Number of terms in the product.  m must be nonnegative.     //
//                                                                            //
//  Return Values:                                                            //
//     If n or m or both are negative,  then 0 is returned.  If m > n, then   //
//     0 is returned.  C(n,m) >= DBL_MAX, then DBL_MAX is returned.           //
//     Otherwise C(n,m) is returned.                                          //
//                                                                            //
//  Example:                                                                  //
//     int m, n;                                                              //
//     long double x;                                                         //
//                                                                            //
//     x = xBinomial_Coefficient( n, m );                                     //
////////////////////////////////////////////////////////////////////////////////
long double xBinomial_Coefficient(int n, int m) {

   long double ln_combination;
   long double combination;
   unsigned long c;

                 // For n < 0 or m < 0 or m > n return 0.0. //

   if ( (n < 0) || (m < 0) || (m > n) ) return 0.0L;

           // If n <= Factorial_Max_Arg(), simply return //
           // the quotient n! / m! (n - m)!.             //

   if ( n <= Factorial_Max_Arg() )
      return xFactorial(n) / ( xFactorial(m) * xFactorial(n-m) );

        // If ln(n!) - ln(m!) - ln((n-m)!) < ln(DBL_MAX) then return //
        // exp(ln(n!) - ln(m!) - ln((n-m)!) otherwise return DBL_MAX //
        // Further since C(n,m) is an integer, if possible correct   //
        // roundoff errors.                                          //

   ln_combination = xLn_Factorial(n) - xLn_Factorial(m) - xLn_Factorial(n-m);
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
