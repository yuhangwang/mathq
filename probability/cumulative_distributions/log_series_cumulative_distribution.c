////////////////////////////////////////////////////////////////////////////////
// File: log_series_cumulative_distribution.c                                 //
// Routine(s):                                                                //
//    Log_Series_Cumulative_Distribution                                      //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                          // required for log()

////////////////////////////////////////////////////////////////////////////////
// double Log_Series_Cumulative_Distribution( int k, double p )               //
//                                                                            //
//  Description:                                                              //
//     A random variable X is said to have a logarthmic series distribution if//
//                       Pr[X = k] = - p^k / [k ln(1-p)],  k = 1,...          //
//     where 0 < p < 1.0.                                                     //
//     The cumulative distribution function of X is Pr[X <= k],               //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 1,...,k <= n and Pr[X = j] is the point      //
//     distribution described above.                                          //
//                                                                            //
//  Arguments:                                                                //
//     int    k   The argument of the logarithmic distribution, k >= 1.       //
//     double p   The shape parameter of a logarithmic distribution,          //
//                0 < p < 1.                                                  //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, pr;                                                          //
//     int k;                                                                 //
//                                                                            //
//     pr = Log_Series_Cumulative_Distribution(k, p);                         //
////////////////////////////////////////////////////////////////////////////////

double Log_Series_Cumulative_Distribution( int k, double p )
{
   double summand = p;
   double sum = summand;
   int i;

   if ( k <= 0) return 0.0;
   for (i = 2; i <= k; i++) {
      summand *= p;
      sum += summand / (double)i;
   }
   return -sum / log(1.0 - p);
}
