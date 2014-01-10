////////////////////////////////////////////////////////////////////////////////
// File: log_series_point_distribution.c                                      //
// Routine(s):                                                                //
//    Log_Series_Point_Distribution                                           //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                          // required for log() and pow()

////////////////////////////////////////////////////////////////////////////////
// double Log_Series_Point_Distribution( int k, double p )                    //
//                                                                            //
//  Description:                                                              //
//     A random variable X is said to have a logarthmic series distribution if//
//                       Pr[X = k] = - p^k / [k ln(1-p)],  k = 1,...          //
//     where 0 < p < 1.0.                                                     //
//     The point distribution of X is Pr[X = k].                              //
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
//     pr = Logarithmic_Series_Point_Distribution(k, p);                      //
////////////////////////////////////////////////////////////////////////////////

double Log_Series_Point_Distribution( int k, double p )
{
   if ( k < 1 ) return 0.0;
      
   return - pow(p, (double) k) / ( (double)k * log(1.0 - p) );
}
