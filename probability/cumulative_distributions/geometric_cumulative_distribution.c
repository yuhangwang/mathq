////////////////////////////////////////////////////////////////////////////////
// File: geometric_cumulative_distribution.c                                  //
// Routine(s):                                                                //
//    Geometric_Cumulative_Distribution                                       //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for pow()

////////////////////////////////////////////////////////////////////////////////
// double Geometric_Cumulative_Distribution( int k, double p )                //
//                                                                            //
//  Description:                                                              //
//     A random variable Y which can assume one of two values {0,1} is said   //
//     to have a Bernoulli distribution if there is a constant p, 0 < p < 1,  //
//     such that Pr[Y = 1] = p and Pr[Y = 0] = 1-p.                           //
//     Let Y(i), i=1,... be a sequence of independent identically distributed //
//     Bernoulli random variables, define Y(0) = 0, and let X be the smallest //
//     integer such that Y(X) = 0 and Y(X + 1) = 1.  Then X is said to have   //
//     a geometric distribution.                                              //
//     The cumulative distribution function of X is Pr[X <= k],               //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 0,...,k <= n and Pr[X = j] is the point      //
//     distribution:                                                          //
//                         Pr[X = j] =  p * (1-p)^j.                          //
//                                                                            //
//  Arguments:                                                                //
//     int    k   The number of 0's which occurred before a 1 occurred.       //
//     double p   The probability of a 1 on a single trial. 0 <= p <= 1.      //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, pr;                                                          //
//     int n,k;                                                               //
//                                                                            //
//     pr = Geometric_Cumulative_Distribution(k, p);                          //
////////////////////////////////////////////////////////////////////////////////

double Geometric_Cumulative_Distribution( int k, double p )
{
   if ( k < 0) return 0.0;
   if ( p == 0.0 ) return 0.0;
   if ( p == 1.0 ) return 1.0;

   return 1.0 - pow((1.0 - p),(double)(k+1));
}
