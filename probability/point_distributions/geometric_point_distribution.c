////////////////////////////////////////////////////////////////////////////////
// File: geometric_point_distribution.c                                       //
// Routine(s):                                                                //
//    Geometric_Point_Distribution                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for pow()

////////////////////////////////////////////////////////////////////////////////
// double Geometric_Point_Distribution( int k, double p )                     //
//                                                                            //
//  Description:                                                              //
//     A random variable Y which can assume one of two values {0,1} is said   //
//     to have a Bernoulli distribution if there is a constant p, 0 < p < 1,  //
//     such that Pr[Y = 1] = p and Pr[Y = 0] = 1-p.                           //
//     Let Y[i], i=1,... be a sequence of independent identically distributed //
//     Bernoulli random variables and let Z[j] = Sum Y[i] where the sum       //
//     extends over i from 1 to j.  Let X = min {j: Z[j] = 1} - 1, then the   //
//     probability distribution of X is                                       //
//                         Pr[X = k] =  p * (1-p)^k                           //
//     where k >= 0. The random variable X is said to have a geometric        //
//     distribution.                                                          //
//     The point distribution of X is Pr[X = k].                              //
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
//     int k;                                                                 //
//                                                                            //
//     pr = Geometric_Point_Distribution(k, p);                               //
////////////////////////////////////////////////////////////////////////////////

double Geometric_Point_Distribution( int k, double p )
{
   if ( k < 0 ) return 0.0;
   if ( p == 0.0 ) return 0.0;
   if ( p == 1.0 ) return (k == 0) ? 1.0 : 0.0;
      
   return p * pow(1.0 - p, (double)k);
}
