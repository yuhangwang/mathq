////////////////////////////////////////////////////////////////////////////////
// File: negative_binomial_point_distribution.c                               //
// Routine(s):                                                                //
//    Negative_Binomial_Point_Distribution                                    //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp(), log()

//                         Externally Defined Routines                        //
double Ln_Gamma_Function( double x );

////////////////////////////////////////////////////////////////////////////////
// double Negative_Binomial_Point_Distribution( int n, int k, double p )      //
//                                                                            //
//  Description:                                                              //
//     A random variable Y which can assume one of two values {0,1} is said   //
//     to have a Bernoulli distribution if there is a constant p, 0 < p < 1,  //
//     such that Pr[Y = 1] = p and Pr[Y = 0] = 1-p.                           //
//     Let Y[i], i=1,... be a sequence of independent identically distributed //
//     Bernoulli random variables and let Z[j] = Sum Y[i] where the sum       //
//     extends over i from 1 to j.  Let X = min {j: Z[j] = n} - n, then the   //
//     probability distribution of X is                                       //
//                  Pr[X = k] = C(n+k-1,k) * p^n * (1-p)^k                    //
//     where n > 0, k >= 0, and C(n,k) = n! / (k! (n-k)!). The random variable//
//     X is said to have a negative binomial distribution.                    //
//     The point distribution of X is Pr[X = k].                              //
//                                                                            //
//     Note that some authors say that X+n has a negative binomial            //
//     distribution.                                                          //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The number of 1's to occur before stopping, n >= 1.         //
//     int    k   The number of 0's which occurred before n 1's occurred.     //
//     double p   The probability of a 1 on a single trial. 0 <= p <= 1.      //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, pr;                                                          //
//     int n,k;                                                               //
//                                                                            //
//     pr = Negative_Binomial_Point_Distribution(n, k, p);                    //
////////////////////////////////////////////////////////////////////////////////

double Negative_Binomial_Point_Distribution( int n, int k, double p )
{
   if ( k < 0 ) return 0.0;
   if ( p == 0.0 ) return 0.0;
   if ( p == 1.0 ) return (k == 0) ? 1.0 : 0.0;
      
   return exp(Ln_Gamma_Function((double)(n + k))
              - Ln_Gamma_Function((double)(k + 1))
              - Ln_Gamma_Function((double)(n))
              + (double)n * log(p) + (double)k * log(1.0 - p));
}
