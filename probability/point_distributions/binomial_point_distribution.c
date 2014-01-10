////////////////////////////////////////////////////////////////////////////////
// File: binomial_point_distribution.c                                        //
// Routine(s):                                                                //
//    Binomial_Point_Distribution                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp(), log()

//                         Externally Defined Routines                        //
double Ln_Gamma_Function(double x);

////////////////////////////////////////////////////////////////////////////////
// double Binomial_Point_Distribution( int n, int k, double p )               //
//                                                                            //
//  Description:                                                              //
//     A random variable Y which can assume one of two values {0,1} is said   //
//     to have a Bernoulli distribution if there is a constant p, 0 < p < 1,  //
//     such that Pr[Y = 1] = p and Pr[Y = 0] = 1-p.                           //
//     If the random variable X is the sum of n independent identically       //
//     distributed Bernoulli random variables Y(i), i = 1,...,n for which     //
//     Pr[Y[i] = 1] = p, then Pr[X = k] = C(n,k) * p^k * (1-p)^(n-k)          //
//     for 0 <= k <= n, where C(n,k) = n! / (k! (n-k)!) and Pr[X = k] = 0 if  //
//     either k < 0 or k > n. The random variable X is said to have a         //
//     binomial distribution.                                                 //
//     The point distribution of X is Pr[X = k],                              //
//                                                                            //
//  Arguments:                                                                //
//     int    n   The total number of trials, n >= 1.                         //
//     int    k   The number of successes (1's), 0 <= k <= n.                 //
//     double p   The probability of a success (1) on a single trial.         //
//                0 <= p <= 1.                                                //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, pr;                                                          //
//     int n,k;                                                               //
//                                                                            //
//     pr = Binomial_Point_Distribution(n, k, p);                             //
////////////////////////////////////////////////////////////////////////////////

double Binomial_Point_Distribution( int n, int k, double p )
{
   if ( k < 0 || k > n) return 0.0;
   if ( p == 0.0 )
      if (k == 0) return 1.0;
      else return 0.0;
   if ( p == 1.0 )
      if ( k == n ) return 1.0;
      else return 0.0;
      
   return exp(Ln_Gamma_Function((double)(n + 1))
              - Ln_Gamma_Function((double)(k + 1))
              - Ln_Gamma_Function((double)(n - k + 1))
              + (double)k * log(p) + (double)(n - k) * log(1.0 - p));
}
