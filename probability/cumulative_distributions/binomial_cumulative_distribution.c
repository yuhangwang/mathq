////////////////////////////////////////////////////////////////////////////////
// File: binomial_cumulative_distribution.c                                   //
// Routine(s):                                                                //
//    Binomial_Cumulative_Distribution                                        //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Beta_Distribution( double x, double a, double b );

////////////////////////////////////////////////////////////////////////////////
// double Binomial_Cumulative_Distribution( int n, int k, double p )          //
//                                                                            //
//  Description:                                                              //
//     A random variable Y which can assume one of two values {0,1} is said   //
//     to have a Bernoulli distribution if there is a constant p, 0 < p < 1,  //
//     such that Pr[Y = 1] = p and Pr[Y = 0] = 1-p.                           //
//     If the random variable X is the sum of n independent identically       //
//     distributed Bernoulli random variables Y(i), i = 1,...,n for which     //
//     Pr[Y(i) = 1] = p, then X is said to have a binomial distribution.      //
//     The cumulative distribution function of X is Pr[X <= k],               //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 0,...,k and Pr[X = j] is the point           //
//     distribution:                                                          //
//                  Pr[X = j] = C(n,j) * p^j * (1-p)^(n-j)                    //
//     for 0 <= j <= n, where C(n,j) = n! / (j! (n-j)!) and Pr[X = j] = 0 if  //
//     either j < 0 or j > n.                                                 //
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
//     pr = Binomial_Cumulative_Distribution(n, k, p);                        //
////////////////////////////////////////////////////////////////////////////////

double Binomial_Cumulative_Distribution( int n, int k, double p )
{
   if ( k < 0) return 0.0;
   if ( k >= n ) return 1.0;
   if ( p == 0.0 ) return 1.0;
   if ( p == 1.0 ) return (k < n) ? 0.0 : 1.0;
   return Beta_Distribution( 1.0 - p, (double)(n-k), (double)(k+1));
}
