////////////////////////////////////////////////////////////////////////////////
// File: negative_binomial_cumulative_distribution.c                          //
// Routine(s):                                                                //
//    Negative_Binomial_Cumulative_Distribution                               //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Beta_Distribution( double x, double a, double b );

////////////////////////////////////////////////////////////////////////////////
// double Negative_Binomial_Cumulative_Distribution( int n, int k, double p ) //
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
//     The cumulative distribution of X is Pr[X <= k],                        //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 0,...,k <= n and Pr[X = j] is the point      //
//     distribution:                                                          //
//                  Pr[X = j] = C(n+j-1,j) * p^n * (1-p)^j.                   //
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
//     pr = Negative_Binomial_Cumulative_Distribution(n, k, p);               //
////////////////////////////////////////////////////////////////////////////////

double Negative_Binomial_Cumulative_Distribution( int n, int k, double p )
{
   if ( k < 0) return 0.0;
   if ( p == 0.0 ) return 0.0;
   if ( p == 1.0 ) return 1.0;

   return Beta_Distribution( p, (double)n, (double)(k+1));
}
