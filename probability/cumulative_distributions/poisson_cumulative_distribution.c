////////////////////////////////////////////////////////////////////////////////
// File: poisson_cumulative_distribution.c                                    //
// Routine(s):                                                                //
//    Poisson_Cumulative_Distribution                                         //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Gamma_Distribution( double x, double a );

////////////////////////////////////////////////////////////////////////////////
// double Poisson_Cumulative_Distribution( int k, double mu )                 //
//                                                                            //
//  Description:                                                              //
//     Classically the Poisson distribution arose as the limit of a Binomial  //
//     distribution, Binomial(n,p), in which n -> infinity as the product     //
//     np -> mu, where mu is a constant.                                      //
//     If X has a Poisson distribution, the point distribution of X is:       //
//                       Pr[X = k] = mu^k e^-mu / k!.                         //
//     The cumulative distribution of X is Pr[X <= k],                        //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//                                   = 1 - Gamma_Distribution(mu,k+1)         //
//     where the sum is over j = 0,...,k <= n, Pr[X = j] is the point         //
//     distribution given above, and Gamma_Distribution(mu,k+1) is the area   //
//     of the gamma distribution with shape parameter k+1 from 0 to mu.       //
//                                                                            //
//  Arguments:                                                                //
//     int    k   The maximum number of events.                               //
//     double mu  The expected number of events.                              //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double mu, pr;                                                         //
//     int n,k;                                                               //
//                                                                            //
//     pr = Poisson_Cumulative_Distribution(k, mu);                           //
////////////////////////////////////////////////////////////////////////////////

double Poisson_Cumulative_Distribution( int k, double mu )
{
   if ( k < 0) return 0.0;

   return 1.0 - Gamma_Distribution( mu, (double)(k+1) );
}
