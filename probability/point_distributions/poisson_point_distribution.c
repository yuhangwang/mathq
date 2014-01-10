////////////////////////////////////////////////////////////////////////////////
// File: poisson_point_distribution.c                                         //
// Routine(s):                                                                //
//    Poisson_Point_Distribution                                              //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp(), log()

//                         Externally Defined Routines                        //
double Ln_Gamma_Function( double x );

////////////////////////////////////////////////////////////////////////////////
// double Poisson_Point_Distribution( int k, double mu )                      //
//                                                                            //
//  Description:                                                              //
//     Classically the Poisson distribution arose as the limit of a Binomial  //
//     distribution, Binomial(n,p), in which n -> infinity as the product     //
//     np -> mu, where mu is a constant.                                      //
//     If X has a Poisson distribution, the point distribution of X is:       //
//                       Pr[X = k] = mu^k e^-mu / k!.                         //
//                                                                            //
//  Arguments:                                                                //
//     int    k   The number of events.                                       //
//     double mu  The expected number of events.                              //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double mu, pr;                                                         //
//     int k;                                                                 //
//                                                                            //
//     pr = Poisson_Point_Distribution(k, mu);                                //
////////////////////////////////////////////////////////////////////////////////

double Poisson_Point_Distribution( int k, double mu )
{
   if ( k < 0 ) return 0.0;
      
   return exp( (double)k * log(mu) - mu - Ln_Gamma_Function((double)(k + 1)) );
}
