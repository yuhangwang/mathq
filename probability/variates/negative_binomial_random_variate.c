////////////////////////////////////////////////////////////////////////////////
// File: negative_binomial_random_variate.c                                   //
// Routine(s):                                                                //
//    Negative_Binomial_Random_Variate                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for log()

//                    Required Externally Defined Routines                    //

extern double Gamma_Random_Variate( double a );
extern int Poisson_Random_Variate( double mu );

////////////////////////////////////////////////////////////////////////////////
// int Negative_Binomial_Random_Variate( int n, double p )                    //
//                                                                            //
//  Description:                                                              //
//     This function returns a Negative Binomial(n,p) distributed random      //
//     variate the algorithm: If Y~Gamma(n) and X~Poisson((1-p)Y/p), then X   //
//     has a Negative Binomial(n,p) distribution.  The number returned        //
//     corresponds to the number of failures which occur before n successes   //
//     occur where p is the probability of a success.                         //
//                                                                            //
//  Arguments:                                                                //
//     int    n                                                               //
//            The total number of successes.                                  //
//     double p                                                               //
//            The probability of a success, 0 < p < 1.                        //
//                                                                            //
//  Return Values:                                                            //
//     A random number distributed as a Negative Binomial(n,p) distribution.  //
//                                                                            //
//  Example:                                                                  //
//     int x, n;                                                              //
//     double p;                                                              //
//                                                                            //
//                  (* Set the probability p, 0 < p < 1 *)                    //
//                (* Set the total number of successes n *)                   //
//                                                                            //
//     x = Negative_Binomial_Random_Variate( n, p );                          //
////////////////////////////////////////////////////////////////////////////////
        
int Negative_Binomial_Random_Variate( int n, double p )
{
   double y = Gamma_Random_Variate( (double)n );

   return Poisson_Random_Variate( (1.0 - p) * y / p ); 
}
