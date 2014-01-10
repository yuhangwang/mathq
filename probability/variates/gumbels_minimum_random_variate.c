////////////////////////////////////////////////////////////////////////////////
// File: gumbels_minimum_random_variate.c                                     //
// Routine(s):                                                                //
//    Gumbels_Minimum_Random_Variate                                          //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                      // required for log()

//                    Required Externally Defined Routines                    //
double Exponential_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Gumbels_Minimum_Random_Variate( void )                              //
//                                                                            //
//  Description:                                                              //
//     This function returns a Gumbel's minimum distributed random variate by //
//     taking the ln of an exponentially distributed random variate.  I.e. if //
//     the random variable e has an exponential distribution 1-exp[-e], then  //
//     X = ln(e) is distributed according to the Gumbel's minimum             //
//     distribution with probability density                                  //
//                       f(x) = exp(x) exp( - exp(x) )                        //
//     and distribution                                                       //
//                        F(x) = 1 - exp( -exp(x) ).                          //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Gumbel's minimum distributed random quantity.                        //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Gumbels_Minimum_Random_Variate();                                  //
////////////////////////////////////////////////////////////////////////////////

double Gumbels_Minimum_Random_Variate( void )
{
   double e = Exponential_Random_Variate();

   return log(e);
}
