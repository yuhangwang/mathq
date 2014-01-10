////////////////////////////////////////////////////////////////////////////////
// File: weibull_random_variate.c                                             //
// Routine(s):                                                                //
//    Weibull_Random_Variate                                                  //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                      // required for pow()

//                    Required Externally Defined Routines                    //
double Exponential_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Weibull_Random_Variate( double a )                                  //
//                                                                            //
//  Description:                                                              //
//     This function returns a Weibull random variate taking the a-th root of //
//     an exponentially distributed random variate.  I.e. if the random       //
//     variable e has an exponential distribution 1-exp[-e], then X = e^(1/a) //
//     is distributed according to the Weibull distribution with              //
//     probability density                                                    //
//                         f(x) = a x^(a-1) exp(-x^a)                         //
//     and distribution                                                       //
//                           F(x) = 1 - exp(-x^a).                            //
//                                                                            //
//  Arguments:                                                                //
//     double a   Shape parameter for the Weibull distribution.               //
//                                                                            //
//  Return Values:                                                            //
//     A Weibull distributed random quantity.                                 //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     double a;                                                              //
//                                                                            //
//     x = Weibull_Random_Variate(a);                                         //
////////////////////////////////////////////////////////////////////////////////

double Weibull_Random_Variate( double a )
{
   double e = Exponential_Random_Variate();

   return ( e == 0.0 ) ?  0.0 : pow(e, 1.0/a);
}
