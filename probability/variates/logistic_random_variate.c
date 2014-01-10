////////////////////////////////////////////////////////////////////////////////
// File: logistic_random_variate.c                                            //
// Routine(s):                                                                //
//    Logistic_Random_Variate                                                 //
////////////////////////////////////////////////////////////////////////////////

#include <math.h> 
#include <float.h>             // required for DBL_MIN

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Logistic_Random_Variate( void )                                     //
//                                                                            //
//  Description:                                                              //
//     This function returns a Logistic random variate: Let u be a U(0,1)     //
//     random variable then x = log u/(1-u) has a logistic distribution       //
//     probability density                                                    //
//                     f(x) = exp(-x) / (1 + exp(-x))^2                       //
//     and distribution                                                       //
//                     F(x) = 1 / (1 + exp(-x)).                              //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Logistic distributed random quantity.                                //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Logistic_Random_Variate();                                         //
////////////////////////////////////////////////////////////////////////////////

double Logistic_Random_Variate( void )
{
   double u = Uniform_0_1_Random_Variate();

   return log (u / (1.0 - u));
}
