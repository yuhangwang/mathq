////////////////////////////////////////////////////////////////////////////////
// File: cauchy_random_variate.c                                              //
// Routine(s):                                                                //
//    Cauchy_Random_Variate                                                   //
////////////////////////////////////////////////////////////////////////////////

#include <float.h>             // required for DBL_MIN

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Cauchy_Random_Variate( void )                                       //
//                                                                            //
//  Description:                                                              //
//     This function returns a Cauchy distributed random variate using the    //
//     following rejection method:  Let u and v be independent U(-1,1) random //
//     variables then if u^2 + v^2 <= 1 then X = u/v has the Cauchy           //
//     distribution with probability density                                  //
//                       f(x) = 1 / [pi * (1 + x * x)]                        //
//     and distribution                                                       //
//                       F(x) = 1/2 + (1/pi) arctan(x).                       //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Cauchy distributed random quantity.                                  //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Cauchy_Random_Variate();                                           //
////////////////////////////////////////////////////////////////////////////////

double Cauchy_Random_Variate( void )
{
   double u, v;
  
//   Continue sampling two independent Uniform(-1,1) random variables u,v     //
//   until u^2 + v^2 <= 1.                                                    //

   do { 
      u = 2.0 * Uniform_0_1_Random_Variate() - 1.0;
      v = 2.0 * Uniform_0_1_Random_Variate() - 1.0;
   } while ( (u * u + v * v) > 1.0 );

//   In the unlikely event that v = 0 in order to avoid division by zero      //
//   return DBL_MAX.  Otherwise return u/v.                                   //

   return ( v == 0.0 ) ?  DBL_MAX : u / v;
}
