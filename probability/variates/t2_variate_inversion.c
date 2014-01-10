////////////////////////////////////////////////////////////////////////////////
// File: t2_variate_inversion.c                                               //
// Routine(s):                                                                //
//    t2_Variate_Inversion                                                    //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for log()
#include <float.h>                            // required for DBL_MAX

//                    Required Externally Defined Routines                    //

extern double Uniform_0_1_Random_Variate( void );

////////////////////////////////////////////////////////////////////////////////
// double t2_Variate_Inversion( void )                                        //
//                                                                            //
//  Description:                                                              //
//     This function returns a random variate distributed as a t-distribution //
//     with two degrees of freedom using the inversion method.  This routine  //
//     would normally be used in the rejection method for generating random   //
//     variates with a more complicated distribution.  The density of the     //
//     t-distribution with two degrees of freedom is:                         //
//                 t2(x) = [ (1 + x^2 / 2) ^ (-3/2) ]/sqrt(8)                 //
//     with corresponding distribution function:                              //
//                    T2(x) = ( 1 + x / sqrt(2 + x^2)) / 2.                   //
//     If a random variable X has an t2 distribution, then T2(X) has a        //
//     uniform (0,1) distribution.  So if U is a uniform(0,1), setting        //
//                      U = ( 1 + X / sqrt(2 + X^2)) / 2                      //
//     and solving for X,                                                     //
//                X = sqrt(2) * [ (U - 1/2) / sqrt(U (1 - U)) ]               //
//     has a t-distribution with 2 degrees of freedom.                        //
//                                                                            //
//  Arguments:                                                                //
//     None                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     A random number with a t-distribution with 2 degrees of freedom.       //
//     The value returned is in the interval [-DBL_MAX, DBL_MAX].             //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = t2_Variate_Inversion();                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                         // required for M_SQRT2 and sqrt()
#include <float.h>                        // required for DBL_MAX

double t2_Variate_Inversion( void )
{
   double u = Uniform_0_1_Random_Variate();
   if (u == 0.0) return -DBL_MAX;
   if (u == 1.0) return DBL_MAX;
   return M_SQRT2 * (u - 0.5) / sqrt( u * (1.0 - u) );
}
