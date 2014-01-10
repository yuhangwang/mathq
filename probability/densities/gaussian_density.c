////////////////////////////////////////////////////////////////////////////////
// File: gaussian_density.c                                                   //
// Routine(s):                                                                //
//    Gaussian_Density                                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for exp(), M_SQRT1_2 and
                               // M_2_SQRTPI

static double const normalization = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Density( double x )                                        //
//                                                                            //
//  Description:                                                              //
//     This function returns the value of the probability density at X = x    //
//     where X has a standard Normal (Gaussian) distribution.                 //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of dPr[X < x]/dx where X ~ N(0,1).                 //
//                                                                            //
//  Return Values:                                                            //
//     The value of the probability density at X = x where X has a Normal     //
//     (Gaussian) distribution with mean 0 and variance 1.                    //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//     pr = Gaussian_Density(x);                                              //
////////////////////////////////////////////////////////////////////////////////

double Gaussian_Density( double x )
{
   return  normalization * exp( - 0.5 * x * x );
}
