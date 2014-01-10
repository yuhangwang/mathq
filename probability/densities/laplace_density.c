////////////////////////////////////////////////////////////////////////////////
// File: laplace_density.c                                                    //
// Routine(s):                                                                //
//    Laplace_Density                                                         //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                              // required for exp(), fabs()

////////////////////////////////////////////////////////////////////////////////
// double Laplace_Density( double x )                                         //
//                                                                            //
//  Description:                                                              //
//     The density of Laplace's distribution is                               //
//                           f(x) = (1/2) exp(-|x|)                           //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Laplace's density.                              //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1/2.                                       //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Laplace_Density(x);                                                //
////////////////////////////////////////////////////////////////////////////////

double Laplace_Density(double x)
{
   return 0.5 * exp(-fabs(x));
}
