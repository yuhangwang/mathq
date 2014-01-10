////////////////////////////////////////////////////////////////////////////////
// File: gumbels_minimum_density.c                                            //
// Routine(s):                                                                //
//    Gumbels_Minimum_Density                                                 //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Gumbels_Minimum_Density( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The density of Gumbel's minimum distribution is                        //
//                        f(x) = exp(x) exp(-exp(x))                          //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Gumbel's minimum density.                       //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and exp(-1) =~ 0.367.                          //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Gumbels_Minimum_Density(x);                                        //
////////////////////////////////////////////////////////////////////////////////

double Gumbels_Minimum_Density(double x)
{
   double temp = exp(x);
   return temp * exp(-temp);
}
