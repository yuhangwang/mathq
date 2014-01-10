////////////////////////////////////////////////////////////////////////////////
// File: gumbels_maximum_density.c                                            //
// Routine(s):                                                                //
//    Gumbels_Maximum_Density                                                 //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Gumbels_Maximum_Density( double x )                                 //
//                                                                            //
//  Description:                                                              //
//     The density of Gumbel's maximum distribution is                        //
//                        f(x) = exp(-x) exp(-exp(-x))                        //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Gumbel's maximum density.                       //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and exp(-1) =~ 0.367.                          //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Gumbels_Maximum_Density(x);                                        //
////////////////////////////////////////////////////////////////////////////////

double Gumbels_Maximum_Density(double x)
{
   double temp = exp(-x);
   return temp * exp(-temp);
}
