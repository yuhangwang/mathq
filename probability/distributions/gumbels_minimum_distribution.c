////////////////////////////////////////////////////////////////////////////////
// File: gumbels_minimum_distribution.c                                       //
// Routine(s):                                                                //
//    Gumbels_Minimum_Distribution                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Gumbels_Minimum_Distribution( double x )                            //
//                                                                            //
//  Description:                                                              //
//     Gumbel's minimum distribution is the integral from -inf to x of the    //
//     density                                                                //
//                        f(x) = exp(x) exp(-exp(x))                          //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = 1 - exp(-exp(x))                   //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Gumbel's minimum distribution Pr[X < x].        //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Gumbels_Minimum_Distribution(x);                                   //
////////////////////////////////////////////////////////////////////////////////

double Gumbels_Minimum_Distribution(double x)
{
   double temp = exp(x);

   return 1.0 - exp(-temp);
}
