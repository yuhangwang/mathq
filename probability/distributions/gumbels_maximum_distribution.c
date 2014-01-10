////////////////////////////////////////////////////////////////////////////////
// File: gumbels_maximum_distribution.c                                       //
// Routine(s):                                                                //
//    Gumbels_Maximum_Distribution                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Gumbels_Maximum_Distribution( double x )                            //
//                                                                            //
//  Description:                                                              //
//     Gumbel's maximum distribution is the integral from -inf to x of the    //
//     density                                                                //
//                        f(x) = exp(-x) exp(-exp(-x))                        //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = exp(-exp(-x))                      //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Gumbel's maximum distribution Pr[X < x].        //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Gumbels_Maximum_Distribution(x);                                   //
////////////////////////////////////////////////////////////////////////////////

double Gumbels_Maximum_Distribution(double x)
{
   double temp = exp(-x);

   return exp(-temp);
}
