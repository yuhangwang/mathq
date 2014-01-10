////////////////////////////////////////////////////////////////////////////////
// File: pareto_distribution.c                                                //
// Routine(s):                                                                //
//    Pareto_Distribution                                                     //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow()

////////////////////////////////////////////////////////////////////////////////
// double Pareto_Distribution( double x, double a )                           //
//                                                                            //
//  Description:                                                              //
//     The Pareto distribution is the integral from -inf to x of the density  //
//                                           0                if x < 1,       //
//                                 f(x) =  a / x^(a+1)        if x >= 1       //
//     i.e.                                                                   //
//                           Pr[X < x] = F(x) = 0             if x < 1.       //
//                           Pr[X < x] = F(x) = 1 - 1 / x^a   if x > 1.       //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Pareto distribution Pr[X < x].              //
//     double a   Shape parameter of the Pareto distribution, a > 0.          //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double a, p, x;                                                        //
//                                                                            //
//     p = Pareto_Distribution(x,a);                                          //
////////////////////////////////////////////////////////////////////////////////

double Pareto_Distribution(double x, double a)
{
   if (x < 1.0) return 0.0;
   return 1.0 - 1.0 / pow(x,a);
}
