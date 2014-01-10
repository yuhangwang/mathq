////////////////////////////////////////////////////////////////////////////////
// File: pareto_density.c                                                     //
// Routine(s):                                                                //
//    Pareto_Density                                                          //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow()

////////////////////////////////////////////////////////////////////////////////
// double Pareto_Density( double x, double a )                                //
//                                                                            //
//  Description:                                                              //
//     The density of the Pareto distribution is                              //
//                                           0                if x < 1,       //
//                                 f(x) =  a / x^(a+1)        if x >= 1       //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Pareto distribution Pr[X < x].              //
//     double a   Shape parameter of the Pareto distribution, a > 0.          //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and a.                                         //
//                                                                            //
//  Example:                                                                  //
//     double a, p, x;                                                        //
//                                                                            //
//     p = Pareto_Density(x,a);                                               //
////////////////////////////////////////////////////////////////////////////////

double Pareto_Density(double x, double a)
{
   if (x < 1.0) return 0.0;
   return a / pow(x,a+1.0);
}
