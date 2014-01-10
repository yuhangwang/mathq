////////////////////////////////////////////////////////////////////////////////
// File: logistic_density.c                                                   //
// Routine(s):                                                                //
//    Logistic_Density                                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                 // required for exp()

////////////////////////////////////////////////////////////////////////////////
// double Logistic_Density( double x )                                        //
//                                                                            //
//  Description:                                                              //
//     The density of the Logistic distribution is                            //
//                       f(x) =  exp(-x) / (1 + exp(-x))^2.                   //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Logistic density.                           //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1/4.                                       //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//                                                                            //
//     p = Logistic_Density(x);                                               //
////////////////////////////////////////////////////////////////////////////////

double Logistic_Density(double x)
{
   double temp = exp(-x);
   double tempp1 = 1.0 + temp;

   return temp / (tempp1 * tempp1);
}
