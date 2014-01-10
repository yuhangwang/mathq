////////////////////////////////////////////////////////////////////////////////
// File: weibull_density.c                                                    //
// Routine(s):                                                                //
//    Weibull_Density                                                         //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow(), exp(), and log()
#include <float.h>                   // required for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// double Weibull_Density( double x, double a )                               //
//                                                                            //
//  Description:                                                              //
//     The density of the Weibull distribution is                             //
//                                           0                if x <= 0,      //
//                              f(x) =  a x^(a-1) exp(-x^a)   if x > 0        //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the Weibull density.                            //
//     double a   Shape parameter of the Weibull distribution, a > 0.         //
//                                                                            //
//  Return Values:                                                            //
//     If a > 1, a real number between 0 and a((a-1)/a)^((a-1)/a)exp(-(a-1)/a)//
//     If a = 1, a real number between 0 and a.                               //
//     If a < 1, a real number between 0 and DBL_MAX.  If f(x) > DBL_MAX, then//
//     DBL_MAX is returned.                                                   //
//                                                                            //
//  Example:                                                                  //
//     double a, p, x;                                                        //
//                                                                            //
//     p = Weibull_Density(x,a);                                              //
////////////////////////////////////////////////////////////////////////////////

double Weibull_Density(double x, double a)
{
   double temp;

   if (x < 0.0) return 0.0;
   if (x == 0.0) {
      if (a > 1.0) return 0.0;
      if (a == 1.0) return 1.0;
      return DBL_MAX;
   }
   temp = log(a) + (a - 1.0) * log(x) - pow(x,a);
   if (temp >= log(DBL_MAX)) return DBL_MAX;
   return exp(temp);
}
