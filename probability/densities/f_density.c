////////////////////////////////////////////////////////////////////////////////
// File: f_density.c                                                          //
// Routine(s):                                                                //
//    F_Density                                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for log() and exp()
#include <float.h>                   // required for DBL_MAX     

//                         Externally Defined Routines                        //
extern double Ln_Beta_Function(double a, double b);

////////////////////////////////////////////////////////////////////////////////
// double F_Density( double x, int v1, int v2 )                               //
//                                                                            //
//  Description:                                                              //
//     The density of the F-distribution is                                   //
//                               0                           if x < 0,        //
//       [v1^(v1/2) v2^(v2/2) / B(v1/2,v2/2)]                                 //
//                 * x^(v1/2-1) * (v2+v1 x)^(-(v1+v2)/2))    if x >= 0,       //
//     where v1 >= 1, v2 >= 1, and B(,) is the (complete) beta function.      //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of the density.                                    //
//     int   v1   The number of degrees of freedom of the numerator of the    //
//                F-test, i.e. (v1/2 - 1) is the exponent of x in the         //
//                density above.                                              //
//     int   v2   The number of degrees of freedom of the denominator of the  //
//                F-test, i.e. (-(v1+v2)/2) is the exponent of the term       //
//                (v2 + v1 x) in the density above.                           //
//                                                                            //
//  Return Values:                                                            //
//     If x <= 0, then 0 is returned, if 0 < x < inf then                     //
//      v1^(v1/2)*(v2^(v2/2)/B(v1/2,v2/2)*x^(v1/2-1)*(v2+v1*x)^(-(v1+v2)/2)   //
//     is returned.                                                           //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    v1, v2;                                                         //
//                                                                            //
//     p = F_Density(x, v1, v2);                                              //
////////////////////////////////////////////////////////////////////////////////

double F_Density( double x, int v1, int v2 )
{
   double ln_density;
   double v12 = (double)v1 / 2.0;
   double v22 = (double)v2 / 2.0;

   if ( x <= 0.0 ) return 0.0;
   ln_density = v12*log((double)v1) + v22 * log((double)v2)
                + (v12 - 1.0) * log(x) - (v12 + v22) * log((double)v2+v1*x)
                - Ln_Beta_Function(v12, v22);
   return exp(ln_density);
}
