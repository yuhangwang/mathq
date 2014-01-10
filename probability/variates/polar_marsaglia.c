////////////////////////////////////////////////////////////////////////////////
// File: polar_marsaglia.c                                                    //
// Routine(s):                                                                //
//    Gaussian_Variate_Polar_Marsaglia                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for sqrt(), log(), and M_PI
#include <float.h>             // required for DBL_MIN
#define M_2PI (M_PI+M_PI)

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Variate_Polar_Marsaglia( void )                            //
//                                                                            //
//  Description:                                                              //
//     This function returns a Gaussian variate with mean 0 and standard      //
//     deviation 1 using the Polar Marsaglia version of the Box-Muller method //
//     for generating a standard Gaussian variate.                            //
//     The Polar-Marsaglia method uses two uniform variates on (-1,1), u and  //
//     v, for which w = u*u + v*v <= 1 to generate two standard Gaussian      //
//     variates:  u * sqrt(-2 ln(w)/w ) and v * sqrt(-2 ln(w)/w ).            //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A random quantity which has a normal (Gaussian) distribution with      //
//     mean 0 and variance 1.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Gaussian_Variate_Polar_Marsaglia( );                               //
////////////////////////////////////////////////////////////////////////////////

double Gaussian_Variate_Polar_Marsaglia( void )
{
   static double next_variate;
   static int use_next_variate = 0;
   double u, v, w;
   double variate;
   
   if (use_next_variate) variate = next_variate;
   else {
      do { 
         u = 2.0 * Uniform_0_1_Random_Variate() - 1.0;
         v = 2.0 * Uniform_0_1_Random_Variate() - 1.0;
         w = u * u + v * v;
      } while ( w > 1.0 );
      if ( w == 0.0 ) w = DBL_MIN;
      w = sqrt(-2.0 * log(w) / w);
      variate = u * w;
      next_variate = v * w;
   }
   use_next_variate = !use_next_variate;
   return variate;
}
