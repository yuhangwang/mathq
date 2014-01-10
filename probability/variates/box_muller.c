////////////////////////////////////////////////////////////////////////////////
// File: box_muller.c                                                         //
// Routine(s):                                                                //
//    Gaussian_Variate_Box_Muller                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for sqrt(), log(), cos(), and sin()
                               // and M_PI
#include <float.h>             // required for DBL_MIN
#define M_2PI (M_PI+M_PI)

//                    Required Externally Defined Routines                    //
        
extern double Uniform_0_1_Random_Variate( void );

////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Variate_Box_Muller( void )                                 //
//                                                                            //
//  Description:                                                              //
//     This function returns a Gaussian variate with mean 0 and standard      //
//     deviation 1 using the Box-Muller method for generating a standard      //
//     Gaussian variate.                                                      //
//     The Box-Muller method uses two independent uniform variates on (0,1),  //
//     u and v, to generate two independent standard Gaussian variates:       //
//     sqrt(-2 ln(u) ) cos(2 pi v) and sqrt(-2 ln(u) ) sin(2 pi v).           //
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
//     x = Gaussian_Variate_Box_Muller();                                     //
////////////////////////////////////////////////////////////////////////////////

double Gaussian_Variate_Box_Muller( void )
{
   static double next_variate;
   static int use_next_variate = 0;
   double u, v;
   double variate;
   
   if (use_next_variate) variate = next_variate;
   else {
      u = Uniform_0_1_Random_Variate();
      if (u == 0.0) u = DBL_MIN;
      u = sqrt(-2.0 * log(u) );
      v = Uniform_0_1_Random_Variate();
      variate = u * cos(M_2PI * v);
      next_variate = u * sin(M_2PI * v);
   }
   use_next_variate = !use_next_variate;
   return variate;
}
