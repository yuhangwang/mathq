////////////////////////////////////////////////////////////////////////////////
// File: sum_12_uniforms.c                                                    //
// Routine(s):                                                                //
//    Gaussian_Variate_Sum_12_Uniforms                                        //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines                    //
        
extern double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Variate_Sum_12_Uniforms( void )                            //
//                                                                            //
//  Description:                                                              //
//     This function returns a Gaussian variate with mean 0 and standard      //
//     deviation 1 by applying the central limit theorem to a sum of          //
//     independent uniform (0,1) distributed random variables.                //
//                                                                            //
//     If Sn is the sum of n independent uniformly (0,1) distributed random   //
//     variables then the random variable                                     //
//                        ( Sn - n/2 ) / sqrt( n/12 )                         //
//     converges in distribution to a standard normal as n becomes large.     //
//                                                                            //
//     This routine choose n = 12.  This routine is not recommended as        //
//     both the Box-Muller method, the Polar-Marsaglia method, and            //
//     Marsaglia's ziggurat method are superior.                              //
//                                                                            //
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
//     x = Gaussian_Variate_Sum_12_Uniforms();                                //
////////////////////////////////////////////////////////////////////////////////

double Gaussian_Variate_Sum_12_Uniforms( void )
{
   double z = 0.0;
   int i;
   
   for (i = 0; i < 12; i++) z += Uniform_0_1_Random_Variate();
   z -= 6.0;

   return z;
}
