////////////////////////////////////////////////////////////////////////////////
// File: uniform_0_1_variate.c                                                //
// Routine(s):                                                                //
//    Uniform_0_1_Random_Variate                                              //
//    Uniform_32_Bits_Random_Variate                                          //
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <time.h>

void Uniform_0_1_Init_Seed( unsigned long seed )
{
   srand((unsigned int)seed);
}

void Uniform_0_1_Init_Time( void )
{
   srand(time(NULL));
}
////////////////////////////////////////////////////////////////////////////////
// double Uniform_0_1_Random_Variate( void )                                  //
//                                                                            //
//  Description:                                                              //
//     This function returns a uniform variate on [0,1].                      //
//                                                                            //
//  Arguments:                                                                //
//     None                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     A random number with a uniform distribution between 0 and 1.           //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Uniform_0_1_Random_Variate();                                      //
////////////////////////////////////////////////////////////////////////////////
        
double Uniform_0_1_Random_Variate( void )
{
   return (double) rand() / ((double) RAND_MAX + 1.0);
}


////////////////////////////////////////////////////////////////////////////////
// unsigned long Uniform_32_Bits_Random_Variate( void )                       //
//                                                                            //
//  Description:                                                              //
//     This function returns a random integer on [0, 2^32].                   //
//                                                                            //
//  Arguments:                                                                //
//     None                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     A uniform random number r, 0 <= r <= 2^32.                             //
//                                                                            //
//  Example:                                                                  //
//     unsigned long x;                                                       //
//                                                                            //
//     x = Uniform_32_Bits_Random_Variate( void )                             //
////////////////////////////////////////////////////////////////////////////////
        
unsigned long Uniform_32_Bits_Random_Variate( void )
{
   unsigned long x = (unsigned long) rand();
   if (rand() < 1073741824) x |= 0x80000000;
   return x;
}
