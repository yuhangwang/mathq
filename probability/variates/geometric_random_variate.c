////////////////////////////////////////////////////////////////////////////////
// File: geometric_random_variate.c                                           //
// Routine(s):                                                                //
//    Geometric_Random_Variate                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for log()

//                    Required Externally Defined Routines                    //

extern double Exponential_Random_Variate( void );

////////////////////////////////////////////////////////////////////////////////
// int Geometric_Random_Variate( double p )                                   //
//                                                                            //
//  Description:                                                              //
//     This function returns a Geometric distributed random variate.          //
//                                                                            //
//  Arguments:                                                                //
//     double p                                                               //
//            The probability of a success, 0 < p < 1.                        //
//                                                                            //
//  Return Values:                                                            //
//     A random number distributed as a Geometric distribution.               //
//                                                                            //
//  Example:                                                                  //
//     int x;                                                                 //
//     double p;                                                              //
//                                                                            //
//                  (* Set the probability p, 0 < p < 1 *)                    //
//                                                                            //
//     x = Geometric_Random_Variate( p );                                     //
////////////////////////////////////////////////////////////////////////////////
        
int Geometric_Random_Variate( double p )
{
   double y = -Exponential_Random_Variate( );

   return (int) (y / log(1.0 - p)); 
}
