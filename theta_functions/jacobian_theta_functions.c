////////////////////////////////////////////////////////////////////////////////
// File: jacobian_theta_functions.c                                           //
// Routine(s):                                                                //
//    Jacobian_Theta_Functions                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for log() and M_PI 

//                         Externally Defined Routines                        //
void Theta_Functions(double nu, double x, double *theta_1, double *theta_2,
                                           double *theta_3, double *theta_4 );

////////////////////////////////////////////////////////////////////////////////
// void Jacobian_Theta_Functions(double z, double q, double *theta_1,         //
//                        double *theta_2, double *theta_3, double *theta_4 ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     The results are not returned in an array in order to avoid the         //
//     potentially confusing theta[0]() = theta_1(), ... .                    //
//     This function calculates the four theta functions as functions of z and//
//     the nome q, 0 <= q < 1.                                                //
//      theta_1(z,q) = 2q^(1/4) Sum[(-1)^j q^(j(j + 1)) sin[(2j+1)z]]         //
//      theta_2(z,q) = 2q^(1/4) Sum[q^(j(j+1)) cos[(2j+1)z]]                  //
//     where the sum extends from j = 0 to inf and                            //
//      theta_3(z,q) = 1 + 2 Sum[q^(j^2) cos[2jz]]                            //
//      theta_4(z,q) = 1 + 2 Sum[(-1)^j q^(j^2) cos[2jz]]                     //
//     where the sum extends from j = 1 to inf.                               //
//                                                                            //
//  Arguments:                                                                //
//     double  z                                                              //
//                The periodic argument of the Theta functions.               //
//     double  q                                                              //
//                The nome, note that 0 <= q < 1.                             //
//     double* theta_1                                                        //
//                The value of theta_1 as calculated above.                   //
//     double* theta_2                                                        //
//                The value of theta_2 as calculated above.                   //
//     double* theta_3                                                        //
//                The value of theta_3 as calculated above.                   //
//     double* theta_4                                                        //
//                The value of theta_4 as calculated above.                   //
//                                                                            //
//  Return Value:                                                             //
//     The values of theta_1, theta_2, theta_3, and theta_4 are returned      //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double z, q;                                                           //
//     double theta_1, theta_2, theta_3, theta_4;                             //
//                                                                            //
//     ( code to initialize z and q )                                         //
//                                                                            //
//     Jacobian_Theta_Functions(z, q, &theta_1, &theta_2, &theta_3, &theta_4);//
////////////////////////////////////////////////////////////////////////////////
                
void Jacobian_Theta_Functions(double z, double q, double *theta_1,
                           double *theta_2, double *theta_3, double *theta_4 )
{
   double nu = z / M_PI;
   double x;

   if ( q == 0.0 ) {
      *theta_1 = 0.0;
      *theta_2 = 0.0;
      *theta_3 = 1.0;
      *theta_4 = 1.0;      
   }
   else {
      x = -log(q) / (M_PI * M_PI);
      Theta_Functions(nu, x, theta_1, theta_2, theta_3, theta_4 );
   }

   return; 
}
