////////////////////////////////////////////////////////////////////////////////
// File: djacobian_theta_functions.c                                          //
// Routine(s):                                                                //
//    DJacobian_Theta_Functions                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for log() and M_PI 

//                         Externally Defined Routines                        //
void DTheta_Functions(double nu, double x, double *dtheta_1, double *dtheta_2,
                                          double *dtheta_3, double *dtheta_4 );

////////////////////////////////////////////////////////////////////////////////
// void DJacobian_Theta_Functions(double z, double q, double *dtheta_1,       //
//                     double *dtheta_2, double *dtheta_3, double *dtheta_4 ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     The results are not returned in an array in order to avoid the         //
//     potentially confusing dtheta[0]() = dtheta_1(), ..., dtheta_4().       //
//     This function calculates the derivatives with respect to z of the      //
//     four theta functions as functions of z and the nome q, 0 <= q < 1.     //
//      dtheta_1(z,q) = 2q^(1/4) Sum[(-1)^j (2j+1) q^(j(j + 1)) cos[(2j+1)z]] //
//      dtheta_2(z,q) = -2q^(1/4) Sum[q^(j(j+1)) (2j+1) sin[(2j+1)z]]         //
//     where the sum extends from j = 0 to inf and                            //
//      dtheta_3(z,q) = 4 Sum[q^(j^2) j sin[2jz]]                             //
//      dtheta_4(z,q) = 4 Sum[(-1)^j q^(j^2) j sin[2jz]]                      //
//     where the sum extends from j = 1 to inf.                               //
//                                                                            //
//  Arguments:                                                                //
//     double  z                                                              //
//                The periodic argument of the Theta functions.               //
//     double  q                                                              //
//                The nome, note that 0 <= q < 1.                             //
//     double* dtheta_1                                                       //
//                The value of dtheta_1 as calculated above.                  //
//     double* dtheta_2                                                       //
//                The value of dtheta_2 as calculated above.                  //
//     double* dtheta_3                                                       //
//                The value of dtheta_3 as calculated above.                  //
//     double* dtheta_4                                                       //
//                The value of dtheta_4 as calculated above.                  //
//                                                                            //
//  Return Value:                                                             //
//     The values of dtheta_1, dtheta_2, dtheta_3, and dtheta_4 are returned  //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double z, q;                                                           //
//     double dtheta_1, dtheta_2, dtheta_3, dtheta_4;                         //
//                                                                            //
//     ( code to initialize z and q )                                         //
//                                                                            //
//     DJacobian_Theta_Functions(z, q, &dtheta_1, &dtheta_2, &dtheta_3,       //
//                                                                &dtheta_4); //
////////////////////////////////////////////////////////////////////////////////
                
void DJacobian_Theta_Functions(double z, double q, double *dtheta_1,
                         double *dtheta_2, double *dtheta_3, double *dtheta_4 )
{
   double nu = z / M_PI;
   double x;

   if ( q == 0.0 ) {
      *dtheta_1 = 0.0;
      *dtheta_2 = 0.0;
      *dtheta_3 = 0.0;
      *dtheta_4 = 0.0; 
      return;     
   }
   else {
      x = -log(q) / (M_PI * M_PI);
      DTheta_Functions(nu, x, dtheta_1, dtheta_2, dtheta_3, dtheta_4 );
   }

   *dtheta_1 /= M_PI;
   *dtheta_2 /= M_PI;
   *dtheta_3 /= M_PI;
   *dtheta_4 /= M_PI;

   return; 
}
