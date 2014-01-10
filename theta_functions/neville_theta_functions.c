////////////////////////////////////////////////////////////////////////////////
// File: neville_theta_functions.c                                            //
// Routine(s):                                                                //
//    Neville_Theta_Functions                                                 //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for exp(), sin(), cos() and M_PI

//                         Externally Defined Routines                        //
void Jacobian_Theta_Functions(double z, double q, double *theta_1,
                          double *theta_2, double *theta_3, double *theta_4 );

////////////////////////////////////////////////////////////////////////////////
// void Neville_Theta_Functions(double u, double k, double K, double r_tau,   //
//       double *theta_s, double *theta_c, double *theta_d, double *theta_n ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     This function calculates the four Neville theta functions as functions //
//     of u and the nome q = exp(-pi * K' / K) = exp (- pi / r_tau).          //
//      theta_s(u,q) = (2*K/pi) theta_1((pi/2K)u,q) / dtheta_1( 0, q),        //
//      theta_c(u,q) = theta_2((pi/2K)u, q) / theta_2(0,q)                    //
//      theta_d(u,q) = theta_3((pi/2K)u, q) / theta_3(0,q)                    //
//      theta_n(u,q) = theta_4((pi/2K)u, q) / theta_4(0,q)                    //
//     where theta_1, theta_2, theta_3, and theta_4 are the Jacobian theta    //
//     functions as calculated by the subroutine Jacobian_Theta_Functions().  //
//                                                                            //
//     Note that:                                                             //
//        theta_2(0,q) = sqrt( 2 k K / pi )                                   //
//        theta_3(0,q) = sqrt( 2 K / pi )                                     //
//        theta_4(0,q) = sqrt( 2 k' K / pi )                                  //
//      and dtheta_1(0,q) = theta_2(0,q) * theta_3(0,q) * theta_4(0,q),       //
//      where q = exp(-pi*K'/K), k' the complementary modulus.                //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The value of the incomplete elliptic integral of the first  //
//                kind with modulus k.                                        //
//     double  k                                                              //
//                The modulus of the elliptic integral of the first kind.     //
//     double  K                                                              //
//                The complete elliptic integral of the first kind with       //
//                modulus k.                                                  //
//     double  r_tau                                                          //
//                The ratio K / K' where K' is the complete elliptic integral //
//                of the first kind with modulus sqrt(1-k*k) unless k = 0.0   //
//                in which case r_tau is 0.0.                                 //
//     double* theta_s                                                        //
//                The value of theta_s as calculated above.                   //
//     double* theta_c                                                        //
//                The value of theta_c as calculated above.                   //
//     double* theta_d                                                        //
//                The value of theta_d as calculated above.                   //
//     double* theta_n                                                        //
//                The value of theta_n as calculated above.                   //
//                                                                            //
//  Return Value:                                                             //
//     The values of theta_s, theta_c, theta_d, and theta_n are returned      //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double u, k, K, r_tau, cK;                                             //
//     double theta_s, theta_c, theta_d, theta_n;                             //
//                                                                            //
//     ( code to initialize k )                                               //
//                                                                            //
//     K = Complete_Elliptic_Integral_First_Kind(k);                          //
//     if ( k == 0.0 ) r_tau = 0.0;                                           //
//     else {                                                                 //
//        cK = Complete_Elliptic_Integral_First_Kind(sqrt(1.0 - k * k));      //
//        r_tau = K / cK;                                                     //
//     }                                                                      //
//                                                                            //
//     (code to select various values of u)                                   //
//        Neville_Theta_Functions( u, k, K, r_tau, &theta_s, &theta_c,        //
//                                                      &theta_d, &theta_n ); //
////////////////////////////////////////////////////////////////////////////////
                
void Neville_Theta_Functions(double u, double k, double K, double r_tau,
          double *theta_s, double *theta_c, double *theta_d, double *theta_n )
{
   double scale = M_PI / (K + K);
   double z = scale * u;
   double q = (k == 0.0 ) ? 0.0 : exp(- M_PI / r_tau);
   double th0[4];
   double thz[4];

   if (q == 0.0) {
      *theta_s = sin(u);
      *theta_c = cos(u);
      *theta_d = 1.0;
      *theta_n = 1.0;
   }
   else {
      Jacobian_Theta_Functions(z, q, &thz[0], &thz[1], &thz[2], &thz[3]);
  
      th0[1] = sqrt(2.0 * fabs(k) * K / M_PI);
      th0[2] = sqrt(2.0 * K / M_PI);
      th0[3] = sqrt(2.0 * sqrt(1.0 - k * k) * K / M_PI);
      th0[0] = th0[1] * th0[2] * th0[3];
      *theta_s = thz[0] / (scale * th0[0]);
      *theta_c = thz[1] / th0[1];
      *theta_d = thz[2] / th0[2];
      *theta_n = thz[3] / th0[3];
   }
   return; 
}
