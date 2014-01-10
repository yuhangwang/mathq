////////////////////////////////////////////////////////////////////////////////
// File: dneville_theta_functions.c                                           //
// Routine(s):                                                                //
//    DNeville_Theta_Functions                                                //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for exp(), sin(), cos() and M_PI

//                         Externally Defined Routines                        //
void Jacobian_Theta_Functions_at_zero(double q, double *theta_1,
                           double *theta_2, double *theta_3, double *theta_4 );
void DJacobian_Theta_Functions(double z, double q, double *dtheta_1,
                        double *dtheta_2, double *dtheta_3, double *dtheta_4 );

////////////////////////////////////////////////////////////////////////////////
// void DNeville_Theta_Functions(double u, double k, double K, double r_tau,  //
//   double *dtheta_s, double *dtheta_c, double *dtheta_d, double *dtheta_n ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     This function calculates the derivatives with respect to u of the      //
//     four Neville theta functions as functions of u and the nome            //
//     q = exp(-pi * K' / K) = exp (- pi / r_tau).                            //
//      dtheta_s(u,q) = dtheta_1((pi/2K)u,q) / dtheta_1( 0, q),               //
//      dtheta_c(u,q) = (pi / 2K) dtheta_2((pi/2K)u, q) / theta_2(0,q)        //
//      dtheta_d(u,q) = (pi / 2K) dtheta_3((pi/2K)u, q) / theta_3(0,q)        //
//      dtheta_n(u,q) = (pi / 2K) dtheta_4((pi/2K)u, q) / theta_4(0,q)        //
//     where theta_1, theta_2, theta_3, and theta_4 are the Jacobian theta    //
//     functions and dtheta_1, dtheta_2, dtheta_3, and dtheta_4 are the       //
//     derivative os the Jacobian theta functions as calculated by the        //
//     subroutine DJacobian_Theta_Functions().                                //
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
//     double* dtheta_s                                                       //
//                The value of dtheta_s as calculated above.                  //
//     double* dtheta_c                                                       //
//                The value of dtheta_c as calculated above.                  //
//     double* dtheta_d                                                       //
//                The value of dtheta_d as calculated above.                  //
//     double* dtheta_n                                                       //
//                The value of dtheta_n as calculated above.                  //
//                                                                            //
//  Return Value:                                                             //
//     The values of dtheta_s, dtheta_c, dtheta_d, and dtheta_n are returned  //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double u, k, K, r_tau, cK;                                             //
//     double dtheta_s, dtheta_c, dtheta_d, dtheta_n;                         //
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
//        DNeville_Theta_Functions( u, k, K, r_tau, &dtheta_s, &dtheta_c,     //
//                                                    &dtheta_d, &dtheta_n ); //
////////////////////////////////////////////////////////////////////////////////
                
void DNeville_Theta_Functions(double u, double k, double K, double r_tau,
       double *dtheta_s, double *dtheta_c, double *dtheta_d, double *dtheta_n )
{
   double coef = M_PI / (K + K);
   double z = coef * u;
   double q = (k == 0.0 ) ? 0.0 : exp(- M_PI / r_tau);
   double th0[4];
   double thz[4];

   if (q == 0.0) {
      *dtheta_s = cos(u);
      *dtheta_c = -sin(u);
      *dtheta_d = 0.0;
      *dtheta_n = 0.0;
   }
   else {
      DJacobian_Theta_Functions(z, q, &thz[0], &thz[1], &thz[2], &thz[3]);
      Jacobian_Theta_Functions_at_zero(q, &th0[0], &th0[1], &th0[2], &th0[3]);
      th0[0] = th0[1] * th0[2] * th0[3];
      *dtheta_s = thz[0] / th0[0];
      *dtheta_c = coef * thz[1] / th0[1];
      *dtheta_d = coef * thz[2] / th0[2];
      *dtheta_n = coef * thz[3] / th0[3];
   }
   return; 
}
