////////////////////////////////////////////////////////////////////////////////
// File: jacobi_sn.c                                                          //
// Routine(s):                                                                //
//    Jacobi_sn                                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                     // required for sin()

//                         Externally Defined Routines                        //
extern double Jacobi_am(double u, char arg,  double x);

////////////////////////////////////////////////////////////////////////////////
// double Jacobi_sn(double u, char arg, double x)                             //
//                                                                            //
//  Description:                                                              //
//     Let F(phi,k) = F(phi \ alpha) = F(phi | m) be Legendre's elliptic      //
//     function of the first kind with modulus k, modular angle alpha where   //
//     k = sin(alpha) or parameter m where m = k^2, i.e.                      //
//        F(phi,k) = Integral(0,phi) dtheta / sqrt(1 - k^2 sin^2(theta))      //
//        F(phi \ alpha) = Integral(0,phi) dtheta /                           //
//                                        sqrt(1 - sin^2(alpha) sin^2(theta)) //
//        F(phi | m) = Integral(0,phi) dtheta / sqrt(1 - m sin^2(theta))      //
//                                                                            //
//     This Jacobi elliptic amplitude function, am, is defined as             //
//               am(u,k) = am(u \ alpha) = am(u | m)  = phi                   //
//     where u = F(phi,k) = F(phi \ alpha) = F(phi | m).                      //
//                                                                            //
//     This function calculates the Jacobi elliptic sine amplitude function,  //
//     sn, which is defined as the sine of the Jacobi elliptic amplitude      //
//     function, am, i.e.                                                     //
//                sn(u,k) = sin(am(u,k)),                                     //
//                sn(u \ alpha) = sin(am(u \ alpha),                          //
//                sn(u | m) = sin(am(u | m).                                  //
//                                                                            //
//     Restrictions:  -1 <= k <= 1 or 0 <= m <= 1.                            //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The first argument of sn(u,x) corresponding to the value of //
//                the elliptic integral of the first kind u = F(am(u,x),x).   //
//     char    arg                                                            //
//                The type of argument of the second argument of sn():        //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the sine amplitude function sn(u,x)  //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.  If 'arg' = 'm', then x must be between 0 and 1   //
//                inclusively and if 'arg' = 'k', then x must be between -1   //
//                and 1 inclusively.                                          //
//                                                                            //
//  Return Value:                                                             //
//     Jacobi's elliptic function sn().                                       //
//                                                                            //
//  Example:                                                                  //
//     double u, x;                                                           //
//     double sn;                                                             //
//     char   arg;                                                            //
//                                                                            //
//     ( code to initialize u, arg, and x )                                   //
//                                                                            //
//     sn = Jacobi_sn( u, arg, x);                                            //
////////////////////////////////////////////////////////////////////////////////

double Jacobi_sn(double u, char arg,  double x)
{
   return sin( Jacobi_am(u, arg, x) );
}
