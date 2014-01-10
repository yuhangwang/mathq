////////////////////////////////////////////////////////////////////////////////
// File: jacobi_cn.c                                                          //
// Routine(s):                                                                //
//    Jacobi_cn                                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                     // required for cos()

//                         Externally Defined Routines                        //
extern double Jacobi_am(double u, char arg,  double x);

////////////////////////////////////////////////////////////////////////////////
// double Jacobi_cn(double u, char arg, double x)                             //
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
//     This function calculates the Jacobi elliptic cosine amplitude          //
//     function, cn, which is defined as the cosine of the Jacobi elliptic    //
//     amplitude function, am, i.e.                                           //
//                cn(u,k) = cos(am(u,k)),                                     //
//                cn(u \ alpha) = cos(am(u \ alpha),                          //
//                cn(u | m) = cos(am(u | m).                                  //
//                                                                            //
//     Restrictions:  -1 <= k <= 1 or 0 <= m <= 1.                            //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The first argument of cn(u,x) corresponding to the value of //
//                the elliptic integral of the first kind u = F(am(u,x),x).   //
//     char    arg                                                            //
//                The type of argument of the second argument of cn():        //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the cosine amplitude function cn(u,x)//
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.  If 'arg' = 'm', then x must be between 0 and 1   //
//                inclusively and if 'arg' = 'k', then x must be between -1   //
//                and 1 inclusively.                                          //
//                                                                            //
//  Return Value:                                                             //
//     Jacobi's elliptic function cn().                                       //
//                                                                            //
//  Example:                                                                  //
//     double u, x;                                                           //
//     double cn;                                                             //
//     char   arg;                                                            //
//                                                                            //
//     ( code to initialize u, arg, and x )                                   //
//                                                                            //
//     cn = Jacobi_cn( u, arg, x);                                            //
////////////////////////////////////////////////////////////////////////////////

double Jacobi_cn(double u, char arg,  double x)
{
   return cos( Jacobi_am(u, arg, x) );
}
