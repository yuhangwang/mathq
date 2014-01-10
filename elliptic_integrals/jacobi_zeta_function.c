////////////////////////////////////////////////////////////////////////////////
// File: jacobi_zeta_function.c                                               //
// Routine(s):                                                                //
//    Jacobi_Zeta_Function                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// double Jacobi_Zeta_Function(double amplitude, char arg, double x )         //
//                                                                            //
//  Description:                                                              //
//     This function calculates Jacobi's zeta function, Z(phi,k),             //
//               Z(phi,k) = E(phi,k) - E(k) * F(phi,k) / K(k),                //
//     where F(phi,k) is Legendre's elliptic integral of the first kind,      //
//     E(phi,k) is Legendre's elliptic integral of the second kind, K(k) is   //
//     the complete elliptic integral of the first kind, and E(k) is the      //
//     complete elliptic integral of the second kind.  The amplitude is the   //
//     amplitude (phi) of both the Legendre elliptic integral of the first    //
//     and second kinds and is subject to the restriction that |phi| <= pi/2. //
//     The modulus is the modulus (k) of all the elliptic integrals and is    //
//     subject to the restriction that |k| <= 1.                              //
//     Note that k = sin(alpha) where alpha is the modular angle and that     //
//     m = k^2 is the parameter of the elliptic integrals, i.e.               //
//                 Z(phi, k) = Z(phi | m) = Z(phi \ alpha), or                //
//                 Z(phi | m) = Z(phi, sqrt(m)) and                           //
//                 Z(phi \ alpha) = Z(phi, sin(alpha)).                       //
//                                                                            //
//  Arguments:                                                                //
//     double  amplitude                                                      //
//                The amplitude of elliptic integrals, |amplitude| <= pi / 2. //
//     char    arg                                                            //
//                The type of argument of the second argument of Z():         //
//                  If arg = 'k', then x = k, the modulus of Z(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                Z(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of Z(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the zeta function Z(phi,x),          //
//                Z(phi \ alpha) or Z(phi | m) corresponding to the value of  //
//                'arg'.  Note that |k| <= 1 or 0 <= m <= 1.                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the zeta function evaluated at amplitude = phi and        //
//     modulus = k (or parameter m = k^2 or modular angle alpha = arcsin(k)). //
//                                                                            //
//  Example:                                                                  //
//     double phi, alpha;                                                     //
//     double zeta;                                                           //
//                                                                            //
//     ( code to initialize phi and alpha )                                   //
//                                                                            //
//     zeta = Jacobi_Zeta_Function( phi,'a',alpha);                           //
//     printf("Jacobi's Zeta function %12.6f\n",zeta);                        //
//     zeta = Jacobi_Zeta_Function( phi,'k',sin(alpha));                      //
//     printf("Jacobi's Zeta function %12.6f\n",zeta);                        //
//     zeta = Jacobi_Zeta_Function( phi,'m',sin(alpha)*sin(alpha));           //
//     printf("Jacobi's Zeta function %12.6f\n",zeta);                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for fabs(), sqrt(), log(), atanl(),
                            // tanl(), sqrtl(), fabsl() and M_PI_2 = pi/2.
#include <float.h>          // required for LDBL_EPSILON, DBL_MAX

double Jacobi_Zeta_Function(double amplitude, char arg, double x )
{
   long double phi;
   long double m;
   long double two_n;
   long double g;
   long double a;
   long double g_old;
   long double a_old;
   long double sum;
   long double tan_2n_phi;
   int sgn_amplitude = (amplitude >= 0.0) ? 1 : -1;

                        // Check for special cases. //

   if ( amplitude == 0.0 ) return 0.0;
   if ( fabs(amplitude) == M_PI_2 ) return 0.0;
   
   if ( x == 0.0 ) return 0.0;

     // Convert modulus, modular angle, or parameter to the parameter. //

   switch (arg) {
      case 'k': m = fabsl((long double) x);
                m *= m;
                break;
      case 'm': m = (long double) x;
                break;
      case 'a': m = sinl((long double)x);
                m *= m;
                break;
      default:  m = fabsl((long double) x);
                m *= m;
   }

                   // Check for special case |k| = m = 1. //

   if ( m == 1.0 ) return sin(amplitude);

           // Perform Arithmetic-Geometric Common Mean Transform. //
   
   phi = (long double) fabs(amplitude);
   two_n = 1.0L;
   g = sqrtl(1.0L - m);
   a = 1.0L;
   sum = 0.0L;

   while (1) {
      tan_2n_phi = tanl( phi );
      two_n += two_n;
      phi += phi 
          - atanl( (a - g) * tan_2n_phi / (a + g * tan_2n_phi * tan_2n_phi) );
      sum += (a - g) * sinl(phi);
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = sqrtl(g_old * a_old);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }
   phi /= two_n;
   return (double) (sgn_amplitude * (0.5L * sum)); 
}
