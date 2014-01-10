////////////////////////////////////////////////////////////////////////////////
// File: heumans_lambda_naught.c                                              //
// Routine(s):                                                                //
//    Heumans_Lambda_Naught                                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>          // required for sin(), sqrt(), fabs(), and M_PI_2

//                         Externally Defined Routines                        //

extern void Complete_Elliptic_Integrals(char arg, double x, double* Fk,
                                                                  double* Ek);
extern void Legendre_Elliptic_Integrals(double amplitude, char arg, double x, 
                                 double* F, double* K, double* E, double* Ek);

////////////////////////////////////////////////////////////////////////////////
// void Heumans_Lambda_Naught(double amplitude, double modular_angle)         //
//                                                                            //
//  Description:                                                              //
//     This function calculates Heuman's lambda naught function               //
//      ^o(phi,alpha) = 2/pi{ K(alpha) * E(phi \ pi/2-alpha)                  //
//                          - [K(alpha) - E(alpha)] * F(phi, \ pi/2-alpha) }, //
//     where F(phi,k) is Legendre's elliptic integral of the first kind,      //
//     E(phi,k) is Legendre's elliptic integral of the second kind, K(k) is   //
//     the complete elliptic integral of the first kind, and E(k) is the      //
//     complete elliptic integral of the second kind.  The amplitude is the   //
//     amplitude (phi) of both the Legendre elliptic integral of the first    //
//     and second kinds and is subject to the restriction that |phi| <= pi/2. //
///    The modular angle is the modular angle of the complete elliptic        //
//     integrals and the complementary modular angle of the Legendre elliptic //
//     integrals.  Note that k = sin(modular_angle) is the modulus of the     //
//     complete elliptic integrals and k' = sqrt(1-k^2) is the modulus of the //
//     incomplete elliptic integrals.                                         //
//                                                                            //
//  Arguments:                                                                //
//     double  amplitude                                                      //
//                The amplitude of elliptic integrals, |amplitude| <= pi / 2. //
//     double  modular_angle                                                  //
//                The parameter, k = sin(modular_angle), is the modulus of the//
//                complete elliptic integrals and the parameter,              //
//                k' = sin(pi/2 - modular_angle), is the modulus of the       //
//                incomplete elliptic integrals, modular_angle <= pi / 2.     //
//                                                                            //
//  Return Value:                                                             //
//     The value of the Heuman's Lambda naught function, given above,         //
//     evaluated at the amplitude = phi and modular angle alpha.              //
//                                                                            //
//  Example:                                                                  //
//     double phi, alpha;                                                     //
//     double lambda;                                                         //
//                                                                            //
//     ( code to initialize phi and alpha )                                   //
//                                                                            //
//     lambda = Heumans_Lambda_Naught( phi, alpha );                          //
//     printf("Heuman's Lambda Naught function %12.6f\n",lambda);             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

double Heumans_Lambda_Naught(double amplitude, double modular_angle )
{
   double k;               // Modulus = sin(modular_angle)
   double ck;              // Complementary Modulus = cos(modular_angle)
   double K, Ek;           // Complete elliptic integrals modulus k
   double cF, cE;          // Elliptic integrals modulus sqrt(1-k^2)
   int sign_amplitude = (amplitude < 0.0) ? -1 : 1;

                        // Check for special cases. //

   if (modular_angle == 0.0) return sin(amplitude);
   if (modular_angle == M_PI_2) return amplitude / M_PI_2; 
   if (amplitude == 0.0) return 0.0;
   if (fabs(amplitude) == M_PI_2) return (double) sign_amplitude;

        // Calculate the complete elliptic integrals with modulus k. //
        // And the elliptic integrals with complementary modulus k'  //
        // and amplitude 'amplitude'.                                //

   k = sin(modular_angle);
   ck = sqrt(1.0 - k * k);
   Legendre_Elliptic_Integrals(amplitude,'k', ck, &cF, &K, &cE, &Ek);
   Complete_Elliptic_Integrals('k', k, &K, &Ek);
   
   return (K * cE - (K - Ek) * cF) / M_PI_2;
}
