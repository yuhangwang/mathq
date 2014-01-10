////////////////////////////////////////////////////////////////////////////////
// File: legendre_elliptic_integrals.c                                        //
// Routine(s):                                                                //
//    Legendre_Elliptic_Integrals                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for fabs(), sqrt(), log(), atanl(),
                            // tanl(), sqrtl(), fabsl() and M_PI_2 = pi/2.
#include <float.h>          // required for LDBL_EPSILON, DBL_MAX

//                         Internally Defined Routines                        //

static void Landen_Transform( long double phi, long double parameter,
                               double *F, double *Fk, double *E, double *Ek);
static void Elliptic_Integrals( long double phi, long double m, double *F,
                                            double *K, double *E, double *Ek);
static void Large_Modulus(double amplitude, long double k, double *F, 
                                            double *K, double *E, double *Ek);

//                         Internally Defined Constants                       //

static const long double PI_4 = 0.7853981633974483096156608458198757L; // pi/4
static const long double PI_2 =  1.5707963267948966192313216916397514L; // pi/2
static const long double PI   =  3.1415926535897932384626433832795029L; // pi

////////////////////////////////////////////////////////////////////////////////
// void Legendre_Elliptic_Integrals(double amplitude, char arg, double x,     //
//                               double* F, double* K, double* E, double* Ek) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the incomplete elliptic integral of the first //
//     Kind, the complete elliptic integral of the first kind, the incomplete //
//     elliptic integral of the second kind, and the complete elliptic        //
//     integral of the second kind.                                           //
//     Legendre's Elliptic Integral of the First Kind, F(phi,k), is the       //
//     integral from 0 to phi of the integrand                                //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) ).                   //
//     If phi = pi / 2, the integral as a function of the modulus k is called //
//     the complete elliptic integral of the first kind.                      //
//     Legendre's Elliptic Integral of the Second Kind, E(phi,k), is the      //
//     integral from 0 to phi of the integrand                                //
//                   sqrt( 1 - k^2 sin^2(theta) ) dtheta .                    //
//     If phi = pi / 2, the integral as a function of the modulus k is called //
//     the complete elliptic integral of the second kind.                     //
//                                                                            //
//     In practice the arguments the elliptic function of the first kind are  //
//     also given as F(phi \ alpha) or F(phi | m) and the elliptic function   //
//     of the second kind are given as E(phi \ alpha) or E(phi | m) where the //
//     angle alpha, called the modular angle, satisfies k = sin(alpha) and    //
//     the argument m = k^2 is simply called the parameter.                   //
//     In terms of these arguments F(phi \ alpha) = F(phi, sin(alpha)),       //
//     F(phi | m) = F(phi, sqrt(m)), E(phi \ alpha) = E(phi, sin(alpha)),     //
//     and E(phi | m) = E(phi sqrt(m)).                                       //
//                                                                            //
//  Arguments:                                                                //
//     double  amplitude                                                      //
//                The upper limit of the integral.                            //
//     char    arg                                                            //
//                The type of argument of the second argument of both F and E.//
//                  If arg = 'k', then x = k, the modulus of F(phi,k) and     //
//                                E(phi,k).                                   //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha) and E(phi \ alpha), where    //
//                                alpha is in radians.                        //
//                  If arg = 'm', then x = m, the parameter of F(phi | m)     //
//                                and E(phi | m).                             //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the elliptic function F(phi,x),      //
//                F(phi \ alpha) or F(phi | m) and E(phi,x), E(phi \ alpha),  //
//                or E(phi | m) corresponding to the value of 'arg'.          //
//     double* F                                                              //
//                The Legendre elliptic integral of the first kind.           //
//                If |modulus| = 1 and |amplitude| = pi / 2, then the integral//
//                is infinite and DBL_MAX is returned if the amplitude = pi/2 //
//                and -DBL_MAX if the amplitude = -pi/2.                      //
//     double* K                                                              //
//                The complete elliptic integral of the first kind.           //
//                If |modulus| = 1, then the integral is infinite and DBL_MAX //
//                returned.                                                   //
//     double* E                                                              //
//                The Legendre elliptic integral of the second kind.          //
//     double* Ek                                                             //
//                The complete elliptic integral of the second kind.          //
//                                                                            //
//  Return Value:                                                             //
//     void.  The values are returned through the argument list.              //
//                                                                            //
//  Example:                                                                  //
//     double f, phi, k;                                                      //
//     double F, K, E, Ek;                                                    //
//                                                                            //
//     ( code to initialize phi and k )                                       //
//                                                                            //
//     Legendre_Elliptic_Integrals( phi, k, &F, &K, &E, &Ek );                //
//     printf("Legendre's elliptic integral of the 1st kind %12.6f\n",F);     //
//     printf("Legendre's elliptic integral of the 2nd kind %12.6f\n",E);     //
//     printf("The complete elliptic integral of the 1st kind %12.6f\n",K);   //
//     printf("The complete elliptic integral of the 2nd kind %12.6f\n",Ek);  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

void Legendre_Elliptic_Integrals(double amplitude, char arg, double x, 
                                 double* F, double* K, double* E, double* Ek)
{
   long double dum;
   long double phi;
   long double k,m;
   int n;
   int sgn_amplitude = (amplitude >= 0.0) ? 1 : -1;

                   // Check for special case: modulus = 0 //

   if ( x == 0.0 ) {
      *F = amplitude;
      *E = amplitude;
      *K = M_PI_2;
      *Ek = M_PI_2;
      return;
   }

 // Convert modulus, modular angle, or parameter to modulus and parameter. //

   switch (arg) {
      case 'k': k = fabsl((long double) x);
                m = k * k;
                break;
      case 'm': m = (long double) x;
                k = sqrtl(fabsl(m));
                break;
      case 'a': k = sinl((long double)x);
                m = k * k;
                break;
      default:  k = fabsl((long double) x);
                m = k * k;
   }

                  // Check for most common case 0 < m < 1. //
                        // Check for special cases. //

   if ( m > 0.0 && m < 1.0 ) {
      Elliptic_Integrals( fabsl((long double) amplitude), m, F, K, E, Ek);
      *F = sgn_amplitude * *F;
      *E = sgn_amplitude * *E;
      return;
   }

          // Check for case m < 0 i.e. a purely imaginary modulus. //
   
   if (m < 0.0 ) {
      phi = PI_2 - fabsl((long double) amplitude);
      Elliptic_Integrals( fabsl(phi), fabsl(m / (1.0L - m)), F, K, E, Ek);
      dum = sqrtl(1.0L - m);
      if (phi > PI_2) {
         *F = (double) (sgn_amplitude * (*K + *F) / dum);
         *E = (double) (sgn_amplitude * (*Ek + *E) * dum);
      }
      else {
         *F = (double) (sgn_amplitude * (*K - *F) / dum);
         *E = (double) (sgn_amplitude * (*Ek - *E) * dum);
      }
      *K /= dum;
      *Ek *= dum;
      return;
   }

                        // Check for the case m = 1. //

   if ( m == 1.0 ) {
      if ( fabs(amplitude) >= M_PI_2) {
         *F = sgn_amplitude * DBL_MAX;
         n = (int) ( (fabs(amplitude) + M_PI_2) / M_PI );
         *E = (double) (n + n) + sin(fabs(amplitude) - n * M_PI);
         *E = sgn_amplitude * *E;
      }
      else {
         x = tan(amplitude);
         *F = sgn_amplitude * (log(fabs(x) + sqrt(1.0 + x * x)));
         *E = sin(amplitude);
      }
      *K = DBL_MAX;
      *Ek = 1.0;
      return;                                            
   }

                   // Check for the remaining case m > 1. //

   Large_Modulus(fabs(amplitude), k, F, K, E, Ek);
   *F *= sgn_amplitude;
   *E *= sgn_amplitude;
   return;
}

////////////////////////////////////////////////////////////////////////////////
// static void Elliptic_Integrals( long double phi, long double m, double *F, //
//                                         double *K, double *E, double *Em)  //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the incomplete elliptic integral of the first  //
//     kind F(phi | m), the complete elliptic integral of the first kind K(m) //
//     = F(pi/2 | m), the incomplete ellptic integral of the second kind      //
//     E(phi | m) and the complete elliptic integral of the second kind       //
//     Em(m) = E(pi/2 | m) for 0 < m < 1.                                     //
//                                                                            //
//  Arguments:                                                                //
//     long double phi                                                        //
//                The upper limit of the integral.                            //
//     long double m                                                          //
//                The parameter of the elliptic integrals of the first and    //
//                second kinds. 0 < m < 1.                                    //
//     double *F                                                              //
//                The incomplete elliptic integral of the first kind,         //
//                F(phi | m).                                                 //
//     double *K                                                              //
//                The complete elliptic integral of the first kind,           //
//                K(m) = F(pi/2 | m).                                         //
//     double *E                                                              //
//                The incomplete elliptic integral of the second kind,        //
//                E(phi | m).                                                 //
//     double *Em                                                             //
//                The complete elliptic integral of the second kind,          //
//                Em(m) = E(pi/2 | m).                                        //
//                                                                            //
//  Return Value:                                                             //
//     The values of the incomplete and complete elliptic integrals of the    //
//     first and second kind for the given parameter m, 0 < m < 1, are        //
//     returned via the argument list.                                        //
//                                                                            //
//  Example:                                                                  //
//     double F,K.E,Em;                                                       //
//     long double phi, m;                                                    //
//                                                                            //
//     ( code to initialize phi and m )                                       //
//                                                                            //
//     Elliptic_Integrals( phi, m, &F, &K, &E, &Em);                          //
//     printf("Legendre's elliptic integral of the 1st kind %12.6f\n",F);     //
//     printf("Legendre's elliptic integral of the 2nd kind %12.6f\n",E);     //
//     printf("The complete elliptic integral of the 1st kind %12.6f\n",K);   //
//     printf("The complete elliptic integral of the 2nd kind %12.6f\n",Em);  //
////////////////////////////////////////////////////////////////////////////////

static void Elliptic_Integrals( long double phi, long double m, double *F,
                                             double *K, double *E, double *Em) 
{
   int n;

   n = (int) ( ( phi + PI_2 ) / PI );
   phi -= n * PI;
   n += n;
   Landen_Transform( fabsl(phi), m, F, K, E, Em);
   if (phi >= 0.0) { *F +=  n * *K; *E +=  n * *Em; }
   else { *F = n * *K - *F; *E = n * *Em - *E; }

   return;
}

////////////////////////////////////////////////////////////////////////////////
// static void Large_Modulus(double amplitude, long double k, double *F,      //
//                                          double *K, double *E, double *Ek) //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the real part of Legendre's Elliptic Integral  //
//     of the Second Kind, Re[ E(amplitude,k) ] for  k > 1.                   //
//                                                                            //
//  Arguments:                                                                //
//     double amplitude                                                       //
//                The upper limit of the integral.                            //
//     long double k                                                          //
//                The modulus of the Elliptic Integral of the second kind.    //
//                k > 1.                                                      //
//                                                                            //
//  Return Value:                                                             //
//     The real part of the value of the incomplete elliptic integral of the  //
//     first and second kinds for the given modulus k, k > 1 and              //
//     and amplitude > 0 and the real part of the complete elliptic integrals //
//     of the first and second kind.                                          //
//                                                                            //
//  Example:                                                                  //
//     double F,K.E,Ek;                                                       //
//     long double phi, k;                                                    //
//                                                                            //
//     ( code to initialize phi > 0 and k > 1 )                               //
//                                                                            //
//     Large_Modulus( phi. k, &F, &K, &E, &Ek);                               //
//     printf("Re[F(phi,k)] = %12.6f where phi = %12.6Lf, k =                 //
//                                                    %12.6f\n",F, phi, k);   //
//     printf("Re[E(phi,k)] = %12.6f where phi = %12.6Lf, k =                 //
//                                                    %12.6f\n",E, phi, k);   //
////////////////////////////////////////////////////////////////////////////////

 static void Large_Modulus(double amplitude, long double k, double *F,
                                             double *K, double *E, double *Ek)
{
   long double m = k * k;
   long double phi = (long double) amplitude;
   long double sin_phi;
   int n;

   n = (int) ( ( phi + PI_2 ) / PI );
   phi -= n * PI;
   n += n;
   
   sin_phi = sinl(phi);
   if ( fabsl(sin_phi) >= 1.0L / k ) 
      if (phi > 0.0L) phi = PI_2;
      else phi = -PI_2;
   else phi = asinl(k * sin_phi);
   Landen_Transform( fabsl(phi), 1.0L/m, F, K, E, Ek);
   *Ek = k * *Ek + (1.0L - m) * *K / k;
   *E = k * *E + (1.0 - m) * *F / k;
   if (phi >= 0.0L) { *F += n * *K; *E +=  n * *Ek; }
   else { *F = n * *K - *F; *E = n * *Ek - *E; }
   *F /= k;
   *K /= k; 
   return;
}

////////////////////////////////////////////////////////////////////////////////
// static void Landen_Transform( long double phi, long double parameter,      //
//                              double *F, double *K, double *E, double *Em)  //
//                                                                            //
//  Description:                                                              //
//     The common mean method, sometimes called the Gauss transform, is a     //
//     variant of the descending Landen transformation in which two sequences //
//     are formed: Setting a[0] = 1 and g[0] = k', the complementary modulus, //
//     a[i] is the arithmetic average and g[i] is the geometric mean of a[i-1]//
//     and g[i-1], i.e. a[i+1] = (a[i] + g[i])/2 and g[i+1] = sqrt(a[i]*g[i]).//
//     The sequences satisfy the inequalities g[0] < g[1] < ... < a[1] < a[0].//
//     Further, lim g[n] = lim a[n].  The sequence of amplitudes is defined by//
//     setting phi[0] = phi, and for i >= 0,                                  //
//       tan(2^(i+1) phi[i+1] - 2^i phi[i]) = (g[i]/a[i]) tan(2^i phi[i]).    //
//     The elliptic integral of the first kind F(phi,k) = lim phi[n]/g[n].    //
//     The elliptic integral of the second kind                               //
//           E(phi,k) = lim (phi[n]/4g[n]) (4 - 2k^2 - Sum(2^j(a[j]-g[j])^2)  //
//                             + (1/2) Sum((a[j] - g[j])sin(2^j+1 phi[j+1]))  //
//     where the limit is as n -> inf and both Sums extend from j = 0 to n.   //
//                                                                            //
//  Arguments:                                                                //
//     long double  phi                                                       //
//                The upper limit of the integral, 0 <= phi <= pi / 2.        //
//     long double  parameter                                                 //
//                The parameter, the modulus squared k^2,                     //
//                0 < parameter < 1.                                          //
//     double *F                                                              //
//                The incomplete elliptic integral of the first kind with     //
//                parameter 'parameter' between 0 and 1 and amplitude between //
//                0 and pi/2.                                                 //
//     double *K                                                              //
//                The complete elliptic integral of the first kind with       //
//                parameter 'parameter' between 0 and 1.                      //
//     double *E                                                              //
//                The incomplete elliptic integral of the second kind with    //
//                parameter 'parameter' between 0 and 1 and amplitude between //
//                0 and pi/2.                                                 //
//     double *Em                                                             //
//                The complete elliptic integral of the first second with     //
//                parameter 'parameter' between 0 and 1.                      //
//                                                                            //
//  Return Value:                                                             //
//     The value of the elliptic integral of the first kind for the given     //
//     amplitude and parameter is returned via F, the value of the complete   //
//     elliptic integral of the first kind for the given parameter is returned//
//     via K, the value of the elliptic integral of the second kind for the   //
//     given amplitude and parameter is returned via E, and the value of the  //
//     complete elliptic integral of the second kind for the given parameter  //
//     is returned via Em.                                                    //
//                                                                            //
//  Example:                                                                  //
//     long double phi, m;                                                    //
//     double F, Fk, E, Em;                                                   //
//                                                                            //
//     ( code to initialize phi and m )                                       //
//                                                                            //
//     Landen_Transform( phi, m, &F, &Fk, &E, &Em );                          //
////////////////////////////////////////////////////////////////////////////////

static void Landen_Transform( long double phi, long double parameter,
                                 double *F, double *K, double *E, double *Em) 
{ 
   long double two_n = 1.0L;
   long double a = 1.0L;
   long double g = sqrtl(1.0L - parameter);
   long double a_old;
   long double g_old;
   long double tan_2n_phi;
   long double sum = 2.0L * (2.0L - parameter);
   long double integral = 0.0L;

   while (1) {
      tan_2n_phi = tanl(two_n * phi );
      sum -= two_n * (a - g) * (a - g);
      two_n += two_n;
      phi -= atanl( (a - g) * tan_2n_phi / (a + g * tan_2n_phi * tan_2n_phi) )
                                                                       / two_n;
      integral += (a - g) * sinl(two_n * phi);
      g_old = g;
      a_old = a;
      a = 0.5L * (g_old + a_old);
      g = sqrtl(g_old * a_old);
      if ( fabsl(a_old - g_old) <= (a_old * LDBL_EPSILON) ) break;
   }
   *F = (double) (phi / g);
   *K = (double) (PI_2 / g);
   *E =  (double) (0.5 * integral + 0.25 * sum * phi / g);
   *Em = (double) ((PI_4 / a) * sum / 2.0L);
}
