////////////////////////////////////////////////////////////////////////////////
// File: inverse_jacobi_elliptic_functions.c                                  //
// Routine(s):                                                                //
//    Inverse_Jacobi_sn                                                       //
//    Inverse_Jacobi_cn                                                       //
//    Inverse_Jacobi_dn                                                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     The incomplete elliptic integral of the first kind, F(phi,k), is the   //
//     integral from 0 to phi of the integrand                                //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) ).                   //
//     The upper limit phi is called the amplitude and the parameter k is     //
//     called the modulus.  This integral is even in k and odd in phi.        //
//     In practise the arguments of the elliptic function of the first kind   //
//     are also given as F(phi \ alpha) or F(phi | m) where the angle alpha,  //
//     called the modular angle, satisfies k = sin(alpha) and the argument    //
//     m = k^2 is simply called the parameter.                                //
//     In terms of these arguments F(phi \ alpha) = F(phi, sin(alpha))        //
//     and F(phi | m) = F(phi, sqrt(m)).                                      //
//     The complementary modulus is k' = sqrt(1 - k^2), similarly the         //
//     complementary parameter is m' = 1 - m.                                 //
//     The complete elliptic integral is K(k) = F(pi/2,k).                    //
//     For fixed modulus, modular angle, or parameter, the elliptic function  //
//     am(u,k):D -> R, called the elliptic amplitude function, assigns to a   //
//     number u in the domain D, which depends on k, the amplitude phi such   //
//     that u = F(phi,k), i.e. u = F(am(u,k),k).                              //
//     For 0 <= |k| < 1, the domain D is all of R.                            //
//                                                                            //
//     The Jacobi elliptic functions are defined in terms am(u,k) below with  //
//     domain specified for 0 < |k| < 1 and range a subset of R.              //
//            sn(u,k) = sin(am(u,k)) : R -> [-1,1]                            //
//            cn(u,k) = cos(am(u,k)) : R -> [-1,1]                            //
//            dn(u,k) = sqrt(1 - k^2 sin^2(am(u,k))) : R -> [k',1]            //
//            cs(u,k) = cn(u,k)/sn(u,k) = cot(am(u,k)) : R\{2n*K} -> R        //
//            ds(u,k) = dn(u,k)/sn(u,k) : R\{2n*K} -> R                       //
//            ns(u,k) = 1/sn(u,k) = csc(am(u,k)) : R\(2n*K} -> R              //
//            sc(u,k) = sn(u,k)/cn(u,k) = tan(am(u,k)) : R\{(2n+1)*K} -> R    //
//            dc(u,k) = dn(u,k)/cn(u,k) : R\{(2n+1)*K} -> R                   //
//            nc(u,k) = 1 / cn(u,k) = sec(am(u,k)) : R\{2n+1)*K} -> R         //
//            sd(u,k) = sn(u,k)/dn(u,k) : R -> [-1/k',1/k']                   //
//            cd(u,k) = cn(u,k)/dn(u,k) : R -> [-1/k',1/k']                   //
//            nd(u,k) = 1/dn(u,k) : R -> [1,1/k'].                            //
//     All the Jacobi elliptic functions are periodic with period either 2K,  //
//     (sc, cs, dn, and nd) or 4K (sn, sd, cn, cd, ns, ds, nc, dc).           //
//                                                                            //
//     If the modulus k = 0, then F(phi,k) = phi, k' = 1, K = pi/2 and        //
//     am(u,0) = u.                                                           //
//            sn(u,0) = sin(u) : R -> [-1,1]                                  //
//            cn(u,0) = cos(u) : R -> [-1,1]                                  //
//            dn(u,0) = 1 : R -> [1,1]                                        //
//            cs(u,0) = cot(u) : R\{n*pi} -> R                                //
//            ds(u,0) = csc(u : R\{n*pi} -> R                                 //
//            ns(u,0) = csc(u) : R\(n*pi} -> R                                //
//            sc(u,0) = tan(u) : R\{(n+1/2)*pi} -> R                          //
//            dc(u,0) = sec(u) : R\{(n+1/2)*pi} -> R                          //
//            nc(u,0) = sec(u) : R\{n+1/2)*pi} -> R                           //
//            sd(u,0) = sin(u) : R -> [-1,1]                                  //
//            cd(u,0) = cos(u) : R -> [-1,1]                                  //
//            nd(u,0) = 1 : R -> [1,1].                                       //
//                                                                            //
//     If the modulus k = 1, then F(phi,k) = ln tan(phi/2 + pi/4), k' = 0,    //
//     K = inf, and am(u,1) = 2*atan(exp(u)) - pi/2.                          //
//            sn(u,1) = tanh(u) : R -> (-1,1)                                 //
//            cn(u,1) = sech(u) : R -> (0,1]                                  //
//            dn(u,1) = sech(u) : R -> (0,1]                                  //
//            cs(u,1) = csch(u) : R\{0} -> R                                  //
//            ds(u,1) = csch(u) : R\{0} -> R                                  //
//            ns(u,1) = coth(u) : R\{0} -> R\[-1,1]                           //
//            sc(u,1) = sinh(u) : R -> R                                      //
//            dc(u,1) = 1       : R -> {1}                                    //
//            nc(u,1) = cosh(u) : R -> [1,inf)                                //
//            sd(u,1) = sinh(u) : R -> R                                      //
//            cd(u,1) = 1       : R -> {1}                                    //
//            nd(u,1) = cosh(u) : R -> [1,inf)                                //
//                                                                            //
//     If the modulus k > 1, then use Jacobi's real transforms:               //
//            sn(u,k) = sn(ku,1/k) / k                                        //
//            cn(u,k) = dn(ku,1/k)                                            //
//            dn(u,k) = cn(ku,1/k).                                           //
//     The remaining elliptic functions are thus:                             //
//            cs(u,k) = k ds(ku,1/k)                                          //
//            ds(u,k) = k cs(ku,1/k)                                          //
//            ns(u,k) = k ns(ku,1/k)                                          //
//            sc(u,k) = sd(ku,1/k) / k                                        //
//            dc(u,k) = cd(ku,1/k)                                            //
//            nc(u,k) = nd(ku,1/k)                                            //
//            sd(u,k) = sc(ku,1/k) / k                                        //
//            cd(u,k) = dc(ku,1/k)                                            //
//            nd(u,k) = nc(ku,1/k).                                           //
//                                                                            //
//     If the modulus is purely imaginary ik, then let k' = sqrt(1 + k^2) and //
//     use the transforms:                                                    //
//            sn(u,ik) = sd(k'u,k/k') / k'                                    //
//            cn(u,ik) = cd(k'u,k/k')                                         //
//            dn(u,ik) = nd(k'u,k/k').                                        //
//     The remaining elliptic functions are thus:                             //
//            cs(u,ik) = k' cs(k'u,k/k')                                      //
//            ds(u,ik) = k' ns(k'u,k/k')                                      //
//            ns(u,ik) = k' ds(k'u,k/k')                                      //
//            sc(u,ik) = sc(k'u,k/k') / k'                                    //
//            dc(u,ik) = nc(k'u,k/k')                                         //
//            nc(u,ik) = dc(k'u,k/k')                                         //
//            sd(u,ik) = sn(k'u,k/k') / k'                                    //
//            cd(u,ik) = cn(k'u,k/k')                                         //
//            nd(u,ik) = dn(k'u,k/k').                                        //
//                                                                            //
//     For real moduli |k| != 1, the inverse elliptic functions are           //
//     multivalued.  The routines below return the principal values of the    //
//     inverse elliptic functions.                                            //
//                                                                            //
//     The inverse Jacobi elliptic functions are given below with domain      //
//     specified for 0 < |k| < 1 and range of the principal value, where      //
//     K is the complete elliptic integral with modulus k and k' is the       //
//     complementary modulus k' = sqrt(1 - k^2).                              //
//      sn^(-1)(x,k) = F(asin(x),k) : [-1,1] -> [-K,K]                        //
//      cn^(-1)(x,k) = F(acos(x),k) : [-1,1] -> [0,2K]                        //
//      dn^(-1)(x,k) = F(asin(sqrt(1 - x*x)/|k|),k) : [k',1] -> [0,K]         //
//                                                                            //
//     The inverse Jacobi elliptic functions are given below with domain      //
//     specified for k = 0 and range of the principal value, where K = pi/2   //
//     is the complete elliptic integral and k' = 1 is the complementary      //
//     modulus.                                                               //
//      sn^(-1)(x,0) = asin(x) : [-1,1] -> [-pi/2,pi/2]                       //
//      cn^(-1)(x,0) = acos(x) : [-1,1] -> [0,pi]                             //
//      dn^(-1)(x,0) = 0 : [1,1] -> [0,0]                                     //
//                                                                            //
//     The inverse Jacobi elliptic functions are given below with domain      //
//     specified for k = 1 and range of the principal value, where K = inf    //
//     is the complete elliptic integral and k' = 0 is the complementary      //
//     modulus.                                                               //
//      sn^(-1)(x,1) = atanh(x) : (-1,1) -> R                                 //
//      cn^(-1)(x,1) = asech(x) : (0,1] -> [0,inf)                            //
//      dn^(-1)(x,1) = asech(x) : (0,1] -> [0,inf)                            //
//                                                                            //
//     The inverse Jacobi elliptic functions are given below with domain      //
//     specified for k > 1 and range of the principal value, where K = K(1/k).//
//      sn^(-1)(x,k) = sn^(-1)(kx,1/k)/k : [-1/k,1/k] -> [-K,K]               //
//      cn^(-1)(x,k) = dn^(-1)(x,1/k)/k : [sqrt(1-(1/k)^2),1] -> [0,K/k]      //
//      dn^(-1)(x,k) = cn^(-1)(x,1/k)/k : [-1,1] -> [0,2K/k]                  //
//                                                                            //
//     The inverse Jacobi elliptic functions are given below with domain      //
//     specified for a purely imaginary modulus and range of the principal    //
//     value, where k' = sqrt(1 + k^2) and  K = K(k/k').                      //
//      sn^(-1)(x,ik) = sd^(-1)(k'x,k/k')/k' : [-1,1] -> [-K/k',K/k']         //
//      cn^(-1)(x,ik) = cd^(-1)(x,k/k')/k' : [-1,1]] -> [0,2K/k']             //
//      dn^(-1)(x,ik) = nd^(-1)(x,k/k')/k' : [1,k'] -> [0,K/k']               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for sqrt(), fabs(), sin(), asin() 
                            // acos(), asinh(), atanh(), and acosh()  
#include <float.h>          // required for DBL_MAX

//                         Externally Defined Routines                        //

extern double Jacobi_am(double x, char arg,  double param);
extern double Legendre_Elliptic_Integral_First_Kind(double amplitude,
                                                     char arg, double param);

//                         Internally Defined Routines                        //

double Inverse_Jacobi_sn(double x, char arg, double param);
double Inverse_Jacobi_cn(double x, char arg, double param);
double Inverse_Jacobi_dn(double x, char arg, double param);
static double Inverse_Jacobi_sd(double x, double k);
static double Inverse_Jacobi_cd(double x, double k);
static double Inverse_Jacobi_nd(double x, double k);
static void parameter(char arg, double param, double *k, double *m);
static double asech(double x);

////////////////////////////////////////////////////////////////////////////////
// double Inverse_Jacobi_sn(double x, char arg, double param)                 //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     sn function for fixed modulus k, modular angle alpha, or parameter m   //
//     at x, i.e. this routine returns u where x = sn(u,k) = sn(u \ alpha)    //
//     = sn(u |m).                                                            //
//     Note that:                                                             //
//                    If m < 1 then -1 <= x <= 1.                             //
//                    If m = 1 then -1 < x < 1.                               //
//                    If m > 1 then -1/k <= x <= 1/k.                         //
//                                                                            //
//     If 0 < |k| < 1, then sn(-1)(x,k) = F(asin(x),k).                       //
//     sn^(-1)(*,k) : [-1,1] -> [-K,K]                                        //
//                                                                            //
//     If k = 0, then sn(-1)(x,0) = asin(x).                                  //
//     sn^(-1)(*,0) : [-1,1] -> [-pi/2,pi/2]                                  //
//                                                                            //
//     If |k| = 1, then sn(-1)(x,1) = atanh(x).                               //
//     sn^(-1)(*,1) : (-1,1) -> R                                             //
//     (If |x| = 1, then sgn(x) * DBL_MAX is returned.)                       //
//                                                                            //
//     If |k| > 1, then sn^(-1)(x,k) = sn^(-1)(kx,1/k)/k.                     //
//     sn^(-1)(*,k) : [-1/|k|,1/|k|] -> [-K(1/k)/|k|,K(1/k)/|k|]              //
//                                                                            //
//     If the modulus is purely imaginary, then                               //
//                    sn^(-1)(x,ik) = sd^(-1)(k'x,k/k')/k'.                   //
//     sn^(-1)(*,ik) : [-1,1] -> [-K(k/k')/k',K(k/k')/k']                     //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of sn(u,k) for some unknown u fixed modulus k,    //
//                modular angle alpha, or parameter m.                        //
//                If m < 1 then -1 <= x <= 1.                                 //
//                If m = 1 then -1 < x < 1.                                   //
//                   If |x| = 1, then sgn(x) DBL_MAX is returned.             //
//                If m > 1 then -1/|k| <= x <= 1/|k|.                         //
//     char    arg                                                            //
//                The type of argument of the second argument of sn():        //
//                  If arg = 'k', then param = k, the modulus of sn(phi,k).   //
//                  If arg = 'a', then param = alpha, the modular angle of    //
//                                sn(phi \ alpha), alpha in radians.          //
//                  If arg = 'm', then param = m, the parameter of            //
//                                sn(phi | m).                                //
//                  The value of arg defaults to 'k'.                         //
//     double  param                                                          //
//                The second argument of the jacobi function sn(u,param).     //
//                'param' may the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//                                                                            //
//  Return Value:                                                             //
//     The principal value of sn^(-1)(x,k) or sn^(-1)(x \ alpha) or           //
//      sn^(-1)(x | m).                                                       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and a )                                         //
//                                                                            //
//     u = Inverse_Jacobi_sn( x, 'a', a);                                     //
//     printf("u %12,6f\n",u);                                                //
//     k = sin(a);                                                            //
//     u = Inverse_Jacobi_sn( x, 'k', k);                                     //
//     printf("u %12,6f\n",u);                                                //
//     u = Inverse_Jacobi_sn( x, 'm', k*k);                                   //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

double Inverse_Jacobi_sn(double x, char arg, double param)
{
   double k;
   double m;
   double kp;
   double phi;

   parameter(arg, param, &k, &m);

                   // Check most common case: 0 < m < 1. //

   if ( m > 0.0 && m < 1.0) {
      phi = asin(x);
      return Legendre_Elliptic_Integral_First_Kind(phi, 'k', k);
   }

                   // Check special cases m = 0 or m = 1 //

   if ( m == 0.0 ) return asin(x);
   if ( m == 1.0 ) { 
      if ( x == -1.0) return -DBL_MAX;
      else {
         if (x == 1.0) return DBL_MAX;
         else return atanh(x);
      }
   }

                  // If m > 1, perform the transformation: //
                  // sn^(-1)(x,k) = sn^(-1)(k*x, 1/k) / k. //

   if (m > 1.0) return Inverse_Jacobi_sn(k * x, 'k', 1.0 / k) / k;

             // Remaining case: m < 0, i.e. k purely imaginary. //
             // Let kp = sqrt(1 + k * k) = sqrt(1 - m), then    //
             // sn^(-1)(x|m) = sd^(-1)(kp*x,|k|/kp) / kp.       //
   
   kp = sqrt(1.0 - m);
   return  Inverse_Jacobi_sd(kp * x, k / kp) / kp;

}

////////////////////////////////////////////////////////////////////////////////
// double Inverse_Jacobi_cn(double x, char arg, double param)                 //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     cn function for fixed modulus k, modular angle alpha, or parameter m   //
//     at x, i.e. this routine returns u where x = cn(u,k) = cn(u \ alpha)    //
//     = cn(u |m).                                                            //
//     Note that:                                                             //
//                    If m < 1 then -1 <= x <= 1.                             //
//                    If m = 1 then 0 < x <= 1.                               //
//                    If m > 1 then sqrt(1-1/k^2) <= x <= 1.                  //
//                                                                            //
//     If 0 < |k| < 1, then cn(-1)(x,k) = F(acos(x),k).                       //
//     cn^(-1)(*,k) : [-1,1] -> [0,2K]                                        //
//                                                                            //
//     If k = 0, then cn(-1)(x,0) = acos(x).                                  //
//     cn^(-1)(*,0) : [-1,1] -> [0,pi]                                        //
//                                                                            //
//     If |k| = 1, then cn(-1)(x,1) = asech(x).                               //
//     cn^(-1)(*,1) : (0,1] -> [0,inf)                                        //
//     (if |x| = 0, then DBL_MAX is returned.)                                //
//                                                                            //
//     If |k| > 1, then cn^(-1)(x,k) = dn^(-1)(x,1/k)/k.                      //
//     cn^(-1)(*,k) : [sqrt(1-1/k^2),1] -> [0,K(1/k)/k]                       //
//                                                                            //
//     If the modulus is purely imaginary, then                               //
//                    cn^(-1)(x,ik) = cd^(-1)(x,k/k')/k'.                     //
//     cn^(-1)(*,k) : [-1,1] -> [0,2K(k/k')/k']                               //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of cn(u,k) for some unknown u fixed modulus k,    //
//                modular angle alpha, or parameter m.                        //
//                If m < 1 then -1 <= x <= 1.                                 //
//                If m = 1 then 0 < x <= 1.                                   //
//                   If |x| = 0, then DBL_MAX is returned.                    //
//                If m > 1 then sqrt(1-1/k^2) <= x <= 1.                      //
//     char    arg                                                            //
//                The type of argument of the second argument of cn():        //
//                  If arg = 'k', then param = k, the modulus of cn(phi,k).   //
//                  If arg = 'a', then param = alpha, the modular angle of    //
//                                cn(phi \ alpha), alpha in radians.          //
//                  If arg = 'm', then param = m, the parameter of            //
//                                cn(phi | m).                                //
//                  The value of arg defaults to 'k'.                         //
//     double  param                                                          //
//                The second argument of the jacobi function cn(u,param).     //
//                'param' may the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//                                                                            //
//  Return Value:                                                             //
//     The value of cn^(-1)(x,k) or cn^(-1)(x \ alpha) or cn^(-1)(x | m).     //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and a )                                         //
//                                                                            //
//     u = Inverse_Jacobi_cn( x, 'a', a);                                     //
//     printf("u %12,6f\n",u);                                                //
//     k = sin(a);                                                            //
//     u = Inverse_Jacobi_cn( x, 'k', k);                                     //
//     printf("u %12,6f\n",u);                                                //
//     u = Inverse_Jacobi_cn( x, 'm', k*k);                                   //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

double Inverse_Jacobi_cn(double x, char arg, double param)
{
   double k;
   double m;
   double kp;
   double phi;

   parameter(arg, param, &k, &m);

                   // Check most common case: 0 < m < 1. //

   if ( m > 0.0 && m < 1.0) {
      phi = acos(x);
      return Legendre_Elliptic_Integral_First_Kind(phi, 'k', k);
   }

                   // Check special cases m = 1 or m = 0 //

   if ( m == 0.0 ) return acos(x);
   if ( m == 1.0 ) return asech(x);

                  // If m > 1, perform the transformation: //
                  // cn^(-1)(x,k) = dn^(-1)(x, 1/k) / k. //

   if (m > 1.0) return Inverse_Jacobi_dn(x, 'k', 1.0 / k) / k;

             // Remaining case: m < 0, i.e. k purely imaginary. //
            // Let kp = sqrt(1 + k * k) = sqrt(1 - m), then     //
             // cn^(-1)(x,m) = cd^(-1)(x,|k|/kp) / kp.          //
   
   kp = sqrt(1.0 - m);
   return  Inverse_Jacobi_cd(x, k / kp) / kp;

}

////////////////////////////////////////////////////////////////////////////////
// double Inverse_Jacobi_dn(double x, char arg, double param)                 //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     dn function for fixed modulus k, modular angle alpha, or parameter m   //
//     at x, i.e. this routine returns u where x = dn(u,k) = dn(u \ alpha)    //
//     = dn(u |m).                                                            //
//     Note that:                                                             //
//                    If m < 0 then 1 <= x <= k'.                             //
//                    If m = 0 then x = 1. (0 is returned regardless of x.)   //
//                    If 0 < m < 1 then k'<= x <= 1.                          //
//                    If m = 1 then 0 < x <= 1.                               //
//                    If m > 1 then -1 <= x <= 1.                             //
//                                                                            //
//     If 0 < |k| < 1, then dn(-1)(x,k) = F(asin(sqrt(1 - x*x)/|k|),k).       //
//     dn^(-1)(*,k) : [k',1] -> [0,K]                                         //
//                                                                            //
//     If k = 0, then by convention dn(-1)(x,0) = 0, dn(u,0) = 1 identically. //
//     dn^(-1)(*,0) : [1,1] -> [0,0]                                          //
//                                                                            //
//     If |k| = 1, then dn(-1)(x,1) = asech(x).                               //
//     dn^(-1)(*,1) : (0,1] -> [0,inf)                                        //
//     (if |x| = 0, then DBL_MAX is returned.)                                //
//                                                                            //
//     If |k| > 1, then dn^(-1)(x,k) = cn^(-1)(x,1/k)/k.                      //
//     dn^(-1)(*,k) : [-1,1] -> [0,2K(1/k)/k]                                 //
//                                                                            //
//     If the modulus is purely imaginary, then                               //
//                    dn^(-1)(x,ik) = nd^(-1)(x,k/k')/k'.                     //
//     dn^(-1)(*,k) : [1,k'] -> [0,K(k/k')/k']                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of dn(u,k) for some unknown u fixed modulus k,    //
//                modular angle alpha, or parameter m.                        //
//                If m < 0 then 1 <= x <= k'.                                 //
//                If m = 0 then x = 1. (0 is returned regardless of x.)       //
//                If 0 < m < 1 then k' <= x <= 1.                             //
//                If m = 1 then 0 < x <= 1.                                   //
//                   If |x| = 0, then DBL_MAX is returned.                    //
//                If m > 1 then -1 <= x <= 1.                                 //
//     char    arg                                                            //
//                The type of argument of the second argument of dn():        //
//                  If arg = 'k', then param = k, the modulus of dn(phi,k).   //
//                  If arg = 'a', then param = alpha, the modular angle of    //
//                                dn(phi \ alpha), alpha in radians.          //
//                  If arg = 'm', then param = m, the parameter of            //
//                                dn(phi | m).                                //
//                  The value of arg defaults to 'k'.                         //
//     double  param                                                          //
//                The second argument of the jacobi function dn(u,param).     //
//                'param' may the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//                                                                            //
//  Return Value:                                                             //
//     The value of dn^(-1)(x,k) or dn^(-1)(x \ alpha) or dn^(-1)(x | m).     //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and a )                                         //
//                                                                            //
//     u = Inverse_Jacobi_dn( x, 'a', a);                                     //
//     printf("u %12,6f\n",u);                                                //
//     k = sin(a);                                                            //
//     u = Inverse_Jacobi_dn( x, 'k', k);                                     //
//     printf("u %12,6f\n",u);                                                //
//     u = Inverse_Jacobi_dn( x, 'm', k*k);                                   //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

double Inverse_Jacobi_dn(double x, char arg, double param)
{
   double k;
   double m;
   double kp;
   double phi;

   parameter(arg, param, &k, &m);

                   // Check most common case: 0 < m < 1. //

   if ( m > 0.0 && m < 1.0) {
      phi = asin(sqrt(1.0 - x*x)/fabs(k));
      return Legendre_Elliptic_Integral_First_Kind(phi, 'k', k);
   }

                   // Check special cases m = 1 or m = 0 //

   if ( m == 0.0 ) return 0.0;
   if ( m == 1.0 ) return asech(x);

                  // If m > 1, perform the transformation: //
                  // dn^(-1)(x,k) = cn^(-1)(x, 1/k) / k. //

   if (m > 1.0) return Inverse_Jacobi_cn(x, 'k', 1.0 / k) / k;

             // Remaining case: m < 0, i.e. k purely imaginary. //
             // Let kp = sqrt(1 + k * k) = sqrt(1 - m), then    //
             // dn^(-1)(x,m) = nd^(-1)(x,|k|/kp) / kp.          //
   
   kp = sqrt(1.0 - m);
   return  Inverse_Jacobi_nd(x, k / kp) / kp;

}


////////////////////////////////////////////////////////////////////////////////
// static double Inverse_Jacobi_sd(double x, double k)                        //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     sd function for fixed modulus k for 0 <= |k| <= 1 at x, i.e. this      //
//     routine returns u where x = sd(u,k).                                   //
//                                                                            //
//     If 0 < |k| < 1, then sd(-1)(x,k) = F(asin(x/sqrt(1+(kx)^2)),k).        //
//     sd^(-1)(*,k) : [-1/k',1/k'] -> [-K,K]                                  //
//                                                                            //
//     If k = 0, then sd(-1)(x,0) = asin(x).                                  //
//     sd^(-1)(*,0) : [-1,1] -> [-pi/2,pi/2]                                  //
//                                                                            //
//     If |k| = 1, then sd(-1)(x,1) = asinh(x).                               //
//     sd^(-1)(*,1) : R -> R                                                  //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of sd(u,k) for some unknown u fixed modulus k,    //
//                where 0 <= |k| <= 1 and  -k' <= x <= k'.                    //
//     double  k                                                              //
//                The modulus of sd(u,k) where 0 <= k <= 1.                   //
//                                                                            //
//  Return Value:                                                             //
//     The value of sd^(-1)(x,k).                                             //
//                                                                            //
//  Example:                                                                  //
//     double u, k;                                                           //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and k )                                         //
//                                                                            //
//     u = Inverse_Jacobi_sd( x, k);                                          //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

static double Inverse_Jacobi_sd(double x, double k)
{
   double m = k * k;
   double phi;

                   // Check most common case: 0 < m < 1. //

   if ( m > 0.0 && m < 1.0) {
      phi = asin( x / sqrt( 1.0 + m * x * x ) );
      return Legendre_Elliptic_Integral_First_Kind(phi, 'k', k);
   }

                   // Check special cases m = 1 or m = 0 //

   if ( m == 0.0 ) return asin(x);
   return asinh(x);
}


////////////////////////////////////////////////////////////////////////////////
// static double Inverse_Jacobi_cd(double x, double k)                        //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     cd function for fixed modulus k for 0 <= |k| <= 1 at x, i.e. this      //
//     routine returns u where x = cd(u,k).                                   //
//                                                                            //
//     If 0 < |k| < 1, then cd(-1)(x,k) = F(acos(x*sqrt(1-k*k)/(1-(kx)^2)),k).//
//     cd^(-1)(*,k) : [-1,1] -> [0,2K]                                        //
//                                                                            //
//     If k = 0, then cd(-1)(x,0) = acos(x).                                  //
//     cd^(-1)(*,0) : [-1,1] -> [0,pi]                                        //
//                                                                            //
//     If |k| = 1, then cd(-1)(x,1) = 0.                                      //
//     cd^(-1)(*,1) : [1,1] -> [0,0]                                          //
//     (cd(u,1) = 1 for all u.)                                               //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of cd(u,k) for some unknown u fixed modulus k,    //
//                where 0 <= |k| <= 1 and  -1 <= x <= 1.                      //
//     double  k                                                              //
//                The modulus of cd(u,k) where 0 <= k <= 1.                   //
//                                                                            //
//  Return Value:                                                             //
//     The value of cd^(-1)(x,k).                                             //
//                                                                            //
//  Example:                                                                  //
//     double u, k;                                                           //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and k )                                         //
//                                                                            //
//     u = Inverse_Jacobi_cd( x, k);                                          //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

static double Inverse_Jacobi_cd(double x, double k)
{
   double m = k * k;
   double phi;

                   // Check most common case: 0 < m < 1. //

   if ( m > 0.0 && m < 1.0) {
      phi = acos( x * sqrt( (1.0 - m)/(1.0 - m*x*x) ) );
      return Legendre_Elliptic_Integral_First_Kind(phi, 'k', k);
   }

                   // Check special cases m = 0 or m = 1 //

   if ( m == 0.0 ) return acos(x);
   return 0.0;

}

////////////////////////////////////////////////////////////////////////////////
// static double Inverse_Jacobi_nd(double x, double k)                        //
//                                                                            //
//  Description:                                                              //
//     This routine returns the principal value of the inverse of the Jacobi  //
//     nd function for fixed modulus k for 0 <= |k| <= 1 at x, i.e. this      //
//     routine returns u where x = nd(u,k).                                   //
//                                                                            //
//     nd^(-1)(x,k) = dn^(-1)(1/x,k).                                         //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The value of nd(u,k) for some unknown u fixed modulus k,    //
//                where 0 <= |k| <= 1 and  -1/k' <= x <= 1/k'.                //
//     double  k                                                              //
//                The modulus of nd(u,k) where 0 <= k <= 1.                   //
//                                                                            //
//  Return Value:                                                             //
//     The value of nd^(-1)(x,k).                                             //
//                                                                            //
//  Example:                                                                  //
//     double u, k;                                                           //
//     double x;                                                              //
//                                                                            //
//     ( code to initialize x and k )                                         //
//                                                                            //
//     u = Inverse_Jacobi_nd( x, 'k', k);                                     //
//     printf("u %12,6f\n",u);                                                //
////////////////////////////////////////////////////////////////////////////////

static double Inverse_Jacobi_nd(double x, double k)
{
   return Inverse_Jacobi_dn(1.0 / x, 'k', k);
}


////////////////////////////////////////////////////////////////////////////////
// static void parameter(char arg, double x, double *k, double *m)            //
//                                                                            //
//  Description:                                                              //
//     This function returns the modulus, k, and the parameter, m, given      //
//     the modular angle, alpha, modulus, k, or parameter m.                  //
//     If the modulus k is purely imaginary, m < 0.                           //
//                          dn(x,ik) = cn(k'x,k/k').                          //
//                                                                            //
//  Arguments:                                                                //
//     char    arg                                                            //
//                The type of argument of the second argument of F():         //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the jacobi functions sn(u,x),        //
//                cn(u,x), or dn(u,x) corresponding to the second argument    //
//                of the elliptic integral of the first kind F(phi,x).        //
//                'x' may the the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//     double* k                                                              //
//                The address of the value the Jacobi elliptic function sn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* m                                                              //
//                The address of the value the Jacobi elliptic function cn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of k and m are returned via the argument list.              //
//                                                                            //
//  Example:                                                                  //
//     char arg;                                                              //
//     double x, k, m;                                                        //
//                                                                            //
//     ( code to initialize arg and x)                                        //
//                                                                            //
//     parameter(arg, x, &k, &m);                                             //
////////////////////////////////////////////////////////////////////////////////

static void parameter(char arg, double x, double *k, double *m)
{
   switch (arg) {
      case 'a': *k = sin(x);
                *m = *k * *k;
                break;
      case 'm': *m = x;
                *k = sqrt(fabs(*m));
                break;
      default:  *k = fabs(x);
                *m = *k * *k;
   }
}

////////////////////////////////////////////////////////////////////////////////
// static double asech(double x)                                              //
//                                                                            //
//  Description:                                                              //
//     This function returns the value u >= 0 such that x = sech(u).          //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of asech(), if x = 0.0, then DBL_MAX is        //
//                returned.                                                   //
//                                                                            //
//  Return Value:                                                             //
//     The value u such that x = sech(u).                                     //
//                                                                            //
//  Example:                                                                  //
//     double x,u;                                                            //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     u = asech(x);                                                          //
////////////////////////////////////////////////////////////////////////////////

static double asech(double x)
{
   if (x == 0.0) return DBL_MAX;
   return acosh(1.0 / x);
}
