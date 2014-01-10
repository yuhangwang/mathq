////////////////////////////////////////////////////////////////////////////////
// File: jacobi_elliptic_functions.c                                          //
// Routine(s):                                                                //
//    Jacobi_sn_cn_dn                                                         //
//    Jacobi_cs_ds_ns                                                         //
//    Jacobi_sc_dc_nc                                                         //
//    Jacobi_sd_cd_nd                                                         //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for sqrt(), sqrtl(), fabs(), fabsl(),
                            // asinl(), sin(), sinl(), cos(), atan(), cosh(), 
                            // tanh() and M_PI_2
#include <float.h>          // required for LDBL_EPSILON

//                         Externally Defined Routines                        //

extern double Jacobi_am(double u, char arg,  double x);

//                         Internally Defined Routines                        //

void Jacobi_sn_cn_dn(double u, char arg, double x, double* sn, double* cn,
                                                                   double* dn);
void Jacobi_cs_ds_ns(double u, char arg, double x, double* cs, double* ds,
                                                                   double* ns);
void Jacobi_sc_dc_nc(double u, char arg, double x, double* sc, double* dc,
                                                                   double* nc);
void Jacobi_sd_cd_nd(double u, char arg, double x, double* sd, double* cd,
                                                                   double* nd);

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sn_cn_dn(double u, char arg, double x, double* sn, double* cn, //
//                                                                double* dn) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sn, cn, and dn //
//     which are defined in terms of the inverse to Legendre's elliptic       //
//     integral of the first kind.  If one sets u to the integral from        //
//     0 to phi, |phi| < pi / 2, of the integrand                             //
//                   dtheta / sqrt( 1 - k^2 sin^2(theta) )                    //
//     then         sn(u,k) = sin(phi),                                       //
//                  cn(u,k) = cos(phi),                                       //
//     and          dn(u,phi) = sqrt( 1 - k^2 sin^2(phi) )                    //
//                                                                            //
//      If the modulus k > 1, then the Jacobi real modulus transforms are     //
//      applied:                                                              //
//                          sn(x,k) = sn(kx,1/k)/k                            //
//                          cn(x,k) = dn(kx,1/k)                              //
//                          dn(x,k) = cn(kx,1/k).                             //
//                                                                            //
//      If the modulus k is purely imaginary, m < 0, then the transforms are  //
//      applied:                                                              //
//                          sn(x,ik) = sn(k'x,k/k')/k'                        //
//                          cn(x,ik) = dn(k'x,k/k')                           //
//                          dn(x,ik) = cn(k'x,k/k').                          //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of am():        //
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
//     double* sn                                                             //
//                The address of the value the Jacobi elliptic function sn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* cn                                                             //
//                The address of the value the Jacobi elliptic function cn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* dn                                                             //
//                The address of the value the Jacobi elliptic function dn    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sn, cn, and dn are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sn, cn, dn;                                                     //
//                                                                            //
//     ( code to initialize u and a )                                         //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sn_cn_dn( u, 'a', a, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sn,cn,dn,u,a);   //
//     Jacobi_sn_cn_dn( u, 'k', k, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sn,cn,dn,u,k);   //
//     Jacobi_sn_cn_dn( u, 'm', m, &sn, &cn, &dn);                            //
//     printf("sn cn dn %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sn,cn,dn,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sn_cn_dn(double u, char arg, double x, double* sn, double* cn,
                                                                    double* dn)
{
   double k;
   double m;
   double phi;
   int n;

   switch (arg) {
      case 'a': k = sin(x);
                m = k * k;
                break;
      case 'm': m = x;
                k = sqrt(fabs(m));
                break;
      default:  k = fabs(x);
                m = k * k;
   }

                   // Check special cases m = 1 or m = 0 //

   if ( m == 1.0 ) {
      *sn = tanh(u);
      *cn = 1.0 / cosh(u);
      *dn = *cn;
      return;
   }
   if ( m == 0.0 ) {
      *sn = sin(u);
      *cn = cos(u);
      *dn = 1.0;
      return;
   }

                  // If m > 1, perform the transformation: //
                  //        sn(u,k) = sn(k*u, 1/k) /k,     //
                  //        cn(u,k) = dn(k*u, 1/k),        //
                  //        dn(u,k) = cn(k*u, 1/k).        //

   if (m > 1.0) {
      Jacobi_sn_cn_dn(k * u, 'k', 1.0 / k, sn, dn, cn);
      *sn /= k;
      return;
   }

                  // If m < 0, perform the transformation: //
//sn(u,m) = [sn(sqrt(1-m)*u,-m/(1-m)) / dn(sqrt(1-m)*u,-m/(1-m))] / sqrt(1-m) //
//cn(u,m) = [cn(sqrt(1-m)*u,-m/(1-m)) / dn(sqrt(1-m)*u,-m/(1-m))]             //
//dn(u,m) = [1.0 / dn(sqrt(1-m)*u,-m/(1-m))]                                  //

   if (m < 0.0) {
      Jacobi_sd_cd_nd(sqrt(1.0 - m) * u, 'm', -m / (1.0 - m), sn, cn, dn);
      *sn /= sqrt(1.0 - m);
      return;
   }
 
   phi = Jacobi_am(u, arg, x);
   *sn = sin(phi);
   *cn = cos(phi);
   *dn = sqrt(1.0 - m * (*sn * *sn));

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_cs_ds_ns(double u, char arg, double x, double* cs, double* ds,
//                                                                double* ns) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: cs, ds, and    //
//     ns which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        cs(u,k) = cn(u,k) / sn(u,k)                         //
//                        ds(u,k) = dn(u,k) / sn(u,k)                         //
//     and                ns(u,k) = 1 / sn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the jacobi functions cs(u,x),        //
//                ds(u,x), or ns(u,x) corresponding to the second argument    //
//                of the elliptic integral of the first kind F(phi,x).        //
//                'x' may the the modulus, modular angle, or parameter        //
//                depending on the value of 'arg'.                            //
//     double* cs                                                             //
//                The address of the value the Jacobi elliptic function cs    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* ds                                                             //
//                The address of the value the Jacobi elliptic function ds    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* ns                                                             //
//                The address of the value the Jacobi elliptic function ns    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of cs, ds, and ns are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double cs, ds, ns;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_cs_ds_ns( u, 'a', a, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",cs,ds,ns,u,a);   //
//     Jacobi_cs_ds_ns( u, 'k', k, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",cs,ds,ns,u,k);   //
//     Jacobi_cs_ds_ns( u, 'm', m, &cs, &ds, &ns);                            //
//     printf("cs ds ns %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",cs,ds,ns,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_cs_ds_ns(double u, char arg, double x, double* cs, double* ds,
                                                                    double* ns)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (sn == 0.0) {
      *cs = DBL_MAX;
      *ds = DBL_MAX;
      *ns = DBL_MAX;
   } else {
      *cs = cn / sn;
      *ds = dn / sn;
      *ns = 1.0 / sn;
   }

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sc_dc_nc(double u, char arg, double x, double* sc, double* dc, //
//                                                                double* nc) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sc, dc, and    //
//     nc which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        sc(u,k) = sn(u,k) / cn(u,k)                         //
//                        dc(u,k) = dn(u,k) / cn(u,k)                         //
//     and                nc(u,k) = 1 / cn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi elliptic functions.              //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the amplitude function am(u,x)       //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.                                                   //
//     double* sc                                                             //
//                The address of the value the Jacobi elliptic function sc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* dc                                                             //
//                The address of the value the Jacobi elliptic function dc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* nc                                                             //
//                The address of the value the Jacobi elliptic function nc    //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sc, dc, and nc are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sc, dc, nc;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sc_dc_nc( u, 'a', a, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sc,dc,nc,u,a);   //
//     Jacobi_sc_dc_nc( u, 'k', k, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sc,dc,nc,u,k);   //
//     Jacobi_sc_dc_nc( u, 'm', m, &sc, &dc, &nc);                            //
//     printf("sc dc nc %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sc,dc,nc,u,m);   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sc_dc_nc(double u, char arg, double x, double* sc,
                                                      double* dc, double* nc)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (cn == 0.0) {
      *sc = DBL_MAX;
      *dc = DBL_MAX;
      *nc = DBL_MAX;
   } else {
      *sc = sn / cn;
      *dc = dn / cn;
      *nc = 1.0 / cn;
   }

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// void Jacobi_sd_cd_nd(double u, char arg, double x, double* sd, double* cd, //
//                                                                double* nd) //
//                                                                            //
//  Description:                                                              //
//     This function calculates the Jacobi elliptic functions: sd, cd, and    //
//     nd which are defined in terms of the Jacobi elliptic functions sn, cn  //
//     and dn (see above)                                                     //
//                        sd(u,k) = sn(u,k) / dn(u,k)                         //
//                        cd(u,k) = cn(u,k) / dn(u,k)                         //
//     and                nd(u,k) = 1 / dn(u,k).                              //
//                                                                            //
//  Arguments:                                                                //
//     double  u                                                              //
//                The argument of the Jacobi's elliptic functions.            //
//     char    arg                                                            //
//                The type of argument of the second argument of a Jacobi     //
//                elliptic function.                                          //
//                  If arg = 'k', then x = k, the modulus of F(phi,k).        //
//                  If arg = 'a', then x = alpha, the modular angle of        //
//                                F(phi \ alpha), alpha in radians.           //
//                  If arg = 'm', then x = m, the parameter of F(phi | m).    //
//                  The value of arg defaults to 'k'.                         //
//     double  x                                                              //
//                The second argument of the amplitude function am(u,x)       //
//                corresponding to the second argument of the elliptic        //
//                integral of the first kind F(phi,x).  'x' may the the       //
//                modulus, modular angle, or parameter depending on the value //
//                of 'arg'.                                                   //
//     double* sd                                                             //
//                The address of the value the Jacobi's elliptic function sd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* sc                                                             //
//                The address of the value the Jacobi's elliptic function cd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//     double* nd                                                             //
//                The address of the value the Jacobi's elliptic function nd  //
//                evaluated at u with modulus | modular angle | parameter x.  //
//                                                                            //
//  Return Value:                                                             //
//     The values of sd, cd, and nd are returned via the argument list.       //
//                                                                            //
//  Example:                                                                  //
//     double u, m, k, a;                                                     //
//     double sd, cd, nd;                                                     //
//                                                                            //
//     ( code to initialize u, and a )                                        //
//                                                                            //
//     k = sin(a);                                                            //
//     m = k * k;                                                             //
//     Jacobi_sd_cd_nd( u, 'a', a, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, angle(radians)       //
//                                                 %12.6f\n",sd,cd,nd,u,a);   //
//     Jacobi_sd_cd_nd( u, 'k', k, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, modulus              //
//                                                 %12.6f\n",sd,cd,nd,u,k);   //
//     Jacobi_sd_cd_nd( u, 'm', m, &sd, &cd, &nd);                            //
//     printf("sd cd nd %12,6f, %12.6f %12.6f, u %12.6f, parameter            //
//                                                 %12.6f\n",sd,cd,nd,u,m);   //
////////////////////////////////////////////////////////////////////////////////
                
void Jacobi_sd_cd_nd(double u, char arg, double x, double* sd,
                                                      double* cd, double* nd)
{
   double sn, cn, dn;

   Jacobi_sn_cn_dn( u, arg, x, &sn, &cn, &dn);

   if (dn == 0.0) {
      *sd = DBL_MAX;
      *cd = DBL_MAX;
      *nd = DBL_MAX;
   } else {
      *sd = sn / dn;
      *cd = cn / dn;
      *nd = 1.0 / dn;
   }

   return; 
}
