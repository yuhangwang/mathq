////////////////////////////////////////////////////////////////////////////////
// File: riemann_zeta_function.c                                              //
// Routine(s):                                                                //
//    Riemann_Zeta_Function                                                   //
//    xRiemann_Zeta_Function                                                  //
//    Riemann_Zeta_Star_Function                                              //
//    xRiemann_Zeta_Star_Function                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                      // required for fabsl(), powl()
#include <float.h>                     // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double Riemann_Zeta_Function(double s);
long double xRiemann_Zeta_Function(long double s);
double Riemann_Zeta_Star_Function(double s);
long double xRiemann_Zeta_Star_Function(long double s);

//                         Externally Defined Routines                        //

extern long double xDirichlet_Eta_Function(long double s);
extern long double xDirichlet_Eta_Star_Function(long double s);

////////////////////////////////////////////////////////////////////////////////
// double Riemann_Zeta_Function( double s )                                   //
//                                                                            //
//  Description:                                                              //
//     The Riemann zeta function is defined as:                               //
//                          zeta(s) = Sum (1/k^s),                            //
//     summed over k = 1,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Riemann zeta function for real s != 1.     //
//     If s = 1, DBL_MAX is returned.                                         //
//     If zeta(s) > DBL_MAX, then DBL_MAX is returned and                     //
//     if zeta(s) < -DBL_MAX, then -DBL_MAX is returned. Note that            //
//     lim zeta(s) = -inf where s->1 from the left and lim zeta(s) = inf      //
//     where s->1 from the right.                                             //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument for the Riemann zeta function.                         //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if zeta(s) > DBL_MAX, then DBL_MAX is returned,               //
//              if zeta(s) < -DBL_MAX, then -DBL_MAX is returned,             //
//              if s = 1, then DBL_MAX is returned but note that s = 1 is     //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     double zeta, s;                                                        //
//     ...                                                                    //
//                         (* User code to set s *)                           //
//                                                                            //
//     zeta = Riemann_Zeta(s);                                                //
////////////////////////////////////////////////////////////////////////////////

double Riemann_Zeta_Function(double s)
{
   long double zeta;

   if (s == 1.0) return DBL_MAX;

   zeta = xRiemann_Zeta_Function((long double)s);

   if (fabsl(zeta) < DBL_MAX) return (double) zeta;
   if (zeta >= DBL_MAX) return DBL_MAX;
   return -DBL_MAX;
}

////////////////////////////////////////////////////////////////////////////////
// long double xRiemann_Zeta_Function( long double s )                        //
//                                                                            //
//  Description:                                                              //
//     The Riemann zeta function is defined as:                               //
//                          zeta(s) = Sum (1/k^s),                            //
//     summed over k = 1,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Riemann zeta function for real s != 1.     //
//     If s = 1, LDBL_MAX is returned.                                        //
//     lim zeta(s) = -inf where s->1 from the left and lim zeta(s) = inf      //
//     where s->1 from the right.                                             //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument for the Riemann zeta function.                         //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if s = 1, then LDBL_MAX is returned but note that s = 1 is    //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     long double zeta, s;                                                   //
//     ...                                                                    //
//                         (* User code to set s *)                           //
//                                                                            //
//     zeta = Riemann_Zeta(s);                                                //
////////////////////////////////////////////////////////////////////////////////
long double xRiemann_Zeta_Function(long double s)
{
   long double zeta;

   if (s == 1.0L) return LDBL_MAX;

   return xDirichlet_Eta_Function(s) / (1.0L - powl(2.0L,1.0L-s));
}


////////////////////////////////////////////////////////////////////////////////
// double Riemann_Zeta_Star_Function(double s)                                //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemann zeta star function defined as      //
//                           zeta*(s) = zeta(s) - 1,                          //
//     where zeta() is the Riemann zeta function for real s != 1.             //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the zeta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta*(s), if zeta*(s) > DBL_MAX, then DBL_MAX is returned,             //
//               if zeta*(s) < -DBL_MAX, then -DBL_MAX is returned,           //
//               if s = 1, then DBL_MAX is returned but note that s = 1 is    //
//               a simple pole of the zeta*() function.                       //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double zeta_star;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the zeta* function. )            //
//                                                                            //
//     zeta_star = Riemann_Zeta_Star_Function(s);                             //
////////////////////////////////////////////////////////////////////////////////

double Riemann_Zeta_Star_Function(double s)
{
   long double x;

   x = xRiemann_Zeta_Star_Function((long double)s);

   if (fabsl(x) < DBL_MAX) return (double) x;
   else return ( x > 0.0L ) ? DBL_MAX : - DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xRiemann_Zeta_Star_Function(long double s)                     //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemannt zeta star function defined as     //
//                           zeta*(s) = zeta(s) - 1,                          //
//     where zeta() is the Riemann zeta function for real s != 1.             //
//                                                                            //
//  Arguments:                                                                //
//     long double s;                                                         //
//        The argument of the zeta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta*(s), if s = 1, then LDBL_MAX is returned but note that s = 1 is   //
//               a simple pole of the zeta*() function.                       //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double zeta_star;                                                 //
//                                                                            //
//     ( User code to set s, the argument to the eta* function. )             //
//                                                                            //
//     zeta_star = xRiemann_Zeta_Star_Function(s);                            //
////////////////////////////////////////////////////////////////////////////////

long double xRiemann_Zeta_Star_Function(long double s)
{
   long double two_s = powl(2.0L,s - 1.0L);
   long double eta = xDirichlet_Eta_Star_Function(s);
   long double temp;

   if (s == 1.0L) return LDBL_MAX;

   if ( s < 0.0L ) {
      two_s = 1.0L / two_s;
      return (eta + two_s) / (1.0L - two_s);
   }
   
   temp = two_s - 1.0L;
   return (two_s * eta + 1.0L) / temp;
}
