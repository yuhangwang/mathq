////////////////////////////////////////////////////////////////////////////////
// File: riemann_zeta_function_int_arg.c                                      //
// Routine(s):                                                                //
//    Riemann_Zeta_Function_int_arg                                           //
//    xRiemann_Zeta_Function_int_arg                                          //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for powl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double      Riemann_Zeta_Function_int_arg(int s);
long double xRiemann_Zeta_Function_int_arg(int s);

//                         Externally Defined Routines                        //

long double xRiemann_Zeta_Function_pos_int_arg(int s);
long double xBernoulli_Number(int n);

////////////////////////////////////////////////////////////////////////////////
// double Riemann_Zeta_Function_int_arg(int s)                                //
//                                                                            //
//  Description:                                                              //
//     The Riemann zeta function is defined as:                               //
//                          zeta(s) = Sum (1/k^s),                            //
//     summed over k = 1,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Riemann zeta function of an integer s for  //
//     s != 1.  If s = 1, DBL_MAX is returned.                                //
//     If zeta(s) > DBL_MAX, then DBL_MAX is returned and                     //
//     if zeta(s) < -DBL_MAX, then -DBL_MAX is returned. Note that            //
//     lim zeta(s) = -inf where s->1 from the left and lim zeta(s) = inf      //
//     where s->1 from the right.                                             //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Riemann zeta function restricted to the         //
//        integers, s != 1.  If s = 1, then DBL_MAX is returned.              //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if zeta(s) > DBL_MAX, then DBL_MAX is returned,               //
//              if zeta(s) < -DBL_MAX, then -DBL_MAX is returned,             //
//              if s = 1, then DBL_MAX is returned but note that s = 1 is     //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double zeta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta function. )     //
//                                                                            //
//     zeta = Riemann_Zeta_Function_int_arg(s);                               //
////////////////////////////////////////////////////////////////////////////////

double Riemann_Zeta_Function_int_arg(int s)
{
   long double x = xRiemann_Zeta_Function_int_arg(s);

   if ( x > DBL_MAX) return DBL_MAX;
   if ( x < -DBL_MAX) return -DBL_MAX;
   return (double) x;
}


////////////////////////////////////////////////////////////////////////////////
// long double xRiemann_Zeta_Function_int_arg(int s)                          //
//                                                                            //
//  Description:                                                              //
//     The Riemann zeta function is defined as:                               //
//                          zeta(s) = Sum (1/k^s),                            //
//     summed over k = 1,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Riemann zeta function of an integer s for  //
//     s != 1.  If s = 1, LDBL_MAX is returned.                               //
//     If zeta(s) > LDBL_MAX, then LDBL_MAX is returned and                   //
//     if zeta(s) < -LDBL_MAX, then -LDBL_MAX is returned. Note that          //
//     lim zeta(s) = -inf where s->1 from the left and lim zeta(s) = inf      //
//     where s->1 from the right.                                             //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Riemann zeta function restricted to the         //
//        integers, s != 1.  If s = 1, then LDBL_MAX is returned.             //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if zeta(s) > LDBL_MAX, then LDBL_MAX is returned,             //
//              if zeta(s) < -LDBL_MAX, then -LDBL_MAX is returned,           //
//              if s = 1, then LDBL_MAX is returned but note that s = 1 is    //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double zeta;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta function. )     //
//                                                                            //
//     zeta = xRiemann_Zeta_Function_int_arg(s);                              //
////////////////////////////////////////////////////////////////////////////////

long double xRiemann_Zeta_Function_int_arg(int s)
{
   long double bn;

   if (s >= 0) return xRiemann_Zeta_Function_pos_int_arg(s);

   bn = xBernoulli_Number(1-s);
   if (fabsl(bn) == LDBL_MAX) return -bn;
   return -bn / (long double)(1-s);
}
