////////////////////////////////////////////////////////////////////////////////
// File: dirichlet_eta_function_int_arg.c                                     //
// Routine(s):                                                                //
//    Dirichlet_Eta_Function_int_arg                                          //
//    xDirichlet_Eta_Function_int_arg                                         //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for powl()
#include <float.h>         // required for DBL_MAX

//                         Internally Defined Routines                        //

double      Dirichlet_Eta_Function_int_arg(int s);
long double xDirichlet_Eta_Function_int_arg(int s);
                                                                                
static long double xDirichlet_Eta_Function_neg_int_arg(int s);

//                         Externally Defined Routines                        //

long double xDirichlet_Eta_Function_pos_int_arg(int s);
long double xBernoulli_Number(int n);

////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Eta_Function_int_arg(int s)                               //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function                               //
//             eta(s) = (1 - 2^(1-s)) zeta(s), where zeta() is the Riemann    //
//     zeta function.  For s > 1,                                             //
//             eta(s) = Sum (-1)^(k-1) (1/k^s), where the sum is summed over  //
//     k = 1,... .  eta(0) = 1/2, eta(1) = ln(2).                             //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the eta function restricted to the integers.        //
//                                                                            //
//  Return Value:                                                             //
//     eta(s), if eta(s) > DBL_MAX, then DBL_MAX is returned,                 //
//             if eta(s) < -DBL_MAX, then -DBL_MAX is returned.               //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double eta;                                                            //
//                                                                            //
//     ( User code to set s, the argument to the eta function. )              //
//                                                                            //
//     eta = Dirichlet_Eta_Function_int_arg(s);                               //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Eta_Function_int_arg(int s)
{
   long double x = xDirichlet_Eta_Function_int_arg(s);

   if ( x > DBL_MAX) return DBL_MAX;
   if ( x < -DBL_MAX) return -DBL_MAX;
   return (double) x;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Eta_Function_int_arg(int s)                         //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function                               //
//        eta(s) = Sum (-1)^(k-1) (1/k^s) = (1 - 2^(1-s)) zeta(s), where the  //
//     sum is summed over k = 1,... and for s > 1 an integer,                 //
//     eta(1) = ln(2).                                                        //
//     For s an integer < 1, the eta function is defined in terms of the      //
//     Riemann zeta function as:                                              //
//                      eta(s) = (1 - 2^(1-s)) zeta(s).                       //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the eta function restricted to the integers.        //
//                                                                            //
//  Return Value:                                                             //
//     eta(s)                                                                 //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double eta;                                                       //
//                                                                            //
//     ( User code to set s, the argument to the eta function. )              //
//                                                                            //
//     eta = xDirichlet_Eta_Function_int_arg(s);                              //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Eta_Function_int_arg(int s)
{
   if (s >= 0) return xDirichlet_Eta_Function_pos_int_arg(s);
   else        return xDirichlet_Eta_Function_neg_int_arg(s);
}

////////////////////////////////////////////////////////////////////////////////
// static long double xDirichlet_Eta_Function_neg_int_arg(int s)              //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function                               //
//             eta(s) = (1 - 2^(1-s)) zeta(s), where zeta() is the Riemann    //
//     zeta function.  For s > 1,                                             //
//             eta(s) = Sum (-1)^(k-1) (1/k^s), where the sum is summed over  //
//     k = 1,... .  eta(0) = 1/2, eta(1) = ln(2).                             //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the eta function restricted to the negative         //
//        integers.                                                           //
//                                                                            //
//  Return Value:                                                             //
//     eta(s), if eta(s) > LDBL_MAX, then LDBL_MAX is returned,               //
//             if eta(s) < -LDBL_MAX, then -LDBL_MAX is returned.             //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double eta;                                                       //
//                                                                            //
//     ( User code to set s, the argument to the eta function. )              //
//                                                                            //
//     eta = xDirichlet_Eta_Function_neg_int_arg(s);                          //
////////////////////////////////////////////////////////////////////////////////

static long double xDirichlet_Eta_Function_neg_int_arg(int s)
{
   long double bn = xBernoulli_Number(1-s);
   long double two_n;
   long double n;
   
   if (fabsl(bn) == 0.0L) return bn;

   n = (long double) (1 - s);
   two_n = powl(2.0L, n) - 1.0L;
   if (bn > 0.0L) {
      if (bn == LDBL_MAX) return LDBL_MAX;
      if (bn >= (LDBL_MAX / two_n) * n) return LDBL_MAX;
      return two_n * bn / n;
   }
   if (bn == -LDBL_MAX) return -LDBL_MAX;
   if (bn <= -(LDBL_MAX / two_n) * n) return -LDBL_MAX;
   return two_n * bn / n;
}
