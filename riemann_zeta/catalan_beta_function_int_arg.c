////////////////////////////////////////////////////////////////////////////////
// File: catalan_beta_function_int_arg.c                                      //
// Routine(s):                                                                //
//    Catalan_Beta_Function_int_arg                                           //
//    xCatalan_Beta_Function_int_arg                                          //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for powl()
#include <float.h>         // required for DBL_MAX

//                         Internally Defined Routines                        //

double      Catalan_Beta_Function_int_arg(int s);
long double xCatalan_Beta_Function_int_arg(int s);
                                                                                
static long double xCatalan_Beta_Function_neg_int_arg(int s);

//                         Externally Defined Routines                        //

extern long double xCatalan_Beta_Function_pos_int_arg(int s);
extern long double xEuler_Number(int n);

////////////////////////////////////////////////////////////////////////////////
// double Catalan_Beta_Function_int_arg(int s)                                //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan Beta function                      //
//        beta(s) = Sum (-1)^k (1/(2k+1)^s) where the sum extends over        //
//     k = 0,... and where s > 1 is an integer,                               //
//     beta(0) = 1 / 2 and beta(1) = pi / 4.                                  //
//     For s < 0, beta(s) = E(-s) / 2, where E(n) is the nth Euler number.    //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the beta function restricted to the integers.       //
//                                                                            //
//  Return Value:                                                             //
//     beta(s), if beta(s) > DBL_MAX, then DBL_MAX is returned,               //
//              if beta(s) < -DBL_MAX, then -DBL_MAX is returned.             //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double beta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the beta function. )             //
//                                                                            //
//     beta = Catalan_Beta_Function_int_arg(s);                               //
////////////////////////////////////////////////////////////////////////////////

double Catalan_Beta_Function_int_arg(int s)
{
   long double x = xCatalan_Beta_Function_int_arg(s);

   if ( x > DBL_MAX) return DBL_MAX;
   if ( x < -DBL_MAX) return -DBL_MAX;
   return (double) x;
}


////////////////////////////////////////////////////////////////////////////////
// long double xCatalan_Beta_Function_int_arg(int s)                          //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan Beta function                      //
//        beta(s) = Sum (-1)^k (1/(2k+1)^s) where the sum extends over        //
//     k = 0,... and where s > 1 is an integer,                               //
//     beta(0) = 1 / 2 and beta(1) = pi / 4.                                  //
//     For s < 0, beta(s) = E(-s) / 2, where E(n) is the nth Euler number.    //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the beta function restricted to the integers.       //
//                                                                            //
//  Return Value:                                                             //
//     beta(s), if beta(s) > LDBL_MAX, then LDBL_MAX is returned,             //
//              if beta(s) < -LDBL_MAX, then -LDBL_MAX is returned.           //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double beta;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the beta function. )             //
//                                                                            //
//     beta = xCatalan_Beta_Function_int_arg(s);                              //
////////////////////////////////////////////////////////////////////////////////

long double xCatalan_Beta_Function_int_arg(int s)
{
   if (s >= 0) return xCatalan_Beta_Function_pos_int_arg(s);
   else        return xCatalan_Beta_Function_neg_int_arg(s);
}


////////////////////////////////////////////////////////////////////////////////
// static long double xCatalan_Beta_Function_neg_int_arg(int s)               //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Beta function for integer s < 0.           //
//                             Beta(s) = E[n] / 2                             //
//     where n = -s.                                                          //
//     If the argument is an even integer less than -1866, +-LDBL_MAX is      //
//     returned.                                                              //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the beta function restricted to the negative        //
//        integers.                                                           //
//                                                                            //
//  Return Value:                                                             //
//     beta(s), if beta(s) > LDBL_MAX, then LDBL_MAX is returned,             //
//              if beta(s) < -LDBL_MAX, then -LDBL_MAX is returned.           //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double beta;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the beta function. )             //
//                                                                            //
//     beta = xCatalan_Beta_Function_neg_int_arg(s);                          //
////////////////////////////////////////////////////////////////////////////////

static long double xCatalan_Beta_Function_neg_int_arg(int s)
{
   long double en = xEuler_Number(-s);
   
   if (fabsl(en) == LDBL_MAX) return en;

   return en / 2.0L;
}
