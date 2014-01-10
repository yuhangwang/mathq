////////////////////////////////////////////////////////////////////////////////
// File: dirichlet_lambda_function.c                                          //
// Routine(s):                                                                //
//    Dirichlet_Lambda_Function                                               //
//    xDirichlet_Lambda_Function                                              //
//    Dirichlet_Lambda_Star_Function                                          //
//    xDirichlet_Lambda_Star_Function                                         //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for fabsl(), powl() 
#include <float.h>         // required for DBL_MAX, LDBL_MAX, LDBL_MAX_EXP

//                         Internally Defined Routines                        //

double Dirichlet_Lambda_Function(double s);
long double xDirichlet_Lambda_Function(long double s);
double Dirichlet_Lambda_Star_Function(double s);
long double xDirichlet_Lambda_Star_Function(long double s);

static long double Lambda_Star_Positive_Arg(long double s);
static long double Lambda_Star_Negative_Arg(long double s);

//                         Externally Defined Routines                        //

long double xRiemann_Zeta_Star_Function(long double s);


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Lambda_Function(double s)                                 //
//                                                                            //
//  Description:                                                              //
//     The Dirichlet lambda function is defined as:                           //
//                       lambda(s) = Sum (1/(2k+1)^s),                        //
//     summed over k = 0,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Dirichlet Lambda function of a real s != 1,//
//     If s = 1, DBL_MAX is returned.                                         //
//     If lambda(s) > DBL_MAX, then DBL_MAX is returned and                   //
//     if lambda(s) < -DBL_MAX, then -DBL_MAX is returned. Note that          //
//     lim lambda(s) = -inf where s->1 from the left and lim lambda(s) = inf  //
//     where s->1 from the right.                                             //
//                                                                            //
//     Note that  lambda(s) = [ zeta(s) + eta(s) ] / 2 , where zeta is the    //
//     Riemann zeta function and eta is the Dirichlet eta function.           //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the Dirichlet lambda function restricted to the     //
//        reals.  If s = 1, then DBL_MAX is returned. If s != 1, then         //
//        lambda(s) = [ zeta(s) + eta(s) ] / 2 is returned.                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda(s), if lambda(s) > DBL_MAX, then DBL_MAX is returned and        //
//                 if lambda(s) < -DBL_MAX, then -DBL_MAX is returned.        //
//                 If s = 1, then DBL_MAX is returned. ( s = 1 is a simple    //
//                 pole. )                                                    //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double lambda;                                                         //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda function. ) //
//                                                                            //
//     lambda = Dirichlet_Lambda_Function(s);                                 //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Lambda_Function(double s)
{
   long double lambda;

   if (s == 1.0) return DBL_MAX;

   lambda = xDirichlet_Lambda_Function((long double) s);
   if (fabsl(lambda) >= DBL_MAX) return (lambda < 0.0L) ? -DBL_MAX : DBL_MAX;
   return (double) lambda;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Lambda_Function(long double s)                      //
//                                                                            //
//  Description:                                                              //
//     The Dirichlet lambda function is defined as:                           //
//                       lambda(s) = Sum (1/(2k+1)^s),                        //
//     summed over k = 0,... where s is a complex number with Re(s) > 1,      //
//     then analytically continued to the rest of the complex plane save for  //
//     a simple pole at s = 1.                                                //
//     This routine calculates the Dirichlet Lambda function of a real s != 1,//
//     If s = 1, LDBL_MAX is returned.                                        //
//     If lambda(s) > LDBL_MAX, then LDBL_MAX is returned and                 //
//     if lambda(s) < -LDBL_MAX, then -LDBL_MAX is returned. Note that        //
//     lim lambda(s) = -inf where s->1 from the left and lim lambda(s) = inf  //
//     where s->1 from the right.                                             //
//                                                                            //
//     Note that  lambda(s) = [ zeta(s) + eta(s) ] / 2 , where zeta is the    //
//     Riemann zeta function and eta is the Dirichlet eta function.           //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the Dirichlet lambda function restricted to the     //
//        reals. If s = 1, then LDBL_MAX is returned. If s != 1, then         //
//        lambda(s) = [ zeta(s) + eta(s) ] / 2 is returned.                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda(s), if lambda(s) > LDBL_MAX, then LDBL_MAX is returned and      //
//                 if lambda(s) < -LDBL_MAX, then -LDBL_MAX is returned.      //
//                 If s = 1, then LDBL_MAX is returned. ( s = 1 is a simple   //
//                 pole. )                                                    //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double lambda;                                                    //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda function. ) //
//                                                                            //
//     lambda = xDirichlet_Lambda_Function(s);                                //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Lambda_Function(long double s)
{
   long double lambda_star;

   if (s == 1.0L) return LDBL_MAX;
   lambda_star = xDirichlet_Lambda_Star_Function(s);
   if (fabsl(lambda_star) == LDBL_MAX) return lambda_star;
   return 1.0L + lambda_star;
}


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Lambda_Star_Function(double s)                            //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet lambda star function defined as  //
//                         lambda*(s) = lambda(s) - 1,                        //
//     where lambda() is the Dirichlet lambda function.                       //
//                                                                            //
//     Note that  lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 , where zeta* is    //
//     the Riemann zeta* function and eta* is the Dirichlet eta* function.    //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        reals.  If s = 1, then DBL_MAX is returned. If s != 1, then         //
//        lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 is returned.                //
//                                                                            //
//  Return Value:                                                             //
//     lambda*(s), if lambda*(s) > DBL_MAX, then DBL_MAX is returned and      //
//                 if lambda*(s) < -DBL_MAX, then -DBL_MAX is returned.       //
//                 If s = 1, then DBL_MAX is returned. ( s = 1 is a simple    //
//                 pole. )                                                    //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double lambda_star;                                                    //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     lambda_star = Dirichlet_Lambda_Star_Function(s);                       //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Lambda_Star_Function(double s)
{
   long double lambda_star;

   if (s == 1.0) return DBL_MAX;
   lambda_star = xDirichlet_Lambda_Star_Function((long double) s);
   if (fabsl(lambda_star) >= DBL_MAX)
      return (lambda_star < 0.0L) ? -DBL_MAX : DBL_MAX;
   return (double) lambda_star;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Lambda_Star_Function(long double s)                 //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet lambda star function defined as  //
//                         lambda*(s) = lambda(s) - 1,                        //
//     where lambda() is the Dirichlet lambda function.                       //
//                                                                            //
//     Note that  lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 , where zeta* is    //
//     the Riemann zeta* function and eta* is the Dirichlet eta* function.    //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        reals.  If s = 1, then LDBL_MAX is returned. If s != 1, then        //
//        lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 is returned.                //
//                                                                            //
//  Return Value:                                                             //
//     lambda*(s), if lambda*(s) > LDBL_MAX, then LDBL_MAX is returned and    //
//                 if lambda*(s) < -LDBL_MAX, then -LDBL_MAX is returned.     //
//                 If s = 1, then LDBL_MAX is returned. ( s = 1 is a simple   //
//                 pole. )                                                    //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double lambda_star;                                               //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     lambda_star = xDirichlet_Lambda_Star_Function(s);                      //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Lambda_Star_Function(long double s)
{
   if (s == 0.0L) return -1.0L;
   if (s > 0.0L) return Lambda_Star_Positive_Arg(s);
   else return Lambda_Star_Negative_Arg(s);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Lambda_Star_Positive_Arg(long double s)                 //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda* function which for s > 1 //
//        lambda(s) = Sum (1/(2k+1)^s), summed over k = 1,... .               //
//        lambda(s) = 0 for s = 0.                                            //
//        lambda(s) = LDBL_MAX for s = 1, note that s = 1 is a simple pole of //
//     the Dirichlet lambda function and lim lambda(s) = -inf where s->1 from //
//     the left and lim lambda(s) = inf where s->1 from the right.            //
//                                                                            //
//     Note that  lambda*(s) = [ ( 2^s - 1 )/ (2^s) ] zeta*(s) - 1 / 2^s,     //
//     where zeta*(s) is Riemann zeta* function.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        non-negative real numbers.  If s = 1, then LDBL_MAX is returned.    //
//                                                                            //
//  Return Value:                                                             //
//     lambda*(s), if s = 1, then LDBL_MAX is returned.                       //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double lambda_star;                                               //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     lambda_star = Lambda_Star_Positive_Arg(s);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Lambda_Star_Positive_Arg(long double s)
{
   long double zeta_star;
   long double two_s;

   if (s == 1.0L) return LDBL_MAX;
   zeta_star = xRiemann_Zeta_Star_Function(s);

   if (s > LDBL_MAX_EXP) return zeta_star;
   two_s = powl(2.0L, s);

   return ( (two_s - 1.0L) * zeta_star - 1.0L ) / two_s;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Lambda_Star_Negative_Arg(long double s)                 //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda* function for s < 0.      //
//                                                                            //
//     Note that  lambda*(s) = [ ( 2^s - 1 )/ (2^s) ] zeta*(s) - 1 / 2^s,     //
//     where zeta*(s) is Riemann zeta* function.                              //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        negative real numbers.  If lambda*(s) > LDBL_MAX, then LDBL_MAX is  //
//        returned and if lambda*(s) < -LDBL_MAX, then -LDBL_MAX is returned. //
//                                                                            //
//  Return Value:                                                             //
//     lambda(s)                                                              //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double lambda_star;                                               //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     lambda_star = Lambda_Star_Negative_Arg(s);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Lambda_Star_Negative_Arg(long double s)
{
   long double zeta_star;
   long double lambda_star;
   long double two_minus_s;

   zeta_star = xRiemann_Zeta_Star_Function(s);
   if (fabsl(zeta_star) == LDBL_MAX) return -zeta_star;
   two_minus_s = powl(2.0L, -s);
   if ( fabsl(zeta_star) > ( LDBL_MAX / (two_minus_s - 1.0L) ) ) 
      return (zeta_star < 0.0L) ? LDBL_MAX : -LDBL_MAX;
   lambda_star = (1.0L - two_minus_s) * zeta_star;
   if ( lambda_star > 0.0L ) return lambda_star - two_minus_s;
   if ( lambda_star < (two_minus_s - LDBL_MAX) ) return -LDBL_MAX;
   return lambda_star - two_minus_s;
}
