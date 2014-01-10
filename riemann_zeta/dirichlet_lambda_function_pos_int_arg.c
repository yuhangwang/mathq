////////////////////////////////////////////////////////////////////////////////
// File: dirichlet_lambda_function_pos_int_arg.c                              //
// Routine(s):                                                                //
//    Dirichlet_Lambda_Function_pos_int_arg                                   //
//    xDirichlet_Lambda_Function_pos_int_arg                                  //
//    Dirichlet_Lambda_Star_Function_int_arg                                  //
//    xDirichlet_Lambda_Star_Function_int_arg                                 //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // 
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double Dirichlet_Lambda_Function_pos_int_arg(int s);
long double xDirichlet_Lambda_Function_pos_int_arg(int s);
double Dirichlet_Lambda_Star_Function_pos_int_arg(int s);
long double xDirichlet_Lambda_Star_Function_pos_int_arg(int s);

//                         Externally Defined Routines                        //

long double xRiemann_Zeta_Star_Function_pos_int_arg(int s);
long double xDirichlet_Eta_Star_Function_pos_int_arg(int s);


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Lambda_Function_pos_int_arg(int s)                        //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda function                  //
//        lambda(s) = Sum (1/(2k+1)^s), summed over k = 0,... for s > 1 an    //
//     integer.                                                               //
//        lambda(s) = 0 for s = 0.                                            //
//        lambda(s) = DBL_MAX for s = 1, note that s = 1 is a simple pole of  //
//     the Dirichlet lambda function and lim lambda(s) = -inf where s->1 from //
//     the left and lim lambda(s) = inf where s->1 from the right.            //
//                                                                            //
//     Note that  lambda(s) = [ zeta(s) + eta(s) ] / 2 , where zeta is the    //
//     Riemann zeta function and eta is the Dirichlet eta function.           //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Dirichlet lambda function restricted to the     //
//        nonnegative integers.  If s = 0, then lambda(0) = 0 is returned.    //
//        If s = 1, then DBL_MAX is returned. If s > 1, then Sum(1/(2k+1)^s)  //
//        summed for k = 0,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda(s), if s = 1, then DBL_MAX is returned ( s = 1 is a simple      //
//     pole).                                                                 //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double lambda;                                                         //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda function. ) //
//                                                                            //
//     if ( s >= 0 ) lambda = Dirichlet_Lambda_Function_pos_int_arg(s);       //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Lambda_Function_pos_int_arg(int s)
{
   if (s == 1) return DBL_MAX;
   return (double) xDirichlet_Lambda_Function_pos_int_arg(s);
}

////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Lambda_Function_pos_int_arg(int s)                  //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda function                  //
//        lambda(s) = Sum (1/(2k+1)^s), summed over k = 0,... for s > 1 an    //
//     integer.                                                               //
//        lambda(s) = 0 for s = 0.                                            //
//        lambda(s) = LDBL_MAX for s = 1, note that s = 1 is a simple pole of //
//     the Dirichlet lambda function and lim lambda(s) = -inf where s->1 from //
//     the left and lim lambda(s) = inf where s->1 from the right.            //
//                                                                            //
//     Note that  lambda(s) = [ zeta(s) + eta(s) ] / 2 , where zeta is the    //
//     Riemann zeta function and eta is the Dirichlet eta function.           //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Dirichlet lambda function restricted to the     //
//        nonnegative integers.  If s = 0, then lambda(0) = 0 is returned.    //
//        If s = 1, then LDBL_MAX is returned. If s > 1, then Sum(1/(2k+1)^s) //
//        summed for k = 0,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda(s), if s = 1, then LDBL_MAX is returned ( s = 1 is a simple     //
//     pole).                                                                 //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double lambda;                                                    //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda function. ) //
//                                                                            //
//     if ( s >= 0 ) lambda = xDirichlet_Lambda_Function_pos_int_arg(s);      //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Lambda_Function_pos_int_arg(int s)
{
   if (s <= 0) return 0.0L;
   if (s == 1) return LDBL_MAX;
   return 1.0L + xDirichlet_Lambda_Star_Function_pos_int_arg(s); 
}


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Lambda_Star_Function_pos_int_arg(int s)                   //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda star function defined for //
//     integer arguments s > 1 as lambda*(s) = lambda(s) - 1, where lambda()  //
//     is the Dirichlet lambda function, i.e. lambda*(s) = Sum (1/(2k+1)^s),  //
//     summed over k = 1,... for s > 1 an integer.                            //
//        lambda(s) = -1 for s = 0.                                           //
//        lambda(s) = DBL_MAX for s = 1, note that s = 1 is a simple pole of  //
//     the Dirichlet lambda function and lim lambda(s) = -inf where s->1 from //
//     the left and lim lambda(s) = inf where s->1 from the right.            //
//                                                                            //
//     Note that  lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 , where zeta* is the//
//     Riemann zeta* function and eta* is the Dirichlet eta* function.        //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        nonnegative integers.  If s = 0, then lambda*(0) = -1 is returned.  //
//        If s = 1, then DBL_MAX is returned. If s > 1, then Sum(1/(2k+1)^s)  //
//        summed for k = 1,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda*(s), if s = 1, then DBL_MAX is returned ( s = 1 is a simple     //
//     pole).                                                                 //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double lambda_star;                                                    //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     if ( s >= 0 )                                                          //
//         lambda_star = Dirichlet_Lambda_Star_Function_pos_int_arg(s);       //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Lambda_Star_Function_pos_int_arg(int s)
{
   if (s == 1) return DBL_MAX;
   if (s == 0) return -1.0;
   return (double) xDirichlet_Lambda_Star_Function_pos_int_arg(s);
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Lambda_Star_Function_pos_int_arg(int s)             //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet Lambda star function defined for //
//     integer arguments s > 1 as lambda*(s) = lambda(s) - 1, where lambda()  //
//     is the Dirichlet lambda function, i.e. lambda*(s) = Sum (1/(2k+1)^s),  //
//     summed over k = 1,... for s > 1 an integer.                            //
//        lambda(s) = -1 for s = 0.                                           //
//        lambda(s) = LDBL_MAX for s = 1, note that s = 1 is a simple pole of //
//     the Dirichlet lambda function and lim lambda(s) = -inf where s->1 from //
//     the left and lim lambda(s) = inf where s->1 from the right.            //
//                                                                            //
//     Note that  lambda*(s) = [ zeta*(s) + eta*(s) ] / 2 , where zeta* is the//
//     Riemann zeta* function and eta* is the Dirichlet eta* function.        //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Dirichlet lambda* function restricted to the    //
//        nonnegative integers.  If s = 0, then lambda*(0) = -1 is returned.  //
//        If s = 1, then LDBL_MAX is returned. If s > 1, then Sum(1/(2k+1)^s) //
//        summed for k = 1,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     lambda*(s), if s = 1, then LDBL_MAX is returned ( s = 1 is a simple    //
//     pole).                                                                 //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double lambda_star;                                                    //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     if ( s >= 0 )                                                          //
//         lambda_star = Dirichlet_Lambda_Star_Function_pos_int_arg(s);       //
//     else printf(" Error argument, %d, is negative\n",s);                   //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double lambda_star;                                               //
//                                                                            //
//     ( User code to set s, the argument to the Dirichlet lambda* function. )//
//                                                                            //
//     if ( s >= 0 )                                                          //
//         lambda_star = xDirichlet_Lambda_Star_Function_pos_int_arg(s);      //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Lambda_Star_Function_pos_int_arg(int s)
{
   if (s == 1) return LDBL_MAX;
   if (s == 0) return -1.0L;
   return ( xRiemann_Zeta_Star_Function_pos_int_arg(s) 
           + xDirichlet_Eta_Star_Function_pos_int_arg(s) ) / 2.0L;
}
