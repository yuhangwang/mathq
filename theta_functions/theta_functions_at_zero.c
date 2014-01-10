////////////////////////////////////////////////////////////////////////////////
// File: theta_functions_at_zero.c                                            //
// Routine(s):                                                                //
//    Theta_Functions_at_zero.c                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for expl(), sinl(), cosl(), sqrtl(),
                            // and M_1_PI cos()

//                         Internally Defined Routines                        //

static void Small_x( long double x, long double theta[] );
static void Large_x( long double x, long double theta[] );

static long double const pi = 3.1415926535897932384626433832795029L;

////////////////////////////////////////////////////////////////////////////////
// void Theta_Functions_at_zero(double x, double *theta_1, double *theta_2,   //
//                                         double *theta_3, double *theta_4 ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     This function calculates the four theta functions as given in Spanier  //
//     and Oldham at nu = 0.                                                  //
//     The results are not returned in an array in order to avoid the         //
//     potentially confusing theta[0]() = theta_1(), ... .                    //
//     For 0 < x < 1/pi,                                                      //
//      theta_1(x) = 0.0                                                      //
//      theta_2(x) = 1/sqrt(pi*x) Sum [(-1)^j exp(-j^2 / x)]                  //
//      theta_3(x) = 1/sqrt(pi*x) Sum [exp(-j^2 / x)]                         //
//      theta_4(x) = 1/sqrt(pi*x) Sum [exp(-(j + 1/2)^2 / x)]                 //
//     where the sum extends from j = -inf to inf.                            //
//     For 1/pi <= x,                                                         //
//      theta_1(x) = 0.0                                                      //
//      theta_2(x) = 2Sum[exp(-(j + 1/2)^2 pi^2 x)                            //
//     where the sum extends from j = 0 to inf and                            //
//      theta_3(x) = 1 + 2Sum[exp(-j^2 pi^2 x)                                //
//      theta_4(x) = 1 + 2Sum[(-1)^j exp(-j^2 pi^2 x)                         //
//     where the sum extends from j = 1 to inf and                            //
//                                                                            //
//  Arguments:                                                                //
//     double  x                                                              //
//                The second argument of the Theta functions, x > 0.          //
//                Note nu = 0 is the periodic first argument.                 //
//     double* theta_1                                                        //
//                The value of theta_1 as calculated above.                   //
//     double* theta_2                                                        //
//                The value of theta_2 as calculated above.                   //
//     double* theta_3                                                        //
//                The value of theta_3 as calculated above.                   //
//     double* theta_4                                                        //
//                The value of theta_4 as calculated above.                   //
//                                                                            //
//  Return Value:                                                             //
//     The values of theta_1, theta_2, theta_3, and theta_4 are returned      //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     double theta_1, theta_2, theta_3, theta_4;                             //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     Theta_Functions_at_zero(x, &theta_1, &theta_2, &theta_3, &theta_4);    //
////////////////////////////////////////////////////////////////////////////////
                
void Theta_Functions_at_zero(double x, double *theta_1, double *theta_2,
                                             double *theta_3, double *theta_4 )
{
   long double theta[4];

 //  If x < 1 / pi, use the series expansion for small x, otherwise use the //
 //  series expansion for large x.                                          //

   if ( x < M_1_PI ) Small_x(x,theta);
   else Large_x(x,theta);

   *theta_1 = 0.0;
   *theta_2 = (double) theta[1];
   *theta_3 = (double) theta[2];
   *theta_4 = (double) theta[3];

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// static void Small_x(long double x, long double theta[] )                   //
//                                                                            //
//  Description:                                                              //
//     For 0 < x < 1/pi,                                                      //
//      theta_1(x) = 0.0                                                      //
//      theta_2(x) = 1/sqrt(pi*x) Sum [(-1)^j exp(-j^2 / x)]                  //
//      theta_3(x) = 1/sqrt(pi*x) Sum [exp(-j^2 / x)]                         //
//      theta_4(x) = 1/sqrt(pi*x) Sum [exp(-(j + 1/2)^2 / x)]                 //
//     where the sum extends from j = -inf to inf.                            //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//                The second argument of the Theta functions, x > 0.          //
//                Note nu = 0 is the periodic first argument.                 //
//     long double theta[]                                                    //
//                The array of values of the theta functions evaluated at     //
//                nu, x, 0 < x < 1/pi: where theta[0] = theta_1(nu,x),...,    //
//                theta_4(nu,x).                                              //
//                                                                            //
//  Return Value:                                                             //
//     The values of theta_1, theta_2, theta_3, and theta_4 are returned      //
//     via the argument list in theta[].                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
static void Small_x( long double x, long double theta[] )
{
   static int const max_j = 6;  //(int)(1.5 + sqrtl(-logl(LDBL_EPSILON)/pi))+1

   long double t_below[4];
   long double t_above[4];
   long double exponents[2];
   long double term[3];
   int phase;
   int n;
   int j;

   phase = (max_j % 2 == 0) ? 1: -1;
   exponents[1] = - max_j;
   exponents[0] = exponents[1] - 0.5L;
   term[0] = expl(- exponents[0] * exponents[0] / x);
   term[1] = expl(- exponents[1] * exponents[1] / x);
   exponents[0] += 1.0L;
   term[2] = expl(- exponents[0] * exponents[0] / x);
   t_below[1] = phase * term[1];
   t_below[2] = term[1];
   t_below[3] = term[2];
   for (j = 1 - max_j; j <= 0; j++) {
      phase = -phase;
      term[0] = term[2];
      exponents[1] += 1.0L;
      term[1] = expl(- exponents[1] * exponents[1] / x);
      exponents[0] += 1.0L;
      term[2] = expl(- exponents[0] * exponents[0] / x);
      t_below[1] += phase * term[1];
      t_below[2] += term[1];
      t_below[3] += term[2];
   }
   phase = (max_j % 2 == 0) ? 1: -1;
   exponents[1] = max_j;
   exponents[0] = exponents[1] + 0.5L;
   term[2] = expl(- exponents[0] * exponents[0] / x);
   term[1] = expl(- exponents[1] * exponents[1] / x);
   exponents[0] -= 1.0L;
   term[0] = expl(- exponents[0] * exponents[0] / x);
   t_above[1] = phase * term[1];
   t_above[2] = term[1];
   t_above[3] = term[2];
   for (j = max_j - 1; j > 0; j--) {
      phase = -phase;
      term[2] = term[0];
      exponents[1] -= 1.0L;
      term[1] = expl(- exponents[1] * exponents[1] / x);
      exponents[0] -= 1.0L;
      term[0] = expl(- exponents[0] * exponents[0] / x);
      t_above[1] += phase * term[1];
      t_above[2] += term[1];
      t_above[3] += term[2];
   }
   x = 1.0L / sqrtl(pi * x);
   theta[1] = x * (t_above[1] + t_below[1]); 
   theta[2] = x * (t_above[2] + t_below[2]); 
   theta[3] = x * (t_above[3] + t_below[3]); 
}

////////////////////////////////////////////////////////////////////////////////
// void Large_x(long double x, long double theta[] )                          //
//                                                                            //
//  Description:                                                              //
//     For x >= 1/pi,                                                         //
//      theta_1(x) = 0.0                                                      //
//      theta_2(x) = 2Sum[exp(-(j + 1/2)^2 pi^2 x)                            //
//     where the sum extends from j = 0 to inf and                            //
//      theta_3(x) = 1 + 2Sum[exp(-j^2 pi^2 x)                                //
//      theta_4(x) = 1 + 2Sum[(-1)^j exp(-j^2 pi^2 x)                         //
//     where the sum extends from j = 1 to inf and                            //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//                The second argument of the Theta functions, x > 0.          //
//                Note nu = 0 is the periodic first argument.                 //
//     long double theta[]                                                    //
//                The array of values of the theta functions evaluated at     //
//                nu, x, x >= 1/pi: where theta[0] = theta_1(nu,x),...,       //
//                theta_4(nu,x).                                              //
//                                                                            //
//  Return Value:                                                             //
//     The values of theta_1, theta_2, theta_3, and theta_4 are returned      //
//     via the argument list in theta[].                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
static void Large_x( long double x, long double theta[] )
{
   static int const max_j = 5;  //(in1G(0.5 + sqrtl(-logl(LDBL_EPSILON)/pi))+1

   long double t_above[4];
   long double exponents[2];
   long double term[2];
   int phase;
   int j;

   x *= (pi * pi);
   phase = (max_j % 2 == 0) ? 1: -1;

   for (j = 0; j < 4; j++) t_above[j] = 0.0L;

   for (j = max_j; j >= 1; j--) {
      exponents[1] = j;
      exponents[0] = exponents[1] + 0.5L;
      term[0] = expl(- exponents[0] * exponents[0] * x);
      term[1] = expl(- exponents[1] * exponents[1] * x);
      t_above[1] += term[0];
      t_above[2] += term[1];
      t_above[3] += phase * term[1];
      phase = -phase;
   }
   t_above[1] += expl(- 0.25L * x);
   theta[1] = t_above[1] + t_above[1]; 
   theta[2] = 1.0L + t_above[2] + t_above[2]; 
   theta[3] = 1.0L + t_above[3] + t_above[3]; 
   return;   
}
