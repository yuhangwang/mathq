////////////////////////////////////////////////////////////////////////////////
// File: dtheta_functions.c                                                   //
// Routine(s):                                                                //
//    DTheta_Functions                                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>           // required for expl(), sinl(), cosl(), sqrtl(),
                            // and M_1_PI cos()

//                         Internally Defined Routines                        //

static void Small_x( long double nu, long double x, long double dtheta[] );
static void Large_x( long double nu, long double x, long double dtheta[] );

static long double const pi = 3.1415926535897932384626433832795029L;

////////////////////////////////////////////////////////////////////////////////
// void DTheta_Functions(double nu, double x, double *dtheta_1,               //
//                     double *dtheta_2, double *dtheta_3, double *dtheta_4 ) //
//                                                                            //
//  Description:                                                              //
//     There is a plethora of notations for the theta functions, you should   //
//     check the expressions below in order to insure that your arguments are //
//     correct.                                                               //
//     This function calculates the derivatives with respect to nu of the     //
//     four theta functions as given in Spanier and Oldham.                   //
//     The results are not returned in an array in order to avoid the         //
//     potentially confusing dtheta[0]() = dtheta_1(), ..., dtheta_4().       //
//     For 0 < x < 1/pi,                                                      //
//      dtheta_1(nu,x) = -2/(x*sqrt(pi*x)) Sum [(-1)^j (nu +j - 1/2)          //
//                                                exp(-(nu - 1/2 + j)^2 / x)] //
//      dtheta_2(nu,x) = -2/(x*sqrt(pi*x)) Sum [(-1)^j (nu + j)               //
//                                                      exp(-(nu + j)^2 / x)] //
//      dtheta_3(nu,x) = -2/(x*sqrt(pi*x)) Sum (nu + j) [exp(-(nu + j)^2 / x)]//
//      dtheta_4(nu,x) = -2/(x*sqrt(pi*x)) Sum (nu + j + 1/2)                 //
//                                               [exp(-(nu + 1/2 + j)^2 / x)] //
//     where the sum extends from j = -inf to inf.                            //
//     For 1/pi <= x,                                                         //
//      dtheta_1(nu,x) = 2 pi Sum[(-1)^j (2j+1) exp(-(j + 1/2)^2 pi^2 x)      //
//                                                          cos[(2j+1)pi nu]] //
//      dtheta_2(nu,x) = -2 pi Sum[(2j+1) exp(-(j + 1/2)^2 pi^2 x)            //
//                                                          sin[(2j+1)pi nu]] //
//     where the sum extends from j = 0 to inf and                            //
//      dtheta_3(nu,x) = -4 pi Sum[j exp(-j^2 pi^2 x) sin[2j pi nu]]          //
//      dtheta_4(nu,x) = -4 pi Sum[(-1)^j j exp(-j^2 pi^2 x) sin[2j pi nu]]   //
//     where the sum extends from j = 1 to inf and                            //
//                                                                            //
//  Arguments:                                                                //
//     double  nu                                                             //
//                The periodic argument of the Theta functions.               //
//     double  x                                                              //
//                The second argument of the Theta functions, x > 0.          //
//     double* dtheta_1                                                       //
//                The value of dtheta_1 as calculated above.                  //
//     double* dtheta_2                                                       //
//                The value of dtheta_2 as calculated above.                  //
//     double* dtheta_3                                                       //
//                The value of dtheta_3 as calculated above.                  //
//     double* dtheta_4                                                       //
//                The value of dtheta_4 as calculated above.                  //
//                                                                            //
//  Return Value:                                                             //
//     The values of dtheta_1, dtheta_2, dtheta_3, and dtheta_4 are returned  //
//     via the argument list.                                                 //
//                                                                            //
//  Example:                                                                  //
//     double x, nu;                                                          //
//     double dtheta_1, dtheta_2, dtheta_3, dtheta_4;                         //
//                                                                            //
//     ( code to initialize x and nu )                                        //
//                                                                            //
//     DTheta_Functions( nu, x, &dtheta_1, &dtheta_2, &dtheta_3, &dtheta_4);  //
////////////////////////////////////////////////////////////////////////////////
                
void DTheta_Functions(double nu, double x, double *dtheta_1, double *dtheta_2,
                                           double *dtheta_3, double *dtheta_4 )
{
   long double dtheta[4];
   int n;

 // Translate nu by a half periods of theta_1 or theta_2 until 0 <= nu < 1. //

   n = (int) fabs(nu);
   if (nu >= 0.0) nu -= (double) n;
   else {
      n++;
      nu += (double) n;
   }

 //  If x < 1 / pi, use the series expansion for small x, otherwise use the //
 //  series expansion for large x.                                          //

   if ( x < M_1_PI ) Small_x(nu,x,dtheta);
   else Large_x(nu,x,dtheta);

  // If nu was translated by an odd number of half periods, change sign of //
  // theta_1 and theta_2. //

   if (n % 2 == 1) {
      dtheta[0] = - dtheta[0];
      dtheta[1] = - dtheta[1];
   }

// Return results. //

   *dtheta_1 = (double) dtheta[0];
   *dtheta_2 = (double) dtheta[1];
   *dtheta_3 = (double) dtheta[2];
   *dtheta_4 = (double) dtheta[3];

   return; 
}

////////////////////////////////////////////////////////////////////////////////
// static void Small_x(long double nu, long double x, long double dtheta[] )  //
//                                                                            //
//  Description:                                                              //
//     For 0 < x < 1/pi,                                                      //
//      dtheta_1(nu,x) = -2/(x*sqrt(pi*x)) Sum [(-1)^j (nu +j - 1/2)          //
//                                                exp(-(nu - 1/2 + j)^2 / x)] //
//      dtheta_2(nu,x) = -2/(x*sqrt(pi*x)) Sum [(-1)^j (nu + j)               //
//                                                      exp(-(nu + j)^2 / x)] //
//      dtheta_3(nu,x) = -2/(x*sqrt(pi*x)) Sum (nu + j) [exp(-(nu + j)^2 / x)]//
//      dtheta_4(nu,x) = -2/(x*sqrt(pi*x)) Sum (nu + j + 1/2)                 //
//                                               [exp(-(nu + 1/2 + j)^2 / x)] //
//     where the sum extends from j = -inf to inf.                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  nu                                                        //
//                The periodic argument of the Theta functions.               //
//     long double x                                                          //
//                The second argument of the Theta functions, x > 0.          //
//     long double dtheta[]                                                   //
//                The array of values of the derivatives of the theta         //
//                functions with respect to nu evaluated at nu, x,            //
//                0 < x < 1/pi: where dtheta[0] = dtheta_1(nu,x),...,         //
//                dtheta_4(nu,x).                                             //
//                                                                            //
//  Return Value:                                                             //
//     The values of dtheta_1, dtheta_2, dtheta_3, and dtheta_4 are returned  //
//     via the argument list in dtheta[].                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
static void Small_x( long double nu, long double x, long double dtheta[] )
{
   static int const max_j = 8;  //(int)(1.5 + sqrtl(-logl(LDBL_EPSILON)/pi))+1

   long double t_below[4];
   long double t_above[4];
   long double exponents[2];
   long double term[3];
   int sign;
   int phase;
   int n;
   int j;

   n = nu;
   sign = (n % 2 == 0) ? 1 : -1;
   nu -= (long double) n;
   
   phase = (max_j % 2 == 0) ? 1: -1;
   exponents[1] = nu - max_j;
   exponents[0] = exponents[1] - 0.5L;
   term[0] = expl(- exponents[0] * exponents[0] / x);
   term[1] = expl(- exponents[1] * exponents[1] / x);
   exponents[0] += 1.0L;
   term[2] = expl(- exponents[0] * exponents[0] / x);
   t_below[0] = phase * term[0];
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
      t_below[0] += phase * (nu +j - 0.5L) * term[0];
      t_below[1] += phase * (nu + j) * term[1];
      t_below[2] += (nu + j) * term[1];
      t_below[3] += (nu + j + 0.5L) * term[2];
   }
   phase = (max_j % 2 == 0) ? 1: -1;
   exponents[1] = nu + max_j;
   exponents[0] = exponents[1] + 0.5L;
   term[2] = expl(- exponents[0] * exponents[0] / x);
   term[1] = expl(- exponents[1] * exponents[1] / x);
   exponents[0] -= 1.0L;
   term[0] = expl(- exponents[0] * exponents[0] / x);
   t_above[0] = phase * term[0];
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
      t_above[0] += phase * (nu +j - 0.5L) * term[0];
      t_above[1] += phase * (nu + j) * term[1];
      t_above[2] += (nu + j) * term[1];
      t_above[3] += (nu + j + 0.5L) * term[2];
   }
   x = -2.0L / (x * sqrtl(pi * x) );
   dtheta[0] = sign * x * (t_above[0] + t_below[0]); 
   dtheta[1] = sign * x * (t_above[1] + t_below[1]); 
   dtheta[2] = x * (t_above[2] + t_below[2]); 
   dtheta[3] = x * (t_above[3] + t_below[3]); 
}

////////////////////////////////////////////////////////////////////////////////
// void Large_x(long double nu, long double x, long double theta[] )          //
//                                                                            //
//  Description:                                                              //
//     For x >= 1/pi,                                                         //
//      dtheta_1(nu,x) = 2*pi Sum[(-1)^j (2j+1) exp(-(j + 1/2)^2 pi^2 x)      //
//                                                           cos[(2j+1)pi nu]]//
//      dtheta_2(nu,x) = -2*pi Sum[(2j+1) exp(-(j + 1/2)^2 pi^2 x)            //
//                                                          sin[(2j+1)pi nu]] //
//     where the sum extends from j = 0 to inf and                            //
//      dtheta_3(nu,x) = -4*pi Sum[j exp(-j^2 pi^2 x) sin[2j pi nu]]          //
//      dtheta_4(nu,x) = -4*pi Sum[(-1)^j j exp(-j^2 pi^2 x) sin[2j pi nu]]   //
//     where the sum extends from j = 1 to inf and                            //
//                                                                            //
//  Arguments:                                                                //
//     long double  nu                                                        //
//                The periodic argument of the Theta functions.               //
//     long double x                                                          //
//                The second argument of the Theta functions, x > 0.          //
//     long double dtheta[]                                                   //
//                The array of values of the derivatives of the theta         //
//                functions with respect to nu evaluated at nu, x,            //
//                x >= 1/pi: where dtheta[0] = dtheta_1(nu,x),...,            //
//                dtheta_4(nu,x).                                             //
//                                                                            //
//  Return Value:                                                             //
//     The values of dtheta_1, dtheta_2, dtheta_3, and dtheta_4 are returned  //
//     via the argument list in dtheta[].                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
                
static void Large_x( long double nu, long double x, long double dtheta[] )
{
   static int const max_j = 5;  //(in1G(0.5 + sqrtl(-logl(LDBL_EPSILON)/pi))+1

   long double t_above[4];
   long double exponents[2];
   long double term[2];
   long double angles[2];
   long double sin_angle;
   int phase;
   int j;

   nu *= pi;
   x *= (pi * pi);
   phase = (max_j % 2 == 0) ? 1: -1;

   for (j = 0; j < 4; j++) t_above[j] = 0.0L;

   for (j = max_j; j >= 1; j--) {
      exponents[1] = j;
      exponents[0] = exponents[1] + 0.5L;
      term[0] = expl(- exponents[0] * exponents[0] * x);
      term[1] = expl(- exponents[1] * exponents[1] * x);
      angles[1] = (j + j) * nu;
      angles[0] = angles[1] + nu;
      t_above[0] += phase * (j + j + 1) * term[0] * cosl(angles[0]);
      t_above[1] += term[0] * (j + j + 1) * sinl(angles[0]);
      sin_angle = sinl(angles[1]);
      t_above[2] += term[1] * j * sin_angle;
      t_above[3] += phase * term[1] * j * sin_angle;
      phase = -phase;
   }
   t_above[0] += expl(- 0.25L * x) * cosl(nu);
   t_above[1] += expl(- 0.25L * x) * sinl(nu);
   dtheta[0] = (pi + pi) * t_above[0]; 
   dtheta[1] = -(pi + pi) * t_above[1]; 
   dtheta[2] = - 4.0 * pi * t_above[2]; 
   dtheta[3] = - 4.0 * pi * t_above[3]; 
   return;   
}
