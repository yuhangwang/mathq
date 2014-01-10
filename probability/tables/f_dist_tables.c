////////////////////////////////////////////////////////////////////////////////
// File: f_dist_tables.c                                                      //
// Routine(s):                                                                //
//    F_Distribution_Tables                                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow(), exp(), and log()
#include <float.h>                   // required for DBL_MAX

//                         Externally Defined Routines                        //
extern double Ln_Beta_Function(double a, double b);
extern double Beta_Distribution(double x, double a, double b);
                                    
////////////////////////////////////////////////////////////////////////////////
// void F_Distribution_Tables( int v1, int v2, double start, double delta,    //
//                int nsteps, double *density, double* distribution_function) //
//                                                                            //
//  Description:                                                              //
//     If X1 and X2 are independent chi-square distributed random variables   //
//     with v1 and v2 degrees of freedom respectively, then the random        //
//     variable F = (X1/v1) / (X2/v2) has an F-distribution with v1 and v2    //
//     degrees of freedom.                                                    //
//                                                                            //
//     The F-distribution is the Pr[F < x] which equals the integral from     //
//     -inf to x of the probability density function                          //
//                               0                           if f < 0,        //
//       [v1^(v1/2) v2^(v2/2) / B(v1/2,v2/2)]                                 //
//                 * f^(v1/2-1) * (v2+v1 f)^(-(v1+v2)/2))    if f >= 0,       //
//     where v1 >= 1, v2 >= 1, and B(,) is the (complete) beta function.      //
//                                                                            //
//     This routine returns the probability density in the array density[]    //
//     where                                                                  //
//            density[i] = f(start + i * delta), i = 0,...,nsteps             //
//     and the distribution function in the array distribution_function[]     //
//     where                                                                  //
//     distribution_function[i] = F(start + i * delta), i = 0,...,nsteps.     //
//     Note the size of the arrays density[] and distribution_function[] must //
//     equal or exceed nsteps + 1.                                            //
//                                                                            //
//  Arguments:                                                                //
//     double a                                                               //
//        A shape parameter, a > 0.  The exponent of x.                       //
//     double b                                                               //
//        A shape parameter, b > 0.  The exponent of (1 - x).                 //
//     double start                                                           //
//        The initial point to start evaluating f(x) and F(x).                //
//     double delta                                                           //
//        The step size between adjacent evaluation points of f(x) and F(x).  //
//     int    nsteps                                                          //
//        The number of steps.                                                //
//     double density[]                                                       //
//        The value of f(start + i * delta), i = 0,...,nsteps.                //
//     double distribution_function[]                                         //
//        The value of F(start + i * delta), i = 0,...,nsteps.                //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     int    v1, v2;                                                         //
//     double start, delta;                                                   //
//     double density[N+1];                                                   //
//     double distribution_function[N+1];                                     //
//                                                                            //
//     F_Distribution_Tables(v1, v2, start, delta, N, density,                //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void F_Distribution_Tables(int v1, int v2,  double start, double delta,
                  int nsteps, double *density, double* distribution_function)
{
   double x = start;
   double ln_density;
   double v12 = (double)v1 / 2.0;
   double v22 = (double)v2 / 2.0;
   double ln_temp;
   double g;
   int i;

   ln_temp = v12 * log((double)v1) + v22 * log((double)v2)
                                                 - Ln_Beta_Function(v12, v22);
   for (i = 0; i <= nsteps; i++) {
      if ( x <= 0.0 ) {
         density[i] = 0.0;
         distribution_function[i] = 0.0;
      } else {
         ln_density = ln_temp + (v12 - 1.0) * log(x)
                                         - (v12 + v22) * log((double)v2+v1*x);
         density[i] = exp(ln_density);
         g = v12 * x;
         distribution_function[i] = Beta_Distribution(g / (v22 + g), v12, v22);
      }
      x += delta;
   }
}
