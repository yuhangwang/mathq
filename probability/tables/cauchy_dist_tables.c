////////////////////////////////////////////////////////////////////////////////
// File: cauchy_dist_tables.c                                                 //
// Routine(s):                                                                //
//    Cauchy_Distribution_Tables                                              //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for atan() and M_PI
                                    
#define PI M_PI

static double const r_pi = 1.0 / PI;

////////////////////////////////////////////////////////////////////////////////
// void Cauchy_Distribution_Tables( double start, double delta, int nsteps,   //
//                           double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The Cauchy distribution is the integral from -inf to x of the density  //
//                          f(x) =  1 / [PI * ( 1 + x^2)],                    //
//     i.e.                                                                   //
//                   Pr[X < x] = F(x) = 1 / 2 + (1 / PI) arctan(x).           //
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
//     double start, delta;                                                   //
//     double density[N+1];                                                   //
//     double distribution_function[N+1];                                     //
//                                                                            //
//     Cauchy_Distribution_Tables(start, delta, N, density,                   //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Cauchy_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   int i;

   for (i = 0; i <= nsteps; i++) {
      density[i] = r_pi / (1.0 + x * x);
      distribution_function[i] = 0.5 + r_pi * atan(x);
      x += delta;
   }
}
