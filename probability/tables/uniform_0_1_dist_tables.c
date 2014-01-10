////////////////////////////////////////////////////////////////////////////////
// File: uniform_0_1_dist_tables.c                                            //
// Routine(s):                                                                //
//    Uniform_0_1_Distribution_Tables                                         //
////////////////////////////////////////////////////////////////////////////////

#include <float.h>                     // required for DBL_MAX

////////////////////////////////////////////////////////////////////////////////
// void Uniform_0_1_Distribution_Tables( double start, double delta,          //
//               int nsteps, double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The Uniform (0,1) distribution is the integral from -inf to x of the   //
//     density                                                                //
//                        f(x) = 1  if 0 < x < 1                              //
//                             = 0  if x < 0 or x > 1                         //
//                        undefined at x = 0 and x = 1.                       //
//     i.e.                                                                   //
//                   Pr[X < x] = F(x) = 0   for x < 0                         //
//                                    = x   for 0 <= x <= 1                   //
//                                    = 1   for x > 1.                        //
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
//     Unform_0_1_Distribution_Tables(start, delta, N, density,               //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Uniform_0_1_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   int i;

   for (i = 0; i <= nsteps; i++) {
      if (x <= 0.0) {
         density[i] = (x == 0.0) ? DBL_MAX : 0.0;
         distribution_function[i] = 0.0;
      } else if (x >= 1.0) {
         density[i] = (x == 1.0) ? -DBL_MAX : 0.0;
         distribution_function[i] = 1.0;
      } else {
         density[i] = 1.0;
         distribution_function[i] = x;
      }
      x += delta;
   }
}
