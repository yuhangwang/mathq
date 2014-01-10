////////////////////////////////////////////////////////////////////////////////
// File: laplace_dist_tables.c                                                //
// Routine(s):                                                                //
//    Laplace_Distribution_Tables                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for exp() and fabs()

////////////////////////////////////////////////////////////////////////////////
// void Laplace_Distribution_Tables( double start, double delta, int nsteps,  //
//                           double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     Laplace's distribution is the integral from -inf to x of the density   //
//                           f(x) = (1/2) exp(-|x|)                           //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = (1/2) exp(x)      if x < 0 and     //
//                                       = 1 - (1/2) exp(-x) if x > 0.        //
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
//     Laplace_Distribution_Tables(start, delta, N, density,                  //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Laplace_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   int i;

   for (i = 0; i <= nsteps; i++) {
      density[i] = 0.5 * exp(-fabs(x));
      distribution_function[i] = (x <= 0.0) ? density[i] : 1.0 - density[i];
      x += delta;
   }
}
