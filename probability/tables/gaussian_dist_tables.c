////////////////////////////////////////////////////////////////////////////////
// File: gaussian_dist_tables.c                                               //
// Routine(s):                                                                //
//    Gaussian_Distribution_Tables                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>              // required for exp(), erf(), M_SQRT1_2 and
                               // M_2_SQRTPI
                                    
static double const normalization = 0.5 * M_SQRT1_2 * M_2_SQRTPI;

////////////////////////////////////////////////////////////////////////////////
// void Gaussian_Distribution_Tables( double start, double delta, int nsteps, //
//                           double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The Gaussian distribution function F(x) is the integral from -inf to x //
//     of the probability density function                                    //
//                     f(x) =  1 / sqrt(2*pi) e^-x^2/2.                       //
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
//     Gaussian_Distribution_Tables(start, delta, N, density,                 //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Gaussian_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   int i;

   for (i = 0; i <= nsteps; i++) {
      density[i] = normalization * exp( - 0.5 * x * x );
      distribution_function[i] = 0.5 * ( 1.0 + erf( M_SQRT1_2 * x ) );
      x += delta;
   }
}
