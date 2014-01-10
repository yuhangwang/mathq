////////////////////////////////////////////////////////////////////////////////
// File: exponential_dist_tables.c                                            //
// Routine(s):                                                                //
//    Exponential_Distribution_Tables                                         //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for exp()
                                    
////////////////////////////////////////////////////////////////////////////////
// void Exponential_Distribution_Tables( double start, double delta,          //
//               int nsteps, double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The exponential distribution is the integral from -inf to x of the     //
//     probability density                                                    //
//                            f(x) =  exp(-x) for x >= 0.                     //
//                                    0       for x < 0.                      //
//     i.e.                                                                   //
//                   Pr[X < x] = F(x) = 1 - exp(-x) for x > 0                 //
//                                      0           for x <= 0.               //
//     This routine returns the probability density in the array density[]    //
//     where                                                                  //
//            density[i] = f(start + i * delta), i = 0,...,nsteps             //
//     and the distribution function in the array distribution_function[]     //
//     where                                                                  //
//     distribution_function[i] = F(start + i * delta), i = 0,...,nsteps.     //
//     Note the size of the arrays density[] and distribution_function[] must //
//     equal or exceed nsteps + 1.                                            //
//                                                                            //
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
//     Exponential_Distribution_Tables(start, delta, N, density,              //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Exponential_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   int i;

   if ( x < 0.0 ) {
      for ( i = 0; i <= nsteps; i++) {
         density[i] = 0.0;
         distribution_function[i] = 0.0;
         x += delta;
         if (x >= 0.0) break;
      }
      for (; i <= nsteps; i++) {
         density[i] = exp(-x);
         distribution_function[i] = 1.0 - density[i];
         x += delta;
      }
   } else { 
      for (i = 0; i <= nsteps; i++) {
         density[i] = exp(-x);
         distribution_function[i] = 1.0 - density[i];
         x += delta;
      }
   }
}
