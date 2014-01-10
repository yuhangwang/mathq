////////////////////////////////////////////////////////////////////////////////
// File: gumbels_maximum_dist_tables.c                                        //
// Routine(s):                                                                //
//    Gumbels_Maximum_Distribution_Tables                                     //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                                   // required for exp()
                                    
////////////////////////////////////////////////////////////////////////////////
// void Gumbels_Maximum_Distribution_Tables( double start, double delta,      //
//               int nsteps, double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     Gumbel's maximum extreme value distribution is the integral from -inf  //
//     to x of the probability density function                               //
//                        f(x) = exp(-x) exp(-exp(-x))                        //
//     i.e.                                                                   //
//                      Pr[X < x] = F(x) = exp(-exp(-x))                      //
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
//     Gumbels_Maximum_Distribution_Tables(start, delta, N, density,          //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Gumbels_Maximum_Distribution_Tables( double start, double delta,
                  int nsteps, double *density, double* distribution_function)
{
   double x = start;
   double temp;
   int i;

   for (i = 0; i <= nsteps; i++) {
      temp = exp(-x);
      distribution_function[i] = exp(-temp);
      density[i] = temp * distribution_function[i];
      x += delta;
   }
}
