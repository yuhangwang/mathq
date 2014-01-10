////////////////////////////////////////////////////////////////////////////////
// File: pareto_dist_tables.c                                                 //
// Routine(s):                                                                //
//    Pareto_Distribution_Tables                                              //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow(), exp(), and log()
#include <float.h>                   // required for DBL_MAX
                                    
////////////////////////////////////////////////////////////////////////////////
// void Pareto_Distribution_Tables( double a, double start, double delta,     //
//               int nsteps, double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The Pareto distribution is the integral from -inf to x of the density  //
//                                           0                if x < 1,       //
//                                 f(x) =  a / x^(a+1)        if x >= 1       //
//     i.e.                                                                   //
//                           Pr[X < x] = F(x) = 0             if x < 1.       //
//                           Pr[X < x] = F(x) = 1 - 1 / x^a   if x > 1.       //
//     The parameter, a > 0, is called the shape parameter.                   //
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
//        The shape parameter, a > 0.                                         //
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
//     Pareto_Distribution_Tables(a, start, delta, N, density,                //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Pareto_Distribution_Tables(double a,  double start, double delta,
                  int nsteps, double *density, double* distribution_function)
{
   double x = start;
   double temp;
   int i;

   for (i = 0; i <= nsteps; i++) {
      if (x < 1.0) {
         density[i] = 0.0;
         distribution_function[i] = 0.0;
      }
      else {
         temp = 1.0 / pow(x,a + 1.0);
         density[i] = a * temp;
         distribution_function[i] = 1.0 - x * temp;
      }
      x += delta;
   }
}
