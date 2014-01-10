////////////////////////////////////////////////////////////////////////////////
// File: t2_dist_tables.c                                                     //
// Routine(s):                                                                //
//    t2_Distribution_Tables                                                  //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for sqrt()
                                    
////////////////////////////////////////////////////////////////////////////////
// void t2_Distribution_Tables( double start, double delta, int nsteps,       //
//                           double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The t2 distribution, Student's t with 2 degrees of freedom, is the     //
//     integral from -inf to x of the density                                 //
//                           f(x) =  1 / (2 + x^2)^(3/2),                     //
//     i.e.                                                                   //
//             Pr[X < x] = F(x) = (1/2) (1 + x / sqrt(2 + x^2)).              //
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
//     t2_Distribution_Tables(start, delta, N, density,                       //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////
#define sqrt_2_over_DBL_EPSILON 9.4906265624251559e+07
#define ONE_over_cbrt_DBL_MIN 3.5553731598732436e+102

void t2_Distribution_Tables( double start, double delta, int nsteps,
                              double *density, double* distribution_function)
{
   double x = start;
   double temp;
   double root_temp;
   int i;

   for (i = 0; i <= nsteps; i++, x += delta) {
      if (fabs(x) >= ONE_over_cbrt_DBL_MIN) {
         density[i] = 0.0;
         distribution_function[i] = (x > 0.0 ? 1.0 : 0.0);
         continue;
      }
      temp = 1.0 / (2.0 + x * x);
      root_temp = sqrt(temp);
      density[i] = temp * root_temp;
      if (fabs(x) > sqrt_2_over_DBL_EPSILON) 
         distribution_function[i] = (x > 0.0 ? 1.0 : 0.0);
      else
         distribution_function[i] = 0.5 * (1.0 + x * root_temp);
   }
}
