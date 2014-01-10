////////////////////////////////////////////////////////////////////////////////
// File: chi_square_dist_tables.c                                             //
// Routine(s):                                                                //
//    Chi_Square_Distribution_Tables                                          //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                     // required for exp() and log().
#include <float.h>                    // required for DBL_MAX

//                         Externally Defined Routines                        //

extern double Gamma_Distribution(double x, double nu);
extern double Ln_Gamma_Function(double nu);
                                    
////////////////////////////////////////////////////////////////////////////////
// void Chi_Square_Distribution_Tables( int n, double start, double delta,    //
//               int nsteps, double *density, double* distribution_function)  //
//                                                                            //
//  Description:                                                              //
//     The Chi Square distribution is the integral from -inf to x of the      //
//     density                                                                //
//                               0                                 if x < 0,  //
//     f(x) = [1 / (2^(n/2) * Gamma(n/2))] * x^(n/2-1) * exp(-x/2) if x >= 0, //
//     where n >= 1 and Gamma() is the gamma function.                        //
//     i.e.                                                                   //
//                         Pr[X < x] = F(x) = 0                    if x < 0.  //
//                                          = Gamma(x/2,n/2)       if x > 0.  //
//     The parameter, n > 0, is called the degrees of freedom.                //
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
//     int    n                                                               //
//        The number of degrees of freedom, n >= 1.                           //
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
//     int n;                                                                 //
//                                                                            //
//     Gamma_Distribution_Tables(n, start, delta, N, density,                 //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Chi_Square_Distribution_Tables(int n,  double start, double delta,
                  int nsteps, double *density, double* distribution_function)
{
   double x = start;
   double temp;
   double n2 = 0.5 * (double) n;
   double x2;
   int i;

   for (i = 0; i <= nsteps; i++) {
      if (x < 0.0) {
         density[i] = 0.0;
         distribution_function[i] = 0.0;
      }
      else if (x == 0.0) {
         distribution_function[i] = 0.0;
         if (n == 1) density[i] = DBL_MAX;
         else if (n == 2) density[i] = 0.5; 
         else density[i] = 0.0;
      } else {
         x2 = 0.5 * x;
         distribution_function[i] = Gamma_Distribution(x2, n2);
         temp = (n2 - 1.0) * log(x2) - x2 - Ln_Gamma_Function(n2);
         density[i] = 0.5 * exp(temp);
      }
      x += delta;
   }
}
