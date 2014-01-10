////////////////////////////////////////////////////////////////////////////////
// File: kumaraswamys_dist_tables.c                                           //
// Routine(s):                                                                //
//    Kumaraswamys_Distribution_Tables                                        //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow(), exp(), and log()
#include <float.h>                   // required for DBL_MAX
                                    
////////////////////////////////////////////////////////////////////////////////
// void Kumaraswamys_Distribution_Tables( double a, double b, double start,   //
//  double delta, int nsteps, double *density, double* distribution_function) //
//                                                                            //
//  Description:                                                              //
//     The Kumarasway's distribution is the integral from -inf to x of the    //
//     density                                                                //
//            f(x) = a b x^(a-1) (1 - x^a)^(b-1) for 0 < x < 1,               //
//                 = 0                           elsewhere                    //
//     i.e.                                                                   //
//            Pr[X < x] = F(x) = 0                 for x <= 0,                //
//                             = 1 - (1 - x^a)^b   for 0 < x < 1,             //
//                             = 1                 for x >= 1.                //
//     The shape parameters, a and b, must be positive.                       //
//     If 0 < a < 1, the probability density function has a discontinuity at  //
//     x = 0 and if 0 < b < 1, it has a discontinuity at x = 1.               //
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
//        A shape parameter, b > 0.  The exponent of (1 - x^a).               //
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
//     double a, b;                                                           //
//     double start, delta;                                                   //
//     double density[N+1];                                                   //
//     double distribution_function[N+1];                                     //
//                                                                            //
//     Kumaraswamys_Distribution_Tables(a, b, start, delta, N, density,       //
//                                                   distribution_function);  //
////////////////////////////////////////////////////////////////////////////////

void Kumaraswamys_Distribution_Tables(double a, double b,  double start,
    double delta, int nsteps, double *density, double* distribution_function)
{
   double x = start;
   double temp;
   double tempa;
   double tempb;
   int i;

   for (i = 0; i <= nsteps; i++) {
      if ( x <= 0.0 ) {
         density[i] = 0.0;
         distribution_function[i] = 0.0;
      } else if (x >= 1.0) {
         density[i] = 0.0;
         distribution_function[i] = 1.0;
      } else {
         tempa = pow(x, a - 1.0);
         temp = 1.0 - tempa * x;
         tempb = pow(temp, b - 1.0);
         density[i] = a * b * tempa * tempb;
         distribution_function[i] = 1.0 - temp * tempb;
      }
      x += delta;
   }
}
