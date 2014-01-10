////////////////////////////////////////////////////////////////////////////////
// File: chi_square_distribution_large_dof.c                                  //
// Routine(s):                                                                //
//    Chi_Square_Distribution_Large_dof                                       //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                         // required for sqrt() and pow()

//                         Externally Defined Routines                        //

extern double Gaussian_Distribution(double x);

////////////////////////////////////////////////////////////////////////////////
// double Chi_Square_Distribution_Large_dof( double x, int n )                //
//                                                                            //
//  Description:                                                              //
//     If X[1], ..., X[n] are independent N[0,1] distributed random variables,//
//     then the random variable Chi^2 = X[1]^2 + ... + X[n]^2 has a Chi Square//
//     distribution with n degrees of freedom.  Mathematically the Chi Square //
//     distribution with n degrees of freedom is equivalent to a Gamma        //
//     distribution with shape parameter n/2 and scale parameter 2.           //
//                                                                            //
//     The Chi-Square distribution is the Pr[Chi^2 < x] which equals the      //
//     integral from -inf to x of the density                                 //
//                               0                              if x < 0,     //
//       [1 / (2^(n/2) * Gamma(n/2))] * x^(n/2-1) * exp(-x/2)   if x >= 0,    //
//     where n >= 1 and Gamma() is the gamma function.                        //
//                                                                            //
//     For a large number of degrees of freedom, n >> 1, the random variable  //
//                Z = [ Chi^2 / n - (1 - 2/9n) ] / sqrt(2/9n)                 //
//     has, approximately, a standard normal distribution, N[0,1].            //
//                                                                            //
//     This routine returns the approximation to Pr[Chi^2 < x] by             //
//     Pr[Z < z] where z = [ x/n - (1-2/9n) ]/sqrt(2/9n).                     //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the chi-square density   //
//                given above.                                                //
//     int    n   The number of degrees of freedom, n should be large.        //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     p = Chi_Square_Distribution_Large_dof(x, n);                           //
////////////////////////////////////////////////////////////////////////////////

double Chi_Square_Distribution_Large_dof(double x, int n)
{
   double var = 2.0 / (9.0 * (double) n);
   double sd  = sqrt(var);
   const double one_third = 1.0 / 3.0;

   if ( x <= 0.0 ) return 0.0;

   return Gaussian_Distribution( (pow(x/n,one_third) - (1.0 - var)) / sd);
}
