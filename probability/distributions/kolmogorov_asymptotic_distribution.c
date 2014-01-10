////////////////////////////////////////////////////////////////////////////////
// File: kolmogorov_asymptotic_distribution.c                                 //
// Routine(s):                                                                //
//    Kolmogorov_Asymptotic_Distribution                                      //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// double Kolmogorov_Asymptotic_Distribution( double dn, int sample_size)     //
//                                                                            //
//  Description:                                                              //
//     This function returns the probability Pr[Dn <= dn] where               //
//           dn = sup {|FN(x) - F(x)| : for all x} where FN                   //
//     is the empirical distribution function                                 //
//              FN(x) = #{sample[i]: sample[i] <= x } / sample_size           //
//     and F is the continuous theoretical distribution function              //
//                            F(x) = Pr[X <= x].                              //
//                                                                            //
//     For values z = dn / sqrt(sample_size) <= 1, the probability is         //
//     calculated using the representation of Jacobi's theta 4 function       //
//     for rapid convergence for small z:                                     //
//            [sqrt(2 pi) / z] Sum[exp(-pi^2 (2k+1)^2 / 8 z^2)],              //
//     where the sum extends for k = 0, ... .                                 //
//                                                                            //
//     For values z = dn / sqrt(sample_size) > 1, the probability is          //
//     calculated using the representation of Jacobi's theta 4 function       //
//     for rapid convergence for large z:                                     //
//                    1 + 2 Sum[(-1)^k exp(-2 k^2 z^2)],                      //
//     where the sum extends for k = 1, ... .                                 //
//                                                                            //
//  Arguments:                                                                //
//     double dn                                                              //
//        The supremum of the absolute difference between the empirical       //
//        distribution function and the continuous theoretical distribution   //
//        function.                                                           //
//     int sample_size                                                        //
//        The number of observations.                                         //
//                                                                            //
//  Return Values:                                                            //
//     For a large sample size, the probability of observing a value for the  //
//     supremum of the difference between the empirical distribution function //
//     and the continuous theoretical distribution function is less than or   //
//     equal to the observed value.                                           //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double dn, pr;                                                         //
//                                                                            //
//     (your code to calculate dn )                                           //
//                                                                            //
//     pr = Kolmogorov_Asymptotic_Distribution(dn, N);                        //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                   // required for fabs(), sqrt(), and exp()
#include <float.h>                  // required for DBL_EPSILON, DBL_MIN

double Kolmogorov_Asymptotic_Distribution( double dn, int sample_size )
{
   const double pi= 3.14159265358979323846;
   double z;
   double x;
   double term;
   double exponent;
   double sn;
   double sqrt_2pi= sqrt(2.0 * pi);           
   double upper_cutoff = sqrt( 0.5 * log(2.0 / DBL_EPSILON) );
   double lower_cutoff = pi / sqrt( 8.0 * fabs(log(DBL_MIN)) );
   int k, n;

   if (dn <= 0.0) return 0.0;
   if (sample_size <= 0.0) return 0.0;
   z = sqrt(sample_size) * dn;
   if ( z > upper_cutoff ) return 1.0;
   if ( z < lower_cutoff ) return 0.0;
   if ( z <= 1.0 ) {
      exponent = pi / z;
      exponent *= (- 0.125 * exponent);
      sn = 0.0;
      for (k = 3; k >= 0; k--) {
         n = k + k + 1;
         sn += exp( n * n * exponent);
      }
      return sqrt_2pi * sn / z;
   } 
   x = -2.0 * z * z;
   sn = 0.0;
   for (k = 4; k >= 1; k--) {
      term  = - exp( k * k * x);
      k--;
      term += exp( k * k * x);
      sn += term;
   }
   return (1.0 - 2.0 * sn);
}
