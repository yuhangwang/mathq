////////////////////////////////////////////////////////////////////////////////
// File: chi_square_distribution.c                                            //
// Routine(s):                                                                //
//    Chi_Square_Distribution                                                 //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //

extern double Gamma_Distribution(double x, double a);

////////////////////////////////////////////////////////////////////////////////
// double Chi_Square_Distribution( double x, int n )                          //
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
//     By making the change of variables: y = x / 2,                          //
//                          Chi^2(x,n) = Gamma(x/2,n/2),                      //
//     where Gamma(x,a) is the Gamma distribution with shape parameter a and  //
//     unit scale parameter.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the density given above. //
//     int    n   The number of degrees of freedom.                           //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     p = Chi_Square_Distribution(x, n);                                     //
////////////////////////////////////////////////////////////////////////////////

double Chi_Square_Distribution(double x, int n)
{
   if ( x <= 0.0 ) return 0.0;

   return Gamma_Distribution( 0.5 * x, 0.5 * (double) n);
}
