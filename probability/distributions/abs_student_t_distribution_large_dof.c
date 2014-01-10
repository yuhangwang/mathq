////////////////////////////////////////////////////////////////////////////////
// File: abs_student_t_distribution_large_dof.c                               //
// Routine(s):                                                                //
//    Absolute_Student_t_Distribution_Large_dof                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                          // required fabs() and sqrt()

//                         Externally Defined Routines                        //
extern double Gaussian_Distribution(double x);

////////////////////////////////////////////////////////////////////////////////
// double Absolute_Student_t_Distribution_Large_dof( double x, int n )        //
//                                                                            //
//  Description:                                                              //
//     The Absolute Student-t distribution is the Pr[|T| < |x|] =             //
//     Pr[-|x| < T < |x|] where T has a Students-t distribution.              //
//     This distribution corresponds to the integral from -|x| to |x| of the  //
//     density                                                                //
//            1/ [n^(1/2)* B(1/2,n/2)] * (1 + x^2/n)^(-(n+1)/2)               //
//     where n >= 1 and B(,) is the (complete) beta function.                 //
//                                                                            //
//     By making the change of variables: g = n / (n + x^2),                  //
//                   |t|(x,n) = 1 - B(n / (n + x^2), n/2, 1/2)                //
//     where B(,,) is the incomplete beta function.                           //
//                                                                            //
//     For a large number of degrees of freedom, n >> 1, the random variable  //
//                   Z = [ T*(1-1/4n) ] / sqrt(1 + T^2/2n)                    //
//     has, approximately, a standard normal distribution, N[0,1], where T    //
//     has a student t distribution with n degrees of freedom.                //
//                                                                            //
//     This routine returns the approximation to Pr[|T| < |x|] by             //
//     2Pr[Z < z]-1 where z = [ x*(1-1/4n) ] / sqrt(1 + x^2/2n).              //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the absolute student t   //
//                density given above.                                        //
//     int    n   The number of degrees of freedom, n should be large.        //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    n;                                                              //
//                                                                            //
//     p = Absolute_Student_t_Distribution_Large_dof(x, n);                   //
////////////////////////////////////////////////////////////////////////////////

double Absolute_Student_t_Distribution_Large_dof(double x, int n)
{
   double t = fabs(x);
   double p = Gaussian_Distribution( t * (1.0 - 0.25/n)/sqrt(1.0 + 0.5*x*x/n));

   return p + p - 1.0; 
}
