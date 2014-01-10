////////////////////////////////////////////////////////////////////////////////
// File: student_t_distribution_large_dof.c                                   //
// Routine(s):                                                                //
//    Student_t__Distribution_Large_dof                                       //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                          // required for sqrt()

//                         Externally Defined Routines                        //
extern double Gaussian_Distribution(double x);

////////////////////////////////////////////////////////////////////////////////
// double Student_t_Distribution_Large_dof( double x, int n )                 //
//                                                                            //
//  Description:                                                              //
//     If X and X2 are independent random variables, X being N(0,1) and       //
//     X2 being Chi_Square with n degrees of freedom, then the random         //
//     variable T = X / sqrt(X2/n) has a Student-t distribution with n        //
//     degrees of freedom.                                                    //
//                                                                            //
//     The Student-t distribution is the Pr[T < x] which equals the integral  //
//     from -inf to x of the density                                          //
//            1/ [n^(1/2)* B(1/2,n/2)] * (1 + x^2/n)^(-(n+1)/2)               //
//     where n >= 1 and B(,) is the (complete) beta function.                 //
//                                                                            //
//     By making the change of variables: g = n / (n + x^2),                  //
//                   t(x,n) = 1 - B(n / (n + x^2), n/2, 1/2) / 2              //
//     where B(,,) is the incomplete beta function.                           //
//                                                                            //
//     For a large number of degrees of freedom, n >> 1, the random variable  //
//                   Z = [ t*(1-1/4n) ] / sqrt(1 + t^2/2n)                    //
//     has, approximately, a standard normal distribution, N[0,1].            //
//                                                                            //
//     This routine returns the approximation to Pr[T < x] by                 //
//     Pr[Z < z] where z = [ x*(1-1/4n) ] / sqrt(1 + x^2/2n).                 //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the student t density    //
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
//     p = Student_t_Distribution_Large_dof(x, n);                            //
////////////////////////////////////////////////////////////////////////////////

double Student_t_Distribution_Large_dof(double x, int n)
{
   return Gaussian_Distribution( x * (1.0 - 0.25/n)/sqrt(1.0 + 0.5*x*x/n));
}
