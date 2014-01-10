////////////////////////////////////////////////////////////////////////////////
// File: f_distribution_large_dofs.c                                          //
// Routine(s):                                                                //
//    F_Distribution_Large_dofs                                               //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                         // required for sqrt() and pow() 

//                         Externally Defined Routines                        //
extern double Gaussian_Distribution(double x);

////////////////////////////////////////////////////////////////////////////////
// double F_Distribution_Large_dofs double x, int v1, int v2 )                //
//                                                                            //
//  Description:                                                              //
//     If X1 and X2 are independent chi-square distributed random variables   //
//     with v1 and v2 degrees of freedom respectively, then the random        //
//     variable F = (X1/v1) / (X2/v2) has an F-distribution with v1 and v2    //
//     degrees of freedom.                                                    //
//                                                                            //
//     The F-distribution is the Pr[F < x] which equals the integral from     //
//     -inf to x of the density                                               //
//                               0                           if f < 0,        //
//       [v1^(v1/2) v2^(v2/2) / B(v1/2,v2/2)]                                 //
//                 * f^(v1/2-1) * (v2+v1 f)^(-(v1+v2)/2))   if f >= 0,        //
//     where v1 >= 1, v2 >= 1, and B(,) is the (complete) beta function.      //
//                                                                            //
//     By making the change of variables: g = v1*f / (v2 + v1*f),             //
//                   F(x,v1,v2) = B(v1*x / (v2 + v1*x), v1/2, v2/2),          //
//     where B(,,) is the incomplete beta function.                           //
//                                                                            //
//     For a large number of numerator degrees of freedom, v1 >> 1, and a     //
//     large number of denominator degrees of freedom, v2 >> 1, the random    //
//     variable                                                               //
//               Z = [ X (1-w2) - (1-w1) ] / sqrt(w1 + x^2 w2)                //
//     where X = F^(1/3), w1 = 2/9v1 and w2 = 2/9v2, has, approximately, a    //
//     standard normal distribution.                                          //
//                                                                            //
//     This routine returns the approximation to Pr[F < f] given by           //
//     Pr[Z < z] where after setting x = f^1/3, w1 = 2/9v1, and w2 = 2/9v2    //
//     z = [x*(1-w2)-(1-w1)]/sqrt(w1+x^2*w2).                                 //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the density given above. //
//     int   v1   The number of degrees of freedom of the numerator of the    //
//                F-test, i.e. (v1/2 - 1) is the exponent of f in the         //
//                integrand above.  Note v1 >> 1.                             //
//     int   v2   The number of degrees of freedom of the denominator of the  //
//                F-test, i.e. (-(v1+v2)/2) is the exponent of the term       //
//                (v2 + v1 f) in the integrand above.  Note v2 >> 1.          //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    v1, v2;                                                         //
//                                                                            //
//     p = F_Distribution_Large_dofs(x, v1, v2);                              //
////////////////////////////////////////////////////////////////////////////////

double F_Distribution_Large_dofs(double f, int v1, int v2)
{
   double w1 = (2.0 / 9.0) / (double) v1;
   double w2 = (2.0 / 9.0) / (double) v2;
   double x;

   if ( f <= 0.0 ) return 0.0;

   x = pow( f, 1.0 / 3.0);
   x = ( x * (1.0 - w2) - (1.0 - w1) ) / sqrt(w1 + x * x * w2);
   return Gaussian_Distribution(x);
}
