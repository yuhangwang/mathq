////////////////////////////////////////////////////////////////////////////////
// File: f_distribution.c                                                     //
// Routine(s):                                                                //
//    F_Distribution                                                          //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
extern double Beta_Distribution(double x, double a, double b);

////////////////////////////////////////////////////////////////////////////////
// double F_Distribution( double x, int v1, int v2 )                          //
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
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the density given above. //
//     int   v1   The number of degrees of freedom of the numerator of the    //
//                F-test, i.e. (v1/2 - 1) is the exponent of f in the         //
//                integrand above.  Note v1 >= 1.                             //
//     int   v2   The number of degrees of freedom of the denominator of the  //
//                F-test, i.e. (-(v1+v2)/2) is the exponent of the term       //
//                (v2 + v1 f) in the integrand above.  Note v2 >= 1.          //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    v1, v2;                                                         //
//                                                                            //
//     p = F_Distribution(x, v1, v2);                                         //
////////////////////////////////////////////////////////////////////////////////

double F_Distribution(double f, int v1, int v2)
{
   double a = (double) v1 / 2.0;
   double b = (double) v2 / 2.0;
   double g = a*f;

   if ( f <= 0.0 ) return 0.0;

   return Beta_Distribution( g / (b + g), a, b);
}
