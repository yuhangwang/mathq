////////////////////////////////////////////////////////////////////////////////
// File: f_distribution_large_denominator_dof.c                               //
// Routine(s):                                                                //
//    F_Distribution_Large_Denominator_dof                                    //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
extern double Chi_Square_Distribution(double x, int n);

////////////////////////////////////////////////////////////////////////////////
// double F_Distribution_Large_Denominator_dof( double x, int v1, int v2 )    //
//                                                                            //
//  Description:                                                              //
//     If X1 and X2 are independent chi-square distributed random variables   //
//     with v1 and v2 degrees of freedom respectively, then the random        //
//     variable F = (X1/v1) / (X2/v2) has an F-distribution with v1 and v2    //
//     degrees of freedom.   The degree of freedom v1 is called the numerator //
//     degree of freedom and the degree of freedom v2 is called the           //
//     denominator degree of freedom.                                         //
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
//     For a large number of denominator degrees of freedom, v2 >> 1, the     //
//     random variable                                                        //
//                              X = [ v1 F ]                                  //
//     has, approximately, a chi-square distribution with v1 degrees of       //
//     freedom and Pr[ F < f ] = Pr[ X < v1*f].                               //
//                                                                            //
//     This routine returns the approximation to Pr[F < f] by                 //
//     Pr[X < x] where x = v1 * f.                                            //
//                                                                            //
//  Arguments:                                                                //
//     double x   The upper limit of the integral of the density given above. //
//     int   v1   The number of degrees of freedom of the numerator of the    //
//                F-test, i.e. (v1/2 - 1) is the exponent of f in the         //
//                integrand above.  Note v1 >= 1.                             //
//     int   v2   The number of degrees of freedom of the denominator of the  //
//                F-test, i.e. (-(v1+v2)/2) is the exponent of the term       //
//                (v2 + v1 f) in the integrand above.  Note v2 is assumed     //
//                to be infinite and is not used.                             //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, x;                                                           //
//     int    v1, v2;                                                         //
//                                                                            //
//     p = F_Distribution_Large_Denominator_dof(x, v1, v2);                   //
////////////////////////////////////////////////////////////////////////////////

double F_Distribution_Large_Denominator_dof(double f, int v1, int v2)
{
   if ( f <= 0.0 ) return 0.0;

   return Chi_Square_Distribution( (double) v1 * f, v1 );
}
