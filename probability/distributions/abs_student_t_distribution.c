////////////////////////////////////////////////////////////////////////////////
// File: abs_student_t_distribution.c                                         //
// Routine(s):                                                                //
//    Absolute_Student_t_Distribution                                         //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
extern double Beta_Distribution(double x, double a, double b);

////////////////////////////////////////////////////////////////////////////////
// double Absolute_Student_t_Distribution( double x, int n )                  //
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
//     p = Absolute_Student_t_Distribution(x, n);                             //
////////////////////////////////////////////////////////////////////////////////

double Absolute_Student_t_Distribution(double x, int n)
{
   double a = (double)n / 2.0;
   double beta = Beta_Distribution( 1.0 / (1.0 + x * x / n), a, 0.5);

   return 1.0 - beta;
}
