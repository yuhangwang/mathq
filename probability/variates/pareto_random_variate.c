////////////////////////////////////////////////////////////////////////////////
// File: pareto_random_variate.c                                              //
// Routine(s):                                                                //
//    Pareto_Random_Variate                                                   //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines                    //

double Gamma_Random_Variate( double a );
double Exponential_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Pareto_Random_Variate( double a )                                   //
//                                                                            //
//  Description:                                                              //
//     This function returns a Pareto random variate: Let e,g be independent  //
//     random variables e an exponentially distribued random variable and g a //
//     gamma(a) random variable, then x = 1 + e / g has a Pareto distribution //
//     with probability density                                               //
//                          f(x) = a / x^(a+1) for x >=1, a >0                //
//     and distribution                                                       //
//                          F(x) = 1 - 1 / x^a.                               //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Pareto distributed random quantity.                                  //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Pareto_Random_Variate();                                           //
////////////////////////////////////////////////////////////////////////////////

double Pareto_Random_Variate( double a )
{
   double e = Exponential_Random_Variate();
   double g = Gamma_Random_Variate(a);

   return 1.0 + e / g;
}
