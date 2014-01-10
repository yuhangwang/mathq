////////////////////////////////////////////////////////////////////////////////
// File: laplace_random_variate.c                                             //
// Routine(s):                                                                //
//    Laplace_Random_Variate                                                  //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
double Exponential_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Laplace_Random_Variate( void )                                      //
//                                                                            //
//  Description:                                                              //
//     This function returns a Laplace random variate: Let u,v be independent //
//     random variables, u a U(0,1) random variable and v an exponentially    //
//     distribued random variable, then x = -v if u < 1/2, x = 0 if u = 1/2,  //
//     and x = v if u > 1/2, has a Laplace distribution with probability      //
//     density                                                                //
//                          f(x) = (1/2) exp(- |x|)                           //
//     and distribution                                                       //
//                          F(x) = (1/2) exp(x)     , x < 0                   //
//                                 1 - (1/2) exp(-x), x >= 0.                 //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Laplace distributed random quantity.                                 //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Laplace_Random_Variate();                                          //
////////////////////////////////////////////////////////////////////////////////

double Laplace_Random_Variate( void )
{
   double u = Uniform_0_1_Random_Variate();
   double e = Exponential_Random_Variate();

   return (u < 0.5) ? -e : (u > 0.5) ? e : 0.0;
}
