////////////////////////////////////////////////////////////////////////////////
// File: kumaraswamys_random_variate.c                                        //
// Routine(s):                                                                //
//    Kumaraswamys_Random_Variate                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for pow()

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// double Kumaraswamys_Random_Variate( double a, double b )                   //
//                                                                            //
//  Description:                                                              //
//     This function returns a Kumaraswamy distributed variate: Let u be a    //
//     U(0,1) random variable then x = (1-u^(1/b))^(1/a) has a Kumaraswamy    //
//     distribution probability density                                       //
//                    f(x) = a b x^(a-1) (1 - x^a)^(b-1)                      //
//     and distribution                                                       //
//                      F(x) = 1 - (1 - x^a)^b.                               //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Kumaraswamy distributed random quantity.                             //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     double a;                                                              //
//     double b;                                                              //
//                                                                            //
//     x = Kumaraswamys_Random_Variate(a,b);                                  //
////////////////////////////////////////////////////////////////////////////////

double Kumaraswamys_Random_Variate( double a, double b )
{
   double u = Uniform_0_1_Random_Variate();

   return pow( 1.0 - pow(u, 1.0 / b), 1.0/a);
}
