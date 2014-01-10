////////////////////////////////////////////////////////////////////////////////
// File: gamma_distribution.c                                                 //
// Routine(s):                                                                //
//    Gamma_Distribution                                                      //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //

extern double Entire_Incomplete_Gamma_Function(double x, double nu);


////////////////////////////////////////////////////////////////////////////////
// double Gamma_Distribution(double x, double nu)                             //
//                                                                            //
//  Description:                                                              //
//     The gamma distribution is defined to be the integral of the gamma      //
//     density which is 0 for x < 0 and x^(nu-1) exp(-nu) for x > 0 where the //
//     parameter nu > 0. The parameter nu is referred to as the shape         //
//     parameter.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double x   Upper limit of the integral of the density given above.     //
//     double nu  The shape parameter of the gamma distribution.              //
//                                                                            //
//  Return Values:                                                            //
//                                                                            //
//  Example:                                                                  //
//     double x, g, nu;                                                       //
//                                                                            //
//     g = Gamma_Distribution( x, nu );                                       //
////////////////////////////////////////////////////////////////////////////////


double Gamma_Distribution(double x, double nu) {
   return  ( x <= 0.0 ) ? 0.0 : Entire_Incomplete_Gamma_Function(x,nu);
}
