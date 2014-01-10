////////////////////////////////////////////////////////////////////////////////
// File: uniform_0_1_distribution.c                                           //
// Routine(s):                                                                //
//    Uniform_0_1_Distribution                                                //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// double Uniform_0_1_Distribution( double x )                                //
//                                                                            //
//  Description:                                                              //
//     This function returns the probability that a random variable with      //
//     a uniform distribution on the interval [0,1] a value less than "x".    //
//                                                                            //
//  Arguments:                                                                //
//     double x   Argument of Pr[X < x] where X ~ U(0,1).                     //
//                                                                            //
//  Return Values:                                                            //
//     The probability of observing a value less than (or equal) to x assuming//
//     a uniform distribution on the interval [0,1].                          //
//                                                                            //
//  Example:                                                                  //
//     double x, pr;                                                          //
//                                                                            //
//     pr = Uniform_0_1_Distribution(x);                                      //
////////////////////////////////////////////////////////////////////////////////

double Uniform_0_1_Distribution( double x )
{
   if (x < 0.0) return 0.0;
   if (x > 1.0) return 1.0;
   return  x;
}
