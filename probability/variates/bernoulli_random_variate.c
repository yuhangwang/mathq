////////////////////////////////////////////////////////////////////////////////
// File: bernoulli_random_variate.c                                           //
// Routine(s):                                                                //
//    Bernoulli_Random_Variate                                                //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// int Bernoulli_Random_Variate( double p )                                   //
//                                                                            //
//  Description:                                                              //
//     This function returns a Bernoulli random variate which assumes the     //
//     value 1 with probability p and the value 0 with probability 1-p.       //
//                                                                            //
//  Arguments:                                                                //
//     double p   The probability that a 1 is returned, 0 <= p <= 1.          //
//                                                                            //
//  Return Values:                                                            //
//     Either a 0 or a 1.                                                     //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//     double p;                                                              //
//                                                                            //
//     x = Bernoulli_Random_Variate( p );                                     //
////////////////////////////////////////////////////////////////////////////////

int Bernoulli_Random_Variate(double p)
{
   return ( Uniform_0_1_Random_Variate() <= p ) ? 1 : 0;
}
