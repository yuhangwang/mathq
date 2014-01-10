////////////////////////////////////////////////////////////////////////////////
// File: gaussian_random_variate.c                                            //
// Routine(s):                                                                //
//    Init_Gaussian_Random_Variate                                            //
//    Gaussian_Random_Variate                                                 //
////////////////////////////////////////////////////////////////////////////////

static double (*rn)(void);

////////////////////////////////////////////////////////////////////////////////
// void Init_Gaussian_Random_Variate( double (*r_generator)(void) )           //
//                                                                            //
//  Description:                                                              //
//     This function saves the pointer to a Gaussian random number            //
//     generator.  Subsequent calls to Gaussian_Random_Variate (below)        //
//     will call this routine.                                                //
//                                                                            //
//  Arguments:                                                                //
//     double (*r_generator)(void)                                            //
//        The user supplied random number generator generating an             //
//        gaussian distributed random number on [-inf,inf].                   //
//                                                                            //
//  Return Values:                                                            //
//     None                                                                   //
//                                                                            //
//  Example:                                                                  //
//     extern double (*gaussrand)(void);                                      //
//     Init_Gaussian_Random_Variate( gaussrand );                             //
////////////////////////////////////////////////////////////////////////////////
        
void Init_Gaussian_Random_Variate( double (*r_generator)(void) )
{
   rn = r_generator;
}


////////////////////////////////////////////////////////////////////////////////
// double Gaussian_Random_Variate( void )                                     //
//                                                                            //
//  Description:                                                              //
//     This function returns an Gaussian distributed random variate on        //
//     [-DBL_MAX, DBL_MAX]                                                    //
//                                                                            //
//  Arguments:                                                                //
//     None                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     A random number with a Gaussian distribution between -DBL_MAX and      //
//     DBL_MAX.                                                               //
//                                                                            //
//  Example:                                                                  //
//     double x;                                                              //
//                                                                            //
//     x = Gaussian_Random_Variate();                                         //
////////////////////////////////////////////////////////////////////////////////
        
double Gaussian_Random_Variate( void ) { return (*rn)(); }
