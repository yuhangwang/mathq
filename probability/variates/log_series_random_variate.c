////////////////////////////////////////////////////////////////////////////////
// File: log_series_random_variate.c                                          //
// Routine(s):                                                                //
//    Log_Series_Random_Variate                                               //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>                               // required for log()

//                    Required Externally Defined Routines                    //

double Uniform_0_1_Random_Variate( void );
        
////////////////////////////////////////////////////////////////////////////////
// int Log_Series_Random_Variate( double p )                                  //
//                                                                            //
//  Description:                                                              //
//     This function returns a Log Series distributed random variate using the//
//     Kemp's procedure which is based on the observation that if U,V are     //
//     independent U(0,1) random variables, then                              //
//                   X = (int)[1 + ln V / ln(1 - (1-p)^U)]                    //
//     has a Log Series distribution with parameter p.   The log series       //
//     distribution with parameter p is:                                      //
//                      Pr[X = x] = p^x / (-x ln(1-p)).                       //
//                                                                            //
//  Arguments:                                                                //
//     double p                                                               //
//            The "shape" parameter, p = 2 Pr[X=2]/Pr[X=1], 0 < p < 1.        //
//                                                                            //
//  Return Values:                                                            //
//     A Log Series distributed random quantity.                              //
//                                                                            //
//  Example:                                                                  //
//     int x;                                                                 //
//     double p;                                                              //
//                                                                            //
//     x = Log_Series_Random_Variate(p);                                      //
////////////////////////////////////////////////////////////////////////////////

int Log_Series_Random_Variate( double p )
{
   double u = Uniform_0_1_Random_Variate();
   double y;
 
   if (u >= p) return 1;
   y = 1.0 - pow(1.0 - p, Uniform_0_1_Random_Variate());
   if (u > y) return 1;
   if (u > y*y) return 2;
   return (int)(1.0 + (log(u) / log(y)));
}
