////////////////////////////////////////////////////////////////////////////////
// File: log_series_dist_tables.c                                             //
// Routine(s):                                                                //
//    Log_Series_Distribution_Tables                                          //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                          // required for log()

////////////////////////////////////////////////////////////////////////////////
// void Log_Series_Distribution_Tables( int size, double p,                   //
//                                          double* pr, double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     This routine returns Pr[X = i + 1] in the user supplied array pr[],    //
//     i = 0, ..., size-1 and returns Pr[X <= i + 1] in the user supplied     //
//     array cumulative[], i = 0,..., size-1.                                 //
//                                                                            //
//  Arguments:                                                                //
//     int    size           The total number of trials, n >= 1.              //
//     double p              The shape parameter of a logarithmic             //
//                           distribution, 0 < p < 1.                         //
//     double pr[]           pr[i] = Pr[X = i], i = 0,...,size-1.             //
//     double cumulative[]   cumulative[i] = Pr[X <= i], i = 0,...,size-1.    //
//                                                                            //
//  Return Values:                                                            //
//     None                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define SIZE                                                           //
//     double p, pr;                                                          //
//     double pr[SIZE], cumulative[SIZE];                                     //
//                                                                            //
//     Log_Series_Distribution_Tables( SIZE, p, pr, cumulative);              //
////////////////////////////////////////////////////////////////////////////////
void Log_Series_Distribution_Tables( int size, double p, double* pr,
                                                           double* cumulative )
{
   int i;
   
   if ( p == 0.0 )
      for (i = 0; i < size; i++) pr[i] = 0.0;
   else {
      pr[0] = - p / log(1.0 - p);
      for (i = 1; i < size; i++) pr[i] = pr[i-1] * p;
      for (i = 1; i < size; i++) pr[i] /= (double)(i + 1);
   }
   cumulative[0] = pr[0];
   for (i = 1; i < size; i++) cumulative[i] = cumulative[i-1] + pr[i];
};
