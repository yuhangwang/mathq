////////////////////////////////////////////////////////////////////////////////
// File: negative_binomial_dist_tables.c                                      //
// Routine(s):                                                                //
//    Negative_Binomial_Distribution_Tables                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp(), log()

////////////////////////////////////////////////////////////////////////////////
// void Negative_Binomial_Distribution_Tables( int n, int size, double p,     //
//                                          double* pr, double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     Let X have the negative binomial distribution NegBinomial(n,p).        //
//     This routine returns Pr[X = i] in the user supplied array pr[], i = 0, //
//     ..., n and returns Pr[X <= i] in the user supplied array cumulative[], //
//     i = 0,..., size-1.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     int    n              The number of 1's to occur before stopping,      //
//                           n >= 1.                                          //
//     int    size           The total number of trials, n >= 1.              //
//     double p              The probability of a success (1) on a single     //
//                           trial. 0 <= p <= 1.                              //
//     double pr[]           pr[i] = Pr[X = i], i = 0,...,size - 1.           //
//     double cumulative[]   cumulative[i] = Pr[X <= i], i = 0,..., size - 1. //
//                                                                            //
//  Return Values:                                                            //
//     None                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define SIZE                                                           //
//     double p, pr;                                                          //
//     double pr[SIZE], cumulative[SIZE];                                     //
//     int n = N;                                                             //
//                                                                            //
//     Negative_Binomial_Distribution_Tables(n, SIZE, p, pr, cumulative);     //
////////////////////////////////////////////////////////////////////////////////
void Negative_Binomial_Distribution_Tables( int n, int size,double p, 
                                                double* pr, double* cumulative )
{
   double q;
   int i;
   
   if ( p == 0.0 ) {
      pr[0] = 1.0;
      for (i = 1; i < size; i++) pr[i] = 0.0;
   }
   else if ( p == 1.0 ) {
      pr[size-1] = 1.0;
      for (i = 0; i < size - 1; i++) pr[i] = 0.0;
   }
   else {
      pr[0] = exp((double)n * log(p));
      q = 1.0 - p;
      for (i = 1; i < size; i++)
         pr[i] = pr[i-1] * q * (double) (n+i-1) / (double)(i);
   }
   cumulative[0] = pr[0];
   for (i = 1; i < size; i++) cumulative[i] = cumulative[i-1] + pr[i];
}
