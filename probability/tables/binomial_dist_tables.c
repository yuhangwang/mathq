////////////////////////////////////////////////////////////////////////////////
// File: binomial_dist_tables.c                                               //
// Routine(s):                                                                //
//    Binomial_Distribution_Tables                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp(), log()

////////////////////////////////////////////////////////////////////////////////
// void Binomial_Distribution_Tables( int n, double p, double* pr,            //
//                                                      double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     Let X have the binomial distribution Binomial(n,p).                    //
//     This routine returns Pr[X = i] in the user supplied array pr[], i = 0, //
//     ..., n and returns Pr[X <= i] in the user supplied array cumulative[], //
//     i = 0,..., n.                                                          //
//                                                                            //
//  Arguments:                                                                //
//     int    n              The total number of trials, n >= 1.              //
//     double p              The probability of a success (1) on a single     //
//                           trial. 0 <= p <= 1.                              //
//     double pr[]           pr[i] = Pr[X = i], i = 0,...,n.                  //
//                           Note: pr must be defined in the calling routine  //
//                           as double pr[N], where N >= n+1.                 //
//     double cumulative[]   cumulative[i] = Pr[X <= i], i = 0,...,n.         //
//                           Note: cumulative must be defined in the calling  //
//                           routine as double cumulative[N], where N >= n+1. //
//                                                                            //
//  Return Values:                                                            //
//     None                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double p, pr;                                                          //
//     double pr[N+1], cumulative[N+1];                                       //
//     int n = N;                                                             //
//                                                                            //
//     Binomial_Distribution_Tables(n, p, pr, cumulative);                    //
////////////////////////////////////////////////////////////////////////////////
void Binomial_Distribution_Tables( int n, double p, double* pr,
                                                           double* cumulative )
{
   double q;
   int i;
   
   if ( p == 0.0 ) {
      pr[0] = 1.0;
      for (i = 1; i <= n; i++) pr[i] = 0.0;
   }
   else if ( p == 1.0 ) {
      pr[n] = 1.0;
      for (i = 0; i < n; i++) pr[i] = 0.0;
   }
   else {
      pr[0] = exp((double)n * log(1.0 - p));
      q = p / (1.0 - p);
      for (i = 0; i < n; i++) pr[i+1] = pr[i] * q * (double) (n-i) / (double)(i+1);
   }
   cumulative[0] = pr[0];
   for (i = 1; i <= n; i++) cumulative[i] = cumulative[i-1] + pr[i];
}
