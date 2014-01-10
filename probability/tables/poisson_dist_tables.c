////////////////////////////////////////////////////////////////////////////////
// File: poisson_dist_tables.c                                                //
// Routine(s):                                                                //
//    Poisson_Distribution_Tables                                             //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                             // required for exp()

////////////////////////////////////////////////////////////////////////////////
// void Poisson_Distribution_Tables( int size, double mu,                     //
//                                          double* pr, double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     Let X have a Poisson(mu) distribution.                                 //
//     This routine returns Pr[X = i] in the user supplied array pr[], i = 0, //
//     ..., n and returns Pr[X <= i] in the user supplied array cumulative[], //
//     i = 0,..., size-1.                                                     //
//                                                                            //
//  Arguments:                                                                //
//     int    size           The max argument - 1 of Pr[]. size >= 1.         //
//                           The dimension of the arrays pr and cumulative.   //
//     double mu             The probability of a success (1) on a single     //
//                           trial. 0 <= p <= 1.                              //
//     double pr[]           pr[i] = Pr[X = i], i = 0,...,size - 1/           //
//                           Note: pr must be defined in the calling routine  //
//                           as double pr[N], where N >= size.                //
//     double cumulative[]   cumulative[i] = Pr[X <= i], i = 0,...,size - 1.  //
//                           Note: cumulative must be defined in the calling  //
//                           routine as double cumulative[N], where N >= size.//
//                                                                            //
//  Return Values:                                                            //
//     None                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define SIZE                                                           //
//     double mu;                                                             //
//     double pr[SIZE], cumulative[SIZE];                                     //
//                                                                            //
//     Poisson_Distribution_Tables(SIZE, mu, pr, cumulative);                 //
////////////////////////////////////////////////////////////////////////////////
void Poisson_Distribution_Tables( int size, double mu, 
                                                double* pr, double* cumulative )
{
   int i;
   
   pr[0] = exp(- mu);
   for (i = 1; i < size; i++) pr[i] = pr[i-1] * mu / (double)(i);

   cumulative[0] = pr[0];
   for (i = 1; i < size; i++) cumulative[i] = cumulative[i-1] + pr[i];
}
