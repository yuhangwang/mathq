////////////////////////////////////////////////////////////////////////////////
// File: geometric_dist_tables.c                                              //
// Routine(s):                                                                //
//    Geometric_Distribution_Tables                                           //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// void Geometric_Distribution_Tables( int size, double p, double* pr,        //
//                                                      double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     Let X have the geometric distribution Geom(p).                         //
//     This routine returns Pr[X = i] in the user supplied array pr[], i = 0, //
//     ..., size - 1 and returns Pr[X <= i] in the user supplied array        //
//     cumulative[], i = 0,..., size-1.                                       //
//                                                                            //
//  Arguments:                                                                //
//     int    size           The max argument - 1 of Pr[]. size >= 1.         //
//                           The dimension of the arrays pr and cumulative.   //
//     double p              The probability of a success (1) on a single     //
//                           trial. 0 <= p <= 1.                              //
//     double pr[]           pr[i] = Pr[X = i], i = 0,...,size - 1.           //
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
//     double p, pr;                                                          //
//     double pr[SIZE], cumulative[SIZE];                                     //
//     int n = N;                                                             //
//                                                                            //
//     Geometric_Distribution_Tables(SIZE, p, pr, cumulative);                //
////////////////////////////////////////////////////////////////////////////////
void Geometric_Distribution_Tables( int size, double p, double* pr, 
                                                          double* cumulative )
{
   double q;
   int i;
   
   if ( p == 0.0 ) {
      for (i = 0; i < size; i++) pr[i] = 0.0;
   }
   else if ( p == 1.0 ) {
      pr[0] = 1.0;
      for (i = 1; i < size; i++) pr[i] = 0.0;
   }
   else {
      pr[0] = p;
      q = 1.0 - p;
      for (i = 1; i < size; i++) pr[i] = pr[i-1] * q;
   }
   cumulative[0] = pr[0];
   for (i = 1; i < size; i++) cumulative[i] = cumulative[i-1] + pr[i];
}
