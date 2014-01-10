////////////////////////////////////////////////////////////////////////////////
// File: hypergeometric_dist_tables.c                                         //
// Routine(s):                                                                //
//    Hypergeometric_Distribution_Tables                                      //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Binomial_Coefficient( int n, int m );

////////////////////////////////////////////////////////////////////////////////
// void Hypergeometric_Distribution_Tables( int n1, int n2, int n,            //
//                                          double* pr, double* cumulative )  //
//                                                                            //
//  Description:                                                              //
//     Given n1 objects labeled 1 and n2 objects labeled 0.  The probability  //
//     by drawing n objects without replacement that k 1's are drawn and      //
//     n - k 0's are drawn is:                                                //
//                    C(n1,k) C(n2,n-k) / C(n1 + n2,n)                        //
//     where C(r,s) = r! / (s! (r-s)!).                                       //
//     Let X be the number of 1's that are drawn. The distribution of X       //
//     is called the hypergeometric distribution.                             //
//     The cumulative distribution of X is Pr[X <= k],                        //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 0,...,k <= n and                             //
//     Pr[X = j] = C(n1,j) C(n2,n-j) / C(n1 + n2,n) if j <= n1 & n-j <= n2    //
//               = 0                                elsewhere                 //
//                                                                            //
//     This routine returns Pr[X = i] in the user supplied array pr[], i = 0, //
//     ..., n and returns Pr[X <= i] in the user supplied array cumulative[], //
//     i = 0,..., n.                                                          //
//                                                                            //
//  Arguments:                                                                //
//     int    n1             The number of objects labeled 1.                 //
//     int    n2             The number of objects labeled 0.                 //
//     int    n              The number of objects drawn without replacement. //
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
//     #define SIZE n+1                                                       //
//     double pr[SIZE], cumulative[SIZE];                                     //
//     int n1, n2, n;                                                         //
//                                                                            //
//     Hypergeometric_Distribution_Tables(int n1, int n2, int n, pr,          //
//                                                              cumulative);  //
////////////////////////////////////////////////////////////////////////////////
void Hypergeometric_Distribution_Tables( int n1, int n2, int n, double* pr, 
                                                          double* cumulative )
{
   double summand;
   int i;
   int k1 = (n <= n2) ? 0 : n - n2;
   int k2 = (n <= n1) ? n : n1;

   for (i = 0; i < k1; i++) pr[i] = 0.0;
   if ( k1 <= k2 )
      pr[k1] = Binomial_Coefficient(n1,k1) * Binomial_Coefficient(n2,n-k1) 
                                              / Binomial_Coefficient(n1+n2,n);
   for (i = k1+1; i <= k2; i++)  
      pr[i] = pr[i-1] * (double)(n1 - i + 1) * (double)(n - i + 1)
                                  / ((double)i * (double)(n2 - n + i));
   for (i = k2+1; i <= n; i++) pr[i] = 0.0;
  
   cumulative[0] = pr[0];
   for (i = 1; i <= n; i++) cumulative[i] = cumulative[i-1] + pr[i];
}
