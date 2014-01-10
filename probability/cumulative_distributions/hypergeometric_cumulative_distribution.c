////////////////////////////////////////////////////////////////////////////////
// File: hypergeometric_cumulative_distribution.c                             //
// Routine(s):                                                                //
//    Hypergeometric_Cumulative_Distribution                                  //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Binomial_Coefficient( int n, int m );

////////////////////////////////////////////////////////////////////////////////
// double Hypergeometric_Cumulative_Distribution( int n1, int n2, int n,      //
//                                                                   int k )  //
//                                                                            //
//  Description:                                                              //
//     Given n1 objects labeled 1 and n2 objects labeled 0.  The probability  //
//     by drawing n objects without replacement that k 1's are drawn and      //
//     n - k 0's are drawn is:                                                //
//                    C(n1,k) C(n2,n-k) / C(n1 + n2,n)                        //
//     where C(r,s) = r! / (s! (r-s)!).                                       //
//     Let X be the number of 1's that are drawn. The distribution of X       //
//     is called the hypergeometric distribution.                             //
//     The cumulative distribution function of X is Pr[X <= k],               //
//                        Pr[X <= k] = Sum Pr[X = j]                          //
//     where the sum is over j = 0,...,k <= n and                             //
//     Pr[X = j] = C(n1,j) C(n2,n-j) / C(n1 + n2,n) if j <= n1 & n-j <= n2    //
//               = 0                                elsewhere.                 //
//                                                                            //
//  Arguments:                                                                //
//     int    n1  The number of objects labeled 1.                            //
//     int    n2  The number of objects labeled 0.                            //
//     int    n   The number of objects drawn without replacement.            //
//                n <= n1 + n2.                                               //
//     int    k   The maximum number of 1's drawn, n-k is the minimum number  //
//                of 0's drawn.                                               //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double p, pr;                                                          //
//     int n1, n2, n, k;                                                      //
//                                                                            //
//     pr = Hypergeometric_Cumulative_Distribution(n1, n2, n, k);             //
////////////////////////////////////////////////////////////////////////////////
double Hypergeometric_Cumulative_Distribution( int n1, int n2, int n, int k )
{
   double pr = 0.0;
   double summand;
   int i;
   int k1 = (n <= n2) ? 0 : n - n2;
   int k2 = (n <= n1) ? n : n1;

   if ( k < k1 ) return 0.0;
   if ( k >= k2 ) return 1.0;
   if ( k2 < k1 ) return 0.0;

   summand = Binomial_Coefficient(n1,k1) * Binomial_Coefficient(n2,n-k1) 
                                             / Binomial_Coefficient(n1+n2,n);
   pr = summand;
   for (i = k1; i < k; i++) {
      summand *= (double)(n1 - i) * (double)(n - i)
                                / ((double)(i + 1) * (double)(n2 - n + i + 1));
      pr += summand;
   }

   return pr;
}
