////////////////////////////////////////////////////////////////////////////////
// File: hypergeometric_point_distribution.c                                  //
// Routine(s):                                                                //
//    Hypergeometric_Point_Distribution                                       //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //
double Binomial_Coefficient( int n, int m );

////////////////////////////////////////////////////////////////////////////////
// double Hypergeometric_Point_Distribution( int n1, int n2, int n, int k )   //
//                                                                            //
//  Description:                                                              //
//     Given n1 objects labeled 1 and n2 objects labeled 0.  The probability  //
//     by drawing n objects without replacement that k 1's are drawn and      //
//     n - k 0's are drawn is:                                                //
//                    C(n1,k) C(n2,n-k) / C(n1 + n2,n)                        //
//     where C(r,s) = r! / (s! (r-s)!).                                       //
//     Let X be the number of 1's that are drawn. The distribution of X       //
//     is called the hypergeometric distribution.                             //
//     The point distribution of X is                                         //
//     Pr[X = k] = C(n1,k) C(n2,n-k) / C(n1 + n2,n) if k <= n1 & n-k <= n2    //
//               = 0                                elsewhere.                 //
//                                                                            //
//  Arguments:                                                                //
//     int    n1  The number of objects labeled 1.                            //
//     int    n2  The number of objects labeled 0.                            //
//     int    n   The number of objects drawn without replacement,            //
//                n <= n1 + n2.                                               //
//     int    k   The number of 1's drawn, n-k is the number of 0's drawn.    //
//                                                                            //
//  Return Values:                                                            //
//     A real number between 0 and 1.                                         //
//                                                                            //
//  Example:                                                                  //
//     double pr;                                                             //
//     int n1,n2,n,k;                                                         //
//                                                                            //
//     pr = Hypergeometric_Point_Distribution(n1, n2, n, k);                  //
////////////////////////////////////////////////////////////////////////////////

double Hypergeometric_Point_Distribution( int n1, int n2, int n, int k )
{
   if ( k < 0 ) return 0.0;
   if ( n > n2 )
      if ( k < n - n2 ) return 0.0;
   if ( k > n ) return 0.0;
   if ( n1 < n )
      if ( k > n1 ) return 0.0;
   
   return Binomial_Coefficient(n1,k) * Binomial_Coefficient(n2,n-k) 
                                          / Binomial_Coefficient(n1 + n2, n);
}
