////////////////////////////////////////////////////////////////////////////////
// File: binomial_random_variate.c                                            //
// Routine(s):                                                                //
//    Binomial_Random_Variate                                                 //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>                    // required for log()

//                    Required Externally Defined Routines                    //

extern double Beta_Random_Variate( double a, double b );
extern double Exponential_Random_Variate( void );

//                    Required Internally Defined Routines                    //

static int Waiting_Time_Variate( int n, double p);

////////////////////////////////////////////////////////////////////////////////
// int Binomial_Random_Variate( int n, double p )                             //
//                                                                            //
//  Description:                                                              //
//     This function returns a Binomial(n,p) distributed random variate using //
//     the recursive algorithm: If Y~Binomial(n,p') and Z given Y has a       //
//     Binomial(n-Y,(p-p')/(1-p')) distribution, where p > p', then X = Y + Z //
//     has a Binomial(n,p) distribution and the waiting time algorithm        //
//     described below.                                                       //
//                                                                            //
//  Arguments:                                                                //
//     int    n                                                               //
//            The total number of trials, n >= 1.                             //
//     double p                                                               //
//            The probability of a success, 0 < p < 1.                        //
//                                                                            //
//  Return Values:                                                            //
//     A random number distributed as a Binomial(n,p) distribution.           //
//                                                                            //
//  Example:                                                                  //
//     int x, n;                                                              //
//     double p;                                                              //
//                                                                            //
//                  (* Set the probability p, 0 < p < 1 *)                    //
//                  (* Set the total number of trials n *)                    //
//                                                                            //
//     x = Binomial_Random_Variate( n, p );                                   //
////////////////////////////////////////////////////////////////////////////////
        
int Binomial_Random_Variate( int n, double p )
{
   double dp;
   int x = 0;
   int i = 0;

   while ( n * p >= 3) {
      i = (int)((n+1) * p);
      dp = Beta_Random_Variate((double)i,(double)(n-i+1));
      if (dp <= p) {
         x += i;
         n -= i;
         p = (p - dp) / (1.0 - dp);
      }
      else {
         n = i - 1;
         p /= dp;
      }
   }
   if (n == 0) return x;
   if (p <= 0.0) return x;
   x += Waiting_Time_Variate(n,p);
   return x;
}


////////////////////////////////////////////////////////////////////////////////
// static int Waiting_Time_Variate( int n, double p )                         //
//                                                                            //
//  Description:                                                              //
//     This function returns a Binomial(n,p) distributed random variate using //
//     the waiting time algorithm:  If Ei are independent identically         //
//     distributed exponentially distributed random variables, then the       //
//     random variable X = min{k: Sum(i=1,k+1)(Ei/(n-i+1) > -ln(1-p)} has a   //
//     binomial(n,p) distribution.                                            //
//                                                                            //
//  Arguments:                                                                //
//     int    n                                                               //
//            The total number of trials, n >= 1.                             //
//     double p                                                               //
//            The probability of a success, 0 < p < 1.                        //
//                                                                            //
//  Return Values:                                                            //
//     A random number distributed as a Binomial(n,p) distribution.           //
//                                                                            //
//  Example:                                                                  //
//     int x, n;                                                              //
//     double p;                                                              //
//                                                                            //
//                  (* Set the probability p, 0 < p < 1 *)                    //
//                  (* Set the total number of trials n *)                    //
//                                                                            //
//     x = Waiting_Time_Variate( n, p );                                      //
////////////////////////////////////////////////////////////////////////////////
static int Waiting_Time_Variate( int n, double p)
{
   double ln = -log(1.0 - p);
   double sum = 0.0;
   int x = 0;
   
   while (sum <= ln) {
      if (x < n)
         sum += (Exponential_Random_Variate() / (double)(n - x));
      else return n;
      x++;
   }
   return x - 1;
}
