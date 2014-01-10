////////////////////////////////////////////////////////////////////////////////
// File: catalan_beta_function.c                                              //
// Routine(s):                                                                //
//    Catalan_Beta_Function                                                   //
//    xCatalan_Beta_Function                                                  //
//    Catalan_Beta_Star_Function                                              //
//    xCatalan_Beta_Star_Function                                             //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for fabsl(), powl(), cosl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double      Catalan_Beta_Function( double s);
long double xCatalan_Beta_Function( long double s);
double      Catalan_Beta_Star_Function( double s);
long double xCatalan_Beta_Star_Function( long double s);

static long double Reflection_Coefficient(long double s);
static long double Alternating_Series_Convergence_Acceleration(long double s);
static long double Sum_Reverse_Order(long double s);

//                         Externally Defined Routines                        //

extern long double xGamma_Function(long double x);


//                         Internally Defined Constants                       //

static const long double pi_2 = 1.57079632679489661923132169163975144L; // pi/2

////////////////////////////////////////////////////////////////////////////////
// double Catalan_Beta_Function(double s)                                     //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan Beta function for real s,          //
//     where for s > 1,                                                       //
//        beta(s) = Sum (-1)^k (1/(2k+1)^s) where the sum extends over        //
//     k = 0,... .                                                            //
//                                                                            //
//     Note that for s < 1, the reflection formula is used.  For s << 0, the  //
//     slope of beta(s) has a large magnitude and small perturbations in the  //
//     argument can create a large absolute error.                            //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the beta function.                                  //
//                                                                            //
//  Return Value:                                                             //
//     beta(s), if beta(s) > DBL_MAX, then DBL_MAX is returned,               //
//              if beta(s) < -DBL_MAX, then -DBL_MAX is returned.             //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double beta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the beta function. )             //
//                                                                            //
//     beta = Catalan_Beta_Function(s);                                       //
////////////////////////////////////////////////////////////////////////////////

double Catalan_Beta_Function(double s)
{
   long double x;

   x = xCatalan_Beta_Function( (long double) s);

   if (fabsl(x) < DBL_MAX) return (double) x;
   else return ( x > 0.0L ) ? DBL_MAX : - DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xCatalan_Beta_Function(long double s)                          //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan Beta function for real s,          //
//     where for s > 1,                                                       //
//        beta(s) = Sum (-1)^k (1/(2k+1)^s) where the sum extends over        //
//     k = 0,... .                                                            //
//                                                                            //
//     Note that for s < 1, the reflection formula is used.  For s << 0, the  //
//     slope of beta(s) has a large magnitude and small perturbations in the  //
//     argument can create a large absolute error.                            //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the beta function.                                  //
//                                                                            //
//  Return Value:                                                             //
//     beta(s)                                                                //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double beta;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the beta function. )             //
//                                                                            //
//     beta = xCatalan_Beta_Function(s);                                      //
////////////////////////////////////////////////////////////////////////////////

long double xCatalan_Beta_Function(long double s)
{
   if (s >= 40.0L) return 1.0L;
   return 1.0L + xCatalan_Beta_Star_Function(s);
}


////////////////////////////////////////////////////////////////////////////////
// double Catalan_Beta_Star_Function(double s)                                //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan beta star function defined as      //
//                          beta*(s) = beta(s) - 1,                           //
//     where beta() is the Catalan beta function.                             //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the beta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     beta*(s), if beta*(s) > DBL_MAX, then DBL_MAX is returned,             //
//               if beta*(s) < -DBL_MAX, then -DBL_MAX is returned.           //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double beta_star;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the beta* function. )            //
//                                                                            //
//     beta_star = Catalan_Beta_Star_Function(s);                             //
////////////////////////////////////////////////////////////////////////////////

double Catalan_Beta_Star_Function(double s)
{
   long double x;

   x = xCatalan_Beta_Star_Function((long double)s);

   if (fabsl(x) < DBL_MAX) return (double) x;
   else return ( x > 0.0L ) ? DBL_MAX : - DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xCatalan_Beta_Star_Function(long double s)                     //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan beta star function defined as      //
//                          beta*(s) = beta(s) - 1,                           //
//     where beta() is the Catalan beta function.                             //
//                                                                            //
//  Arguments:                                                                //
//     long double s;                                                         //
//        The argument of the beta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     beta*(s)                                                               //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double beta_star;                                                 //
//                                                                            //
//     ( User code to set s, the argument to the beta* function. )            //
//                                                                            //
//     beta_star = xCatalan_Beta_Star_Function(s);                            //
////////////////////////////////////////////////////////////////////////////////

long double xCatalan_Beta_Star_Function(long double s)
{
   long double reflection_coefficient;

   if (s >= 0.0L)
      if ( s < 18.0L ) 
         return Alternating_Series_Convergence_Acceleration(s);
      else
         return Sum_Reverse_Order(s);

                       // Use the reflection formula //

   reflection_coefficient = Reflection_Coefficient(s);    
   return reflection_coefficient * xCatalan_Beta_Star_Function(1.0L - s)
          + (reflection_coefficient - 1.0L);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Reflection_Coefficient(long double s)                   //
//                                                                            //
//  Description:                                                              //
//     The reflection formula for the Catalan beta function is:               //
//     beta(s) = { [pi / 2]^(s-1) * Gamma(1-s) * cos(s*pi/2) } *  beta(1-s)   //
//     where s < 0 is a negative real number.                                 //
//                                                                            //
//     This routine returns the coefficient:                                  //
//               { [pi / 2]^(s-1) * Gamma(1-s) * cos(s*pi/2) }                //
//     where s < 0 is a negative real number.                                 //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the beta function, where s < 0.                     //
//                                                                            //
//  Return Value:                                                             //
//     The coefficient c so that beta(s) = c * beta(1-s).                     //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double coef;                                                      //
//                                                                            //
//     ( User code to set s < 0, the argument to the beta function. )         //
//                                                                            //
//     coef = Reflection_Coefficient(s);                                      //
////////////////////////////////////////////////////////////////////////////////

static long double Reflection_Coefficient(long double s)
{
   int k = (int) (s / 4.0L);
   long double v = s - 4.0L * (long double) k;
   long double x = cosl(v * pi_2);
   long double one_s;

   if (fabsl(x) < 1.8L * LDBL_EPSILON) return 0.0L;
   one_s = 1.0L - s;
   x *= xGamma_Function(one_s);
   x /= powl(pi_2,one_s);
   return x;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Alternating_Series_Convergence_Acceleration             //
//                                                            (long double s) //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Beta* function for s <= 18 where the series//
//     converges slowly.                                                      //
//                                                                            //
//     Technique:                                                             //
//     Given an alternating series: Sum (-1)^k a[k] summmed over k = 0,...    //
//     where a[k] > 0 for all k.  If there exits a positive function w(x)     //
//     defined on [0,1] such that a[k] = I[0,1] x^k w(x) dx (the integral     //
//     from 0 to 1 with respect to x of x^k w(x), then                        //
//          Sum (-1)^k a[k] = Sum (-1)^k I[0,1] x^k w(x) dx                   //
//                          = I[0,1] Sum (-1)^k x^k w(x) dx                   //
//                          = I[0,1] w(x) / (1 + x) dx                        //
//                                                                            //
//     Let S = I[0,1] w(x) / (1 + x) dx.                                      //
//     Given a polynomial P(x) of degree n such that P(-1) != 0,              //
//          P(x) = p[0] - p[1] x + p[2] x^2 - ... (-1)^n p[n] x^n,            //
//     define S[n] = I[0,1] (P(-1) - P(x)) w(x) / (1+x) dx / P(-1).           //
//     Then   S[n] = S - I[0,1] P(x) w(x) / (1+x) dx / P(-1)                  //
//     so,   |S[n] - S| <= I[0,1] |P(x)| w(x) / (1+x) dx / |P(x)|             //
//           |S[n] - S| <= (M / |P(-1)|) S, where M = max |P(x)|, x in [0,1]. //
//                                                                            //
//     After a little algebra,                                                //
//               S[n] = (1/P(-1)) Sum Sum p[j] [(-1)^k a[k]],                 //
//     where the first sum is over k = 0,...,n-1 and the second sum is over   //
//     j = k+1,...,n.                                                         //
//     Let d[k] = Sum p[j] summed for j = k+1,...,n for  k = n-1, ..., 0, then//
//        S[n] = (1/c[0]) Sum d[k] (-1)^k a[k] summed k = 0,...,n-1.          //
//                                                                            //
//     The choice for P(x) programmed here is (-1)^n * T*[n](x), where T*[n]  //
//     is the shifted Chebyshev polynomial of the first kind of degree n.     //
//     With this choice of P(x),                                              //
//              P(x) = Sum {( n / (n+j) ) C(n+j,2j) 4^j (-x)^j},              //
//     where the sum is over j = 0,...,n and C(n+j,2j) = (n+j)!/[(2j)!(n-j)!] //
//     and M = 1, P(-1) >= (1/2) (3 + sqrt(8))^n.                             //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the beta* function.                                 //
//                                                                            //
//  Return Value:                                                             //
//     beta*(s)                                                               //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double beta_star;                                                 //
//                                                                            //
//     ( User code to set s >= 18, the argument to the beta* function.)       //
//                                                                            //
//     beta_star = Alternating_Series_Convergence_Acceleration(s);            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static long double Alternating_Series_Convergence_Acceleration(long double s)
{
   static const long double d[29] = {
     1.362725501650887306817e+21L, 1.362725501650887306816e+21L,
     1.362725501650887305248e+21L, 1.362725501650886896000e+21L,
     1.362725501650844334208e+21L, 1.362725501648488235008e+21L,
     1.362725501568066715648e+21L, 1.362725499718371770368e+21L,
     1.362725469310199922688e+21L, 1.362725096810094788608e+21L,
     1.362721590926752350208e+21L, 1.362695647390018306048e+21L,
     1.362542007743905005568e+21L, 1.361803869444099801088e+21L,
     1.358896740140251611136e+21L, 1.349437033675348770816e+21L,
     1.323863206542645919744e+21L, 1.266218975223368122368e+21L,
     1.157712186857668739072e+21L, 9.872015194258554224640e+20L,
     7.640581139674368573440e+20L, 5.220333434317674905600e+20L,
     3.061506212814840135680e+20L, 1.496014168469232680960e+20L,
     5.884825485587356057600e+19L, 1.781624012587768217600e+19L,
     3.882102878793367552000e+18L, 5.404319552844595200000e+17L,
     3.602879701896396800000e+16L
};
   long double term[28];
   long double sum;
   int k;

   for (k = 1; k <= 28; k++) {
      term[k-1] = d[k] * powl((long double)(k+k+1),-s);
      k++;
      term[k-1] = -d[k] * powl((long double)(k+k+1),-s);
   }
   sum = term[27];
   for (k = 26; k >= 0; k--) sum += term[k];
   sum /= d[0];
  
   return -sum;
}


////////////////////////////////////////////////////////////////////////////////
// long double Sum_Reverse_Order(long double s)                               //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Catalan Beta* function for real s,         //
//     where for s > 18 (this requires fewer than 28 calls to powl()),        //
//        beta(s) = Sum (-1)^k (1/(2k+1)^s) where the sum extends over        //
//     k = 1,... .                                                            //
//                                                                            //
//     Since the series is alternating, the absolute error is less than the   //
//     absolute value of the first neglected term.  The number of terms which //
//     are summed is determined by calculating successive lower and upper     //
//     bounds until the two bounds are equal (within truncation errors).      //
//     The final result is then summed starting with the last term and        //
//     terminating with the first term, (this is slightly more accurate).     //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the beta* function.                                 //
//                                                                            //
//  Return Value:                                                             //
//     beta(s)                                                                //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double beta_star;                                                 //
//                                                                            //
//     ( User code to set s >= 18, the argument to the beta* function. )      //
//                                                                            //
//     beta = Sum_Reverse_Order(s);                                           //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static long double Sum_Reverse_Order(long double s)
{
   long double term[30];
   long double lower_bound;
   long double upper_bound;
   long double sum;
   int k;

   lower_bound = -powl(3.0L,-s);
   term[0] = lower_bound;
   term[1] = powl(5.0L,-s);
   upper_bound = term[1] + lower_bound;
   if ( lower_bound == upper_bound) return upper_bound;
   for (k = 3; k < 31; k++) {
      term[k-1] = -powl( (long double) (k + k + 1), -s);
      lower_bound = upper_bound + term[k-1];
      if (lower_bound == upper_bound) break;
      k++;
      term[k-1] = powl( (long double) (k + k + 1), -s);
      upper_bound = lower_bound + term[k-1];
      if (lower_bound == upper_bound) break;
   }
   k -= 2;
   sum = term[k];
   for (k--; k >= 0; k--) sum += term[k];
   return sum;
}
