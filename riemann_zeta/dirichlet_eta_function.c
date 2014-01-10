////////////////////////////////////////////////////////////////////////////////
// File: dirichlet_eta_function.c                                             //
// Routine(s):                                                                //
//    Dirichlet_Eta_Function                                                  //
//    xDirichlet_Eta_Function                                                 //
//    Dirichlet_Eta_Star_Function                                             //
//    xDirichlet_Eta_Star_Function                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for fabsl(), powl(), cosl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double      Dirichlet_Eta_Function( double s);
long double xDirichlet_Eta_Function( long double s);
double      Dirichlet_Eta_Star_Function( double s);
long double xDirichlet_Eta_Star_Function( long double s);

static long double Reflection_Coefficient(long double s);
static long double Alternating_Series_Convergence_Acceleration(long double s);
static long double Sum_Reverse_Order(long double s);

//                         Externally Defined Routines                        //

extern long double xGamma_Function(long double x);


//                         Internally Defined Constants                       //

static const long double pi_2 = 1.57079632679489661923132169163975144L; // pi/2
static const long double pi = 3.14159265358979323846264338327950288L;   // pi


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Eta_Function(double s)                                    //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function for real s,                   //
//                      eta(s) = (1 - 2^(1-s)) zeta(s),                       //
//     where the zeta() is the Riemann zeta function.  For s > 1,             //
//                      eta(s) = Sum (-1)^(k-1) (1/k^s),                      //
//     where the sum is summed over k = 1,... . eta(1) = ln(2).               //
//                                                                            //
//     Note that for s < 1, the reflection formula is used.  For s << 0, the  //
//     slope of eta(s) has a large magnitude and small perturbations of the   //
//     argument can create a large absolute error.                            //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the eta function.                                   //
//                                                                            //
//  Return Value:                                                             //
//     eta(s), if eta(s) > DBL_MAX, then DBL_MAX is returned,                 //
//             if eta(s) < -DBL_MAX, then -DBL_MAX is returned.               //
//                                                                            //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double eta;                                                            //
//                                                                            //
//     ( User code to set s, the argument to the eta function. )              //
//                                                                            //
//     eta = Dirichlet_Eta_Function(s);                                       //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Eta_Function(double s)
{
   long double x;

   x = xDirichlet_Eta_Function( (long double) s);

   if (fabsl(x) < DBL_MAX) return (double) x;
   else return ( x > 0.0L ) ? DBL_MAX : - DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Eta_Function(long double s)                         //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function for real s,                   //
//                      eta(s) = (1 - 2^(1-s)) zeta(s),                       //
//     where the zeta() is the Riemann zeta function.  For s > 1,             //
//                      eta(s) = Sum (-1)^(k-1) (1/k^s),                      //
//     where the sum is summed over k = 1,... . eta(1) = ln(2).               //
//                                                                            //
//     Note that for s < 1, the reflection formula is used.  For s << 0, the  //
//     slope of eta(s) has a large magnitude and small perturbations of the   //
//     argument can create a large absolute error.                            //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the eta function.                                   //
//                                                                            //
//  Return Value:                                                             //
//     eta(s)                                                                 //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double eta;                                                       //
//                                                                            //
//     ( User code to set s, the argument to the eta function. )              //
//                                                                            //
//     eta = xDirichlet_Eta_Function(s);                                      //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Eta_Function(long double s)
{
   if (s >= 64.0L) return 1.0L;
   return 1.0L + xDirichlet_Eta_Star_Function(s);
}


////////////////////////////////////////////////////////////////////////////////
// double Dirichlet_Eta_Star_Function(double s)                               //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet eta star function defined as     //
//                           eta*(s) = eta(s) - 1,                            //
//     where eta() is the Dirichlet eta function.                             //
//                                                                            //
//  Arguments:                                                                //
//     double s                                                               //
//        The argument of the eta* function defined above.                    //
//                                                                            //
//  Return Value:                                                             //
//     eta*(s), if eta*(s) > DBL_MAX, then DBL_MAX is returned,               //
//              if eta*(s) < -DBL_MAX, then -DBL_MAX is returned.             //
//                                                                            //
//  Example:                                                                  //
//     double s;                                                              //
//     double eta_star;                                                       //
//                                                                            //
//     ( User code to set s, the argument to the eta* function. )             //
//                                                                            //
//     eta_star = Dirichlet_Eta_Star_Function(s);                             //
////////////////////////////////////////////////////////////////////////////////

double Dirichlet_Eta_Star_Function(double s)
{
   long double x;

   x = xDirichlet_Eta_Star_Function((long double)s);

   if (fabsl(x) < DBL_MAX) return (double) x;
   else return ( x > 0.0L ) ? DBL_MAX : - DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xDirichlet_Eta_Star_Function(long double s)                    //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Dirichlet eta star function defined as     //
//                           eta*(s) = eta(s) - 1,                            //
//     where eta() is the Dirichlet eta function.                             //
//                                                                            //
//  Arguments:                                                                //
//     long double s;                                                         //
//        The argument of the eta* function defined above.                    //
//                                                                            //
//  Return Value:                                                             //
//     eta*(s)                                                                //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double eta_star;                                                  //
//                                                                            //
//     ( User code to set s, the argument to the eta* function. )             //
//                                                                            //
//     eta_star = xDirichlet_Eta_Star_Function(s);                            //
////////////////////////////////////////////////////////////////////////////////

long double xDirichlet_Eta_Star_Function(long double s)
{
   long double reflection_coefficient;

   if (s >= 0.0L)
      if ( s < 18.0L ) 
         return Alternating_Series_Convergence_Acceleration(s);
      else
         return Sum_Reverse_Order(s);

                       // Use the reflection formula //

   reflection_coefficient = Reflection_Coefficient(s);    
   return reflection_coefficient * xDirichlet_Eta_Star_Function(1.0L - s)
          + (reflection_coefficient - 1.0L);
}


////////////////////////////////////////////////////////////////////////////////
// static long double Reflection_Coefficient(long double s)                   //
//                                                                            //
//  Description:                                                              //
//     The reflection formula for the Dirichlet's eta function is:            //
//     eta(s) = { [2*(1-2^(1-s))/(1-2^s)] * Gamma(1-s)                        //
//                     / (2*pi)^(1-s)*cos((1-s)*pi/2) } *  eta(1-s)           //
//     where s < 0 is a negative real number.                                 //
//                                                                            //
//     This routine returns the coefficient:                                  //
//       [2*(1-2^(1-s))/(1-2^s)] Gamma(1-s) / (2*pi)^(1-s)*cos((1-s)*pi/2)    //
//     where s < 0 is a negative real number.                                 //
//                                                                            //
//  Arguments:                                                                //
//     long double s                                                          //
//        The argument of the eta function, where s < 0.                      //
//                                                                            //
//  Return Value:                                                             //
//     The coefficient c so that eta(s) = c * eta(1-s)                        //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double coef;                                                      //
//                                                                            //
//     ( User code to set s < 0, the argument to the eta function. )          //
//                                                                            //
//     coef = Reflection_Coefficient(s);                                      //
////////////////////////////////////////////////////////////////////////////////

static long double Reflection_Coefficient(long double s)
{
   long double one_s = 1.0L - s;
   int k = (int) (one_s / 4.0L);
   long double v = one_s - 4.0L * (long double) k;
   long double c = cosl(v * pi_2);
   long double x, temp;

   if (fabsl(c) < 1.8L * LDBL_EPSILON) return 0.0L;
   temp = powl(2.0L,one_s);
   x = (1.0L - temp) / (temp - 2.0L);
   x += x;
   x *= c;
   x *= xGamma_Function(one_s);
   x /= powl(pi, one_s);
   return x;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Alternating_Series_Convergence_Acceleration             //
//                                                            (long double s) //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Eta function for s <= 18 where the series  //
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
//        The argument of the eta function.                                   //
//                                                                            //
//  Return Value:                                                             //
//     eta(s)                                                                 //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double eta;                                                       //
//                                                                            //
//     ( User code to set s >= 18, the argument to the eta function.          //
//                                                                            //
//     eta = Alternating_Series_Convergence_Acceleration(s);                  //
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
      term[k-1] = d[k] * powl((long double)(k+1),-s);
      k++;
      term[k-1] = -d[k] * powl((long double)(k+1),-s);
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
//     This routine calculates the Eta function for s >= 18 (this requires    //
//     fewer than 28 calls to powl())                                         //
//                      eta(s) = Sum (-1)^(k-1) (1/k^s),                      //
//     where the sum is summed over k = 1,... .                               //
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
//        The argument of the eta function.                                   //
//                                                                            //
//  Return Value:                                                             //
//     eta(s)                                                                 //
//                                                                            //
//  Example:                                                                  //
//     long double s;                                                         //
//     long double eta;                                                       //
//                                                                            //
//     ( User code to set s >= 18, the argument to the eta function. )        //
//                                                                            //
//     eta = Sum_Reverse_Order(s);                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
static long double Sum_Reverse_Order(long double s)
{
   long double term[30];
   long double lower_bound;
   long double upper_bound;
   long double sum;
   int k;

   lower_bound = -powl(2.0L,-s);
   term[0] = lower_bound;
   term[1] = powl(3.0L,-s);
   upper_bound = term[1] + lower_bound;
   if ( lower_bound == upper_bound) return upper_bound;
   for (k = 4; k < 32; k++) {
      term[k-2] = -powl( (long double) k, -s);
      lower_bound = upper_bound + term[k-2];
      if (lower_bound == upper_bound) break;
      k++;
      term[k-2] = powl( (long double) k, -s);
      upper_bound = lower_bound + term[k-2];
      if (lower_bound == upper_bound) break;
   }
   k -= 2;
   sum = term[k];
   for (k--; k >= 0; k--) sum += term[k];
   return sum;
}
