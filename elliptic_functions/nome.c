////////////////////////////////////////////////////////////////////////////////
// File: nome.c                                                               //
// Routine(s):                                                                //
//    Nome                                                                    //
////////////////////////////////////////////////////////////////////////////////

//                         Externally Defined Routines                        //

extern double Complete_Elliptic_Integral_First_Kind(char arg, double x);

////////////////////////////////////////////////////////////////////////////////
// void Nome(double k)                                                        //
//                                                                            //
//  Description:                                                              //
//     This function calculates the nome q = exp(- pi * K' / K) where K' is   //
//     the complete elliptic integral of the first kind with modulus          //
//     k' = sqrt(1 - k) and K is the complete elliptic integral of the first  //
//     kind with modulus k.                                                   //
//     The modulus, k, must satisfy |k| <= 1.  If k = 0 then the integral K'  //
//     is infinite and q = 0.  If |k| = 1, then the integral K is infinite    //
//     and q = 1.                                                             //
//                                                                            //
//  Arguments:                                                                //
//     double  k                                                              //
//                The the modulus of the complete elliptic integral of the    //
//                first kind, K.                                              //
//                                                                            //
//  Return Value:                                                             //
//     The nome q.  If |k| > 1, then q = 1 is returned.  The user should      //
//     verify that |k| <= 1.                                                  //
//                                                                            //
//  Example:                                                                  //
//     double k;                                                              //
//     double nome;                                                           //
//                                                                            //
//     ( code to initialize k )                                               //
//                                                                            //
//     nome = Nome(k);                                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>        // required for fabs(), exp() and M_PI
#include <float.h>

double Nome(double k) 
{
   double K, cK;

   if ( k == 0.0 ) return 0.0;
   if ( fabs(k) >= 1.0 ) return 1.0;
   K = Complete_Elliptic_Integral_First_Kind('k',k);
   cK = Complete_Elliptic_Integral_First_Kind('m', 1.0 - k * k);
   return exp(-M_PI * cK / K); 
}
