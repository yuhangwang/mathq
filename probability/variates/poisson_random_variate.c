////////////////////////////////////////////////////////////////////////////////
// File: poisson_random_variate.c                                             //
// Routine(s):                                                                //
//    Poisson_Random_Variate                                                  //
////////////////////////////////////////////////////////////////////////////////

//                    Required Externally Defined Routines                    //

double Exponential_Random_Variate( void );
double Gamma_Random_Variate( double a );
int    Binomial_Random_Variate( int n, double p );


//                    Required Internally Defined Routines                    //

static int Inter_Arrival_Time( double mu );

#define THRESHOLD 6

////////////////////////////////////////////////////////////////////////////////
// int Poisson_Random_Variate( double mu )                                    //
//                                                                            //
//  Description:                                                              //
//     This function returns a Poisson distributed random variate using       //
//     the recursive method:  For mu > THRESHOLD, choose an integer n < mu,   //
//     let g be a Gamma_Random_Variate(n) random variate, then if g <= mu,    //
//     then return n + Poisson_Random_Variate(mu - g) and if g > mu, then     //
//     return Binomial_Random_Variate(n-1,mu/g).   For mu <= THRESHOLD, use   //
//     the Inter_Arrival_Time method.                                         //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Values:                                                            //
//     A Poisson distributed random quantity.                                 //
//                                                                            //
//  Example:                                                                  //
//     int x;                                                                 //
//     double mu;                                                             //
//                                                                            //
//     x = Poisson_Random_Variate(mu);                                        //
////////////////////////////////////////////////////////////////////////////////

int Poisson_Random_Variate( double mu )
{
   double g;
   int x;
   int n;
 
   if (mu <= THRESHOLD) return Inter_Arrival_Time(mu);
   n = (int) (0.5 * mu);
   g = Gamma_Random_Variate( (double) n );
   if ( g <= mu ) return n + Poisson_Random_Variate(mu - g);
   return Binomial_Random_Variate(n-1, mu / g);
}


////////////////////////////////////////////////////////////////////////////////
// static int Inter_Arrival_Time( double mu )                                 //
//                                                                            //
//  Description:                                                              //
//     This function returns a Poisson distributed random variate with mean   //
//     mu using the sum of exponential random variates.                       //
//                                                                            //
//  Arguments:                                                                //
//     double mu                                                              //
//            The mean of the Poisson distibution: mu^x exp(-mu) / x!.        //
//                                                                            //
//  Return Values:                                                            //
//     A Poisson distributed random quantity.                                 //
//                                                                            //
//  Example:                                                                  //
//     double mu;                                                             //
//                                                                            //
//     Inter_Arrival_Time( mu );                                              //
////////////////////////////////////////////////////////////////////////////////

static int Inter_Arrival_Time( double mu )
{
   int x = 0;
   double sum = 0.0;
  
   while (sum <= mu) {
      x++;
      sum += Exponential_Random_Variate();
   }
   return x - 1;
}
