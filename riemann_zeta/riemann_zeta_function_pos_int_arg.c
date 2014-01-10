////////////////////////////////////////////////////////////////////////////////
// File: riemann_zeta_function_pos_int_arg.c                                  //
// Routine(s):                                                                //
//    Riemann_Zeta_Function_pos_int_arg                                       //
//    xRiemann_Zeta_Function_pos_int_arg                                      //
//    Riemann_Zeta_Star_Function_int_arg                                      //
//    xRiemann_Zeta_Star_Function_int_arg                                     //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for powl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX
#define CUTOFF 66          // first argument for which zeta*() < LDBL_EPSILON

//                         Internally Defined Routines                        //

double      Riemann_Zeta_Function_pos_int_arg(int s);
long double xRiemann_Zeta_Function_pos_int_arg(int s);
double      Riemann_Zeta_Star_Function_pos_int_arg(int s);
long double xRiemann_Zeta_Star_Function_pos_int_arg(int s);

                 // zeta_star table: zeta_star[0] = zeta*(2) //

static long double zeta_star[] = {
    6.449340668482264364724e-1L,  2.020569031595942853997e-1L,
    8.232323371113819151600e-2L,  3.692775514336992633137e-2L,
    1.734306198444913971452e-2L,  8.349277381922826839798e-3L,
    4.077356197944339378685e-3L,  2.008392826082214417853e-3L,
    9.945751278180853371461e-4L,  4.941886041194645587042e-4L,
    2.460865533080482986443e-4L,  1.227133475784891467815e-4L,
    6.124813505870482926574e-5L,  3.058823630702049355340e-5L,
    1.528225940865187173295e-5L,  7.637197637899762273681e-6L,
    3.817293264999839856479e-6L,  1.908212716553938925660e-6L,
    9.539620338727961131527e-7L,  4.769329867878064631168e-7L,
    2.384505027277329900037e-7L,  1.192199259653110730678e-7L,
    5.960818905125947961244e-8L,  2.980350351465228018606e-8L,
    1.490155482836504123466e-8L,  7.450711789835429491981e-9L,
    3.725334024788457054819e-9L,  1.862659723513049006404e-9L,
    9.313274324196681828718e-10L, 4.656629065033784072989e-10L,
    2.328311833676505492001e-10L, 1.164155017270051977593e-10L,
    5.820772087902700889244e-11L, 2.910385044497099686929e-11L,
    1.455192189104198423593e-11L, 7.275959835057481014521e-12L,
    3.637979547378651190237e-12L, 1.818989650307065947585e-12L,
    9.094947840263889282533e-13L, 4.547473783042154026799e-13L,
    2.273736845824652515227e-13L, 1.136868407680227849349e-13L,
    5.684341987627585609277e-14L, 2.842170976889301855455e-14L,
    1.421085482803160676983e-14L, 7.105427395210852712877e-15L,
    3.552713691337113673298e-15L, 1.776356843579120327473e-15L,
    8.881784210930815903096e-16L, 4.440892103143813364198e-16L,
    2.220446050798041983999e-16L, 1.110223025141066133721e-16L,
    5.551115124845481243724e-17L, 2.775557562136124172582e-17L,
    1.387778780972523276284e-17L, 6.938893904544153697446e-18L,
    3.469446952165922624744e-18L, 1.734723476047576572049e-18L,
    8.673617380119933728342e-19L, 4.336808690020650487497e-19L,
    2.168404344997219785014e-19L, 1.084202172494241406301e-19L,
    5.421010862456645410919e-20L, 2.710505431223468831955e-20L,
    1.355252715610116458149e-20L, 6.776263578045189097995e-21L,
    3.388131789020796818086e-21L, 1.694065894509799165406e-21L,
    8.470329472546998348247e-22L, 4.235164736272833347862e-22L,
    2.117582368136194731844e-22L, 1.058791184068023385227e-22L,
    5.293955920339870323814e-23L, 2.646977960169852961134e-23L,
    1.323488980084899080309e-23L, 6.617444900424404067355e-24L,
    3.308722450212171588947e-24L, 1.654361225106075646230e-24L,
    8.271806125530344403671e-25L, 4.135903062765160926009e-25L,
    2.067951531382576704396e-25L, 1.033975765691287099328e-25L,
    5.169878828456431320410e-26L, 2.584939414228214268128e-26L,
    1.292469707114106670038e-26L, 6.462348535570531803438e-27L,
    3.231174267785265386135e-27L, 1.615587133892632521206e-27L,
    8.077935669463162033159e-28L, 4.038967834731580825622e-28L,
    2.019483917365790349159e-28L, 1.009741958682895153362e-28L,
    5.048709793414475696085e-29L, 2.524354896707237824467e-29L,
    1.262177448353618904375e-29L, 6.310887241768094495683e-30L,
    3.155443620884047239110e-30L, 1.577721810442023616644e-30L,
    7.888609052210118073521e-31L, 3.944304526105059033526e-31L,
    1.972152263052529515685e-31L, 9.860761315262647574833e-32L,
    4.930380657631323786219e-32L, 2.465190328815661892710e-32L,
    1.232595164407830946222e-32L, 6.162975822039154730666e-33L,
    3.081487911019577365185e-33L, 1.540743955509788682543e-33L,
    7.703719777548943412553e-34L, 3.851859888774471706221e-34L,
    1.925929944387235853092e-34L, 9.629649721936179265402e-35L,
    4.814824860968089632681e-35L, 2.407412430484044816333e-35L,
    1.203706215242022408164e-35L, 6.018531076210112040815e-36L,
    3.009265538105056020405e-36L, 1.504632769052528010202e-36L,
    7.523163845262640051005e-37L, 3.761581922631320025502e-37L,
    1.880790961315660012751e-37L, 9.403954806578300063752e-38L,
    4.701977403289150031876e-38L, 2.350988701644575015938e-38L,
    1.175494350822287507969e-38L, 5.877471754111437539844e-39L,
    2.938735877055718769922e-39L, 1.469367938527859384961e-39L,
    7.346839692639296924805e-40L, 3.673419846319648462402e-40L,
    1.836709923159824231201e-40L, 9.183549615799121156006e-41L,
    4.591774807899560578003e-41L, 2.295887403949780289001e-41L,
    1.147943701974890144501e-41L, 5.739718509874450722504e-42L,
    2.869859254937225361252e-42L, 1.434929627468612680626e-42L,
    7.174648137343063403129e-43L, 3.587324068671531701565e-43L,
    1.793662034335765850782e-43L, 8.968310171678829253912e-44L,
    4.484155085839414626956e-44L, 2.242077542919707313478e-44L,
    1.121038771459853656739e-44L, 5.605193857299268283695e-45L,
    2.802596928649634141847e-45L, 1.401298464324817070924e-45L,
    7.006492321624085354619e-46L, 3.503246160812042677309e-46L,
    1.751623080406021338655e-46L, 8.758115402030106693273e-47L,
    4.379057701015053346637e-47L, 2.189528850507526673318e-47L,
    1.094764425253763336659e-47L, 5.473822126268816683296e-48L,
    2.736911063134408341648e-48L, 1.368455531567204170824e-48L,
    6.842277657836020854120e-49L, 3.421138828918010427060e-49L,
    1.710569414459005213530e-49L, 8.552847072295026067650e-50L,
    4.276423536147513033825e-50L, 2.138211768073756516912e-50L,
    1.069105884036878258456e-50L, 5.345529420184391292281e-51L,
    2.672764710092195646141e-51L, 1.336382355046097823070e-51L,
    6.681911775230489115351e-52L, 3.340955887615244557676e-52L,
    1.670477943807622278838e-52L, 8.352389719038111394189e-53L,
    4.176194859519055697095e-53L, 2.088097429759527848547e-53L,
    1.044048714879763924274e-53L, 5.220243574398819621368e-54L,
    2.610121787199409810684e-54L, 1.305060893599704905342e-54L,
    6.525304467998524526710e-55L, 3.262652233999262263355e-55L,
    1.631326116999631131678e-55L, 8.156630584998155658388e-56L,
    4.078315292499077829194e-56L, 2.039157646249538914597e-56L,
    1.019578823124769457298e-56L, 5.097894115623847286492e-57L,
    2.548947057811923643246e-57L, 1.274473528905961821623e-57L,
    6.372367644529809108116e-58L, 3.186183822264904554058e-58L,
    1.593091911132452277029e-58L, 7.965459555662261385144e-59L,
    3.982729777831130692572e-59L, 1.991364888915565346286e-59L,
    9.956824444577826731431e-60L, 4.978412222288913365715e-60L,
    2.489206111144456682858e-60L, 1.244603055572228341429e-60L,
    6.223015277861141707144e-61L };


////////////////////////////////////////////////////////////////////////////////
// double Riemann_Zeta_Function_pos_int_arg(int s)                            //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemann zeta function                      //
//        zeta(s) = Sum (1/k^s), summed over k = 1,... for s > 1 an integer.  //
//        zeta(s) = -(1/2) for s = 0.                                         //
//        zeta(s) = DBL_MAX for s = 1, note that s = 1 is a simple pole of the//
//     Riemann zeta function and lim zeta(s) = -inf where s->1 from the left  //
//     and lim zeta(s) = inf where s->1 from the right.                       //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Riemann zeta function restricted to the         //
//        nonnegative integers.  If s = 0, then zeta(0) = -1/2 is returned.   //
//        If s = 1, then DBL_MAX is returned. If s > 1, then Sum(1/k^s)       //
//        summed for k = 1,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if s = 1, then DBL_MAX is returned but note that s = 1 is     //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double zeta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta function. )     //
//                                                                            //
//     if ( s >= 0 ) zeta = Riemann_Zeta_Function_pos_int_arg(s);             //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

double Riemann_Zeta_Function_pos_int_arg(int s)
{
   if (s < 0) return 0.0;
   if (s > CUTOFF) return 1.0;
   if (s > 1) return (double) xRiemann_Zeta_Function_pos_int_arg(s);
   if (s == 0) return -0.5;
   return DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// long double xRiemann_Zeta_Function_pos_int_arg(int s)                      //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemann zeta function                      //
//        zeta(s) = Sum (1/k^s), summed over k = 1,... for s > 1 an integer.  //
//        zeta(s) = -(1/2) for s = 0.                                         //
//        zeta(s) = LDBL_MAX for s = 1, note that s = 1 is a simple pole of   //
//     the Riemann zeta function and lim zeta(s) = -inf where s->1 from the   //
//     left and lim zeta(s) = inf where s->1 from the right.                  //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the Riemann zeta function restricted to the         //
//        nonnegative integers.  If s = 0, then zeta(0) = -1/2 is returned.   //
//        If s = 1, then LDBL_MAX is returned. If s > 1, then Sum(1/k^s)      //
//        summed for k = 1,... is returned.                                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta(s), if s = 1, then LDBL_MAX is returned but note that s = 1 is    //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double zeta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta function. )     //
//                                                                            //
//     if ( s >= 0 ) zeta = Riemann_Zeta_Function_pos_int_arg(s);             //
//     else printf(" Error argument, %d, is negative\n",s);                   //
////////////////////////////////////////////////////////////////////////////////

long double xRiemann_Zeta_Function_pos_int_arg(int s)
{
   if (s < 0) return 0.0L;
   if (s > CUTOFF) return 1.0L;
   if (s > 1) return 1.0L + xRiemann_Zeta_Star_Function_pos_int_arg(s);
   if (s == 0) return -0.5L;
   return LDBL_MAX;
}

////////////////////////////////////////////////////////////////////////////////
// double Riemann_Zeta_Star_Function_pos_int_arg(int s)                       //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemann zeta star function defined for     //
//     for integer arguments s > 1 as zeta*(s) = zeta(s) - 1, where zeta() is //
//     the Riemann zeta function, i.e. zeta*(s) = Sum (1/k^s), summed over    //
//     k = 2,... for s > 1 an integer.  zeta*(s) has a simple pole at s = 1,  //
//     in which case DBL_MAX is returned.  zeta*(0) = -3/2.                   //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the zeta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta*(s), if s = 1, then DBL_MAX is returned but note that s = 1 is    //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     double zeta;                                                           //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta* function. )    //
//                                                                            //
//     if ( s > 1 ) zeta = Riemann_Zeta_Star_Function_pos_int_arg(s);         //
//     else printf(" Error argument, %d, is less than or equal to 1.\n",s);   //
////////////////////////////////////////////////////////////////////////////////

double Riemann_Zeta_Star_Function_pos_int_arg(int s)
{
   if (s == 0) return -1.5;
   if (s == 1) return DBL_MAX;
   return (double) xRiemann_Zeta_Star_Function_pos_int_arg(s);
}

////////////////////////////////////////////////////////////////////////////////
// long double xRiemann_Zeta_Star_Function_pos_int_arg(int s)                 //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the Riemann zeta star function defined for     //
//     for integer arguments s > 1 as zeta*(s) = zeta(s) - 1, where zeta() is //
//     the Riemann zeta function, i.e. zeta*(s) = Sum (1/k^s), summed over    //
//     k = 2,... for s > 1 an integer.  zeta*(s) has a simple pole at s = 1,  //
//     in which case LDBL_MAX is returned.  zeta*(0) = -3/2.                  //
//                                                                            //
//  Arguments:                                                                //
//     int s                                                                  //
//        The argument of the zeta* function defined above.                   //
//                                                                            //
//  Return Value:                                                             //
//     zeta*(s), if s = 1, then LDBL_MAX is returned but note that s = 1 is   //
//              a simple pole of the zeta() function.                         //
//                                                                            //
//  Example:                                                                  //
//     int    s;                                                              //
//     long double zeta;                                                      //
//                                                                            //
//     ( User code to set s, the argument to the Riemann zeta* function. )    //
//                                                                            //
//     if ( s > 1 ) zeta = xRiemann_Zeta_Star_Function_pos_int_arg(s);        //
//     else printf(" Error argument, %d, is less than or equal to 1.\n",s);   //
////////////////////////////////////////////////////////////////////////////////

long double xRiemann_Zeta_Star_Function_pos_int_arg(int s)
{
   static int max_zeta_star_index = sizeof(zeta_star) / sizeof(long double) - 1;
   int k = s - 2;
   int diff;

   if (s == 0) return -1.5L;
   if (s == 1) return LDBL_MAX;
   if ( k <= max_zeta_star_index ) return zeta_star[k];

           // For large positive integer s, use the approximation //
           //               eta*(s+1) = eta*(s) / 2.              //

   diff = s - max_zeta_star_index - 2;
   return zeta_star[max_zeta_star_index] / powl(2.0L, (long double) diff);
}
