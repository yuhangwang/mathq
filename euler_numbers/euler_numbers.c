////////////////////////////////////////////////////////////////////////////////
// File: euler_numbers.c                                                      //
// Routine(s):                                                                //
//    Euler_Number                                                            //
//    Euler_Number_Sequence                                                   //
//    Euler_Even_Number_Sequence                                              //
//    Max_Euler_Even_Number_Index                                             //
//    xEuler_Number                                                           //
//    xEuler_Number_Sequence                                                  //
//    xEuler_Even_Number_Sequence                                             //
//    xMax_Euler_Even_Number_Index                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     There are several definitions of the Euler numbers, E[n].  The         //
//     definition programed here yields E[0] = 1, E[1] = 0, E[2] = -1,        //
//     E[3] = 0, E[4] = 5, ... .  This is equivalent to using the generating  //
//     function 2 / (exp(x) - exp(-x)) = Sum E[n] x^n / n! where the          //
//     sum extends over n = 0,1,... .                                         //
//     Multiply both sides of generating function by (exp(x) - exp(-x)) / 2   //
//     and express (exp(x) - exp(-2)) / 2 in terms of its power series,       //
//                      1 = Sum (1+(-1)^k)/2 x^k / k! Sum E[j] x^j / j!       //
//                      1 = Sum Sum d[k] E[j] x^(j+k) / (j! k!)               //
//     where d[k] = 0 if k is odd and d[k] = 1 if k is even.                  //
//     Let n = j + 2k, then                                                   //
//                      1 = Sum Sum d[n-j] E[j] x^n / (j! (n - j)!)           //
//     where the first sum extends from n = 0, 1, ..., and the second sum     //
//     extends from j = 0, .., n.  Equating like powers of x on both sides    //
//                      1 = E[0]                                              //
//                      0 = E[1]                                              //
//                      0 = Sum E[2j+1] / [(2j+1)! (n-2j-1)!], for n >= 1,    //
//                        n odd, and the sum extends over j = 0,...,(n-1)/2,  //
//                      0 = Sum E[2j] / [(2j)! (n-2j)!], for n >= 2,          //
//                        n even, and the sum extends over j = 0,...,n/2.     //
//     Therefore                                                              //
//         (1)          E[0] = 1, E[1] = 0                                    //
//                      E[2n+1] = 0, n = 1,2,...                              //
//                      E[2n] = - (2n)! Sum E[2j] / [(2j)! (2n-2j)!] where    //
//                        the sum extends over j = 0, ... n-1, n = 1,2,... .  //
//                                                                            //
//     The even-indexed Euler numbers also satisy the expression              //
//     involving the beta reciprocal sum function:                            //
//                      E[2n] = (-1)^n [ 2 (2n)! (2/pi)^(2n+1) ] beta(2n+1)   //
//     where beta(2n+1) = 1 - (1/3)^(2n+1) + (1/5)^(2n+1) - ... .             //
//     Given that LDBL_EPSILON ~ 10^(-20), the beta function is               //
//     machine computationally equal to 1.0 for 2n >= 40.                     //
//     Therefore for 2n >= 40,                                                //
//         (2)          E[2n] = (-1)^n [ 2 (2n)! (2/pi)^(2n+1) ]              //
//     or                                                                     //
//         (3)          E[2n] = -(2n)(2n-1)(2/pi)^2 E[2n-2].                  //
//                                                                            //
//     Given that LDBL_MAX ~ 10^(4932), using equation (2), the maximum       //
//     even-indexed argument is 1866.  Note that all odd-indexed Euler        //
//     numbers are 0.                                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for fabsl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double      Euler_Number(int n);
void        Euler_Number_Sequence(double E[], int start, int length);
void        Euler_Even_Number_Sequence(double E[], int start, int length);
int         Max_Euler_Even_Number_Index(void);
long double xEuler_Number(int n);
void        xEuler_Number_Sequence(long double E[], int start, int length);
void        xEuler_Even_Number_Sequence(long double E[], int start,
                                                                   int length);
int         xMax_Euler_Even_Number_Index(void);

static long double Iterate_Up(long double en, int n, int m);
static long double Iterate_Down(long double en, int n, int m);

//                         Internally Defined Constants                       //

 // even indexed euler numbers E[2n] = euler_number[n], n = 0,... //

static long double euler_numbers[] = {
         +1.00000000000000000000000000000000000L,
         -1.00000000000000000000000000000000000L,
         +5.00000000000000000000000000000000000L,
         -6.10000000000000000000000000000000000e+1L,
         +1.38500000000000000000000000000000000e+3L,
         -5.05210000000000000000000000000000000e+4L,
         +2.70276500000000000000000000000000000e+6L,
         -1.99360981000000000000000000000000000e+8L,
         +1.93915121450000000000000000000000000e+10L,
         -2.40487967544100000000000000000000000e+12L,
         +3.70371188237525000000000000000000000e+14L,
         -6.93488743931379010000000000000000000e+16L,
         +1.55145341635570869050000000000000000e+19L,
         -4.08707250929312389236100000000000000e+21L,
         +1.25225964140362986546828500000000000e+24L,
         -4.41543893249023104553682821000000000e+26L,
         +1.77519391579539289436664789665000000e+29L,
         -8.07232992358878980621682474532810000e+31L,
         +4.12220603395177021223470796712590450e+34L,
         -2.34895805270431082520178285761989477e+37L,
         +1.48511507181149800178771567814058267e+40L,
         -1.03646227335196121193979573047451860e+43L,
         +7.94757942259759270360804051008807062e+45L,
         -6.66753751668554497743502847477374820e+48L,
         +6.09627864556854215869168574287684315e+51L,
         -6.05328524818862189631438378511164909e+54L,
         +6.50616248668460884771587063408082298e+57L,
         -7.54665993900873909806143256588973674e+60L,
         +9.42032189642024120420228623769058323e+63L,
         -1.26220192518062187199034092372874893e+67L,
         +1.81089114965792304965458077416521587e+70L,
         -2.77571017020715805973669809083715274e+73L,
         +4.53581033300178891747468878715677624e+76L,
         -7.88628420666178941810072074223999042e+79L,
         +1.45618443801396315007150470094942327e+83L,
         -2.85051783223697718732198729556739340e+86L,
         +5.90574720777544365455135032296439571e+89L,
         -1.29297366418786417049760323593869875e+93L,
         +2.98692818328457695093074365221714061e+96L,
         -7.27060171401686414380328065169928185e+99L,
         +1.86229157584126970444824923030431260e+103L,
         -5.01310494081097966129086936788810094e+106L,
         +1.41652557597856259916722069410021670e+110L,
         -4.19664316404024471322573414069418892e+113L,
         +1.30215959052404639812585869133081868e+117L,
         -4.22724068613990906470558992921459310e+120L,
         +1.43432127919765834061336826405785659e+124L,
         -5.08179907245804251645597576430907360e+127L,
         +1.87833293645293026402007579184179893e+131L,
         -7.23653438103385777657187661736782293e+134L,
         +2.90352834666109749705460383476443588e+138L,
         -1.21229373789292182105392954978560988e+142L,
         +5.26306424961699070600224073584236661e+145L,
         -2.37407307193676634703461698760652652e+149L,
         +1.11189009424828230249702335881757893e+153L,
         -5.40307865979529320561911549426347699e+156L,
         +2.72234108557222702137153414458909549e+160L,
         -1.42130105480096698118085204572231882e+164L,
         +7.68426182064690265317095628366647794e+167L,
         -4.29962192543974964281889033648632755e+171L,
         +2.48839157478298716316902455408489408e+175L,
         -1.48875820890620408401048810913362396e+179L,
         +9.20261411885209418840864126560312709e+182L,
         -5.87424445729243560747806550051798443e+186L,
         +3.87013355417592724899726125339465800e+190L,
         -2.63038464627282201918918005755736145e+194L,
         +1.84342186190681643216739318103276967e+198L,
         -1.33150076083199759777989619061195919e+202L,
         +9.90773407946409970275719941594148144e+205L,
         -7.59161615376086554230567716763177264e+209L,
         +5.98738690421595478060934030092899051e+213L,
         -4.85853153680527007166022567445774339e+217L,
         +4.05474737750791455464680535308584710e+221L,
         -3.47892371339090601415585327133292340e+225L,
         +3.06749738825108489449144357479461161e+229L,
         -2.77857404780457414987248665136951661e+233L,
         +2.58465603902711815098815082730837912e+237L,
         -2.46817048046364050455631133967404223e+241L,
         +2.41875397603671333264713788326666700e+245L,
         -2.43169264709107277171036789982532904e+249L,
         +2.50718300057371449601915222347628344e+253L,
         -2.65025200052581375350895159803901660e+257L,
         +2.87130197316667968492991621100369935e+261L,
         -3.18736021623541104699251674698644208e+265L,
         +3.62424164505845624987618515668413679e+269L,
         -4.22000551313026080825687414912160887e+273L,
         +5.03034557853150041609481420707106604e+277L,
         -6.13696178494213385049453688204944205e+281L,
         +7.66062813846337323811799348691311731e+285L,
         -9.78178011283967454892036825005468034e+289L,
         +1.27733166367198064207287773215186928e+294L,
         -1.70535141854472052178024263787253627e+298L,
         +2.32725003482003005917234767874590751e+302L,
         -3.24554745838924695277710327883293385e+306L};

 // Euler numbers even arguments E[20n+198] = xeuler_number[n],//
 // n = 0,...,83. //

static long double xeuler_numbers[] = {
         -3.71689279117523442595544500254463863e+331L,
         -1.06311359653345223864436771012921321e+374L,
         -1.90168660939480511426680883012390551e+417L,
         -1.82360055378359429711155051316921460e+461L,
         -8.23027064165985474131958899398954349e+505L,
         -1.56385817382680568575406775082477747e+551L,
         -1.13601858351201446481652700044863755e+597L,
         -2.89972526042905613264942493451548329e+643L,
         -2.41449160212277660702694362235429277e+690L,
         -6.13923642063764529502742716554716146e+737L,
         -4.49345452991421399120978866459294351e+785L,
         -8.97756427858091603493628046427874307e+833L,
         -4.66649628412561566112495499374042969e+882L,
         -6.04121575288815948930398991671029254e+931L,
         -1.87180552718945404097300923918204788e+981L,
         -1.33828725733423286269574102549431884e+1031L,
         -2.13507148521639810691301131823700866e+1081L,
         -7.36881516180485529387395142592851172e+1131L,
         -5.34629271818275753228776887082132305e+1182L,
         -7.93994931537857379076878224802396551e+1233L,
         -2.35468354083176548121879172953683941e+1285L,
         -1.36256313991614050926850541390176391e+1337L,
         -1.50556372657638989747009698111842025e+1389L,
         -3.11284733054824529224374576953616573e+1441L,
         -1.18158607487418410261222802947146061e+1494L,
         -8.08801726111725843172147531651065154e+1546L,
         -9.81643809823724419884007786352216236e+1599L,
         -2.07910013133804812354853355753212500e+1653L,
         -7.56930837541337892139699683121856588e+1706L,
         -4.66972792066728622296744289609326778e+1760L,
         -4.81610248406603883098616688457328879e+1814L,
         -8.19743050944023143124899502165823523e+1868L,
         -2.27467877743666711876140790482727335e+1923L,
         -1.01708883510992495528701764597271305e+1978L,
         -7.24713502324577739139163393961357841e+2032L,
         -8.14206430760392132452579711261951630e+2087L,
         -1.42777791267301339879573300294254805e+2143L,
         -3.87020326879714541757030010244668726e+2198L,
         -1.60665423922633014521515539843501180e+2254L,
         -1.01242632924908849366960604899985285e+2310L,
         -9.60177614062521427824715338407438326e+2365L,
         -1.35934628187734525852976141636933705e+2422L,
         -2.85024353254759106702939719903467501e+2478L,
         -8.78459548954654778895625926626374915e+2534L,
         -3.95082936134681306929563741719263829e+2591L,
         -2.57476454552965258406462266542008726e+2648L,
         -2.41510717284989111843875921913046006e+2705L,
         -3.23934156259324876775018008688539071e+2762L,
         -6.17404394147989101444878513450330120e+2819L,
         -1.66203723229668229358143280398234422e+2877L,
         -6.28238716467387724217301253134289069e+2934L,
         -3.31559423076960589817743629291310796e+2992L,
         -2.42980470187528459765452596366063711e+3050L,
         -2.45952997288228016987778866037825470e+3108L,
         -3.42116693869848570260313245422990149e+3166L,
         -6.50696101102246931508146345414188986e+3224L,
         -1.68411158327294985991824801902019600e+3283L,
         -5.90367820482716088594963172965098921e+3341L,
         -2.79038893561987423773706900511763896e+3400L,
         -1.77046003842729191681938978065902674e+3459L,
         -1.50152085350928546047713454578423983e+3518L,
         -1.69511069697400132340540764825576045e+3577L,
         -2.53708286823564266907689208268515865e+3636L,
         -5.01462288249826981432317099431657482e+3695L,
         -1.30392553956406453204046059515359276e+3755L,
         -4.44392637758750127722200395339888769e+3814L,
         -1.97793559595110527837198995052984148e+3874L,
         -1.14567558496959386333167877391244865e+3934L,
         -8.60653003813875214291824475170620541e+3993L,
         -8.35722918866828791359937607765858591e+4053L,
         -1.04556714366068261009547015886547413e+4114L,
         -1.68004099070122568539355287290779560e+4174L,
         -3.45639321975421581752867245423113382e+4234L,
         -9.07719145861269666305648782790544439e+4294L,
         -3.03407494116294501445636329338932582e+4355L,
         -1.28706183526776733830096061273568150e+4416L,
         -6.90958658904844127325973020523167394e+4476L,
         -4.68159529384680939138368853648510930e+4537L,
         -3.99265538426709382066133671034331935e+4598L,
         -4.27483540993886283363581320122434618e+4659L,
         -5.73133248005026930526104603547695851e+4720L,
         -9.59807574963659813617059370196913709e+4781L,
         -2.00282376205434874170732124553349950e+4843L,
         -5.19507751303279845445675393607490318e+4904L
};

static long double pi_over_2_squared = 2.46740110027233965470862274996903778L;
static int max_double_index = 2*sizeof(euler_numbers) / sizeof(long double) - 2;
static int max_long_double_index = 1866;

////////////////////////////////////////////////////////////////////////////////
// double Euler_Number(int n)                                                 //
//                                                                            //
//  Description:                                                              //
//     This routine returns the Euler number E[n], where n >= 0.              //
//     If E[n] > DBL_MAX, then DBL_MAX is returned, if E[n] < -DBL_MAX, then  //
//     -DBL_MAX is returned.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     int n                                                                  //
//        The index of the sequence of Euler numbers, n = 0, 1, 2, ....       //
//                                                                            //
//  Return Value:                                                             //
//     E[n]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n;                                                              //
//     double e_num;                                                          //
//                                                                            //
//     ( User code to set n, the index to the sequence of Euler numbers.)     //
//                                                                            //
//     if ( n >= 0 ) e_num = Euler_Number(n);                                 //
//     else printf(" Error argument, %d, is negative\n",n);                   //
////////////////////////////////////////////////////////////////////////////////

double Euler_Number( int n)
{
   long double x;

   x = xEuler_Number(n);
   if (fabsl(x) < DBL_MAX) return (double) x;
   if ( x < 0.0L) return -DBL_MAX;
   return DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// void Euler_Number_Sequence(double e[], int start, int length)              //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Euler numbers                     //
//     e[0] = E[start], e[1] = E[start+1],...,e[length-1] = E[start+length-1].//
//                                                                            //
//  Arguments:                                                                //
//     double e[]                                                             //
//        The array of length at least "length" which upon completion will    //
//        contain the Euler numbers E[start],...,E[start+length-1].           //
//     int start                                                              //
//        The starting index of the sequence of Euler numbers requested.      //
//        "start" must be nonnegative.                                        //
//     int length                                                             //
//        The number of sequential Euler numbers to be returned.              //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Euler numbers is returned in the            //
//     array e[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     double e_num[10];                                                      //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) Euler_Number_Sequence(e_num,start,10);               //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void Euler_Number_Sequence(double *e, int start, int length)
{
   long double x;
   int n;

   for (n = start; n < start+length; n++) {
      x = xEuler_Number(n);
      if (fabsl(x) < DBL_MAX) *e++ = (double) x;
      else *e++ = ( x < 0.0L) ? -DBL_MAX : DBL_MAX;
   }
}


////////////////////////////////////////////////////////////////////////////////
// void Euler_Even_Index_Sequence(double e[], int start, int length)          //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Euler numbers                     //
//     e[0] = E[start], e[1] = E[start+2],...,                                //
//                                         e[length-1] = E[start+2(length-1)].//
//                                                                            //
//  Arguments:                                                                //
//     double e[]                                                             //
//        The array of length at least "length" which upon completion will    //
//        contain the Euler numbers E[start],E[start+2]...,                   //
//        E[start+2(length-1)].                                               //
//     int start                                                              //
//        The starting index of the sequence of Euler numbers requested.      //
//        "start" must be nonnegative and generally an even number.           //
//     int length                                                             //
//        The number of Euler numbers to be returned.                         //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Euler numbers is returned in the            //
//     array e[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     double e_num[10];                                                      //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) Euler_Even_Index_Sequence(e_num,start,10);           //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void Euler_Even_Index_Sequence(double *e, int start, int length)
{
   long double x;
   int n;

   for (n = start; n < start+2*length; n += 2) {
      x = xEuler_Number(n);
      if (fabsl(x) < DBL_MAX) *e++ = (double) x;
      else *e++ = ( x < 0.0L) ? -DBL_MAX : DBL_MAX;
   }
}


////////////////////////////////////////////////////////////////////////////////
// int Max_Euler_Even_Number_Index(void)                                      //
//                                                                            //
//  Description:                                                              //
//     This routine returns the maximum even index of the Euler number        //
//     which is representable as a double.  Note that all Euler numbers       //
//     with an odd index are 0.                                               //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Value:                                                             //
//     The maximum even index of the Euler number representable as a type     //
//     double.                                                                //
//                                                                            //
//  Example:                                                                  //
//     int    big;                                                            //
//                                                                            //
//     big = Max_Euler_Even_Number_Index();                                   //
////////////////////////////////////////////////////////////////////////////////

int Max_Euler_Even_Number_Index(void)
{
   return max_double_index;
}

////////////////////////////////////////////////////////////////////////////////
// long double xEuler_Number(int n)                                           //
//                                                                            //
//  Description:                                                              //
//     This routine returns the Euler number E[n], where n >= 0.              //
//     If n > max_long_double_index, then if n is odd, 0 is returned, if      //
//     n is even and n = 0 (mod 4) then LDBL_MAX is returned and if n is      //
//     even and n = 2 (mod 4) then -LDBL_MAX is returned.                     //
//                                                                            //
//  Arguments:                                                                //
//     int n                                                                  //
//        The index of the sequence of Euler numbers, n = 0, 1, 2, ....       //
//                                                                            //
//  Return Value:                                                             //
//     E[n]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n;                                                              //
//     long double e_num;                                                     //
//                                                                            //
//     ( User code to set n, the index to the sequence of Euler numbers.)     //
//                                                                            //
//     if ( n >= 0 ) e_num = xEuler_Number(n);                                //
//     else printf(" Error argument, %d, is negative\n",n);                   //
////////////////////////////////////////////////////////////////////////////////

long double xEuler_Number( int n)
{
   long double en;
   int m = n >> 1;

   if ( (m + m)  != n) return 0.0L;
   if ( n <= max_double_index ) return euler_numbers[m];
   if (n > max_long_double_index) return (m % 2 == 0) ? LDBL_MAX : -LDBL_MAX;
   m = (n - 188) / 20;
   en = xeuler_numbers[m];
   m = 20 * m + 198;
   if (n >= m) return Iterate_Up(en, m, n);
   return Iterate_Down(en, m, n);
}


////////////////////////////////////////////////////////////////////////////////
// void xEuler_Number_Sequence(long double e[], int start, int length)        //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Euler numbers                     //
//     e[0] = E[start], e[1] = E[start+1],...,e[length-1] = E[start+length-1].//
//                                                                            //
//  Arguments:                                                                //
//     long double e[]                                                        //
//        The array of length at least "length" which upon completion will    //
//        contain the Euler numbers E[start],...,E[start+length-1].           //
//     int start                                                              //
//        The starting index of the sequence of Euler numbers requested.      //
//        "start" must be nonnegative.                                        //
//     int length                                                             //
//        The number of sequential Euler numbers to be returned.              //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Euler numbers is returned in the            //
//     array e[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     long double e_num[10];                                                 //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) xEuler_Number_Sequence(e_num,start,10);              //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void xEuler_Number_Sequence(long double *e, int start, int length)
{
   long double x;
   int n;

   for (n = start; n < start+length; n++) 
      *e++ = xEuler_Number(n);
}


////////////////////////////////////////////////////////////////////////////////
// void xEuler_Even_Index_Sequence(long double e[], int start, int length)    //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Euler numbers                     //
//     e[0] = E[start], e[1] = E[start+2],...,                                //
//                                         e[length-1] = E[start+2(length-1)].//
//                                                                            //
//  Arguments:                                                                //
//     double e[]                                                             //
//        The array of length at least "length" which upon completion will    //
//        contain the Euler numbers E[start],E[start+2]...,                   //
//        E[start+2(length-1)].                                               //
//     int start                                                              //
//        The starting index of the sequence of Euler numbers requested.      //
//        "start" must be nonnegative and generally an even number.           //
//     int length                                                             //
//        The number of Euler numbers to be returned.                         //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Euler numbers is returned in the            //
//     array e[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     long double e_num[10];                                                 //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) xEuler_Even_Index_Sequence(e_num,start,10);          //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void xEuler_Even_Index_Sequence(long double *e, int start, int length)
{
   long double x;
   int n;
   int m;

   for (n = start; n < start+2*length; n += 2) {
      *e++ = xEuler_Number(n);
   }
}


////////////////////////////////////////////////////////////////////////////////
// int xMax_Euler_Even_Number_Index(void)                                     //
//                                                                            //
//  Description:                                                              //
//     This routine returns the maximum even index of the Euler number        //
//     which is representable as a long double.  Note that all Euler          //
//     numbers with an odd index are 0.                                       //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Value:                                                             //
//     The maximum even index of the Euler number representable as a type     //
//     long double.                                                           //
//                                                                            //
//  Example:                                                                  //
//     int    big;                                                            //
//                                                                            //
//     big = xMax_Euler_Even_Number_Index();                                  //
////////////////////////////////////////////////////////////////////////////////

int xMax_Euler_Even_Number_Index(void)
{
   return max_long_double_index;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Iterate_Down(long double en, int n, int m)              //
//                                                                            //
//  Description:                                                              //
//     Starting with E[n] iterate down to E[m], n > m, using the recursion    //
//     formula valid for large even m and n:                                  //
//               E[k] = -[ (pi / 2)^2 / ((k+2) * (k+1)) ] E[k+2],             //
//     for k = n - 2,...,m.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double en                                                         //
//        The nth Euler number E[n].                                          //
//     int n                                                                  //
//        The argument of E[] corresponding the value en, n must be even.     //
//     int m                                                                  //
//        The argument of the desired Euler number E[m], m <= n, m must be    //
//        even.                                                               //
//                                                                            //
//  Return Value:                                                             //
//     E[mm]                                                                  //
//                                                                            //
//  Example:                                                                  //
//     int    n,m;                                                            //
//     long double en,em;                                                     //
//                                                                            //
//     em = Iterate_Down(en,n,m);                                             //
////////////////////////////////////////////////////////////////////////////////

static long double Iterate_Down(long double en, int n, int m)
{
   int k;

   for (k = n - 2; k >= m; k -= 2) {
      en *= -pi_over_2_squared;
      en /= ( (k + 2) * (k + 1) );
   }
   return en;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Iterate_Up(long double en, int n, int m)                //
//                                                                            //
//  Description:                                                              //
//     Starting with E[n] iterate up to E[m], n < m, using the recursion      //
//     formula valid for large even m and n:                                  //
//             E[k] = -[ k * (k-1) / (2/pi)^2 ] E[k-2],                       //
//     for k = n+2, ..., m.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double en                                                         //
//        The nth Euler number E[n].                                          //
//     int n                                                                  //
//        The argument of E[] corresponding the value en, n must be even.     //
//     int m                                                                  //
//        The argument of the desired Euler number E[m], m >= n, m must be    //
//        even.                                                               //
//                                                                            //
//  Return Value:                                                             //
//     E[m]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n,m;                                                            //
//     long double en,em;                                                     //
//                                                                            //
//              (* User code to calculate n, m and en = E[n] *)               //
//                                                                            //
//     em = Iterate_Up(en,n,m);                                               //
////////////////////////////////////////////////////////////////////////////////

static long double Iterate_Up(long double en, int n, int m)
{
   int k;

   for (k = n + 2; k <= m; k += 2) {
      en *= -( k * (k - 1) );
      en /= pi_over_2_squared;
   }
   return en;
}


