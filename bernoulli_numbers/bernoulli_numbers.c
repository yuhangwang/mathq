////////////////////////////////////////////////////////////////////////////////
// File: bernoulli_numbers.c                                                  //
// Routine(s):                                                                //
//    Bernoulli_Number                                                        //
//    Bernoulli_Number_Sequence                                               //
//    Bernoulli_Even_Index_Sequence                                           //
//    Max_Bernoulli_Even_Number_Index                                         //
//    xBernoulli_Number                                                       //
//    xBernoulli_Number_Sequence                                              //
//    xBernoulli_Even_Index_Sequence                                          //
//    xMax_Bernoulli_Even_Number_Index                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     There are several definitions of the Bernoulli numbers, B[n].  The     //
//     definition programed here yields B[0] = 1, B[1] = -1/2, B[2] = 1/6,    //
//     B[3] = 0, B[4] = -1/30, ... .  This is equivalent to using the         //
//     generating function x / (exp(x) - 1) = Sum B[n] x^n / n! where the     //
//     sum extends over n = 0,1,... .                                         //
//     Multiply both sides of generating function by (exp(x) - 1) / x         //
//     and express (exp(x) - 1) / x in terms of its power series,             //
//                      1 = Sum x^k / (k+1)! Sum B[j] x^j / j!                //
//                      1 = Sum Sum B[j] x^(j+k) / (j! (k+1)!)                //
//     Let n = j + k, then                                                    //
//                      1 = Sum Sum B[j] x^n / (j! (n - j + 1)!)              //
//     where the first sum extends from n = 0, 1, ..., and the second sum     //
//     extends from j = 0, .., n.  Equating like powers of x on both sides    //
//                      1 = B[0]                                              //
//                      0 = Sum B[j] / [j! (n-j+1)!], for n >= 1.             //
//     Therefore                                                              //
//         (1)          B[0] = 1                                              //
//                      B[n] = - n! Sum B[j] / [j! (n-j+1)!] for n >= 1       //
//     where the sum extends from j = 0,...,n-1.                              //
//                                                                            //
//     The even-indexed Bernoulli numbers also satisy the expression          //
//     involving the Riemann zeta function:                                   //
//                      B[2n] = (-1)^(n+1) [ 2 (2n)! / (2pi)^(2n) ] zeta(2n)  //
//     where zeta(2n) = 1 + (1/2)^2n + (1/3)^2n + ... .                       //
//     Given that LDBL_EPSILON ~ 10^(-20), the zeta function is               //
//     machine computationally equal to 1.0 for 2n >= 68.                     //
//     Therefore for 2n >= 68,                                                //
//         (2)          B[2n] = (-1)^(n+1) [ 2 (2n)! / (2pi)^(2n) ]           //
//     or                                                                     //
//         (3)          B[2n] = -(2n)(2n-1)/(2pi)^2 B[2n-2].                  //
//                                                                            //
//     Given that LDBL_MAX ~ 10^(4932), using equation (2), the maximum       //
//     even-indexed argument is 2312.  Note that the odd-indexed Bernoulli    //
//     numbers are 0 save B[1] = -1/2.                                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include <math.h>          // required for fabsl()
#include <float.h>         // required for DBL_MAX, LDBL_MAX

//                         Internally Defined Routines                        //

double      Bernoulli_Number( int n);
void        Bernoulli_Number_Sequence(double B[], int start, int length);
void        Bernoulli_Even_Index_Sequence(double B[], int start, int length);
int         Max_Bernoulli_Even_Number_Index(void);
long double xBernoulli_Number(int n);
void        xBernoulli_Number_Sequence(long double B[], int start, int length);
void        xBernoulli_Even_Index_Sequence(long double B[], int start,
                                                                   int length);
int         xMax_Bernoulli_Even_Number_Index(void);

static long double Iterate_Up(long double bn, int n, int m);
static long double Iterate_Down(long double bn, int n, int m);

//                         Internally Defined Constants                       //

 // bernoulli numbers even arguments B[2n] = bernoulli_number[n], n = 0,... //

static long double bernoulli_numbers[] = {
         +1.00000000000000000000000000000000000L,
         +1.66666666666666666666666666666666667e-1L,
         -3.33333333333333333333333333333333333e-2L,
         +2.38095238095238095238095238095238095e-2L,
         -3.33333333333333333333333333333333333e-2L,
         +7.57575757575757575757575757575757576e-2L,
         -2.53113553113553113553113553113553114e-1L,
         +1.16666666666666666666666666666666667e+0L,
         -7.09215686274509803921568627450980392e+0L,
         +5.49711779448621553884711779448621554e+1L,
         -5.29124242424242424242424242424242424e+2L,
         +6.19212318840579710144927536231884058e+3L,
         -8.65802531135531135531135531135531136e+4L,
         +1.42551716666666666666666666666666667e+6L,
         -2.72982310678160919540229885057471264e+7L,
         +6.01580873900642368384303868174835917e+8L,
         -1.51163157670921568627450980392156863e+10L,
         +4.29614643061166666666666666666666667e+11L,
         -1.37116552050883327721590879485616328e+13L,
         +4.88332318973593166666666666666666667e+14L,
         -1.92965793419400681486326681448632668e+16L,
         +8.41693047573682615000553709856035437e+17L,
         -4.03380718540594554130768115942028986e+19L,
         +2.11507486380819916056014539007092199e+21L,
         -1.20866265222965259346027311937082525e+23L,
         +7.50086674607696436685572007575757576e+24L,
         -5.03877810148106891413789303052201258e+26L,
         +3.65287764848181233351104308429711779e+28L,
         -2.84987693024508822262691464329106782e+30L,
         +2.38654274996836276446459819192192150e+32L,
         -2.13999492572253336658107447651910974e+34L,
         +2.05009757234780975699217330956723103e+36L,
         -2.09380059113463784090951852900279702e+38L,
         +2.27526964884635155596492603527692646e+40L,
         -2.62577102862395760473030497361582021e+42L,
         +3.21250821027180325182047923042649852e+44L,
         -4.15982781667947109139170744952623589e+46L,
         +5.69206954820352800238834562191210586e+48L,
         -8.21836294197845756922906534686173330e+50L,
         +1.25029043271669930167323398297028955e+53L,
         -2.00155832332483702749253291988132988e+55L,
         +3.36749829153643742333966769033387530e+57L,
         -5.94709705031354477186604968440515408e+59L,
         +1.10119103236279775595641307904376916e+62L,
         -2.13552595452535011886583850190410657e+64L,
         +4.33288969866411924196166130593792062e+66L,
         -9.18855282416693282262005552155018971e+68L,
         +2.03468967763290744934550279902200201e+71L,
         -4.70038339580357310785752555350060607e+73L,
         +1.13180434454842492706751862577339343e+76L,
         -2.83822495706937069592641563364817647e+78L,
         +7.40642489796788506297508271409209842e+80L,
         -2.00964548027566044834656196727153632e+83L,
         +5.66571700508059414457193460305193570e+85L,
         -1.65845111541362169158237133743199123e+88L,
         +5.03688599504923774192894219151801548e+90L,
         -1.58614682376581863693634015729664388e+93L,
         +5.17567436175456269840732406825071226e+95L,
         -1.74889218402171173396900258776181591e+98L,
         +6.11605199949521852558245252642641678e+100L,
         -2.21227769127078349422883234567129324e+103L,
         +8.27227767987709698542210624599845957e+105L,
         -3.19589251114157095835916343691808149e+108L,
         +1.27500822233877929823100243029266799e+111L,
         -5.25009230867741338994028246245651754e+113L,
         +2.23018178942416252098692981988387281e+116L,
         -9.76845219309552044386335133989802393e+118L,
         +4.40983619784529542722726228748131692e+121L,
         -2.05085708864640888397293377275830155e+124L,
         +9.82144332797912771075729696020975210e+126L,
         -4.84126007982088805087891967099634128e+129L,
         +2.45530888014809826097834674040886904e+132L,
         -1.28069268040847475487825132786017857e+135L,
         +6.86761671046685811921018885984644004e+137L,
         -3.78464685819691046949789954163795568e+140L,
         +2.14261012506652915508713231351482721e+143L,
         -1.24567271371836950070196429616376072e+146L,
         +7.43457875510001525436796683940520613e+148L,
         -4.55357953046417048940633332233212749e+151L,
         +2.86121128168588683453638472510172325e+154L,
         -1.84377235520338697276882026536287855e+157L,
         +1.21811545362210466995013165065995214e+160L,
         -8.24821871853141215484818457296893447e+162L,
         +5.72258779378329433296516498142978616e+165L,
         -4.06685305250591047267679693831158656e+168L,
         +2.95960920646420500628752695815851870e+171L,
         -2.20495225651894575090311752273445985e+174L,
         +1.68125970728895998058311525151360666e+177L,
         -1.31167362135569576486452806355817153e+180L,
         +1.04678940094780380821832853929823090e+183L,
         -8.54328935788337077185982546299082775e+185L,
         +7.12878213224865423522884066771438225e+188L,
         -6.08029314555358993000847118686477458e+191L,
         +5.29967764248499239300942910043247266e+194L,
         -4.71942591687458626443646229013379911e+197L,
         +4.29284137914029810894168296541074669e+200L,
         -3.98767449682322074434477655542938795e+203L,
         +3.78197804193588827138944181161393328e+206L,
         -3.66142336836811912436858082151197349e+209L,
         +3.61760902723728623488554609298914089e+212L,
         -3.64707726451913543621383088655499449e+215L,
         +3.75087554364544090983452410104814189e+218L,
         -3.93458672964390282694891288533713429e+221L,
         +4.20882111481900820046571171111494898e+224L,
         -4.59022962206179186559802940573325591e+227L,
         +5.10317257726295759279198185106496769e+230L,
         -5.78227623036569554015377271242917143e+233L,
         +6.67624821678358810322637794412809363e+236L,
         -7.85353076444504163225916259639312444e+239L,
         +9.41068940670587255245443288258762485e+242L,
         -1.14849338734651839938498599206805593e+246L,
         +1.42729587428487856771416320087122500e+249L,
         -1.80595595869093090142285728117654561e+252L,
         +2.32615353076608052161297985184708876e+255L,
         -3.04957517154995947681942819261542594e+258L,
         +4.06858060764339734424012124124937319e+261L,
         -5.52310313219743616252320044093186392e+264L,
         +7.62772793964343924869949690204961216e+267L,
         -1.07155711196978863132793524001065397e+271L,
         +1.53102008959691884453440916153355334e+274L,
         -2.22448916821798346676602348865048511e+277L,
         +3.28626791906901391668189736436895275e+280L,
         -4.93559289559603449020711938191575963e+283L,
         +7.53495712008325067212266049779283957e+286L,
         -1.16914851545841777278088924731655042e+290L,
         +1.84352614678389394126646201597702232e+293L,
         -2.95368261729680829728014917350525183e+296L,
         +4.80793212775015697668878704043264072e+299L,
         -7.95021250458852528538243631671158693e+302L,
         +1.33527841873546338750122832017820518e+306L};

 // bernoulli numbers even arguments B[20n+270] = xbernoulli_number[n],//
 // n = 0,...,102. //

static long double xbernoulli_numbers[] = {
         +4.13121317607384235973251163948966940e+325L,
         +4.06725630354221225869883600368201604e+358L,
         +1.58852491244122147281469212106982170e+392L,
         +2.25191059133671680915395814672577572e+426L,
         +1.07164338264967557208686546587391661e+461L,
         +1.59752224396858654822751463995972770e+496L,
         +7.01366744280728845244177798142505561e+531L,
         +8.58020723503261785605925064309501976e+567L,
         +2.78229778527875642617754227085498409e+604L,
         +2.28549456528753068146575779851703354e+641L,
         +4.56346231319052136323518242017878446e+678L,
         +2.13274737136019050759574844453691108e+716L,
         +2.25344601183435273327994630683594073e+754L,
         +5.21352427258719957498011735101632252e+792L,
         +2.56421038571922400015654824093410897e+831L,
         +2.60860868193932239358106918827162612e+870L,
         +5.35090808971096424467133422470805781e+909L,
         +2.16117469769979326571518209176467667e+949L,
         +1.68094360014785832214876780698752741e+989L,
         +2.46596231224141873152897352659743310e+1029L,
         +6.69136133257633373813072061684170699e+1069L,
         +3.29737130784864316153222745990138673e+1110L,
         +2.90027368880998769412885765503678326e+1151L,
         +4.47964020731247709293854154677691596e+1192L,
         +1.19641999974741190914414431549965447e+1234L,
         +5.44528951863630699294260477558597778e+1275L,
         +4.16528478663265416856309685061018538e+1317L,
         +5.28504512512583234126389723340519681e+1359L,
         +1.09852143386029963348144968536491412e+1402L,
         +3.69624802581714469084053913210353883e+1444L,
         +1.99061011272471512689500879301421451e+1487L,
         +1.69742181279420379386503220619132270e+1530L,
         +2.26823825575142131255980612298093295e+1573L,
         +4.70323556754388815204962841135454251e+1616L,
         +1.49903614447306459330826068178204826e+1660L,
         +7.27789175914272529417292625836445594e+1703L,
         +5.33592819551240570973377164238950281e+1747L,
         +5.85882068366170855342236377741943082e+1791L,
         +9.55728102705697044611754198378566030e+1835L,
         +2.29850417105056075619235210606259864e+1880L,
         +8.08971060455738243016203150276177139e+1924L,
         +4.13720796521752041053050805386375922e+1969L,
         +3.05345614597816164582345471073790450e+2014L,
         +3.23084795272485636694393980424818620e+2059L,
         +4.86983141221469211917289582228508416e+2104L,
         +1.03923132885160922482233503943089864e+2150L,
         +3.12126806833819945889176493238481974e+2195L,
         +1.31182768398402511174435834778399634e+2241L,
         +7.67252228080694439735866856637964654e+2286L,
         +6.21129738133960687706282445974212906e+2332L,
         +6.92390450563330178846148278663422074e+2378L,
         +1.05744378505391541199102941007672202e+2425L,
         +2.20181928480795405509211770603311317e+2471L,
         +6.22113463665521369604174068513122400e+2517L,
         +2.37425979196319369383757678132139174e+2564L,
         +1.21849507051854920806654411173698559e+2611L,
         +8.37295691983238673049007062562278548e+2657L,
         +7.67132053078199935920009773995131623e+2704L,
         +9.33310084393021656789450800715864493e+2751L,
         +1.50184202300344959033799790094592416e+2799L,
         +3.18412722947632272773220801727926821e+2846L,
         +8.86125834394592566727216453150426569e+2893L,
         +3.22517948915495742349290595788774412e+2941L,
         +1.52975735756834262991256082724328206e+2989L,
         +9.42319952595478495553395998127899248e+3036L,
         +7.51306558344496470479570706050116162e+3084L,
         +7.72773120824010179184551559965944156e+3132L,
         +1.02214342027002972168255108491773037e+3181L,
         +1.73316610732037731038876504765998784e+3229L,
         +3.75589858990025479589481294227571184e+3277L,
         +1.03715034650505289230207763788352270e+3326L,
         +3.63884496365180016808062351190070504e+3374L,
         +1.61751744455702025086491965530118919e+3423L,
         +9.08440236815702565843030025224652660e+3471L,
         +6.42882906445293964054147519865556089e+3520L,
         +5.71753073427761194916291733781074992e+3569L,
         +6.37389859829851340022881911319772874e+3618L,
         +8.88435190810858155115125256646660613e+3667L,
         +1.54453758058234789298017795698410121e+3717L,
         +3.34098802417699522351564081593703704e+3766L,
         +8.97076172060559176230095800755753387e+3815L,
         +2.98303173866272791201688239951587912e+3865L,
         +1.22568255729602707502702153496002615e+3915L,
         +6.20908659583348770719249208717684323e+3964L,
         +3.86958991363574575878627523129665292e+4014L,
         +2.96051276191437626818506412960054931e+4064L,
         +2.77478724685634765115227807646604314e+4114L,
         +3.17956241012382087578705283397501097e+4164L,
         +4.44540885470315644057680807036093474e+4214L,
         +7.56853185420275088133874643207881721e+4264L,
         +1.56614742596747185173656286731874851e+4315L,
         +3.93148221264316732379836632739005868e+4365L,
         +1.19503310506395297988508675434270665e+4416L,
         +4.39052031053386419818620236802663043e+4466L,
         +1.94621918590048211613785506477563525e+4517L,
         +1.03907969660921565101173608760330477e+4568L,
         +6.67026273060300930659501812225273074e+4618L,
         +5.13976436809289046623516243179535059e+4669L,
         +4.74604776985189143857600204752925811e+4720L,
         +5.24330276965152053675952126461515991e+4771L,
         +6.91941751868863603233513125358433165e+4822L,
         +1.08904527435538965461419665176131097e+4874L,
         +2.04110823109132319887750995937125750e+4925L};

static long double b1 = -5.000000000000000000000e-1L;
static long double two_pi_squared =39.4784176043574344753379639995046045L;  
static int max_double_index = 2*sizeof(bernoulli_numbers) / sizeof(long double)
                              - 2;
static int max_long_double_index = 2312;

////////////////////////////////////////////////////////////////////////////////
// double Bernoulli_Number(int n)                                             //
//                                                                            //
//  Description:                                                              //
//     This routine returns the Bernoulli number B[n], where n >= 0.          //
//     If B[n] > DBL_MAX, then DBL_MAX is returned, if B[n] < -DBL_MAX, then  //
//     -DBL_MAX is returned.                                                  //
//                                                                            //
//  Arguments:                                                                //
//     int n                                                                  //
//        The index of the sequence of Bernoulli numbers, n = 0, 1, 2, ....   //
//                                                                            //
//  Return Value:                                                             //
//     B[n]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n;                                                              //
//     double b_num;                                                          //
//                                                                            //
//     ( User code to set n, the index to the sequence of Bernoulli numbers.) //
//                                                                            //
//     if ( n >= 0 ) b_num = Bernoulli_Number(n);                             //
//     else printf(" Error argument, %d, is negative\n",n);                   //
////////////////////////////////////////////////////////////////////////////////

double Bernoulli_Number( int n)
{
   long double x;

   x = xBernoulli_Number(n);
   if (fabsl(x) < DBL_MAX) return (double) x;
   if ( x < 0.0L) return -DBL_MAX;
   return DBL_MAX;
}


////////////////////////////////////////////////////////////////////////////////
// void Bernoulli_Number_Sequence(double b[], int start, int length)          //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Bernoulli numbers                 //
//     b[0] = B[start], b[1] = B[start+1],...,b[length-1] = B[start+length-1].//
//                                                                            //
//  Arguments:                                                                //
//     double b[]                                                             //
//        The array of length at least "length" which upon completion will    //
//        contain the Bernoulli numbers B[start],...,B[start+length-1].       //
//     int start                                                              //
//        The starting index of the sequence of Bernoulli numbers requested.  //
//        "start" must be nonnegative.                                        //
//     int length                                                             //
//        The number of sequential Bernoulli numbers to be returned.          //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Bernoulli numbers is returned in the        //
//     array b[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     double b_num[10];                                                      //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) Bernoulli_Number_Sequence(b_num,start,10);           //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void Bernoulli_Number_Sequence(double *b, int start, int length)
{
   long double x;
   int n;

   for (n = start; n < start+length; n++) {
      x = xBernoulli_Number(n);
      if (fabsl(x) < DBL_MAX) *b++ = (double) x;
      else *b++ = ( x < 0.0L) ? -DBL_MAX : DBL_MAX;
   }
}


////////////////////////////////////////////////////////////////////////////////
// void Bernoulli_Even_Index_Sequence(double b[], int start, int length)      //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Bernoulli numbers                 //
//     b[0] = B[start], b[1] = B[start+2],...,                                //
//                                         b[length-1] = B[start+2(length-1)].//
//                                                                            //
//  Arguments:                                                                //
//     double b[]                                                             //
//        The array of length at least "length" which upon completion will    //
//        contain the Bernoulli numbers B[start],B[start+2]...,               //
//        B[start+2(length-1)].                                               //
//     int start                                                              //
//        The starting index of the sequence of Bernoulli numbers requested.  //
//        "start" must be nonnegative and generally an even number.           //
//     int length                                                             //
//        The number of Bernoulli numbers to be returned.                     //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Bernoulli numbers is returned in the        //
//     array b[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     double b_num[10];                                                      //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) Bernoulli_Even_Index_Sequence(b_num,start,10);       //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void Bernoulli_Even_Index_Sequence(double *b, int start, int length)
{
   long double x;
   int n;

   for (n = start; n < start+2*length; n += 2) {
      x = xBernoulli_Number(n);
      if (fabsl(x) < DBL_MAX) *b++ = (double) x;
      else *b++ = ( x < 0.0L) ? -DBL_MAX : DBL_MAX;
   }
}


////////////////////////////////////////////////////////////////////////////////
// int Max_Bernoulli_Even_Number_Index(void)                                  //
//                                                                            //
//  Description:                                                              //
//     This routine returns the maximum even index of the Bernoulli number    //
//     which is representable as a double.  Note that all Bernoulli numbers   //
//     with an odd index greater than 1, are 0.                               //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Value:                                                             //
//     The maximum even index of the Bernoulli number representable as a type //
//     double.                                                                //
//                                                                            //
//  Example:                                                                  //
//     int    big;                                                            //
//                                                                            //
//     big = Max_Bernoulli_Even_Number_Index();                               //
////////////////////////////////////////////////////////////////////////////////

int Max_Bernoulli_Even_Number_Index(void)
{
   return max_double_index;
}


////////////////////////////////////////////////////////////////////////////////
// long double xBernoulli_Number(int n)                                       //
//                                                                            //
//  Description:                                                              //
//     This routine returns the Bernoulli number B[n], where n >= 0.          //
//     If n > max_long_double_index, then if n is odd, 0 is returned, if      //
//     n is even and n = 0 (mod 4) then -LDBL_MAX is returned and if n is     //
//     even and n = 2 (mod 4) then LDBL_MAX is returned.                      //
//                                                                            //
//  Arguments:                                                                //
//     int n                                                                  //
//        The index of the sequence of Bernoulli numbers, n = 0, 1, 2, ....   //
//                                                                            //
//  Return Value:                                                             //
//     B[n]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n;                                                              //
//     long double b_num;                                                     //
//                                                                            //
//     ( User code to set n, the index to the sequence of Bernoulli numbers.) //
//                                                                            //
//     if ( n >= 0 ) b_num = xBernoulli_Number(n);                            //
//     else printf(" Error argument, %d, is negative\n",n);                   //
////////////////////////////////////////////////////////////////////////////////

long double xBernoulli_Number( int n)
{
   long double bn;
   int m = n >> 1;

   if ( (m + m)  != n) return (n == 1) ? b1 : 0.0L;
   if ( n <= max_double_index ) return bernoulli_numbers[m];
   if (n > max_long_double_index) return (m % 2 == 1) ? LDBL_MAX : -LDBL_MAX;
   m = (n - 260) / 20;
   bn = xbernoulli_numbers[m];
   m = 20 * m + 270;
   if (n >= m) return Iterate_Up(bn, m, n);
   return Iterate_Down(bn, m, n);
}


////////////////////////////////////////////////////////////////////////////////
// void xBernoulli_Number_Sequence(long double b[], int start, int length)    //
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Bernoulli numbers                 //
//     b[0] = B[start], b[1] = B[start+1],...,b[length-1] = B[start+length-1].//
//                                                                            //
//  Arguments:                                                                //
//     long double b[]                                                        //
//        The array of length at least "length" which upon completion will    //
//        contain the Bernoulli numbers B[start],...,B[start+length-1].       //
//     int start                                                              //
//        The starting index of the sequence of Bernoulli numbers requested.  //
//        "start" must be nonnegative.                                        //
//     int length                                                             //
//        The number of sequential Bernoulli numbers to be returned.          //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Bernoulli numbers is returned in the        //
//     array b[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     long double b_num[10];                                                 //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) xBernoulli_Number_Sequence(b_num,start,10);          //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void xBernoulli_Number_Sequence(long double *b, int start, int length)
{
   int n;

   for (n = start; n < start+length; n++) 
      *b++ = xBernoulli_Number(n);
}


////////////////////////////////////////////////////////////////////////////////
// void xBernoulli_Even_Index_Sequence(long double b[], int start, int length)//
//                                                                            //
//  Description:                                                              //
//     This routine returns the sequence of Bernoulli numbers                 //
//     b[0] = B[start], b[1] = B[start+2],...,                                //
//                                         b[length-1] = B[start+2(length-1)].//
//                                                                            //
//  Arguments:                                                                //
//     long double b[]                                                        //
//        The array of length at least "length" which upon completion will    //
//        contain the Bernoulli numbers B[start],B[start+2]...,               //
//        B[start+2(length-1)].                                               //
//     int start                                                              //
//        The starting index of the sequence of Bernoulli numbers requested.  //
//        "start" must be nonnegative and generally an even number.           //
//     int length                                                             //
//        The number of Bernoulli numbers to be returned.                     //
//        "length" must be positive.                                          //
//                                                                            //
//  Return Value:                                                             //
//     type void, the sequence of Bernoulli numbers is returned in the        //
//     array b[].                                                             //
//                                                                            //
//  Example:                                                                  //
//     int    start;                                                          //
//     long double b_num[10];                                                 //
//                                                                            //
//     ( User code to set start.)                                             //
//                                                                            //
//     if ( start >= 0 ) xBernoulli_Even_Index_Sequence(b_num,start,10);      //
//     else printf(" Error argument, %d, is negative\n",start);               //
////////////////////////////////////////////////////////////////////////////////

void xBernoulli_Even_Index_Sequence(long double *b, int start, int length)
{
   long double x;
   int n;
   int m;

   for (n = start; n < start+2*length; n += 2) {
      *b++ = xBernoulli_Number(n);
   }
}


////////////////////////////////////////////////////////////////////////////////
// int xMax_Bernoulli_Even_Number_Index(void)                                 //
//                                                                            //
//  Description:                                                              //
//     This routine returns the maximum even index of the Bernoulli number    //
//     which is representable as a long double.  Note that all Bernoulli      //
//     numbers with an odd index greater than 1, are 0.                       //
//                                                                            //
//  Arguments:                                                                //
//                                                                            //
//  Return Value:                                                             //
//     The maximum even index of the Bernoulli number representable as a type //
//     long double.                                                           //
//                                                                            //
//  Example:                                                                  //
//     int    big;                                                            //
//                                                                            //
//     big = xMax_Bernoulli_Even_Number_Index();                              //
////////////////////////////////////////////////////////////////////////////////

int xMax_Bernoulli_Even_Number_Index(void)
{
   return max_long_double_index;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Iterate_Down(long double bn, int n, int m)              //
//                                                                            //
//  Description:                                                              //
//     Starting with B[n] iterate down to B[m], n > m, using the recursion    //
//     formula valid for large even m and n:                                  //
//               B[k] = -[ (2 pi)^2 / ((k+2) * (k+1)) ] B[k+2],               //
//     for k = n - 2,...,m.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double bn                                                         //
//        The nth Bernoulli number B[n].                                      //
//     int n                                                                  //
//        The argument of B[] corresponding the value bn, n must be even.     //
//     int m                                                                  //
//        The argument of the desired Bernoulli number B[m], m <= n, m must be//
//        even.                                                               //
//                                                                            //
//  Return Value:                                                             //
//     B[mm]                                                                  //
//                                                                            //
//  Example:                                                                  //
//     int    n,m;                                                            //
//     long double bn,bm;                                                     //
//                                                                            //
//     bm = Iterate_Down(bn,n,m);                                             //
////////////////////////////////////////////////////////////////////////////////

static long double Iterate_Down(long double bn, int n, int m)
{
   int k;

   for (k = n - 2; k >= m; k -= 2) {
      bn *= -two_pi_squared;
      bn /= ( (k + 2) * (k + 1) );
   }
   return bn;
}


////////////////////////////////////////////////////////////////////////////////
// static long double Iterate_Up(long double bn, int n, int m)                //
//                                                                            //
//  Description:                                                              //
//     Starting with B[n] iterate up to B[m], n < m, using the recursion      //
//     formula valid for large even m and n:                                  //
//             B[k] = -[ k * (k-1) / (2 pi)^2 ] B[k-2],                       //
//     for k = n+2, ..., m.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     long double bn                                                         //
//        The nth Bernoulli number B[n].                                      //
//     int n                                                                  //
//        The argument of B[] corresponding the value bn, n must be even.     //
//     int m                                                                  //
//        The argument of the desired Bernoulli number B[m], m >= n, m must be//
//        even.                                                               //
//                                                                            //
//  Return Value:                                                             //
//     B[m]                                                                   //
//                                                                            //
//  Example:                                                                  //
//     int    n,m;                                                            //
//     long double bn,bm;                                                     //
//                                                                            //
//              (* User code to calculate n, m and bn = B[n] *)               //
//                                                                            //
//     bm = Iterate_Up(bn,n,m);                                               //
////////////////////////////////////////////////////////////////////////////////

static long double Iterate_Up(long double bn, int n, int m)
{
   int k;

   for (k = n + 2; k <= m; k += 2) {
      bn *= -( k * (k - 1) );
      bn /= two_pi_squared;
   }
   return bn;
}


