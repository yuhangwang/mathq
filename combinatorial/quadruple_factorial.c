////////////////////////////////////////////////////////////////////////////////
// File: quadruple_factorial.c                                                //
// Routine(s):                                                                //
//    Quadruple_Factorial                                                     //
//    xQuadruple_Factorial                                                    //
//    Quadruple_Factorial_Max_Arg                                             //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//     The quadruple factorial n!!!! is defined as:                           //
//                   n!!!! = n (n - 4) ... 4,  if n = 0 (mod 4)               //
//                   n!!!! = n (n - 4) ... 3,  if n = 3 (mod 4)               //
//                   n!!!! = n (n - 4) ... 2,  if n = 2 (mod 4)               //
//                   n!!!! = n (n - 4) ... 1,  if n = 1 (mod 4),              //
//     where by convention (-3)!!!! = (-2)!!!! = (-1)!!!! = 0!!!! = 1.        //
//                                                                            //
//     The functions Quadruple_Factorial(n) and xQuadruple_Factorial(n) return//
//     n!!!!  for -3 <= n <= Quadruple_Factorial_Max_Arg().                   //
//     If n > Quadruple_Factorial_Max_Arg(), then DBL_MAX is returned and     //
//     The function Quadruple_Factorial_Max_Arg() returns the maximum argument//
//     to the functions Quadruple_Factorial() and xQuadruple_Factorial().     //
////////////////////////////////////////////////////////////////////////////////

#include <float.h>                         // required for DBL_MAX

//                         Internally Defined Routines                        //

double Quadruple_Factorial(int n);
long double xQuadruple_Factorial(int n);
int Quadruple_Factorial_Max_Arg( void );

//                         Internally Defined Constants                       //

static long double const quadruple_factorials[] = {
       1.000000000000000000000e+0L,          //   0!!!!
       1.000000000000000000000e+0L,          //   1!!!!
       2.000000000000000000000e+0L,          //   2!!!!
       3.000000000000000000000e+0L,          //   3!!!!
       4.000000000000000000000e+0L,          //   4!!!!
       5.000000000000000000000e+0L,          //   5!!!!
       1.200000000000000000000e+1L,          //   6!!!!
       2.100000000000000000000e+1L,          //   7!!!!
       3.200000000000000000000e+1L,          //   8!!!!
       4.500000000000000000000e+1L,          //   9!!!!
       1.200000000000000000000e+2L,          //  10!!!!
       2.310000000000000000000e+2L,          //  11!!!!
       3.840000000000000000000e+2L,          //  12!!!!
       5.850000000000000000000e+2L,          //  13!!!!
       1.680000000000000000000e+3L,          //  14!!!!
       3.465000000000000000000e+3L,          //  15!!!!
       6.144000000000000000000e+3L,          //  16!!!!
       9.945000000000000000000e+3L,          //  17!!!!
       3.024000000000000000000e+4L,          //  18!!!!
       6.583500000000000000000e+4L,          //  19!!!!
       1.228800000000000000000e+5L,          //  20!!!!
       2.088450000000000000000e+5L,          //  21!!!!
       6.652800000000000000000e+5L,          //  22!!!!
       1.514205000000000000000e+6L,          //  23!!!!
       2.949120000000000000000e+6L,          //  24!!!!
       5.221125000000000000000e+6L,          //  25!!!!
       1.729728000000000000000e+7L,          //  26!!!!
       4.088353500000000000000e+7L,          //  27!!!!
       8.257536000000000000000e+7L,          //  28!!!!
       1.514126250000000000000e+8L,          //  29!!!!
       5.189184000000000000000e+8L,          //  30!!!!
       1.267389585000000000000e+9L,          //  31!!!!
       2.642411520000000000000e+9L,          //  32!!!!
       4.996616625000000000000e+9L,          //  33!!!!
       1.764322560000000000000e+10L,         //  34!!!!
       4.435863547500000000000e+10L,         //  35!!!!
       9.512681472000000000000e+10L,         //  36!!!!
       1.848748151250000000000e+11L,         //  37!!!!
       6.704425728000000000000e+11L,         //  38!!!!
       1.729986783525000000000e+12L,         //  39!!!!
       3.805072588800000000000e+12L,         //  40!!!!
       7.579867420125000000000e+12L,         //  41!!!!
       2.815858805760000000000e+13L,         //  42!!!!
       7.438943169157500000000e+13L,         //  43!!!!
       1.674231939072000000000e+14L,         //  44!!!!
       3.410940339056250000000e+14L,         //  45!!!!
       1.295295050649600000000e+15L,         //  46!!!!
       3.496303289504025000000e+15L,         //  47!!!!
       8.036313307545600000000e+15L,         //  48!!!!
       1.671360766137562500000e+16L,         //  49!!!!
       6.476475253248000000000e+16L,         //  50!!!!
       1.783114677647052750000e+17L,         //  51!!!!
       4.178882919923712000000e+17L,         //  52!!!!
       8.858212060529081250000e+17L,         //  53!!!!
       3.497296636753920000000e+18L,         //  54!!!!
       9.807130727058790125000e+18L,         //  55!!!!
       2.340174435157278720000e+19L,         //  56!!!!
       5.049180874501576312500e+19L,         //  57!!!!
       2.028432049317273600000e+20L,         //  58!!!!
       5.786207128964686173750e+20L,         //  59!!!!
       1.404104661094367232000e+21L,         //  60!!!!
       3.080000333445961550625e+21L,         //  61!!!!
       1.257627870576709632000e+22L,         //  62!!!!
       3.645310491247752289463e+22L,         //  63!!!!
       8.986269831003950284800e+22L,         //  64!!!!
       2.002000216739875007906e+23L,         //  65!!!!
       8.300343945806283571200e+23L,         //  66!!!!
       2.442358029135994033940e+24L,         //  67!!!!
       6.110663485082686193664e+24L,         //  68!!!!
       1.381380149550513755455e+25L,         //  69!!!!
       5.810240762064398499840e+25L,         //  70!!!!
       1.734074200686555764097e+26L,         //  71!!!!
       4.399677709259534059438e+26L,         //  72!!!!
       1.008407509171875041482e+27L,         //  73!!!!
       4.299578163927654889882e+27L,         //  74!!!!
       1.300555650514916823073e+28L,         //  75!!!!
       3.343755059037245885173e+28L,         //  76!!!!
       7.764737820623437819414e+28L,         //  77!!!!
       3.353670967863570814108e+29L,         //  78!!!!
       1.027438963906784290228e+30L,         //  79!!!!
       2.675004047229796708138e+30L,         //  80!!!!
       6.289437634704984633726e+30L,         //  81!!!!
       2.750010193648128067568e+31L,         //  82!!!!
       8.527743400426309608890e+31L,         //  83!!!!
       2.247003399673029234836e+32L,         //  84!!!!
       5.346021989499236938667e+32L,         //  85!!!!
       2.365008766537390138109e+33L,         //  86!!!!
       7.419136758370889359734e+33L,         //  87!!!!
       1.977362991712265726656e+34L,         //  88!!!!
       4.757959570654320875413e+34L,         //  89!!!!
       2.128507889883651124298e+35L,         //  90!!!!
       6.751414450117509317358e+35L,         //  91!!!!
       1.819173952375284468523e+36L,         //  92!!!!
       4.424902400708518414134e+36L,         //  93!!!!
       2.000797416490632056840e+37L,         //  94!!!!
       6.413843727611633851490e+37L,         //  95!!!!
       1.746406994280273089782e+38L,         //  96!!!!
       4.292155328687262861710e+38L,         //  97!!!!
       1.960781468160819415703e+39L,         //  98!!!!
       6.349705290335517512975e+39L,         //  99!!!!
       1.746406994280273089782e+40L,         // 100!!!!
       4.335076881974135490328e+40L,         // 101!!!!
       1.999997097524035804017e+41L,         // 102!!!!
       6.540196449045583038364e+41L,         // 103!!!!
       1.816263274051484013374e+42L,         // 104!!!!
       4.551830726072842264844e+42L,         // 105!!!!
       2.119996923375477952258e+43L,         // 106!!!!
       6.998010200478773851050e+43L,         // 107!!!!
       1.961564335975602734444e+44L,         // 108!!!!
       4.961495491419398068680e+44L,         // 109!!!!
       2.331996615713025747484e+45L,         // 110!!!!
       7.767791322531438974665e+45L,         // 111!!!!
       2.196952056292675062577e+46L,         // 112!!!!
       5.606489905303919817608e+46L,         // 113!!!!
       2.658476141912849352132e+47L,         // 114!!!!
       8.932960020911154820865e+47L,         // 115!!!!
       2.548464385299503072589e+48L,         // 116!!!!
       6.559593189205586186602e+48L,         // 117!!!!
       3.137001847457162235516e+49L,         // 118!!!!
       1.063022242488427423683e+50L,         // 119!!!!
       3.058157262359403687107e+50L,         // 120!!!!
       7.937107758938759285788e+50L,         // 121!!!!
       3.827142253897737927329e+51L,         // 122!!!!
       1.307517358260765731130e+52L,         // 123!!!!
       3.792115005325660572013e+52L,         // 124!!!!
       9.921384698673449107235e+52L,         // 125!!!!
       4.822199239911149788435e+53L,         // 126!!!!
       1.660547044991172478535e+54L,         // 127!!!!
       4.853907206816845532176e+54L,         // 128!!!!
       1.279858626128874934833e+55L,         // 129!!!!
       6.268859011884494724965e+55L,         // 130!!!!
       2.175316628938435946881e+56L,         // 131!!!!
       6.407157512998236102473e+56L,         // 132!!!!
       1.702211972751403663328e+57L,         // 133!!!!
       8.400271075925222931453e+57L,         // 134!!!!
       2.936677449066888528289e+58L,         // 135!!!!
       8.713734217677601099363e+58L,         // 136!!!!
       2.332030402669423018760e+59L,         // 137!!!!
       1.159237408477680764541e+60L,         // 138!!!!
       4.081981654202975054322e+60L,         // 139!!!!
       1.219922790474864153911e+61L,         // 140!!!!
       3.288162867763886456451e+61L,         // 141!!!!
       1.646117120038306685648e+62L,         // 142!!!!
       5.837233765510254327681e+62L,         // 143!!!!
       1.756688818283804381632e+63L,         // 144!!!!
       4.767836158257635361854e+63L,         // 145!!!!
       2.403330995255927761045e+64L,         // 146!!!!
       8.580733635300073861691e+64L,         // 147!!!!
       2.599899451060030484815e+65L,         // 148!!!!
       7.104075875803876689163e+65L,         // 149!!!!
       3.604996492883891641568e+66L,         // 150!!!!
       1.295690778930311153115e+67L,         // 151!!!!
       3.951847165611246336918e+67L,         // 152!!!!
       1.086923608997993133442e+68L,         // 153!!!!
       5.551694599041193128015e+68L,         // 154!!!!
       2.008320707341982287329e+69L,         // 155!!!!
       6.164881578353544285593e+69L,         // 156!!!!
       1.706470066126849219504e+70L,         // 157!!!!
       8.771677466485085142264e+70L,         // 158!!!!
       3.193229924673751836853e+71L,         // 159!!!!
       9.863810525365670856948e+71L,         // 160!!!!
       2.747416806464227243401e+72L,         // 161!!!!
       1.421011749570583793047e+73L,         // 162!!!!
       5.204964777218215494070e+73L,         // 163!!!!
       1.617664926159970020540e+74L,         // 164!!!!
       4.533237730665974951612e+74L,         // 165!!!!
       2.358879504287169096458e+75L,         // 166!!!!
       8.692291177954419875097e+75L,         // 167!!!!
       2.717677075948749634506e+76L,         // 168!!!!
       7.661171764825497668224e+76L,         // 169!!!!
       4.010095157288187463978e+77L,         // 170!!!!
       1.486381791430205798642e+78L,         // 171!!!!
       4.674404570631849371351e+78L,         // 172!!!!
       1.325382715314811096603e+79L,         // 173!!!!
       6.977565573681446187321e+79L,         // 174!!!!
       2.601168135002860147623e+80L,         // 175!!!!
       8.226952044312054893578e+80L,         // 176!!!!
       2.345927406107215640987e+81L,         // 177!!!!
       1.242006672115297421343e+82L,         // 178!!!!
       4.656090961655119664245e+82L,         // 179!!!!
       1.480851367976169880844e+83L,         // 180!!!!
       4.246128605054060310186e+83L,         // 181!!!!
       2.260452143249841306845e+84L,         // 182!!!!
       8.520646459828868985568e+84L,         // 183!!!!
       2.724766517076152580753e+85L,         // 184!!!!
       7.855337919350011573845e+85L,         // 185!!!!
       4.204440986444704830731e+86L,         // 186!!!!
       1.593360887987998500301e+87L,         // 187!!!!
       5.122561052103166851816e+87L,         // 188!!!!
       1.484658866757152187457e+88L,         // 189!!!!
       7.988437874244939178389e+88L,         // 190!!!!
       3.043319296057077135575e+89L,         // 191!!!!
       9.835317220038080355486e+89L,         // 192!!!!
       2.865391612841303721791e+90L,         // 193!!!!
       1.549756947603518200607e+91L,         // 194!!!!
       5.934472627311300414372e+91L,         // 195!!!!
       1.927722175127463749675e+92L,         // 196!!!!
       5.644821477297368331929e+92L,         // 197!!!!
       3.068518756254966037203e+93L,         // 198!!!!
       1.180960052834948782460e+94L,         // 199!!!!
       3.855444350254927499350e+94L,         // 200!!!!
       1.134609116936771034718e+95L,         // 201!!!!
       6.198407887635031395150e+95L,         // 202!!!!
       2.397348907254946028394e+96L,         // 203!!!!
       7.865106474520052098675e+96L,         // 204!!!!
       2.325948689720380621171e+97L,         // 205!!!!
       1.276872024852816467401e+98L,         // 206!!!!
       4.962512238017738278775e+98L,         // 207!!!!
       1.635942146700170836524e+99L,         // 208!!!!
       4.861232761515595498248e+99L,         // 209!!!!
       2.681431252190914581542e+100L,        // 210!!!!
       1.047090082221742776822e+101L,        // 211!!!!
       3.468197351004362173432e+101L,        // 212!!!!
       1.035442578202821841127e+102L,        // 213!!!!
       5.738262879688557204499e+102L,        // 214!!!!
       2.251243676776746970166e+103L,        // 215!!!!
       7.491306278169422294612e+103L,        // 216!!!!
       2.246910394700123395245e+104L,        // 217!!!!
       1.250941307772105470581e+105L,        // 218!!!!
       4.930223652141075864664e+105L,        // 219!!!!
       1.648087381197272904815e+106L,        // 220!!!!
       4.965671972287272703492e+106L,        // 221!!!!
       2.777089703254074144689e+107L,        // 222!!!!
       1.099439874427459917820e+108L,        // 223!!!!
       3.691715733881891306785e+108L,        // 224!!!!
       1.117276193764636358286e+109L,        // 225!!!!
       6.276222729354207566998e+109L,        // 226!!!!
       2.495728514950334013452e+110L,        // 227!!!!
       8.417111873250712179470e+110L,        // 228!!!!
       2.558562483721017260474e+111L,        // 229!!!!
       1.443531227751467740410e+112L,        // 230!!!!
       5.765132869535271571073e+112L,        // 231!!!!
       1.952769954594165225637e+113L,        // 232!!!!
       5.961450587069970216905e+113L,        // 233!!!!
       3.377863072938434512558e+114L,        // 234!!!!
       1.354806224340788819202e+115L,        // 235!!!!
       4.608537092842229932503e+115L,        // 236!!!!
       1.412863789135582941406e+116L,        // 237!!!!
       8.039314113593474139889e+116L,        // 238!!!!
       3.237986876174485277893e+117L,        // 239!!!!
       1.106048902282135183801e+118L,        // 240!!!!
       3.405001731816754888790e+118L,        // 241!!!!
       1.945514015489620741853e+119L,        // 242!!!!
       7.868308109103999225281e+119L,        // 243!!!!
       2.698759321568409848474e+120L,        // 244!!!!
       8.342254242951049477535e+120L,        // 245!!!!
       4.785964478104467024959e+121L,        // 246!!!!
       1.943472102948687808644e+122L,        // 247!!!!
       6.692923117489656424215e+122L,        // 248!!!!
       2.077221306494811319906e+123L,        // 249!!!!
       1.196491119526116756240e+124L,        // 250!!!!
       4.878114978401206399697e+124L,        // 251!!!!
       1.686616625607393418902e+125L,        // 252!!!!
       5.255369905431872639362e+125L,        // 253!!!!
       3.039087443596336560849e+126L,        // 254!!!!
       1.243919319492307631923e+127L,        // 255!!!!
       4.317738561554927152390e+127L,        // 256!!!!
       1.350630065695991268316e+128L,        // 257!!!!
       7.840845604478548326990e+128L,        // 258!!!!
       3.221751037485076766680e+129L,        // 259!!!!
       1.122612026004281059621e+130L,        // 260!!!!
       3.525144471466537210305e+130L,        // 261!!!!
       2.054301548373379661671e+131L,        // 262!!!!
       8.473205228585751896369e+131L,        // 263!!!!
       2.963695748651301997400e+132L,        // 264!!!!
       9.341632849386323607309e+132L,        // 265!!!!
       5.464442118673189900046e+133L,        // 266!!!!
       2.262345796032395756330e+134L,        // 267!!!!
       7.942704606385489353033e+134L,        // 268!!!!
       2.512899236484921050366e+135L,        // 269!!!!
       1.475399372041761273012e+136L,        // 270!!!!
       6.130957107247792499655e+136L,        // 271!!!!
       2.160415652936853104025e+137L,        // 272!!!!
       6.860214915603834467499e+137L,        // 273!!!!
       4.042594279394425888054e+138L,        // 274!!!!
       1.686013204493142937405e+139L,        // 275!!!!
       5.962747202105714567109e+139L,        // 276!!!!
       1.900279531622262147497e+140L,        // 277!!!!
       1.123841209671650396879e+141L,        // 278!!!!
       4.703976840535868795361e+141L,        // 279!!!!
       1.669569216589600078790e+142L,        // 280!!!!
       5.339785483858556634467e+142L,        // 281!!!!
       3.169232211274054119199e+143L,        // 282!!!!
       1.331225445871650869087e+144L,        // 283!!!!
       4.741576575114464223765e+144L,        // 284!!!!
       1.521838862899688640823e+145L,        // 285!!!!
       9.064004124243794780908e+145L,        // 286!!!!
       3.820617029651637994280e+146L,        // 287!!!!
       1.365574053632965696444e+147L,        // 288!!!!
       4.398114313780100171979e+147L,        // 289!!!!
       2.628561196030700486463e+148L,        // 290!!!!
       1.111799555628626656335e+149L,        // 291!!!!
       3.987476236608259833617e+149L,        // 292!!!!
       1.288647493937569350390e+150L,        // 293!!!!
       7.727969916330259430202e+150L,        // 294!!!!
       3.279808689104448636190e+151L,        // 295!!!!
       1.180292966036044910751e+152L,        // 296!!!!
       3.827283056994580970658e+152L,        // 297!!!!
       2.302935035066417310200e+153L,        // 298!!!!
       9.806627980422301422207e+153L,        // 299!!!!
       3.540878898108134732252e+154L,        // 300!!!!
       1.152012200155368872168e+155L,        // 301!!!!
       6.954863805900580276805e+155L,        // 302!!!!
       2.971408278067957330929e+156L,        // 303!!!!
       1.076427185024872958605e+157L,        // 304!!!!
       3.513637210473875060113e+157L,        // 305!!!!
       2.128188324605577564702e+158L,        // 306!!!!
       9.122223413668629005951e+158L,        // 307!!!!
       3.315395729876608712502e+159L,        // 308!!!!
       1.085713898036427393575e+160L,        // 309!!!!
       6.597383806277290450577e+160L,        // 310!!!!
       2.837011481650943620851e+161L,        // 311!!!!
       1.034403467721501918301e+162L,        // 312!!!!
       3.398284500854017741889e+162L,        // 313!!!!
       2.071578515171069201481e+163L,        // 314!!!!
       8.936586167200472405680e+163L,        // 315!!!!
       3.268714957999946061830e+164L,        // 316!!!!
       1.077256186770723624179e+165L,        // 317!!!!
       6.587619678244000060710e+165L,        // 318!!!!
       2.850770987336950697412e+166L,        // 319!!!!
       1.045988786559982739786e+167L,        // 320!!!!
       3.457992359534022833614e+167L,        // 321!!!!
       2.121213536394568019549e+168L,        // 322!!!!
       9.207990289098350752640e+168L,        // 323!!!!
       3.389003668454344076906e+169L,        // 324!!!!
       1.123847516848557420925e+170L,        // 325!!!!
       6.915156128646291743729e+170L,        // 326!!!!
       3.011012824535160696113e+171L,        // 327!!!!
       1.111593203253024857225e+172L,        // 328!!!!
       3.697458330431753914842e+172L,        // 329!!!!
       2.282001522453276275431e+173L,        // 330!!!!
       9.966452449211381904135e+173L,        // 331!!!!
       3.690489434800042525987e+174L,        // 332!!!!
       1.231253624033774053642e+175L,        // 333!!!!
       7.621885084993942759938e+175L,        // 334!!!!
       3.338761570485812937885e+176L,        // 335!!!!
       1.240004450092814288732e+177L,        // 336!!!!
       4.149324712993818560775e+177L,        // 337!!!!
       2.576197158727952652859e+178L,        // 338!!!!
       1.131840172394690585943e+179L,        // 339!!!!
       4.216015130315568581688e+179L,        // 340!!!!
       1.414919727130892129224e+180L,        // 341!!!!
       8.810594282849598072778e+180L,        // 342!!!!
       3.882211791313788709785e+181L,        // 343!!!!
       1.450309204828555592101e+182L,        // 344!!!!
       4.881473058601577845823e+182L,        // 345!!!!
       3.048465621865960933181e+183L,        // 346!!!!
       1.347127491585884682295e+184L,        // 347!!!!
       5.047076032803373460510e+184L,        // 348!!!!
       1.703634097451950668192e+185L,        // 349!!!!
       1.066962967653086326613e+186L,        // 350!!!!
       4.728417495466455234857e+186L,        // 351!!!!
       1.776570763546787458100e+187L,        // 352!!!!
       6.013828364005385858719e+187L,        // 353!!!!
       3.777048905491925596211e+188L,        // 354!!!!
       1.678588210890591608374e+189L,        // 355!!!!
       6.324591918226563350835e+189L,        // 356!!!!
       2.146936725949922751563e+190L,        // 357!!!!
       1.352183508166109363444e+191L,        // 358!!!!
       6.026131677097223874063e+191L,        // 359!!!!
       2.276853090561562806300e+192L,        // 360!!!!
       7.750441580679221133141e+192L,        // 361!!!!
       4.894904299561315895666e+193L,        // 362!!!!
       2.187485798786292266285e+194L,        // 363!!!!
       8.287745249644088614934e+194L,        // 364!!!!
       2.828911176947915713597e+195L,        // 365!!!!
       1.791534973639441617814e+196L,        // 366!!!!
       8.028072881545692617266e+196L,        // 367!!!!
       3.049890251869024610296e+197L,        // 368!!!!
       1.043868224293780898317e+198L,        // 369!!!!
       6.628679402465933985911e+198L,        // 370!!!!
       2.978415039053451961006e+199L,        // 371!!!!
       1.134559173695277155030e+200L,        // 372!!!!
       3.893628476615802750723e+200L,        // 373!!!!
       2.479126096522259310731e+201L,        // 374!!!!
       1.116905639645044485377e+202L,        // 375!!!!
       4.265942493094242102913e+202L,        // 376!!!!
       1.467897935684157637023e+203L,        // 377!!!!
       9.371096644854140194562e+203L,        // 378!!!!
       4.233072374254718599579e+204L,        // 379!!!!
       1.621058147375811999107e+205L,        // 380!!!!
       5.592691134956640597056e+205L,        // 381!!!!
       3.579758918334281554323e+206L,        // 382!!!!
       1.621266719339557223639e+207L,        // 383!!!!
       6.224863285923118076570e+207L,        // 384!!!!
       2.153186086958306629866e+208L,        // 385!!!!
       1.381786942477032679969e+209L,        // 386!!!!
       6.274302203844086455482e+209L,        // 387!!!!
       2.415246954938169813709e+210L,        // 388!!!!
       8.375893878267812790181e+210L,        // 389!!!!
       5.388969075660427451878e+211L,        // 390!!!!
       2.453252161703037804094e+212L,        // 391!!!!
       9.467768063357625669740e+212L,        // 392!!!!
       3.291726294159250426541e+213L,        // 393!!!!
       2.123253815810208416040e+214L,        // 394!!!!
       9.690346038726999326169e+214L,        // 395!!!!
       3.749236153089619765217e+215L,        // 396!!!!
       1.306815338781222419337e+216L,        // 397!!!!
       8.450550186924629495838e+216L,        // 398!!!!
       3.866448069452072731142e+217L,        // 399!!!!
       1.499694461235847906087e+218L,        // 400!!!!
       5.240329508512701901540e+218L,        // 401!!!!
       3.397121175143701057327e+219L,        // 402!!!!
       1.558178571989185310650e+220L,        // 403!!!!
       6.058765623392825540591e+220L,        // 404!!!!
       2.122333450947644270124e+221L,        // 405!!!!
       1.379231197108342629275e+222L,        // 406!!!!
       6.341786787995984214346e+222L,        // 407!!!!
       2.471976374344272820561e+223L,        // 408!!!!
       8.680343814375865064807e+223L,        // 409!!!!
       5.654847908144204780026e+224L,        // 410!!!!
       2.606474369866349512096e+225L,        // 411!!!!
       1.018454266229840402071e+226L,        // 412!!!!
       3.584981995337232271765e+226L,        // 413!!!!
       2.341107033971700778931e+227L,        // 414!!!!
       1.081686863494535047520e+228L,        // 415!!!!
       4.236769747516136072616e+228L,        // 416!!!!
       1.494937492055625857326e+229L,        // 417!!!!
       9.785827402001709255931e+229L,        // 418!!!!
       4.532267958042101849108e+230L,        // 419!!!!
       1.779443293956777150499e+231L,        // 420!!!!
       6.293686841554184859343e+231L,        // 421!!!!
       4.129619163644721306003e+232L,        // 422!!!!
       1.917149346251809082173e+233L,        // 423!!!!
       7.544839566376735118115e+233L,        // 424!!!!
       2.674816907660528565221e+234L,        // 425!!!!
       1.759217763712651276357e+235L,        // 426!!!!
       8.186227708495224780878e+235L,        // 427!!!!
       3.229191334409242630553e+236L,        // 428!!!!
       1.147496453386366754480e+237L,        // 429!!!!
       7.564636383964400488336e+237L,        // 430!!!!
       3.528264142361441880558e+238L,        // 431!!!!
       1.395010656464792816399e+239L,        // 432!!!!
       4.968659643162968046897e+239L,        // 433!!!!
       3.283052190640549811938e+240L,        // 434!!!!
       1.534794901927227218043e+241L,        // 435!!!!
       6.082246462186496679499e+241L,        // 436!!!!
       2.171304264062217036494e+242L,        // 437!!!!
       1.437976859500560817629e+243L,        // 438!!!!
       6.737749619460527487208e+243L,        // 439!!!!
       2.676188443362058538980e+244L,        // 440!!!!
       9.575451804514377130938e+244L,        // 441!!!!
       6.355857718992478813919e+245L,        // 442!!!!
       2.984823081421013676833e+246L,        // 443!!!!
       1.188227668852753991307e+247L,        // 444!!!!
       4.261076053008897823268e+247L,        // 445!!!!
       2.834712542670645551008e+248L,        // 446!!!!
       1.334215917395193113544e+249L,        // 447!!!!
       5.323259956460337881055e+249L,        // 448!!!!
       1.913223147800995122647e+250L,        // 449!!!!
       1.275620644201790497954e+251L,        // 450!!!!
       6.017313787452320942086e+251L,        // 451!!!!
       2.406113500320072722237e+252L,        // 452!!!!
       8.666900859538507905592e+252L,        // 453!!!!
       5.791317724676128860709e+253L,        // 454!!!!
       2.737877773290806028649e+254L,        // 455!!!!
       1.097187756145953161340e+255L,        // 456!!!!
       3.960773692809098112855e+255L,        // 457!!!!
       2.652423517901667018205e+256L,        // 458!!!!
       1.256685897940479967150e+257L,        // 459!!!!
       5.047063678271384542164e+257L,        // 460!!!!
       1.825916672384994230026e+258L,        // 461!!!!
       1.225419665270570162411e+259L,        // 462!!!!
       5.818455707464422247904e+259L,        // 463!!!!
       2.341837546717922427564e+260L,        // 464!!!!
       8.490512526590223169622e+260L,        // 465!!!!
       5.710455640160856956834e+261L,        // 466!!!!
       2.717218815385885189771e+262L,        // 467!!!!
       1.095979971863987696100e+263L,        // 468!!!!
       3.982050374970814666553e+263L,        // 469!!!!
       2.683914150875602769712e+264L,        // 470!!!!
       1.279810062046751924382e+265L,        // 471!!!!
       5.173025467198021925592e+265L,        // 472!!!!
       1.883509827361195337280e+266L,        // 473!!!!
       1.272175307515035712843e+267L,        // 474!!!!
       6.079097794722071640815e+267L,        // 475!!!!
       2.462360122386258436582e+268L,        // 476!!!!
       8.984341876512901758823e+268L,        // 477!!!!
       6.080997969921870707392e+269L,        // 478!!!!
       2.911887843671872315951e+270L,        // 479!!!!
       1.181932858745404049559e+271L,        // 480!!!!
       4.321468442602705745994e+271L,        // 481!!!!
       2.931041021502341680963e+272L,        // 482!!!!
       1.406441828493514328604e+273L,        // 483!!!!
       5.720555036327755599867e+273L,        // 484!!!!
       2.095912194662312286807e+274L,        // 485!!!!
       1.424485936450138056948e+275L,        // 486!!!!
       6.849371704763414780302e+275L,        // 487!!!!
       2.791630857727944732735e+276L,        // 488!!!!
       1.024901063189870708249e+277L,        // 489!!!!
       6.979981088605676479045e+277L,        // 490!!!!
       3.363041507038836657128e+278L,        // 491!!!!
       1.373482382002148808506e+279L,        // 492!!!!
       5.052762241526062591666e+279L,        // 493!!!!
       3.448110657771204180648e+280L,        // 494!!!!
       1.664705545984224145279e+281L,        // 495!!!!
       6.812472614730658090188e+281L,        // 496!!!!
       2.511222834038453108058e+282L,        // 497!!!!
       1.717159107570059681963e+283L,        // 498!!!!
       8.306880674461278484940e+283L,        // 499!!!!
       3.406236307365329045094e+284L,        // 500!!!!
       1.258122639853265007137e+285L,        // 501!!!!
       8.620138720001699603453e+285L,        // 502!!!!
       4.178360979254023077925e+286L,        // 503!!!!
       1.716743098912125838727e+287L,        // 504!!!!
       6.353519331258988286042e+287L,        // 505!!!!
       4.361790192320859999347e+288L,        // 506!!!!
       2.118429016481789700508e+289L,        // 507!!!!
       8.721054942473599260735e+289L,        // 508!!!!
       3.233941339610825037595e+290L,        // 509!!!!
       2.224512998083638599667e+291L,        // 510!!!!
       1.082517227422194536960e+292L,        // 511!!!!
       4.465180130546482821496e+292L,        // 512!!!!
       1.659011907220353244286e+293L,        // 513!!!!
       1.143399681014990240229e+294L,        // 514!!!!
       5.574963721224301865342e+294L,        // 515!!!!
       2.304032947361985135892e+295L,        // 516!!!!
       8.577091560329226272961e+295L,        // 517!!!!
       5.922810347657649444386e+296L,        // 518!!!!
       2.893406171315412668112e+297L,        // 519!!!!
       1.198097132628232270664e+298L,        // 520!!!!
       4.468664702931526888213e+298L,        // 521!!!!
       3.091707001477293009969e+299L,        // 522!!!!
       1.513251427597960825423e+300L,        // 523!!!!
       6.278028974971937098279e+300L,        // 524!!!!
       2.346048969039051616312e+301L,        // 525!!!!
       1.626237882777056123244e+302L,        // 526!!!!
       7.974835023441253549978e+302L,        // 527!!!!
       3.314799298785182787891e+303L,        // 528!!!!
       1.241059904621658305029e+304L,        // 529!!!!
       8.619060778718397453192e+304L,        // 530!!!!
       4.234637397447305635038e+305L,        // 531!!!!
       1.763473226953717243158e+306L,        // 532!!!!
       6.614849291633438765804e+306L,        // 533!!!!
       4.602578455835624240005e+307L,        // 534!!!!
                                        };

static const int N = sizeof(quadruple_factorials) / sizeof(long double);

////////////////////////////////////////////////////////////////////////////////
// double Quadruple_Factorial( int n )                                        //
//                                                                            //
//  Description:                                                              //
//     This function computes n!!!! for -3<=n<=Quadruple_Factorial_Max_Arg(). //
//     If n > Quadruple_Factorial_Max_Arg(), then DBL_MAX is returned and if  //
//     n < -3, then 0 is returned.                                            //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Quadruple Factorial function.               //
//                                                                            //
//  Return Values:                                                            //
//     If n < -3, then 0 is returned.  If n > Quadruple_Factorial_Max_Arg(),  //
//     then DBL_MAX is returned.  If -3 <= n <= Quadruple_Factorial_Max_Arg(),//
//     then n!!!! is returned.                                                //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     double x;                                                              //
//                                                                            //
//     x = Quadruple_Factorial( n );                                          //
////////////////////////////////////////////////////////////////////////////////
double Quadruple_Factorial(int n) {

               // For an argument n < -3, return 0.0 //

   if ( n < -3 ) return 0.0;


           // For a large postive argument (n >= N) return DBL_MAX //

   if ( n >= N ) return DBL_MAX;

                         // Otherwise return n!!!!. //

   if (n <= 0) return 1.0;
   return (double) quadruple_factorials[n];
}


////////////////////////////////////////////////////////////////////////////////
// long double xQuadruple_Factorial( int n )                                  //
//                                                                            //
//  Description:                                                              //
//     This function computes n!!!! for -3<=n<=Quadruple_Factorial_Max_Arg(). //
//     If n > Quadruple_Factorial_Max_Arg(), then DBL_MAX is returned and if  //
//     n < -3, then 0 is returned.                                            //
//                                                                            //
//  Arguments:                                                                //
//     int    n   Argument of the Quadruple Factorial function.               //
//                                                                            //
//  Return Values:                                                            //
//     If n < -3, then 0 is returned.  If n > Quadruple_Factorial_Max_Arg(),  //
//     then DBL_MAX is returned.  If -3 <= n <= Quadruple_Factorial_Max_Arg(),//
//     then n!!!! is returned.                                                //
//                                                                            //
//  Example:                                                                  //
//     int n;                                                                 //
//     long double x;                                                         //
//                                                                            //
//     x = xQuadruple_Factorial( n );                                         //
////////////////////////////////////////////////////////////////////////////////
long double xQuadruple_Factorial(int n) {

               // For an argument n < -3, return 0.0 //

   if ( n < -3 ) return 0.0L;


           // For a large postive argument (n >= N) return DBL_MAX //

   if ( n >= N ) return DBL_MAX;

                         // Otherwise return n!!!!. //

   if (n <= 0) return 1.0L;
   return quadruple_factorials[n];
}


////////////////////////////////////////////////////////////////////////////////
// int Quadruple_Factorial_Max_Arg( void )                                    //
//                                                                            //
//  Description:                                                              //
//     This function returns the maximum argument of the Quadruple Factorial  //
//     function for which a positive number < DBL_MAX is returned.            //
//                                                                            //
//  Arguments:                                                                //
//     none                                                                   //
//                                                                            //
//  Return Values:                                                            //
//     N-1                                                                    //
//                                                                            //
//  Example:                                                                  //
//     int x;                                                                 //
//                                                                            //
//     x = Quadruple_Factorial_Max_Arg();                                     //
////////////////////////////////////////////////////////////////////////////////

int Quadruple_Factorial_Max_Arg( void ) { return N - 1; }
