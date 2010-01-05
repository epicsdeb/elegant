/* GSL
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Jorma Olavi Tähtinen, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 */


#include "gsl_code.h"

gsl_error_handler_t * gsl_error_handler = NULL;
gsl_stream_handler_t * gsl_stream_handler = NULL;
FILE * gsl_stream = NULL ;



static double am21_data[37] = {
  0.0065809191761485,
  0.0023675984685722,
  0.0001324741670371,
  0.0000157600904043,
  0.0000027529702663,
  0.0000006102679017,
  0.0000001595088468,
  0.0000000471033947,
  0.0000000152933871,
  0.0000000053590722,
  0.0000000020000910,
  0.0000000007872292,
  0.0000000003243103,
  0.0000000001390106,
  0.0000000000617011,
  0.0000000000282491,
  0.0000000000132979,
  0.0000000000064188,
  0.0000000000031697,
  0.0000000000015981,
  0.0000000000008213,
  0.0000000000004296,
  0.0000000000002284,
  0.0000000000001232,
  0.0000000000000675,
  0.0000000000000374,
  0.0000000000000210,
  0.0000000000000119,
  0.0000000000000068,
  0.0000000000000039,
  0.0000000000000023,
  0.0000000000000013,
  0.0000000000000008,
  0.0000000000000005,
  0.0000000000000003,
  0.0000000000000001,
  0.0000000000000001
};
static cheb_series am21_cs = {
  am21_data,
  36,
  -1, 1,
  20
};
static double ath1_data[36] = {
  -0.07125837815669365,
  -0.00590471979831451,
  -0.00012114544069499,
  -0.00000988608542270,
  -0.00000138084097352,
  -0.00000026142640172,
  -0.00000006050432589,
  -0.00000001618436223,
  -0.00000000483464911,
  -0.00000000157655272,
  -0.00000000055231518,
  -0.00000000020545441,
  -0.00000000008043412,
  -0.00000000003291252,
  -0.00000000001399875,
  -0.00000000000616151,
  -0.00000000000279614,
  -0.00000000000130428,
  -0.00000000000062373,
  -0.00000000000030512,
  -0.00000000000015239,
  -0.00000000000007758,
  -0.00000000000004020,
  -0.00000000000002117,
  -0.00000000000001132,
  -0.00000000000000614,
  -0.00000000000000337,
  -0.00000000000000188,
  -0.00000000000000105,
  -0.00000000000000060,
  -0.00000000000000034,
  -0.00000000000000020,
  -0.00000000000000011,
  -0.00000000000000007,
  -0.00000000000000004,
  -0.00000000000000002
};
static cheb_series ath1_cs = {
  ath1_data,
  35,
  -1, 1,
  15
};
static double am22_data[33] = {
 -0.01562844480625341,
  0.00778336445239681,
  0.00086705777047718,
  0.00015696627315611,
  0.00003563962571432,
  0.00000924598335425,
  0.00000262110161850,
  0.00000079188221651,
  0.00000025104152792,
  0.00000008265223206,
  0.00000002805711662,
  0.00000000976821090,
  0.00000000347407923,
  0.00000000125828132,
  0.00000000046298826,
  0.00000000017272825,
  0.00000000006523192,
  0.00000000002490471,
  0.00000000000960156,
  0.00000000000373448,
  0.00000000000146417,
  0.00000000000057826,
  0.00000000000022991,
  0.00000000000009197,
  0.00000000000003700,
  0.00000000000001496,
  0.00000000000000608,
  0.00000000000000248,
  0.00000000000000101,
  0.00000000000000041,
  0.00000000000000017,
  0.00000000000000007,
  0.00000000000000002
};
static cheb_series am22_cs = {
  am22_data,
  32,
  -1, 1,
  15
};

static double ath2_data[32] = {
   0.00440527345871877,
  -0.03042919452318455,
  -0.00138565328377179,
  -0.00018044439089549,
  -0.00003380847108327,
  -0.00000767818353522,
  -0.00000196783944371,
  -0.00000054837271158,
  -0.00000016254615505,
  -0.00000005053049981,
  -0.00000001631580701,
  -0.00000000543420411,
  -0.00000000185739855,
  -0.00000000064895120,
  -0.00000000023105948,
  -0.00000000008363282,
  -0.00000000003071196,
  -0.00000000001142367,
  -0.00000000000429811,
  -0.00000000000163389,
  -0.00000000000062693,
  -0.00000000000024260,
  -0.00000000000009461,
  -0.00000000000003716,
  -0.00000000000001469,
  -0.00000000000000584,
  -0.00000000000000233,
  -0.00000000000000093,
  -0.00000000000000037,
  -0.00000000000000015,
  -0.00000000000000006,
  -0.00000000000000002
};
static cheb_series ath2_cs = {
  ath2_data,
  31,
  -1, 1,
  16
};
static double ai_data_f[9] = {
  -0.03797135849666999750,
   0.05919188853726363857,
   0.00098629280577279975,
   0.00000684884381907656,
   0.00000002594202596219,
   0.00000000006176612774,
   0.00000000000010092454,
   0.00000000000000012014,
   0.00000000000000000010
};
static cheb_series aif_cs = {
  ai_data_f,
  8,
  -1, 1,
  8
};
static double ai_data_g[8] = {
   0.01815236558116127,
   0.02157256316601076,
   0.00025678356987483,
   0.00000142652141197,
   0.00000000457211492,
   0.00000000000952517,
   0.00000000000001392,
   0.00000000000000001
};
static cheb_series aig_cs = {
  ai_data_g,
  7,
  -1, 1,
  7
};
static double data_aip[36] = {
 -0.0187519297793867540198,
 -0.0091443848250055004725,
  0.0009010457337825074652,
 -0.0001394184127221491507,
  0.0000273815815785209370,
 -0.0000062750421119959424,
  0.0000016064844184831521,
 -0.0000004476392158510354,
  0.0000001334635874651668,
 -0.0000000420735334263215,
  0.0000000139021990246364,
 -0.0000000047831848068048,
  0.0000000017047897907465,
 -0.0000000006268389576018,
  0.0000000002369824276612,
 -0.0000000000918641139267,
  0.0000000000364278543037,
 -0.0000000000147475551725,
  0.0000000000060851006556,
 -0.0000000000025552772234,
  0.0000000000010906187250,
 -0.0000000000004725870319,
  0.0000000000002076969064,
 -0.0000000000000924976214,
  0.0000000000000417096723,
 -0.0000000000000190299093,
  0.0000000000000087790676,
 -0.0000000000000040927557,
  0.0000000000000019271068,
 -0.0000000000000009160199,
  0.0000000000000004393567,
 -0.0000000000000002125503,
  0.0000000000000001036735,
 -0.0000000000000000509642,
  0.0000000000000000252377,
 -0.0000000000000000125793
};
static cheb_series aip_cs = {
  data_aip,
  35,
  -1, 1,
  17
};
static double cos_data[11] = {
  0.165391825637921473505668118136,
 -0.00084852883845000173671196530195,
 -0.000210086507222940730213625768083,
  1.16582269619760204299639757584e-6,
  1.43319375856259870334412701165e-7,
 -7.4770883429007141617951330184e-10,
 -6.0969994944584252706997438007e-11,
  2.90748249201909353949854872638e-13,
  1.77126739876261435667156490461e-14,
 -7.6896421502815579078577263149e-17,
 -3.7363121133079412079201377318e-18
};
static cheb_series cos_cs = {
  cos_data,
  10,
  -1, 1,
  10
};
static double sin_data[12] = {
  -0.3295190160663511504173,
   0.0025374284671667991990,
   0.0006261928782647355874,
  -4.6495547521854042157541e-06,
  -5.6917531549379706526677e-07,
   3.7283335140973803627866e-09,
   3.0267376484747473727186e-10,
  -1.7400875016436622322022e-12,
  -1.0554678305790849834462e-13,
   5.3701981409132410797062e-16,
   2.5984137983099020336115e-17,
  -1.1821555255364833468288e-19
};
static cheb_series sin_cs = {
  sin_data,
  11,
  -1, 1,
  11
};
static double data_bif[9] = {
  -0.01673021647198664948,
   0.10252335834249445610,
   0.00170830925073815165,
   0.00001186254546774468,
   0.00000004493290701779,
   0.00000000010698207143,
   0.00000000000017480643,
   0.00000000000000020810,
   0.00000000000000000018
};
static cheb_series bif_cs = {
  data_bif,
  8,
  -1, 1,
  8
};
static double data_big[8] = {
   0.02246622324857452,
   0.03736477545301955,
   0.00044476218957212,
   0.00000247080756363,
   0.00000000791913533,
   0.00000000001649807,
   0.00000000000002411,
   0.00000000000000002
};
static cheb_series big_cs = {
  data_big,
  7,
  -1, 1,
  7
};
static double data_bif2[10] = {
  0.0998457269381604100,
  0.4786249778630055380,
  0.0251552119604330118,
  0.0005820693885232645,
  0.0000074997659644377,
  0.0000000613460287034,
  0.0000000003462753885,
  0.0000000000014288910,
  0.0000000000000044962,
  0.0000000000000000111
};
static cheb_series bif2_cs = {
  data_bif2,
  9,
  -1, 1,
  9
};
static double data_big2[10] = {
  0.033305662145514340,
  0.161309215123197068,
  0.0063190073096134286,
  0.0001187904568162517,
  0.0000013045345886200,
  0.0000000093741259955,
  0.0000000000474580188,
  0.0000000000001783107,
  0.0000000000000005167,
  0.0000000000000000011
};
static cheb_series big2_cs = {
  data_big2,
  9,
  -1, 1,
  9
};
static double data_bip[24] = {
  -0.08322047477943447,
   0.01146118927371174,
   0.00042896440718911,
  -0.00014906639379950,
  -0.00001307659726787,
   0.00000632759839610,
  -0.00000042226696982,
  -0.00000019147186298,
   0.00000006453106284,
  -0.00000000784485467,
  -0.00000000096077216,
   0.00000000070004713,
  -0.00000000017731789,
   0.00000000002272089,
   0.00000000000165404,
  -0.00000000000185171,
   0.00000000000059576,
  -0.00000000000012194,
   0.00000000000001334,
   0.00000000000000172,
  -0.00000000000000145,
   0.00000000000000049,
  -0.00000000000000011,
   0.00000000000000001
};
static cheb_series bip_cs = {
  data_bip,
  23,
  -1, 1,
  14
};
static double data_bip2[29] = {    
  -0.113596737585988679,
   0.0041381473947881595,
   0.0001353470622119332,
   0.0000104273166530153,
   0.0000013474954767849,
   0.0000001696537405438,
  -0.0000000100965008656,
  -0.0000000167291194937,
  -0.0000000045815364485,
   0.0000000003736681366,
   0.0000000005766930320,
   0.0000000000621812650,
  -0.0000000000632941202,
  -0.0000000000149150479,
   0.0000000000078896213,
   0.0000000000024960513,
  -0.0000000000012130075,
  -0.0000000000003740493,
   0.0000000000002237727,
   0.0000000000000474902,
  -0.0000000000000452616,
  -0.0000000000000030172,
   0.0000000000000091058,
  -0.0000000000000009814,
  -0.0000000000000016429,
   0.0000000000000005533,
   0.0000000000000002175,
  -0.0000000000000001737,
  -0.0000000000000000010
};
static cheb_series bip2_cs = {
  data_bip2,
  28,
  -1, 1,
  10
};
/*************** DERIV **************/
static double aif_data[8] = {
   0.10527461226531408809,
   0.01183613628152997844,
   0.00012328104173225664,
   0.00000062261225638140,
   0.00000000185298887844,
   0.00000000000363328873,
   0.00000000000000504622,
   0.00000000000000000522
};
static cheb_series aif_cs_der = {
  aif_data,
  7,
  -1, 1,
  7
};
static double aig_data[9] = {
   0.021233878150918666852,
   0.086315930335214406752,
   0.001797594720383231358,
   0.000014265499875550693,
   0.000000059437995283683,
   0.000000000152403366479,
   0.000000000000264587660,
   0.000000000000000331562,
   0.000000000000000000314
};
static cheb_series aig_cs_der = {
  aig_data,
  8,
  -1, 1,
  8
};
static double an20_data[16] = {
    0.0126732217145738027,
   -0.0005212847072615621,
   -0.0000052672111140370,
   -0.0000001628202185026,
   -0.0000000090991442687,
   -0.0000000007438647126,
   -0.0000000000795494752,
   -0.0000000000104050944,
   -0.0000000000015932426,
   -0.0000000000002770648,
   -0.0000000000000535343,
   -0.0000000000000113062,
   -0.0000000000000025772,
   -0.0000000000000006278,
   -0.0000000000000001621,
   -0.0000000000000000441
};
static cheb_series an20_cs_der = {
  an20_data,
  15,
  -1, 1,
  8
};
static double aph0_data[15] = {
 -0.0855849241130933257,
  0.0011214378867065261,
  0.0000042721029353664,
  0.0000000817607381483,
  0.0000000033907645000,
  0.0000000002253264423,
  0.0000000000206284209,
  0.0000000000023858763,
  0.0000000000003301618,
  0.0000000000000527010,
  0.0000000000000094555,
  0.0000000000000018709,
  0.0000000000000004024,
  0.0000000000000000930,
  0.0000000000000000229
};
static cheb_series aph0_cs_der = {
  aph0_data,
  14,
  -1, 1,
  7
};
static double an21_data[24] = {
    0.0198313155263169394,
   -0.0029376249067087533,
   -0.0001136260695958196,
   -0.0000100554451087156,
   -0.0000013048787116563,
   -0.0000002123881993151,
   -0.0000000402270833384,
   -0.0000000084996745953,
   -0.0000000019514839426,
   -0.0000000004783865344,
   -0.0000000001236733992,
   -0.0000000000334137486,
   -0.0000000000093702824,
   -0.0000000000027130128,
   -0.0000000000008075954,
   -0.0000000000002463214,
   -0.0000000000000767656,
   -0.0000000000000243883,
   -0.0000000000000078831,
   -0.0000000000000025882,
   -0.0000000000000008619,
   -0.0000000000000002908,
   -0.0000000000000000993,
   -0.0000000000000000343
};
static cheb_series an21_cs_der = {
  an21_data,
  23,
  -1, 1,
  12
};
static double aph1_data[22] = {
  -0.1024172908077571694,
   0.0071697275146591248,
   0.0001209959363122329,
   0.0000073361512841220,
   0.0000007535382954272,
   0.0000001041478171741,
   0.0000000174358728519,
   0.0000000033399795033,
   0.0000000007073075174,
   0.0000000001619187515,
   0.0000000000394539982,
   0.0000000000101192282,
   0.0000000000027092778,
   0.0000000000007523806,
   0.0000000000002156369,
   0.0000000000000635283,
   0.0000000000000191757,
   0.0000000000000059143,
   0.0000000000000018597,
   0.0000000000000005950,
   0.0000000000000001934,
   0.0000000000000000638
};
static cheb_series aph1_cs_der = {
  aph1_data,
  21,
  -1, 1,
  10
};
static double an22_data[33] = {
    0.0537418629629794329,
   -0.0126661435859883193,
   -0.0011924334106593007,
   -0.0002032327627275655,
   -0.0000446468963075164,
   -0.0000113359036053123,
   -0.0000031641352378546,
   -0.0000009446708886149,
   -0.0000002966562236472,
   -0.0000000969118892024,
   -0.0000000326822538653,
   -0.0000000113144618964,
   -0.0000000040042691002,
   -0.0000000014440333684,
   -0.0000000005292853746,
   -0.0000000001967763374,
   -0.0000000000740800096,
   -0.0000000000282016314,
   -0.0000000000108440066,
   -0.0000000000042074801,
   -0.0000000000016459150,
   -0.0000000000006486827,
   -0.0000000000002574095,
   -0.0000000000001027889,
   -0.0000000000000412846,
   -0.0000000000000166711,
   -0.0000000000000067657,
   -0.0000000000000027585,
   -0.0000000000000011296,
   -0.0000000000000004645,
   -0.0000000000000001917,
   -0.0000000000000000794,
   -0.0000000000000000330
};
static cheb_series an22_cs_der = {
  an22_data,
  32,
  -1, 1,
  18
};
static double aph2_data[32] = {
   -0.2057088719781465107,
    0.0422196961357771922,
    0.0020482560511207275,
    0.0002607800735165006,
    0.0000474824268004729,
    0.0000105102756431612,
    0.0000026353534014668,
    0.0000007208824863499,
    0.0000002103236664473,
    0.0000000644975634555,
    0.0000000205802377264,
    0.0000000067836273921,
    0.0000000022974015284,
    0.0000000007961306765,
    0.0000000002813860610,
    0.0000000001011749057,
    0.0000000000369306738,
    0.0000000000136615066,
    0.0000000000051142751,
    0.0000000000019351689,
    0.0000000000007393607,
    0.0000000000002849792,
    0.0000000000001107281,
    0.0000000000000433412,
    0.0000000000000170801,
    0.0000000000000067733,
    0.0000000000000027017,
    0.0000000000000010835,
    0.0000000000000004367,
    0.0000000000000001769,
    0.0000000000000000719,
    0.0000000000000000294
};
static cheb_series aph2_cs_der = {
  aph2_data,
  31,
  -1, 1,
  16
};
static double aip1_data[25] = {
    0.0358865097808301538,
    0.0114668575627764899,
   -0.0007592073583861400,
    0.0000869517610893841,
   -0.0000128237294298592,
    0.0000022062695681038,
   -0.0000004222295185921,
    0.0000000874686415726,
   -0.0000000192773588418,
    0.0000000044668460054,
   -0.0000000010790108052,
    0.0000000002700029447,
   -0.0000000000696480108,
    0.0000000000184489907,
   -0.0000000000050027817,
    0.0000000000013852243,
   -0.0000000000003908218,
    0.0000000000001121536,
   -0.0000000000000326862,
    0.0000000000000096619,
   -0.0000000000000028935,
    0.0000000000000008770,
   -0.0000000000000002688,
    0.0000000000000000832,
   -0.0000000000000000260
};
static cheb_series aip1_cs_der = {
  aip1_data,
  24,
  -1, 1,
  14
};
static double aip2_data[15] = {
    0.0065457691989713757,
    0.0023833724120774592,
   -0.0000430700770220586,
    0.0000015629125858629,
   -0.0000000815417186163,
    0.0000000054103738057,
   -0.0000000004284130883,
    0.0000000000389497963,
   -0.0000000000039623161,
    0.0000000000004428184,
   -0.0000000000000536297,
    0.0000000000000069650,
   -0.0000000000000009620,
    0.0000000000000001403,
   -0.0000000000000000215
};
static cheb_series aip2_cs_der = {
  aip2_data,
  14,
  -1, 1,
  9
};
static double bif_data[8] = {
   0.1153536790828570243,
   0.0205007894049192875,
   0.0002135290278902876,
   0.0000010783960614677,
   0.0000000032094708833,
   0.0000000000062930407,
   0.0000000000000087403,
   0.0000000000000000090
};
static cheb_series bif_cs_der = {
  bif_data,
  7,
  -1, 1,
  7
};
static double big_data[9] = {
   -0.097196440416443537390,
    0.149503576843167066571,
    0.003113525387121326042,
    0.000024708570579821297,
    0.000000102949627731379,
    0.000000000263970373987,
    0.000000000000458279271,
    0.000000000000000574283,
    0.000000000000000000544
};
static cheb_series big_cs_der = {
  big_data,
  8,
  -1, 1,
  8
};
static double bif2_data[10] = {
   0.323493987603522033521,
   0.086297871535563559139,
   0.002994025552655397426,
   0.000051430528364661637,
   0.000000525840250036811,
   0.000000003561751373958,
   0.000000000017146864007,
   0.000000000000061663520,
   0.000000000000000171911,
   0.000000000000000000382
};
static cheb_series bif2_cs_der = {
  bif2_data,
  9,
  -1, 1,
  9
};
static double big2_data[10] = {
   1.6062999463621294578,
   0.7449088819876088652,
   0.0470138738610277380,
   0.0012284422062548239,
   0.0000173222412256624,
   0.0000001521901652368,
   0.0000000009113560249,
   0.0000000000039547918,
   0.0000000000000130017,
   0.0000000000000000335
};
static cheb_series big2_cs_der = {
  big2_data,
  9,
  -1, 1,
  9
};
static double bip1_data[24] = {
   -0.1729187351079553719,
   -0.0149358492984694364,
   -0.0005471104951678566,
    0.0001537966292958408,
    0.0000154353476192179,
   -0.0000065434113851906,
    0.0000003728082407879,
    0.0000002072078388189,
   -0.0000000658173336470,
    0.0000000074926746354,
    0.0000000011101336884,
   -0.0000000007265140553,
    0.0000000001782723560,
   -0.0000000000217346352,
   -0.0000000000020302035,
    0.0000000000019311827,
   -0.0000000000006044953,
    0.0000000000001209450,
   -0.0000000000000125109,
   -0.0000000000000019917,
    0.0000000000000015154,
   -0.0000000000000004977,
    0.0000000000000001155,
   -0.0000000000000000186
};
static cheb_series bip1_cs_der = {
  bip1_data,
  23,
  -1, 1,
  13
};
static double bip2_data[29] = {
    -0.13269705443526630495,
    -0.00568443626045977481,
    -0.00015643601119611610,
    -0.00001136737203679562,
    -0.00000143464350991284,
    -0.00000018098531185164,
     0.00000000926177343611,
     0.00000001710005490721,
     0.00000000476698163504,
    -0.00000000035195022023,
    -0.00000000058890614316,
    -0.00000000006678499608,
     0.00000000006395565102,
     0.00000000001554529427,
    -0.00000000000792397000,
    -0.00000000000258326243,
     0.00000000000121655048,
     0.00000000000038707207,
    -0.00000000000022487045,
    -0.00000000000004953477,
     0.00000000000004563782,
     0.00000000000000332998,
    -0.00000000000000921750,
     0.00000000000000094157,
     0.00000000000000167154,
    -0.00000000000000055134,
    -0.00000000000000022369,
     0.00000000000000017487,
     0.00000000000000000207
};
static cheb_series bip2_cs_der = {
  bip2_data,
  28,
  -1, 1,
  14
};




gsl_complex gsl_complex_mul_real (gsl_complex a, double x) {
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x * GSL_REAL (a), x * GSL_IMAG (a));
  return z;
}
gsl_complex gsl_complex_rect (double x, double y) {
  gsl_complex z;
  GSL_SET_COMPLEX (&z, x, y);
  return z;
}
double gsl_sf_airy_Ai(const double x, gsl_mode_t mode) {
  EVAL_RESULT(gsl_sf_airy_Ai_e(x, mode, &result));
}
double gsl_sf_airy_Bi(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_e(x, mode, &result));
}

unsigned int GSL_MODE_PREC(gsl_mode_t mt);
unsigned int GSL_MODE_PREC(gsl_mode_t mt){ 
  return  (mt & (unsigned int)7); 
}
static int cheb_eval_mode_e(const cheb_series * cs,
                 const double x,
                 gsl_mode_t mode,
                 gsl_sf_result * result) {
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  int eval_order;

  if(GSL_MODE_PREC(mode) == GSL_PREC_DOUBLE)
    eval_order = cs->order;
  else
    eval_order = cs->order_sp;

  for(j = eval_order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    dd = temp;
  }

  result->val = y*d - dd + 0.5 * cs->c[0];
  result->err = GSL_DBL_EPSILON * fabs(result->val) + fabs(cs->c[eval_order]);
  return GSL_SUCCESS;
}
static int airy_mod_phase(const double x, gsl_mode_t mode, gsl_sf_result * mod, gsl_sf_result * phase) {
  gsl_sf_result result_m;
  gsl_sf_result result_p;
  double m, p;
  double sqx;

  if(x < -2.0) {
    double z = 16.0/(x*x*x) + 1.0;
    cheb_eval_mode_e(&am21_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath1_cs, z, mode, &result_p);
  }
  else if(x <= -1.0) {
    double z = (16.0/(x*x*x) + 9.0)/7.0;
    cheb_eval_mode_e(&am22_cs, z, mode, &result_m);
    cheb_eval_mode_e(&ath2_cs, z, mode, &result_p);
  }
  else {
    mod->val = 0.0;
    mod->err = 0.0;
    phase->val = 0.0;
    phase->err = 0.0;
    GSL_ERROR ("x is greater than 1.0", GSL_EDOM);
  }

  m =  0.3125 + result_m.val;
  p = -0.625  + result_p.val;

  sqx = sqrt(-x);

  mod->val   = sqrt(m/sqx);
  mod->err  = fabs(mod->val) * (GSL_DBL_EPSILON + fabs(result_m.err/result_m.val));
  phase->val = M_PI_4 - x*sqx * p;
  phase->err = fabs(phase->val) * (GSL_DBL_EPSILON + fabs(result_p.err/result_p.val));

  return GSL_SUCCESS;
}
static int airy_aie(const double x, gsl_mode_t mode, gsl_sf_result * result) {
  double sqx = sqrt(x);
  double z   = 2.0/(x*sqx) - 1.0;
  double y   = sqrt(sqx);
  gsl_sf_result result_c;
  cheb_eval_mode_e(&aip_cs, z, mode, &result_c);
  result->val = (0.28125 + result_c.val)/y;
  result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  return GSL_SUCCESS;
}
int gsl_sf_airy_Ai_e(const double x, const gsl_mode_t mode, gsl_sf_result * result) {
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result cos_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_cos = gsl_sf_cos_err_e(theta.val, theta.err, &cos_result);
    result->val  = mod.val * cos_result.val;
    result->err  = fabs(mod.val * cos_result.err) + fabs(cos_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_cos);
  }
  else if(x <= 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs, z, mode, &result_c1);
    result->val  = 0.375 + (result_c0.val - x*(0.25 + result_c1.val));
    result->err  = result_c0.err + fabs(x*result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    double x32 = x * sqrt(x);
    double s   = exp(-2.0*x32/3.0);
    gsl_sf_result result_aie;
    int stat_aie = airy_aie(x, mode, &result_aie);
    result->val  = result_aie.val * s;
    result->err  = result_aie.err * s + result->val * x32 * GSL_DBL_EPSILON;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    CHECK_UNDERFLOW(result);
    return stat_aie;
  }
}
static int airy_bie(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  const double ATR =  8.7506905708484345;
  const double BTR = -2.0938363213560543;

  if(x < 4.0) {
    double sqx = sqrt(x);
    double z   = ATR/(x*sqx) + BTR;
    double y   = sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }
  else {
    double sqx = sqrt(x);
    double z   = 16.0/(x*sqx) - 1.0;
    double y   = sqrt(sqx);
    gsl_sf_result result_c;
    cheb_eval_mode_e(&bip2_cs, z, mode, &result_c);
    result->val = (0.625 + result_c.val)/y;
    result->err = result_c.err/y + GSL_DBL_EPSILON * fabs(result->val);
  }

  return GSL_SUCCESS;
}
int gsl_sf_airy_Bi_e(const double x, gsl_mode_t mode, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */
  if(x < -1.0) {
    gsl_sf_result mod;
    gsl_sf_result theta;
    gsl_sf_result sin_result;
    int stat_mp  = airy_mod_phase(x, mode, &mod, &theta);
    int stat_sin = gsl_sf_sin_err_e(theta.val, theta.err, &sin_result);
    result->val  = mod.val * sin_result.val;
    result->err  = fabs(mod.val * sin_result.err) + fabs(sin_result.val * mod.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_ERROR_SELECT_2(stat_mp, stat_sin);
  }
  else if(x < 1.0) {
    const double z = x*x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big_cs, z, mode, &result_c1);
    result->val  = 0.625 + result_c0.val + x*(0.4375 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x <= 2.0) {
    const double z = (2.0*x*x*x - 9.0)/7.0;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs, z, mode, &result_c1);
    result->val  = 1.125 + result_c0.val + x*(0.625 + result_c1.val);
    result->err  = result_c0.err + fabs(x * result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double y = 2.0*x*sqrt(x)/3.0;
    const double s = exp(y);

    if(y > GSL_LOG_DBL_MAX - 1.0) {
      OVERFLOW_ERROR(result);
    }
    else {
      gsl_sf_result result_bie;
      int stat_bie = airy_bie(x, mode, &result_bie);
      result->val  = result_bie.val * s;
      result->err  = result_bie.err * s + fabs(1.5*y * (GSL_DBL_EPSILON * result->val));
      result->err += GSL_DBL_EPSILON * fabs(result->val);
      return stat_bie;
    }
  }
}
void gsl_error (const char * reason, const char * file, int line, int gsl_errno) {
  if (gsl_error_handler) 
    {
      (*gsl_error_handler) (reason, file, line, gsl_errno);
      return ;
    }

  gsl_stream_printf ("ERROR", file, line, reason);

  fflush (stdout);
  fprintf (stderr, "Default GSL error handler invoked.\n");
  fflush (stderr);

  abort ();
}
void gsl_stream_printf (const char *label, const char *file, int line, 
                   const char *reason) {
  if (gsl_stream == NULL)
    {
      gsl_stream = stderr;
    }
  if (gsl_stream_handler)
    {
      (*gsl_stream_handler) (label, file, line, reason);
      return;
    }
  fprintf (gsl_stream, "gsl: %s:%d: %s: %s\n", file, line, label, reason);

}
int gsl_sf_cos_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_c = gsl_sf_cos_e(x, result);
  result->err += fabs(sin(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_c;
}
int gsl_sf_sin_err_e(const double x, const double dx, gsl_sf_result * result)
{
  int stat_s = gsl_sf_sin_e(x, result);
  result->err += fabs(cos(x) * dx);
  result->err += GSL_DBL_EPSILON * fabs(result->val);
  return stat_s;
}

static int cheb_eval_e(const cheb_series * cs,
                              const double x,
                              gsl_sf_result * result) {
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

int gsl_sf_cos_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = 1.0 - 0.5*x2;
      result->err = fabs(x2*x2/12.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = 1.0;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }

      if(octant > 1) {
        sgn_result = -sgn_result;
      }

      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}
int gsl_sf_sin_e(double x, gsl_sf_result * result) {
  /* CHECK_POINTER(result) */

  {
    const double P1 = 7.85398125648498535156e-1;
    const double P2 = 3.77489470793079817668e-8;
    const double P3 = 2.69515142907905952645e-15;

    const double sgn_x = GSL_SIGN(x);
    const double abs_x = fabs(x);

    if(abs_x < GSL_ROOT4_DBL_EPSILON) {
      const double x2 = x*x;
      result->val = x * (1.0 - x2/6.0);
      result->err = fabs(x*x2*x2 / 100.0);
      return GSL_SUCCESS;
    }
    else {
      double sgn_result = sgn_x;
      double y = floor(abs_x/(0.25*M_PI));
      int octant = y - ldexp(floor(ldexp(y,-3)),3);
      int stat_cs;
      double z;

      if(GSL_IS_ODD(octant)) {
        octant += 1;
        octant &= 07;
        y += 1.0;
      }

      if(octant > 3) {
        octant -= 4;
        sgn_result = -sgn_result;
      }
      
      z = ((abs_x - y * P1) - y * P2) - y * P3;

      if(octant == 0) {
        gsl_sf_result sin_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&sin_cs, t, &sin_cs_result);
        result->val = z * (1.0 + z*z * sin_cs_result.val);
      }
      else { /* octant == 2 */
        gsl_sf_result cos_cs_result;
        const double t = 8.0*fabs(z)/M_PI - 1.0;
        stat_cs = cheb_eval_e(&cos_cs, t, &cos_cs_result);
        result->val = 1.0 - 0.5*z*z * (1.0 - z*z * cos_cs_result.val);
      }

      result->val *= sgn_result;

      if(abs_x > 1.0/GSL_DBL_EPSILON) {
        result->err = fabs(result->val);
      }
      else if(abs_x > 100.0/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * abs_x * GSL_DBL_EPSILON * fabs(result->val);
      }
      else if(abs_x > 0.1/GSL_SQRT_DBL_EPSILON) {
        result->err = 2.0 * GSL_SQRT_DBL_EPSILON * fabs(result->val);
      }
      else {
        result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      }

      return stat_cs;
    }
  }
}

static int airy_deriv_mod_phase(const double x, gsl_mode_t mode,
                     gsl_sf_result * ampl, gsl_sf_result * phi) {
  const double pi34 = 2.356194490192344928847;
  gsl_sf_result result_a;
  gsl_sf_result result_p;
  double a, p;
  double sqx;
  double x32;

  if(x <= -4.0) {
    double z = 128.0/(x*x*x) + 1.0;
    cheb_eval_mode_e(&an20_cs_der, z, mode, &result_a);
    cheb_eval_mode_e(&aph0_cs_der, z, mode, &result_p);
  }
  else if(x <= -2.0) {
    double z = (128.0/(x*x*x) + 9.0) / 7.0;
    cheb_eval_mode_e(&an21_cs_der, z, mode, &result_a);
    cheb_eval_mode_e(&aph1_cs_der, z, mode, &result_p);
  }
  else if(x <= -1.0) {
    double z = (16.0/(x*x*x) + 9.0) / 7.0;
    cheb_eval_mode_e(&an22_cs_der, z, mode, &result_a);
    cheb_eval_mode_e(&aph2_cs_der, z, mode, &result_p);
  }
  else {
    ampl->val = 0.0;
    ampl->err = 0.0;
    phi->val  = 0.0;
    phi->err  = 0.0;
    GSL_ERROR ("x is greater than 1.0", GSL_EDOM);
  }

  a =  0.3125 + result_a.val;
  p = -0.625  + result_p.val;
 
  sqx = sqrt(-x);
  x32   = x*sqx;

  ampl->val = sqrt(a * sqx);
  ampl->err = fabs(ampl->val) * (GSL_DBL_EPSILON + fabs(result_a.err/result_a.val));
  phi->val  = pi34 - x * sqx * p;
  phi->err = fabs(phi->val) * (GSL_DBL_EPSILON + fabs(result_p.err/result_p.val));

  return GSL_SUCCESS;
}

double gsl_sf_airy_Ai_deriv(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Ai_deriv_e(x, mode, &result));
}
int gsl_sf_airy_Ai_deriv_e(const double x, gsl_mode_t mode, 
                           gsl_sf_result * result) {
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int status_ap = airy_deriv_mod_phase(x, mode, &a, &p);
    double c    = cos(p.val);
    result->val  = a.val * c;
    result->err  = fabs(result->val * p.err) + fabs(c * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const double x3 = x*x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&aif_cs_der, x3, mode, &result_c1);
    cheb_eval_mode_e(&aig_cs_der, x3, mode, &result_c2);
    result->val  = (x*x*(0.125 + result_c1.val) - result_c2.val) - 0.25;
    result->err  = fabs(x*x*result_c1.err) + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x*x*x < 9.0/4.0 * GSL_LOG_DBL_MIN*GSL_LOG_DBL_MIN) {
    gsl_sf_result result_aps;
    const double arg = -2.0*x*sqrt(x)/3.0;
    const int stat_a = gsl_sf_airy_Ai_deriv_scaled_e(x, mode, &result_aps);
    const int stat_e = gsl_sf_exp_mult_err_e(arg, 1.5*fabs(arg*GSL_DBL_EPSILON),
                                                result_aps.val, result_aps.err,
                                                result);
    return GSL_ERROR_SELECT_2(stat_e, stat_a);
  }
  else {
    UNDERFLOW_ERROR(result);
  }
}

int gsl_sf_airy_Ai_deriv_scaled_e(const double x, gsl_mode_t mode, gsl_sf_result * result) {
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int status_ap = airy_deriv_mod_phase(x, mode, &a, &p);
    double c    = cos(p.val);
    result->val  = a.val * c;
    result->err  = fabs(result->val * p.err) + fabs(c * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x <= 1.0) {
    const double x3 = x*x*x;
    const double x2 = x*x;
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&aif_cs_der, x3, mode, &result_c0);
    cheb_eval_mode_e(&aig_cs_der, x3, mode, &result_c1);

    result->val  = (x2*(0.125 + result_c0.val) - result_c1.val) - 0.25;
    result->err  = fabs(x2*result_c0.val) + result_c1.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > GSL_ROOT3_DBL_EPSILON*GSL_ROOT3_DBL_EPSILON) {
      /* scale only if x is positive */
      double s = exp(2.0*x*sqrt(x)/3.0);
      result->val *= s;
      result->err *= s;
    }

    return GSL_SUCCESS;
  }
  else if(x <= 4.0) {
    const double sqrtx = sqrt(x);
    const double z = (16.0/(x*sqrtx) - 9.0)/7.0;
    const double s = sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&aip1_cs_der, z, mode, &result_c0);
    result->val  = -(0.28125 + result_c0.val) * s;
    result->err  = result_c0.err * s;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sqrtx = sqrt(x);
    const double z = 16.0/(x*sqrtx) - 1.0;
    const double s = sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&aip2_cs_der, z, mode, &result_c0);
    result->val  = -(0.28125 + result_c0.val) * s;
    result->err  = result_c0.err * s;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
int gsl_sf_exp_mult_err_e(const double x, const double dx,
                             const double y, const double dy,
                             gsl_sf_result * result)
{
  const double ay  = fabs(y);

  if(y == 0.0) {
    result->val = 0.0;
    result->err = fabs(dy * exp(x));
    return GSL_SUCCESS;
  }
  else if(   ( x < 0.5*GSL_LOG_DBL_MAX   &&   x > 0.5*GSL_LOG_DBL_MIN)
          && (ay < 0.8*GSL_SQRT_DBL_MAX  &&  ay > 1.2*GSL_SQRT_DBL_MIN)
    ) {
    double ex = exp(x);
    result->val  = y * ex;
    result->err  = ex * (fabs(dy) + fabs(y*dx));
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double ly  = log(ay);
    const double lnr = x + ly;

    if(lnr > GSL_LOG_DBL_MAX - 0.01) {
      OVERFLOW_ERROR(result);
    }
    else if(lnr < GSL_LOG_DBL_MIN + 0.01) {
      UNDERFLOW_ERROR(result);
    }
    else {
      const double sy  = GSL_SIGN(y);
      const double M   = floor(x);
      const double N   = floor(ly);
      const double a   = x  - M;
      const double b   = ly - N;
      const double eMN = exp(M+N);
      const double eab = exp(a+b);
      result->val  = sy * eMN * eab;
      result->err  = eMN * eab * 2.0*GSL_DBL_EPSILON;
      result->err += eMN * eab * fabs(dy/y);
      result->err += eMN * eab * fabs(dx);
      return GSL_SUCCESS;
    }
  }
}
double gsl_sf_airy_Bi_deriv_scaled(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_deriv_scaled_e(x, mode, &result));
}
double gsl_sf_airy_Bi_deriv(const double x, gsl_mode_t mode)
{
  EVAL_RESULT(gsl_sf_airy_Bi_deriv_e(x, mode, &result));
}
int gsl_sf_airy_Bi_deriv_e(const double x, gsl_mode_t mode, 
                           gsl_sf_result * result){
  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int status_ap = airy_deriv_mod_phase(x, mode, &a, &p);
    double s    = sin(p.val);
    result->val  = a.val * s;
    result->err  = fabs(result->val * p.err) + fabs(s * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const double x3 = x*x*x;
    const double x2 = x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif_cs_der, x3, mode, &result_c1);
    cheb_eval_mode_e(&big_cs_der, x3, mode, &result_c2);
    result->val  = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x2 * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 2.0) {
    const double z = (2.0*x*x*x - 9.0) / 7.0;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif2_cs_der, z, mode, &result_c1);
    cheb_eval_mode_e(&big2_cs_der, z, mode, &result_c2);
    result->val  = x*x * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x*x * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < GSL_ROOT3_DBL_MAX*GSL_ROOT3_DBL_MAX) {
    gsl_sf_result result_bps;
    const double arg = 2.0*(x*sqrt(x)/3.0);
    int stat_b = gsl_sf_airy_Bi_deriv_scaled_e(x, mode, &result_bps);
    int stat_e = gsl_sf_exp_mult_err_e(arg, 1.5*fabs(arg*GSL_DBL_EPSILON),
                                          result_bps.val, result_bps.err,
                                          result);
    return GSL_ERROR_SELECT_2(stat_e, stat_b);
  }
  else {
    OVERFLOW_ERROR(result);
  }
}
int gsl_sf_airy_Bi_deriv_scaled_e(const double x, gsl_mode_t mode, 
                                  gsl_sf_result * result) {
  const double atr =  8.7506905708484345;   /* 16./(sqrt(8)-1) */
  const double btr = -2.0938363213560543;   /* -(sqrt(8)+1)/(sqrt(8)-1) */

  /* CHECK_POINTER(result) */

  if(x < -1.0) {
    gsl_sf_result a;
    gsl_sf_result p;
    int status_ap = airy_deriv_mod_phase(x, mode, &a, &p);
    double s     = sin(p.val);
    result->val  = a.val * s;
    result->err  = fabs(result->val * p.err) + fabs(s * a.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return status_ap;
  }
  else if(x < 1.0) {
    const double x3 = x*x*x;
    const double x2 = x*x;
    gsl_sf_result result_c1;
    gsl_sf_result result_c2;
    cheb_eval_mode_e(&bif_cs_der, x3, mode, &result_c1);
    cheb_eval_mode_e(&big_cs_der, x3, mode, &result_c2);
    result->val  = x2 * (result_c1.val + 0.25) + result_c2.val + 0.5;
    result->err  = x2 * result_c1.err + result_c2.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);

    if(x > GSL_ROOT3_DBL_EPSILON*GSL_ROOT3_DBL_EPSILON) {
      /* scale only if x is positive */
      const double s = exp(-2.0*x*sqrt(x)/3.0);
      result->val *= s;
      result->err *= s;
    }

    return GSL_SUCCESS;
  }
  else if(x < 2.0) {
    const double z = (2.0*x*x*x - 9.0) / 7.0;
    const double s = exp(-2.0*x*sqrt(x)/3.0);
    gsl_sf_result result_c0;
    gsl_sf_result result_c1;
    cheb_eval_mode_e(&bif2_cs_der, z, mode, &result_c0);
    cheb_eval_mode_e(&big2_cs_der, z, mode, &result_c1);
    result->val  = s * (x*x * (0.25 + result_c0.val) + 0.5 + result_c1.val);
    result->err  = s * (x*x * result_c0.err + result_c1.err);
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x < 4.0) {
    const double sqrtx = sqrt(x);
    const double z = atr/(x*sqrtx) + btr;
    const double s = sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&bip1_cs_der, z, mode, &result_c0);
    result->val  = s * (0.625 + result_c0.val);
    result->err  = s * result_c0.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else {
    const double sqrtx = sqrt(x);
    const double z = 16.0/(x*sqrtx) - 1.0;
    const double s = sqrt(sqrtx);
    gsl_sf_result result_c0;
    cheb_eval_mode_e(&bip2_cs_der, z, mode, &result_c0);
    result->val  = s * (0.625 + result_c0.val);
    result->err  = s * result_c0.err;
    result->err += GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
}
gsl_complex gsl_complex_add (gsl_complex a, gsl_complex b) {
  double ar = GSL_REAL (a), ai = GSL_IMAG (a);
  double br = GSL_REAL (b), bi = GSL_IMAG (b);

  gsl_complex z;
  GSL_SET_COMPLEX (&z, ar + br, ai + bi);
  return z;
}

