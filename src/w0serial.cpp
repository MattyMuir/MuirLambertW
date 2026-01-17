#include <cmath>

static inline double AddEm(double x)
{
    static constexpr double emHigh = 0.36787944117144232160;
    static constexpr double emLow = -1.2428753672788363168e-17;

    return (x + emHigh) + emLow;
}

// [EM, -0.2]
static inline double Approx1(double x)
{
    // Rational approximation coefficients for algorithm 1, index 1, order 6/6
    static constexpr double P[] = {
        -0.999999999999999998781454,-2.218532275103755425409341,1.09193891725289689417007,5.676164924352752669529278,4.333697663393533268812102,1.069483295881317351836729,0.06209728006746276596365885
    };
    static constexpr double Q[] = {
        1,4.550176256700878418784205,7.705264281250672121011024,5.98062498820306285618724,2.106065724108968678807767,0.2831366693201161636277512,0.008433340110008441063524488
    };

    double t = sqrt(AddEm(x));

    double numer = P[6];
    for (size_t i = 0; i < 6; i++)
        numer = numer * t + P[5 - i];

    double denom = Q[6];
    for (size_t i = 0; i < 6; i++)
        denom = denom * t + Q[5 - i];

    return numer / denom;
}

// [-0.2, 0.29]
static inline double Approx2(double x)
{
    // Rational approximation coefficients for algorithm 1, index 2, order 8/7
    static constexpr double P[] = {
        4.69506e-15,-1.3032001e-13,-14.62923727738790490843,-112.0311635820189938392,-319.2626531507523886724,-414.3459534818585609992,-237.7255507343262504594,-47.11464158808262121342,-1.
    };
    static constexpr double Q[] = {
        29.258474554755502910434,263.07362657044956799772,923.4585739869901129263,1592.8905039245497144734,1389.990974014904455648,561.71139116540236682278,79.0696738395877184621,1.
    };

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * x + P[7 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * x + Q[6 - i];

    return x / (numer / denom + 1.0 + x);
}

// [0.29, 4]
static inline double Approx3(double x)
{
    // Rational approximation coefficients for algorithm 1, index 3, order 9/8
    static constexpr double P[] = {
        -440.6731001967694598229,-1678.822471546733977823,-504.690210255024449673,4255.73984382946080597,5102.289142949515499772,1480.65590448134185755,-246.4351557243390048574,-138.0298996206123389102,-12.78796032476732884865,-0.2448638517438543484516
    };
    static constexpr double Q[] = {
        1762.6924155049607541829,12003.366679806332641726,30978.09793714236012718,38473.957611820841411168,24318.929508956407826684,7691.197450762467632937,1128.2558170022085444327,65.173884390065698879591,1.
    };

    double numer = P[9];
    for (size_t i = 0; i < 9; i++)
        numer = numer * x + P[8 - i];

    double denom = Q[8];
    for (size_t i = 0; i < 8; i++)
        denom = denom * x + Q[7 - i];

    return x * 0.25 + 0.25 + numer / denom;
}

// [4, 35]
static inline double Approx4(double x)
{
    // Rational approximation coefficients for algorithm 1, index 4, order 7/7
    static constexpr double P[] = {
        371352.6890273707190189,764816.1289134941936569,1.590788999468928463328e6,272562.4960524054950652,-143824.95152375100549954,-20840.30840375205222913,-206.8709250331294554843,7.110575694956416395276
    };
    static constexpr double Q[] = {
        -191680.7497864744754626,-430288.66700330475150672,-840876.8188997657662196,-550849.70537404649021723,-105906.43155790512238976,-4989.754151281373681581,19.133153468520517672461,1.
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 1.9013671875 + numer / denom;
}

// [35, 400]
static inline double Approx5(double x)
{
    // Rational approximation coefficients for algorithm 1, index 5, order 7/7
    static constexpr double P[] = {
        -1.792881327664732774125672,-17.61464238625069531482495,-9.46468505228302440381521,-0.098507032250218364378744,0.252684589808521844673978,0.0185401997044637489930649,0.000320124383079910459803588,1.00596762044525170188858e-6
    };
    static constexpr double Q[] = {
        1.,5.634133203100493766531353,5.701846569650484362210726,1.62801389384192401912878,0.1587858940217511099807044,0.005384586308983632735195283,0.00005485573193546679034405104,1.000471535629899062777137e-7
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 2.599609375 + numer / denom;
}

// [400, 7225]
static inline double Approx6(double x)
{
    // Rational approximation coefficients for algorithm 1, index 6, order 7/7
    static constexpr double P[] = {
        -5.005623424227718509862604,-1.726465528305267763227662,-0.0930465095885677277099738,0.00345136751278117870639794,0.00023925068112691758179999,3.176064297510035854535192e-6,1.150279068130744295450373e-8,7.97525056173957603661583e-12
    };
    static constexpr double Q[] = {
        1.,0.6238687213326071358573221,0.09597659788456336449777539,0.005083872187709157313884376,0.0001019391699294342050514312,7.534616420138851318126195e-7,1.734589555390420938538666e-9,7.324840574097317633321993e-13
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 4.4892578125 + numer / denom;
}

// [7225, 115600]
static inline double Approx7(double x)
{
    // Rational approximation coefficients for algorithm 1, index 7, order 7/7
    static constexpr double P[] = {
        -6.80850916801115975512343,-0.937716606754753009171354,-0.019651274882844411834533,0.000073241085858857224564555,2.58534440356849950813133e-6,1.012062876528937260891374e-8,1.012302477751943323736401e-11,1.884509485121327887975332e-15
    };
    static constexpr double Q[] = {
        1.,0.2548503903981785486641877,0.01304323120742714196353544,0.0002089306300337398523361538,1.200229358904857318594755e-6,2.456451693171267057402716e-9,1.530533070068020571046581e-12,1.722166998754853949457752e-16
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 6.9462890625 + numer / denom;
}

// [115600, 1960000]
static inline double Approx8(double x)
{
    // Rational approximation coefficients for algorithm 1, index 8, order 7/7
    static constexpr double P[] = {
        -7.612093884341071008010845,-0.274802792714947248733489,-0.00145745680927150862083744,1.35999973363915456373421e-6,1.19196933211232836882692e-8,1.155229958994611342409266e-11,2.8532804068730851773796e-15,1.308453166623186115267132e-19
    };
    static constexpr double Q[] = {
        1.,0.07003385475422981402026354,0.0009228901402300251218971713,3.730855899369000609795312e-6,5.359968694684899993961814e-9,2.729613387105260895868582e-12,4.218623564427737659800048e-16,1.175343889812314173725193e-20
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 9.4150390625 + numer / denom;
}

// [1960000, 33640000]
static inline double Approx9(double x)
{
    // Rational approximation coefficients for algorithm 1, index 9, order 7/7
    static constexpr double P[] = {
        -8.027947460271808480269167,-0.0712181706141578815996098,-0.000092002807890987479538579,2.0851384228527376946667e-8,4.43163923372767059192999e-11,1.040659663685645368914985e-14,6.22354864600744979061285e-19,6.90567068355775459798299e-24
    };
    static constexpr double Q[] = {
        1.,0.01756139045708885178643158,0.00005678555412565004401422836,5.598424227113294734502591e-8,1.956262231361565891562691e-11,2.419677496364215364069321e-15,9.075611079244018685539657e-20,6.135139434681698755303782e-25
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 12.0029296875 + numer / denom;
}

// [33640000, 625000000]
static inline double Approx10(double x)
{
    // Rational approximation coefficients for algorithm 1, index 10, order 7/7
    static constexpr double P[] = {
        -8.24084249505589349742764,-0.0173640630156215691060427,-5.26165875796207160722595e-6,3.1501370303779933098242e-10,1.47472083420334075784529e-13,8.13083425756224047368241e-18,1.144915848510323444807518e-22,2.99026302450106453943181e-28
    };
    static constexpr double Q[] = {
        1.,0.004228645260524294370495129,3.260233706589085539503716e-6,7.637479264024251424078956e-10,6.328507995588061328180158e-14,1.853386497214566930267771e-18,1.643953680433643000010534e-23,2.625816506484210095089451e-29
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 14.646484375 + numer / denom;
}

// [625000000, 1.2321e10]
static inline double Approx11(double x)
{
    // Rational approximation coefficients for algorithm 1, index 11, order 7/7
    static constexpr double P[] = {
        -8.3751159851081905254953,-0.00403793316970308814262656,-2.77539139953009958918235e-7,4.110797040983118598181e-12,4.182908162231595424257e-16,5.23980765061453367142234e-21,1.67992781689060350287687e-26,9.98781622714724674922145e-33
    };
    static constexpr double Q[] = {
        1.,0.000975717886293310054664437,1.72701179442394066917732e-7,9.269739055908127656955358e-12,1.757663902508001236141819e-16,1.176735106854517192845353e-21,2.384034220840144864074194e-27,8.692224140447244647516727e-34
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 17.396484375 + numer / denom;
}

// [1.2321e10, 2.45025e11]
static inline double Approx12(double x)
{
    // Rational approximation coefficients for algorithm 1, index 12, order 7/7
    static constexpr double P[] = {
        -8.47904179929216848163146,-0.00091920288370704052148185,-1.41801707419537036755195e-8,4.7815945158179396523238e-14,1.08443748230076127789288e-18,3.04937459252979197979618e-24,2.19533456384229817212162e-30,2.93063227712280141389727e-37
    };
    static constexpr double Q[] = {
        1.,0.0002201916057056314181424858,8.773797493718583179902109e-9,1.059355403941572180156795e-13,4.516712623402576062514378e-19,6.797879437499254827837037e-25,3.095660300044424090738299e-31,2.537025865290478330568713e-38
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 20.2265625 + numer / denom;
}

// [2.45025e11, 5.0625e12]
static inline double Approx13(double x)
{
    // Rational approximation coefficients for algorithm 1, index 13, order 7/7
    static constexpr double P[] = {
        -8.54252631199205009492691,-0.000205663328564006721917119,-7.0089605019484734334422e-10,5.5087419705783248519153e-16,2.69267545543350657250446e-21,1.6751649979420392588839e-27,2.67173520967342607675829e-34,7.90033312046015246921766e-42
    };
    static constexpr double Q[] = {
        1.,0.00004909267637852783469934642,4.35322411173719211708329e-10,1.16865385223501684593348e-15,1.107125349746267653108508e-21,3.700216644642381656609618e-28,3.7399451062026616043743e-35,6.800307233759839460468346e-43
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 23.0849609375 + numer / denom;
}

// [5.0625e12, 1.0609e14]
static inline double Approx14(double x)
{
    // Rational approximation coefficients for algorithm 1, index 14, order 7/7
    static constexpr double P[] = {
        -8.59542220139790648743983,-0.0000453592241753706977451745,-3.3809176602658933735261e-11,5.9390755599427177205496e-18,6.2878487229374374461311e-24,8.5613017313268683778171e-31,2.99005646037029901845583e-38,1.93604238789417158143534e-46
    };
    static constexpr double Q[] = {
        1.,0.00001078368206327979150689663,2.098243755937184929063929e-11,1.235481492793517848680097e-17,2.566406645809936923832783e-24,1.880327185598498590960048e-31,4.165515027248826798458935e-39,1.659928415472216380340524e-47
    };

    double t = sqrt(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return 25.994140625 + numer / denom;
}

// [1.0609e14, 5.9209720277e47]
static inline double Approx15(double x)
{
    // Rational approximation coefficients for algorithm 1, index 15, order 7/6
    static constexpr double P[] = {
        0.7432598296406762110144502,-0.5952760274551788492906215,-0.0668761716145358047390643,-0.0014678068208685779597395,-4.085792684429854816067e-6,4.28068948312128597983008e-8,1.0199481223789291568498e-10,5.81829262478192878862e-15
    };
    static constexpr double Q[] = {
        1.,0.3214298794122304285047461,0.01961181540054920446844377,0.0002951873837001442714238258,2.559254266891386662646856e-7,-8.188206603754071240509357e-9,-1.334261865436649758928393e-11
    };

    double t = log(x);

    double numer = P[7];
    for (size_t i = 0; i < 7; i++)
        numer = numer * t + P[6 - i];

    double denom = Q[6];
    for (size_t i = 0; i < 6; i++)
        denom = denom * t + Q[5 - i];

    return t + numer / denom;
}

// [5.9209720277e47, DBL_MAX]
static inline double Approx16(double x)
{
    // Rational approximation coefficients for algorithm 1, index 16, order 8/7
    static constexpr double P[] = {
        0.3654395923466255825690192,-0.5171047620893658168387926,-0.0337328495559685441962204,-0.00054759067987015982542188,-3.0334321402991061053707e-6,-6.125076497093456058971e-9,-4.179206569160681749502e-12,-7.08445931807102286001e-16,-3.6926689468640405451e-21
    };
    static constexpr double Q[] = {
        1.,0.220533378532386649632492,0.008943830909731700163501178,0.0001125724579079293206308385,5.180133726226615474763108e-7,8.910908966615756121574719e-10,5.177785017654958872034204e-13,7.165956863781358128164073e-17
    };

    double t = log(x);

    double numer = P[8];
    for (size_t i = 0; i < 8; i++)
        numer = numer * t + P[7 - i];

    double denom = Q[7];
    for (size_t i = 0; i < 7; i++)
        denom = denom * t + Q[6 - i];

    return t + numer / denom;
}

double MuirW0(double x)
{
    if (x == INFINITY) return INFINITY;

    if (x < 1960000)
    {
        if (x < 35)
        {
            if (x < 0.29)
            {
                if (x < -0.2) return Approx1(x);
                return Approx2(x);
            }

            if (x < 4) return Approx3(x);
            return Approx4(x);
        }

        if (x < 7225)
        {
            if (x < 400) return Approx5(x);
            return Approx6(x);
        }

        if (x < 115600) return Approx7(x);
        return Approx8(x);
    }

    if (x < 2.45025e11)
    {
        if (x < 625000000)
        {
            if (x < 33640000) return Approx9(x);
            return Approx10(x);
        }

        if (x < 1.2321e10) return Approx11(x);
        return Approx12(x);
    }

    if (x < 1.0609e14)
    {
        if (x < 5.0625e12) return Approx13(x);
        return Approx14(x);
    }

    if (x < 5.9209720277e47) return Approx15(x);
    return Approx16(x);
}