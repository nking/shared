package thirdparty.net.oelen.polsolve.pzeros;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.List;

/**
 *
 * @author woelen
 */
public class TestPolynomials {
    // Test polynomials can be specified by means of a two-letter string and an
    // array of numerical values:
    //   RR : Roots Real, all numbers are real roots. Numbers separated by ; stand for
    //        a complex conjugate pair. For example, the polynomial specified as
    //        {"RR", "1", "-2", "3;4", "5"} is the polynomial (z-1)(z+2)(z-3-4i)(z-3+4i)(z-5)
    //   CR : Coefs Real, numbers are coefficients. For example, the polynomial specified as
    //        {"CR", "1", "2", "3", "0", "-5"} is the polynomial 1 + 2*z + 3*z^2 - 5*z^4
    //   RC : Roots Complex, number strings stand for real roots, Numbers, separated by ;
    //        stand for a single complex root. For example, the polynomial, specified as
    //        {"RC", "1;1", "2", "3;-2", "-3") is the polynomial (z-1-i)(z-2)(z-3+2i)(z+3)
    //   CC : Coefs Complex. Numbers are coefficients. For example, the polynomial specified as
    //        {"CC", "1;1", "0", "-2;-4", "8", "0;5"} is 1+i + (-2-4i)*z^2 + 8*z^3 + 5i*z^4
    // 
    // You can add new polynomials, just by name, and they will automatically be recognized
    // by the getPolynomials() method, given below. This property allows very easy extension
    // of the test suite for the polynomial solvers.
    
    
    // A seemingly random polynomial, which is a hard one for many solvers.
    // It also makes things quite difficult for the PZeros solver!
    static 
    public String[] WOELEN0 = {"RR", "4.52", "2.37", "-4.68",
        "2.25",
        "2.19",
        "-3.69",
        "-0.92",
        "-1.15",
        "-1.49",
        "-0.42",
        "-4.64",
        "2.00",
        "-2.54",
        "-3.36",
        "-3.50",
        "-3.85",
        "0.06",
        "4.50",
        "2.20",
        "-3.46",
        "4.48",
        "-3.41",
        "-4.15",
        "3.08",
        "-3.96;1.13",
        "-1.32;1.77",
        "3.44;0.32"};
    
    // Just a simple polynomial, which tests the functionality of the
    // software. If this does not give good results, then the software
    // to be tested is broken.
    static
    public String[] WOELEN1 = {"RR", "3;4", "1", "2", "4;5", "-3;4", "-6;7"};
    
    // The following polynomial has 23 roots, all of them exact multiples of
    // 0.01 with a deviation of at most 1 ppb. This polynomial is ill-conditioned.
    static
    public String[] WOELEN2 = {"CR",
        "-4.045012939518347816847920585748e10",
        "-1.888087406087532166014161339279e11",
        "-3.295125720566197014641938085859e11",
        "-2.250119152904938422148671897994e11",
        "4.942583110950258363757447196915e10",
        "1.850106964655267571073048912017e11",
        "1.149520963469206315579723066639e11",
        "9951381603.182799281310727107729",
        "-2.596480853186198601189459542701e10",
        "-1.623619660766524990709394417409e10",
        "-3289311362.176895609539780253830",
        "1045774648.353365050544923692971",
        "924637650.8764884154837085833074",
        "267920883.0027485899226404130771",
        "12186397.19954828202744702780974",
        "-19714940.83023376923476613354756",
        "-8334224.243558339590157089007706",
        "-1542880.590688567989365630921971",
        "16799.64105418819220397255403594",
        "90981.74341258074200094916007031",
        "25730.30188928914633538149677721",
        "3734.644325747884227019977959842",
        "291.4318471130768686177267609720",
        "9.698231185127349718300138192717"
    };
    
    
    // The 10-th degree Wilkinson polynomial. It is quite ill-conditioned,
    // but all polynomials solvers should have no problems with this when
    // standard 53 bit arithmetic precision is used.
    static 
    public String[] WILK10 = {"RR", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    
    // Another extremely ill-conditioned polynomial, which was analysed by
    // Wilkinson in the 1960's already. The root equal to 20 has a sensitivity
    // for changes in the coefficients which is 20^19 times as large as the
    // sensitivity of the root equal to 1. Many polynomial solvers, which
    // only have 53 bit precision available have a hard time with this.
    static 
    public String[] WILKINSON = {"RR", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"};

    static 
    public String[] JT3 = {"RR", "1", "0.1", "0.01", "0.001", "0.0001", "0.00001", "1e-6", "1e-7", "1e-8", "1e-9", "1e-10"};
    
    static
    public String[] JT4 = {"RR", "0.1", "0.1", "0.1", "0.5", "0.6", "0.7"};
    
    static
    public String[] JT5 = {"RR", "0.1", "0.1", "0.1", "0.1", "0.2", "0.2", "0.2", "0.3", "0.3", "0.4"};
    
    static
    public String[] JT6 = {"RR", "0.1", "1.001", "0.998", "1.00002", "0.99999"};
    
    static
    public String[] JT7a = {"RR", "0.001", "0.01", "0.1", "0.1;0.00001", "1", "10"};
    
    static
    public String[] JT7b = {"RR", "0.001", "0.01", "0.1", "0.1;-0.000001", "0.1;0.000001", "1", "10"};
    
    static
    public String[] JT7c = {"RR", "0.001", "0.01", "0.1", "0.1;-0.0000001", "0.1;0.0000001", "1", "10"};
    
    static
    public String[] JT8 = {"RR", "1.0002", "1.0002", "1.0002", "1.0002", "1.0002", "3.9998", "4.0002"};

    
    /* JT9 polynomial roots:
        -9.5105651629515357211       -3.090169943749474241
        -9.5105651629515357211        3.090169943749474241
        -5.8778525229247312917       -8.090169943749474241
        -5.8778525229247312917        8.090169943749474241
        -0.095105651629515357212      0.03090169943749474241
        -0.095105651629515357212     -0.03090169943749474241
        -0.058778525229247312917      0.08090169943749474241
        -0.058778525229247312917     -0.08090169943749474241
         0                          -10
         0                           10
         0                           -0.1
         0                            0.1
         0.058778525229247312917      0.08090169943749474241
         0.058778525229247312917     -0.08090169943749474241
         0.095105651629515357212     -0.03090169943749474241
         0.095105651629515357212      0.03090169943749474241
         5.8778525229247312917        8.090169943749474241
         5.8778525229247312917       -8.090169943749474241
         9.5105651629515357211       -3.090169943749474241
         9.5105651629515357211        3.090169943749474241
    */
    static
    public String[] JT9 = {"CR", "1",                     "0", "0", "0", "0", "0", "0", "0", "0", "0",
                                 "9999999999.9999999999", "0", "0", "0", "0", "0", "0", "0", "0", "0",
                                 "1"};    // (z¹⁰ - 1/10¹⁰)(z¹⁰ + 10¹⁰)

    
    /* MIGNOTTE1 polynomial roots:
        -1.6416435693506390758        -0.28966326727303548058
        -1.6416435693506390758         0.28966326727303548058
        -1.4435031163886131792        -0.83405213461823124326
        -1.4435031163886131792         0.83405213461823124326
        -1.0711208704931316813        -1.2778419997014787862
        -1.0711208704931316813         1.2778419997014787862
        -0.56941162198766528743        1.567505254453898159
        -0.56941162198766528743       -1.567505254453898159
        -0.01                          1E-22
        -0.01                         -1E-22
         0.0011110930358771479945      1.6681042376166902171
         0.0011110930358771479945     -1.6681042376166902171
         0.57163381651625665665        1.5675052312172674459
         0.57163381651625665665       -1.5675052312172674459
         1.0733430864364517267        -1.2778419640993868033
         1.0733430864364517267         1.2778419640993868033
         1.4457253566860692764        -0.83405210330876849706
         1.4457253566860692764         0.83405210330876849706
         1.643865825545394416         -0.28966325490757405807
         1.643865825545394416          0.28966325490757405807
    */
    static
    public String[] MIGNOTTE1 = {"CR", "1", "200", "10000", "0", "0", "0", "0", "0", "0", "0",
                                       "0",   "0",     "0", "0", "0", "0", "0", "0", "0", "0",
                                       "1"};   // x²⁰ + (100x+1)²    
    
    /* MIGNOTTE2 polynomial roots:
        -2.2521645747763735138          0
        -2.099961758981865955          -0.81421651543673016602
        -2.099961758981865955           0.81421651543673016602
        -1.6639091433471375182         -1.5184685738528918752
        -1.6639091433471375182          1.5184685738528918752
        -1.0028980451623218692         -2.0176430300764448855
        -1.0028980451623218692          2.0176430300764448855
        -0.20620167103887976323         2.244323609267870081
        -0.20620167103887976323        -2.244323609267870081
        -0.01                           0
        -0.01                           0
        -0.01                           0
         0.61858172170157911977        -2.167895847573263897
         0.61858172170157911977         2.167895847573263897
         1.3600605645541593242         -1.7986817427491601583
         1.3600605645541593242          1.7986817427491601583
         1.9180940269885225814         -1.1865457085548345649
         1.9180940269885225814          1.1865457085548345649
         2.2173165926741308372         -0.41416010784257236963
         2.2173165926741308372          0.41416010784257236963
    */
    static
    public String[] MIGNOTTE2 = {"CR", "1", "300", "30000", "1000000", "0", "0", "0", "0", "0", "0",
                                       "0",   "0",     "0",       "0", "0", "0", "0", "0", "0", "0",
                                       "1"};   // x²⁰ + (100x+1)³
    
    
    /* MIGNOTTE3 polynomial roots:
        -2.613287409914964481          -0.59695706055660041996
        -2.613287409914964481           0.59695706055660041996
        -2.0952693780634052963         -1.6726365104206022639
        -2.0952693780634052963          1.6726365104206022639
        -1.1618331320385936599         -2.4170297764504124058
        -1.1618331320385936599          2.4170297764504124058
        -0.010000000000023207944       -4.0197338439365768428E-14
        -0.010000000000023207944        4.0197338439365768428E-14
        -0.01                           0
        -0.01                           0
        -0.01                           0
        -0.0099999999999535841117       0
         0.002142839926418955442        2.6827006448820861499
         0.002142839926418955442       -2.6827006448820861499
         1.1661188248551808005         -2.4170297495298379895
         1.1661188248551808005          2.4170297495298379895
         2.0995551000106903299         -1.672636476850064341
         2.0995551000106903299          1.672636476850064341
         2.6175731552246733513         -0.59695704561590193993
         2.6175731552246733513          0.59695704561590193993    
    */
    static 
    public String[] MIGNOTTE3 = {"CR", "1", "600", "150000", "20000000", "1500000000", "60000000000",
                                       "1000000000000", "0",   "0",     "0",       "0", "0", "0", "0", "0", "0", "0",
                                       "1", "300", "30000", "1000000"};   // (x¹⁷+(100x+1)³)(100x+1)³
    
    static
    public String[] KAMENY1 = {"CC","9e-24", "-6e-12", "1", "0", "0", "0", "0", "0;1e-6"}; // (z-3c²)² + i*cx⁷  (with c = 1e-6)
    
    static
    public String[] KAMENY2 = {"CC","9", "0", "-6e12", "0", "1e24", "0", "0", "0", "0", "0;1e12"}; // (c²z²-3)² + i*c²x⁹  (with c = 1e6)
    
    
    
    /* KAMENY3 polynomial roots:
        -251.18864315095800633                   0
         -77.62155952763026605                -238.89459588805661795
         -77.62155952763026605                 238.89459588805661795
          -0.17320508075688772935E-5             0
          -0.17320508075688772935E-5             0
           0.17320508075688772935E-5             0
           0.17320508075688772935E-5             0
         203.21588110310926922                -147.64497998748985988
         203.21588110310926922                 147.64497998748985988
    */
    static
    public String[] KAMENY3 = {"CR","9", "0", "-6e12", "0", "1e24", "0", "0", "0", "0", "1e12"}; // (c²z²-3)² + c²x⁹  (with c = 1e6)
    
    
    /* KAMENY4 polynomial roots:
        -1E8                            0
        -1E8                            0
        -0.11892071150027210667E-5      0
        -0.11892071150027210667E-5      0
         0                             -0.11892071150027210667E-5
         0                             -0.11892071150027210667E-5
         0                              0.11892071150027210667E-5
         0                              0.11892071150027210667E-5
         0.11892071150027210667E-5      0
         0.11892071150027210667E-5      0
         5E7                           -8.6602540378443864676E7
         5E7                           -8.6602540378443864676E7
         5E7                            8.6602540378443864676E7
         5E7                            8.6602540378443864676E7
    */
    static
    public String[] KAMENY4 = {"CR", "4", "0", "0", "0", "-4e24", "0", "0", "4", "1e48", "0", "0", "2e24", "0", "0", "1"};

    
    
    
    
    // The code below is for retrieving the polynomials, specified above
    // in this file. It can be used for obtaining the names and the coefficients
    // of the polynomials.
    
    public static void getPolynomials(List<String> nameList, List<String[]> polList) {
        nameList.clear();
        polList.clear();    
        try {
            Field[] fields = TestPolynomials.class.getFields();
            for (Field f : fields) {
                Class<?> clazz = f.getType();
                if (!clazz.isArray()) {
                    continue;
                }
                clazz = clazz.getComponentType();
                if (clazz != String.class) {
                    continue;
                }
                int modifiers = f.getModifiers();
                if (!Modifier.isStatic(modifiers)) {
                    continue;
                }
                
                nameList.add(f.getName());
                polList.add((String[])f.get(null));
            }
        }
        catch (IllegalAccessException ignore) { /* never occurs */ }
    }
    
    
    
    public static boolean hasRealCoefs(String[] pol) {
        return pol != null && (pol[0].equals("RR") || pol[0].equals("CR"));
    }
}
