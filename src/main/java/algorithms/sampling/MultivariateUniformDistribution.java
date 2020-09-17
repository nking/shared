package algorithms.sampling;

import java.security.NoSuchAlgorithmException;

/**
 * create random vectors uniformly distributed on the surface of the unit n-1 
 * sphere.
 * @author nichole
 */
public class MultivariateUniformDistribution {
    
    /*
     U(0, 1) is 
            standard uniform (rectangular) distribution with a=0 and b=1
            <pre>
            PDF = 1/(b-1)
            CDF = (x-a)/(b-a) = x/1 for unit standard
            inverseCDF = a + p*(b-a)
            mean = (a+b)/2 = median
            variance = ((b-a)^2)/12
            information content = log_2(b)
            </pre>
    */
    
    /**
     * random uniform points in an n-dimensional uniform unit standard ball.
     * @param nDimensions
     * @return
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandardSample(int nDimensions) throws NoSuchAlgorithmException {
        
        //also see https://www.sciencedirect.com/science/article/pii/S0047259X10001211
        /*Ondecompositionalalgorithmsforuniformsamplingfromn-spheresandn-balls
        Harman, and Lacko
        Journal of Multivariate Analysis 101 (2010) 2297–2304
        */
        
        // generate nDimensions individ random standard normal N(0,1)        
        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(nDimensions);
        
        // then normalization of each is the sqrt(sum of squares of all)
        double sqSumSq = 0;
        for (double v : u) {
            sqSumSq += (v*v);
        }
        
        for (int i = 0; i < u.length; ++i) {
            u[i] /= sqSumSq;
        }
        
        /*
        // the sas blog uses a radius factor too:
        // r = radius * u##(1/d);
        //   which means r = radius * U(0,1) of dimension nDimensions to the power (1/d)
        //      so a vector u_n of (1)^(1/d) is all 1's(??), then multiplying by radius
        //      which is 1 here
       
        https://blogs.sas.com/content/iml/2016/04/06/generate-points-uniformly-in-ball.html
        NOTE: https://support.sas.com/rnd/app/iml/matrix.html
        multiplication operator, elementwise:   #
        multiplication operator, matrix:   *
        power operator, elementwise:   ##
        power operator, matrix:   **
        
        proc iml;
        call randseed(12345);
        d = 3;                               // dimension = number of variables 
        N = 1000;                            // sample size = number of obs    
        radius = 2;                          // radius of circle 
        Y = randfun(N // d, "Normal");       // Y ~ MVN(0, I(d)) 
        u = randfun(N, "Uniform");           // U ~ U(0,1)       
        r = radius * u##(1/d);               // r proportional to d_th root of U 
        X = r # Y / sqrt(Y[,##]);            // Y[,##] is sum of squares for each row 
        // X contains N random uniformly distributed points in the d-bal
        */
        
        return u;
    }
    
    /*
    from https://rdrr.io/cran/NonNorMvtDist/man/MvtUniform.html
    
    */
}
