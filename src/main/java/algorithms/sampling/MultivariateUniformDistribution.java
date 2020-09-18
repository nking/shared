package algorithms.sampling;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

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
     * generate random uniform points on an n-dimensional unit standard hypersphere.
     * NOTE: to generate multiple samples, invoke the overloaded method multiple times.
     * <pre>
     * paper containing references:
     * "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     * https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     * </pre>
     * @param d number of dimensions.
     * @return vector of size d+1 
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandardOnNSphere(int d) 
        throws NoSuchAlgorithmException {
        
         SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        return generateUnitStandard(d, rand, true);
    }
    
    /**
     * generate random uniform points in an n-dimensional unit standard ball.
     * NOTE: to generate multiple samples, invoke the overloaded method multiple times.
     * <pre>
     * paper containing references:
     * "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     * https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     * </pre>
     * @param d number of dimensions.
     * @return vector of size d+1 
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandardInNBall(int d) 
        throws NoSuchAlgorithmException {
        
         SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        return generateUnitStandard(d, rand, false);
    }
    
    /**
     * generate random uniform points on an n-dimensional unit standard hypersphere
     * (for onSurface=true) or in an n-dimensional unit standard ball.
     * NOTE: to generate multiple samples, invoke the overloaded method multiple times.
     * <pre>
     * paper containing references:
     * "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     * https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     * </pre>
     * @param d number of dimensions.
     * @param rand
     * @param onSurface
     * @return vector of size d+1 for onSurface == true, else vector of size d
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandard(int d, SecureRandom rand, boolean onSurface) 
        throws NoSuchAlgorithmException {
        
        if (onSurface) {
            d++;
        }
        
        /*
        also see the SAS blog for similar algorithm of d dimensions and N samples:
        https://blogs.sas.com/content/iml/2016/04/06/generate-points-uniformly-in-ball.html
        radius = 2;                          // radius of circle 
        Y = randfun(N // d, "Normal");       // Y ~ MVN(0, I(d)) 
        u = randfun(N, "Uniform");           // U ~ U(0,1)       
        r = radius * u##(1/d);               // r proportional to d_th root of U 
        X = r # Y / sqrt(Y[,##]);            // Y[,##] is sum of squares for each row 
        // X contains N random uniformly distributed points in the d-bal
        */
        
        // note: Voelker et al. 2017 demonstrate with a python library called Nengo
       
        // generate d dimensions individ random standard normal N(0,1)  
        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, d);
        
        // then normalization of each is the sqrt(sum of squares of all)
        double norm = 0;
        for (double v : u) {
            norm += (v*v);
        }
        norm = 1./Math.sqrt(norm);
        
        if (!onSurface) {
            // to avoid bunching at the center of hypersphere
            double b = Math.pow(rand.nextDouble(), 1./d);
        
            norm *= b;
        }
        
        for (int i = 0; i < u.length; ++i) {
            u[i] *= norm;
        }
                
        return u;
    }
    
}
