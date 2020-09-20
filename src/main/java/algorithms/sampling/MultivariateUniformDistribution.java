package algorithms.sampling;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 * create random vectors uniformly distributed on or within unit standard n-dimensional
  hyper-sphere.
  
 from http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
 ‘uniform random sampling‘ is also known as ‘(continuous) simple random sampling‘, 
 is also known as ‘uniform sampling‘.
 For discrete probabilities, this means that all possible elements of 
 have an equal probability of being selected. 
 For continuous probabilities, this means that the likelihood of an element 
 falling in any subinterval is directly proportional to the length of the 
 subinterval.
 
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
     * generate random uniform points on an n-dimensional unit standard hyper-sphere.
     * NOTE: the result is one sample of all dimensions.
     <pre>
     paper containing references:
     "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     </pre>
     <pre>
     from http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
     for both the n-sphere and the n-ball, the n signifies how many degrees 
     of freedom it has:
       A unit n-dimensional sphere is defined such that:
          S^n = {x is member of the real set of dimension (n+1) : |x|=1}
       A unit n-dimensional ball is defined such that:
          B^n = {x is member of the real set of dimension n : |x|.lte.1}
       So the perimeter of a circle is a 1-sphere, 
          and the interior of the circle (a disk) is a 2-ball.
     </pre>
     * Muller / Marsaglia (‘Normalized Gaussians’).
     * @param d number of dimensions.
     * @return vector of size d
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
     * NOTE: the result is one sample of all dimensions.
     <pre>
     paper containing references:
     "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     </pre>
     <pre>
     from http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
     for both the n-sphere and the n-ball, the n signifies how many degrees 
     of freedom it has:
       A unit n-dimensional sphere is defined such that:
          S^n = {x is member of the real set of dimension (n+1) : |x|=1}
       A unit n-dimensional ball is defined such that:
          B^n = {x is member of the real set of dimension n : |x|.lte.1}
       So the perimeter of a circle is a 1-sphere, 
          and the interior of the circle (a disk) is a 2-ball.
     </pre>
     * Muller / Marsaglia (‘Normalized Gaussians’).
     * @param d number of dimensions.
     * @return vector of size d
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
     NOTE: the result is one sample of all dimensions.
     Muller / Marsaglia (‘Normalized Gaussians’).
     <pre>
     paper containing references:
     "Efficiently sampling vectors and coordinates from the n-sphere and n-ball"
    Aaron R. Voelker, Jan Gosmann, Terrence C. Stewart Centre for Theoretical Neuroscience – Technical Report
    January 4, 2017
     https://www.researchgate.net/profile/Aaron_Voelker/publication/312056739_Efficiently_sampling_vectors_and_coordinates_from_the_n-sphere_and_n-ball/links/586d3e1208ae6eb871bce1d6/Efficiently-sampling-vectors-and-coordinates-from-the-n-sphere-and-n-ball.pdf
     </pre>
     <pre>
     from http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
     for both the n-sphere and the n-ball, the n signifies how many degrees 
     of freedom it has:
       A unit n-dimensional sphere is defined such that:
          S^n = {x is member of the real set of dimension (n+1) : |x|=1}
       A unit n-dimensional ball is defined such that:
          B^n = {x is member of the real set of dimension n : |x|.lte.1}
       So the perimeter of a circle is a 1-sphere, 
          and the interior of the circle (a disk) is a 2-ball.
     </pre>
     * Muller / Marsaglia (‘Normalized Gaussians’).
     * @param d number of dimensions.
     * @param rand
     * @param onSurface
     * @return vector of size d
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandard(int d, SecureRandom rand, boolean onSurface) 
        throws NoSuchAlgorithmException {
        
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
        //     that is freely available, but also proprietary.
       
        // generate d dimensions individ random standard normal N(0,1)  
        double[] u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, d);
        /*double[] u = new double[d];
        for (int i = 0; i < d; ++i) {
            u[i] = rand.nextDouble();
        }*/
        
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
