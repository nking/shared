package algorithms.sampling;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Arrays;

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
    
    /*
    https://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
    In general, the volume of the n-ball in n-dimensional Euclidean space, and 
    the surface area of the n-sphere in (n + 1)-dimensional Euclidean space, 
    of radius R, are proportional to the nth power of the radius, R (with 
    different constants of proportionality that vary with n). 
    
    let R = radius.
    the volume of the n-ball:
        V_n(R) = V_n*R^n
        for unit-radius: V_n = V_n(1)
    the surface area of the n-sphere
        S_n(R) = S_n*R^n 
        for unit-radius: S_n = S_n(1)
    
    Hausdorff measure is an outer measure, it is the number of points in a set.
    
    n general, the volume, in n-dimensional Euclidean space, of the unit n-ball, 
    is given by:
        V_n = (pi^(n/2))/(gamma(0.5*n + 1) = (pi^(n/2))/((n/2)!)
    
        S_{n-1} = (2*pi^(n/2))/(gamma(0.5*n)
    
    S_{n+1} = 2*pi*V_n
    
    coordinate system in an n-dimensional Euclidean space have 
        a radial coordinate r, 
        and n − 1 angular coordinates φ1, φ2, ... φn−1, 
            where the angles φ1, φ2, ... φn−2 range over [0,π] radians 
            (or over [0,180] degrees) 
            and φn−1 ranges over [0,2π) radians (or over [0,360) degrees). 
        If xi are the Cartesian coordinates, then we may compute 
        x1, ... xn from r, φ1, ... φn−1 with:
            x_1 = r*cos(φ1)
            x_2 = r*sin(φ1)*cos(φ2)
            x_3 = r*sin(φ1)*cos(φ2)
            ...
            x_{n-1} = r*sin(φ1)...sin(φ_{n-2})*cos(φ_{n-1})
            x_{n}   = r*sin(φ1)...sin(φ_{n-2})*sin(φ_{n-1})
    
        r = sqrt( x_{n}^2 + x_{n-1}^2 + ... x_{2}^2 + x_{1}^2 )
    
        φ1 = arccot( (x_1)^2 / sqrt( x_{n}^2 + x_{n-1}^2 + ... x_{2}^2 ))
           = arccos( (x_1)^2 / sqrt( x_{n}^2 + x_{n-1}^2 + ... x_{2}^2 + x_{1}^2 ))
    
        φ2 = arccot( (x_2)^2 / sqrt( x_{n}^2 + x_{n-1}^2 + ... x_{3}^2 ))
           = arccos( (x_2)^2 / sqrt( x_{n}^2 + x_{n-1}^2 + ... x_{3}^2 + x_{2}^2 ))
        ...
        φ_{n-2} = arccot( (x_{n-2})^2 / sqrt( x_{n}^2 + x_{n-1}^2))
                = arccos( (x_{n-2})^2 / sqrt( (x_{n}^2 + x_{n-1}^2 + x_{n-2}^2 ))
        φ_{n-1} = 2*arccot( (x_{n-1}^2 + sqrt( x_{n}^2 + x_{n-1}^2)) / x_n )
                = for x_n >= 0 : arccos( (x_{n-1})^2 / sqrt( (x_{n}^2 + x_{n-1}^2 ))
                = for x_n < 0 : 2*pi - arccos( (x_{n-1})^2 / sqrt( (x_{n}^2 + x_{n-1}^2 ))
    
       For the above, if there is a range where x_{k+1}, ... xn are zero 
          and there is a x_k != 0, then φk = 0 when x_k>0 and φk = pi when x_k>0.
    
       For special cases where the inverse transform is not unique; 
          φk for any k will be ambiguous whenever all of xk, xk+1, ... xn are zero; 
          in this case φk may be chosen to be zero.
    
       If ever need orthogonal:
          Polyspherical coordinates...points in Rn may be expressed by taking 
          the ray starting at the origin and passing through z ∈ Real^{n − 1}, 
           rotating it towards the first basis vector by θ, and traveling a 
           distance r along the ray. (Repeating this decomposition eventually 
           leads to the standard spherical coordinate system.)
          ... The possible polyspherical coordinate systems correspond to 
          binary trees with n leaves. Each non-leaf node in the tree corresponds 
          to a splitting and determines an angular coordinate. For instance, 
          the root of the tree represents Real^n, and its immediate children 
          represent the first splitting into Real^p and Real^q.  
          Leaf nodes correspond to Cartesian coordinates for S^{n−1}.
          ...
          special orthogonal group: A splitting Real^n = Real^p × Real^q 
              determines a subgroup (subgroup p times subgroup q is in subgroup n).
              There are special factors relating the VOLUME measure on Real^n 
              and the AREA measure on S^{n − 1} as products. 
              There is one factor for each angle, and the factor is determined 
              by the tree.
              If p == 1 and q > 1, the factor has a denominator that is 
              uses the beta function 0.5*Beta(p/2, q/2).
              can invet it to get theta which is defined in terms of
              polyspherical coordinates.
          see https://en.wikipedia.org/wiki/N-sphere#Polyspherical_coordinates
    */
    
    /*
    interpoint spacings:
    
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
     * @param rand
     * @return vector of size d
     * @throws NoSuchAlgorithmException 
     */
    public static double[] generateUnitStandardOnNSphere(int d, SecureRandom rand) 
        throws NoSuchAlgorithmException {
        
        return _generateUnitStandard(d+1, rand, true);
    }
    
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
        
        return _generateUnitStandard(d+1, rand, true);
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
    public static double[] generateUnitStandardInNBall(int d, SecureRandom rand) 
        throws NoSuchAlgorithmException {
        
        return _generateUnitStandard(d, rand, false);
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
        
        return _generateUnitStandard(d, rand, false);
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
     * @param onSurface when true the method generates points on the 
     * unit hyper-sphere that and assumes that d has already been incremented to d+1,
     * else when false the method generates points within the unit hyper-sphere
     * that is within the d-ball.
     * @return vector of length d-1 if onSurface is true, else returns a vector
     * of length d.
     * @throws NoSuchAlgorithmException 
     */
    public static double[] _generateUnitStandard(int d, SecureRandom rand, boolean onSurface) 
        throws NoSuchAlgorithmException {
        
        /*
        From 2010 Harman and Lacko:
        "On decomposition alalgorithms for uniform sampling from n-spheres and n-balls"
         Journal o fMultivariate Analysis 101 (2010) 2297–2304.
         the multivariate normal distribution with independent standardized 
         components is radially symmetric, i.e., it is invariant under orthogonal
         rotations.  Therefore, if Y=N_n(0_n, I_n), then S_n = Y/‖Y‖ 
         has the uniform distribution on the unit n-sphere, which was noted 
         already by Muller[23].  Moreover, multiplying S_n by U^(1/n), where
         U has the uniform distribution on the unit interval (0,1), we obtain
         the uniform distribution on the unit n-ball; see, e.g., Section 3.29 in [12]."
        */
        
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
        //     that is freely available, but also proprietary.  see license in
        //     comments at end of this file.
       
        double[] u = new double[d];
                
        u = UnivariateNormalDistribution.randomSampleOfUnitStandard(rand, d);
        
        // then normalization of each is the sqrt(sum of squares of all)
        double norm = 0;
        for (double v : u) {
            norm += (v * v);
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
        
        if (onSurface) {
            return Arrays.copyOf(u, u.length - 1);
        }
                
        return u;
    }
    
    // TODO: if needed, add method to Uniformly sampling coordinates from the n-sphere
    //  from Section 2.3 of Voelker, Gosman, Stewart 2017
    // To sample coordinates from the unit n-sphere (i.e., uniform points from 
    //  the sphere projected onto an arbitrary unit vector) we could simply 
    //  modify §2.1 to return only a single element 
    //  (which is _generateUnitStandard with onSurface = true
    
    /*
    on subject of Uniformly sampling coordinates from the n-sphere
    
    
    https://github.com/nengo/nengo-extras/issues/34
    by user https://github.com/tcstewar
    
    sqrt(1-betaincinv((d+1)/2, 1/2, x+1))
    (where x is the original intercept, uniformly distributed between -1 and 1, and d is the dimensionality)
    
    def transform(self, x):
        sign = 1
        if x > 0:
            x = -x
            sign = -1
        return sign * np.sqrt(1-scipy.special.betaincinv((self.dimensions+1)/2.0, 0.5, x+1))
    
    def sample(self, n, d=None, rng=np.random):
        // where self.base=nengo.dists.Uniform(-1, 1)
        s = self.base.sample(n=n, d=d, rng=rng)
        for i in range(len(s)):
            s[i] = self.transform(s[i])
        return s
    
    -------------
    when have an implementation of betaincinv,
    could pursue this for cases of very large number of dimensions,
    else, just sampling 1 item from the uniform sphere sampling above is more
    efficient.
    
    Nengo license:
    https://github.com/nengo/nengo/blob/11b809e486c8b8591320ff15ac991157b0b7716f/LICENSE.rst
    Copyright (c) 2013-2020 Applied Brain Research

Nengo is made available under a proprietary license that permits using, copying, 
    sharing, and making derivative works from Nengo and its source code for any 
    non-commercial purpose, as long as the above copyright notice and this 
    permission notice are included in all copies or substantial portions of 
    the software.
If you would like to use Nengo commercially, licenses can be purchased 
    from Applied Brain Research. Please contact info@appliedbrainresearch.com 
    for more information.
    
These are details from Nengo issue above and source code below:
    
    nengo.dists.CosineSimilarity(dist.dimensions + 2).ppf(np.linspace(0, 1, n))
    
    looking at
    https://github.com/nengo/nengo/blob/master/nengo/dists.py
    and numpydocs
    
    numpy.linspace(start, stop, num=50, endpoint=True, retstep=False, dtype=None, axis=0)[source]
        Return evenly spaced numbers over a specified interval.
        Returns num evenly spaced samples, calculated over the interval [start, stop].
        The endpoint of the interval can optionally be excluded.
    
    => base/parent class SqrtBeta receives dist.dimensions + 2
       which assigns m=1=subdimensions, and n=(dist.dimensions + 2-m)=dimensions
    => CosineSimilarity.ppf receives
       the results from linspace: n numbers between 0 and 1
       (where n is not the same variable as SqrtBeta local n, presumably).
       --> y = n numbers from 0 to 1.
       --> x = super().ppf(abs(y * 2 - 1))
           states that SqrtBeta.ppf receives y*2 - 1
           SqrtBeta local y is Cumulative probabilities in [0, 1].
           sb.y is y*2 - 1 in sb.ppf
           sq_x = betaincinv(self.m / 2.0, self.n / 2.0, sb.y)
           return np.sqrt(sq_x) which is Evaluation points ``x`` in [0, 1] such that ``P(X <= x) = y``.
              
    cython_special.betaincinv?
    
    class CosineSimilarity(SubvectorLength):
    def ppf(self, y):
        x = super(CosineSimilarity, self).ppf(abs(y*2 - 1))
        return np.where(y > 0.5, x, -x)
    def ppf(self, y):
        x = super().ppf(abs(y * 2 - 1))
        return np.where(y > 0.5, x, -x)
    
    
    class SubvectorLength(SqrtBeta):
    """Distribution of the length of a subvectors of a unit vector.
    Parameters
    ----------
    dimensions : int
        Dimensionality of the complete unit vector.
    subdimensions : int, optional
        Dimensionality of the subvector.
    """
    def __init__(self, dimensions, subdimensions=1):
        super().__init__(dimensions - subdimensions, subdimensions)
    @property
    def dimensions(self):
        return self.n + self.m
    @property
    def subdimensions(self):
        return self.m

    
    class SqrtBeta(Distribution):
    def ppf(self, y):
        """Percent point function (inverse cumulative distribution).
        .. note:: Requires SciPy.
        Parameters
        ----------
        y : array_like
            Cumulative probabilities in [0, 1].
        Returns
        -------
        ppf : array_like
            Evaluation points ``x`` in [0, 1] such that ``P(X <= x) = y``.
        """
        from scipy.special import betaincinv  # pylint: disable=import-outside-toplevel

        sq_x = betaincinv(self.m / 2.0, self.n / 2.0, y)
        return np.sqrt(sq_x)
    
    x would be 
    
    */
   
}
