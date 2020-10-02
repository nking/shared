package algorithms.statistics;

import thirdparty.smile.math.special.Beta;

/**
 *  functions useful for testing uniformity of points on the n-sphere
 * from the paper:
 * <pre>
 * Measuring spatial uniformity with the hypersphere chord length distribution
           2020, Panagiotis Sidiropoulos
        https://arxiv.org/pdf/2004.05692.pdf
        arXiv:2004.05692v1 [cs.CG] 12 Apr 2020
 * </pre>
 * @author nichole
 */
public class HypersphereChordLength {
    
    /**
     * M points selected uniformly and independently from the surface of a
        n-dimensional hypersphere of radius R.
     * @param d pairwise euclidean distances
     * @param r radius of the n-dimensional hypersphere.  e.g.use r=1 for unit standard.
     * @param n the number of dimensions of the hypersphere
     * @return 
     */
    public static double[] pdf(double[] d, double r, int n) {
        
        /*
        term1 = r^2 * beta( (n-1)/2, 1/2 )
        
        term2 = ( (d^2/r^2) - (d^4/(4*r^4)) )
        
        f_n(d) = (d/term1) * ( term2^((n-3)/2) )
        */
        
        double[] pdf = new double[d.length];
        
        double rsq = r*r;
        double rq = rsq * rsq;
        double term1 = rsq * Beta.beta((n-1)/2., 0.5);
        double term2, dsq, dq;
        
        for (int i = 0; i < d.length; ++i) {
            dsq = d[i]*d[i];
            dq = dsq * dsq;
            term2 = (dsq/rsq) - (0.25*dq/rq);
            pdf[i] = (d[i]/term1) * Math.pow(term2, (n-3)/2.);
        }
        
        return pdf;
    }
    
    /**
     * M points selected uniformly and independently from the surface of a
        n-dimensional hypersphere of radius R.
     * @param d pairwise euclidean distances
     * @param r radius of the n-dimensional hypersphere.  e.g.use r=1 for unit standard.
     * @param n the number of dimensions of the hypersphere
     * @return 
     */
    public static double[] cdf(double[] d, double r, int n) {
        
        /*
        x = ( (d^2/r^2) - (d^4/(4*r^4)) )
        
        a = (n-1)/2
        b = 1/2
        
        I is regularized incomplete beta function
        
        for d < r*sqrt(2):
            P(D<=d) = F_n(d) = (1/2) * I(a, b, x)
        for d >= r*sqrt(2)
            P(D<=d) = F_n(d) = 1. - (1/2) * I(a, b, x)
        */
        
        double[] cdf = new double[d.length];
        
        double rsqrt2 = r * Math.sqrt(2.);
        double rsq = r*r;
        double rq = rsq * rsq;
        double dsq, dq, x;
        
        for (int i = 0; i < d.length; ++i) {
            dsq = d[i]*d[i];
            dq = dsq * dsq;
            x = (dsq/rsq) - (0.25*dq/rq);
            
            cdf[i] = 0.5 * Beta.regularizedIncompleteBetaFunction((n-1)/2., 0.5, x);
            
            if (d[i] >= rsqrt2) {
                cdf[i] = 1. - cdf[i];
            }
        }
        
        return cdf;
    }
    
    /**
     *
     * @param k the k-th moment
     * @param r radius of the n-dimensional hypersphere.  e.g.use r=1 for unit standard.
     * @param n the number of dimensions of the hypersphere
     * @return 
     */
    public static double momentAboutOrigin(int k, double r, int n) {
        
        /*
        term1 = 2^(k+n-2)
        term2 = B((n-1)/2, 1/2)
        term3 = B( (k+n-1)/2, (n-1)/2)
        term4 = r^k
        
        E(D^k) = (term1/term2) * term3 * term4
        */
        
        // could use 1 << (k + n - 2) instead for power of 2:
        double term1 = Math.pow(2, k + n - 2);
        double term2 = Beta.beta((n-1)/2., 0.5);
        double term3 = Beta.beta((k+n-1)/2., (n-1)/2.);
        double term4 = Math.pow(r, k);
        
        return (term1/term2) * term3 * term4;
    }
    
    /**
     *
     * @param r radius of the n-dimensional hypersphere.  e.g.use r=1 for unit standard.
     * @param n the number of dimensions of the hypersphere
     * @return 
     */
    public static double meanOfChordLengthDistribution(double r, int n) {
        
        /*
        term1 = (gamma(n/2))^2
        term2 = gamma(n-(1/2))
        term3 = sqrt(pi)
        term4 = r*2^(n-1)
        
        mu = term4*(term1/(term2 * term3))
        */
        
        double term1 = Gamma.lanczosGamma9(n/2.);
        term1 *= term1;
        
        double term2 = Gamma.lanczosGamma9(n - (1./2.));
        
        double term3 = Math.sqrt(Math.PI);
        
        // could use 1 << (n-1) instead for power of 2:
        double term4 = r * Math.pow(2, n - 1.);
        
        return term4*(term1/(term2 * term3));
    }
    
    /**
     *
     * @param r radius of the n-dimensional hypersphere.  e.g.use r=1 for unit standard.
     * @param n the number of dimensions of the hypersphere
     * @return 
     */
    public static double varianceOfChordLengthDistribution(double r, int n) {
        
        /*
        term1 = (gamma(n/2))^4
        term2 = (gamma(n-(1/2)))^2
        term3 = pi
        term4 = 2^(2*n-2)
        
        sigma^2 = (r^2)*(2 - ((term1*term4)/(term3*term2)))
        */
        
        double term1 = Gamma.lanczosGamma9(n/2.);
        term1 = Math.pow(term1, 4.);
        
        double term2 = Gamma.lanczosGamma9(n - (1./2.));
        term2 *= term2;
        
        double term3 = Math.PI;
        
        // could use 1 << (2*n-2) instead for power of 2:
        double term4 = Math.pow(2, 2.*n-2);
        
        return (r*r)*(2. - ((term1*term4)/(term3*term2)));
    }
    
    /*
    IV. HYPERSPHERE CHORD LENGTH DISTRIBUTION AS A UNIFORMITY MEASURE
       As already mentioned, some of the most interesting properties of the 
    hypersphere chord length distribution arise from the fact that this is the 
    limit distribution of the distances of uniformly selected hypersphere points. 
    To summarise, the hypersphere chord length distribution is the limit 
    distribution of 3 (related but distinct) distributions:
        1) M point pairs:
           The distance distribution of 
           M point-pairs ||pi − p′i ||, i = 1, 2, ..., M , 
           if the 2M relevant points pi and p′i are independently selected 
           from a uniform random distribution.
        2) M(M − 1)/2 point pairs:
           The intra-distance distribution of a set of 
           M points pi,i=1,2,...,M,
           if the M relevant points are independently selected 
           from a uniform random distribution, 
           before the M(M − 1)/2 pairwise distances 
           ||pi − pj||,i,j = 1,2,...,M,i ̸= j are estimated.
        3) M - 1 point pairs:
           The distance distribution of a set of 
           M − 1 points pi,i = 1,2,...,M − 1 
           from a fixed point p0 if both pi and p0 are selected 
           from a uniform random distribution 
           before the M − 1 pairwise distances 
           ||pi − p0 ||, i = 1, 2, ..., M − 1 are estimated.
    The second and the third distribution allow the hypersphere chord length 
    distribution to be used as an uniformity measure, as is described in the 
    current section.
    
    to quantify the ”uniformity” of an input point distribution on a N-sphere, 
    the L1 distance is used.  
    
        L_1(g) = integral_{x=0_to_2} ( | g_N(x)- f_N(x)| * dx )

           where f_N(x) is hypersphere chord length distribution, the pdf.
               and
           g_N is intra-distance distribution of the input point distribution 
               using models above from (1), (2), or (3).
           x is the chord length d, so the integration is from d=0 to d=2 (probably 2*r).
    
    Initially, uniform pointsets of size M′ (where M′ is much smaller than M) 
    are generated on the N-sphere and L1 values are sorted, before acquiring 
    the α%-largest L1 value and finally extrapolating for pointsets of size M. 
    This value is the threshold with which the input distance distribution 
    L1(g) is compared to determine whether it is uniform or not.
    
    Elaborating on this idea, based on the computational cost of iteratively 
    estimating pairwise distances, L1 can be used to qualitatively assess 
    whether an N-dimensional point sample S (consisting of M points, and having 
    an intra-distribution g) originates from a uniform N-sphere (or 
    N-hemisphere) distribution following on of the three following approaches:
        • If S dimension N and sample size M imply a non-prohibitive 
          computational cost, then Q uniformly distributed point sets of size M 
          and dimension N are randomly generated and L1 is estimated for all of 
          them (as well as for S). If the α%-largest L1 value of Q is smaller 
          than L1(g) then S can be declared as non-uniform with confidence 
          (100 − α)%.
    */
    
    
}
