package algorithms.statistics;

import algorithms.SubsetChooser;
import algorithms.matrix.MatrixUtil;
import algorithms.misc.Distances;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import algorithms.util.PairInt;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
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
     * @return vector of the probabilities
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
    
    /**
    <pre>
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
    
        L_1(g) = integral_{x=0_to_2} ( | g_N(x) - f_N(x)| * dx )

           where f_N(x) is hypersphere chord length distribution == the pdf,
               and
           g_N is intra-distance distribution of the input point distribution 
               using model (2) above.
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
    </pre>
     * @param x data points in format [nSamples][nDimensions] that are to be
     * tested as uniformly distributed on an nDimensions hypersphere
     * (note, there may be some confusion in definitions as the literature using 
     * von Mises–Fisher distributions for the hypersphere have a multiplier
     * that includes a dimension that is then multiplied by an n-sphere of dimension n-1).
     * @param m the number of points to choose from x.length in calculating L1
     * @return the calculated L1 statistic from eqn (26) of the paper.
    */
    public static double calcL1UniformityStatistic(double[][] x, int m) throws NoSuchAlgorithmException {
        
        if (m > x.length) {
            throw new IllegalArgumentException("m must be less tna x.length");
        }
        
        // g_N:
        // -- from the n-sphere points to be tested, choose m of the points.
        // -- the M(M − 1)/2 point pairs will be constructed 
        //    and each distance is the euclidean distance 
        //    ||pi − pj||,i,j = 1,2,...,M,i ̸= j.
        // looks like after the pairwise distances are calculated, they 
        //    need to be ordered (from 0 to 2*r with r=1)
        //    before used in the L_1(g) calculation.
        
        // f_N(x)
        //   x is the chord length d, so the integration is from d=0 to d=2 (probably 2*r).
        //double[] pdf(double[] d, double r, int n)
        
        /*
        L_1(g) = integral_{x=0_to_2} ( | g_N(x)- f_N(x)| * dx )

           where f_N(x) is hypersphere chord length distribution == the pdf,
               and
           g_N is intra-distance distribution of the input point distribution 
               using model (2) above.
           x is the chord length d, so the integration is from d=0 to d=2 (probably 2*r).
        */
        //L_1(g) = integral_{x=0_to_2} ( | g_N(x)- f_N(x)| * dx )
        
        int[] xMIdx = chooseM(m, x.length);

        double[] d1 = chooseMCalcPairwiseDistances(x, xMIdx);
        
        int[] indexes1 = MiscSorter.mergeSortIncreasing(d1);
        
        //double[] d2 = calculateDistancesFromOrigin(x, xMIdx);
        
        //int[] indexes2 = MiscSorter.mergeSortIncreasing(d2);
        
        double[] p2 = pdf(d1, 1, x[0].length);
        
        // the details needed to create the integral seem to be these:
        //
        // calculate p1 for each item in d2 by:
        //   -- create a histogram from d1 and normalize it so that the
        //      histogram counts sum to 1
        //   -- the normalized count for the bin holding each item of d2
        //      is the value of g_N(x) to be subtracted from the corresponding item in p2
     
        double[] minMaxD1 = MiscMath0.getMinMax(d1);
        //double[] minMaxD2 = MiscMath0.getMinMax(d2);
        double minD = minMaxD1[0];
        double maxD = minMaxD1[1];
     
        HistogramHolder h1 = Histogram.calculateScottsHistogram(d1);
        float[] p1 = h1.getYHistFloat();
        float[] p1x = h1.getXHist();
        float[] minMaxP1 = MiscMath0.getMinMax(p1);
        float p1Norm = 1.f/minMaxP1[1];
        for (int i = 0; i < p1.length; ++i) {
            p1[i] *= p1Norm;
        }
        
        try {
            String str = h1.plotHistogram("p1", "p1_hist");
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(d1, p2, null, null, "p2");
            String str2 = plotter.writeFile("p2_hist");
        } catch (IOException ex) {
            Logger.getLogger(HypersphereChordLength.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        float binWidth1 = p1x[1] - p1x[0];
        float p1Min = p1x[0] - (binWidth1/2.f);
        float p1Max = p1x[p1x.length - 1] + (binWidth1/2.f);
        
        double l1Sum = 0;
        double d;
        double diff;
        int bin1;
        for (int i = 0; i < d1.length; ++i) {
            d = d1[i];
            
            if (d < 0 || d > 2) {
                int z = 0;
            }
            
            bin1 = (int) ((d - p1Min)/binWidth1);
            if ((bin1 < 0) || (bin1 >= p1x.length)) {
                throw new IllegalStateException("bin1 is outside of bounds of histogram");
            }
            diff = p1[bin1] - p2[i];
            l1Sum += (diff*diff);
        }
        l1Sum = Math.sqrt(l1Sum);
        
        /*
            https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm

            chisq stat for deg of freedom

            For an upper-tail one-sided test, find the column corresponding to 
                1-α in the table containing upper-tail critical 
                 and reject the null hypothesis if the test statistic is 
                 greater than the tabled value.

            Upper-tail critical values of chi-square distribution with ν degrees of freedom
                            Probability less than the critical value
               ν           0.90      0.95     0.975      0.99     0.999

               1          2.706     3.841     5.024     6.635    10.828
               2          4.605     5.991     7.378     9.210    13.816
               3          6.251     7.815     9.348    11.345    16.266
               4          7.779     9.488    11.143    13.277    18.467
               5          9.236    11.070    12.833    15.086    20.515
               6         10.645    12.592    14.449    16.812    22.458
               7         12.017    14.067    16.013    18.475    24.322
               8         13.362    15.507    17.535    20.090    26.125
               9         14.684    16.919    19.023    21.666    27.877
              10         15.987    18.307    20.483    23.209    29.588
            */
        
        return l1Sum;
    }
    
    /**
     * randomly choose m numbers from 0 to n-1, inclusive.
     * @param m
     * @param n
     * @return
     * @throws NoSuchAlgorithmException 
     */
    static int[] chooseM(int m, int n) throws NoSuchAlgorithmException {
                        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        TIntSet s = new TIntHashSet();
        int[] mr = new int[m];
        int i;
        int rd;
        for (i = 0; i < m; ++i) {
            rd = rand.nextInt(n);
            while (s.contains(rd)) {
                rd = rand.nextInt(n);
            }
            mr[i] = rd;
            s.add(mr[i]);
        }
        return mr;
    }

    static double[] chooseMCalcPairwiseDistances(double[][] x, int[] idx) throws NoSuchAlgorithmException {
        
        if (idx.length > x.length) {
            throw new IllegalArgumentException("m cannot be larger than x.length");
        }
        
        int m = idx.length;
                
        int np = (int) MiscMath0.computeNDivKTimesNMinusKExact(m, 2);
        
        double[] d = new double[np];
                
        SubsetChooser chooser = new SubsetChooser(m, 2);
        
        int[] selectedIndexes = new int[2];
        
        int c = 0;
        while (chooser.getNextSubset(selectedIndexes) != -1) {
            int idx0 = idx[selectedIndexes[0]];
            int idx1 = idx[selectedIndexes[1]];
            d[c] = Math.sqrt(Distances.calcEuclideanSquared(x[idx0], x[idx1]));
            c++;
        }
        
        return d;
    }

    /**
     * calculate the distance of each point from the origin.
     * assumes points are already w.r.t origin as they are being tested for
     * uniformity on the unit sphere
     * @param x points on the unit sphere
     * @param idx subset of indexes of x to calculate the distances for.
     * @return 
     */
    private static double[] calculateDistancesFromOrigin(double[][] x, int[] idx) {
                
        int n = idx.length;
        int nd = x[0].length;

        // 
        double[] d = new double[n];
        double s;   
        int i, j;
        for (i = 0; i < n; ++i) {
            s = 0;
            for (j = 0; j < nd; ++j) {
                s += (x[idx[i]][j]*x[idx[i]][j]);
            }
            d[i] = Math.sqrt(s);
        }
        
        return d;
    }
    
    
}
