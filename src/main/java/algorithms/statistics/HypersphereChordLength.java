package algorithms.statistics;

import algorithms.SubsetChooser;
import algorithms.misc.Distances;
import algorithms.misc.Histogram;
import algorithms.misc.HistogramHolder;
import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.security.SecureRandom;
import java.util.Arrays;
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
    
    public static enum POINT_DISTRIBUTION_TYPE {
        INTRA_DISTANCE_2, POINT_DISTANCE_3
    }
    
    public static class NonUniformityStats {
        double oneMinusAlpha;
        double oneMinusAlphaCritVal;
        double l1MaxSphere;
        double l1MaxX;
        boolean isConsistentWithNonUniform;
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("isConsistentWithNonUniform=").append(isConsistentWithNonUniform);
            sb.append(String.format("\n1-alpha=%.5f, crit. val.=%.4e \nL1(sphere)=%.4f \nL1(X)=%.4f", 
                oneMinusAlpha, oneMinusAlphaCritVal, l1MaxSphere, l1MaxX));
            return sb.toString();
        }
    }
    
    /**
    test to identify non-uniform spatial distribution, assessing 
    whether the points span the hypersphere 
    in a way that is compatible with the uniform distribution.
    <pre>
    see notes belowx in method comments for calcL1UniformityStatistic
    </pre>
    * TODO: determine the number of iterations to calculate L1's in some manner.
    * It's currently set to a value of 100 internally.
    */
    public static NonUniformityStats calcConfidenceOfNonUniformity(double[][] x, int m,
        POINT_DISTRIBUTION_TYPE type, SecureRandom rand) {
        
        if (m > x.length) {
            throw new IllegalArgumentException("m must be less than x.length");
        }
        
        int nDimensions = x[0].length;
        int nSamples = x.length;
        
        // number of L1 statistics
        int q = 100;
        
        double[] l1Sphere = new double[q];
        double[] l1X = new double[q];
        double[][] s;
        
        for (int i = 0; i < q; ++i) {
            
            l1X[i] = calcL1UniformityStatistic(x, m, type, rand);
            
            s = MultivariateUniformDistribution
                .generateUnitStandardNSphereWithRejection(nDimensions, nSamples, rand, true);
            
            l1Sphere[i] = calcL1UniformityStatistic(s, m, type, rand);
        }
        
        double[] minMaxSphere = MiscMath0.getMinMax(l1Sphere);
        double[] minMaxX = MiscMath0.getMinMax(l1X);        
        double[] avgAndStDevX = MiscMath0.getAvgAndStDev(l1X);
        double[] avgAndStDevSphere = MiscMath0.getAvgAndStDev(l1Sphere);
        
        //1-alpha for hypersphere chord distribution
        double alphaCV = findCVForAlpha95Percent(nDimensions);
        
        System.out.printf("c.v. for 1-alpha=%.4e\n", alphaCV);
        System.out.printf("S: min, max L1=%.4e : %.4e,  m=%.4e, stDev=%.4e\n", 
            minMaxSphere[0], minMaxSphere[1],
            avgAndStDevSphere[0], avgAndStDevSphere[1]);
        System.out.printf("X: min, max L1=%.4e : %.4e,  m=%.4e,, stDev=%.4e\n", 
            minMaxX[0], minMaxX[1],
            avgAndStDevX[0], avgAndStDevX[1]);
        
        // Type I error rejects a null hypothesis as false when it is actually true.
        // Type II error accepts a null hypothesis as true when it is actually false.

        // alpha=0.05 (probability of Type I error)
        // confidence level or a confidence coefficient, (1 - α)100% = 95%
        // confidence interval is interval in x capturing 95% of area 
        //     under curve, e.g. mu +- 2*sigma/sqrt(n)
        
        NonUniformityStats stats = new NonUniformityStats();
        stats.l1MaxSphere = minMaxSphere[1];
        stats.l1MaxX = minMaxX[1];
        stats.oneMinusAlpha = 0.95;
        stats.oneMinusAlphaCritVal = alphaCV;
        
        //TODO: follow up on this:
        
        // see Section IV, page 9:
        //If the α%-largest L1 value of Q is smaller than L1(g) then S can be 
        //  declared as non-uniform with confidence (100 − α)%.
        /*if (stats.l1MaxSphere > stats.l1MaxX) {
            stats.isConsistentWithNonUniform = true;
        } else {
            stats.isConsistentWithNonUniform = false;
        }*/
        //  min stats.l1MaxX minL1 > (stats.l1MaxSphere avgL1 + c*StDevSphereL1)
        //       implies not the same distributions
        if (minMaxX[0] > (avgAndStDevSphere[0] + alphaCV*avgAndStDevSphere[1])) {
            stats.isConsistentWithNonUniform = true;
        } else {
            stats.isConsistentWithNonUniform = false;
        }
        
        return stats;
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
     * @param type type of distance distribution to create internally using 
     * methods outlined about form Section IV. of 2020 Sidiropoulos.
     * @param rand
     * @return the calculated L1 statistic from eqn (26) of the paper.
    */
    public static double calcL1UniformityStatistic(double[][] x, int m,
        POINT_DISTRIBUTION_TYPE type, SecureRandom rand) {
        
        if (m > x.length) {
            throw new IllegalArgumentException("m must be less tna x.length");
        }
        
        double[] d1 = createOrderedDistanceDistribution(x, m, type, rand);
        
        double l1Sum = calcL1UniformityStatistic(d1, x[0].length);
        
        return l1Sum;
    }
    
    static double[] createOrderedDistanceDistribution(double[][] x, int m,
        POINT_DISTRIBUTION_TYPE type, SecureRandom rand) {
        
        if (m > x.length) {
            throw new IllegalArgumentException("m must be less tna x.length");
        }
        
        double[] d1;

        if (type == null || type.equals(POINT_DISTRIBUTION_TYPE.INTRA_DISTANCE_2)) {

            int[] xMIdx = chooseM(m, x.length, rand);

            d1 = calcPairwiseDistances(x, xMIdx);
            
        } else {
            // uses POINT_DISTRIBUTION_TYPE.POINT_DISTANCE_3
            //    resulting length of d is
            
            int[] xMIdx = chooseM(m + 1, x.length, rand);
            
            int pt1Idx = xMIdx[m];
            xMIdx = Arrays.copyOfRange(xMIdx, 0, m);

            d1 = calcPairwiseDistances(x, xMIdx, pt1Idx);
            
            assert(d1.length == xMIdx.length);
        }
        
        int[] indexes1 = MiscSorter.mergeSortIncreasing(d1);
        
        return d1;
    }
    
    /**
    calculate the uniformity measure based upon the point distances
    (presumable measured by methods (2) or (3) of Sective IV. of the 
    2020 Sidiropoulos paper.
    <pre>
    to quantify the ”uniformity” of an input point distribution on a N-sphere, 
    the L1 distance is used.  
        
        L_1(g) = integral_{x=0_to_2} ( | g_N(x) - f_N(x)| * dx )

           where f_N(x) is hypersphere chord length distribution == the pdf,
               and
           g_N is intra-distance distribution of the input point distribution 
               using model (2) above.
           x is the chord length d, so the integration is from d=0 to d=2 (presumably 2*r with r=1).
    </pre>
     * @param d array of point distances, sorted by non-decreasing order.
     * @param nDimensions
     * @return the calculated L1 statistic from eqn (26) of the paper.
    */
    public static double calcL1UniformityStatistic(double[] d, int nDimensions) {
        
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
        
        double[] p2 = pdf(d, 1, nDimensions);
        
        double[] minMaxD1 = MiscMath0.getMinMax(d);
        //double[] minMaxD2 = MiscMath0.getMinMax(d2);
        double minD = minMaxD1[0];
        double maxD = minMaxD1[1];
     
        HistogramHolder h1 = Histogram.calculateScottsHistogram(d);
        float[] p1 = h1.getYHistFloat();
        float[] p1x = h1.getXHist();
        float[] minMaxP1 = MiscMath0.getMinMax(p1);
        float p1Norm = 1.f/minMaxP1[1];
        for (int i = 0; i < p1.length; ++i) {
            p1[i] *= p1Norm;
        }
        
        /*try {
            PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
            plotter.addPlot(p1x, p1, p1x, p1, "p1");
            String str = plotter.writeFile("p1_hist");
            
            plotter = new PolygonAndPointPlotter();
            plotter.addPlot(d, p2, null, null, "p2");
            String str2 = plotter.writeFile("p2_hist");
        } catch (IOException ex) {
            Logger.getLogger(HypersphereChordLength.class.getName()).log(Level.SEVERE, null, ex);
        }*/
        
        float binWidth1 = p1x[1] - p1x[0];
        float p1Min = p1x[0] - (binWidth1/2.f);
        float p1Max = p1x[p1x.length - 1] + (binWidth1/2.f);
        
        double l1Sum = 0;
        double dist;
        double diff;
        int bin1;
        for (int i = 0; i < d.length; ++i) {
            dist = d[i];
            
            if (dist < 0 || dist > 2) {
                int z = 0;
            }
            
            bin1 = (int) ((dist - p1Min)/binWidth1);
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
     * @param rand
     * @return
     */
    static int[] chooseM(int m, int n, SecureRandom rand) {
             
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

    static double[] calcPairwiseDistances(double[][] x, int[] idx) {
        
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
    
    static double[] calcPairwiseDistances(double[][] x, int[] idx, int p0Idx) {
        
        if (idx.length > x.length) {
            throw new IllegalArgumentException("m cannot be larger than x.length");
        }
        
        int m = idx.length;
                
        double[] d = new double[idx.length];
               
        for (int i = 0; i < idx.length; ++i) {
            d[i] = Math.sqrt(Distances.calcEuclideanSquared(x[p0Idx], x[idx[i]]));
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
    
    /**
     * a rough critical value for which alpha is 95% for a unit
     * radius nDimension hypersphere.
     * @param nDimensions the number of dimensions of the hypersphere
     * @return  rough critical value for 95% quantile 
     */
    public static double findCVForAlpha95Percent(int nDimensions) {
        
        // quick look at finding critical values using the CDF.
        // no inverse function, so "trial-and-error".
        //    For r = 1.  n = [2:10:+1,15:50:+5,60:100:+10]
        //       find d's where alpha=0.95
        
        double r = 1;
        
        TDoubleList dList = new TDoubleArrayList();
        if (nDimensions < 2) {
            //fine resolution between 1.9 and 2.0
            for (double k = 1.99; k <= 2.0; k+=0.001) {
                dList.add(k);
            }
        } else if (nDimensions < 5) {
            //fine resolution between 1.9 and 2.0
            for (double k = 1.8; k <= 2.0; k+=0.005) {
                dList.add(k);
            }
        } else if (nDimensions < 8) {
            //fine resolution between 1.76 and 1.88
            for (double k = 1.76; k <= 1.88; k+=0.01) {
                dList.add(k);
            }
        } else if (nDimensions <= 10) {
            //fine resolution between 1.76 and 1.88
            for (double k = 1.74; k <= 1.88; k+=0.01) {
                dList.add(k);
            }
        }  else if (nDimensions < 100) {
            //fine resolution between 1.5 and ?
            for (double k = 1.5; k <= 1.75; k+=0.01) {
                dList.add(k);
            }
        } else if (nDimensions < 800) {
            //fine resolution between 1.5 and ?
            for (double k = 1.45; k <= 1.6; k+=0.005) {
                dList.add(k);
            }
        } else if (nDimensions <= 1000) {
            //fine resolution between 1.5 and ?
            for (double k = 1.45; k <= 1.455; k+=0.001) {
                dList.add(k);
            }
        } else if (nDimensions <= 2000) {
            //fine resolution between 1.5 and ?
            for (double k = 1.44; k <= 1.456; k+=0.001) {
                dList.add(k);
            }
        } else if (nDimensions <= 3050) {
            //fine resolution between 1.5 and ?
            for (double k = 1.435; k <= 1.45; k+=0.001) {
                dList.add(k);
            }
        } else if (nDimensions <= 5000) {
            //fine resolution between 1.5 and ?
            for (double k = 1.42; k <= 1.45; k+=0.001) {
                dList.add(k);
            }
        }
        double[] ds = dList.toArray();
        
        int n, idx;
        double[] cdf;
        double tol = 1e-2;
        
        // binary search between d=1.5*r and d=2*R for n <= 100
        cdf = HypersphereChordLength.cdf(ds, r, nDimensions);
        //System.out.println(FormatArray.toString(cdf, "%.9f"));

        // Type I error rejects a null hypothesis that is actually true.
        // Type II error accepts a null hypothesis that is actually false.

        // alpha=0.05 (probability of Type I error)
        // confidence level or a confidence coefficient, (1 - α)100% = 95%
        idx = CDFRandomSelect.binarySearchForNearest(cdf, 0.95, tol);
        // confidence interval is interval in x capturing 95% of area 
        //     under curve, e.g. N(0,1): mu +- 1.96*sigma/sqrt(n)
        
        return ds[idx];
    }
}
