package algorithms.dimensionReduction;

import algorithms.statistics.BruteForce;
import algorithms.statistics.MultivariateDistance;
import algorithms.dimensionReduction.PrincipalComponents.PCAStats;
import algorithms.matrix.MatrixUtil;
import algorithms.statistics.Standardization;
import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.EVD;

/**
 *
 * @author nichole
 */
public class PrincipalComponentsTest extends TestCase {
    
    public PrincipalComponentsTest(String testName) {
        super(testName);
    }
    
    public void testPCAZC() throws Exception {

        System.out.println("testPCAZC");
        
        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670
        
        double[][] x = readPlaces();
        
        int i, j;
        
        //NOTE: to match the results of the psu tutorial, follow section
        //    Example 11-3: Place Rated (after Standardization)

        /*
        double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        System.out.printf("mean x=\n%s\n", FormatArray.toString(mean, "%.5e"));
        System.out.flush();
         */

        double[] mean = MatrixUtil.columnMeans(x);

        double[][] xZC = Standardization.zeroCenterMean(x);

        /*
        double[] mean = MatrixUtil.rowMeans(x);
        double m;
        for (i = 0 ; i < x.length; ++i) {
            m = mean[i];
            for (j = 0 ; j < x[i].length; ++j) {
                x[i][j] -= m;
            }
        }*/

        /*
        Step 1: Examine the eigenvalues to determine how many principal 
        components should be considered:
        
        Table 1. Eigenvalues and the proportion of variation explained by the 
        principal components.

            Component	Eigenvalue	Proportion	Cumulative
                1	3.2978	0.3664	0.3664
                2	1.2136	0.1348	0.5013
                3	1.1055	0.1228	0.6241
                4	0.9073	0.1008	0.7249
                5	0.8606	0.0956	0.8205
                6	0.5622	0.0625	0.8830
                7	0.4838	0.0538	0.9368
                8	0.3181	0.0353	0.9721
                9	0.2511	0.0279	1.0000	
            The first principal component explains about 37% of the variation. 
            Furthermore, the first four principal components explain 72%, while 
            the first five principal components explain 82% of the variation.         
        */

        System.out.printf("xZX is [%d X %x]\n", xZC.length, xZC[0].length);
        
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(xZC, 3);

        System.out.printf("mean of A =\n%s\n", FormatArray.toString(mean, "%.4e"));
        System.out.printf("pA=\n%s\n", FormatArray.toString(stats.principalAxes, "%.4e"));
        System.out.printf("pC=\n%s\n", FormatArray.toString(stats.principalComponents, "%.4e"));

        double[] expectedEig = new double[]{0.3775, 0.0511, 0.0279, 0.0230, 0.0168, 0.0120, 0.0085, 0.0039, 0.0018};
        double[] expectedFracEig = new double[]{0.7227, 0.0977, 0.0535, 0.0440, 0.0321, 0.0229, 0.0162, 0.0075, 0.0034};
        double[] expectedCumulativeFracEig = new double[]{0.7227, 0.8204, 0.8739, 0.9178, 0.9500, 0.9728,
            0.9890, 0.9966, 1.00};
        double[] expectedPA0 = new double[]{
                0.03551, 0.0933, 0.4078, 0.1004, 0.1501, 0.0321, 0.8743, 0.1590, 0.0195
        };

        double tol = 1E-3;
        assertEquals(expectedEig.length, stats.eigenValues.length);
        for (i = 0; i < expectedEig.length; ++i) {
            assertTrue(Math.abs(expectedEig[i] - stats.eigenValues[i]) < tol);
        }


        /*
        correlation between stats.pC and the original variables:

        	            Principal Component
            Variable	1	2	3
            --------   --  --  --
            Climate	0.190	0.017	0.207
            Housing	0.544	0.020	0.204
            Health	0.782	-0.605	0.144
            Crime	0.365	0.294	0.585
            Transportation	0.585	0.085	0.234
            Education	0.394	-0.273	0.027
            Arts	0.985	0.126	-0.111
            Recreation	0.520	0.402	0.519
            Economy	0.142	0.150	0.239
         */
        // projected [297 X 9]
        // X [297 X 9]
        // correlation between projected and X
        int k;
        double[][] dCor = MultivariateDistance.fastDCor(xZC, stats.principalComponents);

        // correlation and fast distance correlation are similar, except dCor ~ abs(cor)
        System.out.printf("dCor=\n%s\n", FormatArray.toString(dCor, "%.4e"));

        double[][] expectedCorr = new double[9][];
        expectedCorr[0] = new double[]{0.190,	0.017,	0.207};
        expectedCorr[1] = new double[]{0.544,	0.020,	0.204};
        expectedCorr[2] = new double[]{0.782,	-0.605,	0.144};
        expectedCorr[3] = new double[]{0.365,	0.294,	0.585};
        expectedCorr[4] = new double[]{0.585,	0.085,	0.234};
        expectedCorr[5] = new double[]{0.394,	-0.273,	0.027};
        expectedCorr[6] = new double[]{0.985,	0.126,	-0.111};
        expectedCorr[7] = new double[]{0.520,	0.402,	0.519};
        expectedCorr[8] = new double[]{0.142,	0.150,	0.239};

        double diff, e;
        assertEquals(expectedCorr.length, dCor.length);
        int n10 = 0;
        int n20 = 0;
        int n30 = 0;
        int n40 = 0;
        for (i = 0; i < expectedCorr.length; ++i) {
            assertEquals(expectedCorr[i].length, dCor[i].length);
            for (j = 0; j < expectedCorr[i].length; ++j) {
                // these are roughly equivalent... diff is 10-20% of expected
                e = Math.abs(expectedCorr[i][j]);
                diff = Math.abs(e - dCor[i][j]);
                if (diff <= 0.1*e) {
                    ++n10;
                } else
                if (diff <= 0.2*e) {
                    ++n20;
                } else
                if (diff <= 0.3*e) {
                    ++n30;
                } else
                if (diff <= 0.4*e) {
                    ++n40;
                }
                //System.out.printf("%.3e %.3e => %.3e => %.3e\n", expectedCorr[i][j], dCor[i][j], diff, 0.1*e);
            }
        }
        System.out.printf("out of %d: %d within 10%% of expected, then %d within 20%% of expected, " +
                        "then %d within 30%% of expected, then %d within 40%% of expected\n",
                expectedCorr.length*expectedCorr[0].length, n10, n20, n30, n40);

        double[][] xR = PrincipalComponents.reconstruct(mean,stats);

        int n = x.length;
        int p = stats.nComponents;
        k = x[0].length;
        double ssd = 0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < p; ++j) {
                diff = x[i][j] - xR[i][j];
                ssd += (diff * diff);
                System.out.printf("x[%d][%d] = %.4e   diff = %.4e\n", i, j, x[i][j], diff);
            }
        }
        System.out.printf("ssd=%.4e\n", ssd);
        System.out.flush();

        // quick look at the CUR method
        boolean useCUR = true;
        PCAStats stats2 = PrincipalComponents.calcPrincipalComponents(xZC, stats.nComponents, useCUR);
        System.out.printf("principal axes (SVD)=\n%s\n", FormatArray.toString(stats.principalAxes, "%.4e"));
        System.out.printf("principal axes (CUR)=\n%s\n", FormatArray.toString(stats2.principalAxes, "%.4e"));
        //System.out.printf("principal projections (=u_p*s)=\n%s\n", FormatArray.toString(pC, "%.4e"));
        System.out.flush();

        // quick look at the number of components for 85%
        boolean useEVProp = true;
        double prop = 0.85;
        PCAStats stats3 = PrincipalComponents.calcPrincipalComponents(xZC, useEVProp, prop);
        assertEquals(3, stats3.nComponents);

    }

    public void testPCAUSN() throws Exception {

        System.out.println("testPCAUSN");

        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670

        double[][] x = readPlaces();

        int i, j;

        //NOTE: to match the results of the psu tutorial, follow section
        //    Example 11-5

        // correlation values furthest from 0 (positive or negative)
        // are the most strongly correlated.
        // for this table just printed, |correlation| >= 0.5 are significant.

        // principal component scores = principal components
        // X * pa^T
        //  [nx9] * [9x3] *= [nX3]
        //double[][] y = MatrixUtil.multiply(x, MatrixUtil.transpose(stats.principalAxes));
        //System.out.printf("principal component scores=\n%s\n", FormatArray.toString(y, "%.4e"));

        // variance-covariance(standardized data) == correlation(unstandardized data).
        //                                        == correlation(zero mean centered data).
        // therefore, pca using the standardized data == pca using the correlation matrix.
        // eigen of cov(standardized) == eigen of cor(unstandardized)

        System.out.printf("\ntestPCA using unit normal standardization\n");

        double[][] x2 = readPlaces();
        double[][] cor = BruteForce.correlation(x);

        double[] outMean = new double[x2[0].length];
        double[] outS = new double[x2[0].length];
        double[][] z = Standardization.standardUnitNormalization(x2, outMean, outS);
        double[][] cov = BruteForce.covariance(z);

        EVD evdCov = EVD.factorize(new DenseMatrix(cov));
        EVD evdCor = EVD.factorize(new DenseMatrix(cor));
        System.out.printf("EVD(cov(z))=\n%s\n", FormatArray.toString(evdCov.getRealEigenvalues(), "%.3e"));
        System.out.printf("EVD(cor(x2s))=\n%s\n", FormatArray.toString(evdCor.getRealEigenvalues(), "%.3e"));

        /*
        cov = (1/(n-1)) * X^T*X where X has been zero-centered
        cor_i_j = cov(X_i, X_j) / (sigma_i * sigma_j)
        */

        System.out.println("unit standardized Z:");
        PCAStats stats2 = PrincipalComponents.calcPrincipalComponents(z, 2);
        System.out.printf("Z pA=\n%s\n", FormatArray.toString(stats2.principalAxes, "%.4e"));
        System.out.printf("Z pC=\n%s\n", FormatArray.toString(stats2.principalComponents, "%.4e"));

    }
    
    public void testPCA2() throws Exception {
        
        double[][] x = new double[4][2];
        x[0] = new double[]{1, 2};
        x[1] = new double[]{2, 1};
        x[2] = new double[]{3, 4};
        x[3] = new double[]{4, 3};

        double[][] xT = MatrixUtil.transpose(x);

        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xT[0], xT[1], null, null, "X");
        
        int i, j;
        
        double[][] xZC = Standardization.zeroCenterMean(x);

        double[][] xZCT = MatrixUtil.transpose(xZC);
        plot.addPlot(xZCT[0], xZCT[1], null, null, "X-m");

        double n = x.length;
        int p = 1;
        // x is 4 X 2.
        // p = 1
        // U is 4x4,  U_p is 4x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // U_p * diag(s) = [4x1][1x1]=[4x1]
        // X * V_p = [4x2]*[2x1]=[4x1]
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(xZC, p);

        plot.writeFile("_pca_");
    }

    public void testReconstruction() throws Exception {

        System.out.println("testReconstruction");
        
        double d2r = Math.PI/180.;
        double angle, dx, dy;
        double[][] x = new double[4][2];
       
        angle = -30;//-45;
        double x1 = 3./Math.sqrt(2);
        double y1 = 1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt1 = dx;
        double yt1 = dy;
        
        y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt2 = dx;
        double yt2 = dy;
        
        x1=7./Math.sqrt(2); y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt3 = dx;
        double yt3 = dy;
        
        y1*=-1.;
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        double xt4 = dx;
        double yt4 = dy;
       
        x[0] = new double[]{xt1, yt1};
        x[1] = new double[]{xt2, yt2};
        x[2] = new double[]{xt3, yt3};
        x[3] = new double[]{xt4, yt4};
        
        int i, j;
        
        System.out.printf("x0:\n%s\n", FormatArray.toString(x, "%11.3e"));
        System.out.flush();

        double[] mean = MatrixUtil.columnMeans(x);

        double[][] xZC = Standardization.zeroCenterMean(x);

        System.out.printf("mean = \n%s\n", FormatArray.toString(mean, "%11.3e"));
        System.out.flush();
        
        double n = x.length;
        
        // U is 2x2,  U_p is 2x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // x is nx2
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(xZC, 2);
                
        double[][] b = stats.principalComponents;

        double[][] xR = PrincipalComponents.reconstruct(mean, stats);

        int p = stats.nComponents;
        int k = x[0].length;
        double diff;
        double ssd = 0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j < p; ++j) {
                diff = x[i][j] - xR[i][j];
                ssd += (diff * diff);
                System.out.printf("x[%d][%d] = %.4e   diff = %.4e\n", i, j, x[i][j], diff);
            }
        }
        System.out.printf("ssd=%.4e\n", ssd);
        System.out.flush();

    }
    
    /**
     * get places data from a PSU statistics tutorial.  the
     * array has n=329 and dimensions=9.
     * the dimensions are 
     *     climate housing health crime trans educate arts recreate econ
     @return
     * @throws IOException 
     */
    private double[][] readPlaces() throws IOException {
        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670
        
        double[][] x = new double[329][9];
        for (int i = 0; i < 329; ++i) {
            x[i] = new double[9];
        }
        
        String path = ResourceFinder.findFileInTestResources("places.txt");

        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(path)));
            
            String line = null;

            String pattern = "^(\\d+)";
            for (int i = 0; i < 9; ++i) {
                pattern = pattern + "\\s+(\\d+)";
            }
            pattern = pattern + "$";
            Pattern p = Pattern.compile(pattern);
            Matcher m = null;
            line = in.readLine();
            do {
                line = line.trim();
                //521  6200  237  923 4031 2757   996 1405 7633   1
                //521  6200  237  923 4031 2757   996 1405 7633   1
                //575  8138 1656  886 4883 2438  5564 2632 4350   2
                //468  7339  618  970 2531 2560   237  859 5250   3
                m = p.matcher(line);
                if (m.matches()) {
                    for (int c = 1; c <= 9; ++c) {
                        String s = m.group(c);
                        x[count][c-1] = Integer.valueOf(s);
                        x[count][c-1] = Math.log10(x[count][c-1]);
                    }
                    count++;
                }
                line = in.readLine();
            } while (line != null);
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }

        return x;
    }
}
