package algorithms.pca;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.Standardization;
import algorithms.pca.PrincipalComponents.PCAStats;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class PrincipalComponentsTest extends TestCase {
    
    public PrincipalComponentsTest(String testName) {
        super(testName);
    }
    
    public void testPCA() throws Exception {
        
        // from:
        // https://online.stat.psu.edu/stat505/book/export/html/670
        
        double[][] x = readPlaces();
        
        int i, j;
        
        //NOTE: to match the results of the psu tutorial, follow section
        //    Example 11-3: Place Rated (after Standardization)
        
        double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        System.out.println("mean x=");
        for (i = 0; i < mean.length; ++i) {
            System.out.printf("%11.3e  ", mean[i]);
        }
        System.out.println();
        System.out.flush();
        
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
        
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 3);
        
        double[][] b = MatrixUtil.multiply(x, stats.principalDirections);
        b = MatrixUtil.multiply(b, stats.vTP);
        System.out.println("b = x * pr.dir * v^T_p=");
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                System.out.printf("%11.3e  ", b[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        double[][] xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                //xMinusB[i][j] -= mean[j];
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.println("x:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        System.out.println("(x - (m + B*a))^2:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", xMinusB[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
    }
    
    public void testPCA2() throws Exception {
        
        double[][] x = new double[4][2];
        x[0] = new double[]{1, 2};
        x[1] = new double[]{2, 1};
        x[2] = new double[]{3, 4};
        x[3] = new double[]{4, 3};
        
        int i, j;
        
        double[] mean = MatrixUtil.mean(x);
        System.out.println("mean x=");
        for (i = 0; i < mean.length; ++i) {
            System.out.printf("%11.3e  ", mean[i]);
        }
        System.out.println();
        System.out.flush();
        
        /*double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        System.out.println("mean x=");
        for (i = 0; i < mean.length; ++i) {
            System.out.printf("%11.3e  ", mean[i]);
        }
        System.out.println();
        System.out.flush();
        */
        
        double n = x.length;
        
        // U is 2x2,  U_p is 2x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // x is nx2
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 1);
                
        double[][] b = MatrixUtil.multiply(x, stats.principalDirections);
        b = MatrixUtil.multiply(b, stats.vTP);
        System.out.println("b = x * pr.dir * v^T_p=");
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                System.out.printf("%11.3e  ", b[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        double[][] xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                //xMinusB[i][j] -= mean[j];
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.println("x:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        System.out.println("(x - (m + B*a))^2:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", xMinusB[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
    }
    
    public void testReconstruction() throws Exception {
        
        double d2r = Math.PI/180.;
        double angle, dx, dy;
        double[][] x = new double[4][2];
       
        angle = -30;//-45;
        double x1=3./Math.sqrt(2); double y1=1./Math.sqrt(2);
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
        
        System.out.println("x0:");
        for (i = 0; i < x.length; ++i) {
            for (j = 0; j < x[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        double[] mean = new double[x[0].length];
        double[] stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        System.out.println("mean x=");
        for (i = 0; i < mean.length; ++i) {
            System.out.printf("%11.3e  ", mean[i]);
        }
        System.out.println();
        System.out.println("stdev x=");
        for (i = 0; i < stDev.length; ++i) {
            System.out.printf("%11.3e  ", stDev[i]);
        }
        System.out.println();
        System.out.flush();
        
        double n = x.length;
        
        // U is 2x2,  U_p is 2x1
        // V is 2x2, V_P is 2x1, V^T_p is 1x2
        // x is nx2
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 1);
                
        double[][] b = MatrixUtil.multiply(x, stats.principalDirections);
        b = MatrixUtil.multiply(b, stats.vTP);
        System.out.println("b = x * pr.dir * v^T_p=");
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                System.out.printf("%11.3e  ", b[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        double[][] xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                //xMinusB[i][j] -= mean[j];
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.println("x:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        System.out.println("(x - (m + B*a))^2:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", xMinusB[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        // =====================
        System.out.println("=== new points in same reference frame ===");
        
        angle = -30;//-45;
        x1=10./Math.sqrt(2); y1=1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt1 = dx;
        yt1 = dy;
        
        x1=10./Math.sqrt(2); y1=-1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt2 = dx;
        yt2 = dy;
        
        x1=15./Math.sqrt(2); y1=1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt3 = dx;
        yt3 = dy;
        
        x1=15./Math.sqrt(2); y1=-1./Math.sqrt(2);
        dx= x1*Math.cos(angle*d2r) + y1*Math.sin(angle*d2r);
        dy= -x1*Math.sin(angle*d2r) + y1*Math.cos(angle*d2r);
        xt4 = dx;
        yt4 = dy;
        
        x[0] = new double[]{xt1, yt1};
        x[1] = new double[]{xt2, yt2};
        x[2] = new double[]{xt3, yt3};
        x[3] = new double[]{xt4, yt4};
        
        System.out.println("*x0:");
        for (i = 0; i < x.length; ++i) {
            for (j = 0; j < x[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        mean = new double[x[0].length];
        stDev = new double[x[0].length];
        x = Standardization.standardUnitNormalization(x, mean, stDev);
        
        b = MatrixUtil.multiply(x, stats.principalDirections);
        b = MatrixUtil.multiply(b, stats.vTP);
        
        double[][] bReconstruction = PrincipalComponents.reconstruct(x, stats);
        for (i = 0; i < bReconstruction.length; ++i) {
            for (j = 0; j < bReconstruction[i].length; ++j) {
                bReconstruction[i][j] *= stDev[j];
                bReconstruction[i][j] += mean[j];
            }
        }
        
        System.out.println("*b = x * pr.dir * v^T_p=");
        for (i = 0; i < b.length; ++i) {
            for (j = 0; j < b[i].length; ++j) {
                System.out.printf("%11.3e  ", b[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        System.out.println("*B =");
        for (i = 0; i < bReconstruction.length; ++i) {
            for (j = 0; j < bReconstruction[i].length; ++j) {
                System.out.printf("%11.3e  ", bReconstruction[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
        xMinusB = MatrixUtil.copy(x);
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                //xMinusB[i][j] -= mean[j];
                xMinusB[i][j] -= b[i][0];
            }
        }
        System.out.println("*x:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", x[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        System.out.println("*(x - (m + B*a))^2:");
        for (i = 0; i < xMinusB.length; ++i) {
            for (j = 0; j < xMinusB[i].length; ++j) {
                System.out.printf("%11.3e  ", xMinusB[i][j]);
            }
            System.out.println();
        }
        System.out.flush();
        
    }
    
    /**
     * get places data from a PSU statistics tutorial.  the
     * array has n=329 and dimensions=9.
     * the dimensions are 
     *     climate housing health crime trans educate arts recreate econ
     * @return
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
