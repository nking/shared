package algorithms.pca;

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
        
        /*
        Step 1: Examine the eigenvalues to determine how many principal 
        components should be considered:
        
        Table 1. Eigenvalues and the proportion of variation explained by the 
        principal components.

            Component	Eigenvalue	Proportion	Cumulative
            1	0.3775	0.7227	0.7227
            2	0.0511	0.0977	0.8204
            3	0.0279	0.0535	0.8739
            4	0.0230	0.0440	0.9178
            5	0.0168	0.0321	0.9500
            6	0.0120	0.0229	0.9728
            7	0.0085	0.0162	0.9890
            8	0.0039	0.0075	0.9966
            9	0.0018	0.0034	1.0000
            Total	0.5225	 	 
            If you take all of these eigenvalues and add them up, then you get the total variance of 0.5223.

            The proportion of variation explained by each eigenvalue is given in the third column. For example, 0.3775 divided by the 0.5223 equals 0.7227, or, about 72% of the variation is explained by this first eigenvalue. The cumulative percentage explained is obtained by adding the successive proportions of variation explained to obtain the running total. For instance, 0.7227 plus 0.0977 equals 0.8204, and so forth. Therefore, about 82% of the variation is explained by the first two eigenvalues together.
        */
        
        PCAStats stats = PrincipalComponents.calcPrincipalComponents(x, 3);
        
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
