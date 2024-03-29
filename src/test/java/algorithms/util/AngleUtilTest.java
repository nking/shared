package algorithms.util;

import algorithms.util.AngleUtil;
import algorithms.util.ResourceFinder;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class AngleUtilTest extends TestCase {
    
    public AngleUtilTest() {
    }
    
    public void testSubtract() throws Exception {
        
        AngleUtil instance = new AngleUtil();
                        
        BufferedReader br = null;
        FileReader reader = null;
        
        //1 1 1.394857 0.343024 5.231353 (8.000000, -45.000000, 56.000000, -20.000000)
        
        String pattern = "(\\d)\\s(\\d)\\s([\\-0-9\\.]+)\\s([\\-0-9\\.]+)\\s" +
            "([\\-0-9\\.]+)\\s" 
            + "\\(" 
            + "([\\-\\.0-9]+),\\s([\\-\\.0-9]+),\\s([\\-\\.0-9]+),\\s([\\-\\.0-9]+)" 
            + "\\)";
        Pattern p = Pattern.compile(pattern);
        Matcher m = null;
        
        String filePath1 = ResourceFinder.findFileInTestResources("angles.tsv");
        try {
            reader = new FileReader(new File(filePath1));
            br = new BufferedReader(reader);
                        
            String line = br.readLine();
            while (line != null) {                
                m = p.matcher(line);
                if (m.matches()) {
                    double diffX1 = Double.valueOf(m.group(6)).doubleValue();
                    double diffY1 = Double.valueOf(m.group(7)).doubleValue();
                    double diffX2 = Double.valueOf(m.group(8)).doubleValue();
                    double diffY2 = Double.valueOf(m.group(9)).doubleValue();
                    
                    double expected = (2 * Math.PI) 
                        - Double.valueOf(m.group(5)).doubleValue();
                    
                    double result = instance.subtract(diffX1, diffY1, diffX2, 
                        diffY2);
                    
                    if (result < 0) {
                        result = (2 * Math.PI) + result;
                    }
                    
                    //System.out.println("line=" + line);
                    //System.out.println("expected=" + expected + " result=" 
                    //    + result);
                    
                    assertTrue(Math.abs(expected - result) < 0.1);
                }
                
                line = br.readLine();
            }
            
        } finally {
            if (reader != null) {
                reader.close();
            }
            if (br != null) {
                br.close();
            }
        }
        
    }
    
    public void testPolarAngleCCW() throws Exception {
        
        /*
                  +Y
                 90
        QII       |       QI
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QIII     |       QIV 
                 270
        */
        double x, y, expected, result;
        
        x = 10;
        y = 10;
        expected = 45.*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = 10;
        expected = (45. + 90.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = -10;
        expected = (45. + 180.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = 10;
        y = -10;
        expected = (45. + 270.)*Math.PI/180.;
        result = AngleUtil.polarAngleCCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
    }
    
    public void testPolarAngleCW() throws Exception {
        
         /*
                  +Y
                 270
        QIII      |       QIV
                  |     
                  |
     180-------------------- +X  0, 360
                  |   
                  |      
         QII      |       QI 
                  90
        */
        double x, y, expected, result;
        
        x = 10;
        y = 10;
        expected = 315.*Math.PI/180.;
        result = AngleUtil.polarAngleCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = 10;
        expected = (45. + 180.)*Math.PI/180.;
        result = AngleUtil.polarAngleCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = -10;
        y = -10;
        expected = (45. + 90.)*Math.PI/180.;
        result = AngleUtil.polarAngleCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        x = 10;
        y = -10;
        expected = 45.*Math.PI/180.;
        result = AngleUtil.polarAngleCW(x, y);
        assertTrue(Math.abs(result - expected) < 0.1);
    }
    
    public void testAngleAddition() throws Exception {
        
        boolean useRadians = false;
        double sum = AngleUtil.calcAngleAddition(0, 350, useRadians);
        
        assertTrue(Math.abs(sum - 710) < 0.1);
        
        double avg = AngleUtil.getAngleAverage(0, 350, useRadians);
        assertTrue(Math.abs(avg - 355) < 0.1);
        
        // --- same test using doubles ---
        useRadians = false;
        sum = AngleUtil.calcAngleAddition(0., 350., useRadians);
        
        assertTrue(Math.abs(sum - 710.) < 0.1);
        
        avg = AngleUtil.getAngleAverage(0., 350., useRadians);
        assertTrue(Math.abs(avg - 355.) < 0.1);
        
        
    }
    
    public void testCalculateAverageWithQuadrantCorrections() throws Exception {
        
        boolean useRadians = false;
        double[] angles = new double[]{0, 350, 340};
        double expected = (350. + 340. + 360.)/3.;
        double result = AngleUtil.calculateAverageWithQuadrantCorrections(angles,
            useRadians);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        useRadians = false;
        angles = new double[]{350, 10, 20};
        expected = (350. + 370. + 380.)/3.;
        expected -= 360;
        result = AngleUtil.calculateAverageWithQuadrantCorrections(angles,
            useRadians);
        assertTrue(Math.abs(result - expected) < 0.1);
        
        //-----
        useRadians = true;
        angles = new double[]{0, 350.*Math.PI/180., 340.*Math.PI/180.};
        double e0 = 2.*Math.PI;
        double e1 = 350.*Math.PI/180.;
        double e2 = 340.*Math.PI/180.;
        expected = (e0 + e1 + e2)/3.;
        result = AngleUtil.calculateAverageWithQuadrantCorrections(angles,
            useRadians);
        assertTrue(Math.abs(result - expected) < 0.01);
        
        useRadians = true;
        angles = new double[]{350.*Math.PI/180., 10.*Math.PI/180., 20.*Math.PI/180.};
        expected = ((350.*Math.PI/180.) + (370.*Math.PI/180.) + (380.*Math.PI/180.))/3.;
        expected -= (2.*Math.PI);
        result = AngleUtil.calculateAverageWithQuadrantCorrections(angles,
            useRadians);
        assertTrue(Math.abs(result - expected) < 0.01);
    }
    
}
