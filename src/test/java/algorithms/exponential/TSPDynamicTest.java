package algorithms.exponential;

import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPDynamicTest extends TestCase {
    
    public TSPDynamicTest() {
    }
    
    public void test1() {
        System.out.println("test1");
         // Create adjacency matrix
        int n = 6;
        double[][] distanceMatrix = new double[n][n];
        for (double[] row : distanceMatrix) {
            java.util.Arrays.fill(row, 10000);
        }
        distanceMatrix[5][0] = 10;
        distanceMatrix[1][5] = 12;
        distanceMatrix[4][1] = 2;
        distanceMatrix[2][4] = 4;
        distanceMatrix[3][2] = 6;
        distanceMatrix[0][3] = 8;
        int startNode = 0;
        
        int[] expectedTour = new int[]{0, 3, 2, 4, 1, 5, 0};
        int expectedCost = 42;
        TSPDynamic solver = new TSPDynamic(startNode, distanceMatrix);

        List<Integer> tour = solver.getTour();
        // Prints: [0, 3, 2, 4, 1, 5, 0]
        System.out.println("Tour: " + tour);

        double cost = solver.getTourCost();
        // Print: 42.0
        System.out.println("Tour cost: " + solver.getTourCost());

        assertEquals(expectedTour.length, tour.size());
        assertTrue(Math.abs(expectedCost-cost) < 1e-17);
        
        for (int i = 0; i < expectedTour.length; ++i) {
            assertEquals(expectedTour[i], tour.get(i).intValue());
        }        
        
    }
    
    public void test2() {
        System.out.println("test2");
         // https://www.baeldung.com/cs/tsp-dynamic-programming
         
        //A B C D E
        //0 1 2 3 4
        
        int n = 5;
        double[][] distanceMatrix = new double[n][n];
        int inf = 1000;
        distanceMatrix[0] = new double[]{inf, 12, 10, 19, 8};
        distanceMatrix[1] = new double[]{12, inf, 3, 7, 6};
        distanceMatrix[2] = new double[]{10, 3, inf, 2, 20};
        distanceMatrix[3] = new double[]{19, 7, 2, inf, 4};
        distanceMatrix[4] = new double[]{8, 6, 20, 4, inf};
        
        int startNode = 0;
        
        //A B C D E
        //0 1 2 3 4
        //A, E, D, C, B -> 0, 4, 3, 2, 1, 0
        int[] expectedTour = new int[]{0, 4, 3, 2, 1, 0};
        int expectedCost = 29;//8+4+2+3+12
        TSPDynamic solver = new TSPDynamic(startNode, distanceMatrix);

        for (int ii = 0; ii < 2; ++ii) {
            if (ii == 0) {
                solver.solveIteratively();
            } else {
                solver.solveRecursively();
            }

            List<Integer> tour = solver.getTour();
            System.out.println("Tour: " + tour);

            double cost = solver.getTourCost();
            System.out.println("Tour cost: " + solver.getTourCost());

            assertEquals(expectedTour.length, tour.size());
            assertTrue(Math.abs(expectedCost - cost) < 1e-17);

            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour.get(i).intValue());
            }
        }
        
    }
}
