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
        
        int[] expectedTour0 = new int[]{0, 3, 2, 4, 1, 5, 0};
        int[] expectedTour1 = new int[]{1, 5, 0, 3, 2, 4, 1};
        int[] expectedTour;
        int startNode;
        int expectedCost = 42;
        TSPDynamic solver;
        
        for (int ii = 0; ii < 4; ++ii) {
            if ((ii & 1) == 0) { // even ii
                startNode = 0;
                solver = new TSPDynamic(startNode, distanceMatrix);
                solver.solveIteratively();
                expectedTour = expectedTour0;
            } else {
                startNode = 1;
                solver = new TSPDynamic(startNode, distanceMatrix);
                solver.solveRecursively();
                expectedTour = expectedTour1;
            }
            solver.solveIteratively();
            List<Integer> tour = solver.getTour();
            // Prints: [0, 3, 2, 4, 1, 5, 0]
            System.out.println("Tour: " + tour);

            double cost = solver.getTourCost();
            // Print: 42.0
            System.out.println("Tour cost: " + solver.getTourCost());

            assertEquals(expectedTour.length, tour.size());
            assertTrue(Math.abs(expectedCost - cost) < 1e-17);

            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour.get(i).intValue());
            }
        }
    }
    
    public void test2() {
        System.out.println("test2");
         // https://www.baeldung.com/cs/tsp-dynamic-programming
         
        //A B C D E
        //0 1 2 3 4
        
        int n = 5;
        double[][] distanceMatrix = new double[n][n];
        int zed = 0;
        distanceMatrix[0] = new double[]{zed, 12, 10, 19, 8};
        distanceMatrix[1] = new double[]{12, zed, 3, 7, 6};
        distanceMatrix[2] = new double[]{10, 3, zed, 2, 20};
        distanceMatrix[3] = new double[]{19, 7, 2, zed, 4};
        distanceMatrix[4] = new double[]{8, 6, 20, 4, zed};
                
        //A B C D E
        //0 1 2 3 4
        //A, E, D, C, B -> 0, 4, 3, 2, 1, 0
        int[] expectedTour0 = new int[]{0, 4, 3, 2, 1, 0};
        int[] expectedTour1 = new int[]{1, 0, 4, 3, 2, 1};
        int[] expectedTour;
        int startNode;
        int expectedCost = 29;//8+4+2+3+12
        TSPDynamic solver;

        for (int ii = 0; ii < 4; ++ii) {
            if ((ii & 1) == 0) { // even ii
                startNode = 0;
                solver = new TSPDynamic(startNode, distanceMatrix);
                solver.solveIteratively();
                expectedTour = expectedTour0;
            } else {
                startNode = 1;
                solver = new TSPDynamic(startNode, distanceMatrix);
                solver.solveRecursively();
                expectedTour = expectedTour1;
            }

            List<Integer> tour = solver.getTour();
            System.out.println("Tour: " + tour);

            double cost = solver.getTourCost();
            System.out.println("Tour cost: " + solver.getTourCost());
/*
            assertEquals(expectedTour.length, tour.size());
            assertTrue(Math.abs(expectedCost - cost) < 1e-17);

            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour.get(i).intValue());
            }*/
        }
        
    }
}
