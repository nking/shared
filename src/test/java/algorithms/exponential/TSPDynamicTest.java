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
        for (int i = 0; i < n; ++i) {
            distanceMatrix[i][i] = 0;
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
        TSPOptimalDynamic solver;
        TSPGreedy solver2;
        List<Integer> tour, tour2;
        double cost, cost2;
        
        int ni = 2;
        for (int ii = 0; ii < ni; ++ii) {
            if (ii < ni/2) {
                startNode = 0;
                expectedTour = expectedTour0;
            } else {
                startNode = 1;
                expectedTour = expectedTour1;
            }
            solver = new TSPOptimalDynamic(startNode, distanceMatrix);
            solver.solveIteratively();
            System.out.printf("startNode=%d  solveIteratively: ", startNode);
            tour = solver.getTour();
            System.out.println("Tour: " + tour);
            cost = solver.getTourCost();
            System.out.printf("Tour cost: %.1f, expected=%d\n", cost, expectedCost);

            solver2 = new TSPGreedy(startNode, distanceMatrix);
            solver2.solveRecursively();
            // symmetric distance matrix, so cycle in other direction also a correct answer
            System.out.printf("greedy startNode=%d  solveRecursively: ", startNode);
            tour2 = solver2.getTour();
            System.out.println("greedy Tour: " + tour2);
            cost2 = solver2.getTourCost();
            System.out.printf("greedy Tour cost: %.1f, expected=%d\n", cost2, expectedCost);

            assertEquals(expectedTour.length, tour.size());
            assertTrue(Math.abs(expectedCost - cost) < 1e-17);
            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour.get(i).intValue());
            }
            
            assertEquals(expectedTour.length, tour2.size());
            assertTrue(Math.abs(expectedCost - cost2) < 1e-17);
            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour2.get(i).intValue());
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
        int[] expectedTour1 = new int[]{1, 2, 3, 4, 0, 1};
        int[] expectedTour;
        int startNode;
        int expectedCost = 29;//8+4+2+3+12
        TSPOptimalDynamic solver;
        TSPGreedy solver2;
        List<Integer> tour, tour2;
        double cost, cost2;
        
        int ni = 2;
        for (int ii = 0; ii < ni; ++ii) {
            if (ii < ni/2) {
                startNode = 0;
                expectedTour = expectedTour0;
            } else {
                startNode = 1;
                expectedTour = expectedTour1;
            }
            solver = new TSPOptimalDynamic(startNode, distanceMatrix);
            solver.solveIteratively();
            System.out.printf("startNode=%d  solveIteratively: ", startNode);
            tour = solver.getTour();
            System.out.println("Tour: " + tour);
            cost = solver.getTourCost();
            System.out.printf("Tour cost: %.1f, expected=%d\n", cost, expectedCost);

            solver2 = new TSPGreedy(startNode, distanceMatrix);
            solver2.solveRecursively();
            // symmetric distance matrix, so cycle in other direction also a correct answer
            System.out.printf("greedy startNode=%d  solveRecursively: ", startNode);
            tour2 = solver2.getTour();
            System.out.println("greedy Tour: " + tour2);
            cost2 = solver2.getTourCost();
            System.out.printf("greedy Tour cost: %.1f, expected=%d\n", cost2, expectedCost);

            assertEquals(expectedTour.length, tour.size());
            assertTrue(Math.abs(expectedCost - cost) < 1e-17);
            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour.get(i).intValue());
            }
            
            assertEquals(expectedTour.length, tour2.size());
            assertTrue(Math.abs(expectedCost - cost2) < 1e-17);
            for (int i = 0; i < expectedTour.length; ++i) {
                assertEquals(expectedTour[i], tour2.get(i).intValue());
            }
        }
        
    }
}
