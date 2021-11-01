package algorithms.tsp;

import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TSPCorrectedInterviewBitTest extends TestCase {
    
    public TSPCorrectedInterviewBitTest() {
    }
    
    public void testGreedyFail() {
        int n = 4;
        double[][] x = new double[4][];
        x[0] = new double[]{0, 0};
        x[1] = new double[]{0, 1};
        x[2] = new double[]{2, 0};
        x[3] = new double[]{3, 1};
        
        // optimal is 0, 1, 3, 2, 0 = 7.41
        // greedy is 0, 1, 2, 3, 0  = 7.81
        int i, j;
        double[] xi, xj;
        double xd, yd;
        double[][] dist = new double[4][];
        for (i = 0; i < n; ++i) {
            dist[i] = new double[4];
            xi = x[i];
            for (j = 0; j < n; ++j) {
                xj = x[j];
                xd = (xi[0] - xj[0]);
                yd = (xi[1] - xj[1]);
                dist[i][j] = Math.sqrt(xd*xd + yd*yd);
            }
        }
        
        int startNode = 0;
        int[] expectedGreedy = new int[]{0, 1, 2, 3, 0};
        int[] expectedOptimal = new int[]{0, 1, 3, 2, 0};
        expectedOptimal = new int[]{0, 2, 3, 1, 0};
        double expectedOptimalCost = 7.41;
        double expectedGreedyCost = 7.81;
    
        List<Integer> tour, tour2;
        double cost, cost2;
        
        TSPOptimalDynamic solver = new TSPOptimalDynamic(startNode, dist);
        solver.solveIteratively();        
        tour = solver.getTour();
        cost = solver.getTourCost();
        
        System.out.printf("startNode=%d  optimal: ", startNode);
        System.out.printf("optimal tour=%s, expected=%s\n ", tour.toString(), Arrays.toString(expectedOptimal));
        System.out.printf("Tour cost: %.3f, expected=%.3f\n", cost, expectedOptimalCost);
        
        TSPGreedy solver2 = new TSPGreedy(startNode, dist);
        solver2.solveIteratively();        
        tour2 = solver2.getTour();
        cost2 = solver2.getTourCost();
        System.out.printf("startNode=%d  greedy: ", startNode);
        System.out.printf("greedy tour=%s, expected=%s\n ", tour2.toString(), Arrays.toString(expectedGreedy));
        System.out.printf("Tour cost: %.3f, expected=%.3f\n", cost2, expectedGreedyCost);
        
        assertEquals(expectedOptimal.length, tour.size());
        assertEquals(expectedGreedy.length, tour2.size());
        for (i = 0; i < n; ++i) {
            assertEquals(expectedOptimal[i], tour.get(i).intValue());
            assertEquals(expectedGreedy[i], tour2.get(i).intValue());
        }
        assertTrue(Math.abs(expectedOptimalCost - cost) < 1.e-2);
        assertTrue(Math.abs(expectedGreedyCost - cost2) < 1.e-2);
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
