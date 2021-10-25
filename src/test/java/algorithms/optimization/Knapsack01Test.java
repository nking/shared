package algorithms.optimization;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class Knapsack01Test extends TestCase {
    
    public Knapsack01Test(String testName) {
        super(testName);
    }
    
    public void test011() {

        // values
        int[] ratingsOfObjectOfInterest = new int[]{30, 14, 16, 9};

        // weights
        int[] distances = new int[]{6,  3,  4, 2};

        int capacity = 10;

        int expSumValue = 46;

        int result = Knapsack01.solveDynamically0(ratingsOfObjectOfInterest, distances, capacity);

        assertTrue(result == expSumValue);
        
        result = Knapsack01.approxDynamically(ratingsOfObjectOfInterest, distances, capacity);

        assertTrue(result == expSumValue);
        
        /*
        memo[5][11]
              0     1     2     3     4     5     6     7     8     9     10
        0     0     0     0     0     0     0     0     0     0     0     0
        1     0     0     0     0     0     0    30    30    30    30    30
        2     0     0     0    14    14    14    30    30    30    44    44
        3     0     0     0    14    16    16    30    30    30    44    46
        4     0     0     9    14    16    23    30    30    39    44    46
        
        memo[5][6] 
            0     1     2     3     4     5
        0   0     0     0     0     0     0
        1   0     0     0    30    30    30
        2   0     0    14    30    30    44
        3   0     0    16    30    30    46
        4   0     9    16    30    39    46
        */
    }
    
    public void test012() {

        int[] weights = new int[] {
            9, 13, 153, 50, 15, 68, 27, 39, 23, 52,
            11,32,24,48,73,42,43,22,7,18,4,30 
        };
        int[] values = new int[] {
            150,35,200,160,60,45,60,40,30,10,70,30,
            15,10,40,70,75,80,20,12,50,10
        };

        int result;
        int expSumValue = 1030;

        result = Knapsack01.solveDynamically0(values, weights, 400);
        
        System.out.println("result: " +result);

        assertTrue(result == expSumValue);
        
        result = Knapsack01.approxDynamically(values, weights, 400);

        System.out.println("result approx: " +result);
                
    }
        
}
