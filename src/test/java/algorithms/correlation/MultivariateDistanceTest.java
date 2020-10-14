package algorithms.correlation;

import junit.framework.TestCase;

public class MultivariateDistanceTest extends TestCase {
    
    public MultivariateDistanceTest(String testName) {
        super(testName);
    }

  
    public void testEfficientDCov() {
        /*
        Distance Based Independence Tests in Secion 2.3 of paper
        "A Statistically and Numerically Efficient Independence Test Based On
        Random Projections and Distance Covariance"
        2017, Huang and Hu
        
        (n * covsq / term3) > (InverseNormalCDF(1 - (alpha_s/2))^2
             where alpha_s is .gt. 0 and .lt. 0.215
                 e.g. alpha_s = 0.05
        
        For multivariate test, see Section 3.2 of same paper,
        which is summarized in algorithm 2 in the appendix.
        
        test w/:
             p = q = 10 and the number of Monte Carlo iterations in RPDC is K = 50
        
        // Type I error rejects a null hypothesis as false when it is actually true.
        // Type II error accepts a null hypothesis as true when it is actually false.
        // Type 3? Asking the wrong question, making the right decision
              for the wrong reason, etc.

        */
    }


}
