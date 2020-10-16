package algorithms.correlation;

import java.security.SecureRandom;
import junit.framework.TestCase;

public class MultivariateDistanceTest extends TestCase {
    
    private SecureRandom rand = null;
    
    public MultivariateDistanceTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp(); 
        
        rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //seed = 1318454033349663L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
    }
    
  
    /*
    Distance Based Independence Tests in Secion 5 of paper
    "A Statistically and Numerically Efficient Independence Test Based On
    Random Projections and Distance Covariance"
    2017, Huang and Hu
    https://arxiv.org/pdf/1701.06054.pdf


    test w/:
         αs = 0.05
         p = q = 10 and 

         can use K = 50 Monte Carlo iterations for the independent cases 

    N=400
    sampleSize: 100, 500, 1000, 5000, 10000 ?

    EXAMPLE 5.1 (fig 1)
     test that it finds independence:
        X vector of dimension 10, with each entry drawn as U(0,1), and no final normalization
            (i.e. within 10 dimensional hypercube)
        Y vector = U(0,1)^2 for each entry
        (estimator converges to 0 for all K)

    EXAMPLE 5.2 (fig 2)
    test dependence only in the first 2 variables:
        same X as above.
    Y: y0=x0^2, y1=x1^2, then for remaining entries draw U(0,1)^2 
        (estimator converges to 0 as sample size n grows in the independent case; and 
         converges to some nonzero number as the number of the Monte Carlo 
         iterations K grows in the dependent case)

    EXAMPLE 5.3 (fig 3)
    test that it finds independence:
        X vector of dimension p, with each entry drawn as U(0,1), and no final normalization
            (i.e. within 10 dimensional hypercube)
        Y vector of dimension q where q != p
          vector = U(0,1)^2 for each entry
        (estimator converges to 0 quickly)

    EXAMPLE 5.4 (fig 4, Table 1)
    test dependence in a small number of entries in X and Y
    (which means that the dependency structure between X and Y is 
    low-dimensional though X or Y could be of high dimensions.)
        X vector of dimension p, with each entry drawn as U(0,1).
        Y vector of dimension q, 
             with 1st 5 entires being X^2 entries.
             remaining entries being U(0,1)^2
        (estimator power decreases as p and q increase for fixed nSample size)
        (it may fail to detect the low dimensional dependency in high 
        dimensional data. A possible remedy for this issue is performing 
        dimension reduction before applying the proposed method.)


    EXAMPLE 5.1
    We set the dimension of the data to be p = q = 10. 
    We generate random vectors X ∈ R10 and Y ∈ R10 from the standard 
    multivariate normal distribution N (0, I10). 
    The joint distribution of (X, Y ) is also normal and we have 
    Cor(Xi, Yi) = ρ,i = 1,...,10, and the rest correlation are all 0. 
    We set the value of ρ to be 0 and 0.1 to represent independent and 
    correlated scenarios, respectively. The sample size n is set to be from 100 to 1500 with an increment of 100.


    EXAMPLE 5.5 (fig 5)
    We set the dimension of data to be p = q = 10. 
    We generate random vector X ∈ R10 from the standard multivariate 
    normal distribution N (0, I10 ). Let the i-th entry of Y be 
    Yi = log(Xi2) + εi, i = 1, . . . , q, where εi’s are independent 
    random errors, εi ∼ N (0, σ2). 
    We set the value of σ to be 1 and 3 to represent low and high 
    noise ratios, respectively. In the σ = 1 case, the sample size n 
    is from 100 to 1000 with an increment 20; 
    and in the σ = 3 case, the sample size n is from 100 to 4000 
    with an increment 100.

    
    EXAMPLE 5.6. (fig 6)
    We set the dimension of data to be p = q = 10. 
    We generate random vector X ∈ R10 from the standard multivariate 
    normal distribution N(0,I10). Let the i-th entry of Y be 
    Yi = log(Xi2) + ǫi,i = 1,...,q, where ǫi’s
    are independent random errors, ǫi ∼ N(0,σ2).
    We set the value of σ to be 1 and 3 to represent low and high noise ratios, 
    respectively. In the σ = 1 case, the sample size n is from 100 to 1000 
    with an increment 20; and in the σ = 3 case, the sample size n is 
    from 100 to 4000 with an increment 100.

    
    EXAMPLE 5.7. (fig 7)
    We set the dimension of data to be p = q = 10. We generate random 
    vector Xt ∈ R10, t = 1, . . . , n, from the standard multivariate 
    normal distributionN(0,I10).Letthei-thentryofYt be
    Yt,i =log(Z2 )+ǫt,i,t=1,...,T t,i
    and Yt,i = log(X2 ) + ǫt,i,t = T + 1,...,n, 
    where Zt i.i.d. ∼ N(0,I10) and ǫt,i’s are independent random errors, 
    ǫt,i ∼ N(0,1). We set the value of T to be 0.5n and 0.8n to 
    represent early and late dependency transition, respectively. 
    In the early change case, the sample size n is from 500 to 2000 
    with an increment 100; and in the late change case, the sample 
    size n is from 500 to 4000 with an increment 100.


    // Type I error rejects a null hypothesis as false when it is actually true.
    // Type II error accepts a null hypothesis as true when it is actually false.
    // Type 3? Asking the wrong question, making the right decision
          for the wrong reason, etc.
    */

    public void test1() {
        
        // x, y independent
        
        double tolCov = 1.e-3;
        
        double alpha = 0.05;
        int p = 10;
        int q = 10;
        
        int k = 50; // for independent X and Y
        
        int[] nSamples = new int[]{100, 500, 1000, 5000
            //, 10000
        };
        
        double[][] x, y;
        
        int nIter = 1;//400;
        
        double[] dcov = new double[nSamples.length];
        boolean indep1, indep2;
        
        int n = 0;
        for (int nSample : nSamples) {
            x = createWithinUnitHypercube(p, nSample);
            y = createWithinUnitHypercube(q, nSample);
            for (int i = 0; i < y.length; ++i) {
                for (int j = 0; j < y[i].length; ++j) {
                    y[i][j] *= y[i][j];
                }
            }
            System.out.printf("nSample=%d, k=%d, nIter=%d alpha=%.4e\n", 
                nSample, k, nIter, alpha);
            System.out.flush();
            
            dcov[n] = MultivariateDistance.efficientDCov(x, y, k, rand);
            //estimator should converge to 0 for all k
            
            indep1 = MultivariateDistance.areIndependent1(x, y, k, nIter, alpha, rand);
            indep2 = MultivariateDistance.areIndependent2(x, y, k, alpha, rand);
            
            System.out.printf("cov=%.4e indep=(%b, %b)\n", dcov[n], indep1, indep2);
            System.out.flush();

            assertTrue(indep1);
            assertTrue(indep2);
            
            ++n;
        }
        
        assertTrue(Math.abs(dcov[nSamples.length - 1]) < tolCov);
        
    }

    private double[][] createWithinUnitHypercube(int nCols, int nRows) {
        double[][] a = new double[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                a[i][j] = rand.nextDouble();
            }
        }
        return a;
    }

}
