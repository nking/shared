package algorithms.statistics;

import java.security.SecureRandom;
import java.util.logging.Logger;

import algorithms.misc.MiscMath0;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

public class DerivGEVTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testDerivWRTX() throws Exception {
        
        double k = 1.80f;
        double sigma = 0.85f;
        double mu = 0.441f;
        double yConst = 1;
        
        double[] xp = new double[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = i/xp.length;
        }
       
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new double[0], 
            new double[0], new double[0], new double[0]);
        double[] yGEV = gev.generateCurve(xp, mu, sigma, k);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        plotter.addPlot(xp, yGEV, null, null, null, null, "");
        String filePath = plotter.writeFile();
       
        
        // generate a curve, and look at several points to see where
        //   the slope is > 0,        == 0 ,    and < 0
        //                x=0.02    x=0.0333    x=0.7
        
        /*
        for (int i = 0; i < xp.length; i++) {
            log.fine( xp[i] + ":" + DerivGEV.derivWRTX(yConst, mu, k, sigma, xp[i]));
        }
        */
        
        double d = DerivGEV.derivWRTX(yConst, mu, sigma, k, xp[0]);
        assertTrue( d > 0);
        
        // top of curve where slope should be near zero
        d = DerivGEV.derivWRTX(yConst, mu, sigma, k, 0.04278); // d=-0.000589
        assertTrue( Math.abs(d) < 0.01);
        
        // slope should be near -1 at x near 0.4
        d = DerivGEV.derivWRTX(yConst, mu, sigma, k, 0.4f);
        assertTrue( Math.abs(d + 1) < 0.1);
        
        d = DerivGEV.derivWRTX(yConst, mu, sigma, k, 0.7f);
        assertTrue( d < 0);
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100; i++) {
            double x = sr.nextDouble()*sr.nextInt(120000);
            double deriv = DerivGEV.derivWRTX(yConst, mu, sigma, k, x);
            log.fine( x + ":" + deriv);
        }
        
        /*
        double s = 0.85f;
        s = 1.9f;
        
        k = 1.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=1");
        filePath = plotter.writeFile();
        
        k = 1.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=1.5");
        filePath = plotter.writeFile();
        
        k = 2.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=2");
        filePath = plotter.writeFile();
        
        k = 2.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=2.5");
        filePath = plotter.writeFile();
        
        k = 3.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=3");
        filePath = plotter.writeFile();
        
        k = 3.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=3.5");
        filePath = plotter.writeFile();
        
        k = 4.0f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "k=4");
        filePath = plotter.writeFile();
        
        k = 4.5f;
        sigma = s;
        mu = 0.441f;
        yConst = 1;
        yGEV = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        plotter.addPlot(xp, yGEV, null, null, null, null, "");
        filePath = plotter.writeFile();
        */
    }
    
    public void testDerivWRTK() throws Exception {
        
        // k is the shape parameter
        
        double k = 1.80f;
        double sigma = 0.85f;
        double mu = 0.441f;
        double yConst = 1;
        double xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        double[] xp = new double[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = i/xp.length;
        }
        
        /*
        // ====== exploring known GEV curve and expected residual suggested for a change in k ====
        double minChiSqSum = double.MAX_VALUE;
        int minChiSqSumIdx = 0;
        double[] yg0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        double[] yg0e = new double[yg0.length];
        Arrays.fill(yg0e, 0.03f);
        double s2  = 0.85f; //1.07f;
        double mu2 = 0.441f;//0.59f;
        double k2  = 1.80f - 1.7f; //2.75f;
        double avgResidK = 0.f;
        int yMaxModelIdx = MiscMath.findYMaxIndex(GeneralizedExtremeValue.generateNormalizedCurve(xp, mu2, s2, k2));
        log.fine("yMaxModelIdx=" + yMaxModelIdx);
        double[] suggestedK = new double[xp.length];
        for (int i = 0; i < xp.length; i++) {            
            Double d = DerivGEV.derivWRTK(yConst, mu2, s2, k2, xp[i]);
            double deltaK = 0.0001f;
            Double d2 = DerivGEV.derivWRTK(yConst, mu2, s2, (k2 - deltaK), xp[i]);
            Double dd = (d2-d)/deltaK;            
            double preconditionedResidual = d2/dd;
            double chiSqSum = DerivGEV.chiSqSum(mu2, s2, (k2 + preconditionedResidual), xp, yg0, yg0e);
            double chiSqSum2 = DerivGEV.chiSqSum(mu2, s2, (k2 - preconditionedResidual), xp, yg0, yg0e);
            log.fine( String.format("x[%d]=%4.3f  (d/dk=%4.5f, d2/dkdk=%4.5f) ==> (+%4.4f  chiSqSum=%4.4f) (-%4.4f  chiSqSum=%4.4f)", 
                i, xp[i], d, dd, preconditionedResidual, chiSqSum, preconditionedResidual, chiSqSum2));
            avgResidK += preconditionedResidual;
            suggestedK[i] = (chiSqSum < chiSqSum2) ? preconditionedResidual : -1.f*preconditionedResidual;
            if (suggestedK[i] < minChiSqSum) {
                minChiSqSum =  suggestedK[i];
                minChiSqSumIdx = i;
            }
        }
        avgResidK /= xp.length;
        double stdDevSuggestedK = 0;
        for (int i = 0; i < xp.length; i++) {
            stdDevSuggestedK += Math.pow((suggestedK[i] - avgResidK), 2);
        }
        stdDevSuggestedK =  (Math.sqrt(stdDevSuggestedK/(xp.length - 1.0f)));//N-1 because had to calculate mean from the data
        double avgKWithoutOutliers = 0;
        int countWithoutOutliers = 0;
        for (int i = 0; i < xp.length; i++) {
            double diff = suggestedK[i] - avgResidK;
            if (Math.abs(diff - avgResidK) < 3.*stdDevSuggestedK) {
                avgKWithoutOutliers += suggestedK[i];
                countWithoutOutliers++;
            } else {
                log.fine("  exclude " + suggestedK[i]);
            }
        }
        avgKWithoutOutliers /= countWithoutOutliers;
        log.fine("===> avg d/dk/d2/dkdk = " + avgKWithoutOutliers 
            + "  d/dk/d2/dkdk at minChiSqSum = " + suggestedK[minChiSqSumIdx]
            + "  d/dk/d2/dkdk at yMax=" + suggestedK[yMaxModelIdx]);
        // results show it's still the best approach to take the minChiSqSum solution from the min change in d/dk
        //    the problem with this approach is we are actually needing to know what deltaK to test chisqsum for
        //    and the curve gives us d/dk modified by d2/dkdk as the suggested step in either direction that
        //    will create a change in the gev.
        //    the correct answer would test both + modified d/dk and - d/dk
        //  is the idx at position +1 past model y max always the best answer for d/dk?
        // ===> looks like a good pattern would be to get the modified d/dk at position yMaxIdx + 1
        //      and test whether adding that or subtracting it improves the chi sq sum and take the best answer.
        */
        
        double factor = 0.0001f;
        
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new double[0], 
            new double[0], new double[0], new double[0]);
        double[] yGEV0 = gev.generateCurve(xp, mu, sigma,k - (k*factor));
        double[] yGEV1 = gev.generateCurve(xp, mu, sigma, k);
        double[] yGEV2 = gev.generateCurve(xp, mu, sigma,k + (k*factor));
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
     
        
        plotter.addPlot(xp, yGEV0, null, null, null, null, "k-deltaK");
        plotter.addPlot(xp, yGEV1, null, null, null, null, "k=" + k);
        plotter.addPlot(xp, yGEV2, null, null, null, null, "k+deltaK");
        String filePath = plotter.writeFile();
       
        //double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (k*factor) * i;
            double d = (k + delta);
                        
            // looks like the deriv increases with increasing shape k for this region of the curve
            
            double deriv = DerivGEV.derivWRTK(yConst, mu, sigma, d, xPoint);
            
            log.fine("k=" + d + "  derivWRTK=" + deriv);

            // these are all nearly the same value since x is the same and delta k is small
            //assertTrue(last < deriv.doubleValue());
            
            double deriv2 = DerivGEV.estimateDerivUsingDeltaK(mu, sigma, k, xPoint);
            
            log.fine("  compare derivWRTK to  estimateDerivUsingDeltaK " + deriv + " : " + deriv2);
            
            //last = deriv.doubleValue();
        }
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            double x = sr.nextDouble()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextDouble()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTK(yConst, mu, sigma, k, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaK(mu, sigma, k, x);
            assertNotNull(deriv2);
                        
            double m = sr.nextDouble()*sr.nextInt((int)Math.ceil(x));
            double s = sr.nextDouble()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextDouble()*sr.nextInt(10);
            }
            double k_ = sr.nextDouble()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextDouble()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTK(yConst, m, s, k_, x);
            assertNotNull(deriv);

            //log.info( x + ":" + deriv);
        }
        
    }

    public void testDerivWRTSigma() throws Exception {
        
        log.info("testDerivWRTSigma");
        
        // sigma is the scale factor
        
        double k = 1.80f;
        double sigma = 0.85f;
        double mu = 0.441f;
        double yConst = 1;
        double xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        double[] xp = new double[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = i/xp.length;
        }

        double factor = 0.0001f;
        
        double[] yGEV0 = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, (sigma - (sigma*factor)), k);
        double[] yGEV1 = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, sigma, k);
        double[] yGEV2 = GeneralizedExtremeValue.generateNormalizedCurve(xp, mu, (sigma + (sigma*factor)), k);
                    
        double[] yGEV = GeneralizedExtremeValue.genCurve(xp, mu, sigma, k);
        int idx = MiscMath0.findYMaxIndex(yGEV);
        assertTrue(idx > -1);
        assertTrue(idx < yGEV.length);
        //idx = 12;
        double yFactor = yGEV[idx];
        double xForYFactor = xp[idx];
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);
        plotter.addPlot(xp, yGEV0, null, null, null, null, "sigma-delta");
        plotter.addPlot(xp, yGEV1, null, null, null, null, "sigma=" + sigma);
        plotter.addPlot(xp, yGEV2, null, null, null, null, "sigma+delta");
        String filePath = plotter.writeFile();
       
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (sigma*factor) * i;
            double d =  (sigma + delta);
                        
            // looks like the deriv increases with increasing scale sigma for this region of the curve
            
            Double deriv = DerivGEV.derivWRTSigma(yConst, mu, d, k, xPoint);
            if (deriv != null) {
                log.fine("sigma=" + d + "  derivWRTSigma=" + deriv);
            }
            
            assertNotNull(deriv);
            
            assertTrue(last < deriv.doubleValue());
            
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaSigma(mu, sigma, k, xPoint);
            
            log.fine("  compare derivWRTSigma to  estimateDerivUsingDeltaSigma " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();           
        }
        log.fine("==> yFactor=" + yFactor + " delta=" + sigma*factor);
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            double x = sr.nextDouble()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextDouble()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTSigma(yConst, mu, sigma, k, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaSigma(mu, sigma, k, x);
            assertNotNull(deriv2);
            //log.fine( x + ":" + deriv);
            
            double m = sr.nextDouble()*sr.nextInt((int)Math.ceil(x));
            double s = sr.nextDouble()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextDouble()*sr.nextInt(10);
            }
            double k_ = sr.nextDouble()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextDouble()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTSigma(yConst, m, s, k_, x);
            assertNotNull(deriv);
        }
    }

    public void testDerivWRTMu() throws Exception {
        
        log.info("testDerivWRTMu");
        
        // mu is the location variable
        
        double k = 1.80f;
        double sigma = 0.85f;
        double mu = 0.441f;
        double yConst = 1;
        double xPoint = 0.4f; // point at which slope is approx -1 for dy/dx
        
        double[] xp = new double[30];

        for (int i = 0; i < xp.length; i++) {
            xp[i] = i/xp.length;
        }

        double factor = 0.0001f;
        
        GeneralizedExtremeValue gev = new GeneralizedExtremeValue(new double[0], 
            new double[0], new double[0], new double[0]);
                    
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(0.f, 1.0f, 0f, 1.3f);        
        
        double last = Integer.MIN_VALUE;
        for (int i = 1; i < 10; i++) {
            double delta = (mu*factor) * i;
            double d = (mu + delta);
                                    
            Double deriv = DerivGEV.derivWRTMu(yConst, d, sigma, k, xPoint);
            if (deriv != null) {
                log.info("mu=" + d + "  derivWRTMu=" + deriv);
            }
            double[] yGEV1 = gev.generateCurve(xp, d, sigma, k);
            plotter.addPlot(xp, yGEV1, null, null, null, null, "mu=" + d);
            
            assertNotNull(deriv);
            
            //assertTrue(last < deriv.doubleValue());
            

            Double deriv2 = DerivGEV.estimateDerivUsingDeltaMu(mu, sigma, k, xPoint);
            
            log.fine("  compare derivWRTMu to  estimateDerivUsingDeltaMu " + deriv + " : " + deriv2);
            
            last = deriv.doubleValue();       
        }
        String filePath = plotter.writeFile();
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.currentTimeMillis();
        sr.setSeed( seed);
        log.info("seed=" + seed);
        
        for (int i = 0; i < 100000; i++) {
            double x = sr.nextDouble()*sr.nextInt(120000);
            while (x == 0) {
                x = sr.nextDouble()*sr.nextInt(120000);
            }
            Double deriv = DerivGEV.derivWRTMu(yConst, mu, sigma, k, x);
            assertNotNull(deriv);
            Double deriv2 = DerivGEV.estimateDerivUsingDeltaMu(mu, sigma, k, x);
            assertNotNull(deriv2);
            //log.fine( x + ":" + deriv);
            
            double m = sr.nextDouble()*sr.nextInt((int)Math.ceil(x));
            double s = sr.nextDouble()*sr.nextInt(10);
            while (s == 0) {
                s = sr.nextDouble()*sr.nextInt(10);
            }
            double k_ = sr.nextDouble()*sr.nextInt(10);
            while (k_ == 0) {
                k_ = sr.nextDouble()*sr.nextInt(10);
            }
            deriv = DerivGEV.derivWRTMu(yConst, m, s, k_, x);
            assertNotNull(deriv);
        }
    }

}
