package thirdparty.dlib.optimization;

/**
 ported to java from dlib
 Copyright (C) 2008  Davis E. King (davis@dlib.net)
 License: Boost Software License   See LICENSE.txt for the full license.
  
 */
public class ObjectiveDeltaStopStrategy {
    
    private double minDelta = 1.e-7;
    private boolean beenUsed = false;
    private int maxIter = 0;
    private double prevFunctValue = 0;
    private boolean verbose = false;
    
    private int curIter = 0;
    
    /**
     *
     @param eps
     */
    public ObjectiveDeltaStopStrategy(double eps) {
        this.minDelta = eps;
    }
    
    /**
     *
     @param eps
     @param maxIter
     */
    public ObjectiveDeltaStopStrategy(double eps, int maxIter) {
        this.minDelta = eps;
        this.maxIter = maxIter;
    }
    
    /**
     *
     @return
     */
    public ObjectiveDeltaStopStrategy beVerbose() {
        this.verbose = true;
        return this;
    }
    
    //x is coefficients for polynomial function, f_value, g

    /**
     *
     @param x
     @param fValue
     @param g
     @return
     */
    public boolean shouldContinueSearch(double[] x, double fValue,
        double[] g) {

        if (verbose) {
            System.out.println("iteration: " + curIter 
                + "   objective: " + fValue);
        }

        ++curIter;
        
        if (beenUsed) {
            // Check if we have hit the max allowable number of iterations.  (but only
            // check if _max_iter is enabled (i.e. not 0)).
            if (maxIter != 0 && curIter > maxIter) {
                return false;
            }

            // check if the function change was too small
            double diff = Math.abs(fValue - prevFunctValue);
            if (diff < minDelta) {
                return false;
            }
        }
        
        beenUsed = true;
        prevFunctValue = fValue;
        return true;
    }
}
