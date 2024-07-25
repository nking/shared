package algorithms.signalProcessing;

import algorithms.misc.MiscMath0;

/**
  The code below is from project
     http://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
     
 * @author nichole
 */
public class Util {
        
    public static float[] createGaussian(float sigma, float mu) {
        
        float normalization = (float)(sigma * Math.sqrt(2.f * Math.PI));
        
        float hwi = estimateHWZI(sigma, 0.001f);
        if (hwi < 0) {
            hwi *= -1.f;
        }
        int halfWidthInPixels = (int)Math.ceil(hwi);
        
        int start = -1*halfWidthInPixels;
        int stopExcl = halfWidthInPixels + 1;
        
        float d, dsq;
       
        int nPoints = stopExcl - start;
        
        float[] yPoints = new float[nPoints];
        int count = 0;
        for (int i = start; i < stopExcl; i++) {
            
            float x = i;
            d = (x - mu);
            dsq = d*d;
            
            float y = (float) Math.exp(-1.f * dsq/(2.f * sigma * sigma));
            
            yPoints[count] = y/normalization;
                        
            count++;
        }
                                
        return yPoints;
    }
    
    public static float estimateHWZI(float sigma, float fractionMax) {
        
        /*
           exp( (-(x0 - mu)^2)/2o~^2) / exp( (-(xcenter - mu)^2)/2o~^2) = fractionMax
        
           where xcenter = mu is the center of the gaussian
        
           exp( (-(x0 - mu)^2)/2o~^2) / 1 = fractionMax
        
           (-(x0 - mu)^2)/2o~^2) = ln(fractionMax)
        
           -(x0 - mu)^2 = -1*(2o~^2 * ln(fractionMax))
           
           x0 - mu = math.sqrt(-1*(2o~^2 * ln(fractionMax)))
        
           x0 = mu + o~ * math.sqrt(-1*(2 * ln(fractionMax)))
        
        */

        float x0 = (float)(sigma * Math.sqrt(-2. * Math.log(fractionMax)));
                
        return x0;
    }

    /**
     * given a histogram, measure the FWHM.
     * For a gaussian normal distribution, FWHM = 2*sqrt(2*ln(2)) * standard deviation
     * @param x
     * @return
     */
    public static float measureFWHM(float[] x) {
                
        int xMaxIdx = MiscMath0.findYMaxIndex(x);
        float xMax = x[xMaxIdx];
                
        int xHalfMax0Idx = -1;
        int xHalfMax1Idx = -1;
        for (int i = xMaxIdx; i > -1; i--) {
            if (x[i] > (xMax / 2.)) {
                xHalfMax0Idx = i;
            } else {
                break;
            }
        }
        for (int i = xMaxIdx - 1; i < x.length; i++) {
            if (x[i] > (xMax / 2.)) {
                xHalfMax1Idx = i;
            } else {
                break;
            }
        }
        
        //System.out.println("max=" + xMax + "\n" +
        //   " " + xHalfMax0Idx + " : " + xMaxIdx + " : " + xHalfMax1Idx);
        
        if ((xHalfMax1Idx == -1) || (xHalfMax0Idx == -1)) {
            return -1;
        }
        return xHalfMax1Idx - xHalfMax0Idx;
    }
    
}
