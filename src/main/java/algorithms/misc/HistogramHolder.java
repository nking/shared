package algorithms.misc;

import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;

/**
   first implemented in projects
     https://github.com/nking/two-point-correlation
     w/ Copyright (c) 2013-2015 Nichole King
     http://nking.github.io/two-point-correlation/
     using The MIT License (MIT)
     and
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class HistogramHolder {

    protected float[] xHist = null;
    protected int[] yHist = null;
    protected float[] yHistFloat = null;
    protected float[] yErrors = null;
    protected float[] xErrors = null;
    
    public float[] getHistArea(float maxXToUse, int nPartitions) {
        
        if (yHistFloat == null) {
            return null;
        }
        
        double[] area = new double[nPartitions];
        
        float binSize = maxXToUse/(float)nPartitions;
        
        // trapezoidal rule for area under the curve
        
        for (int i = 0; i < (xHist.length - 1); ++i) {
            
            float yTerm = yHistFloat[i + 1] + yHistFloat[i];
            float xLen = xHist[i + 1] - xHist[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            
            float x = xHist[i];
            
            int partition = (int)(x/binSize);
            
            if (partition > (nPartitions - 1)) {
                partition = nPartitions - 1;
            }
            
            area[partition] += (yTerm * xLen);
        }
                
        double sum = 0;
        for (int i = 0; i < nPartitions; ++i) {
            sum += area[i]; // area should be multiplied by 0.5, but that's not needed for normalization
        }
        
        float[] frac = new float[nPartitions];
        for (int i = 0; i < nPartitions; ++i) {
            frac[i] = (float)(area[i]/sum);
        }
        
        return frac;
    }

    /**
     * integrate the area of the histogram and the area of the restricted
     * range of the histogram (from x0 to x1, inclusive) and then return
     * the value of the restricted range over the total.
     * @param x0
     * @param x1
     * @return 
     */
    public float getHistAreaFractionOfTotal(float x0, float x1) {
        
        if (yHistFloat == null) {
            return 0;
        }
        
        float sumTot = 0;
        float sumR = 0;
        
        // trapezoidal rule for area under the curve
        
        for (int i = 0; i < (xHist.length - 1); ++i) {
            
            float yTerm = yHistFloat[i + 1] + yHistFloat[i];
            float xLen = xHist[i + 1] - xHist[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            
            float v = (yTerm * xLen);
            
            float x = xHist[i];
            
            sumTot += v;
            
            if ((x >= x0) && (x <= x1)) {
                sumR += v;
            }
        }
        
        float frac = sumR / sumTot;
        
        return frac;
    }
    
    /**
     * integrate the area under the curve of the histogram.
     * @return 
     */
    public float getHistArea() {
        
        if (yHistFloat == null) {
            return 0;
        }
        
        float sumTot = 0;
        
        // trapezoidal rule for area under the curve
        
        for (int i = 0; i < (xHist.length - 1); ++i) {
            
            float yTerm = yHistFloat[i + 1] + yHistFloat[i];
            float xLen = xHist[i + 1] - xHist[i];
            if (xLen < 0) {
                xLen *= -1;
            }
            
            float v = (yTerm * xLen);
                        
            sumTot += v;
        }
        
        sumTot *= 0.5;
        
        return sumTot;
    }

    public int calculateHalfYMaxIndexPastYMax() {

        if (yHistFloat == null) {
            return -1;
        }

        int yMaxIndex = MiscMath0.findYMaxIndex(yHistFloat);

        int halfMaxIndex = -1;
        float halfMax = yHistFloat[yMaxIndex]/2.0f;

        for (int i = yMaxIndex; i < yHistFloat.length; i++) {
            if (halfMax <= yHistFloat[i]) {
                halfMaxIndex = i;
            }
        }
        return halfMaxIndex;
    }
    
    public String plotHistogram(String label, 
        long outputFileNumber) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath0.findMin(yh);
        int yMaxIdx = MiscMath0.findYMaxIndex(yh);
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath0.findMin(xh);
        float xMax = MiscMath0.findMax(xh);        
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, 1.1f*yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileNumber);
    }
    
    public String plotHistogram(String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath0.findMin(yh);
        int yMaxIdx = MiscMath0.findYMaxIndex(yh);
        if (yMaxIdx == -1) {
            return null;
        }
        float yMax = yh[yMaxIdx];
        
        float xMin = MiscMath0.findMin(xh);
        float xMax = MiscMath0.findMax(xh);        
                
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, 1.1f*yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileSuffix);
    }
    
    public String plotHistogram(float xMin, float xMax, String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float yMin = MiscMath0.findMin(yh);
        int yMaxIdx = MiscMath0.findYMaxIndex(yh);
        if (yMaxIdx == -1) {
            return null;
        }
        float yMax = yh[yMaxIdx];
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, 1.1f*yMax,
            xh, yh, xh, yh, label);

        return plotter.writeFile(outputFileSuffix);
    }
    
    public String plotLogHistogram(String label, 
        String outputFileSuffix) throws IOException {
                
        float[] xh = xHist;
        float[] yh = yHistFloat;
        
        float[] yLogH = new float[yh.length];
        for (int i = 0; i < yh.length; ++i) {
            yLogH[i] = (float)Math.log(yh[i]/Math.log(10));
        }
        
        float yMin = MiscMath0.findMin(yLogH);
        int yMaxIdx = MiscMath0.findYMaxIndex(yLogH);
        float yMax = yLogH[yMaxIdx];
        
        float xMin = MiscMath0.findMin(xh);
        float xMax = MiscMath0.findMax(xh);
                        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();

        plotter.addPlot(
            xMin, xMax, yMin, 1.1f*yMax,
            xh, yLogH, xh, yLogH, label);

        return plotter.writeFile(outputFileSuffix);
    }

    /**
     * @return the xHist
     */
    public float[] getXHist() {
        return xHist;
    }

    /**
     * @return the yHist
     */
    public int[] getYHist() {
        return yHist;
    }

    /**
     * @return the yHistFloat
     */
    public float[] getYHistFloat() {
        return yHistFloat;
    }

    /**
     * @return the yErrors
     */
    public float[] getYErrors() {
        return yErrors;
    }

    /**
     * @return the xErrors
     */
    public float[] getXErrors() {
        return xErrors;
    }

    /**
     * @param xHist the xHist to set
     */
    public void setXHist(float[] xHist) {
        this.xHist = xHist;
    }

    /**
     * @param yHist the yHist to set
     */
    public void setYHist(int[] yHist) {
        this.yHist = yHist;
    }

    /**
     * @param yHistFloat the yHistFloat to set
     */
    public void setYHistFloat(float[] yHistFloat) {
        this.yHistFloat = yHistFloat;
    }
    
    /**
     * @param yHistInt the yHistFloat to set
     */
    public void setYHistFloat(int[] yHistInt) {
        this.yHistFloat = new float[yHistInt.length];
        for (int i = 0; i < yHistInt.length; ++i) {
            yHistFloat[i] = yHistInt[i];
        }
    }

    /**
     * @param yErrors the yErrors to set
     */
    public void setYErrors(float[] yErrors) {
        this.yErrors = yErrors;
    }

    /**
     * @param xErrors the xErrors to set
     */
    public void setXErrors(float[] xErrors) {
        this.xErrors = xErrors;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        
        sb.append("histogram=[");
        
        for (int i = 0; i < xHist.length; i++) {
            
            sb.append("(").append(xHist[i]).append(", ");
            
            if (yHist != null) {
                sb.append(yHist[i]);
            } else {
                sb.append(yHistFloat[i]);
            }
            sb.append(") ");
            
        }
        
        sb.append("]%n");
        
        if (xErrors != null) {
            
            sb.append("histogram errors=[");
            
            for (int i = 0; i < xErrors.length; i++) {

                sb.append("(").append(xErrors[i]).append(", ");

                sb.append(yErrors[i]).append(") ");
            }
            
            sb.append("]%n");
        }      
        
        return sb.toString();
    }
}
