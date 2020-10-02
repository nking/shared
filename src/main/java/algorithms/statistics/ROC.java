package algorithms.statistics;

import algorithms.misc.MiscSorter;
import algorithms.util.PairFloatArray;
import java.util.Arrays;

/**
 * calculates ROC plot points and area under the curve given a set of test
 * examples with f-scores, and positive or negative labels.
 * 
 * @author nichole
 */
public class ROC {
   
    /**
     * calculates for 2-classes, the ROC plot points and area under the curve 
     * given a set of test examples with f-scores, and positive or negative labels.
     * The algorithms are from “An introduction to ROC analysis” by Fawcett, 2006, 
     * Pattern Recognition Letters 27, 61–874
       https://www.math.ucdavis.edu/~saito/data/roc/fawcett-roc.pdf

       Runtime is O(n*log_2(n)).
       
     * @param fScores
     * @param labels
     * @return 
     */
    public static ROCResults calcAUCAndPoints(double[] fScores, boolean[] labels) {
        
        fScores = Arrays.copyOf(fScores, fScores.length);
        
        // sort by decreasing fScores
        int[] indexes = MiscSorter.mergeSortDecreasing(fScores);
        
        int[] nPCounts = countLabels(labels);
        float neg = nPCounts[0];
        float pos = nPCounts[1];
        
        double tol = 1.e-15;
        float x, y;
        float fp = 0;
        float tp = 0;
        float fpPrev = 0;
        float tpPrev = 0;
        
        PairFloatArray pts = new PairFloatArray(fScores.length + 1);
        
        double fPrev = Double.NEGATIVE_INFINITY;
        
        double auc = 0;
        
        double diff;
        int i = 0;
        while (i < fScores.length) {
            diff = Math.abs(fPrev - fScores[i]);
            if (diff > tol) {
                x = fp/neg;
                y = tp/pos;
                pts.add(x, y);
                
                auc += (trapezoidalArea(fp, fpPrev, tp, tpPrev));
                
                fPrev = fScores[i];
                fpPrev = fp;
                tpPrev = tp;
            }
            if (labels[indexes[i]]) {
                tp++;
            } else {
                fp++;
            }
            ++i;
        }
        x = fp / neg;
        y = tp / pos;
        pts.add(x, y);
        
        auc += (trapezoidalArea(fp, fpPrev, tp, tpPrev));
        auc /= (neg * pos);
        
        ROCResults roc = new ROCResults();
        roc.sortedFScores = fScores;
        roc.indexes = indexes;
        roc.pts = pts;
        roc.auc = auc;
        
        return roc;
    }
    
    /**
     Not currently implemented, but could be upon need.
     Following the Hand and Till 2001 summary in Section 9 of 
        “An introduction to ROC analysis” by Fawcett, 2006, 
     Pattern Recognition Letters 27, 61–874
       https://www.math.ucdavis.edu/~saito/data/roc/fawcett-roc.pdf
    */
    public static double multiClassAUCHandTill2001() {
        /*
        Hand and Till (2001) derive a multi-class generalization of the AUC 
        that is insensitive to class distribution and error costs.  
        The derivation is too detailed to summarize here, but it is based 
        upon the fact that the AUC is equivalent to the probability that the 
        classifier will rank a randomly chosen positive instance higher than 
        a randomly chosen negative instance. From this probabilistic form, 
        they derive a formulation that measures the unweighted pairwise 
        discriminability of classes.  n is the number of classes and AUC(ci,cj) 
        is the area under the two-class ROC curve involving classes ci and cj.  
        The summation is calculated over all pairs of distinct classes, 
        irrespective of order. There are |C|*(|C| - 1)/2 such pairs, so the 
        time complexity of their measure is O(|C|^2*n*log_2(n)).
        
        nClasses = |C|.
        
        Their measure, which they call M, is equivalent to:
            AUC_total= (1/(nClasses*(nClasses-1)) 
                         * summation_{ci,cj}( AUC(ci, cj) )
        */
        throw new UnsupportedOperationException("not implemented, but could be upon need");
    }
    
    private static double trapezoidalArea(double x1, double x2, double y1, 
        double y2) {
        
        double base = Math.abs(x1 - x2);
        double height = (y1 + y2) / 2.;
        
        return base * height;
    }

    private static int[] countLabels(boolean[] labels) {
        int[] npc = new int[2];
        for (int i = 0; i < labels.length; ++i) {
            if (labels[i]) {
                npc[1]++;
            } else {
                npc[0]++;
            }
        }
        return npc;
    }
    
    public static class ROCResults {
        double[] sortedFScores;
        int[] indexes;
        PairFloatArray pts;
        double auc;
    }
}
