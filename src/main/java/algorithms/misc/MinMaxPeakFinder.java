package algorithms.misc;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;

/**
 * find peaks in sequential data by finding the minima above  threshold and
 * the maxima and then the maxima which are a gain factor above either
 * adjacent minima.
 * The success of the method depends upon a reasonable lowThreshold and
 * gain.
 * 
 * NOTE: methods such as those in MedianSmooth could be used to find a
 * windowed mean and max (and hence, standard deviation, but they are
 * dependent upon the size of the window.  A peak that is a wide gradual
 * plateau and a narrow window might miss the peak.
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
 * 
 * @author nichole
 */
public class MinMaxPeakFinder {
    
    /**
     * find the indexes of the array values which are local maxima 
     * whose values are above a calculated lowThreshold and a 
     * fixed factor 2.5 above
     * one of the adjacent minima.
     @param values
     @return 
     */
    public int[] findPeaks(float[] values) {
                
        float mean3Percent = calculateMeanOfSmallest(values, 0.03f);
        
        return findPeaks(values, mean3Percent, 2.5f);
    }
    
    /**
     * copies and sorts values, then takes the fraction*values.length
     * first indexes and calculates the mean for them and returns it.
     * 
     @param values
     @param fraction fraction of sorted values to determine the mean of, i.e. range is [0, fraction*n)
     @return 
     */
    public float calculateMeanOfSmallest(float[] values, float fraction) {
        
        float[] a = Arrays.copyOf(values, values.length);
        
        Arrays.sort(a);
        
        float mean = 0;
        int end = Math.round(fraction * values.length);
        if (end == 0) {
            end = 1;
        }
        for (int i = 0; i < end; ++i) {
            mean += a[i];
        }
        mean /= (float)end;
        
        return mean;
    }
    
    /**
     * copies and sorts values, then takes the fraction*values.length
     * first indexes and calculates the mean for them and returns it.
     * 
     @param values
     @param fraction fraction of sorted values to determine the mean of, i.e. range is [0, fraction*n)
     @return 
     */
    public float calculateMeanOfSmallest(int[] values, float fraction) {
        
        int[] a = Arrays.copyOf(values, values.length);
        
        Arrays.sort(a);
        
        float mean = 0;
        int end = Math.round(fraction * values.length);
        if (end == 0) {
            end = 1;
        }
        for (int i = 0; i < end; ++i) {
            mean += a[i];
        }
        mean /= (float)end;
        
        return mean;
    }
    
    /**
     * find the indexes of the array values which are local maxima whose
     * values are above lowThreshold and a factor factorAboveMin above
     * one of the adjacent minima.
     * 
     @param values
     @param lowThreshold
     @param factorAboveMin
     @return an array with indexes of peak values.  if none were found an empty
     * array is returned.
     */
    public int[] findPeaks(float[] values, float lowThreshold, 
        float factorAboveMin) {
        
        int[] minMaxIdxs = findMinimaMaxima(values, lowThreshold);
        
        if (minMaxIdxs.length == 0) {
            return new int[0];
        } else if (minMaxIdxs.length == 1) {
            if (minMaxIdxs[0] < 0) {
                return new int[0];
            }
            return minMaxIdxs;
        } 
        
        TIntList peaks = new TIntArrayList(minMaxIdxs.length/2);

        // choose candidates from minMaxIndexes that are
        //     >= factorAboveMin for one adjacent minima
        for (int ii = 0; ii < minMaxIdxs.length; ii++) {

            int idx = minMaxIdxs[ii];

            if (idx > -1) {
                // this is a maxima

                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIdxs[iii];
                    if (idx2 < 0) {
                        float compare = values[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (values[idx] >= lowThreshold && 
                            values[idx] >= factorAboveMin * compare) {
                            
                            peaks.add(idx);
                            found = true;
                        }
                        break;
                    }
                }
                
                if (found && (ii > 0)) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIdxs.length; iii++) {
                    int idx2 = minMaxIdxs[iii];
                    if (idx2 < 0) {
                        float compare = values[-1*idx2];
                        if (compare < lowThreshold) {
                            // avoids divide by very small number sometimes
                            compare = lowThreshold;
                        }
                        if (values[idx] >= lowThreshold 
                            && values[idx] >= factorAboveMin * compare) {
                            
                            peaks.add(idx);
                        }
                        
                        break;
                    }
                }
            }
        }

        return peaks.toArray(new int[peaks.size()]);
    }
    
    /**
     * find the minima above lowThreshold and the maxima in values
     * and return their indexes.  the negative values are -1*index for a minima
     * while positive indexes are the indexes of maxima.
     * 
     @param values
     @param lowThreshold
     @return 
     */
    public int[] findMinimaMaxima(float[] values, float lowThreshold) {
        
        TIntList minMaxIdxs = new TIntArrayList();
        
        float lastK = values[0];
        boolean incr = true;
        for (int ii = 1; ii < values.length; ii++) {

            float currentK = values[ii];

            if ((currentK < lastK) && incr) {
                if (values[ii - 1] > lowThreshold) {
                    minMaxIdxs.add(ii - 1);
                }
                incr = false;
            } else if ((currentK > lastK) && !incr) {
                // values below outputLowThreshold[0] are handled by
                // callers.  TODO: redesign the caller and this method
                // to not need to understand peculiarities of the data.
                minMaxIdxs.add(-1*(ii - 1));
                incr = true;
            }

            lastK = currentK;
        }

        if (incr) {
            // add the last point
             minMaxIdxs.add(values.length - 1);
        }
        
        return minMaxIdxs.toArray(new int[minMaxIdxs.size()]);
    }
}
