/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package algorithms.imageProcessing;

import algorithms.FixedSizeSortedVector;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import gnu.trove.list.TIntList;
import java.util.Arrays;

/**
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
   * 
 * @author nichole
 */
public class Filters {
    
    /**
     * @author nichole
     * @param img
     * @param size
     * @return
     */
    public static float[][] maximumFilter(float[][] img, int size) {

        int nRows = img.length;
        int nCols = img[0].length;

        // return_value = out
        float[][] out = new float[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            out[i] = new float[nCols];
        }

        // have adapted median window algorithm for this:
        StatsInSlidingWindow maxWindow = new StatsInSlidingWindow();
        maxWindow.calculateMaximum(img, out, size, size);

        return out;
    }
    
       /**
     Find peaks in an image as coordinate list
     Peaks are the local maxima in a region of `2 * min_distance + 1`
     (i.e. peaks are separated by at least `min_distance`).
     If peaks are flat (i.e. multiple adjacent pixels have identical
     intensities), the coordinates of all such pixels are returned.
     If both `threshold_abs` and `threshold_rel` are provided, the maximum
     of the two is chosen as the minimum intensity threshold of peaks.

      adapted from
     https://github.com/scikit-image/scikit-image/blob/92a38515ac7222aab5e606f9de46caf5f503a7bd/skimage/feature/peak.py

     The implementation below is adapted from the scipy implementation which has
     * the following copyright:

     https://github.com/scikit-image/scikit-image/blob/master/LICENSE.txt

    -- begin scipy, skimage copyright ---
    Unless otherwise specified by LICENSE.txt files in individual
    directories, all code is

    Copyright (C) 2011, the scikit-image team
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in
        the documentation and/or other materials provided with the
        distribution.
     3. Neither the name of skimage nor the names of its contributors may be
        used to endorse or promote products derived from this software without
        specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
    IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
    INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    -- end scipy, skimage copyright ---
    */
    public void peakLocalMax(float[][] img, int minDistance,
        float thresholdRel,
        TIntList outputKeypoints0, TIntList outputKeypoints1) {

        int excludeBorder = minDistance;
        int numPeaks = Integer.MAX_VALUE;
        //int numPeaksPerLabel = Integer.MAX_VALUE;
        
        /*
        The peak local maximum function returns the coordinates of local peaks
        (maxima) in an image. A maximum filter is used for finding local maxima.
        This operation dilates the original image. After comparison of the dilated
        and original image, this function returns the coordinates or a mask of the
        peaks where the dilated image equals the original image.
        */

        int nRows = img.length;
        int nCols = img[0].length;

        //# Non maximum filter
        int size = 2 * minDistance + 1;
        float[][] imageMax = maximumFilter(img, size);
        assert(nRows == imageMax.length);
        assert(nCols == imageMax[0].length);
        //mask = image == image_max

        //debugPrint("before shift imageMax=", imageMax);

        // a fudge to match results of scipy which must store same windows at
        // locations shifted by minDistance or so in x and y from the
        // beginning of the sliding window
        if (minDistance != 0) {
            applyShift(imageMax, minDistance, nRows, nCols);
        }
        
        // 1's where same, else 0's
        int[][] mask = new int[nRows][nCols];
        for (int i = 0; i < nRows; ++i) {
            mask[i] = new int[nCols];
            for (int j = 0; j < nCols; ++j) {
                if (img[i][j] == imageMax[i][j]) {
                    mask[i][j] = 1;
                }
            }
        }
        
        // exclude border
        for (int i = 0; i < nRows; ++i) {
            if ((i < excludeBorder) || (i > (nRows - 1 - excludeBorder))){
                Arrays.fill(mask[i], 0);
            } else {
                Arrays.fill(mask[i], 0, excludeBorder, 0);
                Arrays.fill(mask[i], nCols - excludeBorder, nCols, 0);
            }
        }


        // find top peak candidates above a threshold.
        // TODO: should this use mask so excluding borders?
        float thresholdAbs = MiscMath0.findMin(img);
        float thresholdMax = thresholdRel * MiscMath0.findMax(img);
        thresholdAbs = Math.max(thresholdAbs, thresholdMax);
        
        // mask &= image > 0.1
        for (int i = 0; i < nRows; ++i) {
            for (int j = 0; j < nCols; ++j) {
                if (imageMax[i][j] > thresholdAbs) {
                    mask[i][j] &= 1;
                } else {
                    mask[i][j] = 0;
                }
            }
        }
        
        //TODO: should num_peaks be this.nKeypoints?  re-read paper...
        if (numPeaks == Integer.MAX_VALUE) {
            // find non-zero pixels in mask
            float[] values = new float[nRows * nCols];
            int[] pixIdxs = new int[values.length];
            int count = 0;
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        values[count] = img[i][j];
                        //(row * width) + col
                        pixIdxs[count] = (j * nRows) + i;
                        count++;
                    }
                }
            }
            values = Arrays.copyOf(values, count);
            pixIdxs = Arrays.copyOf(pixIdxs, count);
            MiscSorter.sortByDecr(values, pixIdxs);
            
            for (int i = 0; i < values.length; ++i) {
                int pixIdx = pixIdxs[i];
                int jj = pixIdx/nRows;
                int ii = pixIdx - (jj * nRows);
                outputKeypoints0.add(ii);
                outputKeypoints1.add(jj);
            }
        } else {
            //need to sort to keep top numPeaks
            FixedSizeSortedVector<Pix> vec = new
                FixedSizeSortedVector<Pix>(numPeaks, Pix.class);
            for (int i = 0; i < mask.length; ++i) {
                for (int j = 0; j < mask[i].length; ++j) {
                    if (mask[i][j] > 0.f) {
                        Pix pix = new Pix(i, j, Float.valueOf(img[i][j]));
                        vec.add(pix);
                    }
                }
            }
            for (int i = 0; i < vec.getNumberOfItems(); ++i) {
                Pix pix = vec.getArray()[i];
                outputKeypoints0.add(pix.i);
                outputKeypoints1.add(pix.j);
            }
        }
    }
    
    private static class Pix implements Comparable<Pix> {

        public final int i;
        public final int j;
        public final Float value;
        public Pix(int i, int j, Float v) {
            this.i = i;
            this.j = j;
            this.value = v;
        }
        @Override
        public int compareTo(Pix other) {
            // changed for a descending sort
            return other.value.compareTo(this.value);
        }

    }
   
    private static void applyShift(float[][] imageMax, int minDistance, int nRows,
        int nCols) {

        for (int i = 0; i < nRows; ++i) {
            System.arraycopy(imageMax[i], 0, imageMax[i], minDistance,
                nCols - minDistance);
            for (int j = 0; j < minDistance; ++j) {
                imageMax[i][j] = 0;
            }
            for (int j = (nCols - minDistance) - 1; j < nCols; ++j) {
                imageMax[i][j] = 0;
            }
        }
        for (int j = 0; j < nCols; ++j) {
            for (int i = (nRows - minDistance) - 1; i >= minDistance; --i) {
                imageMax[i][j] = imageMax[i - minDistance][j];
            }
            for (int i = 0; i < minDistance; ++i) {
                imageMax[i][j] = 0;
            }
            for (int i = (nRows - minDistance) - 1; i < nRows; ++i) {
                imageMax[i][j] = 0;
            }
        }
    }
}
