package thirdparty.fpetitjean.dba;

/**
 * Dynamic Time Warping (DTW) Barycenter Averaging
 * is used to globally aligning temporal sequences of different speeds and measure their similarity.
 * 
 * code from publication
 * <pre>
 * "A global averaging method for dynamic time warping, with applications to clustering"
    Petitjean, Ketterlin, & Gancarski 
    Pattern recognition, 2011 - Elsevier
 * </pre>
 * 
 * *****************************************************************************
 * Copyright (C) 2018 Francois PETITJEAN
 *
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * source code is available at:
 * https://github.com/fpetitjean/DBA/blob/master/DBA.java
 * 
 * license is available at:
 * https://github.com/fpetitjean/DBA/blob/master/LICENSE
 *****************************************************************************
 */
import algorithms.util.FormatArray;
import static java.lang.Math.sqrt;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * This toy class show the use of DBA.
 *
 * @author Francois Petitjean
 */
public class DynamicTimeWarpingBarycenterAveraging {

    static final long serialVersionUID = 1L;

    private final static int NIL = -1;
    private final static int DIAGONAL = 0;
    private final static int LEFT = 1;
    private final static int UP = 2;

    // encoding that 0 => (-1,-1), 1 => (0,-1), 2=> (-1,0)
    private final static int[] moveI = {-1, 0, -1};
    private final static int[] moveJ = {-1, -1, 0};

    /**
     * Performs the DBA averaging by first finding the median over a sample,
     * then doing n iterations of the update.
     *
     @param sequences set of sequences to average
     @param nIterations the number of iterations to run it for (default 15)
     @return 
     */
    public static double[] performDBA(double[][] sequences, int nIterations) {

        int maxLength = 0;
        for (int i = 0; i < sequences.length; i++) {
            // java allows jagged arrays
            maxLength = Math.max(maxLength, sequences[i].length);
        }
        double[][] costMatrix = new double[maxLength][maxLength];
        int[][] pathMatrix = new int[maxLength][maxLength];
        
        //NLK: since calc medoid is at most 50*(sequence.length * sequence.length * sequences[i].length)
        //     which is at most 50*(N^3), could replace this with geometric median using
        //     DTW for the distance metric.
        //     see algorithms.optimization.GeometricMedian
        int medoidIndex = approximateMedoidIndex(sequences, costMatrix);
        
        double[] center = Arrays.copyOf(sequences[medoidIndex], sequences[medoidIndex].length);

        for (int i = 0; i < nIterations; i++) {
            System.out.println("center=" + Arrays.toString(center));
            center = dBAUpdate(center, sequences, costMatrix, pathMatrix);
        }
        return center;
    }

    /**
     * Performs the DBA averaging by first finding the median over a sample,
     * then doing n iterations of the update
     *
     @param sequences set of sequences to average
     @return 
     */
    public static double[] performDBA(double[][] sequences) {
        return performDBA(sequences, 15);
    }

    /**
     * 
     @param sequences input array of rows of sequences
     @param mat output array to hold cost matrices
     @return 
     */
    private static int approximateMedoidIndex(double[][] sequences, double[][] mat) {
        /*
         * we are finding the medoid, as this can take a bit of time, if
         * there is more than 50 time series, we sample 50 as possible
         * medoid candidates
         */
        ArrayList<Integer> allIndices = new ArrayList<>();
        for (int i = 0; i < sequences.length; i++) {
            allIndices.add(i);
        }
        Collections.shuffle(allIndices);
        ArrayList<Integer> medianIndices = new ArrayList<>();
        for (int i = 0; i < sequences.length && i < 50; i++) {
            medianIndices.add(allIndices.get(i));
        }

        int indexMedoid = -1;
        double lowestSoS = Double.MAX_VALUE;

        //NLK: runtime complexity is at most 50*(sequence.length * sequence.length * sequences[i].length)
        for (int medianCandidateIndex : medianIndices) {
            double[] possibleMedoid = sequences[medianCandidateIndex];
            // NLK: runtime complexity is sequence.length * sequence.length * sequences[i].length
            double tmpSoS = sumOfSquares(possibleMedoid, sequences, mat);
            if (tmpSoS < lowestSoS) {
                indexMedoid = medianCandidateIndex;
                lowestSoS = tmpSoS;
            }
        }
        return indexMedoid;
    }

    private static double sumOfSquares(double[] sequence, double[][] sequences, double[][] mat) {
        double sos = 0.0;
        // NLK: runtime complexity is sequence.length * sequence.length * sequences[i].length
        for (int i = 0; i < sequences.length; i++) {
            // NLK: runtime complexity is S.length * T.length
            double dist = DTW(sequence, sequences[i], mat);
            sos += dist * dist;
        }
        return sos;
    }

    /**
     * calculate dynamic time warping similarity measure.
     @param S
     @param T
     @param costMatrix
     @return 
     */
    public static double DTW(double[] S, double[] T, double[][] costMatrix) {
        int i, j;
        costMatrix[0][0] = squaredDistance(S[0], T[0]);
        for (i = 1; i < S.length; i++) {
            costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(S[i], T[0]);
        }
        for (j = 1; j < T.length; j++) {
            costMatrix[0][j] = costMatrix[0][j - 1] + squaredDistance(S[0], T[j]);
        }
        for (i = 1; i < S.length; i++) {
            for (j = 1; j < T.length; j++) {
                costMatrix[i][j] = Min3(costMatrix[i - 1][j - 1], costMatrix[i][j - 1], costMatrix[i - 1][j])
                        + squaredDistance(S[i], T[j]);
            }
        }

        return sqrt(costMatrix[S.length - 1][T.length - 1]);
    }

    private static double[] dBAUpdate(double[] C, double[][] sequences, double[][] costMatrix, int[][] pathMatrix) {
        double[] updatedMean = new double[C.length];
        int[] nElementsForMean = new int[C.length];

        int i, j, move;
        double res = 0.0;
        int centerLength = C.length;
        int seqLength;

        for (double[] T : sequences) {
            seqLength = T.length;

            costMatrix[0][0] = squaredDistance(C[0], T[0]);
            pathMatrix[0][0] = DynamicTimeWarpingBarycenterAveraging.NIL;

            for (i = 1; i < centerLength; i++) {
                costMatrix[i][0] = costMatrix[i - 1][0] + squaredDistance(C[i], T[0]);
                pathMatrix[i][0] = DynamicTimeWarpingBarycenterAveraging.UP;
            }
            for (j = 1; j < seqLength; j++) {
                costMatrix[0][j] = costMatrix[0][j - 1] + squaredDistance(T[j], C[0]);
                pathMatrix[0][j] = DynamicTimeWarpingBarycenterAveraging.LEFT;
            }

            for (i = 1; i < centerLength; i++) {
                for (j = 1; j < seqLength; j++) {
                    double diag = costMatrix[i - 1][j - 1], left = costMatrix[i][j - 1], top = costMatrix[i - 1][j];
                    if (diag <= left) {
                        if (diag <= top) {
                            res = diag;
                            move = DIAGONAL;
                        } else {
                            res = top;
                            move = UP;
                        }
                    } else {
                        if (left <= top) {
                            res = left;
                            move = LEFT;
                        } else {
                            res = top;
                            move = UP;
                        }
                    }

                    pathMatrix[i][j] = move;
                    res = costMatrix[i + moveI[move]][j + moveJ[move]];
                    costMatrix[i][j] = res + squaredDistance(C[i], T[j]);
                }
            }

            i = centerLength - 1;
            j = seqLength - 1;
            
            //NLK: if making a version of this algorithm to use geometric median instead of geometric mean
            //  used here, would need to change this too (though might revise the iteration over each sequence too).
            while (pathMatrix[i][j] != DynamicTimeWarpingBarycenterAveraging.NIL) {
                updatedMean[i] += T[j];
                nElementsForMean[i]++;
                move = pathMatrix[i][j];
                i += moveI[move];
                j += moveJ[move];
            }
            /*assert (i != 0 || j != 0);
            updatedMean[i] += T[j];
            nElementsForMean[i]++;*/
            if (i != 0 || j != 0) {
                updatedMean[i] += T[j];
                nElementsForMean[i]++;
            }
        }

        for (int t = 0; t < centerLength; t++) {
            updatedMean[t] /= nElementsForMean[t];
        }

        return updatedMean;

    }

    private static double Min3(final double a, final double b, final double c) {
        if (a < b) {
            if (a < c) {
                return a;
            } else {
                return c;
            }
        } else {
            if (b < c) {
                return b;
            } else {
                return c;
            }
        }
    }

    private static int ArgMin3(final double a, final double b, final double c) {
        if (a < b) {
            if (a < c) {
                return 0;
            } else {
                return 2;
            }
        } else {
            if (b < c) {
                return 1;
            } else {
                return 2;
            }
        }
    }

    private static double squaredDistance(double a, double b) {
        double diff = a - b;
        return diff * diff;
    }

}
