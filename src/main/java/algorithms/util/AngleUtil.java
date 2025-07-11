package algorithms.util;

import algorithms.sort.MultiArrayMergeSort;
import algorithms.QuickSort;
import java.util.List;

/**
 *
 * @author nichole
 */
public class AngleUtil {

    /**
    calculates the difference of angles between pairs of points in set1 and
    set 2 ( angle of (diffX1, diffY1) - angle of (diffX2, diffY2).
    * diffX1, diffX2 are the difference
    * in x and y between 2 points in set 1, and diffX2, diffY2 are the
    * difference for the same matched point pair in set2.  The reference
    * frame is polar clockwise.
    * For example, (diffX1, diffY1) being (1, 4) and (diffX2, diffY2) being
    * (3.536, 2.12) leads to angles 284 minus 329 = -45.  Note that to
    * transform (diffX1, diffY1) to (diffX2, diffY2) one would apply
    * -1*result of this method.
    *
    * <pre>
    *       clockwise -->
              +Y        V
              270
           III | IV
       180 --------- 0   +X
           II  |  I
               90
    * </pre>
    * Note that the subtraction is from closest angles.
    * For example, if point 1 is in qI and point 2 is in qIV, 360 is added to
    * point 1's angle, then the equation is (2*PI + point1 angle) - (point 2 angle).

     * @param diffX1
     * @param diffY1
     * @param diffX2
     * @param diffY2
     * @return
     */
    public double subtract(double diffX1, double diffY1, double diffX2,
        double diffY2) {

        /*  clockwise -->
              +Y        V
              270
           III | IV
       180 --------- 0   +X
           II  |  I
               90
        */

        double theta1 = polarAngleCW(diffX1, diffY1);

        double theta2 = polarAngleCW(diffX2, diffY2);
        
        double twoPi = 2. * Math.PI;

        if (theta1 > theta2) {
            // add a phase to next value if it's closer to current with addition
            if ((theta1 - theta2) > Math.abs(theta1 - (theta2 + twoPi))) {
                theta2 += twoPi;
            }
        } else if (theta2 > theta1) {
            // add a phase to next value if it's closer to current with addition
            if ((theta2 - theta1) > Math.abs(theta2 - (theta1 + twoPi))) {
                theta1 += twoPi;
            }
        }

        return theta1 - theta2;
    }

    /**
    calculates the polar theta in radians given x and y w.r.t. origin.  theta increases
    * in value in a counter clockwise directionCCW (CCW).  range of returned
    * values is 0 to 2*pi.

     * @param x
     * @param y
     * @return angle in range [0, 2*Math.PI)
     */
    public static double polarAngleCCW(double x, double y) {

        /*
                  +Y
                 90
        QII       |       QI
                  |
                  |
     180-------------------- +X  0, 360
                  |
                  |
         QIII     |       QIV
                 270


        Math.atan2 angles are:
                 90
           135    |    45
                  |
        180 ---------------  0
                  |
          -135    |   -45
                 -90
        so, for d < 0, need 360+d
        */

        double theta = Math.atan2(y, x);

        if (theta < 0) {
            theta += 2. * Math.PI;
        }

        return theta;
    }

    /**
    calculates the polar theta in radians given x and y w.r.t. origin.  theta increases
    * in value in a clockwise directionCCW (CW).

     * @param x
     * @param y
     * @return polar angle in radians
     */
    public static double polarAngleCW(double x, double y) {

        /*
                  +Y
                 270
        QIII      |       QIV
                  |
                  |
     180-------------------- +X  0, 360
                  |
                  |
         QII      |       QI
                  90

        Math.atan2 angles are:
                 90
           135    |    45
                  |
        180 ---------------  0
                  |
          -135    |   -45
                 -90
        so need -1*d, then for d < 0, need 360+d
        */

        double theta = -1 * Math.atan2(y, x);

        if (theta < 0) {
            theta += 2. * Math.PI;
        }

        return theta;
    }

    /**
     * get the difference in angles, making correction for quadrants.
     * @param rotDegrees0
     * @param rotDegrees1
     * @return 
     */
    public static float getAngleDifference(float rotDegrees0, float rotDegrees1) {

        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */

        double twoPi = 360;

        if (rotDegrees0 > rotDegrees1) {
            // add a phase to next value if it's closer to current with addition
            if ((rotDegrees0 - rotDegrees1) > Math.abs(rotDegrees0 - (rotDegrees1 + twoPi))) {
                rotDegrees1 += twoPi;
            }
        } else if (rotDegrees1 > rotDegrees0) {
            // add a phase to next value if it's closer to current with addition
            if ((rotDegrees1 - rotDegrees0) > Math.abs(rotDegrees1 - (rotDegrees0 + twoPi))) {
                rotDegrees0 += twoPi;
            }
        }

        return rotDegrees0 - rotDegrees1;
    }

    public static float getAngleAverageInDegrees(float rotDegrees0, float rotDegrees1) {

        double angleAvg = getAngleAverage(rotDegrees0, rotDegrees1, false);

        return (float) angleAvg;
    }

    public static double getAngleAverageInRadians(double rotation0, double rotation1) {

        double angleAvg = getAngleAverage(rotation0, rotation1, true);

        return angleAvg;
    }

    /**
     * given twoPi in degrees or in radians, return the angle average.
     * @param rot0
     * @param rot1
     * @param useRadians
     * @return
     */
    protected static double getAngleAverage(double rot0, double rot1, boolean useRadians) {

         /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
        double angleAvg = calcAngleAddition(rot0, rot1, useRadians)/2.;

        return angleAvg;
    }

    /**
     * given twoPi in degrees or in radians, return the angle sum corrected to
     * the larger angle frame, e.g. 0 + 350 = 710.
     * @param theta1
     * @param theta2
     * @param useRadians
     * @return
     */
    public static double calcAngleAddition(double theta1, double theta2,
        boolean useRadians) {

        if (theta1 < 0 || theta2 < 0) {
            throw new IllegalArgumentException(
                "theta1 and theta2 cannot be negative numbers");
        }

        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */

        double twoPi = useRadians ? 2. * Math.PI : 360;

        if (theta1 > theta2) {
            // add a phase to next value if it's closer to current with addition
            if ((theta1 - theta2) > Math.abs(theta1 - (theta2 + twoPi))) {
                theta2 += twoPi;
            }
        } else if (theta2 > theta1) {
            // add a phase to next value if it's closer to current with addition
            if ((theta2 - theta1) > Math.abs(theta2 - (theta1 + twoPi))) {
                theta1 += twoPi;
            }
        }

        return theta1 + theta2;
    }

    /**
     * given angles in degrees or in radians, return the angle sum corrected to
     * the larger angle frame, e.g. 0 + 350 = 710.  It always chooses the
     * smallest difference in quadrant space
     * between the two angles (adding a cycle if needed) before adding.
     * @param theta1
     * @param theta2
     * @param useRadians
     * @param outputQuadrantCorrected quadrant corrected values of
     * [rotationInDegrees0, rotationInDegrees1]
     * @return
     */
    public static double calcAngleAddition(double theta1, double theta2,
        boolean useRadians, double[] outputQuadrantCorrected) {

        if (theta1 < 0 || theta2 < 0) {
            throw new IllegalArgumentException(
                "theta1 and theta2 cannot be negative numbers");
        }

        /*
                  270
               III | IV
          180  --------- 0
               II  |  I
                   90
        */
        double twoPi = useRadians ? 2. * Math.PI : 360;
        
        if (theta1 > theta2) {
            // add a phase to next value if it's closer to current with addition
            if ((theta1 - theta2) > Math.abs(theta1 - (theta2 + twoPi))) {
                theta2 += twoPi;
            }
        } else if (theta2 > theta1) {
            // add a phase to next value if it's closer to current with addition
            if ((theta2 - theta1) > Math.abs(theta2 - (theta1 + twoPi))) {
                theta1 += twoPi;
            }
        }
        
        outputQuadrantCorrected[0] = theta1;
        outputQuadrantCorrected[1] = theta2;
        
        return theta1 + theta2;
    }

    /**
     * calculate the angular average of rotDegrees0 and rotDegrees1 with
     * corrections for quadrants if necessary.  For example, 0 averaged with
     * 350 is 355.  Any changes in the given angles for quadrant math
     * are seen in the populated outputQuadrantCorrected.
     * @param outputQuadrantCorrected quadrant corrected values of
     * [rotationInDegrees0, rotationInDegrees1]
     * @return
    */
    protected static double getAngleAverageInDegrees(
        double rotationInDegrees0, double rotationInDegrees1,
        double[] outputQuadrantCorrected) {

        boolean useRadians = false;

        double sum = calcAngleAddition(rotationInDegrees0,
            rotationInDegrees1, useRadians, outputQuadrantCorrected);

        return sum/2.;
    }

    /**
     calculate the average of the angles using quadrant corrections (note,
     * angles is modified by descending sort).
     Runtime complexity is O(N*lg2N) for sort plus O(N).
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.

     * @param angles angles in units of degrees.  Note that this array
     * is modified by use here so pass a copy if it should not be.
     * @param useRadians
     * @return
     */
    public static double calculateAverageWithQuadrantCorrections(
        double[] angles, boolean useRadians) {

        //TODO: because 360 is usually the max value, if N is > 80, can
        // use CountingSort here for better performance.
        // CountingSort is O(maxValue) while MergeSort is O(N*lg_2(N))
        QuickSort.descendingSort(angles);

        double twoPI = useRadians ? 2. * Math.PI : 360;
        
        // visit in order of large values to low to store quadrant corrections
        double avg = 0;
        for (int i = 0; i < angles.length; ++i) {

            double angle1 = angles[i];

            if (i == 0) {

                avg = angle1;

            } else {

                double angle0 = avg;

                // add a phase to next value if it's closer to current with addition
                if ((angle0 - angle1) > Math.abs(angle0 - (angle1 + twoPI))) {
                    angle1 += twoPI;
                    if (useRadians) {
                        angles[i] = angle1;
                    } else {
                        angles[i] = Math.round(angle1);
                    }
                }

                avg = (angle0 + angle1)/2.f;
            }
        }
        
        avg = 0;
        for (double a : angles) {
            avg += a;
        }
        avg /= (double)angles.length;
        
        if (avg > twoPI) {
            avg -= twoPI;
        }

        return avg;
    }

    /**
     Runtime complexity is O(N*lg2N) for sort plus O(N).
     * Note, may change to use CountingSort for N > 80 and max(angles) ~ 360
     * in the future.

     */
    public static void correctForQuadrants(List<Double> thetas,
        boolean useRadians) {

        if (thetas == null || thetas.size() < 2) {
            return;
        }

        double[] angles = new double[thetas.size()];
        int[] indexes = new int[thetas.size()];

        for (int i = 0; i < thetas.size(); ++i) {
            angles[i] = thetas.get(i).doubleValue();
            indexes[i] = i;
        }

        MultiArrayMergeSort.sortByDecr(angles, indexes);

        double twoPi = useRadians ? 2. * Math.PI : 360;
        
        double avg = 0;
        for (int i = 0; i < angles.length; ++i) {
            double angle1 = angles[i];
            if (i == 0) {
                avg = angle1;
            } else {
                double angle0 = avg;
                if ((angle0 - angle1) > Math.abs(angle0 - (angle1 + twoPi))) {
                    angle1 += twoPi;
                    angles[i] = Math.round(angle1);
                }
                avg = (angle0 + angle1)/2.f;
            }
        }

        for (int i = 0; i < angles.length; ++i) {
            int idx = indexes[i];
            double v = angles[i];
            thetas.set(idx, v);
        }
    }

}
