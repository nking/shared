package algorithms.compGeometry;

import algorithms.compGeometry.convexHull.GrahamScanPairInt;
import algorithms.compGeometry.convexHull.GrahamScanTooFewPointsException;
import algorithms.util.AngleUtil;
import algorithms.util.PairIntWithIndex;
import algorithms.util.PairIntArray;
import algorithms.util.PairFloatArray;
import algorithms.util.PairInt;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class MiscellaneousCurveHelper {

    private Logger log = Logger.getLogger(this.getClass().getName());

    // choosing a minimum size empirically from looking at edges in tests
    private static int minLedgeWidth = 4;

    protected static final int[] eightNeighborsX =
        new int[]{-1, -1,  0,  1, 1, 1, 0, -1};
    protected static final int[] eightNeighborsY =
        new int[]{ 0, -1, -1, -1, 0, 1, 1,  1};

    /**
     * determine whether the closed curve points are ordered in a counter clockwise
     * manner 
     * calculating the cross product between adjacent edges in sequence around
     * the polygon to determine if there are fewer that are positive (CCW)
     * than negative (CW).
     * The given closedCurve cannot intersect itself or have holes in it.
     * NOTE: the answer returns true if the points are ordered in CW manner, but
     * if one needs the answer w.r.t. viewing an image which has y increasing
     * downward, need the opposite of the return here.
     *
     * @param closedCurve
     * @return
     */
    public boolean curveIsOrderedClockwise(PairIntArray closedCurve) {

        if (closedCurve.getN() < 2) {
            return false;
        }

        int nNeg = 0;
        int n = closedCurve.getN();

        for (int i = 0; i < n; i++) {

            long xm1, ym1, x, y, xp1, yp1;

            if (i == 0) {
                xm1 = closedCurve.getX(closedCurve.getN() - 1);
                ym1 = closedCurve.getY(closedCurve.getN() - 1);
                xp1 = closedCurve.getX(i + 1);
                yp1 = closedCurve.getY(i + 1);
            } else if (i == (closedCurve.getN() - 1)) {
                xm1 = closedCurve.getX(i - 1);
                ym1 = closedCurve.getY(i - 1);
                xp1 = closedCurve.getX(0);
                yp1 = closedCurve.getY(0);
            } else {
                xm1 = closedCurve.getX(i - 1);
                ym1 = closedCurve.getY(i - 1);
                xp1 = closedCurve.getX(i + 1);
                yp1 = closedCurve.getY(i + 1);
            }
            x = closedCurve.getX(i);
            y = closedCurve.getY(i);

            long dxmxm1 = (x - xm1);
            long dymym1 = (y - ym1);
            long dxp1mx = (xp1 - x);
            long dyp1my = (yp1 - y);

            //(xi - xi-1) * (yi+1 - yi) - (yi - yi-1) * (xi+1 - xi)
            long crossProduct = (dxmxm1 * dyp1my) - (dymym1 * dxp1mx);

            if (crossProduct < 0) {
                // clockwise when crossProduct is negative
                nNeg++;
            }
        }

        int nPos = n - nNeg;//n - 2 - nNeg;

        //log.info(closedCurve.toString());
        //log.info("n=" + n + " nNegative=" + nNeg + " nPositive=" + nPos);

        return ((n > 2) && (nNeg >= nPos)) || (nNeg > nPos);
     }

    /**
     * determine whether the closed curve points are ordered in a counter clockwise
     * manner by first computing the convex hull then
     * calculating the cross product between adjacent edges in sequence around
     * the polygon to determine if there are fewer that are positive (CCW)
     * than negative (CW).
     * The given closedCurve cannot intersect itself or have holes in it.
     * NOTE: the answer returns true if the points are ordered in CW manner, but
     * if one needs the answer w.r.t. viewing an image which has y increasing
     * downward, need the opposite of the return here.
     *
     * @param closedCurve
     * @return
     */
    public boolean curveIsOrderedClockwise2(PairIntArray closedCurve) {

        if (closedCurve.getN() < 2) {
            return false;
        } else if (closedCurve.getN() < 4) {
            return curveIsOrderedClockwise(closedCurve);
        }
        
        int n = closedCurve.getN();
        
        PairIntWithIndex[] p = new PairIntWithIndex[n];
        for (int i = 0; i < n; ++i) {
            p[i] = new PairIntWithIndex(closedCurve.getX(i), closedCurve.getY(i),  i);
        }
        
        GrahamScanPairInt<PairIntWithIndex> scan = new GrahamScanPairInt<PairIntWithIndex>();
        try {
            scan.computeHull(p);
            
            // hull returns points in clockwise order
            
            n = scan.getHull().size() - 1;
            //PairIntArray hull = new PairIntArray(n);
            //List<Integer> hullCurveIndexes = new ArrayList<Integer>();
            //int[] deltaIndexes = new int[n];
            
            // nPos or nNeg might be 1 and then other n-2 if there is wrap-around
            int nNeg = 0;
            int nPos = 0;
            for (int i = 0; i < n; ++i) {
                
                PairIntWithIndex p0 = scan.getHull().get(i);
                
                //hull.add(Math.round(p0.getX()), Math.round(p0.getY()));
                //hullCurveIndexes.add(Integer.valueOf(p0.getPixIndex()));
                
                // for CW input, expect these to be + numbers
                int deltaIndex = scan.getHull().get(i + 1).getPixIndex() - p0.getPixIndex();
                if (deltaIndex > 0) {
                    nPos++;
                } else {
                    nNeg++;
                }
            }
            
            //boolean isCW = curveIsOrderedClockwise(hull);
            //assert(isCW);
            
            if (nPos > nNeg) {
                return true;
            }
            
            return false;
            
        } catch (GrahamScanTooFewPointsException ex) {
            return curveIsOrderedClockwise(closedCurve);
        }
    }

    /**
         *           Y
         *          90
         *     135   |    +45
         *           |
         *   180------------ 0   X
         *           |
         *    225    |   315
         *          270
         * 
         *  -45    90    45          y/x
                -  |  +
            0 -----|----- 0
                +  |  -
            45    90    -45
        
     * @param angle360
     * @return 
     */
    private int convert360To45RefFrame(int angle360) {
        
        if ((angle360 < 23) || (angle360 > 337)) {
            return 0;
        } else if ((angle360 > 157) && (angle360 < 203)) {
            return 0;
        } else  if (((angle360 > 22) && (angle360 < 68)) || 
            ((angle360 > 202) && (angle360 < 248))) {
            // in range of +45 or +225
            return 45;
        } else if (((angle360 > 67) && (angle360 < 113)) || 
            ((angle360 > 247) && (angle360 < 293))) {
            // in range of +90 or +270
            return 90;
        } else { //if (((t > 112) && (t < 158)) || ((t > 292) && (t < 338))) {
            // in range of +135 or +315
            return -45;
        }
    }
    
    /**
     * return true if correlation shows that the 2 curves are adjacent
     * to one another.  Note that the method needs the points within the
     * curves to be ordered in a similar manner and for the endpoints of the
     * curves to be accurate.  If a point in the middle of the curve is
     * the first or last point, it may prevent comparison of it with another
     * edge's endpoints.
     *
     * @param curve0
     * @param curve1
     * @param correlationOffset offset of where the shorter curve starts
     *  with respect to the longer.  For example, an offset of -2 means that
     * the first 2 points in the shorter curve are outside of the longer curve,
     * but the next point in the longer curve is adjacent to the shorter.
     * Another example: if offset is +2, the first pixel in the shorter curve
     * is adjacent to the third pixel in the longer curve.  NOTE: the offset
     * is only useful if this method returns true;
     * @return
     */
    protected boolean correlation(PairIntArray curve0, PairIntArray curve1,
        int[] correlationOffset) {

        correlationOffset[0] = Integer.MAX_VALUE;

        //TODO: look at string matching algorithms to explore improvements here

        PairIntArray shorter, longer;
        if (curve0.getN() <= curve1.getN()) {
            shorter = curve0;
            longer = curve1;
        } else {
            shorter = curve1;
            longer = curve0;
        }

        /*
        len0 = 5; len1 = 11;
         #####
             +++++++++++
          #####
             +++++++++++
           #####
             +++++++++++
            #####
             +++++++++++
             #####
             +++++++++++

             #####
             +++++++++++

                       #####
             +++++++++++

        ccs = sqrt(sumsqdiff)/nOverlapping if nOverlapping > 0.

        if (css <= 1 pix * nOverlapping) {
            store as a possible adjacent curve
        }
        compare possible adjacent curves for the smallest css, and store that
        offset in correlationOffset and return true, else false
        */

        double cSSMin = Double.MAX_VALUE;
        int cSSMinOffset = Integer.MAX_VALUE;
        int cSSMinNOverlapping = 0;

        double sqrtTwo = Math.sqrt(2);

        for (int i = 0; i < (longer.getN() + shorter.getN() - 1); i++) {
            //siIdx is first index in shorter for comparison
            //sfIdx is last index in shorter for comparison
            //liIdx is first index of longer for comparison
            //lfIdx is last index of longer for comparison
            int siIdx, sfIdx, liIdx, lfIdx, offset;
            if (i < shorter.getN()) {
                /*
                 #####
                     +++++++++++ i=0
                  #####
                     +++++++++++
                   #####
                     +++++++++++
                    #####
                     +++++++++++
                     #####
                     +++++++++++ i=4
                */
                //sfIdx is inclusive endpoint
                sfIdx = shorter.getN() - 1;
                siIdx = sfIdx - i;
                liIdx = 0;
                //lfIdx is inclusive endpoint
                lfIdx = sfIdx - siIdx;
                offset = i - sfIdx;

            } else if (i < longer.getN() ) {

                /*
                      #####
                     +++++++++++  i=5

                       #####
                     +++++++++++

                        #####
                     +++++++++++

                         #####
                     +++++++++++

                          #####
                     +++++++++++

                           #####
                     +++++++++++ i = 10
                */
                //sfIdx is inclusive endpoint
                sfIdx = shorter.getN() - 1;
                siIdx = 0;
                liIdx = i - shorter.getN() + 1;
                //lfIdx is inclusive endpoint
                lfIdx = liIdx + shorter.getN() - 1;
                offset = i - sfIdx;

            } else {
                /*
                            #####
                     +++++++++++ i = 12
                     01234567890
                             #####
                     +++++++++++

                              #####
                     +++++++++++

                               #####
                     +++++++++++ i=15
                     01234567890
                */
                liIdx = i - shorter.getN() + 1;
                //sfIdx is inclusive endpoint
                sfIdx = longer.getN() - liIdx - 1;
                siIdx = 0;
                //lfIdx is inclusive endpoint
                lfIdx = liIdx + (sfIdx - siIdx);
                offset = liIdx;

            }

            int nOverLapping = (sfIdx - siIdx) + 1;

            if ((sfIdx - siIdx) != (lfIdx - liIdx)) {
                throw new IllegalStateException(
                    "sample ranges are not correct");
            }

            double sumSq = 0;

            int s = siIdx;
            int l = liIdx;
            while (s <= sfIdx) {
                int xs = shorter.getX(s);
                int xl = longer.getX(l);
                int dx = xs - xl;
                int ys = shorter.getY(s);
                int yl = longer.getY(l);
                int dy = ys - yl;
                sumSq += ((dx*dx) + (dy*dy));
                s++;
                l++;
            }

            double tmp = Math.sqrt(sumSq/nOverLapping);

            // assuming adjacent pixel has distance of sqrt(2) at the most
            if (tmp <= sqrtTwo) {

                if ((tmp < cSSMin) ||
                (tmp == cSSMin && (nOverLapping > cSSMinNOverlapping))
                ) {

                    cSSMin = tmp;

                    cSSMinOffset = offset;

                    cSSMinNOverlapping = nOverLapping;
                }
            }
        }

        if (cSSMin < Double.MAX_VALUE) {

            correlationOffset[0] = cSSMinOffset;

            return true;
        }

        return false;
    }

    public double[] calculateXYCentroids(PairIntArray xy, float[] weights) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {
            double x1 = xy.getX(i);
            xc += (weights[i] * x1);

            double y1 = xy.getY(i);
            yc += (weights[i] * y1);
        }

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(PairIntArray xy) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {

            xc += xy.getX(i);

            yc += xy.getY(i);
        }

        xc /= (double)xy.getN();

        yc /= (double)xy.getN();

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(PairFloatArray xy) {

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < xy.getN(); i++) {

            xc += xy.getX(i);

            yc += xy.getY(i);
        }

        xc /= (double)xy.getN();

        yc /= (double)xy.getN();

        return new double[]{xc, yc};
    }

    public double[] calculateXYCentroids(List<PairIntArray> xyList) {

        double xc = 0;
        double yc = 0;

        for (PairIntArray points : xyList) {

            double[] xycen = calculateXYCentroids(points);

            xc += xycen[0];
            yc += xycen[1];
        }

        xc /= (double)xyList.size();
        yc /= (double)xyList.size();

        return new double[]{xc, yc};
    }

    public PairInt calculateXYCentroids2(Collection<PairInt> points) {
        
        double[] xyCen = calculateXYCentroids(points);
        
        PairInt p = new PairInt((int)Math.round(xyCen[0]), (int)Math.round(xyCen[1]));
        
        return p;
    }
    
    public double[] calculateXYCentroids(Collection<PairInt> points) {

        double xc = 0;
        double yc = 0;

        for (PairInt p : points) {

           int x = p.getX();
           int y = p.getY();

            xc += x;
            yc += y;
        }

        xc /= (double)(points.size());

        yc /= (double)(points.size());

        return new double[]{xc, yc};
    }

    public int[] calculateRoundedXYCentroids(Set<PairInt> points) {

        double[] xyCen = calculateXYCentroids(points);

        int[] out = new int[2];
        out[0] = (int)Math.round(xyCen[0]);
        out[1] = (int)Math.round(xyCen[1]);

        return out;
    }
    
    public int[] calculateRoundedXYCentroids(TIntSet pixelIdxs, int imgWidth) {

        double xc = 0;
        double yc = 0;

        TIntIterator iter = pixelIdxs.iterator();
        while (iter.hasNext()) {

            int pixIdx = iter.next();
            
            int y = pixIdx/imgWidth;
            int x = pixIdx - (y * imgWidth);
            
            xc += x;
            yc += y;
        }

        xc /= (double)pixelIdxs.size();
        yc /= (double)pixelIdxs.size();
        
        int[] out = new int[2];
        out[0] = (int)Math.round(xc);
        out[1] = (int)Math.round(yc);

        return out;
    }

    public double[] calculateXYCentroids(float[] x, float[] y) {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < x.length; i++) {

            xc += x[i];

            yc += y[i];
        }

        xc /= (double)(x.length);

        yc /= (double)(x.length);

        return new double[]{xc, yc};
    }
    
    public double[] calculateXYCentroids(int[] x, int[] y) {

        if (x == null) {
            throw new IllegalArgumentException("x cannot be null");
        }
        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        }
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }

        double xc = 0;
        double yc = 0;

        for (int i = 0; i < x.length; i++) {

            xc += x[i];

            yc += y[i];
        }

        xc /= (double)(x.length);

        yc /= (double)(x.length);

        return new double[]{xc, yc};
    }

    /**
     * does removing the point at idx create a gap between it's neighboring
     * pixels?  this uses the simplest test of only the points at idx-1
     * and idx+1.
     *
     * @param edge
     * @param idx
     * @return
     */
    public boolean doesDisconnect(PairIntArray edge, int idx) {

        // test for endpoints first
        if (idx == 0) {

            if (edge.getN() < 3) {
                return true;
            }

            // does this point currently connect to the last point?
            float diffX = edge.getX(idx) - edge.getX(edge.getN() - 1);
            if (diffX < 0) {
                diffX *= -1;
            }
            float diffY = edge.getY(idx) - edge.getY(edge.getN() - 1);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (((diffX < 2) && (diffY < 2))) {
                // this is connected to the last point in the edge
                // check to see if lastPoint and idx + 1 are adjacent
                diffX = edge.getX(idx + 1) - edge.getX(edge.getN() - 1);
                if (diffX < 0) {
                    diffX *= -1;
                }
                if (diffX > 1) {
                    return true;
                }

                diffY = edge.getY(idx + 1) - edge.getY(edge.getN() - 1);
                if (diffY < 0) {
                    diffY *= -1;
                }
                if (diffY > 1) {
                    return true;
                }
            }
            return false;
        }

        if (idx == (edge.getN() - 1)) {

            if (edge.getN() < 3) {
                return true;
            }

            // does this point currently connect to the first point?
            float diffX = edge.getX(idx) - edge.getX(0);
            if (diffX < 0) {
                diffX *= -1;
            }
            float diffY = edge.getY(idx) - edge.getY(0);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (((diffX < 2) && (diffY < 2))) {
                // this is connected to the first point in the edge
                // check to see if lastPoint - 1 and first point are adjacent
                diffX = edge.getX(idx - 1) - edge.getX(0);
                if (diffX < 0) {
                    diffX *= -1;
                }
                if (diffX > 1) {
                    return true;
                }

                diffY = edge.getY(idx - 1) - edge.getY(0);
                if (diffY < 0) {
                    diffY *= -1;
                }
                if (diffY > 1) {
                    return true;
                }
            }
            return false;
        }

        if ((idx + 1) < edge.getN()) {
            float diffX = edge.getX(idx - 1) - edge.getX(idx + 1);
            if (diffX < 0) {
                diffX *= -1;
            }
            if (diffX > 1) {
                return true;
            }

            float diffY = edge.getY(idx - 1) - edge.getY(idx + 1);
            if (diffY < 0) {
                diffY *= -1;
            }
            if (diffY > 1) {
                return true;
            }

            return false;
        }

        return false;
    }

    public double distanceFromPointToALine(float lineX0, float lineY0,
        float lineX1, float lineY1, float xP, float yP) {

        /*
        en.wikipedia.org/wiki/Distance_from_a_point_to_a_line

        for the edge, we have the 2 points (lineX0, lineY0) and (lineX1, lineY1)

        distance between that edge and a point (xP, yP) is

        defining diffX = lineX1 - lineX0
                 diffY = lineY1 - lineY0;

        d =
           ( diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0 )
           ( --------------------------------------------------- )
           (         (diffX*diffX + diffY*diffY)^0.5             )
        )
        */

        float diffX = lineX1 - lineX0;
        float diffY = lineY1 - lineY0;

        if (diffY == 0) {
            // horizontal line
            return Math.abs(yP - lineY0);
        } else if (diffX == 0) {
            // vertical line
            return Math.abs(xP - lineX0);
        }

        double pt1 = Math.abs(diffY*xP - diffX*yP - lineX0*lineY1 + lineX1*lineY0);

        double pt2 = Math.sqrt(diffX*diffX + diffY*diffY);

        double dist = pt1/pt2;

        return dist;
    }

    public void sortByX(PairIntArray curve) {
        if (curve.getN() < 2) {
            return;
        }
        sortByX(curve, 0, curve.getN() - 1);
    }

    private void sortByX(PairIntArray curve, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = partitionByX(curve, idxLo, idxHi);
            sortByX(curve, idxLo, idxMid - 1);
            sortByX(curve, idxMid + 1, idxHi);
        }
    }

    private int partitionByX(PairIntArray curve, int idxLo, int idxHi) {

        int x = curve.getX(idxHi);  //for comparison
        int store = idxLo - 1;      //store to swap after pivot

        for (int i = idxLo; i < idxHi; i++) {
            if (curve.getX(i) <= x) {
                store++;
                int swapX = curve.getX(store);
                int swapY = curve.getY(store);
                curve.set(store, curve.getX(i), curve.getY(i));
                curve.set(i, swapX, swapY);
            }
        }
        store++;

        int swapX = curve.getX(store);
        int swapY = curve.getY(store);
        curve.set(store, curve.getX(idxHi), curve.getY(idxHi));
        curve.set(idxHi, swapX, swapY);

        return store;
    }

    protected void findNeighbors(int x, int y, Set<PairInt> outputNeighbors,
        Set<PairInt> points, Set<PairInt> tmpAddedPoints,
        Set<PairInt> tmpRemovedPoints, int imageWidth, int imageHeight) {

        outputNeighbors.clear();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (tmpRemovedPoints.contains(p2)) {
                continue;
            }
            if (tmpAddedPoints.contains(p2) || points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    public Set<PairInt> findNeighbors(int x, int y, Set<PairInt> points) {

        Set<PairInt> neighbors = new HashSet<PairInt>();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                neighbors.add(p2);
            }
        }

        return neighbors;
    }
        
    public void findNeighbors(int x, int y, Set<PairInt> points, 
        Set<PairInt> excludePoints, int[] dxs, int[] dys, 
        Set<PairInt> outputNeighbors) {
        
        outputNeighbors.clear();
        
        for (int i = 0; i < dxs.length; i++) {
            
            int x2 = x + dxs[i];
            int y2 = y + dys[i];
            
            PairInt p2 = new PairInt(x2, y2);
            
            if (excludePoints.contains(p2)) {
                continue;
            }
            if (points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    public void findNeighbors(int x, int y, Set<PairInt> outputNeighbors,
        Set<PairInt> points, Set<PairInt> excludePoints, int imageWidth, int imageHeight) {

        outputNeighbors.clear();

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (excludePoints.contains(p2)) {
                continue;
            }
            if (points.contains(p2)) {
                outputNeighbors.add(p2);
            }
        }
    }

    /**
     * iterate through points, counting the number of pixels on the image
     * boundaries, and return true if the number reaches numberOfPixels.
     * @param numberOfPixels the number of pixels for which to return true
     * if they are on the image boundaries.
     * @param points
     * @param imageWidth
     * @param imageHeight
     * @return
     */
    public boolean hasNumberOfPixelsOnImageBoundaries(int numberOfPixels,
        Set<PairInt> points, int imageWidth, int imageHeight) {

        int n = 0;

        for (PairInt p : points) {

            int x = p.getX();
            int y = p.getY();

            if ((x == 0) || (y == 0) || (x == (imageWidth - 1)) ||
                (y == (imageHeight - 1))) {

                n++;

                if (n == numberOfPixels) {
                    return true;
                }
            }
        }

        return (n >= numberOfPixels);
    }

    public int countNeighbors(int x, int y, Set<PairInt> points, int imageWidth,
        int imageHeight) {

        int nn = 0;

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            if ((x2 < 0) || (x2 > (imageWidth - 1)) || (y2 < 0) ||
                (y2 > (imageHeight - 1))) {
                continue;
            }

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                nn++;
            }
        }

        return nn;
    }
    
    public int countNeighbors(int x, int y, Set<PairInt> points) {

        int nn = 0;

        for (int i = 0; i < eightNeighborsX.length; i++) {

            int x2 = x + eightNeighborsX[i];
            int y2 = y + eightNeighborsY[i];

            PairInt p2 = new PairInt(x2, y2);

            if (points.contains(p2)) {
                nn++;
            }
        }

        return nn;
    }

    public boolean isAdjacent(PairIntArray edge, int idx1, int idx2) {
        
        if (idx2 < 0) {
            return false;
        }

        int x1 = edge.getX(idx1);
        int y1 = edge.getY(idx1);

        int x2 = edge.getX(idx2);
        int y2 = edge.getY(idx2);

        int diffX = Math.abs(x1 - x2);
        int diffY = Math.abs(y1 - y2);

        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }

        return false;
    }
    
    public boolean isAdjacent(PairIntArray edge, int idx, int x, int y) {
        
        if (idx < 0 || idx > (edge.getN() - 1)) {
            return false;
        }

        int x1 = edge.getX(idx);
        int y1 = edge.getY(idx);

        int diffX = Math.abs(x1 - x);
        int diffY = Math.abs(y1 - y);

        if ((diffX < 2) && (diffY < 2)) {
            return true;
        }

        return false;
    }
    
    public boolean isAdjacent(PairIntArray edge, int idx1, int idx2,
        float spacingBetweenPoints) {

        int x1 = edge.getX(idx1);
        int y1 = edge.getY(idx1);

        int x2 = edge.getX(idx2);
        int y2 = edge.getY(idx2);

        int diffX = Math.abs(x1 - x2);
        int diffY = Math.abs(y1 - y2);
        
        float dist = (float)Math.sqrt(diffX*diffX + diffY*diffY);

        if (dist <= spacingBetweenPoints) {
            return true;
        }

        return false;
    }

    /**
     * given 3 counter-clockwise ordered points on a curve, calculate the angle 
     * along the curve at the middle point, its directionCCW is from p0 to p1.
     * <pre>
     * For example:
     * 
     * 135 degrees
     *       .---
     *       | .
     *           p2   
     *             p1
     *                p0
     * </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public double calculateAngleAtMidpoint(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        /*
        given the points (x1, y1) (x2, y2) and (x3, y3), 
        calculates the angle at the midpoint (x2, y2) for the path along
        the points.
        */
        
        double theta1 = AngleUtil.polarAngleCCW(x2 - x1, y2 - y1);
        
        double theta2 = AngleUtil.polarAngleCCW(x3 - x2, y3 - y2);
        
        double theta = AngleUtil.getAngleAverageInRadians(theta1, theta2);
                
        return theta;
    }
    
    /**
     * given 3 counter-clockwise ordered points on a curve, calculate the angle 
     * tangent to the curve at the middle point - its directionCCW follows
     * the right hand rule.
     * <pre>
     * For example:
     *                  45 degrees
     *               __
     *               . |
     *       p2    .
     *          p1
     *             p0
     * </pre>
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param x3
     * @param y3
     * @return
     */
    public double calculateAngleTangentToMidpoint(int x1, int y1, 
        int x2, int y2, int x3, int y3) {

        double theta = calculateAngleAtMidpoint(x1, y1, x2, y2, x3, y3);
               
        double thetaMinus90 = theta - Math.PI/2.;
        if (thetaMinus90 < 0) {
            thetaMinus90 += (2.*Math.PI);
        }
        
        return thetaMinus90;
    }

    /**
     * given theta and the point (xp, yp), determine which directionCCW and hence
     * polar angle (clockwise) is perpendicular away from the centroid.
     * The reference point (xm, ym) is the point from which theta was also
     * calculated, which is probably the point for kMaxIdx.  The points are also
     * checked to make sure they aren't in the points set.
     *
     * @param theta
     * @param xp
     * @param yp
     * @param xm
     * @param ym
     * @param centroidXY
     * @param points
     * @return
     */
    public double calculatePerpendicularAngleAwayFromCentroid(
        double theta, int xp, int yp, int xm, int ym, double[] centroidXY,
        Set<PairInt> points) {

        /*
        rotate the point (xm, ym) around (xp, yp) 90 degrees and -90 degrees.
        The rotated point which is furthest from the centroid is the
        directionCCW of the vector pointing away from the centroid.
        */

        /*
        math.cos(math.pi/2) = 0
        math.sin(math.pi/2) = 1
        math.sin(-math.pi/2) = -1

        double xr = centroidX + ((y - centroidY) * sine(angle)));
        double yr = centroidY + ((-(x - centroidX) * sine(angle)))
        */

        int xmRot90 = xp + (ym - yp);
        int ymRot90 = yp + (-(xm - xp));

        int xmRotNegative90 = xp  - (ym - yp);
        int ymRotNegative90 = yp + (xm - xp);

        boolean rot90IsInPoints = points.contains(
            new PairInt(Math.round(xmRot90), Math.round(ymRot90)));

        boolean rotNegative90IsInPoints = points.contains(
            new PairInt(Math.round(xmRotNegative90),
            Math.round(ymRotNegative90)));

        double distSqRot90 = (xmRot90 - centroidXY[0]) * (xmRot90 - centroidXY[0])
            + (ymRot90 - centroidXY[1]) * (ymRot90 - centroidXY[1]);

        double distSqRotNegative90 =
            (xmRotNegative90 - centroidXY[0]) * (xmRotNegative90 - centroidXY[0])
            + (ymRotNegative90 - centroidXY[1]) * (ymRotNegative90 - centroidXY[1]);

        double perp = theta;

        if (distSqRot90 > distSqRotNegative90) {
            perp += Math.PI/2.;
        } else if (distSqRot90 > distSqRotNegative90) {
            if (rot90IsInPoints && !rotNegative90IsInPoints) {
                perp -= Math.PI/2.;
            } else if (!rot90IsInPoints && rotNegative90IsInPoints) {
                perp += Math.PI/2.;
            } else {
                throw new IllegalStateException("Error in algorithm:" +
                " consider changing the test 90 and -90 points so that" +
                " one will always be in points set.");
            }
        } else {
            perp -= Math.PI/2.;
        }

        if (perp >= 2*Math.PI) {
            perp = perp - 2*Math.PI;
        } else if (perp < 0) {
            perp += 2*Math.PI;
        }

        return perp;
    }

    public double calculateArea(PairIntArray closedCurve) {
        
        int n = closedCurve.getN();
        
        double sum = 0;
        
        for (int i = 0; i < (n - 1); ++i) {
            
            double t = 0.5 * (closedCurve.getY(i + 1) + closedCurve.getY(i)) *
                (closedCurve.getX(i + 1) - closedCurve.getX(i));
            
            sum += t;
        }
        
        sum += ((closedCurve.getY(0) + closedCurve.getY(n - 1)) *
                (closedCurve.getX(0) - closedCurve.getX(n -1)));
        
        return sum;
    }

    public PairIntArray createContiguousCircle(float radius) {

        int shift = (int)Math.ceil(radius);
    
        return createContiguousCircle(radius, shift, shift);
    }
    
    public PairIntArray createContiguousCircle(float radius, int xShift, 
        int yShift) {
        
        // for a change in y to be at least 1 pixel, theta would be:
        //   theta = asin(1/r)
        double theta = Math.asin(1./radius);
        
        Set<PairInt> added = new HashSet<PairInt>();
        
        PairIntArray circle = new PairIntArray();
        double t = 0;
        int x, y;
        while (t <= Math.PI/2.) {
            x = xShift + (int)Math.round(radius * Math.cos(t));
            y = yShift + (int)Math.round(radius * Math.sin(t));
            PairInt p = new PairInt(x, y);
            t += theta;
            if (added.contains(p)) {
                continue;
            }
            circle.add(x, y);
            added.add(p);
        }
        int n90 = circle.getN();
        for (int i = (n90 - 1); i > -1; --i) {
            x = xShift -1 * (circle.getX(i) - xShift);
            y = circle.getY(i);
            PairInt p = new PairInt(x, y);
            if (added.contains(p)) {
                continue;
            }
            circle.add(x, y);
            added.add(p);
        }
        int n180 = circle.getN();
        for (int i = (n180 - 1); i > -1; --i) {
            x = circle.getX(i);
            y = yShift -1 * (circle.getY(i) - yShift);
            PairInt p = new PairInt(x, y);
            if (added.contains(p)) {
                continue;
            }
            circle.add(x, y);
            added.add(p);
        }
        
        return circle;
    }
}
