package algorithms.compGeometry.convexHull;

/**
 *adapted from 
 * https://code.google.com/p/two-point-correlation/source/browse/src/test/java/algorithms/compGeometry/convexHull/
 * under MIT License (MIT), Nichole King 2013

 then moved to project:
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 */
public class GrahamScanTooFewPointsException extends Exception {

    protected static final long serialVersionUID = 12345678;

    GrahamScanTooFewPointsException(String errorMessage) {
        super(errorMessage);
    }

}
