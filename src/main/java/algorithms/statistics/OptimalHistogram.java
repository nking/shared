package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class OptimalHistogram {
    
    // implementation of 
    // "Optimal Histograms With Quality Guarantees",
    // 1998, Jagadish et al., Proceedings of the 24th VLDB Conference.
    //http://infolab.stanford.edu/~datar/courses/cs361a/papers/vopt.pdf
    
    // also see curvature project for fast integer multi-dimensional histograms:
    // src/algorithms/imageProcessing/OtsuThresholding.java
    
    // KDE's
    // for gaussian smoothing: h = 1.06 * sigma * n^(-1/5)
        //     or, let A=min(standard deviation, interquantile range/1.34).
        //      then h = 0.9 * A * n^(-1/5)
        //
}
