package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class Beta {
    
    /**
     * implemented as gamma(alpha)*gamma(beta)/gamma(alpha + beta)
     * 
     @param alpha alpha parameter
     @param beta beta parameter
     @return the value for beta(alpha, beta)
     */
    public static double beta(double alpha, double beta) {
        double ga = Gamma.lanczosGamma9(alpha);
        double gb = Gamma.lanczosGamma9(beta);
        double gab = Gamma.lanczosGamma9(alpha + beta);
        return (ga * gb)/gab;
    }
}
