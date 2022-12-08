package algorithms.misc;

/**
 *
 * @author nichole
 */
public class EpanechnikovKernel implements IKernel {

    /**
     *
     @param x
     @param xTilde
     @param h
     @return
     */
    public double kernel(double x, double[] xTilde, double h) {
        int nS = xTilde.length;
        double s5 = Math.sqrt(5);
        double z;
        double sum = 0;
        for (int i = 0; i < nS; i++) {
            z = (x - xTilde[i])/h;
            if (Math.abs(z) < s5) {
                sum += (1. - (z*z/5.));
            }
        }
        sum *= (0.75/s5);
        //(0.75 * (1. - (x*x/5.)) )/ s5
        return sum/nS;
    }

    /**
     * derivative of the kernel w.r.t. x
     @param x
     @return
     */
    public double dKdx(double x) {
        return 3.*x/10.;
    }
    /**
     * 2nd derivative of the kernel w.r.t. x
     @param x
     @return
     */
    public double d2Kdx2(double x) {
        return 3./10.;
    }
}
