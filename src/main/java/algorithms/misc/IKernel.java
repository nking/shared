package algorithms.misc;

/**
 *
 * @author nichole
 */
public interface IKernel {

    /**
     *
     @param x
     @param xTilde
     @param h
     @return
     */
    public double kernel(double x, double[] xTilde, double h);

    /**
     *
     @param x
     @return
     */
    public double dKdx(double x);

    /**
     * 2nd derivative of the kernel w.r.t. x
     @param x
     @return
     */
    public double d2Kdx2(double x);

}
