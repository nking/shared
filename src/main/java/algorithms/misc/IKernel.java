package algorithms.misc;

public interface IKernel {

    public double kernel(double x, double[] xTilde, double h);

    public double dKdx(double x);

    /**
     * 2nd derivative of the kernel w.r.t. x
     * @param x
     * @return
     */
    public double d2Kdx2(double x);

}
