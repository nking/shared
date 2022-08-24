package algorithms.misc;

public class GaussianKernel implements IKernel {

    /**
     * calculate (1./(h * n * Math.sqrt(2. * Math.PI)) * summation over i of exp(-0.5*( (x - xTilde[i])/h)^2 )
     * runtime complexity is O(xTilde.length).
     * @param x the grid point
     * @param xTilde the observed data
     * @param h
     * @return
     */
    public double kernel(double x, double[] xTilde, double h) {
        int nS = xTilde.length;
        double ch = 1./(Math.sqrt(2. * Math.PI));
        double z;
        double sum = 0;
        for (int i = 0; i < nS; i++) {
            z = (x - xTilde[i])/h;
            sum += Math.exp(-0.5 * z * z);
        }
        sum *= ch/(h*nS);
        return sum;
    }

    /**
     * calculate (1./(Math.sqrt(2. * Math.PI)) * exp(-0.5*( z^2 )
     * @param z
     * @return
     */
    public double kernel(double z) {
        double ch = 1./(Math.sqrt(2. * Math.PI));
        return ch * Math.exp(-0.5 * z * z);
    }

    /**
     * derivative of the kernel w.r.t. x
     * @param x
     * @return
     */
    public double dKdx(double x) {
        double c = (1./Math.sqrt(2. * Math.PI));
        double e = Math.exp(-x*x/2.);
        return -x * c * e;
    }
    /**
     * 2nd derivative of the kernel w.r.t. x
     * @param x
     * @return
     */
    public double d2Kdx2(double x) {
        double c = (1./Math.sqrt(2. * Math.PI));
        double e = Math.exp(-x*x/2.);
        return  (-x*x - 1.) * c * e;
    }
}
