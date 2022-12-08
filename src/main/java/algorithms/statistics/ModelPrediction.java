package algorithms.statistics;

/**
 *
 * @author nichole
 */
public class ModelPrediction {

    /**
     *
     */
    protected double[] yFit;

    /**
     *
     */
    protected double[] yErr;

    /**
     *
     @return
     */
    public double[] getYErr() {
        return yErr;
    }

    /**
     *
     @param yErr
     */
    public void setYErr(double[] yErr) {
        this.yErr = yErr;
    }

    /**
     *
     @return
     */
    public double[] getYFit() {
        return yFit;
    }

    /**
     *
     @param yFit
     */
    public void setYFit(double[] yFit) {
        this.yFit = yFit;
    }
}
