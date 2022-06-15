package algorithms.statistics;

public class ModelPrediction {
    protected double[] yFit;
    protected double[] yErr;

    public double[] getYErr() {
        return yErr;
    }

    public void setYErr(double[] yErr) {
        this.yErr = yErr;
    }

    public double[] getYFit() {
        return yFit;
    }

    public void setYFit(double[] yFit) {
        this.yFit = yFit;
    }
}
