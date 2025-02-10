package algorithms.statistics;

import junit.framework.TestCase;
import no.uib.cipr.matrix.NotConvergedException;

public class CovarianceTest extends TestCase {

    public void test() {
        double[] x = new double[]{0,2,4,5,7,8,10};
        double[] y = new double[]{4,3,1,1,2,1,0};
        double expCov = -4.12;
        double expCor = - 0.86;
        double tol = 1E-2;

        double[] xc = Standardization.zeroCenterMean(x);
        double[] yc = Standardization.zeroCenterMean(y);

        double[] outXMean = new double[1];
        double[] outYMean = new double[1];
        double[] outXStdv = new double[1];
        double[] outYStdv = new double[1];
        double[] xs = Standardization.standardUnitNormalization(x, 1, outXMean, outXStdv);
        double[] ys = Standardization.standardUnitNormalization(y, 1, outYMean, outYStdv);

        double corPearson = Covariance.correlationPearson(x, y);
        double corrUnCentered = Covariance.correlationSample(x, y,
                Covariance.STAND_TYPE.UNCENTERED);

        assertTrue(Math.abs(corPearson - expCor) < tol);
        assertTrue(Math.abs(corrUnCentered - expCor) < tol);

        double cov = Covariance.covarianceSample(x, y, Covariance.STAND_TYPE.UNCENTERED);
        assertTrue(Math.abs(cov - expCov) < tol);
        cov = Covariance.covarianceSample(xc, yc, Covariance.STAND_TYPE.MEAN_CENTERED);
        assertTrue(Math.abs(cov - expCov) < tol);
        cov = Covariance.covarianceSample(xs, ys, Covariance.STAND_TYPE.UNIT_STAND_MEAN0_STD1);
        cov *= (outXStdv[0] * outYStdv[0]);
        assertTrue(Math.abs(cov - expCov) < tol);
    }

    public void test11() throws NotConvergedException {

        double[] x = new double[]{13,14,15,16,16,18};
        double[] y = new double[]{8,10,15,20,27,30};
        double[] xc = Standardization.zeroCenterMean(x);
        double[] yc = Standardization.zeroCenterMean(y);

        double[] outXMean = new double[1];
        double[] outYMean = new double[1];
        double[] outXStdv = new double[1];
        double[] outYStdv = new double[1];

        double expLinearCor = 0.948;

        double[] xs = Standardization.standardUnitNormalization(x, 1, outXMean, outXStdv);
        double[] ys = Standardization.standardUnitNormalization(y, 1, outYMean, outYStdv);

        double corPearson = Covariance.correlationPearson(x, y);
        double corrUnCentered = Covariance.correlationSample(x, y,
                Covariance.STAND_TYPE.UNCENTERED);
        double corrMeanCentered = Covariance.correlationSample(xc, yc,
                Covariance.STAND_TYPE.MEAN_CENTERED);
        double corrUnit = Covariance.correlationSample(xs, ys,
                Covariance.STAND_TYPE.UNIT_STAND_MEAN0_STD1);

        assertTrue(Math.abs(corPearson - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrUnCentered - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrMeanCentered - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrUnit - 0) < 1E-3);
    }

    public void testCorrelationCovarianceDependent() throws NotConvergedException {

        double[] x = new double[]{1, 2, 3, 5, 8};
        double[] y = new double[]{0.11, 0.12, 0.13, 0.15, 0.18};
        double[] xc = Standardization.zeroCenterMean(x);
        double[] yc = Standardization.zeroCenterMean(y);

        double[] outXMean = new double[1];
        double[] outYMean = new double[1];
        double[] outXStdv = new double[1];
        double[] outYStdv = new double[1];
        /*
        double[] standardUnitNormalization(double[] data,
        int nDimensions, double[] outputMean, double[] outputStandardDeviation)
         */
        double[] xs = Standardization.standardUnitNormalization(x, 1, outXMean, outXStdv);
        double[] ys = Standardization.standardUnitNormalization(y, 1, outYMean, outYStdv);

        // y = 0.10 + 0.01 x
        double expLinearCor = 1.0;
        double expLinearCov = expLinearCor*(outXStdv[0] * outYStdv[0]);
        double cosineSim = 0.92081;

        double corPearson = Covariance.correlationPearson(x, y);
        double corrUnCentered = Covariance.correlationSample(x, y,
                Covariance.STAND_TYPE.UNCENTERED);
        double corrMeanCentered = Covariance.correlationSample(xc, yc,
                Covariance.STAND_TYPE.MEAN_CENTERED);
        double corrUnit = Covariance.correlationSample(xs, ys,
                Covariance.STAND_TYPE.UNIT_STAND_MEAN0_STD1);

        double covUncentered = Covariance.covarianceSample(x, y, Covariance.STAND_TYPE.UNCENTERED);
        double covMeanCentered = Covariance.covarianceSample(xc, yc, Covariance.STAND_TYPE.MEAN_CENTERED);
        double covUnit = Covariance.covarianceSample(xs, ys, Covariance.STAND_TYPE.UNIT_STAND_MEAN0_STD1);
        covUnit *= (outXStdv[0] * outYStdv[0]);

        assertTrue(Math.abs(corPearson - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrUnCentered - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrMeanCentered - expLinearCor) < 1E-3);
        assertTrue(Math.abs(corrUnit - 0) < 1E-3);

        double tol = 1E-2;
        assertTrue(Math.abs(covUncentered - expLinearCov) < tol);
        assertTrue(Math.abs(covMeanCentered - expLinearCov) < tol);
        assertTrue(Math.abs(covUnit - expLinearCov) < tol);

    }
}
