package algorithms.statistics;

/*
methods for the critical values of the chi-squared distribution to look up
the chi-squared statistic test or limiting probability (=p-value)
for the degrees of freedom.

Methods and 2 chi-squared inverse survivial function tables
are given from NIST
https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm

and 1988, Lin, J.T., The Statistician 37, 3-5

NOTE: The cumulative distribution function and its inverse (quantile function) 
of the chi-square distribution do not have closed form so the pre-calculated
tables, or approximations are needed.

 * @author nichole
 */
public class ChiSquaredCriticalValues {
    
    /*
    For a two-sided test:
        find the column corresponding to 1-α/2 in the table for
        upper-tail critical values 
        and reject the null hypothesis if the test statistic is greater than the 
        table value. 

        Similarly, find the column corresponding to α/2 in the table for 
        lower-tail critical values and reject the null hypothesis if the test 
        statistic is less than the table value.
    
    For a one-sided upper-tail test:
        find the column corresponding to 1-α in the table containing 
        upper-tail critical and reject the null hypothesis if the test 
        statistic is greater than the table value.
    
    For a one-sided lower-tail test:
        find the column corresponding to α in the lower-tail critical values 
        table and reject the null hypothesis if the computed test statistic 
        is less than the tabled value.
    */
    
    public static enum PROB_UT {
        Z_90, Z_95, Z_975, Z_99, Z_999;
    }
    public static enum PROB_LT {
        Z_10, Z_05, Z_025, Z_01, Z_001;
    }
    
    /**
     * a table of upper-tail critical values of the chi-squared distribution
     * where columns are the limit probabilities and rows are the 
     * degrees of freedom.  
     * The values are a member of the variate X.
     * A quantile is a range within the set of X.
     * One can reject the null hypothesis if the test statistic is greater 
     * than the tabled value.
     */
    static double[][] upperTail = new double[100][];
    static {                      //Z_90, Z_95,  Z_975,  Z_99,  Z_999;
        upperTail[0] = new double[]{2.706, 3.841, 5.024, 6.635, 10.828};
        upperTail[1] = new double[]{4.605, 5.991, 7.378, 9.210, 13.816};
        upperTail[2] = new double[]{6.251, 7.815, 9.348, 11.345, 16.266};
        upperTail[3] = new double[]{7.779, 9.488, 11.143, 13.277, 18.467};
        upperTail[4] = new double[]{9.236, 11.070, 12.833, 15.086, 20.515};
        upperTail[5] = new double[]{10.645, 12.592, 14.449, 16.812, 22.458};
        upperTail[6] = new double[]{12.017, 14.067, 16.013, 18.475, 24.322};
        upperTail[7] = new double[]{13.362, 15.507, 17.535, 20.090, 26.125};
        upperTail[8] = new double[]{14.684, 16.919, 19.023, 21.666, 27.877};
        upperTail[9] = new double[]{15.987, 18.307, 20.483, 23.209, 29.588};
        upperTail[10] = new double[]{17.275, 19.675, 21.920, 24.725, 31.264};
        upperTail[11] = new double[]{18.549, 21.026, 23.337, 26.217, 32.910};
        upperTail[12] = new double[]{19.812, 22.362, 24.736, 27.688, 34.528};
        upperTail[13] = new double[]{21.064, 23.685, 26.119, 29.141, 36.123};
        upperTail[14] = new double[]{22.307, 24.996, 27.488, 30.578, 37.697};
        upperTail[15] = new double[]{23.542, 26.296, 28.845, 32.000, 39.252};
        upperTail[16] = new double[]{24.769, 27.587, 30.191, 33.409, 40.790};
        upperTail[17] = new double[]{25.989, 28.869, 31.526, 34.805, 42.312};
        upperTail[18] = new double[]{27.204, 30.144, 32.852, 36.191, 43.820};
        upperTail[19] = new double[]{28.412, 31.410, 34.170, 37.566, 45.315};
        upperTail[20] = new double[]{29.615, 32.671, 35.479, 38.932, 46.797};
        upperTail[21] = new double[]{30.813, 33.924, 36.781, 40.289, 48.268};
        upperTail[22] = new double[]{32.007, 35.172, 38.076, 41.638, 49.728};
        upperTail[23] = new double[]{33.196, 36.415, 39.364, 42.980, 51.179};
        upperTail[24] = new double[]{34.382, 37.652, 40.646, 44.314, 52.620};
        upperTail[25] = new double[]{35.563, 38.885, 41.923, 45.642, 54.052};
        upperTail[26] = new double[]{36.741, 40.113, 43.195, 46.963, 55.476};
        upperTail[27] = new double[]{37.916, 41.337, 44.461, 48.278, 56.892};
        upperTail[28] = new double[]{39.087, 42.557, 45.722, 49.588, 58.301};
        upperTail[29] = new double[]{40.256, 43.773, 46.979, 50.892, 59.703};
        upperTail[30] = new double[]{41.422, 44.985, 48.232, 52.191, 61.098};
        upperTail[31] = new double[]{42.585, 46.194, 49.480, 53.486, 62.487};
        upperTail[32] = new double[]{43.745, 47.400, 50.725, 54.776, 63.870};
        upperTail[33] = new double[]{44.903, 48.602, 51.966, 56.061, 65.247};
        upperTail[34] = new double[]{46.059, 49.802, 53.203, 57.342, 66.619};
        upperTail[35] = new double[]{47.212, 50.998, 54.437, 58.619, 67.985};
        upperTail[36] = new double[]{48.363, 52.192, 55.668, 59.893, 69.347};
        upperTail[37] = new double[]{49.513, 53.384, 56.896, 61.162, 70.703};
        upperTail[38] = new double[]{50.660, 54.572, 58.120, 62.428, 72.055};
        upperTail[39] = new double[]{51.805, 55.758, 59.342, 63.691, 73.402};
        upperTail[40] = new double[]{52.949, 56.942, 60.561, 64.950, 74.745};
        upperTail[41] = new double[]{54.090, 58.124, 61.777, 66.206, 76.084};
        upperTail[42] = new double[]{55.230, 59.304, 62.990, 67.459, 77.419};
        upperTail[43] = new double[]{56.369, 60.481, 64.201, 68.710, 78.750};
        upperTail[44] = new double[]{57.505, 61.656, 65.410, 69.957, 80.077};
        upperTail[45] = new double[]{58.641, 62.830, 66.617, 71.201, 81.400};
        upperTail[46] = new double[]{59.774, 64.001, 67.821, 72.443, 82.720};
        upperTail[47] = new double[]{60.907, 65.171, 69.023, 73.683, 84.037};
        upperTail[48] = new double[]{62.038, 66.339, 70.222, 74.919, 85.351};
        upperTail[49] = new double[]{63.167, 67.505, 71.420, 76.154, 86.661};
        upperTail[50] = new double[]{64.295, 68.669, 72.616, 77.386, 87.968};
        upperTail[51] = new double[]{65.422, 69.832, 73.810, 78.616, 89.272};
        upperTail[52] = new double[]{66.548, 70.993, 75.002, 79.843, 90.573};
        upperTail[53] = new double[]{67.673, 72.153, 76.192, 81.069, 91.872};
        upperTail[54] = new double[]{68.796, 73.311, 77.380, 82.292, 93.168};
        upperTail[55] = new double[]{69.919, 74.468, 78.567, 83.513, 94.461};
        upperTail[56] = new double[]{71.040, 75.624, 79.752, 84.733, 95.751};
        upperTail[57] = new double[]{72.160, 76.778, 80.936, 85.950, 97.039};
        upperTail[58] = new double[]{73.279, 77.931, 82.117, 87.166, 98.324};
        upperTail[59] = new double[]{74.397, 79.082, 83.298, 88.379, 99.607};
        upperTail[60] = new double[]{75.514, 80.232, 84.476, 89.591, 100.888};
        upperTail[61] = new double[]{76.630, 81.381, 85.654, 90.802, 102.166};
        upperTail[62] = new double[]{77.745, 82.529, 86.830, 92.010, 103.442};
        upperTail[63] = new double[]{78.860, 83.675, 88.004, 93.217, 104.716};
        upperTail[64] = new double[]{79.973, 84.821, 89.177, 94.422, 105.988};
        upperTail[65] = new double[]{81.085, 85.965, 90.349, 95.626, 107.258};
        upperTail[66] = new double[]{82.197, 87.108, 91.519, 96.828, 108.526};
        upperTail[67] = new double[]{83.308, 88.250, 92.689, 98.028, 109.791};
        upperTail[68] = new double[]{84.418, 89.391, 93.856, 99.228, 111.055};
        upperTail[69] = new double[]{85.527, 90.531, 95.023, 100.425, 112.317};
        upperTail[70] = new double[]{86.635, 91.670, 96.189, 101.621, 113.577};
        upperTail[71] = new double[]{87.743, 92.808, 97.353, 102.816, 114.835};
        upperTail[72] = new double[]{88.850, 93.945, 98.516, 104.010, 116.092};
        upperTail[73] = new double[]{89.956, 95.081, 99.678, 105.202, 117.346};
        upperTail[74] = new double[]{91.061, 96.217, 100.839, 106.393, 118.599};
        upperTail[75] = new double[]{92.166, 97.351, 101.999, 107.583, 119.850};
        upperTail[76] = new double[]{93.270, 98.484, 103.158, 108.771, 121.100};
        upperTail[77] = new double[]{94.374, 99.617, 104.316, 109.958, 122.348};
        upperTail[78] = new double[]{95.476, 100.749, 105.473, 111.144, 123.594};
        upperTail[79] = new double[]{96.578, 101.879, 106.629, 112.329, 124.839};
        upperTail[80] = new double[]{97.680, 103.010, 107.783, 113.512, 126.083};
        upperTail[81] = new double[]{98.780, 104.139, 108.937, 114.695, 127.324};
        upperTail[82] = new double[]{99.880, 105.267, 110.090, 115.876, 128.565};
        upperTail[83] = new double[]{100.980, 106.395, 111.242, 117.057, 129.804};
        upperTail[84] = new double[]{102.079, 107.522, 112.393, 118.236, 131.041};
        upperTail[85] = new double[]{103.177, 108.648, 113.544, 119.414, 132.277};
        upperTail[86] = new double[]{104.275, 109.773, 114.693, 120.591, 133.512};
        upperTail[87] = new double[]{105.372, 110.898, 115.841, 121.767, 134.746};
        upperTail[88] = new double[]{106.469, 112.022, 116.989, 122.942, 135.978};
        upperTail[89] = new double[]{107.565, 113.145, 118.136, 124.116, 137.208};
        upperTail[90] = new double[]{108.661, 114.268, 119.282, 125.289, 138.438};
        upperTail[91] = new double[]{109.756, 115.390, 120.427, 126.462, 139.666};
        upperTail[92] = new double[]{110.850, 116.511, 121.571, 127.633, 140.893};
        upperTail[93] = new double[]{111.944, 117.632, 122.715, 128.803, 142.119};
        upperTail[94] = new double[]{113.038, 118.752, 123.858, 129.973, 143.344};
        upperTail[95] = new double[]{114.131, 119.871, 125.000, 131.141, 144.567};
        upperTail[96] = new double[]{115.223, 120.990, 126.141, 132.309, 145.789};
        upperTail[97] = new double[]{116.315, 122.108, 127.282, 133.476, 147.010};
        upperTail[98] = new double[]{117.407, 123.225, 128.422, 134.642, 148.230};
        upperTail[99] = new double[]{118.498, 124.342, 129.561, 135.807, 149.449};
    }
    
    /**
     * find the upper-tail critical value to compare to a test statistic.
     * One can reject the null hypothesis if the computed test statistic 
        is greater than the table value.
     * @param probability the probability to look-up in the table.  
     * e.g. for a one-sided test, this probability is the same as 1 minus the significance 
     * level alpha.  For a two-sided test, this probability is the same as the significance 
     * level 1 minus alpha/2.
     * @param degreesOfFreedom
     * @return 
     */
    public static double upperTailStatistic(PROB_UT probability, int degreesOfFreedom) {
        if (degreesOfFreedom < 1 || degreesOfFreedom > 100) {
            throw new IllegalArgumentException("degreesOfFreedom must be in range [1,100], inclusive");
        }
        degreesOfFreedom--;
        switch (probability) {
            case Z_90:
                return upperTail[degreesOfFreedom][0];
            case Z_95:
                return upperTail[degreesOfFreedom][1];
            case Z_975:
                return upperTail[degreesOfFreedom][2];
            case Z_99:
                return upperTail[degreesOfFreedom][3];
            case Z_999:
                return upperTail[degreesOfFreedom][4];
        }
        return Double.NaN;
    }
    
      /**
     * a table of lower-tail critical values of the chi-squared distribution
     * where columns are the limit probabilities and rows are the 
     * degrees of freedom (+1 as first row is 0).  
     * The values are a member of the variate X.
     * A quantile is a range within the set of X.
     * One can reject the null hypothesis if the computed test statistic 
        is less than the tabled value.
     */
    static double[][] lowerTail = new double[100][];
    static {                      //Z_10, Z_05, Z_025, Z_01, Z_001;
        lowerTail[0] = new double[]{.016, .004, .001, .000, .000};
        lowerTail[1] = new double[]{.211, .103, .051, .020, .002};
        lowerTail[2] = new double[]{.584, .352, .216, .115, .024};
        lowerTail[3] = new double[]{1.064, .711, .484, .297, .091};
        lowerTail[4] = new double[]{1.610, 1.145, .831, .554, .210};
        lowerTail[5] = new double[]{2.204, 1.635, 1.237, .872, .381};
        lowerTail[6] = new double[]{2.833, 2.167, 1.690, 1.239, .598};
        lowerTail[7] = new double[]{3.490, 2.733, 2.180, 1.646, .857};
        lowerTail[8] = new double[]{4.168, 3.325, 2.700, 2.088, 1.152};
        lowerTail[9] = new double[]{4.865, 3.940, 3.247, 2.558, 1.479};
        lowerTail[10] = new double[]{5.578, 4.575, 3.816, 3.053, 1.834};
        lowerTail[11] = new double[]{6.304, 5.226, 4.404, 3.571, 2.214};
        lowerTail[12] = new double[]{7.042, 5.892, 5.009, 4.107, 2.617};
        lowerTail[13] = new double[]{7.790, 6.571, 5.629, 4.660, 3.041};
        lowerTail[14] = new double[]{8.547, 7.261, 6.262, 5.229, 3.483};
        lowerTail[15] = new double[]{9.312, 7.962, 6.908, 5.812, 3.942};
        lowerTail[16] = new double[]{10.085, 8.672, 7.564, 6.408, 4.416};
        lowerTail[17] = new double[]{10.865, 9.390, 8.231, 7.015, 4.905};
        lowerTail[18] = new double[]{11.651, 10.117, 8.907, 7.633, 5.407};
        lowerTail[19] = new double[]{12.443, 10.851, 9.591, 8.260, 5.921};
        lowerTail[20] = new double[]{13.240, 11.591, 10.283, 8.897, 6.447};
        lowerTail[21] = new double[]{14.041, 12.338, 10.982, 9.542, 6.983};
        lowerTail[22] = new double[]{14.848, 13.091, 11.689, 10.196, 7.529};
        lowerTail[23] = new double[]{15.659, 13.848, 12.401, 10.856, 8.085};
        lowerTail[24] = new double[]{16.473, 14.611, 13.120, 11.524, 8.649};
        lowerTail[25] = new double[]{17.292, 15.379, 13.844, 12.198, 9.222};
        lowerTail[26] = new double[]{18.114, 16.151, 14.573, 12.879, 9.803};
        lowerTail[27] = new double[]{18.939, 16.928, 15.308, 13.565, 10.391};
        lowerTail[28] = new double[]{19.768, 17.708, 16.047, 14.256, 10.986};
        lowerTail[29] = new double[]{20.599, 18.493, 16.791, 14.953, 11.588};
        lowerTail[30] = new double[]{21.434, 19.281, 17.539, 15.655, 12.196};
        lowerTail[31] = new double[]{22.271, 20.072, 18.291, 16.362, 12.811};
        lowerTail[32] = new double[]{23.110, 20.867, 19.047, 17.074, 13.431};
        lowerTail[33] = new double[]{23.952, 21.664, 19.806, 17.789, 14.057};
        lowerTail[34] = new double[]{24.797, 22.465, 20.569, 18.509, 14.688};
        lowerTail[35] = new double[]{25.643, 23.269, 21.336, 19.233, 15.324};
        lowerTail[36] = new double[]{26.492, 24.075, 22.106, 19.960, 15.965};
        lowerTail[37] = new double[]{27.343, 24.884, 22.878, 20.691, 16.611};
        lowerTail[38] = new double[]{28.196, 25.695, 23.654, 21.426, 17.262};
        lowerTail[39] = new double[]{29.051, 26.509, 24.433, 22.164, 17.916};
        lowerTail[40] = new double[]{29.907, 27.326, 25.215, 22.906, 18.575};
        lowerTail[41] = new double[]{30.765, 28.144, 25.999, 23.650, 19.239};
        lowerTail[42] = new double[]{31.625, 28.965, 26.785, 24.398, 19.906};
        lowerTail[43] = new double[]{32.487, 29.787, 27.575, 25.148, 20.576};
        lowerTail[44] = new double[]{33.350, 30.612, 28.366, 25.901, 21.251};
        lowerTail[45] = new double[]{34.215, 31.439, 29.160, 26.657, 21.929};
        lowerTail[46] = new double[]{35.081, 32.268, 29.956, 27.416, 22.610};
        lowerTail[47] = new double[]{35.949, 33.098, 30.755, 28.177, 23.295};
        lowerTail[48] = new double[]{36.818, 33.930, 31.555, 28.941, 23.983};
        lowerTail[49] = new double[]{37.689, 34.764, 32.357, 29.707, 24.674};
        lowerTail[50] = new double[]{38.560, 35.600, 33.162, 30.475, 25.368};
        lowerTail[51] = new double[]{39.433, 36.437, 33.968, 31.246, 26.065};
        lowerTail[52] = new double[]{40.308, 37.276, 34.776, 32.018, 26.765};
        lowerTail[53] = new double[]{41.183, 38.116, 35.586, 32.793, 27.468};
        lowerTail[54] = new double[]{42.060, 38.958, 36.398, 33.570, 28.173};
        lowerTail[55] = new double[]{42.937, 39.801, 37.212, 34.350, 28.881};
        lowerTail[56] = new double[]{43.816, 40.646, 38.027, 35.131, 29.592};
        lowerTail[57] = new double[]{44.696, 41.492, 38.844, 35.913, 30.305};
        lowerTail[58] = new double[]{45.577, 42.339, 39.662, 36.698, 31.020};
        lowerTail[59] = new double[]{46.459, 43.188, 40.482, 37.485, 31.738};
        lowerTail[60] = new double[]{47.342, 44.038, 41.303, 38.273, 32.459};
        lowerTail[61] = new double[]{48.226, 44.889, 42.126, 39.063, 33.181};
        lowerTail[62] = new double[]{49.111, 45.741, 42.950, 39.855, 33.906};
        lowerTail[63] = new double[]{49.996, 46.595, 43.776, 40.649, 34.633};
        lowerTail[64] = new double[]{50.883, 47.450, 44.603, 41.444, 35.362};
        lowerTail[65] = new double[]{51.770, 48.305, 45.431, 42.240, 36.093};
        lowerTail[66] = new double[]{52.659, 49.162, 46.261, 43.038, 36.826};
        lowerTail[67] = new double[]{53.548, 50.020, 47.092, 43.838, 37.561};
        lowerTail[68] = new double[]{54.438, 50.879, 47.924, 44.639, 38.298};
        lowerTail[69] = new double[]{55.329, 51.739, 48.758, 45.442, 39.036};
        lowerTail[70] = new double[]{56.221, 52.600, 49.592, 46.246, 39.777};
        lowerTail[71] = new double[]{57.113, 53.462, 50.428, 47.051, 40.519};
        lowerTail[72] = new double[]{58.006, 54.325, 51.265, 47.858, 41.264};
        lowerTail[73] = new double[]{58.900, 55.189, 52.103, 48.666, 42.010};
        lowerTail[74] = new double[]{59.795, 56.054, 52.942, 49.475, 42.757};
        lowerTail[75] = new double[]{60.690, 56.920, 53.782, 50.286, 43.507};
        lowerTail[76] = new double[]{61.586, 57.786, 54.623, 51.097, 44.258};
        lowerTail[77] = new double[]{62.483, 58.654, 55.466, 51.910, 45.010};
        lowerTail[78] = new double[]{63.380, 59.522, 56.309, 52.725, 45.764};
        lowerTail[79] = new double[]{64.278, 60.391, 57.153, 53.540, 46.520};
        lowerTail[80] = new double[]{65.176, 61.261, 57.998, 54.357, 47.277};
        lowerTail[81] = new double[]{66.076, 62.132, 58.845, 55.174, 48.036};
        lowerTail[82] = new double[]{66.976, 63.004, 59.692, 55.993, 48.796};
        lowerTail[83] = new double[]{67.876, 63.876, 60.540, 56.813, 49.557};
        lowerTail[84] = new double[]{68.777, 64.749, 61.389, 57.634, 50.320};
        lowerTail[85] = new double[]{69.679, 65.623, 62.239, 58.456, 51.085};
        lowerTail[86] = new double[]{70.581, 66.498, 63.089, 59.279, 51.850};
        lowerTail[87] = new double[]{71.484, 67.373, 63.941, 60.103, 52.617};
        lowerTail[88] = new double[]{72.387, 68.249, 64.793, 60.928, 53.386};
        lowerTail[89] = new double[]{73.291, 69.126, 65.647, 61.754, 54.155};
        lowerTail[90] = new double[]{74.196, 70.003, 66.501, 62.581, 54.926};
        lowerTail[91] = new double[]{75.100, 70.882, 67.356, 63.409, 55.698};
        lowerTail[92] = new double[]{76.006, 71.760, 68.211, 64.238, 56.472};
        lowerTail[93] = new double[]{76.912, 72.640, 69.068, 65.068, 57.246};
        lowerTail[94] = new double[]{77.818, 73.520, 69.925, 65.898, 58.022};
        lowerTail[95] = new double[]{78.725, 74.401, 70.783, 66.730, 58.799};
        lowerTail[96] = new double[]{79.633, 75.282, 71.642, 67.562, 59.577};
        lowerTail[97] = new double[]{80.541, 76.164, 72.501, 68.396, 60.356};
        lowerTail[98] = new double[]{81.449, 77.046, 73.361, 69.230, 61.137};
        lowerTail[99] = new double[]{82.358, 77.929, 74.222, 70.065, 61.918};
    }
   
    /**
     * find the lower-tail critical value to compare to a test statistic.
     * One can reject the null hypothesis if the computed test statistic 
        is less than the table value.
     * @param probability the probability to look-up in the table.  
     * e.g. for a one-sided test, this probability is the same as the significance 
     * level alpha.  For a two-sided test, this probability is the same as the significance 
     * level alpha/2.
     * @param degreesOfFreedom
     * @return 
     */
    public static double lowerTailStatistic(PROB_LT probability, int degreesOfFreedom) {
        if (degreesOfFreedom < 1 || degreesOfFreedom > 100) {
            throw new IllegalArgumentException("degreesOfFreedom must be in range [1,100], inclusive");
        }
        degreesOfFreedom--;
        switch (probability) {
            case Z_10:
                return lowerTail[degreesOfFreedom][0];
            case Z_05:
                return lowerTail[degreesOfFreedom][1];
            case Z_025:
                return lowerTail[degreesOfFreedom][2];
            case Z_01:
                return lowerTail[degreesOfFreedom][3];
            case Z_001:
                return lowerTail[degreesOfFreedom][4];
        }
        return Double.NaN;
    }
    
    /**
     * Calculate the p-value of a chi-squared statistic.
     from  
        "Approximating the cumulative chi-squared distribution and
        it's inverse"  1988, Lin, J.T., The Statistician 37, 3-5
        https://www.jstor.org/stable/2348373?seq=1#metadata_info_tab_contents
        
     * @param chisqStat
     * @param degreesOfFreedom
     * @return 
     */
    public static double approxPValueLin(double chisqStat, int degreesOfFreedom) {
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException("degreesOfFreedom must larger than 0");
        }
        double p;
        double a1 = -0.9911;
        double b1 = 0.8055;
        double a2 = -0.6763;
        double b2 = -1.2451;
        
        //EQN (1):

        double z = Math.sqrt(chisqStat) - Math.sqrt(degreesOfFreedom);

        if (z <= 0) {
            p = 1. - 0.5 * Math.exp(b1 * z + a1 * z * z);
        } else {
            p = 0.5 * Math.exp(b2 * z + a2 * z * z);
        }
        
        System.out.printf("Lin: chisq_t==%12.4f, df=%d, => p=%9.4f, 1-p=%9.4f%n",
            chisqStat, degreesOfFreedom, p, (1 - p));
        System.out.flush();
        
        return 1. - p;
    }
    
    /**
     * Calculate the chi-squared statistic given the upper tail p-value.
     from  
        "Approximating the cumulative chi-squared distribution and
        it's inverse"  1988, Lin, J.T., The Statistician 37, 3-5
        https://www.jstor.org/stable/2348373?seq=1#metadata_info_tab_contents
        
     * @param p the p-value OR (1-p)
     * @param degreesOfFreedom
     * @return 
     */
    public static double approxChiSqStatLin(double p, int degreesOfFreedom) {
        if (degreesOfFreedom < 1) {
            throw new IllegalArgumentException("degreesOfFreedom must larger than 0");
        }
        double chisqStat;
        double a1 = -0.9911;
        double b1 = 0.8055;
        double a2 = -0.6763;
        double b2 = -1.2451;
        double c1, c2;
        double z;

        //EQN (2):
        if (p >= 0.5) {
            c1 = -Math.log(2 * (1. - p));
            z = (-b1 + Math.sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1);
        } else {
            c2 = -Math.log(2 * p);
            z = (-b2 - Math.sqrt(b2 * b2 - 4 * a2 * c2)) / (2 * a2);
        }
        chisqStat = Math.pow(z + Math.sqrt(degreesOfFreedom), 2.);

        System.out.printf("Lin: df=%d, p=%9.4f => chisq=%12.4f%n",
            degreesOfFreedom, p, chisqStat);
        System.out.flush();
        
        return chisqStat;
    }
    
    /*
     * Calculate the p-value of a chi-squared statistic.
     * "Exploring How to Simply Approximate the P-value of a Chi-squared Statistic"
     * 2018, Beh, E., Australian Journal of Statistics, Vol. 47 No. 3
     * https://doi.org/10.17713/ajs.v47i3.757
     * @param chisqStat
     * @param degreesOfFreedom
     * @return 
    public static double approxUpperTailPValueBeh(double chisqStat, int degreesOfFreedom) {
        
        double p;
        
        double a =-1.37266; 
        double b = 1.06807; 
        double c = 2.13161; 
        double d =-0.0458;
        
        double numer = Math.sqrt(chisqStat) - (a + b*Math.sqrt(degreesOfFreedom));
        double denom = c - d*Math.sqrt(degreesOfFreedom);
        
        // P(x^2 > X^2) ~  0.1 * ( numer/denom)^2
        p = (numer/denom);
        p *= p;
        p *= 0.1;
        
        System.out.printf("Beh: chisq_t==%12.4f, df=%d, => p=%9.4f, 1-p=%9.4f%n",
            chisqStat, degreesOfFreedom, p, (1 - p));
        System.out.flush();
        
        return p;
    }
    */
    
}
