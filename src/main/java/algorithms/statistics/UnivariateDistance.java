package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import java.util.Arrays;

/**
 * calculates the square of the unbiased distance covariance and correlation 
   between univariate vectors x and y as
   a weighted  distance between the joint characteristic function and 
   the product of marginal distributions following algorithms in
   <pre>
      "A fast algorithm for computing distance correlation"
      2019 Chaudhuri and Hu, Computational Statistics And Data Analysis,
      Volume 135, July 2019, Pages 15-24.
      https://arxiv.org/pdf/1810.11332.pdf
   </pre>
   Runtime is O(n * lg_2(n)) where n is the number of points in x which is
   the same as the number in y

 * TODO: consider implementing Brownian Distance Covariance.
 * TODO: add notes for Hilbert-Schmidt independence measure (HSIC) - Lasso.
 
 * NOTE: can use this within feature screening: 
     Li, R., Zhong, W., and Zhu, L. (2012). 
     Feature screening via distance correlation learning. 
     Journal of the American Statistical Association, 107(499):1129–1139

 * @author nichole
 */
public class UnivariateDistance {
    
    /**
     * calculates the unbiased distance covariance between univariate vectors x and y as
     * "a weighted  distance between the joint characteristic function and 
     * the product of marginal distributions; 
     * it is 0 if and only if two random vectors  and  are independent. 
     * This measure can detect the presence of a dependence structure when the 
     * sample size is large enough."
     * computational cost is O(n log(n)) and memory cost is O(n).
     * This algorithm is ported from the Matlab code in
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri and Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * https://arxiv.org/pdf/1810.11332.pdf
     * 
     * The small corrections to make the covariance unbiased are from 
     * equation 2.5 of “A Statistically And Numerically Efficient Independence 
       Test Based On Random Projections And Distance Covariance”, 
       2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
       See also Appendix C of Shen and Vogelstein, 
       "The Chi-Square Test of Distance Correlation"
     * 
     * Runtime is O(n * lg_2(n)) where n is the number of points in x which is
     * the same as the number in y.
     * 
     * NOTE: redundant points are possible in the rankings because "ties" are handled
     * in the algorithm.  This is one advantage over the similar
     * algorithm of Huo and Szekely (2016).
     *
     @param x sample of univariate observations of a variable
     @param y second sample of univariate observations (can be of another variable).  the length of y must be the
          same as the length of x.
     @return covariance of X and Y along with intermediate data.  
     */
    public static DCov fastDCov(double[] x, double[] y) {

        if (x.length != y.length) {
            throw new IllegalArgumentException("length of x must equal length of y");
        }

        int n = x.length;

        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);
        
        /*
        System.out.printf("X=%s\n", FormatArray.toString(x, "%.3f"));
        System.out.printf("Y=%s\n", FormatArray.toString(y, "%.3f"));
        */

        //"the algorithm essentially consists of two sorting steps. 
        //    First we sort X and calculate ai· for i = 1,...,n. 
        //    Then we sort Y and calculate D and all bi· for i = 1,...,n."
        // x and y sorted:
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(x, y);
        
        /*
        System.out.printf("sortedX=%s\n", FormatArray.toString(x, "%.3f"));
        System.out.printf("sortedY=%s\n", FormatArray.toString(y, "%.3f"));
        System.out.printf("indexes=%s\n", Arrays.toString(indexes));
        System.out.println("//v = [ x  y  x.∗y ];  [n][3]");
        System.out.println("//csumv is cumulative sum of v along columns; [n + 1][3]");
        */

        double[] si = MiscMath0.cumulativeSum(x);
        assert (si.length == n);

        double s = si[n - 1];

        //NOTE: in matlab, the .' is transpose without conjugation
        //NOTE: in matlab, the .* is multiplying same element positions by one another
        //NOTE: in matlab, vector .' * vector is the outer product
        //%a_x is the vector of row sums of distance matrix of x
        //a_x = (−(n −2 ): 2: n ) . ’ . ∗ x + ( s − 2∗ s i ) ;
        //    a_i= (2*i − n)*x_i + (s_n − 2*s_i).
        double[] a_x = new double[n];
        for (int i = 1; i <= n; ++i) {
            a_x[i - 1] = ((2. * i - n) * x[i - 1]) + (s - 2. * si[i - 1]);
        }

        //%Compute Frobenius inner product between the
        //%distance matrices while finding indexes corresponding
        //%to sorted y using merge−sort.
        //%Weight vectors for building the 1st term in Eqn (9)
        //v = [ x y x.∗y ];  [n][3]
        //nw = size( v, 2 );
        //   [ x y x .∗ y ] is column 0 = x, column 1 = y, column 2 = x.*y
        //   matlab .* is multiplying same element positions by one another
        int nw = 3;

        double[][] v = new double[n][nw];
        for (int row = 0; row < n; ++row) {
            v[row][0] = x[row];
            v[row][1] = y[row];
            v[row][2] = x[row] * y[row];
        }

        //% The columns of idx are buffers to store sort indices and output buffer
        //%of merge-sort
        //%we alternate between them and avoid unnecessary copying
        int[][] idx = new int[n][2];
        for (int i = 0; i < n; ++i) {
            idx[i] = new int[2];
            idx[i][0] = i + 1;
        }

        //% iv1, iv2, iv3, iv4 are used to store sum of weights
        //% On output, for 1 ≤ j ≤ n
        // [n][1]
        // iv1(j) = summation_{i<j,y_i<y_j}( 1 )
        // iv2(j) = summation_{i<j,y_i<y_j}( x_i )
        // iv3(j) = summation_{i<j,y_i<y_j}( y_i )
        // iv4(j) = summation_{i<j,y_i<y_j}( x_i * y_i )
        double[] iv1 = new double[n];
        double[] iv2 = new double[n];
        double[] iv3 = new double[n];
        double[] iv4 = new double[n];

        //NOTE: in matlab, the .' is transpose without conjugation
        //NOTE: in matlab, the .* is multiplying same element positions by one another
        //NOTE: in matlab, vector .' * vector is the outer product
        int i, j, k, kf, st1, st2, e1, e2, idx1, idx2, gap, z, zz;
        int[] idx_r = new int[n];
        double[][] csumv = new double[n + 1][nw];
        double[][] tempv = new double[n][nw];
        for (i = 0; i < n; ++i) {
            tempv[i] = new double[nw];
            csumv[i] = new double[nw];
        }
        csumv[n] = new double[nw];
        double[][] tempcs;

        //% The merge sort loop .
        i = 1;
        int r = 1;
        int idx_s = 2;
        while (i < n) {
            gap = 2 * i;
            k = 0;

            //idx_r = idx(:, r);
            for (z = 0; z < idx.length; ++z) {
                idx_r[z] = idx[z][r - 1];
            }

            // csumv = [ zeros(1, nw); 
            // cumsum( v(idx_r, :) ) ] ;
            // dimensions: row0:  nw columns
            //             n rows of nw columns
            //tempv [n][nw];
            for (z = 0; z < n; ++z) {
                for (zz = 0; zz < nw; ++zz) {
                    tempv[z][zz] = v[idx_r[z] - 1][zz];
                }
            }
            //cumulative sum along columns: tempcs [n][nw]
            tempcs = MiscMath0.cumulativeSumAlongColumns(tempv);
            //csumv [n+1][nw];
            Arrays.fill(csumv[0], 0);
            for (z = 0; z < n; ++z) {
                for (zz = 0; zz < nw; ++zz) {
                    csumv[z + 1][zz] = tempcs[z][zz];
                }
            }

            /*
            System.out.printf("\ni=%d, r=%d, gap=%d, idx_s=%d\n", i, r, gap, idx_s);
            System.out.printf("  idx[0]=%s\n", Arrays.toString(idx[0]));
            System.out.printf("  idx[1]=%s\n", Arrays.toString(idx[1]));
            System.out.printf("  iv1=%s\n", FormatArray.toString(iv1, "%.3f"));
            System.out.printf("  iv2=%s\n", FormatArray.toString(iv2, "%.3f"));
            System.out.printf("  iv3=%s\n", FormatArray.toString(iv3, "%.3f"));
            System.out.printf("  iv4=%s\n", FormatArray.toString(iv4, "%.3f"));
            for (z = 0; z < n + 1; ++z) {
                System.out.printf("  csumv[%d]=%s\n", z, FormatArray.toString(csumv[z], "%.3f"));
            }
            */

            //for j = 1:gap:n;
            for (j = 1; j <= n; j += gap) {
                // starts and ends of comparison intervals
                st1 = j;
                e1 = Math.min(st1 + i - 1, n);
                st2 = j + i;
                e2 = Math.min(st2 + i - 1, n);

                /*
                System.out.printf("\n**j=%d, st1=%d, e1=%d, st2=%d, e2=%d\n",
                        j, st1, e1, st2, e2);
                */

                while ((st1 <= e1) && (st2 <= e2)) {
                    k++;
                    //starting indexes from latest sorted order
                    idx1 = idx_r[st1 - 1];
                    idx2 = idx_r[st2 - 1];

                    /*
                    System.out.printf("\n****k=%d, idx1=%d, idx2=%d, st1=%d, e1=%d, st2=%d, e2=%d\n",
                            k, idx1, idx2, st1, e1, st2, e2);
                    System.out.printf("    idx[0]=%s\n", Arrays.toString(idx[0]));
                    System.out.printf("    idx[1]=%s\n", Arrays.toString(idx[1]));
                    System.out.printf("    iv1=%s\n", FormatArray.toString(iv1, "%.3f"));
                    System.out.printf("    iv2=%s\n", FormatArray.toString(iv2, "%.3f"));
                    System.out.printf("    iv3=%s\n", FormatArray.toString(iv3, "%.3f"));
                    System.out.printf("    iv4=%s\n", FormatArray.toString(iv4, "%.3f"));
                    System.out.printf("    (y[%d] >= y[%d]) = %b\n", idx1 - 1, idx2 - 1, (y[idx1 - 1] >= y[idx2 - 1]));
                    */

                    if (y[idx1 - 1] >= y[idx2 - 1]) {
                        //System.out.printf("    *idx[%d][%d] = %d\n", k - 1, idx_s - 1, idx1);
                        idx[k - 1][idx_s - 1] = idx1;
                        st1++;
                    } else {

                        //System.out.printf("    *idx[%d][%d] =%d,  INCR iv1:iv4\n", k - 1, idx_s - 1, idx2);
                        idx[k - 1][idx_s - 1] = idx2;
                        st2++;

                        // an inversion pair is y[idx1], y[idx2]
                        //   which is (csumv[e1 + 1 - 1][1] - csumv[st1 - 1][1]), (csumv[e2 + 1 - 1][1] - csumv[st2 - 1][1])

                        iv1[idx2 - 1] += (e1 - st1 + 1);
                        // use of csumv is similar to use of Summed Area Tables
                        iv2[idx2 - 1] += (csumv[e1 + 1 - 1][0] - csumv[st1 - 1][0]);
                        iv3[idx2 - 1] += (csumv[e1 + 1 - 1][1] - csumv[st1 - 1][1]);
                        iv4[idx2 - 1] += (csumv[e1 + 1 - 1][2] - csumv[st1 - 1][2]);
                    } // end if-else

                    //System.out.printf("    (st1 <= e1) = %b, (st2 <= e2) = %b\n",
                    //        (st1 <= e1), (st2 <= e2));

                } // end while (( st1 <= e1 ) && ( st2 <= e2 ) )

                if (st1 <= e1) {
                    kf = k + e1 - st1 + 1;
                    // note: idx is int[n][2];  idx_r is int[n]
                    //idx( (k+1):kf, s ) = idx_r( st1 : e1, : );
                    int c = st1;
                    for (z = (k + 1); z <= kf; ++z) {
                        idx[z - 1][idx_s - 1] = idx_r[c - 1];
                        //System.out.printf("    *idx[%d][%d]=idx_r[%d] = %d, \n", z - 1, idx_s - 1, c - 1, idx_r[c - 1]);
                        c++;
                    }
                    k = kf;
                } else if (st2 <= e2) {
                    kf = k + e2 - st2 + 1;
                    //idx( ( k+1):kf, s ) = idx_r( st2 : e2, : );
                    int c = st2;
                    for (z = (k + 1); z <= kf; ++z) {
                        idx[z - 1][idx_s - 1] = idx_r[c - 1];
                        //System.out.printf("    *idx[%d][%d]=idx_r[%d] = %d, \n", z - 1, idx_s - 1, c - 1, idx_r[c - 1]);
                        c++;
                    }
                    k = kf;
                }

            } // end for j=

            i = gap;
            r = 3 - r;
            idx_s = 3 - idx_s;
        }

        //% d is the Frobenius inner product of the distance matrices
        // The second term of Eqn (9):
        //covterm = n∗( x − mean(x) ) .’ ∗ ( y − mean(y) );
        //              [1][n] * [n][1] = [1][1]
        double[] mx = MiscMath0.mean(x, 1);
        assert (mx.length == 1);
        double[] my = MiscMath0.mean(y, 1);
        assert (my.length == 1);
        double[] xz = Arrays.copyOf(x, x.length);
        double[] yz = Arrays.copyOf(y, y.length);
        for (z = 0; z < n; ++z) {
            xz[z] -= mx[0];
            yz[z] -= my[0];
        }
        double covtermXY = n * MatrixUtil.innerProduct(xz, yz);

        //v is double[n][nw]; v = [ x y x.∗y ];
        //c1 = iv1 .’ ∗ v (:, 3 );
        //c2 = sum( iv4 );
        //c3 = iv2 .’ ∗ y;
        //c4 = iv3 .’ ∗ x;
        double[] v3 = new double[n];
        double c2 = 0;
        for (z = 0; z < n; ++z) {
            v3[z] = v[z][2]; //x.*y
            c2 += iv4[z];    //sum(inversion pairs of x*y)
        }

        double c1 = MatrixUtil.innerProduct(iv1, v3); //sum(x.*y)
                                                      // c2 is sum(inversion pairs  x dot y)
        double c3 = MatrixUtil.innerProduct(iv2, y);  //(sum(inversion pairs  x)) dot y
        double c4 = MatrixUtil.innerProduct(iv3, x);  //(sum(inversion pairs  y)) dot x

        // d = 4∗( ( c1 + c2 ) − ( c3 + c4 ) ) − 2∗ covterm;
        final double d = (4. * ((c1 + c2) - (c3 + c4))) - 2. * covtermXY;

        double[] b_y = _calcB(y, idx, n, r);

        // extract the indexes for sorting y w.r.t. original indexes
        int[] indexesToSortY = new int[n];
        int c = 0;
        for (z = n; z >= 1; z--) {
            indexesToSortY[c] = indexes[ idx[z - 1][r - 1] - 1 ];
            ++c;
        }

        //NOTE: minor edits to nsq, ncb, and nq to make the covariance
        // unbiased following equation 2.5 of
        // “A Statistically And Numerically Efficient Independence Test Based On 
        // Random Projections And Distance Covariance”, 
        // 2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
        //  See also Appendix C of Shen and ogelstein, 
        //  "The Chi-Square Test of Distance Correlation"
        
        //%covsq equals V^2_n(x, y) the square of the distance covariance
        //%between x and y
        double nsq = n * (n - 3.);//(double)(n*n);
        double ncb = nsq * (n - 2.);//nsq*(double)(n);
        double nq = ncb * (n - 1.);//ncb*(double)(n);
        //Eqn (3):
        //term1 = d / nsq;
        //term2 = 2∗ ( a_x .’ ∗ b_y ) / ncb;
        //term3 = sum( a_x ) ∗ sum( b_y ) / nq;
        //covsq = ( term1 + term3 ) − term2;
        double term1 = d / nsq;
        double term2 = (2. / ncb) * (MatrixUtil.innerProduct(a_x, b_y));
        double a_dot_dot = 0;
        double b_dot_dot = 0;
        for (z = 0; z < n; ++z) {
            a_dot_dot += a_x[z];
            b_dot_dot += b_y[z];
        }

        double term3 = (a_dot_dot * b_dot_dot) / nq;
        double covsq = (term1 + term3) - term2;
        
        /*       
        System.out.printf("iv1=%s\n", FormatArray.toString(iv1, "%.3f"));
        System.out.printf("iv2=%s\n", FormatArray.toString(iv2, "%.3f"));
        System.out.printf("iv3=%s\n", FormatArray.toString(iv3, "%.3f"));
        System.out.printf("iv4=%s\n", FormatArray.toString(iv4, "%.3f"));
        */
        
        /*
        TODO: add test here:
        
        Distance Based Independence Tests in Section 2.3 of paper
        "A Statistically and Numerically Efficient Independence Test Based On
        Random Projections and Distance Covariance"
        2017, Huang and Hu
        
        (n * covsq / term3) > (InverseNormalCDF(1 - (alpha_s/2))^2
             where alpha_s is .gt. 0 and .lt. 0.215
             e.g. alpha_s = 0.05        
        */
                    
        DCov dcov = new DCov();
        dcov.covsq = covsq;
        dcov.d = d;
        dcov.indexesToSortX = indexes;
        dcov.indexesToSortY = indexesToSortY;
        dcov.aDotDot = a_dot_dot;
        dcov.bDotDot = b_dot_dot;
        dcov.ai = a_x;
        dcov.bi = b_y;
        dcov.iv1 = iv1;
        dcov.iv2 = iv2;
        dcov.iv3 = iv3;
        dcov.iv4 = iv4;

        return dcov;
    }
    
    /**
     * calculate the covariance matrix for a using the fast algorithm of
     * <pre>
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri and Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * </pre>
     @param a an mxn matrix of data with a sample of a variable per row.
     * <pre>
     * e.g.  row 0 holds a sample of variable x, row 1 holds a sample of variable y, 
     * etc. so that the array format is [nVariables][nSample]
     *      a[0] = new double[]{x_0, x_1, x_2, ... x_n}
     *      a[1] = new double[]{y_0, y_1, y_2, ... y_n}
     *      ...
     * </pre>
     @return the covariance matrix as a double array of size [a.length][a.length]
     */
    /*public static double[][] multivariateCovariance(double[][] a) {
    
        from paper:
            "For multivariate random variables,
            random projection approach can be adopted (Huang and Huo, 2017), which
            depends on calculation of distance covariance for univariate variables."  

        throw new UnsupportedOperationException("implementation is not yet finished.");
    }*/
      
        /*
         matlab code from 
         "A fast algorithm for computing distance correlation" by Chaudhuri and Hu, 2018
    
        % The merge sort loop .  
        i = 1; r = 1; s = 2;
        while i < n;
           gap = 2∗i;
           k = 0;
           idx_r = idx(:, r);
           csumv = [ zeros(1, nw); cumsum( v(idx_r, :) ) ] ;
           for j = 1:gap:n;
              st1 = j; e1 = min( st1 + i − 1, n );
              st2 = j + i; e2 = min( st2 + i − 1, n );
              while ( st1 <= e1 ) && ( st2 <= e2 );
                 k = k +1;
                 idx1 = idx_r( st1 );
                 idx2 = idx_r( st2 );
                 if y( idx1 ) >= y(idx2);
                    idx(k, s) = idx1;
                    st1 = st1 + 1;
                 else
                    idx(k, s ) = idx2;
                    st2 = st2 + 1;
                    iv1( idx2, 1) = iv1( idx2 ) + e1 − st1 +1;
                    iv2( idx2 ) = iv2( idx2 ) + ( csumv( e1+1, 1) − csumv ( st1, 1 ) ) ;
                    iv3( idx2 ) = iv3( idx2 ) + ( csumv( e1+1, 2) − csumv ( st1, 2 ) ) ;
                    iv4( idx2 ) = iv4( idx2 ) + ( csumv( e1+1, 3) − csumv ( st1, 3 ) ) ;
                 end;
              end;
              if st1 <= e1;
                 kf = k + e1 − st1 + 1;
                 idx( (k+1):kf, s ) = idx_r( st1 : e1, : );
                 k = kf ;
              elseif st2 <=e2;
                 kf = k + e2 − st2 + 1;
                 idx( ( k+1):kf, s ) = idx_r( st2 : e2, : );
                 k = kf;
              end;
           end;
           //% d is the Frobenius inner product of the distance matrices
           covterm = n∗( x − mean(x) ) .’ ∗ ( y − mean(y) );
           c1 = iv1 .’ ∗ v (:, 3 );
           c2 = sum( iv4 );
           c3 = iv2 .’ ∗ y;
           c4 = iv3 .’ ∗ x;
           d = 4∗( ( c1 + c2 ) − ( c3 + c4 ) ) − 2∗ covterm;
           //% b_y is the vector of row sums of distance matrix of y
           ySorted = y ( idx( n : −1: 1, r ) );
           si = cumsum( ySorted );
           s = si( n );
           b_y = zeros( n , 1 );
           b_y( idx( n : −1: 1, r ) ) = (−(n −2 ): 2: n ) .’ .∗ ySorted + ( s − 2∗ s i );
           //%covsq equals V^2_n(x, y) the square of the distance covariance
           //%between x and y
           nsq = n∗n;
           ncb = nsq∗n;
           nq = ncb∗n;
           term1 = d / nsq;
           term2 = 2∗ ( a_x .’ ∗ b_y ) / ncb;
           term3 = sum( a_x ) ∗ sum( b_y ) / nq;
           covsq = ( term1 + term3 ) − term2;
        */  
    
    /**
     * checks the descending sort algorithm for having ported the code from 1-based array indexes to 0-based indexes.
     * "A fast algorithm for computing distance correlation" by Chaudhuri and Hu, 2018
     @param x
     @param y
     @return 2 dimensional array of size double[2][x.length] holding the sorted x and y in each row, respectively
     */
    static double[][] _sortCheck(double[] x, double[] y) {
        
        if (x.length != y.length) {
            throw new IllegalArgumentException("length of x must equal length of y");
        }
        
        int n = x.length;
        
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);
        
        // x and y sorted:
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(x, y);

        //% The columns of idx are buffers to store sort indices and output buffer
        //%of merge-sort
        //%we alternate between them and avoid unnecessary copying
        int[][] idx = new int[n][2];
        for (int i = 0; i < n; ++i) {
            idx[i] = new int[2];
            idx[i][0] = i+1;
        }
        
        //NOTE: in matlab, the .' is transpose without conjugation
        //NOTE: in matlab, the .* is multiplying same element positions by one another
        //NOTE: in matlab, vector .' * vector is the outer product
        
        int i, j, k, kf, st1, st2, e1, e2, idx1, idx2, gap, z, zz;
        int[] idx_r = new int[n];
                                
        //% The merge sort loop .
        i = 1;
        int r = 1; 
        // s = 2 in matlab code, was repurposed use of var s which was an integer type due to integer array x's cumulative sum.
        int idx_s = 2;
        while (i < n) {
           gap = 2*i; 
           k = 0;
           
           //idx_r = idx(:, r);
           for (z = 0; z < idx.length; ++z) {
               idx_r[z] = idx[z][r-1];
           }
           
           for (j = 1; j <= n; j+=gap) {
              st1 = j; 
              e1 = Math.min(st1 + i - 1, n);
              st2 = j + i; 
              e2 = Math.min(st2 + i - 1, n);
             // System.out.printf("i=%d, j=%d, gap=%d, st1=%d, st2=%d, e1=%d, e2=%d\n", i, j, gap, st1, st2, e1, e2);
              while (( st1 <= e1 ) && ( st2 <= e2 ) ) {
                 k++;
                 idx1 = idx_r[st1-1];
                 idx2 = idx_r[st2-1];
                 //System.out.printf("k=%d, idx1=%d, idx2=%d, st1=%d, st2=%d, y[idx1-1]=%.3e, y[idx2-1]=%.3e\n",
                 //        k, idx1, idx2, st1, st2, y[idx1-1], y[idx2-1]);
                 //System.out.printf("compare y[%d] to y[%d]\n", idx1-1, idx2-1);
                 if (y[idx1-1] >= y[idx2-1]) {                 
                    idx[k-1][idx_s-1] = idx1;
                    st1++;
                 } else {
                    idx[k-1][idx_s-1] = idx2;
                    st2++;
                 } // end if-else
              } // end while
              if (st1 <= e1) {
                  kf = k + e1 - st1 + 1;
                  // note: idx is int[n][2];  idx_r is int[n]
                  //idx( (k+1):kf, s ) = idx_r( st1 : e1, : );
                  int c = st1;
                  for (z = (k+1); z <= kf; ++z) {
                      idx[z-1][idx_s-1] = idx_r[c-1];
                      c++;
                  }
                  k = kf;
              } else if (st2 <= e2) {
                  kf = k + e2 - st2 + 1;
                  //idx( ( k+1):kf, s ) = idx_r( st2 : e2, : );
                  int c = st2;
                  for (z = k+1; z <= kf; ++z) {
                      idx[z-1][idx_s-1] = idx_r[c-1];
                      c++;
                  }
                  k = kf;
              }
           } // end for j=
           
           i = gap;
           r = 3 - r; 
           idx_s = 3 - idx_s;  //2;  3-2=1 // 3-1=2
           
        }

       // System.out.printf("indexes=\n%s\n", Arrays.toString(indexes));
       // System.out.printf("idx=\n%s\n", FormatArray.toString(idx, "%d"));
       // System.out.printf("last col used in idx was %d\n", (idx_s == 2) ? /*1-1*/0 : /*2-1*/1);
       // System.out.printf("i=%d, r=%d, idx_s=%d\n", i, r, idx_s);
        //System.out.printf("idx_r=\n%s\n", FormatArray.toString(idx_r, "%d"));
       // System.out.printf("after x=\n%s\n", FormatArray.toString(x, "%.3e"));
       // System.out.printf("after y=\n%s\n", FormatArray.toString(y, "%.3e"));

        // the indexes for the sorting of y after sorted by x:
        // lastCol = (idx_s == 2) ? /*1-1*/0 : /*2-1*/1;
        // idx[*][lastCol]

        double[][] sorted = new double[2][x.length];
        sorted[0] = x;
        sorted[1] = y;
        return sorted;
    }

    private static double[] _calcB(double[] y, int[][] idx, int n, int r) {

        /*
        System.out.printf("\n\n*_calcB\n");
        System.out.printf("idx[0]=%s\n", Arrays.toString(idx[0]));
        System.out.printf("idx[1]=%s\n", Arrays.toString(idx[1]));
        */

        double[] ySorted = new double[n];

        //% b_y is the vector of row sums of distance matrix of y
        // ySorted = y ( idx( n : −1: 1, r ) );
        int c = 0;
        int z;

        for (z = n; z >= 1; z--) {
            ySorted[c] = y[idx[z - 1][r - 1] - 1];
            //System.out.printf("ySorted[%d]=y[idx[%d][%d]-1]=y[%d]=%.3f\n",
            //        c, z - 1, r - 1, idx[z - 1][r - 1] - 1, y[idx[z - 1][r - 1] - 1]);
            c++;
        }

        //si = cumsum( ySorted );
        double[] si = MiscMath0.cumulativeSum(ySorted);
        double s = si[n - 1];

        //b_y = zeros( n , 1 );
        double[] b_y = new double[n];

        //b_y(idx(n : −1: 1, r)) = (−(n−2): 2: n) .’ .∗ ySorted + (s − 2∗ si);
        c = 0;
        double cc = -(n - 2.);
        for (z = n; z >= 1; z--) {
            b_y[idx[z - 1][r - 1] - 1] = (cc * ySorted[c]) + (s - (2. * si[c]));

            /*
            System.out.printf("(s - (2.*si[%d])) = (%.3f - (2. * %.3f)) = %.3f\n",
                    c, s, si[c], (s - (2. * si[c])));
            System.out.printf("b_y[idx[%d][%d]-1] => b_y[%d] = (%.0f * ySorted[%d]) + (s - (2.*si[%d]) = %.3f\n",
                    z - 1, r - 1, idx[z - 1][r - 1] - 1, cc, c, c, (cc * ySorted[c]) + (s - (2. * si[c])));
            */

            c++;
            cc += 2;
        }

        return b_y;
    }

    /**
     * calculates the distance covariance between univariate vectors x and y as
     * "a weighted distance between the joint characteristic function and the
     * product of marginal distributions; it is 0 if and only if two random
     * vectors and are independent. This measure can detect the presence of a
     * dependence structure when the sample size is large enough."
     *
     * This algorithm is ported from the Matlab code in "A fast
     * algorithm for computing distance correlation" 2019 Chaudhuri and Hu,
     * Computational Statistics And Data Analysis, Volume 135, July 2019, Pages
     * 15-24. https://arxiv.org/pdf/1810.11332.pdf
     *
     * Runtime is O(n*log_2(n)) where n is the number of points in x which is
     * the same as the number in y.
     *
     * NOTE: redundant points are possible in the rankings as "ties" are handled
     * in the algorithm. This is one advantage over the similar algorithm of Huo
     * and Szekely (2016).
     *
     * NOTE: that this method follows Algorithm 1 in the paper which stops at
     * the intermediate steps to show the merge steps clearly.
     * The complete algorithm is in the method fastDCov().
     *
     @param x
     @param y
     @return
     */
    static DCov _univariateCovariance(double[] x, double[] y) {

        if (x.length != y.length) {
            throw new IllegalArgumentException("length of x must equal length of y");
        }

        int n = x.length;

        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);

        int[][] idx = new int[2][n];
        idx[0] = new int[n];
        idx[1] = new int[n];
        for (int i = 0; i < n; ++i) {
            idx[0][i] = i + 1;
        }

        int[] idx_r = new int[n];

        double[] csumT = new double[n + 1];
        double[] d = new double[n];

        //the merge sort is descending order
        int gap, k, kf, z, j, st1, e1, st2, e2, idx1, idx2;
        int i = 1;
        int r = 1;
        int idx_s = 2;
        while (i < n) {
            gap = 2 * i;
            k = 0;

            System.arraycopy(idx[r - 1], 0, idx_r, 0, n);

            //csumT = cusum(y[idx r]);
            //csumT = (0, csumT );
            csumT[0] = 0;
            for (z = 0; z < n; ++z) {
                csumT[z + 1] = y[idx_r[z] - 1];
            }
            for (z = 2; z <= n; ++z) {
                csumT[z] += csumT[z - 1];
            }

            j = 1;

            while (j < n) {

                st1 = j;
                e1 = Math.min(st1 + i - 1, n);
                st2 = j + i;
                e2 = Math.min(st2 + i - 1, n);
                while ((st1 <= e1) && (st2 <= e2)) {
                    k++;
                    idx1 = idx_r[st1 - 1];
                    idx2 = idx_r[st2 - 1];
                    //NOTE: we're comparing x here to compare to the paper's results in a unit test.
                    // in fastDCov() we compare y because x is already sorted
                    if (x[idx1 - 1] >= x[idx2 - 1]) {
                        idx[idx_s - 1][k - 1] = idx1;
                        st1++;
                    } else {
                        idx[idx_s - 1][k - 1] = idx2;
                        st2++;
                        // d[] here is similar to iv3 in method fastDCov
                        //   which is the paper's Appendix Matlab code.
                        // iv3(j) = summation_{i<j,y_i<y_j}( y_i )
                        d[idx2 - 1] += (csumT[e1 + 1 - 1] - csumT[st1 - 1]);
                    } // end if-else
                } // end while
                if (st1 <= e1) {
                    kf = k + e1 - st1 + 1;
                    int c = st1;
                    for (z = (k + 1); z <= kf; ++z) {
                        idx[idx_s - 1][z - 1] = idx_r[c - 1];
                        c++;
                    }
                    k = kf;
                } else if (st2 <= e2) {
                    kf = k + e2 - st2 + 1;
                    int c = st2;
                    for (z = (k + 1); z <= kf; ++z) {
                        idx[idx_s - 1][z - 1] = idx_r[c - 1];
                        c++;
                    }
                    k = kf;
                }
                j += gap;

            } // end while j

            i = gap;
            r = 3 - r;
            idx_s = 3 - idx_s;
        } // end while i

        int[] indexes = new int[n];
        double[] sortedD = new double[n];
        for (z = 0; z < n; ++z) {
            indexes[z] = idx[r - 1][z] - 1;
            sortedD[z] = d[indexes[z]];
        }

        DCov dcov = new DCov();
        dcov.iv3 = sortedD;
        dcov.indexesToSortX = indexes;

        return dcov;
    }

    /**
     * calculates the distance covariance between univariate vectors x and y as
     * "a weighted  distance between the joint characteristic function and 
     * the product of marginal distributions; 
     * it is 0 if and only if two random vectors  and  are independent. 
     * This measure can detect the presence of a dependence structure when the 
     * sample size is large enough."
     * 
     * This algorithm is ported from the Matlab code in
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri and Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * 
     * Runtime is O(n * lg_2(n)) where n is the number of points in x which is
     * the same as the number in y.
     * 
     * NOTE: redundant points are possible in the rankings as "ties" are handled
     * in the algorithm.  This is one advantage over the similar
     * algorithm of Huo and Szekely (2016).
     * 
     * runtime is  O(n*log_2(n)).
     * 
     @param x sample of univariate observations of a variable
     @param y second sample of univariate observations (can be of another variable)
     @return the distance correlation of X and Y.  also returns intermediate data.
     * NOTE that the distance correlation is exactly zero if and only if two 
     * random variables are independent.
     */
    public static DCor fastDCor(double[] x, double[] y) {
        DCor dcor = new DCor();
        double tol = 1e-15;
        dcor.covXXSq = fastDCov(x, x);
        if (dcor.covXXSq.covsq < tol) {
            dcor.corSq = 0;
            System.err.println("cov(XX) is 0, so correlation is 0");
            return dcor;
        }
        dcor.covYYSq = fastDCov(y, y);
        if (dcor.covYYSq.covsq < tol) {
            dcor.corSq = 0;
            System.err.println("cov(YY) is 0, so correlation is 0");
            return dcor;
        }
        dcor.covXYSq = fastDCov(x, y);
        if (dcor.covXYSq.covsq < tol) {
            dcor.corSq = 0;
            System.err.println("cov(XY) is 0, so correlation is 0");
            return dcor;
        }
        dcor.corSq = dcor.covXYSq.covsq/Math.sqrt(dcor.covXXSq.covsq * dcor.covYYSq.covsq);
        return dcor;
    }
  
    /**
     *
     */
    public static class DCov {

        /**
         *
         */
        public double covsq;

        /**
         *
         */
        public double d;
        int[] indexesToSortX;
        int[] indexesToSortY;
        
        /**
         * this is a.. of eqn (2) of the paper.
         * a_.. = summation_{i:1,n}( a_i )
         * where a_i = summation_{j:1,n}( a_i_j )
         * where a_i_j = |x_i - x_j|_p
         * and |·|_p is the euclidean norm in real space with dimension p
         */
        public double aDotDot;
        
        /**
         * this is a_i of eqn (2) of the paper.
         * where a_i = summation_{j:1,n}( a_i_j )
         * where a_i_j = |x_i - x_j|_p
         * where p is the number of dimensions of vector x
         * and |·|_p is the euclidean norm in real space with dimension p
         */
        public double[] ai;
        
        /**
         <pre>
         this is b.. of eqn (2) of the paper.
         b_.. = summation_{i:1,n}( b_i )
         where b_i = summation_{j:1,n}( b_i_j )
         where b_i_j = |y_i - y_j|_q
         where q is the number of dimensions of vector y
         and |·|_q is the euclidean norm in real space with dimension q
         </pre>
         */
        public double bDotDot;
        
        /**
         <pre>
         this is b_i of eqn (2) of the paper.
         where b_i = summation_{j:1,n}( b_i_j )
         where b_i_j = |y_i - y_j|_q
         where q is the number of dimensions of vector y
         and |·|_q is the euclidean norm in real space with dimension q
         </pre>
         */
        public double[] bi;
        
        /**
        iv1(j) = summation_{i .lt. j, y_i .lt. y_j}( 1 )
        */
        public double[] iv1;
        
        /**
        iv2(j) = summation_{i .lt. j, y_i .lt. y_j}( x_i )
        */
        public double[] iv2;
        
        /**
        iv3(j) = summation_{i .lt. j, y_i .lt. y_j}( y_i )
        */
        public double[] iv3;
        
        /**
        iv4(j) = summation_{i .lt. j, y_i .lt. y_j}( x_i * y_i )
        */
        public double[] iv4;
                
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("d=").append(d).append("\n");
            sb.append("covsq=").append(covsq).append("\n");
            if (indexesToSortX != null) {
                sb.append("indexesToSortX=").append(Arrays.toString(indexesToSortX)).append("\n");
            }
            if (indexesToSortY != null) {
                sb.append("indexesToSortY=").append(Arrays.toString(indexesToSortY)).append("\n");
            }
            if (ai != null) {
                sb.append("ai=").append(FormatArray.toString(ai, "%.3f")).append("\n");
                sb.append("a..=").append(aDotDot).append("\n");
            }
            if (bi != null) {
                sb.append("bi=").append(FormatArray.toString(bi, "%.3f")).append("\n");
                sb.append("b..=").append(bDotDot).append("\n");
            }
            if (iv1 != null) {
                sb.append("iv1=").append(Arrays.toString(iv1)).append("\n");
            }
            if (iv2 != null) {
                sb.append("iv2=").append(Arrays.toString(iv2)).append("\n");
            }
            if (iv3 != null) {
                sb.append("iv3=").append(Arrays.toString(iv3)).append("\n");
            }
            if (iv4 != null) {
                sb.append("iv4=").append(Arrays.toString(iv4)).append("\n");
            }
                        
            return sb.toString();
        }
        
        /**
         *
         @return
         */
        public String toString2() {
            StringBuilder sb = new StringBuilder();
            sb.append("d=").append(d).append("\n");
            sb.append("covsq=").append(covsq).append("\n");
            if (ai != null) {
                sb.append("a..=").append(aDotDot).append("\n");
            }
            if (bi != null) {
                sb.append("b..=").append(bDotDot).append("\n");
            }
                        
            return sb.toString();
        }
    }
  
    /**
     *
     */
    public static class DCor {

        /**
         *
         */
        public DCov covXYSq;

        /**
         *
         */
        public DCov covXXSq;

        /**
         *
         */
        public DCov covYYSq;

        /**
         *
         */
        public double corSq;
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("corSq=").append(corSq).append("\n");
            if (covXYSq != null) {
                sb.append("covXYSq: ").append(covXYSq.toString());
                sb.append("\n");
            }
            if (covXXSq != null) {
                sb.append("covXXSq: ").append(covXXSq.toString());
                sb.append("\n");
            }
            if (covYYSq != null) {
                sb.append("covYYSq: ").append(covYYSq.toString());
                sb.append("\n");
            }            
            return sb.toString();
        }
    }
}
