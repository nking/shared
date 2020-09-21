package algorithms.correlation;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import java.util.Arrays;

/**
 * calculates the distance covariance between univariate vectors x and y as
     "a weighted  distance between the joint characteristic function and 
     the product of marginal distributions.
    
     * The method covariance() is a port of the Matlab code from
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * https://arxiv.org/pdf/1810.11332.pdf
     * 
     * Runtime is O(n * lg_2(n)) where n is the number of points in x which is
     * the same as the number in y
 * 
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
     * calculates the distance covariance between univariate vectors x and y as
     * "a weighted  distance between the joint characteristic function and 
     * the product of marginal distributions; 
     * it is 0 if and only if two random vectors  and  are independent. 
     * This measure can detect the presence of a dependence structure when the 
     * sample size is large enough."
     * 
     * This algorithm is an implementation/port of the Matlab code from
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
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
     * @param x sample of univariate observations of a variable
     * @param y second sample of univariate observations (can be of another variable)
     * @return 
     */
    public static DCov fastDcov(double[] x, double[] y) {
                  
        if (x.length != y.length) {
            throw new IllegalArgumentException("length of x must equal length of y");
        }
        
        int n = x.length;
        
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);
        
        //"the algorithm essentially consists of two sorting steps. 
        //    First we sort X and calculate ai· for i = 1,...,n. 
        //    Then we sort Y and calculate D and all bi· for i = 1,...,n."
                
        // x and y sorted:
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(x, y);
        
        double[] si = MiscMath0.cumulativeSum(x);
        assert(si.length == n);
        
        double s = si[n - 1];
        
        //NOTE: in matlab, the .' is transpose without conjugation
        //NOTE: in matlab, the .* is multiplying same element positions by one another
        //NOTE: in matlab, vector .' * vector is the outer product
        
        //%a_x is the vector of row sums of distance matrix of x
        //a_x = (−(n −2 ): 2: n ) . ’ . ∗ x + ( s − 2∗ s i ) ;
        //    a_i= (2*i − n)*x_i + (s_n − 2*s_i).
        double[] a_x = new double[n];
        for (int i = 1; i <= n; ++i) {
            a_x[i-1] = ((2.*i - n) * x[i-1]) + (s - 2.*si[i-1]);
        }
        
        //%Compute Frobenius inner product between the
        //%distance matrices while finding indexes corresponding
        //%to sorted y using merge−sort.
        
        //%Weight vectors for building the 1st term in Eqn (9)
        //v = [ x y x.∗y ];
        //nw = size( v, 2 );
        
        //   [ x y x .∗ y ] is column 0 = x, column 1 = y, column 2 = x.*y
        //   matlab .* is multiplying same element positions by one another
        
        int nw = 3;
        
        double[][] v = new double[n][nw];
        for (int row = 0; row < n; ++row) {
            v[row][0] = x[row];
            v[row][1] = y[row];
            v[row][2] = x[row]*y[row];
        }
        
        //% The columns of idx are buffers to store sort indices and output buffer
        //%of merge-sort
        //%we alternate between them and avoid uneccessary copying
        int[][] idx = new int[n][2];
        for (int i = 0; i < n; ++i) {
            idx[i] = new int[2];
            idx[i][0] = i+1;
        }
        
        //% iv1, iv2, iv3, iv4 are used to store sum of weights
        //% On output, for 1 ≤ j ≤ n
        // [n][1]
        int[] iv1 = new int[n];
        int[] iv2 = new int[n];
        int[] iv3 = new int[n];
        int[] iv4 = new int[n];
        
        //NOTE: in matlab, the .' is transpose without conjugation
        //NOTE: in matlab, the .* is multiplying same element positions by one another
        //NOTE: in matlab, vector .' * vector is the outer product
        
        int i, j, k, kf, st1, st2, e1, e2, idx1, idx2, gap, z, zz;
        int[] idx_r = new int[n];
        double[][] csumv = new double[n+1][nw];
        double[][] tempv = new double[n][nw];
        for (i = 0; i < n; ++i){
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
           gap = 2*i; 
           k = 0;
           
           //idx_r = idx(:, r);
           for (z = 0; z < idx.length; ++z) {
               idx_r[z] = idx[z][r-1];
           }
           
           // csumv = [ zeros(1, nw); 
           // cumsum( v(idx_r, :) ) ] ;
           
           // dimensions: row0:  nw columns
           //             n rows of nw columns
           for (z = 0; z < n; ++z) {
               for (zz = 0; zz < nw; ++zz) {
                   tempv[z][zz] = v[ idx_r[z] -1][zz];
               }
           }
           //[n][nw]
           tempcs = MiscMath0.cumulativeSumAlongColumns(tempv);
           Arrays.fill(csumv[0], 0);
           for (z = 0; z < n; ++z) {
               for (zz = 0; zz < nw; ++zz) {
                   csumv[z + 1][zz] = tempcs[z][zz];
               }
           }
           
           //for j = 1:gap:n;
           for (j = 1; j <= n; j+=gap) {
              st1 = j;
              e1 = Math.min(st1 + i - 1, n);
              st2 = j + i;
              e2 = Math.min(st2 + i - 1, n);
              while (( st1 <= e1 ) && ( st2 <= e2 ) ) {
                 k++;
                 idx1 = idx_r[st1-1];
                 idx2 = idx_r[st2-1];
                 if (y[idx1-1] >= y[idx2-1]) {                 
                    idx[k-1][idx_s-1] = idx1;
                    st1++;
                 } else {
                    idx[k-1][idx_s-1] = idx2;
                    st2++;
                    iv1[idx2-1] += (e1 - st1 + 1);
                    iv2[idx2-1] += (csumv[e1+1-1][0] - csumv[st1-1][0]);
                    iv3[idx2-1] += (csumv[e1+1-1][1] - csumv[st1-1][1]);
                    iv4[idx2-1] += (csumv[e1+1-1][2] - csumv[st1-1][2]);
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
                  for (z = (k+1); z <= kf; ++z) {
                      idx[z-1][idx_s-1] = idx_r[c-1];
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
       assert(mx.length == 1);
       double[] my = MiscMath0.mean(y, 1);
       assert(my.length == 1);
       double[] xz = Arrays.copyOf(x, x.length);
       double[] yz = Arrays.copyOf(y, y.length);
       for (z = 0; z < n; ++z) {
           xz[z] -= mx[0];
           yz[z] -= my[0];
       }
       double covtermXY = n * MatrixUtil.dot(xz, yz);
       
       //v is double[n][nw]; v = [ x y x.∗y ];
       //c1 = iv1 .’ ∗ v (:, 3 );
       //c2 = sum( iv4 );
       //c3 = iv2 .’ ∗ y;
       //c4 = iv3 .’ ∗ x;
       
       double[] v3 = new double[n];
       double c2 = 0;
       for (z = 0; z < n; ++z) {
           v3[z] = v[z][2];
           c2 += iv4[z];
       }
       
       double c1 = MatrixUtil.dot(iv1, v3);
       double c3 = MatrixUtil.dot(iv2, y);
       double c4 = MatrixUtil.dot(iv3, x);
       
       // d = 4∗( ( c1 + c2 ) − ( c3 + c4 ) ) − 2∗ covterm;
       double d = (4.*( (c1 + c2) - (c3 + c4))) - 2.*covtermXY;
       
       double[] b_y = _calcB(y, idx, n, r);
      
//System.out.printf("a_x=%s\n", Arrays.toString(a_x));
//System.out.printf("b_y=%s\n", Arrays.toString(b_y));

       //NOTE: minor edits to nsq, ncb, and nq following equation 2.5 of
       // “A Statistically And Numerically Efficient Independence Test Based On 
       // Random Projections And Distance Covariance”, 
       // 2017, Cheng Huang, And Xiaoming Huo, Annals of Statistics
       
       //%covsq equals V^2_n(x, y) the square of the distance covariance
       //%between x and y
       double nsq = n*(n-3.);//(double)(n*n);
       double ncb = nsq*(n-2.);//nsq*(double)(n);
       double nq = ncb*(n-1.);//ncb*(double)(n);
       //Eqn (3):
       //term1 = d / nsq;
       //term2 = 2∗ ( a_x .’ ∗ b_y ) / ncb;
       //term3 = sum( a_x ) ∗ sum( b_y ) / nq;
       //covsq = ( term1 + term3 ) − term2;
       double term1 = d / nsq;
       double term2 = (2./ncb) * (MatrixUtil.dot(a_x, b_y));
       double term3A = 0;
       double term3B = 0;
       for (z = 0; z < n; ++z) {
           term3A += a_x[z];
           term3B += b_y[z];
       }
       double term3 = (term3A * term3B) / nq;
       double covsq =  (term1 + term3) - term2;
           
       DCov dcov = new DCov();
       dcov.covsq = covsq;
       dcov.d = d;
       dcov.indexes = indexes;
       dcov.sortedX = x;
       dcov.sortedY = y;
       
       return dcov;
    }
    
    /**
     * calculate the covariance matrix for a using the fast algorithm of
     * <pre>
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * </pre>
     * @param a an mxn matrix of data with a sample of a variable per row.
     * <pre>
     * e.g.  row 0 holds a sample of variable x, row 1 holds a sample of variable y, 
     * etc. so that the array format is [nVariables][nSample]
     *      a[0] = new double[]{x_0, x_1, x_2, ... x_n}
     *      a[1] = new double[]{y_0, y_1, y_2, ... y_n}
     *      ...
     * </pre>
     * @return the covariance matrix as a double array of size [a.length][a.length]
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
     * checks the sort algorithm for having ported the code from 1-based array indexes to 0-based indexes.
     * "A fast algorithm for computing distance correlation" by Chaudhuri and Hu, 2018
     * @param x
     * @param y
     * @return 2 dimensional array of size double[2][x.length] holding the sorted x and y in each row, respectively
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
         
        int nw = 3;
        
        double[][] v = new double[n][nw];
        v[0] = Arrays.copyOf(x, n);
        v[1] = Arrays.copyOf(y, n);
        v[2] = Arrays.copyOf(x, n);
        for (int row = 0; row < n; ++row) {
            v[2][row] *= y[row];
        }
        
        //% The columns of idx are buffers to store sort indices and output buffer
        //%of merge-sort
        //%we alternate between them and avoid uneccessary copying
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
              while (( st1 <= e1 ) && ( st2 <= e2 ) ) {
                 k++;
                 idx1 = idx_r[st1-1];
                 idx2 = idx_r[st2-1];
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
        
        double[][] sorted = new double[2][x.length];
        sorted[0] = x;
        sorted[1] = y;
        return sorted;
    }

    private static double[] _calcB(double[] y, int[][] idx, int n, int r) {
        
       double[] ySorted = new double[n];

       //% b_y is the vector of row sums of distance matrix of y
       // ySorted = y ( idx( n : −1: 1, r ) );
       int c = 0;
       int z;
       for (z = n; z >=1; z--) {
           ySorted[c] = y[ idx[z-1][r-1] -1];
           c++;
       }
       
       //si = cumsum( ySorted );
       double[] si = MiscMath0.cumulativeSum(ySorted);
       double s = si[n-1];

       //b_y = zeros( n , 1 );
       double[] b_y = new double[n];

       //b_y(idx(n : −1: 1, r)) = (−(n−2): 2: n) .’ .∗ ySorted + (s − 2∗ si);
       c = 0;
       double cc = -(n-2.);
       for (z = n; z >=1; z--) {
           b_y[ idx[z-1][r-1] -1] = (cc * ySorted[c]) + (s - (2.*si[c]));
           c++;
           cc += 2;
       }

       return b_y;
    }
    
    /**
     * calculates the distance covariance between univariate vectors x and y as
     * "a weighted  distance between the joint characteristic function and 
     * the product of marginal distributions; 
     * it is 0 if and only if two random vectors  and  are independent. 
     * This measure can detect the presence of a dependence structure when the 
     * sample size is large enough."
     * 
     * This algorithm is an implementation/port of the Matlab code from
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
     * Volume 135, July 2019, Pages 15-24.
     * 
     * Runtime is O(n*log_2(n)) where n is the number of points in x which is
     * the same as the number in y.
     * 
     * NOTE: redundant points are possible in the rankings as "ties" are handled
     * in the algorithm.  This is one advantage over the similar
     * algorithm of Huo and Szekely (2016).
     * 
     * NOTE: that this method follows Algorithm 1 in the paper which stops at
     * the intermediate steps to show the merge steps clearly.
     * The complete algorithm is present as fastDcov.
     * 
     * @param x
     * @param y
     * @return 
     */
    public static DCov _univariateCovariance(double[] x, double[] y) {

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
                    if (x[idx1 - 1] >= x[idx2 - 1]) {
                        idx[idx_s - 1][k - 1] = idx1;
                        st1++;
                    } else {
                        idx[idx_s - 1][k - 1] = idx2;
                        st2++;
                        // similar to iv3 in method univariateCovariance2):
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

        double[] sortedX = new double[n];
        double[] sortedY = new double[n];
        int[] indexes = new int[n];
        double[] sortedD = new double[n];
        for (z = 0; z < n; ++z) {
            indexes[z] = idx[r - 1][z] - 1;
            sortedX[z] = x[indexes[z]];
            sortedY[z] = y[indexes[z]];
            sortedD[z] = d[indexes[z]];
        }

        DCov dcov = new DCov();
        dcov.dcov = sortedD;
        dcov.indexes = indexes;
        dcov.sortedX = sortedX;
        dcov.sortedY = sortedY;

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
     * This algorithm is an implementation/port of the Matlab code from
     * "A fast algorithm for computing distance correlation"
     * 2019 Chaudhuri & Hu, Computational Statistics And Data Analysis,
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
     * @param x sample of univariate observations of a variable
     * @param y second sample of univariate observations (can be of another variable)
     * @return 
     */
    public static DCor fastDcor(double[] x, double[] y) {
        DCor dcor = new DCor();
        double tol = 1e-15;
        dcor.covXXSq = fastDcov(x, x);
        if (dcor.covXXSq.covsq < tol) {
            return dcor;
        }
        dcor.covYYSq = fastDcov(y, y);
        if (dcor.covYYSq.covsq < tol) {
            return dcor;
        }
        dcor.covXYSq = fastDcov(x, y);
        dcor.corSq = dcor.covXYSq.covsq/Math.sqrt(dcor.covXXSq.covsq * dcor.covYYSq.covsq);
        return dcor;
    }
  
    public static class DCov {
        public double covsq;
        public double d;
        int[] indexes;
        double[] sortedX;
        double[] sortedY;
        double[] dcov;
        
        public double[] getDCov() {
            return dcov;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            if (indexes != null) {
                sb.append("indexes=").append(Arrays.toString(indexes)).append("\n");
            }
            if (sortedX != null) {
                sb.append("sortedX=").append(Arrays.toString(sortedX)).append("\n");
            }
            if (sortedY != null) {
                sb.append("sortedY=").append(Arrays.toString(sortedY)).append("\n");
            }
            if (dcov != null) {
                sb.append("dcov=").append(Arrays.toString(dcov)).append("\n");
            }
            sb.append("d=").append(d).append("\n");
            sb.append("covsq=").append(covsq).append("\n");
            
            return sb.toString();
        }
    }
  
    public static class DCor {
        public DCov covXYSq;
        public DCov covXXSq;
        public DCov covYYSq;
        public double corSq;
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("corSq=").append(corSq).append("\n");
            if (covXYSq != null) {
                sb.append("covXYSq: ").append(covXYSq.toString());
            }
            if (covXXSq != null) {
                sb.append("covXXSq: ").append(covXXSq.toString());
            }
            if (covYYSq != null) {
                sb.append("covYYSq: ").append(covYYSq.toString());
            }            
            return sb.toString();
        }
    }
}
