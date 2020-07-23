package algorithms.correlation;

import algorithms.misc.MiscMath0;
import algorithms.misc.MiscSorter;
import java.util.Arrays;

/**
 *
 * implementation of Chaudhuri & Hu 2019.
 * 
 * TODO: consider implementing Brownian Distance Covariance.
 * TODO: add notes for Hilbert-Schmidt independence measure (HSIC) - Lasso.
 * 
 * @author nichole
 */
public class Distance {
    
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
     * @param x
     * @param y
     * @return 
     */
    public static DCov covariance(double[] x, double[] y) {
                  
        if (x.length != y.length) {
            throw new IllegalArgumentException("length of x must equal length of y");
        }
        
        int n = x.length;
        
        x = Arrays.copyOf(x, x.length);
        y = Arrays.copyOf(y, y.length);
        
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
        double[] b_y = new double[n];
        
        //%Compute Frobenius inner product between the
        //%distance matrices while finding indexes corresponding
        //%to sorted y using merge−sort.
        
        //%Weight vectors
        //v = [ x y x .∗ y ];
        //nw = size( v, 2 );
        
        //   [ x y x .∗ y ] is column 0 = x, column 1 = y, column 2 = x.*y
        //   matlab .* is multiplying same element positions by one another
        
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
        
        //% iv1, iv2, iv3, iv4 are used to store sum of weights
        //% On output, for 1 ≤ j ≤ n
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
        
        double[][] covterm = new double[n][n];
        double[][] c1 = new double[n][n];
        double c2;
        double[][] c3 = new double[n][n];
        double[][] c4 = new double[n][n];
        double[][] d = new double[n][n];
        for (i = 0; i < n; ++i){
            covterm[i] = new double[n];
            c1[i] = new double[n];
            c3[i] = new double[n];
            c4[i] = new double[n];
            d[i] = new double[n];
        } 
                               
        //% The merge sort loop .
        i = 1;
        int r = 1; 
        int idx_s = 2;
        while (i < n) {
           gap = 2*i; 
           k = 0;
           
           //idx_r = idx(:, r);
           //NOTE: it looks like idx_r will never contain a 0 as that is out-of-bounds index for array v
           for (z = 0; z < idx.length; ++z) {
               idx_r[z] = idx[z][r-1];
           }
           
           // csumv = [ zeros(1, nw); cumsum( v(idx_r, :) ) ] ;
           
           // csumv = [ zeros(1, nw); cumsum( v(idx_r, :) ) ] ;
           // dimensions: row0:  nw columns
           //             n rows of nw columns
           for (z = 0; z < n; ++z) {
               for (zz = 0; zz < nw; ++zz) {
                   tempv[z][zz] = v[ idx_r[z] -1][zz];
               }
           }
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
                    iv1[idx2-1] += e1 - st1 + 1;
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
        
System.out.println("iv1=" + Arrays.toString(iv1)); 
System.out.println("iv2=" + Arrays.toString(iv2)); 
System.out.println("iv3=" + Arrays.toString(iv3)); 
System.out.println("iv4=" + Arrays.toString(iv4)); 
        
       //% d is the Frobenius inner product of the distance matrices
       //covterm = n∗( x − mean(x) ) .’ ∗ ( y − mean(y) );
       //   vector .' * is outer product:  nX1 * 1Xn = nxn
       double[] mx = MiscMath0.mean(x, 1);
       assert(mx.length == 1);
       double[] my = MiscMath0.mean(y, 1);
       assert(my.length == 1);
       for (z = 0; z < n; ++z) {
           for (zz = 0; zz < n; ++zz) {
               covterm[z][zz] = n * ( (x[z] - mx[0]) * (y[zz] - my[0]) );
           }
       }

       //NOTE v is double[n][nw];
       //c1 = iv1 .’ ∗ v (:, 3 );
       //c2 = sum( iv4 );
       //c3 = iv2 .’ ∗ y;
       //c4 = iv3 .’ ∗ x;
       c2 = 0;
       for (z = 0; z < n; ++z) {
           c2 += iv4[z];
           for (zz = 0; zz < n; ++zz) {
               c1[z][zz] = iv1[z] * v[zz][2];
               c3[z][zz] = iv2[z] * y[zz];
               c4[z][zz] = iv3[z] * x[zz];
           }
       }
       // d = 4∗( ( c1 + c2 ) − ( c3 + c4 ) ) − 2∗ covterm;
       // c1, c3, c4, covterm are nxn
       // c2 is 1x1
       for (z = 0; z < n; ++z) {
           for (zz = 0; zz < n; ++zz) {
               d[z][zz] = (4.*((c1[z][zz] + c2) - (c3[z][zz] + c4[z][zz]))) - 2.*covterm[z][zz];
           }
       }
       
       double[] ySorted = new double[n];

       //% b_y is the vector of row sums of distance matrix of y
       // ySorted = y ( idx( n : −1: 1, r ) );
       int c = 0;
       for (z = n; z >=1; z--) {
           ySorted[c] = y[ idx[z-1][r-1] -1];
           c++;
       }

       //si = cumsum( ySorted );
       si = MiscMath0.cumulativeSum(ySorted);
       s = si[n-1];

       //b_y = zeros( n , 1 );
       Arrays.fill(b_y, 0);

       //b_y(idx(n : −1: 1, r)) = (−(n−2): 2: n) .’ .∗ ySorted + (s − 2∗ si);
       c = 0;
       double cc = -(n-2.);
System.out.printf("cc=%.4f", cc);    
       for (z = n; z >=1; z--) {
           b_y[ idx[z-1][r-1] -1] = (cc * ySorted[c]) + (s - (2.*si[c]));
           c++;
           cc += 2;
       }
System.out.printf("to %.4f\n", cc); 

System.out.printf("b_y=%s\n", Arrays.toString(b_y));

       //%covsq equals V^2_n(x, y) the square of the distance covariance
       //%between x and y
       double nsq = (double)(n*n);
       double ncb = (double)(nsq*n);
       double nq = (double)(ncb*n);
       //term1 = d / nsq;
       //term2 = 2∗ ( a_x .’ ∗ b_y ) / ncb;
       //term3 = sum( a_x ) ∗ sum( b_y ) / nq;
       //covsq = ( term1 + term3 ) − term2;
       double[][] term1 = new double[n][];
       double[][] term2 = new double[n][n];
       double term3A = 0;
       double term3B = 0;
       for (z = 0; z < n; ++z) {
           term1[z] = Arrays.copyOf(d[z], n);
           term2[z] = new double[n];
           term3A += a_x[z];
           term3B += b_y[z];
           for (zz = 0; zz < n; ++zz) {
               term1[z][zz] /= nsq;
               term2[z][zz] = (2./ncb) * (a_x[z] * b_y[zz]);
           }
       }
       double term3 = (term3A * term3B) / nq;
       double[][] covsq = new double[n][n];
       for (z = 0; z < n; ++z) {
           for (zz = 0; zz < n; ++zz) {
               covsq[z][zz] = (term1[z][zz] + term3) - term2[z][zz];
           }
       }
       
       // for debugging only:
       System.out.println("d=");
       for (z = 0; z < n; ++z) {
           for (zz = 0; zz < n; ++zz) {
               System.out.printf("  %.3f", d[z][zz]);
           }
           System.out.printf("\n");
       }
       
       DCov dcov = new DCov();
       dcov.covsq = covsq;
       dcov.d = d;
       dcov.indexes = indexes;
       dcov.sortedX = x;
       dcov.sortedX = y;
           
       return dcov;
    }
      
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
    
    public static class DCov {
        double[][] covsq;
        double[][] d;
        int[] indexes;
        double[] sortedX;
        double[] sortedY;
        
        double[] d2;
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
     * Runtime is O(n^2) where n is the number of points in x which is
     * the same as the number in y.
     * 
     * NOTE: redundant points are possible in the rankings as "ties" are handled
     * in the algorithm.  This is one advantage over the similar
     * algorithm of Huo and Szekely (2016).
     * 
     * @param x
     * @param y
     * @return 
     */
    public static DCov _bruteForceCovariance(double[] x, double[] y) {
                  
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
            idx[0][i] = i+1;
        }
        
        int[] idx_r = new int[n];
        
        double[] csumT = new double[n+1];
        double[] d = new double[n];
        
        int gap, k, kf, z, j, st1, e1, st2, e2, idx1, idx2;
        int i = 1;
        int r = 1;
        int idx_s = 2;
        while (i < n) {
           gap = 2*i; 
           k = 0;
           
           System.arraycopy(idx[r-1], 0, idx_r, 0, n);
      
           //csumT = cusum(y[idx r]); 
           //csumT = (0, csumT );
           csumT[0] = 0;
           for (z = 0; z < n; ++z) {
               csumT[z + 1] = x[ idx_r[z] -1];
           }
           for (z = 2; z <= n; ++z) {
               csumT[z] += csumT[z-1];
           }
           
           j = 1;
        
           while (j < n) {
               
              st1 = j;
              e1 = Math.min(st1 + i - 1, n);
              st2 = j + i;
              e2 = Math.min(st2 + i - 1, n);
              while (( st1 <= e1 ) && ( st2 <= e2 ) ) {
                 k++;
                 idx1 = idx_r[st1-1];
                 idx2 = idx_r[st2-1];
                 if (x[idx1-1] >= x[idx2-1]) {                 
                    idx[idx_s-1][k-1] = idx1;
                    st1++;
                 } else {
                    idx[idx_s-1][k-1] = idx2;
                    st2++;

                    d[idx2-1] += (csumT[e1 + 1 -1] - csumT[st1 -1]);
                 } // end if-else
              } // end while
              if (st1 <= e1) {
                  kf = k + e1 - st1 + 1;
                  int c = st1;
                  for (z = (k+1); z <= kf; ++z) {
                      idx[idx_s-1][z-1] = idx_r[c-1];
                      c++;
                  }
                  k = kf;
              } else if (st2 <= e2) {
                  kf = k + e2 - st2 + 1;
                  int c = st2;
                  for (z = (k+1); z <= kf; ++z) {
                      idx[idx_s-1][z-1] = idx_r[c-1];
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
        
        System.out.printf("sortedX=%s\n", Arrays.toString(sortedX));
        System.out.printf("sortedY=%s\n", Arrays.toString(sortedY));
        System.out.printf("indexes=%s\n", Arrays.toString(indexes));
        System.out.printf("sortedD=%s\n", Arrays.toString(sortedD));

       DCov dcov = new DCov();
       dcov.d2 = sortedD;
       dcov.indexes = indexes;
       dcov.sortedX = sortedX;
       dcov.sortedY = sortedY;
       
       return dcov;
    }
}
