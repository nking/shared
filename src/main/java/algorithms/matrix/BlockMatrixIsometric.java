package algorithms.matrix;

/**
 * a class to hold a matrix that is partitioned into blocks of the same size.
 * @author nichole
 */
public class BlockMatrixIsometric {
    /**
     * the matrix
     */
    private final double[][] a;
    /**
     * the blocks' sizes in the row dimension
     */
    private final int bSize0; 
    
    /**
     * the blocks' sizes in the column dimension
     */
    private final int bSize1;
    
    /**
     * constructor.  The matrix a is not copied in.   This is a 
     * pass-by-reference constructor.
     * @param a
     * @param bSize0 block size along the first dimension (rows)
     * @param bSize1 block size along the second dimension (columns)
     */
    public BlockMatrixIsometric(double[][] a, int bSize0, int bSize1) {
        int n = a.length;
        int m = a[0].length;
        if ((n % bSize0) != 0) {
            throw new IllegalArgumentException("a.length must be evenly divided"
                    + " by bSize0");
        }
        if ((m % bSize1) != 0) {
            throw new IllegalArgumentException("a[0].length must be evenly divided"
                    + " by bSize1");
        }
        this.a = a;
        this.bSize0 = bSize0;
        this.bSize1 = bSize1;
    }
    
    /**
     * set values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * to the given block b.  b must be the same size as the block size
     * of this instance.
     * @param b a given block of values to replace a block of this instance with.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void setBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= getA().length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= getA()[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        for (int i = 0; i < this.bSize0; ++i) {
            System.arraycopy(b[i], 0, getA()[start0 + i], start1, bSize1);
        }
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     * @return contents of (blockNumber0, blockNumber1) of matrix a in this instance.
     */
    public double[][] getBlock(int blockNumber0, int blockNumber1) {
        double[][] b = MatrixUtil.zeros(bSize0, bSize1);
        getBlock(b, blockNumber0, blockNumber1);
        return b;
    }
    
    /**
     * get values in internal matrix a for block number (blockNumber0, blockNumber1) 
     * and set them into the given block b.
     * @param b the output block for the contents of (blockNumber0, blockNumber1) 
     * of matrix a in this instance.
     * @param blockNumber0 the number of the block along the first dimension,
     * that is, rows.
     * @param blockNumber1 the number of the block along the second dimension,
     * that is, columns.
     */
    public void getBlock(double[][] b, int blockNumber0, int blockNumber1) {
        if (b.length != bSize0 || b[0].length != bSize1) {
            throw new IllegalArgumentException("b size must equal the block size"
                    + " within this instance");
        }
        if ((blockNumber0 * bSize0) >= getA().length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber0 is out of bounds");
        }
        if ((blockNumber1 * bSize1) >= getA()[0].length) {
            throw new ArrayIndexOutOfBoundsException("blockNumber1 is out of bounds");
        }
        int start0 = blockNumber0 * bSize0;
        int start1 = blockNumber1 * bSize1;
        for (int i = 0; i < this.bSize0; ++i) {
            System.arraycopy(getA()[start0 + i], start1, b[i], 0, bSize1);
        }
    }

    /**
     * get the entire matrix of this instance.  Note that this is a reference
     * to the instance and is not a copy of it, so modifications of the
     * returned matrix affect this instance too.
     * @return the a
     */
    public double[][] getA() {
        return a;
    }

    /**
     * get the size of dimension 0 of the blocks in this instance
     * @return the bSize0
     */
    public int getBlockSize0() {
        return bSize0;
    }

    /**
     * get the size of dimension 1 of the blocks in this instance
     * @return the bSize1
     */
    public int getBlockSize1() {
        return bSize1;
    }
    
}
