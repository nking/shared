package algorithms.puzzles;

import java.util.*;

/**
 * A backtracking solver that in practice, completes in less time than 9! by recalculating the ordering of indexes.
 * It could be improved.
 *
 * NOTE: the code expects a solvable board.

 @author Nichole King

 * Another approach to this problem is
 * SudokuDLXShavanKaul which uses the Dancing Links datastructure of Knuth's
 and an exact cover of constraints approach..

 Here are some brief notes on what would be faster with a Dancing Links representation and
 exact cover matrix(es) rather than the sets and maps that I'm using below.

 The Knuth repesentation of the exact cover problem reformats a problem into a binary matrix where,
 for any column, only 1 row can be set for that column.

 The Soduko problem constraints would need such a matrix for the number '1',
 and a similar matrix for the number '2', etc.
 where columns are column numbers 0:8 of the board and rows are row numbers 0:8 of the board.

 The 3x3 sections would each need a matrix, but the columns would be numbers 1-9 (or mapped to 0:8)
 and the rows would be the internal indexes of the 3x3 section mapped to values 0:8.

 So the total number of exact set format style matrices would be:
 9 matrices of 9x9, 1 for each number
 9 matrices of 9x9, 1 for each 3x3 section of the board.
 = 18 matrices of 9x9 for a dense representation.

 Dancing Links makes a sparse representation of each of those matrices by using linked lists
 with each node having up, down, left, right nodes where up, down are along columns and left, right
 are along rows.  Each row is a circularly doubly linked list and so is each column.

 those matrices hold candidate values, that is, remaining possible values.

 The matrices, combined into one matrix, can be formatted as described by
 https://www.ocf.berkeley.edu/~jchu/publicportal/sudoku/sudoku.paper.html
 He also gives an example of the initial matrix for a 4x4 suduko before
 modifying it for a given puzzle board:
 https://www.ocf.berkeley.edu/~jchu/publicportal/sudoku/4x4.dlx.64x64.xls

 The column titles are:
 row 1				row 2				row 3				row 4				row 1				row 2				row 3				row 4				col 1				col 2				co3				col 4				region 1				region 2				region 3				region 4
 each having 4 columns in the title.
 The rows are:
 (cand,r,c)
 (1,1,1)
 (2,1,1)
 (3,1,1)
 (4,1,1)
 (1,1,2)
 (2,1,2)
 (3,1,2)
 (4,1,2)
 etc... 64 rows and 64 columns where cand is a possible value.

 The sparse representation, if used here, would make the Section methods
 addToRemaining(), removeFromRemaining(), populateMinRemaining()
 faster.

 The SudukoDLXShavanKaul use of Dancing Links (DLX) uses one matrix for all constraints
 and uses an encoding to store the initial grid information.

 A brief comparison of test results shows that the SudukoDLXShavanKaul code is about
 5 times faster than this SudukoBackTracking code implemented below.
 */
public class SudokuBackTracking {

    protected static int ONE = (int)'1';
    protected boolean solved = false;
    protected char[][] board = null;
    boolean reordRecur = false;
    long nIter = 0;

    protected static class Section {
        int rowSection; // range [0,2]
        int colSection;  // range [0,2]
        //int rowOffset; // offset from 0, used to lookup
        //int colOffset; // offset from 0, used to lookup
        // keys are board row number - rowOffset, similarly for cols.
        // values are the possible numbers present in that column
        Map<Integer, Map<Integer, Set<Integer>>> remainingPossibleValues = new HashMap<>();
        Set<Integer> remaining = new HashSet<Integer>();
        int[][] values = new int[3][];
        public Section(char[][] board, int rowSection, int colSection) {
            if (colSection < 0  || colSection > 2) {
                throw new IllegalArgumentException(
                        "colSection must be >= 0 and <= 2");
            }
            if (rowSection < 0  || rowSection > 2) {
                throw new IllegalArgumentException(
                        "rowSection must be >= 0 and <= 2");
            }
            this.colSection = colSection;
            this.rowSection = rowSection;

            // initialize data structures and fill
            int i, j, k;
            int colOffset = colSection*3;
            int rowOffset = rowSection*3;
            int d, r, c;

            for (i = 0; i < 9; ++i) {
                remaining.add(i);
            }
            for (i = 0; i < 3; ++i) {
                values[i] = new int[3];
                Arrays.fill(values[i], -1);
            }
            for (i = 0; i < 3; ++i) {
                remainingPossibleValues.put(i, new HashMap<Integer, Set<Integer>>());
                for (j = 0; j < 3; ++j) {
                    remainingPossibleValues.get(i).put(j, new HashSet<Integer>());
                }
            }
            for (i = rowOffset; i < rowOffset + 3; ++i) {
                for (j = colOffset; j < colOffset + 3; ++j) {
                    if (board[i][j] == '.') continue;
                    d = (int)board[i][j] - ONE;
                    c = j - colOffset;
                    r = i - rowOffset;
                    values[r][c] = d;
                    remaining.remove(d);
                }
            }
            for (int v : remaining) {
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        remainingPossibleValues.get(i).get(j).add(v);
                    }
                }
            }
            for (i = rowOffset; i < rowOffset + 3; ++i) {
                for (j = colOffset; j < colOffset + 3; ++j) {
                    if (board[i][j] == '.') continue;
                    d = (int)board[i][j] - ONE;
                    c = j - colOffset;
                    r = i - rowOffset;
                    remainingPossibleValues.get(r).get(c).clear();
                }
            }

            // for all rows for the columns in our section, rm remaining
            for (i = 0; i < board.length; ++i) {
                for (j = colOffset; j < colOffset + 3; ++j) {
                    if (board[i][j] == '.') continue;
                    d = (int)board[i][j] - ONE;
                    c = j - colOffset;
                    for (k = 0; k < 3; ++k) {
                        remainingPossibleValues.get(k).get(c).remove(d);
                    }
                }
            }
            // for all cols for the rows in our section, rm remaining
            for (i = rowOffset; i < rowOffset + 3; ++i) {
                r = i - rowOffset;
                for (j = 0; j < board[i].length; ++j) {
                    if (board[i][j] == '.') continue;
                    d = (int)board[i][j] - ONE;
                    for (k = 0; k < 3; ++k) {
                        remainingPossibleValues.get(r).get(k).remove(d);
                    }
                }
            }
        }

        public boolean[] removeFromRemaining(int row, int col, int value) {
            int c = col - colSection*3;
            int r = row - rowSection*3;//
            assert(c < 3 && c > -1);
            assert(r < 3 && r > -1);
            boolean[] b = new boolean[1+3+3];
            values[r][c] = value;

            b[0] = remaining.remove(value);

            int i;
            // remove from all [i][c]
            for (i = 0; i < 3; ++i) {
                b[1 + i] = remainingPossibleValues.get(i).get(c).remove(value);
            }
            // remove from all [r][i]
            for (i = 0; i < 3; ++i) {
                b[4 + i] = remainingPossibleValues.get(r).get(i).remove(value);
            }
            return b;
        }
        public void addToRemaining(int row, int col, int value, boolean[] b) {
            int c = col - colSection*3;
            int r = row - rowSection*3;
            assert(c < 3 && c > -1);
            assert(r < 3 && r > -1);
            values[r][c] = -1;
            if (b[0])  remaining.add(value);
            // remove from all [i][c]
            int i;
            for (i = 0; i < 3; ++i) {
                if (b[1+i]) remainingPossibleValues.get(i).get(c).add(value);
            }
            // remove from all [r][i]
            for (i = 0; i < 3; ++i) {
                if (b[4+i]) remainingPossibleValues.get(r).get(i).add(value);
            }
        }
        protected static String print(Collection<Integer> a) {
            StringBuilder sb = new StringBuilder();
            for (int x : a) {
                if (sb.length() > 0) sb.append(",");
                sb.append(String.format("%3d", x));
            }
            return sb.toString();
        }

        /**
         *
         * @param row row number in 9x9 reference frame
         * @param col col number in 9x9 reference frame
         * @param sections
         * @return
         */
        public static int[] intersectionOfRemaining(int row, int col,
                                                    Section[][] sections) {
            int rowSect = row/3;
            int colSect = col/3;
            int rowInRowSect = row - rowSect*3;
            int colInColSect = col - colSect*3;
            Set<Integer> remaining = sections[rowSect][colSect].remainingPossibleValues.get(rowInRowSect).get(colInColSect);
            // this can happen when an assignment of same value occurs in same section,
            // but neither the row nor col were same as these.
            // so remove and from remainingPossibleValues not in remaining.
            //assert(sections[rowSect][colSect].remaining.size() >= remaining.size());
            int[] out = new int[remaining.size()];
            int i = 0;
            for (int v : remaining) {
                if (!sections[rowSect][colSect].remaining.contains(v)) continue;
                out[i] = v;
                ++i;
            }
            if (i < out.length) {
                out = Arrays.copyOf(out, i);
            }
            return out;
        }

        /**
         @param out map is populated with the method.  key=sectIdx,
         value = min number of remaining values for the cell
         calculated from min num values from row and col.
         */
        public void populateMinRemaining(Map<Integer, Integer> out) {
            int i, j, sectIdx, sectRow, sectCol;
            for (i = 0; i < 3; ++i) { // row
                for (j = 0; j < 3; ++j) { //col
                    if (values[i][j] > -1)  continue;
                    sectRow = rowSection*3 + i;
                    sectCol = colSection*3 + j;
                    sectIdx = (sectRow * 9) + sectCol;
                    assert(sectIdx < 81);
                    Set<Integer> remaining2 = remainingPossibleValues.get(i).get(j);
                    int n = 0;
                    for (int nr : remaining2) {
                        if (remaining.contains(nr)) ++n;
                    }
                    //assert(remaining.size() >= remaining2.size());
                    out.put(sectIdx, n);
                }
            }
        }
    }
    public static class Sections {
        boolean valid = true;
        Section[][] sections = new Section[3][];
        public Sections(char[][] board) {
            int i, j;
            for (i = 0; i < 3; ++i) {
                sections[i] = new Section[3];
                for (j = 0; j < 3; ++j) {
                    sections[i][j] = new Section(board, i, j);
                }
            }
        }

        public boolean[] removeSameRowOtherCols(int row, int col, int v) {
            int rowSect = row/3;
            int colSect = col/3;
            int i,j;
            boolean[] b = new boolean[6];
            int idx = 0;
            int r = row - rowSect*3;
            for (i = 0; i < 3; ++i) {// colSect
                if (i == colSect) continue;
                for (j = 0; j < 3; ++j) {
                    b[idx] = sections[rowSect][i].remainingPossibleValues.get(r).get(j).remove(v);
                    ++idx;
                }
            }
            return b;
        }
        public void addSameRowOtherCols(int row, int col, int v, boolean[] b) {
            int rowSect = row/3;
            int colSect = col/3;
            int i,j;
            int idx = 0;
            int r = row - rowSect*3;
            for (i = 0; i < 3; ++i) {// colSect
                if (i == colSect) continue;
                for (j = 0; j < 3; ++j) {
                    if (b[idx]) sections[rowSect][i].remainingPossibleValues.get(r).get(j).add(v);
                    ++idx;
                }
            }
        }
        public boolean[] removeSameColOtherRows(int row, int col, int v) {
            int rowSect = row/3;
            int colSect = col/3;
            int i,j;
            boolean[] b = new boolean[6];
            int idx = 0;
            int c = col - colSect*3;
            for (i = 0; i < 3; ++i) {//rowSect
                if (i == rowSect) continue;
                for (j = 0; j < 3; ++j) {
                    b[idx] = sections[i][colSect].remainingPossibleValues.get(j).get(c).remove(v);
                    ++idx;
                }
            }
            return b;
        }
        public void addSameColOtherRows(int row, int col, int v, boolean[] b) {
            int rowSect = row/3;
            int colSect = col/3;
            int i,j;
            int idx = 0;
            int c = col - colSect*3;
            for (i = 0; i < 3; ++i) {//rowSect
                if (i == rowSect) continue;
                for (j = 0; j < 3; ++j) {
                    if (b[idx]) sections[i][colSect].remainingPossibleValues.get(j).get(c).add(v);
                    ++idx;
                }
            }
        }
        protected void populateValues(Set<Integer> set) {
            set.clear();
            for (int i = 0; i < 9; ++i) set.add(i);
        }
        public boolean valuesAreValidSolution() {
            int i, j, k;
            Set<Integer> expected = new HashSet<>();

            // for each column, check that all rows have 0-8 disjointly
            for (int col = 0; col < 9; ++col) {
                int colSect = col / 3;
                int colInSect = col - 3*colSect;
                populateValues(expected);
                for (int row = 0; row < 9; ++row) {
                    int rowSect = row / 3;
                    int rowInSect = row - 3*rowSect;
                    int v = sections[rowSect][colSect].values[rowInSect][colInSect];
                    boolean bv = expected.remove(v);
                    if (!bv) return false;
                }
            }

            // for each row, check that all cols have 0-8 disjointly
            for (int row = 0; row < 9; ++row) {
                int rowSect = row / 3;
                int rowInSect = row - 3*rowSect;
                populateValues(expected);
                for (int col = 0; col < 9; ++col) {
                    int colSect = col / 3;
                    int colInSect = col - 3*colSect;
                    int v = sections[rowSect][colSect].values[rowInSect][colInSect];
                    boolean bv = expected.remove(v);
                    if (!bv) return false;
                }
            }

            // for each section, check that values 0-8 are present
            for (int rowSect = 0; rowSect < 3; ++rowSect) {
                for (int colSect = 0; colSect < 3; ++colSect) {
                    populateValues(expected);
                    for (int r = 0; r < 3; ++r) {
                        for (int c = 0; c < 3; ++c) {
                            int v = sections[rowSect][colSect].values[r][c];
                            boolean bv = expected.remove(v);
                            if (!bv) return false;
                        }
                    }
                }
            }
            return true;
        }
    }

    /**
     * given a solvable board, modifies the board to hold the solution.
     * @param board input incomplete but solvable board.
     * empty cells should have a '.'.
     *              the numbers 1 through 9 are represented as their characters.
     */
    public void solveSudoku(char[][] board) {
        solved = false;
        this.board = board;

        nIter = 0;

        Sections sections = new Sections(board);

        int[] indexes = calcOrder(sections.sections)[0];

        recursion(0, indexes, sections);

        assert(solved);

        System.out.printf("nIter=%d\n", nIter);
    }

    protected void recursion(int idx, int[] indexes, Sections sections) {
        ++nIter;
        if (solved) return;
        if (!sections.valid) return;
        if (idx == indexes.length) {
            //System.out.printf("found solution at nIter=%d.  sections is valid=%b\n", nIter,
            //        sections.valuesAreValidSolution());
            solved = true;

            int row, col;
            Section sect;
            for (int rowSect = 0; rowSect < 3; ++rowSect) {
                for (int colSect = 0; colSect < 3; ++colSect) {
                    sect = sections.sections[rowSect][colSect];
                    for (int rowInRowSect = 0; rowInRowSect < 3; ++rowInRowSect) {
                        row = rowInRowSect + rowSect*3;
                        for (int colInRowSect = 0; colInRowSect < 3; ++colInRowSect) {
                            col = colInRowSect + colSect*3;
                            board[row][col] = (char)(ONE + sect.values[rowInRowSect][colInRowSect]);
                        }
                    }
                }
            }
            return;
        }

        int sectIdx = indexes[idx];
        int row = calcRowIdx(sectIdx);
        int col = calcColIdx(sectIdx);
        int rowSect = row/3;
        int colSect = col/3;

        int[] possible = Section.intersectionOfRemaining(row, col, sections.sections);
        if (possible.length == 0) {
            sections.valid=false;
            return;
        } else if (possible.length > 1) {
            // recalc order to see if can improve/reduce number of posssibilities
            int[][] indexesAndN = calcOrder(sections.sections);
            if (indexesAndN[1][0] == 1) {
                int[] indexes2 = new int[indexes.length];
                // copy indexes[0:idx-1] into indexes3
                System.arraycopy(indexes, 0, indexes2, 0, idx);
                System.arraycopy(indexesAndN[0], 0, indexes2, idx, indexesAndN[0].length);
                reordRecur = true;
                recursion(idx, indexes2, sections);
                return;
            }
        }

        int i, v;
        for (i = 0; i < possible.length; ++i) {
            if (solved) return;
            v = possible[i];

            boolean[] b = sections.sections[rowSect][colSect].removeFromRemaining(row, col, v);

            boolean[] b1 = sections.removeSameRowOtherCols(row, col, v);

            boolean[] b2 = sections.removeSameColOtherRows(row, col, v);

            sections.valid = true;

            recursion(idx+1, indexes, sections);

            // undo changes just made
            sections.sections[rowSect][colSect].addToRemaining(row, col, v, b);

            sections.addSameRowOtherCols(row, col, v, b1);
            sections.addSameColOtherRows(row, col, v, b2);
        }
    }

    protected int[][] calcOrder(Section[][] sections) {
        int i, j;
        // key = sectIdx, value = number of possible values
        Map<Integer, Integer> sectIdxNR = new HashMap<>();
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                sections[i][j].populateMinRemaining(sectIdxNR);
            }
        }

        int[] indexes = new int[sectIdxNR.size()];
        int[] nVals = new int[sectIdxNR.size()];

        i = 0;
        for (int sectIdx : sectIdxNR.keySet()) {
            indexes[i] = sectIdx;
            nVals[i] = sectIdxNR.get(sectIdx);
            ++i;
        }

        //counting sort, O(N) where N = sectIdxNR.size(),
        // that is the number of unfilled cells
        int[] sortedIdxs = sort(nVals);
        int[][] out = new int[2][sortedIdxs.length];
        for (i = 0; i < sortedIdxs.length; ++i) {
            out[0][i] = indexes[sortedIdxs[i]];
            out[1][i] = nVals[i];
        }

        return out;
    }

    protected int calcColIdx(int sectIdx) {
        return sectIdx % 9;
    }
    protected int calcRowIdx(int sectIdx) {
        return sectIdx / 9;
    }
    protected int calcSectIdx(int row, int col) {
        return (row * 9) + col;
    }

    // counting sort, O(n)
    protected int[] sort(int[] a) {
        int[] c = new int[a.length + 1];
        int i;
        for (i = 0; i < a.length; ++i) {
            c[a[i]]++;
        }
        for (i = 1; i < c.length; ++i) {
            c[i] += c[i-1];
        }
        int[] b = new int[a.length];
        int[] indexes = new int[a.length];
        int aI;
        for (i = a.length - 1; i > -1; --i) {
            aI = a[i];
            c[aI]--;
            b[c[aI]] = aI;
            indexes[c[aI]] = i;
        }
        System.arraycopy(b, 0, a, 0, a.length);
        return indexes;
    }

}
