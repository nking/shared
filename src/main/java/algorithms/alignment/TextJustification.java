package algorithms.alignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TextJustification {

    protected int counter = 0;


    /**
     * given a string of words and a column width, minimize the amount of space added between words where the
     * number of spaces between words is at least 1.
     * returns the total number of added spaces and fill output list with the text justified to line lengths <= width.
     *
     * @param words list of words to consecutively justify to line lengths <= width
     * @param width widht of column to place consecutive words and spaces between them
     * @param output and output array to fill with justified text
     * @return
     */
    public int justify(String[] words, int width, List<String> output) {

        counter = 0;

        //key = wordIdx, value=map with key=position in row, value=penalty
        Map<Integer, BestSoln> memo = new HashMap<>();
        List<List<Integer>> rows = new ArrayList<>();

        rows.add(new ArrayList<>());
        rows.get(0).add(0);
        int currRSum = words[0].length();

        BestSoln best = recurse(0, 1, 0, currRSum, width, words, rows, memo);

        output.clear();
        // transform rows into words
        for (int[] row : best.rows) {
            StringBuilder sb = new StringBuilder();
            for (int i : row) {
                if (sb.length() > 0) {
                    sb.append(" ");
                }
                sb.append(words[i]);
            }
            output.add(sb.toString());
        }

        System.out.printf("nRecursions=%d, words.length=%d\n", counter, words.length);

        return best.cost;
    }

    protected BestSoln recurse(int rIdx, int wIdx, int prevPSum, int currRSum, final int width,
                               String[] words, List<List<Integer>> rows, Map<Integer, BestSoln> memo) {
        ++counter;

        if (wIdx >= words.length) {
            return new BestSoln(rows, prevPSum + width - currRSum);
        }

        if (memo.containsKey(wIdx)) {
            return memo.get(wIdx);
        }

        assert(currRSum > 0);
        String w = words[wIdx];
        int cost0 = w.length() + 1;

        // add to current row if can
        BestSoln p0 = new BestSoln(null, Integer.MAX_VALUE);
        if ((width - (currRSum + cost0)) >= 0) {
            if (rows.size() <= rIdx) {
                rows.add(new ArrayList<>());
            }
            rows.get(rIdx).add(wIdx);
            p0 = recurse(rIdx, wIdx + 1, prevPSum, currRSum + cost0, width, words, rows, memo);
            // remove wIdx so can use it on p1 recursion:
            rows.get(rIdx).remove(rows.get(rIdx).size() - 1);
        }

        // posIdx = 0 for new row
        if (rows.size() <= rIdx + 1) {
            rows.add(new ArrayList<>());
        }
        rows.get(rIdx + 1).add(wIdx);
        int prevPSum1 = prevPSum + (width - currRSum);
        int cost1 = cost0 - 1;
        BestSoln p1 = recurse(rIdx + 1, wIdx + 1, prevPSum1, cost1, width, words, rows, memo);
        rows.get(rIdx + 1).remove(rows.get(rIdx + 1).size() - 1);

        if (p0.cost <= p1.cost) {
            memo.put(wIdx, p0);
            return p0;
        } else {
            memo.put(wIdx, p1);
            return p1;
        }
    }

    protected class BestSoln {
        int[][] rows;
        final int cost;
        public BestSoln(List<List<Integer>> rows, int penalty) {
            this.cost = penalty;
            if (rows != null) {
                this.rows = new int[rows.size()][];
                for (int i = 0; i < rows.size(); ++i) {
                    List<Integer> row = rows.get(i);
                    this.rows[i] = new int[row.size()];
                    for (int j = 0; j < row.size(); ++j) {
                        this.rows[i][j] = row.get(j);
                    }
                }
            } else {
                this.rows = null;
            }
        }
    }
}
