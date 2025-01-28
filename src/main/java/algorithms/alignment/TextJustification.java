package algorithms.alignment;

import java.util.*;

public class TextJustification {

    protected int counter = 0;


    /**

     * @param output an empty output array to fill with justified text
     * @return total number of added spaces
     */
    public int justify(String[] words, int width, List<String> output) {

        counter = 0;

        //key = wordIdx, value=map with key=position in row, value=penalty
        Map<Integer, BestSoln> memo = new HashMap<>();
        List<List<Integer>> rows = new ArrayList<>();

	// initialize with first word
        rows.add(new ArrayList<>());
        rows.get(0).add(0);
        int currRSum = words[0].length();

	//(row index, word index, int prevPSum, int currRSum,...
        BestSoln best = recursion(0, 1, 0, currRSum, width, words, rows, memo);

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

    protected BestSoln recursion(int rIdx, int wIdx, int prevPSum, int currRSum, final int width,
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
        int cost0 = w.length() + 1; // the 1 is counting space before the word
        if (wIdx == 0) --cost0;

        // add to current row if can
        BestSoln p0 = new BestSoln(null, Integer.MAX_VALUE);
        if ((width - (currRSum + cost0)) >= 0) {
            if (rows.size() <= rIdx) {
                rows.add(new ArrayList<>());
            }
            rows.get(rIdx).add(wIdx);
            p0 = recursion(rIdx, wIdx + 1, prevPSum, currRSum + cost0, width, words, rows, memo);
            //backtracking: remove wIdx so can use it on p1 recursion:
            rows.get(rIdx).remove(rows.get(rIdx).size() - 1);
        }

        // posIdx = 0 for new row
        if (rows.size() <= rIdx + 1) {
            rows.add(new ArrayList<>());
        }
        rows.get(rIdx + 1).add(wIdx);
        int prevPSum1 = prevPSum + (width - currRSum);
        int cost1 = cost0 - 1;
        BestSoln p1 = recursion(rIdx + 1, wIdx + 1, prevPSum1, cost1, width, words, rows, memo);
	//backtracking: remove so can try another permutation
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

    public static class Best {
        int added;
        List<String> lines;
        public Best(int added, List<StringBuilder> lines) {
            this.added = added;
            this.lines = new ArrayList<>();
            for (StringBuilder sb : lines) {
                this.lines.add(sb.toString());
            }
        }
    }

    /**
     * given a string of words and a column width, minimize the amount of space added between words where the
     * number of spaces between words is at least 1.
     * returns the total number of added spaces at ends of lines, and the words partitioned into lines.
     *
     * This version is recursive, like justify(), but is more readable.
     *
     * @param words list of words to consecutively justify to line lengths LEQ width
     * @param width width of column to place consecutive words and spaces between them
     * @return the total number of added spaces at ends of lines, and the words partitioned into lines.
     */
    public static Best justify3(String[] words, int width) {
        Map<Integer, Best> memo = new HashMap<>();
        List<StringBuilder> lines = new ArrayList<>();
        StringBuilder line = new StringBuilder(words[0]);
        lines.add(line);
        int endSpace = width - words[0].length();
        return r3(1, 0,
                words, width,
                memo, lines, line);
    }

    // 1 space is added before each word, except when it is first on a line.
    // line is already added to lines.
    private static Best r3(int wIdx, int added,
                          String[] words, int width,
                          Map<Integer, Best> memo, List<StringBuilder> lines, StringBuilder line) {
        if (wIdx == words.length) {
            // add end of current line
            added += (width - line.length());
            return new Best(added, lines);
        }
        if (memo.containsKey(wIdx)) {
            return memo.get(wIdx);
        }

        // we add a space before the word
        int wLen = 1 + words[wIdx].length();

        Best solnIncl = null;

        // include if can
        if (line.length() + wLen <= width) {
            line.append(" ").append(words[wIdx]);
            solnIncl = r3(wIdx+1, added, words, width,
                    memo, lines, line);

            //back track, undo update to line
            line.delete(line.length() - wLen, line.length());
        }

        // exclude by putting it on the next line
        int diff = width - line.length();

        StringBuilder line2 = new StringBuilder().append(words[wIdx]);
        lines.add(line2);
        Best solnExcl = r3(wIdx+1, added + diff,
                words, width, memo,
                lines, line2);

        // back track, remove line
        lines.remove(lines.size() - 1);

        if (solnIncl == null) {
            memo.put(wIdx, solnExcl);
        } else {
            // consider how to handle ties if want to enumerate
            if (solnIncl.added <= solnExcl.added) {
                memo.put(wIdx, solnIncl);
            } else {
                memo.put(wIdx, solnExcl);
            }
        }

        return memo.get(wIdx);
    }
}
