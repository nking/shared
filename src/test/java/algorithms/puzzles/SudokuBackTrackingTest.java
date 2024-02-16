package algorithms.puzzles;

import junit.framework.TestCase;

import java.util.*;

public class SudokuBackTrackingTest extends TestCase {

    public void test0() throws InterruptedException {

        char[][] board;
        char[][] expected;

        //board = getBoard1();
        //expected = getSoln1();
        //board = getBoard3();
        //expected = getSoln3();

        //board = getBoard2();
        //expected = getSoln2();

        board = getBoard4();
        expected = getSoln4();

        System.out.printf("board=\n");
        for (int i = 0; i < board.length; ++i) {
            for (int j = 0; j < board[i].length; ++j) {
                System.out.printf("%3s", board[i][j]);
            }
            System.out.printf("\n");
        }
        System.out.printf("expected=\n");
        for (int i = 0; i < board.length; ++i) {
            for (int j = 0; j < board[i].length; ++j) {
                System.out.printf("%3s", expected[i][j]);
            }
            System.out.printf("\n");
        }

        SudokuBackTracking sudoku = new SudokuBackTracking();
        sudoku.solveSudoku(board);

        System.out.printf("result=\n");
        for (int i = 0; i < board.length; ++i) {
            for (int j = 0; j < board[i].length; ++j) {
                System.out.printf("%3s", board[i][j]);
            }
            System.out.printf("\n");
        }
        for (int i = 0; i < board.length; ++i) {
            for (int j = 0; j < board[i].length; ++j) {
                if (expected[i][j] != board[i][j]) {
                    System.out.printf("Error at [%d][%d]\n", i, j);
                }
                assertTrue(expected[i][j] == board[i][j]);
            }
        }
    }

    protected char[][] getBoard1() {
        char[][] b = new char[9][];
        b[0]= new char[]{'5','3','.','.','7','.','.','.','.'};
        b[1]= new char[]{'6','.','.','1','9','5','.','.','.'};
        b[2]= new char[]{'.','9','8','.','.','.','.','6','.'};
        b[3]= new char[]{'8','.','.','.','6','.','.','.','3'};
        b[4]= new char[]{'4','.','.','8','.','3','.','.','1'};
        b[5]= new char[]{'7','.','.','.','2','.','.','.','6'};
        b[6]= new char[]{'.','6','.','.','.','.','2','8','.'};
        b[7]= new char[]{'.','.','.','4','1','9','.','.','5'};
        b[8]= new char[]{'.','.','.','.','8','.','.','7','9'};
        return b;
    }
    protected char[][] getSoln1() {
        char[][] b = new char[9][];
        b[0]=new char[]{'5','3','4','6','7','8','9','1','2'};
        b[1]=new char[]{'6','7','2','1','9','5','3','4','8'};
        b[2]=new char[]{'1','9','8','3','4','2','5','6','7'};
        b[3]=new char[]{'8','5','9','7','6','1','4','2','3'};
        b[4]=new char[]{'4','2','6','8','5','3','7','9','1'};
        b[5]=new char[]{'7','1','3','9','2','4','8','5','6'};
        b[6]=new char[]{'9','6','1','5','3','7','2','8','4'};
        b[7]=new char[]{'2','8','7','4','1','9','6','3','5'};
        b[8]=new char[]{'3','4','5','2','8','6','1','7','9'};
        return b;
    }

    public static char[][] getBoard2() {
        char[][] b = new char[9][];
        b[0]=new char[]{'.','.','9','7','4','8','.','.','.'};
        b[1]=new char[]{'7','.','.','.','.','.','.','.','.'};
        b[2]=new char[]{'.','2','.','1','.','9','.','.','.'};
        b[3]=new char[]{'.','.','7','.','.','.','2','4','.'};
        b[4]=new char[]{'.','6','4','.','1','.','5','9','.'};
        b[5]=new char[]{'.','9','8','.','.','.','3','.','.'};
        b[6]=new char[]{'.','.','.','8','.','3','.','2','.'};
        b[7]=new char[]{'.','.','.','.','.','.','.','.','6'};
        b[8]=new char[]{'.','.','.','2','7','5','9','.','.'};
        return b;
    }
    public static char[][] getSoln2() {
        char[][] b = new char[9][];
        b[0]=new char[]{'5','1','9','7','4','8','6','3','2'};
        b[1]=new char[]{'7','8','3','6','5','2','4','1','9'};
        b[2]=new char[]{'4','2','6','1','3','9','8','7','5'};
        b[3]=new char[]{'3','5','7','9','8','6','2','4','1'};
        b[4]=new char[]{'2','6','4','3','1','7','5','9','8'};
        b[5]=new char[]{'1','9','8','5','2','4','3','6','7'};
        b[6]=new char[]{'9','7','5','8','6','3','1','2','4'};
        b[7]=new char[]{'8','3','2','4','9','1','7','5','6'};
        b[8]=new char[]{'6','4','1','2','7','5','9','8','3'};
        return b;
    }

    public static char[][] getBoard3() {
        char[][] b = new char[9][];
        b[0]=new char[]{'4','.','.','.','9','6','.','.','8'};
        b[1]=new char[]{'.','5','9','.','2','4','.','.','6'};
        b[2]=new char[]{'.','6','.','3','.','.','.','9','4'};
        b[3]=new char[]{'.','.','2','.','.','.','.','6','.'};
        b[4]=new char[]{'6','8','.','.','.','.','4','5','1'};
        b[5]=new char[]{'.','7','.','.','.','.','.','8','.'};
        b[6]=new char[]{'8','1','5','4','.','.','6','2','7'};
        b[7]=new char[]{'7','.','.','.','.','.','8','.','.'};
        b[8]=new char[]{'2','.','.','.','6','8','.','1','5'};
        return b;
    }
    public static char[][] getSoln3() {
        char[][] b = new char[9][];
        b[0]=new char[]{'4','2','7','1','9','6','5','3','8'};
        b[1]=new char[]{'3','5','9','8','2','4','1','7','6'};
        b[2]=new char[]{'1','6','8','3','5','7','2','9','4'};
        b[3]=new char[]{'9','4','2','5','8','1','7','6','3'};
        b[4]=new char[]{'6','8','3','9','7','2','4','5','1'};
        b[5]=new char[]{'5','7','1','6','4','3','9','8','2'};
        b[6]=new char[]{'8','1','5','4','3','9','6','2','7'};
        b[7]=new char[]{'7','3','6','2','1','5','8','4','9'};
        b[8]=new char[]{'2','9','4','7','6','8','3','1','5'};
        return b;
    }

    public static char[][] getBoard4() {
        char[][] b = new char[9][];
        b[0]=new char[]{'1','2','.','4','.','.','3','.','.'};
        b[1]=new char[]{'3','.','.','.','1','.','.','5','.'};
        b[2]=new char[]{'.','.','6','.','.','.','1','.','.'};
        b[3]=new char[]{'7','.','.','.','9','.','.','.','.'};
        b[4]=new char[]{'.','4','.','6','.','3','.','.','.'};
        b[5]=new char[]{'.','.','3','.','.','2','.','.','.'};
        b[6]=new char[]{'5','.','.','.','8','.','7','.','.'};
        b[7]=new char[]{'.','.','7','.','.','.','.','.','5'};
        b[8]=new char[]{'.','.','.','.','.','.','.','9','8'};
        return b;
    }
    public static char[][] getSoln4() {
        char[][] b = new char[9][9];
        b[0]=new char[]{'1','2','8','4','6','5','3','7','9'};
        b[1]=new char[]{'3','7','4','2','1','9','8','5','6'};
        b[2]=new char[]{'9','5','6','8','3','7','1','4','2'};
        b[3]=new char[]{'7','6','5','1','9','8','4','2','3'};
        b[4]=new char[]{'2','4','9','6','7','3','5','8','1'};
        b[5]=new char[]{'8','1','3','5','4','2','9','6','7'};
        b[6]=new char[]{'5','9','2','3','8','6','7','1','4'};
        b[7]=new char[]{'4','8','7','9','2','1','6','3','5'};
        b[8]=new char[]{'6','3','1','7','5','4','2','9','8'};
        return b;
    }

}
