package algorithms.puzzles;

import junit.framework.TestCase;

public class SudokuDLXShavanKaulTest extends TestCase {


    public void test0() {

        char[][] board;
        char[][] expected;

        board = SudokuBackTrackingTest.getBoard2();
        expected = SudokuBackTrackingTest.getSoln2();

        board = SudokuBackTrackingTest.getBoard3();
        expected = SudokuBackTrackingTest.getSoln3();

        board = SudokuBackTrackingTest.getBoard4();
        expected = SudokuBackTrackingTest.getSoln4();

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

        int size = 3;
        SudokuDLXShavanKaul s = new SudokuDLXShavanKaul(size);
        s.read(board);
        s.solve();

        s.print();
    }

}
