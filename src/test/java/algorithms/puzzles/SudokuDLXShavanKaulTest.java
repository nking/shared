package algorithms.puzzles;

import junit.framework.TestCase;

public class SudokuDLXShavanKaulTest extends TestCase {

    public void test0() {

        boolean doPrint = false;
        char[][] board;
        char[][] expected;
        for (int ii = 4; ii < 5; ++ii) {

            switch (ii) {
                case 0:
                    board = SudokuBackTrackingTest.getBoard1();
                    expected = SudokuBackTrackingTest.getSoln1();
                    break;
                case 1:
                    board = SudokuBackTrackingTest.getBoard2();
                    expected = SudokuBackTrackingTest.getSoln2();
                    break;
                case 2:
                    board = SudokuBackTrackingTest.getBoard3();
                    expected = SudokuBackTrackingTest.getSoln3();
                    break;
                case 3:
                    board = SudokuBackTrackingTest.getBoard4();
                    expected = SudokuBackTrackingTest.getSoln4();
                    break;
                case 4:
                    board = SudokuBackTrackingTest.getBoard5();
                    expected = SudokuBackTrackingTest.getSoln5();
                    break;
                default:
                    throw new IllegalArgumentException("i not recognized");
            }

            if (doPrint) {
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
            }

            int size = 3;
            SudokuDLXShavanKaul s = new SudokuDLXShavanKaul(size);
            s.read(board);
            s.solve();
            //s.print();;

            board = s.getSolution();

            if (doPrint) {
                System.out.printf("result=\n");
                for (int i = 0; i < board.length; ++i) {
                    for (int j = 0; j < board[i].length; ++j) {
                        System.out.printf("%3s", board[i][j]);
                    }
                    System.out.printf("\n");
                }
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

    }

}
