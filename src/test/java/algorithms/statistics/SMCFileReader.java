package algorithms.statistics;

import algorithms.util.ResourceFinder;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class SMCFileReader {

    public static double[] readDiffFile(String fileName) throws IOException {

        String path = ResourceFinder.findFileInTestResources(fileName);
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        TDoubleList d = new TDoubleArrayList();
        try {
            in = new BufferedReader(new FileReader(new File(path)));
            String line = in.readLine();
            while (line != null) {
                d.add(Double.parseDouble(line));
                line = in.readLine();
            }
        } finally {
            if (in != null) {
                in.close();
            }
            if (reader != null) {
                reader.close();
            }
        }
        double[] X = d.toArray();
        return X;
    }
}
