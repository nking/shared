package algorithms.optimization;

import algorithms.util.ResourceFinder;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import static junit.framework.Assert.assertEquals;

/**
 *
 * @author nichole
 */
public class ORSetCoverReader {
    
    public static TIntSet getSCPDatasetCoverSet(String filename) {
        if (filename.startsWith("scpe1")) {
            // all sets have weight of 1
            return null;
        } else if (filename.startsWith("scp41")) {
            // set weights are not all 1
            return null;
        } else {
            throw new IllegalArgumentException("filename " + filename 
            + " isn't one of the choices.");
        }
    }
    
    public static double getSCPDatasetZ(String filename) {
        if (filename.startsWith("scpe1")) {
            // all sets have weight of 1
            return 5;
        } else if (filename.startsWith("scp41")) {
            // set weights are not all 1
            return 429;
        } else {
            throw new IllegalArgumentException("filename " + filename 
            + " isn't one of the choices.");
        }
    }
    
    /**
     * get [nSets, nItems]
     * @param filename
     * @return 
     */
    public static int[] getSCPDatasetNumberOfRowsCols(String filename) {
        if (filename.startsWith("scpe1")) {
            // all sets have weight of 1
            return new int[]{50, 500};
        } else if (filename.startsWith("scp41")) {
            // set weights are not all 1
            return new int[]{200, 1000};
        } else {
            throw new IllegalArgumentException("filename " + filename 
            + " isn't one of the choices.");
        }
    }
    public static void getSCPDataset(String filename, List<TIntSet> sets, 
        double[] weights) throws IOException {
    
        String filePath = ResourceFinder.findFileInTestResources(filename);
            
        readFile(filePath, sets, weights);
                
    }
    
    private static int readFile(String filePath, List<TIntSet> sets, double[] weights) 
        throws IOException {
        
        FileReader reader = null;
        BufferedReader in = null;
        
        int nSets, nItems;

        int col, row, i, ii, nSet;
        TIntSet set;
        
        String space = "\\s+";
        
        /*
        mItems nSets
        weights(nSets of them over many lines)
        number_of_items_in_set_0
        items_in_set_0
        ...
        number_of_items_in_set_nSet-1
        items_in_set_nSet-1
        */
        try {
            in = new BufferedReader(new FileReader(new File(filePath)));
            
            String line = in.readLine().trim();
            String[] lineElements = line.split(space);
            nSets = Integer.valueOf(lineElements[0]); //200
            nItems = Integer.valueOf(lineElements[1]); //1000
            assertEquals(nItems, weights.length);
            
            // read each weight
            for (i = 0; i < nItems;) {
                // read the set weights
                line = in.readLine().trim();
                lineElements = line.split(space);
                for (String element : lineElements) {
                    weights[i] = Double.valueOf(element);
                    i++;
                }
            }
            
            int nU;
            
            // read the number of items in each set followed by the items in each set
            line = in.readLine();
            while (line != null) {
                
                line = line.trim();
                if (line.length()==0) {
                    break;
                }
                
                nU = Integer.valueOf(line);
                set = new TIntHashSet();
                i = 0;
                while (i < nU) {
                    line = in.readLine().trim();
                    lineElements = line.split(space);
                    for (String element : lineElements) {
                        ii = Integer.valueOf(element) - 1;// subtr 1 for zero-based indexes
                        set.add(ii);
                        i++;
                    }
                }
                sets.add(set);
                
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
        assertEquals(nSets, sets.size());
        return nItems;
    }

}
