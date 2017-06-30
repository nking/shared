package algorithms.misc;

import gnu.trove.list.TIntList;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;

/**
 * miscellaneous methods for counting the occurrence
 * of numbers
 * 
 * @author nichole
 */
public class Frequency {
    
    /**
     * create sorted frequency arrays of input values.
     * 
     * runtime complexity is O(N*log_2(N))
     * 
     * @param input
     * @param outputValues
     * @param outputCounts
     * @param excludeZeros 
     */
    public void calcFrequency(int[][] input, TIntList outputValues,
        TIntList outputCounts, boolean excludeZeros) {
        
        TIntIntMap valueCountMap = new TIntIntHashMap();
        for (int i = 0; i < input.length; ++i) {
            for (int j = 0; j < input[i].length; ++j) {
       
                int v = input[i][j];
                if (v == 0 && excludeZeros) {
                    continue;
                }
                //NOTE: trove maps return value 0 when key isn'r present
                int c = valueCountMap.get(v);
                valueCountMap.put(v, c + 1);
            }
        }
        outputValues.clear();
        outputValues.add(valueCountMap.keys());
        outputValues.sort();
    
        for (int i = 0;i < outputValues.size(); ++i) {
            int v = outputValues.get(i);
            outputCounts.add(valueCountMap.get(v));
        }
    }
}
