package algorithms.shortestPaths;

import algorithms.util.PairInt;
import algorithms.util.PairIntWithIndex;
import algorithms.util.ResourceFinder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class StanfordMilesReader {

    public final int nCities = 128;
    public final String[] names = new String[nCities];
    // coordinates of the cities as integers
    public final PairIntWithIndex[] coords = new PairIntWithIndex[nCities];
    // map with key=index of city1, value = map w/ key=index of city2, value = distance between cities.  undirectod graph
    public final Map<Integer, Map<Integer, Integer>> distMap = new HashMap<>();

    public void loadFile() throws IOException {

        String path = ResourceFinder.findFileInTestResources("miles.dat");
        FileReader reader = null;
        BufferedReader in = null;
        int count = 0;
        try {
            in = new BufferedReader(new FileReader(new File(path)));

            String line = null;

            String pattern = "^([a-zA-Z\\s-]+,\\s[a-zA-Z]{2})\\[(\\d+),(\\d+)\\](\\d+)$";
            Pattern p = Pattern.compile(pattern);

            Matcher m = null;
            line = in.readLine();
            while (line != null) {
                if (line.startsWith("*")) {
                    // is a comment
                    line = in.readLine();
                    continue;
                }
                //Youngstown, OH[4110,8065]115436
                // group 0 is the whole line
                // group 1 is the city, state
                // group 2 is latitude times 100
                // group 3 is longitude times 100
                // group 4 is population
                m = p.matcher(line);
                if (m.matches()) {
                    String name = m.group(1);
                    names[count] = name;
                    coords[count] = new PairIntWithIndex(Integer.valueOf(m.group(2)), Integer.valueOf(m.group(3)), count);
                    count++;
                } else {
                    int count2 = count - 2;
                    assert(count2 >= 0);
                    assert(count - 1 >= 0);
                    Map<Integer, Integer> d = new HashMap<Integer, Integer>();
                    distMap.put(count - 1, d);
                    String[] edges = line.split("\\s");
                    for (int i = 0; i < edges.length; ++i) {
                        assert(count2 >= 0);
                        d.put(count2, Integer.valueOf(edges[i]));
                        count2--;
                    }
                }
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

        // make the reverse entries
        Map<Integer, Map<Integer, Integer>> revMap = new HashMap<>();
        Iterator<Map.Entry<Integer, Map<Integer, Integer>>> distMapIter = distMap.entrySet().iterator();
        while (distMapIter.hasNext()) {
            Map.Entry<Integer, Map<Integer, Integer>> entry = distMapIter.next();
            Integer key1 = entry.getKey();
            Map<Integer, Integer> map1 = entry.getValue();
            Iterator<Map.Entry<Integer, Integer>> iter2 = map1.entrySet().iterator();
            while (iter2.hasNext()) {
                Map.Entry<Integer, Integer> entry2 = iter2.next();
                Integer key2 = entry2.getKey();
                Integer dist = entry2.getValue();

                // make entry for key2, map2  and add key1 to map2
                Map<Integer, Integer> map2 = revMap.get(key2);
                if (map2 == null) {
                    map2 = new HashMap<Integer, Integer>();
                    revMap.put(key2, map2);
                }
                map2.put(key1, dist);
            }
        }
        distMap.putAll(revMap);
    }

    public void fillWithCoordinates(long[] outLat, long[] outLng) {
        if (outLat.length != nCities || outLng.length != nCities) {
            throw new IllegalArgumentException("outLat and outLng must be length " + nCities);
        }
        for (int i = 0; i < coords.length; ++i) {
            outLat[i] = coords[i].getX();
            outLng[i] = coords[i].getY();
        }
    }
}
