package algorithms.alignment;

import junit.framework.TestCase;
import org.w3c.dom.Text;

import java.util.ArrayList;
import java.util.List;

public class TextJustificationTest extends TestCase {

    public void test1() {
        String w = "Google’s mission is to organize the world's information and make it universally accessible and " +
                "useful. That's why Search makes it easy to discover a broad range of information from a wide variety " +
                "of sources.";
        String[] words = w.split(" ");
        int width = 30;
        int expCost = 71;

        TextJustification.Best best3 = TextJustification.justify3(words, width);
        // count number of spaces at end of lines
        int endCount = 0;
        int cWIdx = 0;
        for (String line : best3.lines) {
            String[] ws = line.split(" ");
            for (String wd : ws) {
                assertEquals(words[cWIdx], wd);
                ++cWIdx;
            }
            assert(line.length() <= width);
            //System.out.printf("len=%d, width=%d\n    [%s]\n", line.length(), width, line);
            endCount += (width - line.length());
        }
        assertEquals(words.length, cWIdx);
        assertEquals(expCost, endCount);
        assertEquals(expCost, best3.added);
    }
    public void test0() {

        // Reference: "Dynamo: Amazon’s Highly Available Key-value Store"
        // DeCandia et al. 2007 (SOSP’07)
        String w = "Reliability at massive scale is one of the biggest challenges we face at Amazon.com, one of the " +
                "largest e-commerce operations in the world; even the slightest outage has significant financial " +
                "consequences and impacts customer trust. The Amazon.com platform, which provides services for many " +
                "web sites worldwide, is implemented on top of an infrastructure of tens of thousands of servers and " +
                "network components located in many datacenters around the world. At this scale, small and large " +
                "components fail continuously and the way persistent state is managed in the face of these failures " +
                "drives the reliability and scalability of the software systems.";

        String[] words = w.split("\\s+");

        final int width = 72;

        int expCost = 80;

        TextJustification t = new TextJustification();

        List<String> out = new ArrayList<>();

        int cost = t.justify(words, width, out);
        int cost0 = 0;
        for (String ri : out) {
            //System.out.println(ri);
            cost0 += (width - ri.length());
        }
        assertEquals(expCost, cost0);
        assertEquals(cost0, cost);
        //System.out.println("cost0=" + cost0);
        //System.out.println();

        //System.out.println("words length=" + words.length);

        TextJustification.Best best3 = TextJustification.justify3(words, width);
        // count number of spaces at end of lines
        int endCount = 0;
        int cWIdx = 0;
        for (String line : best3.lines) {
            String[] ws = line.split(" ");
            for (String wd : ws) {
                assertEquals(words[cWIdx], wd);
                ++cWIdx;
            }
            assert(line.length() <= width);
            //System.out.printf("len=%d, width=%d\n    [%s]\n", line.length(), width, line);
            endCount += (width - line.length());
        }
        assertEquals(words.length, cWIdx);
        assertEquals(expCost, endCount);
        assertEquals(expCost, best3.added);
    }

    /*
    a fun look at other missions while at it.

    Microsoft's mission is to empower people and organizations to achieve more. This mission is focused on providing
    people with the tools they need to make a positive impact.

    Meta's mission is to help people connect and build communities, and to create the technology that enables this.
    Meta's mission statement is "Building the future of human connection and the technology that makes it possible".

    The Cold Spring Harbor Laboratory (CSHL) mission is to improve human quality of life by advancing biomedical
    research and education. CSHL is a private, not-for-profit organization on Long Island, New York.    ...
    */

}
