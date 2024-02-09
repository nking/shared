package algorithms.alignment;

import junit.framework.TestCase;
import org.w3c.dom.Text;

import java.util.ArrayList;
import java.util.List;

public class TextJustificationTest extends TestCase {

    public void test0() {

        // REference: "Dynamo: Amazon’s Highly Available Key-value Store"
        // DeCandia et al. 2007 (SOSP’07)
        String w = "Reliability at massive scale is one of the biggest challenges we face at Amazon.com, one of the largest e-commerce operations in the world; even the slightest outage has significant financial consequences and impacts customer trust. The Amazon.com platform, which provides services for many web sites worldwide, is implemented on top of an infrastructure of tens of thousands of servers and network components located in many datacenters around the world. At this scale, small and large components fail continuously and the way persistent state is managed in the face of these failures drives the reliability and scalability of the software systems.";

        String[] words = w.split("\\s+");

        int width = 72;

        TextJustification t = new TextJustification();

        List<String> out = new ArrayList<>();

        int cost = t.justify(words, width, out);
        int cost0 = 0;
        for (String ri : out) {
            System.out.println(ri);
            cost0 += (width - ri.length());
        }
        assertEquals(cost0, cost);
        System.out.println("cost0=" + cost0);
        System.out.println();

        System.out.println("words length=" + words.length);

    }

}
