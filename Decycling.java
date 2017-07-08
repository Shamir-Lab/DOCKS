import gurobi.GRBException;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Vector;

public class Decycling {
    
    public static void main(String[] argv) throws IOException, InterruptedException, GRBException {
            if (argv.length != 8 && argv.length != 7)
                    System.out.println("Usage: java -jar decycling.jar <outputfile> <k> <alphabet>");
            
            // Parse input
            String outfile = argv[0];
            int k = Integer.parseInt(argv[1]);
            String alphabet = argv[2];
            int alphabetSize = alphabet.length(); 
            PrintWriter writer = new PrintWriter(outfile);

            // Find k-mer hitting set and write to file
            Vector<Integer> d = decycling(k, alphabetSize, alphabet); // Generate decycling set
            System.out.println("Finished decycling " + d.size());
            deBruijn G = new deBruijn(k, alphabetSize, alphabet);
            for (int i = 0; i < d.size(); i++) {
            	writer.println(G.getEdgeLabel(d.get(i)));
            }

            writer.close();
    }
    
    public static Vector<Integer> decycling(int k, int alphabetSize, String alphabet) {

        Vector<Integer> d_set = new Vector<Integer>();

        // go through the cycles using fkm
        int []a = new int[k+1]; for (int i = 0; i <= k; i++) a[i] = 0;
        int i = 1; a[0] = -1;
          
        do {
        	if (k % (i) == 0) {
        		int d_node = getDecyclingEdge(a, i, d_set.size()+1, k, alphabetSize);
                d_set.add(d_node);
        	}
        	for (i = k; a[i]==alphabetSize-1; i--);
        	if (i != 0) {
        		a[i] = a[i] + 1;
        		for (int j = i+1; j <= k; j++) a[j] = a[j-i];
        	}
        } while (i != 0);
        return d_set;
    }
    
    public static int getDecyclingEdge(int[] a, int i, int d, int k, int alphabetSize) {
        double u = Math.PI * 2.0 / (double)k;
        int q;
        if (i < k)
        	q = k;
        else {
        	for (q = 1;;q++) {
        		double s = 0;
                for (int l =1; l <= k; l++) {
                	s += a[l]*Math.sin(((l-1+k-q)%k)*u);
                }
                if (s < .0001) break;
            }
            for (q++;q<k+k;q++) {
            	double s = 0;
                for (int l = 1; l <=k; l++)
                	s += a[l]*Math.sin(((l-1+k+k-q)%k)*u);
                if (s>=0.0001) break;
            }
            if (q>k) q-=k;
        }
        int code = 0;
        for (int j = q+1; j <=k; j++) code = alphabetSize*code + a[j];
        for (int j = 1; j <=q; j++) code =alphabetSize*code + a[j];

        return code;
    }
}
