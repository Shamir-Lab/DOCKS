import gurobi.GRBException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Vector;

public class DOCKS {

      
        public static void main(String[] argv) throws IOException, InterruptedException, GRBException {
                if (argv.length != 8 && argv.length != 7)
                        System.out.println("Usage: java -jar decycling.jar <outputfile> <inputfile> <k> <L> <alphabet> <ILP time limit> <0/1/2/3-random/greedyL/greedyAny/none> <x-optional with any>");
                
                // Parse input
                String outfile = argv[0];
                String infile = argv[1];
                int k = Integer.parseInt(argv[2]);
                int L = Integer.parseInt(argv[3]);
                String alphabet = argv[4];
                int alphabetSize = alphabet.length(); 
                Vector<Integer> e = new Vector<Integer>();
                PrintWriter writer = new PrintWriter(outfile);
                int time = Integer.parseInt(argv[5]);
                int random = Integer.parseInt(argv[6]);
                int x = 1;
                if (argv.length > 7) x = Integer.parseInt(argv[7]);

                // Find k-mer hitting set
                Vector<Integer> d = readKmers(infile, alphabet, alphabetSize);
                System.out.println("Finished reading decycling " + d.size());
                d.addAll(e);
                pathCover(k, L, d, d.size() > 0 ? 1 : 0, random, x, alphabetSize, time, alphabet, argv[0], writer); // Find L-path cover
                writer.flush();

                writer.close();
        }    
       
        public static int getInt(String a, HashMap<Character, Integer> map, int alphabetSize) {
        	int rc = 0;
        	for (int i = 0; i < a.length(); i++)
        	        rc += Math.pow(alphabetSize, a.length()-i-1) * map.get(a.charAt(i));
        	return rc;
        }
        
     
        public static Vector<Integer> readKmers(String file, String alphabet, int alphabetSize) throws IOException {
        	Vector<Integer> rc = new Vector<Integer>();
        	HashMap<Character, Integer> map = new HashMap<Character, Integer>();
        	for (int i = 0; i < alphabet.length(); i++) map.put(alphabet.charAt(i), i);
        	if (new File(file).exists()) {
        	BufferedReader br = new BufferedReader(new FileReader(file));
        	    for(String line; (line = br.readLine()) != null; ) {
        	        rc.add(getInt(line, map, alphabetSize));
        	    }
            	br.close();
        	}
        	
        	return rc;
        }
                
        public static int pathCover(int k, int L, Vector<Integer> d, int decycling, int random, int x, int alphabetSize, int time,
                        String alphabet, String file, PrintWriter writer) throws IOException, InterruptedException,GRBException {

        	// Remove decycling set
            int dsize = d.size();
            deBruijn G = new deBruijn(k, alphabetSize, alphabet);

            for (int i = 0; i < d.size() && decycling == 1; i++) {
            	G.removeEdge(d.get(i));
            }
                
            // Remove L-path cover
            int rsize = 0;

            if (random == 0)
            	rsize = G.removeRandomRemaining(L, 0, writer);
            else if (random == 1)
            	rsize = G.removeRemaining(L, writer);
            else if (random == 4)
            	rsize = G.removeMemoryRemaining(L, random, writer);
            else if (random == 2)
            	rsize = G.removeAnyRemaining(L, x, writer);
            writer.flush();
                        
            if (time > 0) {
            	PrintWriter out2 = new PrintWriter(file+".ilp"); 
                rsize = deBruijn.removeILP(k, L, alphabetSize, time, d, file, out2); out2.flush(); out2.close();}
                String output = (k + "\t" + L + "\t" + dsize + "\t" + rsize + "\t" + 
                		(rsize+dsize) + "\t" + (double)(rsize+dsize) / (double)dsize);
                System.out.println(output);
            return rsize+dsize;
        }
}