import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;
import java.util.Arrays;

import gurobi.*;

public class deBruijn {
        
	public static int alphabetSize;
        public static String alphabet;
        public static HashMap<Character, Integer> mapAlphabet;
        public double [][] D;
        public double [][] F;
        public double [] hit;
        public int []sort;
        public int k;
        public int bits = 8;
        public int length = 0;
        public int numVertices = 0;
        public int numEdges = 0;
        public int pow, pow2, pow3;
        public int pow_1;
        public int pow_mask;
        public int log;
        public static byte[] E;
 
        public static int removeILP(int k, int L, int alphabetSize, int time, Vector<Integer> decycling, String file, 
        		PrintWriter out) throws IOException, GRBException {
            
            GRBEnv env = new GRBEnv("decycline_"+k+"_"+L+"_"+alphabetSize+".log");
            env.set(GRB.DoubleParam.TimeLimit, time);
            GRBModel model = new GRBModel(env);
            
            // Create variables
            int nume = (int)Math.pow(alphabetSize, k);
            int numv = (int)Math.pow(alphabetSize, k-1);
            GRBVar[] vertices = model.addVars(nume, GRB.BINARY); // Indicator variables if vertex is removed
            GRBVar[] Lvertices = model.addVars(numv, GRB.CONTINUOUS);
            model.update();
           
            // Objective
            GRBLinExpr expr = new GRBLinExpr();
            double[] coeffs = new double[nume];
            Arrays.fill(coeffs, 1.0);
            expr.addTerms(coeffs, vertices);
            model.setObjective(expr, GRB.MINIMIZE);
            
            // Set feasible solution
        	if (new File(file).exists()) {
            	BufferedReader br = new BufferedReader(new FileReader(file));
            	for(String line; (line = br.readLine()) != null; ) {
            		expr = new GRBLinExpr();
                    vertices[getInt(line, alphabetSize)].set(GRB.DoubleAttr.Start, 1.0);
            	}
        	}
        	
            // Constrains
            for (int i = 0; i < numv; i++) {
                        
            	// All paths are shorter than L-k
                expr = new GRBLinExpr();
                expr.addTerm(1.0, Lvertices[i]);
                model.addConstr(expr, GRB.LESS_EQUAL, L-k, " L" + i);
                expr.addTerm(1.0, Lvertices[i]);
                model.addConstr(expr, GRB.GREATER_EQUAL, 0, " L" + i);
                        
                // Make sure Lv >= Lu + 1 - (L-k+1) xv
                for (int j = 0; j < alphabetSize; j++) {
                	expr = new GRBLinExpr();
                    String u = getChar(j) + getString(i, k-1, alphabetSize);
                    int indexe = getInt(u, alphabetSize);
                    int indexv = getInt(u.substring(0, k-1), alphabetSize);
                    expr.addTerm(1.0, Lvertices[indexv]); expr.addTerm(-1.0, Lvertices[i]); expr.addTerm(k-1-L, vertices[indexe]);
                    model.addConstr(expr, GRB.LESS_EQUAL, -1.0, "L"+getString(i, k, alphabetSize)+j);
                }
            }
                
            // Add decycling
            for (int i = 0; i < decycling.size(); i++) {
            	expr = new GRBLinExpr();
                expr.addTerm(1.0, vertices[decycling.get(i)]);
                model.addConstr(expr, GRB.EQUAL, 1.0, "decycling"+i);
            }

            model.update();
            model.optimize(); 
            
            double[] xvals = model.get(GRB.DoubleAttr.X, model.getVars());
            int count = 0;
            for (int i = 0; i < nume; i++) {
                if (xvals[i] > 0) {
                	out.write(getString(i, k, alphabetSize) + "\n");
                    count++;
                }	
            }
            
            // Dispose of model and environment
            model.dispose();
            env.dispose();                   
            
            return count - decycling.size();
        }

        // Find a random edge that has not been deleted
        private int randomEdge(Random R) {
        	int index = -1;
            do {
            	index = R.nextInt(E.length);
            } while (!(E[index] == 1));
            return index;
        }
        
        // Remove random edges till no path of length = L-k+1 exists
        public int removeRandomRemaining(int L, int seed, PrintWriter writer) {
                Random R = new Random(seed);
                boolean length;
                int index = -1;
                int num = 0;
                sort = topologicalSort();
                do {
                	length = maxLength() >= L-k+1;
                    index = randomEdge(R);
                    if (index > -1 && length) {
                    	removeEdge(index);
                        writer.println(getEdgeLabel(index));
                        num++;
                    }
            } while (index > -1 && length);
            return num;
        }
        
        // Remove vertices with greatest hitting (all paths) number
        // till no path of length = L-k+1 exists
        public int removeAnyRemaining(int L, int x, PrintWriter writer) { 
            int []index;
            int num = 0;
            int l = L-k+1;
            sort = topologicalSort();
            D = new double[1][numVertices];         
            F = new double[1][numVertices];

            while(maxLength() >= l) {
            	calcAnyPaths();
                index = calcHitAny(x);

                for (int i = 0; i < index.length; i++) {
                	removeEdge(index[i]); 
                	writer.println(getEdgeLabel(index[i]));
                	num++;
                }
            } 
            return num;
        }
        
        public void calcAnyPaths() {
        	if (alphabetSize == 4)
            calcAnyPaths4();
        	else calcAnyPathsK();
        }

        // Calculate hitting (all paths) numbers for all vertices
        // assuming alphabet size is 4 (relevant to DNA)
        public double calcAnyPaths4() {

            for (int i = 0; i < pow; i++) {F[0][i] = D[0][i] = 1;}
          
            for (int i = 0; i < sort.length; i++) {
            	D[0][sort[i]] += E[sort[i]]*D[0][(sort[i] >> 2)] + E[sort[i] + pow]*D[0][((sort[i] + pow) >> 2)] + 
            			E[sort[i] + pow2]*D[0][((sort[i] + pow2) >> 2)] + E[sort[i] + pow3]*D[0][((sort[i] + pow3) >> 2)];
           	 	int ind = (sort[sort.length-i-1] * 4); 
           	 	F[0][sort[sort.length-i-1]] += E[ind]*F[0][ind & pow_mask] + E[ind+1]*F[0][(ind+1) & pow_mask] + E[ind+2]*F[0][(ind+2) & pow_mask] + 
             		E[ind+3]*F[0][(ind+3) & pow_mask];
            }
            
            return 0;
        }
        
        // Calculate hitting (all paths) numbers for all vertices
        public double calcAnyPathsK() {

            for (int i = 0; i < pow; i++) {F[0][i] = D[0][i] = 1;}
            
            for (int i = 0; i < sort.length; i++) {
            	for (int l=0; l < alphabetSize; l++)
            	D[0][sort[i]] += E[sort[i] + l*pow]*D[0][((sort[i] + l*pow) / alphabetSize)];
            	int ind = (sort[sort.length-i-1] * alphabetSize); 
            	for (int l=0; l < alphabetSize; l++)
            		F[0][sort[sort.length-i-1]] += E[ind+l]*F[0][(ind+l) % pow];
            }
            
            return 0;
        }
        
        // Find vertex with largest hitting (any path) number
        public int[] calcHitAny(int x) {
        	double[] max = new double[x];
        	int[] imax = new int[x];
            for (int i = 0; i < length; i++) {
                double num = E[i]*F[0][i % pow]*D[0][i / alphabetSize];
                int j = 0;
                while (j < x && num <= max[j]) j++;
                if (j < x && num > max[j]) {
                	max[j] = num; imax[j] = i;}
            }

            return imax;
        }
    
        // Remove vertices with greatest hitting (l-long paths) number
        // till no path of length = L-k+1 exists
        public int removeRemaining(int L, PrintWriter writer) { 
            int index = -1;
            int num = 0;
            int l = L-k+1;
            D = new double[l+1][numVertices];
            F = new double[l+1][numVertices];
            
            while (calcPaths(l) && (index = calcHit(l)) >= 0) {
                    if (index > -1) {
                            removeEdge(index);
                            writer.println(getEdgeLabel(index));
                            num++;
                   }
            }
            sort = topologicalSort();
            System.out.println("MAX length = " + maxLength());
            return num;
        }
    
        // Remove vertices with greatest hitting (l-long paths) number
        // till no path of length = L-k+1 exists
        // use less memory at the expense of more running time
        public int removeMemoryRemaining(int L, int l, PrintWriter writer) throws InterruptedException { 
        	int index = -1;
        	boolean length;
        	int num = 0;
        	D = new double[2][numVertices];
        	F = new double[2][numVertices];

        	do {
        		length = maxLength() >= L-k+1;
        		if (length) {
                calcMemoryHit(L-k+1, l);
        		index = findMax();
                if (index > -1 && length) {
                        removeEdge(index);
                        writer.println(getEdgeLabel(index));
                        num++;
                	}
        		}
        	} while (index > -1 && length);
        	return num;
        }

        // Find edge hitting max L-long paths
        public boolean calcPaths(int L) {

        	if (alphabetSize == 4)
        		calcPaths4(L);
        	else {
        		calcPathsK(L, log);
        	}
        return true;
        }

     // Find vertex with largest hitting (l-long path) number
    public int calcHit(int L) {
    	double max = 0; int imax = -1;

        for (int i = 0; i < length; i++) {

            double num = 0;
            for (int j = (1-E[i]) * L; j < L; j++) {
            	num = num + F[j][i % pow]*D[L-j-1][i / alphabetSize];
            }

            if (num > max) {max = num; imax = i;}
        }

        return imax;
    }

    // Calc hitting number for all vertices trading memory for runtime
    public void calcMemoryHit(int L, int l) throws InterruptedException {
    	// Initialize
    	for (int i = 0; i < E.length; i++) hit[i] = 0;

    	for (int j = 0; j < L-1; j++) {
    		System.out.print("\t" + j + "\t");
    		calcMemoryFinishPaths(j, l);

   		for (int i = 0; i < E.length; i++)
                if (E[i] == 1) {
                        int source = i % pow;
                        int target = i / alphabetSize;

                       	hit[i] += F[0][source]*D[0][target];
                }
        }
    }
    
    // Find vertex with maximum hitting number
    public int findMax() {

    	double max = 0;
    	int maxi = -1;
        for (int i = 0; i < E.length; i++) {

            	if (hit[i]*E[i] > max) {
                	max = hit[i]; maxi = i;
                }
        }
        return maxi;
    }

    // Calc number of l-long paths ending at each vertex
    // trading runtime for memory
    	public double calcMemoryFinishPaths(int L, int t) throws InterruptedException {

    	Arrays.fill(F[0], 1);
    	for (int j = 1; j <= L; j++) {

    	Arrays.fill(F[1], 0);
        	for (int i = 0; i < F[0].length; i++) {

        		int ind = i * alphabetSize;
                for (int l = 0; l < alphabetSize; l++) {
                	if (E[ind] == 1) F[1][i] = F[1][i] + D[0][ind % pow];
                    ind += 1;
                }
            }

            F[0] = Arrays.copyOf(F[1], F[0].length);
    	}
    
    return 0;
    	}

    	// Calculate hitting (l-long paths) numbers for all vertices
        // assuming alphabet size is 4 (relevant to DNA)
    	public double calcPaths4(int L) {
    		for (int i = 0; i < pow; i++) {F[0][i] = D[0][i] = 1;}
        
    		for (int j = 1; j <= L; j++) {
    			for (int i = 0; i < pow; i++) {
    				int ind = (i * 4);                       
                    F[j][i] = E[ind]*F[j-1][ind & pow_mask] + E[ind+1]*F[j-1][(ind+1) & pow_mask] + E[ind+2]*F[j-1][(ind+2) & pow_mask] + 
                    		E[ind+3]*F[j-1][(ind+3) & pow_mask];
                    D[j][i] = E[i]*D[j-1][(i >> 2)] + E[i + pow]*D[j-1][((i + pow) >> 2)] + 
                    		E[i + pow2]*D[j-1][((i + pow2) >> 2)] + E[i + pow3]*D[j-1][((i + pow3) >> 2)];
                }
    		}
        
        return 0;
    	}

    	// Calculate hitting (l-long paths) numbers for all vertices
    	public double calcPathsK(int L, int log) {
    		for (int i = 0; i < pow; i++) {F[0][i] = D[0][i] = 1;}
    
    		for (int j = 1; j <= L; j++) {
    			for (int i = 0; i < pow; i++) {
    				int ind = (i * alphabetSize); F[j][i] = 0; D[j][i] = 0;
                    for (int l = 0; l < alphabetSize; l++){
                    	F[j][i] += E[ind+l]*F[j-1][(ind+l) % pow];
                    	D[j][i] += E[i + l*pow]*D[j-1][((i + l*pow) / alphabetSize/*>> log*/)];
                    }
    			}
    		}

    		return 0;
    	}
 
    	public void removeEdge(int i) {
    		if (E[i] == 1) numEdges--;
    		E[i] = 0;
    	}

    	// Initialize a de Bruijn object
    	public deBruijn(int _k, int _alphabetSize, String _alphabet) {
    		alphabetSize = _alphabetSize; k = _k; alphabet = _alphabet;
    		generateGraph(k);
    		pow = (int)Math.pow(alphabetSize, k-1);
    		pow2 = 2*pow; pow3 = 3*pow;
    		pow_mask = pow-1;
    		pow_1 = (int)Math.pow(alphabetSize, k-2);
    		log = (int)(Math.log(alphabetSize)/Math.log(2));
    		mapAlphabet = new HashMap<Character, Integer>();
    		for (int i = 0; i < alphabetSize; i++)
    			mapAlphabet.put(alphabet.charAt(i), i);
    	}

    	// Initialize a de Bruijn graph
    	public void generateGraph(int k) {
    		int num = (int)Math.pow(alphabetSize, k);
    		E = new byte[num];
    		length = num;
    		Arrays.fill(E, (byte)1);
    		numEdges = num;
    		numVertices = num / alphabetSize;
    	}

    	public String getEdgeLabel(int i) {
    		return getString(i, k, alphabetSize);
    	}

    	public static String getString(int a, int k, int alphabetSize) {
    		String rc = "";
    		for (int i = 0; i < k; i++) {
    			rc = getChar(a % alphabetSize) + rc;
    			a = a / alphabetSize;
    		}
    		return rc;
    	}

    	public static int getInt(String a, int alphabetSize) {
    		int rc = 0;
    		int k = a.length();
    		for (int i = 0; i < a.length(); i++)
    			rc += Math.pow(alphabetSize, k-i-1) * getInt(a.charAt(i));

    		return rc;
    	}

    	public static char getChar(int i) {
    		return alphabet.charAt(i);
    	}

    	public static int getInt(char a) {
    		return mapAlphabet.get(a);
    	}

    	public int[] getAdj(int v) {
    		int count = 0; int[] adj = new int[alphabetSize];
    		int pow = (int)Math.pow(alphabetSize, k-1);
    		for (int i = 0; i < alphabetSize; i++) {
    			int ind = v + i * pow; 
    			if (E[ind] == 1) 
    				adj[count++] = ind / alphabetSize;
    		}
    		int []rc = new int[count];
    		for (int i =0; i < count; i++) {
    			rc[i] = adj[i]; 
    		}
    		return rc;              
    	}

    	private int dfs(boolean[] used, boolean[] finished, int[] res, int ind, int u) {
    		used[u] = true; boolean cycle = false;

    		for (int v : getAdj(u)) {
    			// Return true for cycle
    			if (used[v] && !finished[v]) cycle = true;
    			if (!used[v]) {
    				ind = dfs(used, finished, res, ind, v);
    				cycle = cycle || (ind == -1);
    			}
    		}
    		finished[u] = true;

    		res[ind] = u;
    		if (cycle) return -1;
    		else return ind+1;
    	}

    	public int[] topologicalSort() {
    		int n = numVertices;

    		boolean[] used = new boolean[n];
    		boolean[] finished = new boolean[n];
    		int[] res = new int[n];
    		int ind = 0;
    		for (int i = 0; i < n; i++) {

    			if (!used[i]) {
    				ind = dfs(used, finished, res, ind, i);

    				if (ind == -1) return null;
    			}
    		}

    		int[] rc = new int[res.length];
    		for (int i = 0; i < rc.length; i++)
    			rc[i] = res[res.length-i-1];
    		return res;
    	}

    	public int maxLength() {

    		int[] depth = new int[sort.length];

    		int max = -1;
    		for (int i = 0; i < sort.length; i++) {
              
    			int maxdepth = -1;
    			for (int j = 0; j < alphabetSize; j++) {
    				int edgeind = sort[i] + j*pow;
    				int vertexind = edgeind / alphabetSize; // >> log;
    				if (depth[vertexind] > maxdepth && E[edgeind] == 1) maxdepth = depth[vertexind];
    			}

              depth[sort[i]] = maxdepth + 1;
              if (depth[sort[i]] > max) {max = depth[sort[i]];}
    		}
    		return max;
    	}
}
