import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class BranchAndBound {

	ArrayList<Integer> initialRoute, optimumRoute;
    int nodes = 0;
    static int cost = 0;
    static int optimumCost = Integer.MAX_VALUE; //Set cost of optimal solution to very high value 
    int[][] DisMatrix;
	int[][] OrderMatrix;
	int size;
    int sourceCity;
    String result = new String();
    String optResult;//Store tour and cost
    String tracefile ;
    double endtime;


    ArrayList<Integer> origProb = new ArrayList<Integer>();

    
	/**
	 * Constructor
	 */
	public BranchAndBound(TSP data){

		


		DisMatrix = data.DisMatrix;
		OrderMatrix = data.OrderMatrix;
		size= DisMatrix.length; 
		
		optResult = new String();//Store tour and cost

		sourceCity = 0;//First node by default
		//Create original list of vertices
		for(int x=0;x<size;x++){
			origProb.add(x);
    	}
		this.tracefile=data.tracefile;
		endtime = data.cutoff + Math.round(System.currentTimeMillis()/10.0) / 100.0;//Set cutoff
	}
	
	
    /**
     * Sorting by lowerBound
     */
    public static <K, V extends Comparable<? super V>> Map<K, V>  sortByValue( Map<K, V> map ){
	    List<Map.Entry<K, V>> list =   new LinkedList<Map.Entry<K, V>>( map.entrySet() );
	    Collections.sort( list, new Comparator<Map.Entry<K, V>>()
	    {
	        public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
	        {
	            return (o1.getValue()).compareTo( o2.getValue() );
	        }
	    } );
	
	    Map<K, V> result = new LinkedHashMap<K,V>();
	    for (Map.Entry<K, V> entry : list)
	    {
	        result.put( entry.getKey(), entry.getValue() );
	    }
	    return result;
    }
    
    /**
     * Search for new solution
     * @param from node where we start the search.
     * @param route followed route for arriving to node "from".
     */
    public String search (int from, ArrayList<Integer> route, long starttime) {
    	//System.out.println("Entered"+endtime);
    	double ct=Math.round(System.currentTimeMillis()/10.0) / 100.0;
    	if( endtime - ct <= 0){
    		//System.out.println("Exiting"+ct);
    		return optResult;
    	}

    	//Initialise new subproblems for each possible edge
        Map<Integer, Integer> ProbSet = new HashMap<Integer,Integer>();//to node,lowerBound -- sorted on lowerBound
        
    	//Add original source
    	if(route.size()==0){
    		route.add(from);//Add source node
    	}

    	// Solution found:All nodes visited - Set new optimal
    	if (route.size() == size) {
        	newOpt(route,from, starttime);
        }
        else {
        	//Identify nodes yet to be explored
        	ArrayList<Integer> remNodes = (ArrayList<Integer>) origProb.clone();
        	remNodes.removeAll(route);
        	ProbSet = (Map<Integer, Integer>) identifySubProblems(from,remNodes);
	        ProbSet=sortByValue(ProbSet);
	        if (ProbSet.size()!=0){
	        	calcOpt(ProbSet,from,route, starttime);
	        }
	    }
		return optResult;
        
    }
    
    /**
     * Find sub problems
     */
    HashMap identifySubProblems(int from, ArrayList<Integer> remNodes){
    	HashMap<Integer, Integer> ProbSet = new HashMap<Integer,Integer>();//to node,lowerBound -- sorted on lowerBound
    	for (int to:remNodes){ //Node not visited
	    	// update the route's partial cost
			cost += DisMatrix[from][to];//Check cost on Adding new node "to"
			int lowB =cost;
			if (cost < optimumCost) { //check with upper bound
	            ArrayList<Integer> remN = (ArrayList)remNodes.clone();
	            remN.remove(remN.indexOf(to));
	
	            if(remN.size()>0){//Not leaf node
	            	//Check lower bound
	        		lowB = calcLowerBound(cost, from,to,remN);
	            }
	        	if(lowB < optimumCost){//Prune rest
	        			ProbSet.put(to, lowB);//Add to subproblems
	        	}
			}
        // update the route's cost (back to the previous value) - backtrack
        cost -= DisMatrix[from][to];
        }
    	return ProbSet;
    	
    }
    
    /**
     * Find new opt by depth first search
     */
     void calcOpt(Map<Integer, Integer> PS,int from,ArrayList route, long starttime){
    	 
    	 for(int to:PS.keySet()){//For each sub problem
         	 cost += DisMatrix[from][to];//Check cost on Adding new edge
             ArrayList<Integer> modRoute = (ArrayList)route.clone(); 
             modRoute.add(to);
             search(to,modRoute, starttime);
             cost -= DisMatrix[from][to];//Check cost on Adding new edge
         }
     }
     
     /**
      * Add new opt
      */
     
     void newOpt(ArrayList route,int from, long starttime){
    	 
    	 // update the route's cost
         cost += DisMatrix[from][sourceCity];//b to a
         if (cost < optimumCost) {
        	 
        	 try {
        		 optResult="";
	        	 BufferedWriter outTrace = new BufferedWriter(new FileWriter(tracefile,true));//Append mode
	        	 
	        	 optimumCost = cost;
	             optimumRoute = (ArrayList<Integer>)route.clone();
	             optimumRoute.add(sourceCity);
	             for (int i : optimumRoute){
	            	 int v = i+1;
	            	 optResult+= "v"+v+", ";
	             }
	             optResult+= "Cost:"+optimumCost;//Store cost as last element
	             
	             //Write timestamp and quality to trace file
	             result = optimumRoute.toString() + ", Cost:"+optimumCost + "\n";
	             //System.out.println(result);
		long currtime = System.currentTimeMillis();
		double timestamp = (currtime- starttime) /1000.0 ;
		outTrace.write(String.format("%.2f", timestamp) +" "+optimumCost+ "\n");
	             //outTrace.write((new Date()).toString()+" "+optimumCost+ "\n");
			 	 outTrace.close();
			 	
             } catch (IOException e) {
 				// TODO Auto-generated catch block
 				e.printStackTrace();
 			}
             
         }
         // update the route's cost (back to the previous value) - backtrack
         cost -= DisMatrix[from][sourceCity];
         
     }
    
   /**
    * Returns lower bound for a partial solution a-route-b
    * @param partialCost
    * @param route
    * @return
    */
    int calcLowerBound(int partialCost, int from, int to, ArrayList<Integer> remNodes){
    	int c = partialCost;
    	int shortestEdgeCostS = Integer.MAX_VALUE; 
    	int shortestEdgeCostT = Integer.MAX_VALUE;
    	int src= from;//a
    	int trgt = to;//b
    	for(int i:remNodes){//Shortest edge out of src
    		int edgeCostS = DisMatrix[src][i];
			if(shortestEdgeCostS> edgeCostS){
				shortestEdgeCostS = edgeCostS; //Save as shortest edge
			}
			int edgeCostT = DisMatrix[trgt][i];
			if(shortestEdgeCostT> edgeCostT){
				shortestEdgeCostT = edgeCostT; //Save as shortest edge
			}
    	}
    	
    	c += (shortestEdgeCostS>shortestEdgeCostT)?shortestEdgeCostT:shortestEdgeCostS;//Add costs of shortest outgoing edges
    	
    	int noOfRemNodes = remNodes.size();
    	
    	Graph g = new Graph(noOfRemNodes,noOfRemNodes*noOfRemNodes);
    	int m=0,n=0;
    	for(int i : remNodes){
    		for(int j:remNodes){
    			g.addEdge(m,n++,DisMatrix[i][j]);	
    		}
    		n=0;
    		m++;
    	}

    	PrimsVertices p=new PrimsVertices(g);
    	int[] x = p.primMST(g, 0);
    	c+=x[0];//Add cost of MST

    	return c;
    	
    }

    
}
