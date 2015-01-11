import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

/*
 * Class to store details of edges
 *
 */
class Edge implements Comparable<Edge>
{
    int source, destination;
    int weight;

    /**
     * Constructor
     * @param src
     * @param dest
     * @param weight
     */
    public Edge(int src, int dest, int weight)
    {
        this.source = src;
        this.destination = dest;
        this.weight = weight;
    }
    
    /**
     * Source node
     * @return int
     */
    public int getSource()
    {
        return source;
    }
    
    /**
     * Destination Node
     * @return int
     */
    public int getDestination()
    {
        return destination;
    }
    
    /**
     * Edge weight
     * @return int
     */
    public int getWeight()
    {
        return weight;
    }
    
    /**
     * Adjacent vertex
     * @param v
     * @return int
     */
    public int other(int v){
    	int v1=this.source;
    	int v2=this.destination;
    	if(v==v1){
    		return v2;
    	}
    	else if(v==v2){
    		return v1;
    	}
    	else{
    		return -1;//not adjacent
    	}
    }
    
    /**
     * Comparing Edges 
     */
    public int compareTo(Edge edge)
    {
        // For comparisons during insertions and deletions
        int cost1 = this.weight;
        int cost2 = edge.weight;
        int from1 = this.source;
        int from2 = edge.source;
        int to1   = this.destination;
        int to2   = edge.destination;

        if (cost1==cost2 && from1==from2 && to1==to2)
          return(0);
        else 
          return (cost1>cost2)?1:-1;
    }
}

/**
 * Class to store graph structure
 * @author Abhinaya
 *
 */
class Graph{
	TreeSet<Edge> edges = new TreeSet<Edge>(); //Ordered set of edges
	int noOfNodes,noOfEdges;
	ArrayList<Edge> adjEdges[];

	ArrayList<Integer> nodes=new ArrayList<Integer>();//A list of vertices
	
	static int ApproxTour[];
	/**
	 * Constructor
	 * @param nodes
	 * @param edges
	 */
	Graph(int nodes, int edges){
		this.noOfNodes=nodes;
		this.noOfEdges=edges;
		ApproxTour = new int[noOfNodes+1];
		adjEdges=new ArrayList[noOfNodes];
	}
	
	/**
	 * Add an edge to graph
	 * @param from
	 * @param to
	 * @param cost
	 */
	public void addEdge(int from, int to, int cost) {
		Edge e=new Edge(from, to, cost);
		edges.add(e);  // Update set
		if(adjEdges[from] == null){
			adjEdges[from]= new ArrayList<Edge>();
		}
		if(adjEdges[to] == null){
			adjEdges[to]= new ArrayList<Edge>();
		}
		if(!adjEdges[from].contains(e))
			adjEdges[from].add(e);
		if(!adjEdges[to].contains(e))
			adjEdges[to].add(e);
		if(!nodes.contains(from)){
			nodes.add(from);
		}
		if(!nodes.contains(to)){
			nodes.add(to);
		}
	}
	
	/**
	 * Get set of edges
	 * @return TreeSet
	 */
	public TreeSet<Edge> getEdges() {
		return edges;
	}
	
	/**
	 * Get set of adjacent edges
	 * @param v
	 * @return ArrayList
	 */
	public ArrayList<Edge> getAdjEdges(int v) {
		return adjEdges[v];
	}	

	/**
	 * Get list of nodes
	 * @return ArrayList
	 */
	public ArrayList<Integer> getNodes() {
		return nodes;
	}
	
	/**
	 * Get edge weight
	 * @param from
	 * @param to
	 * @return cost
	 */
	public int getCost(int from , int to){
		for(Edge e:getEdges()){
			if((e.source==from && e.destination ==to) || (e.source==to && e.destination == from))
				return e.weight;
		}
		return 0;
	}
	

	/**
	 * Calculate TSP tour
	 * @param G
	 * @return int
	 */
	public static int[] APPROX_TSP_TOUR(Graph G){
		PrimsVertices p=new PrimsVertices(G);
		ApproxTour = p.primMST(G,1);
		return ApproxTour;
	}
}

/**
 * Prims Algorithm
 * @author Abhinaya
 *
 */
class PrimsVertices{

	int distance[];
	int visited[];
	Edge MST[];
	IndexMinPQ<Integer> pq ;
	int cost;
	int noOfNodes =0;
	int prev=-1;
	Graph net;
	int flag=1;
	int tour[];
	static int index =1;
	//Constructor
	PrimsVertices(Graph g){

		index=1;
		net =g;
		cost=0;
		noOfNodes = g.noOfNodes;
		pq = new IndexMinPQ<Integer>(noOfNodes);
		MST=new Edge[noOfNodes];
		tour = new int[noOfNodes+1];
		
	}
	
	/**
	 * Create MST
	 * @param g
	 * @param type = 1 for tour and 0 for MST
	 * @return cost
	 */
	int[] primMST(Graph g,int type){
		distance =new int[noOfNodes];
		visited =new int[noOfNodes];
		int i=0;//nodes visited
		for(int x : g.getNodes()){
			visited[x]=0;//false
			distance[x]=Integer.MAX_VALUE;
		}

		HashMap<Integer,ArrayList<Edge>> tree = new HashMap<Integer,ArrayList<Edge>>(); //root, all destinations
		
		int r= g.getNodes().get(0);//root
		
		for (int x : g.getNodes()) {
	        tree.put(x, new ArrayList<Edge>());//Create new arraylist for each node    
			if(visited[x]!=1){//Visit node if not visited
	            	prim(g,x);
	        }
		}
		
		//Calculate tour cost
		if(type==1){//tour
			//Add edges to tree
			for (i=0;i<MST.length;i++) {
				if(MST[i]!=null){
				 tree.get(MST[i].source).add(MST[i]);
				}
			 }
			
			cost =preorderTreeWalk(tree,r);
			tour[noOfNodes]=cost;
			return tour;	
		}
		else{//Return MST cost
			cost=0;
//			if(MST.length>0)
//				System.out.print("MST:");
			for (i=0;i<MST.length;i++) {
				if(MST[i]!=null){
					//System.out.print(MST[i].getSource()+"-"+MST[i].getDestination()+",");
				 	cost+=MST[i].getWeight();
				}
			 }
			// System.out.println("Cost="+cost);
			 int[] x = {cost};
			 return x;
		}
	}
		
	/**
	 * Prims Algorithm
	 * @param g
	 * @param vertex
	 */
	void prim(Graph g, int x) {
			distance[x] = 0;
	        pq.insert(x, distance[x]);//distance of node from itself is 0
	        while (!pq.isEmpty()) {
	            int v = pq.delMin();
	            scan(g,v);
	        }
	}
	
	/**
	 * Scan graph for unvisited nodes
	 * @param g
	 * @param v1
	 */
	void scan(Graph g, int v1) {
		        visited[v1] = 1;
		        for (Edge e : g.getAdjEdges(v1)) {
		            int v2 = e.other(v1);
		            if (visited[v2]==1) continue;         // node has been visited
		            if (e.weight < distance[v2]) {
		            	distance[v2] = e.weight;
		            	if(v2==e.source)
		            		MST[v2]=new Edge(e.destination,e.source,e.weight);
		            	else
		            		MST[v2] = e;
		                if (pq.contains(v2)) pq.decreaseKey(v2, distance[v2]);
		                else                pq.insert(v2, distance[v2]);
		            }
		        }
	}
	
	/**
	 * Detect tour based on MST
	 * @param tree
	 * @param root
	 * @return cost
	 */
	int preorderTreeWalk(HashMap<Integer,ArrayList<Edge>> tree,int root){
			 //for each node traverse
			 ArrayList<Edge> temp = new ArrayList<Edge>();
			 
			 temp = tree.get(root);
			 //int vertex= root+1;
			 //System.out.print("Tour: v"+vertex);
			 tour[0]=root+1;
			 tree.remove(root);
			 cost = 0;
			 
			 cost=getChild(tree,temp,cost);
			 return cost;	 
	}
		 
	/**
	 * Visit child nodes recursively
	 * @param tree
	 * @param temp
	 * @param cost
	 * @return cost
	 */
	int getChild(HashMap<Integer,ArrayList<Edge>> tree, ArrayList<Edge> temp,int cost){
			 int i;
			 int visited=0;
			 ArrayList<Edge> child = new ArrayList<Edge>();
			 
			 for(i=0;i<temp.size();i++){
				 visited = temp.get(i).destination;
				 
				if(i==0){
						 cost+= temp.get(i).weight;
				}
			    else{
			    	if(prev!=-1){
			    		 cost+= net.getCost(prev,temp.get(i).destination);
			    		
			    	}
			    	else{
			    		 cost+= net.getCost(temp.get(i-1).destination,temp.get(i).destination);
			    	}
			    	
				 }
				 //Print tour to file
				 tour[index++]=visited+1;
				 child = tree.get(visited);
				 flag++;
				 tree.remove(visited);
				 if(child.size()==0){
					 prev=visited;
				 }
  				 cost = getChild(tree,child,cost);
					
				 
			 }
			 if(flag==noOfNodes && temp.size()!=0){
				 flag=0;
				 cost+=net.getCost(net.getNodes().get(0),visited);
				 tour[noOfNodes-1]=visited+1;
			 }
			 
			 return cost;
			 
	 }
		 
}


