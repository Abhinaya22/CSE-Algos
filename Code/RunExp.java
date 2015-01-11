import java.io.*;
import java.util.ArrayList;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/*
This is an example of how your experiments should look like.
Feel free to use and modify the code below, or write your own experimental code, as long as it produces the desired output.

This code assumes you are running your program from the folder where the graph text instances are.

Elias Khalil for CSE 6140 - Fall 2014
*/
public class RunExp{
	public enum Type {BnB , Approx , Heur, LS1 , LS2};
	public static void main(String[] args) {

		TSP Tsp = new TSP();

        String input, output, SolFile, SolTrace="";
        
		//Initialise input parameters
        String filename = args[0];//"berlin52";//{"berlin52", "burma14", "ch150", "gr202", "kroA100", "ulysses16"};Integer.parseInt(args[0]);
        double cutOff = Double.parseDouble(args[1]);//in seconds
        String method = args[2];//"BnB";//BnB | Approx | Heur| LS1 | LS2 //Integer.parseInt(args[2]);
		int seed = Integer.parseInt(args[3]);
		output = "./output/"+filename+"_"+method+"_"+cutOff;
		input = "./DATA/"+filename+".tsp";
		Type m;
		
		//Thread
		ExecutorService executor = Executors.newSingleThreadExecutor();
        Future<String> future;

		
		
		if(method.equalsIgnoreCase("BnB")){
			SolFile = output+".sol";
			SolTrace = output+".trace";
			m = Type.BnB;
		}
		else if(method.equalsIgnoreCase("LS1")||method.equalsIgnoreCase("LS2")){
				SolFile = output+"_"+seed+".sol";
				SolTrace = output+"_"+seed+".trace";
				m = method.equalsIgnoreCase("LS1")?Type.LS1:Type.LS2;
				/*
				 *  each line has two values (comma-separated):
					1)    A timestamp in seconds (double)
					2)    Quality of the best found solution at that point in time (integer). 
					Note that to produce these lines, you should record every time a new improved solution is found.
				 * */
		}
		else{
			SolFile = output+".sol";
			m = method.equalsIgnoreCase("Approx")?Type.Approx:Type.Heur;
			/*
			 * File format:
				1.     line 1: quality of best solution found (integer)
				2.     line 2: list of vertex IDs of the tour (comma-separated): v_1,v_2,v_2,...,v_n,v_1
			 */
		}
		
		try {
			BufferedWriter outFile = new BufferedWriter(new FileWriter(SolFile));
			
			//Parse Input
			TSP tspData = Tsp.parseTSP(input,SolTrace,cutOff);
			int q; 
			double start, end, total;
			System.out.println("method:"+m);
			switch(m){
			
			case Approx:
				System.out.println("Entered");
				//MST 2 approximation - Prims algorithm
				//Creating graph
				Graph G = parseEdges(tspData);
				
				start = System.currentTimeMillis();
				int[] tour = Graph.APPROX_TSP_TOUR(G);
				end = System.currentTimeMillis();
				
				
				total = (end-start)/1000;
				System.out.println("Running time:"+total);
				
				//Write Quality to file
				q = (tour[tour.length-1]);//TO DO: Check formula
				outFile.write("Quality:"+q+"\nList Of Vertices on Tour:");
				//Write tour to file
				for(int i=0;i<tour.length-1;i++){
					outFile.write("v"+tour[i]+", ");
				}
				outFile.write("v"+tour[0]);
				break;
				
			case Heur:
				
				//Call Greedy Heuristic
				start = System.currentTimeMillis();
				TSP heur = Tsp.Heur();
				end = System.currentTimeMillis();
	
				total = (end-start)/1000;
				System.out.println("Running time:"+total);
				//Write Quality to file
				q = (heur.HeurCost);//TO DO: Make Integer
				outFile.write("Quality:"+q+"\nList Of Vertices on Tour:");
				//Write tour to file
				for(int i =0; i < heur.dim; i++) {
					outFile.write("v"+(heur.HeurTour[i]+1)+", ");
				}
				outFile.write("v"+(heur.HeurTour[0]+1));
				break;
				
			case LS1:
				
				//Call Hill Climbing algorithm
				start = System.currentTimeMillis();
				TSP hill = Tsp.Hill2(seed);
				end = System.currentTimeMillis();
				
				total = (end-start)/1000;
				System.out.println("Running time:"+total);
				//Write Quality to file
				q = (hill.HillCost2);//TO DO: Check formula
				outFile.write("Quality:"+q+"\nList Of Vertices on Tour:");
				//Write tour to file
				for(int i =0; i < hill.dim; i++) {
					outFile.write("v"+(hill.HillTour2[i]+1)+", ");
				}
				outFile.write("v"+(hill.HillTour2[0]+1));
				break;
				
			case LS2:
		        //Call Simulated Annealing
				start = System.currentTimeMillis();
				TSP sa = Tsp.SA(seed);
				end = System.currentTimeMillis();
				
				total = (end-start)/1000;
				System.out.println("Running time:"+total);
				//Write Quality to file
				q = (sa.SACost);//TO DO: Check formula
				outFile.write("Quality:"+q+"\nList Of Vertices on Tour:");
				//Write tour to file
				for(int i =0; i < sa.dim; i++) {
					outFile.write("v"+(sa.SATour[i]+1)+", ");
				}
				outFile.write("v"+(sa.SATour[0]+1));
				
				break;
				
			case BnB:
				//future =  executor.submit(new BranchAndBound(tspData));
				//Call Branch and Bound algorithm 
				BranchAndBound bnb = new BranchAndBound(tspData);
				int src = 0;
		    	ArrayList<Integer> route = new ArrayList<Integer>();
				//start = System.nanoTime();
		    	
				//try {
					start = System.currentTimeMillis();
				
				
					long starttime = System.currentTimeMillis();
					String bnbSol=bnb.search(src, route, starttime);

					end = System.currentTimeMillis();
					total = (end-start)/1000;
					System.out.println("Running time:"+total);
					//bnbSol = future.get((long) cutOff, TimeUnit.SECONDS);
					int i = bnbSol.indexOf(" Cost:");
					//Write Quality to file
					q = Integer.parseInt(bnbSol.substring(i+6));//TO DO: Check formula
					outFile.write("Quality:"+q+"\nList Of Vertices on Tour:");
					//Write tour to file
					outFile.write(bnbSol.substring(0,i-1));//Tour
					
//				} catch (InterruptedException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				} catch (ExecutionException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				} catch (TimeoutException e) {
//					// TODO Auto-generated catch block
//					e.printStackTrace();
//				}
		        //end = System.nanoTime();
		        //total = (end-start)/1000000;
		        
				
				
				break;
			
			}
						
			
				
			//Branch and Bound:
			//Initialise problem
//			TimerTask timerTask = new MyTimer("BnB",tspData);
//			// running timer task as daemon thread
//			Timer timer = new Timer(true);
//			timer.scheduleAtFixedRate(timerTask, 0, 10 * 1000);
//
//			System.out.println("Branch and Bound begins! :" + new Date());
//
//			// cancel after sometime
//			        try {
//			            Thread.sleep(20000);
//			        } catch (InterruptedException e) {
//			            e.printStackTrace();
//			        }
//			        timer.cancel();
//			        System.out.println("TimerTask cancelled! :" + new Date());
//			        try {
//			        	Thread.sleep(30000);
//			        } catch (InterruptedException e) {
//			            e.printStackTrace();
//			        }

			
	    	
	        
	        
	        outFile.close();
	        
			
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		//}//end of 1 for
	}//end of main


	
	 static Graph parseEdges(TSP tsp) throws IOException{
			
			//initialize graph
			Graph g=new Graph(tsp.dim,tsp.dim*tsp.dim);
			//Read each vertex
			for(int i=0;i<tsp.dim;i++)
				for(int j=0;j<tsp.dim;j++)
					g.addEdge(i,j,tsp.DisMatrix[i][j]);
			
			return g;
	}


}//end of RunExp


