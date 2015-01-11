import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Arrays;
import java.util.concurrent.*; 
import java.util.Random;





public class TSP{

		int type= 0;
		int dim = 0;

		int OptCost = 0;
		int HeurCost = 0;
		int MSTCost = 0;
        int SACost = 0;
		int HillCost = 0;
		int HillCost2 = 0;


		int[][] DisMatrix;
		int[][] OrderMatrix;

		int[] HeurTour;
        int[] SATour;
		
		int[] MSTTour;
		int[] HillTour;
		int[] HillTour2;


		
		public String tracefile = "";
		public double cutoff;		
		public TSP() {}


/////////////////////////////////////////////////// Beginning of parseTSP /////////////////////////////////////////////////////////////






		public TSP parseTSP(String s,String tracefile,double cutOff){
			
			this.tracefile = tracefile;
			this.cutoff = cutOff;
		
			try
			{
				//I only put the one filename here 
				// the filename must be an input

				FileReader reader = new FileReader(s);
				BufferedReader br = new BufferedReader(reader);
				String row;
				//what we may need
				
				
				
				while((row = br.readLine()) != null )
				{
					//Get the dimension
					if(row.indexOf("DIMENSION") != -1)
					{
						String[] s1 = row.split(" ");
						this.dim = Integer.parseInt(s1[1]);
						
						System.out.println("dimension:"+this.dim);
					}
					
					//Get the type of data
					if(row.indexOf("EDGE_WEIGHT_TYPE") != -1){
						if(row.indexOf("EUC_2D") != -1){
							this.type = 1;			
						} else{
							this.type = 2;
						}
					}
					
					//Get the optimal cost
					if(row.indexOf("OPTIMAL_COST") != -1)
					{
						String[] s2 = row.split(" ");
						this.OptCost = Integer.parseInt(s2[1]);
						System.out.println("OptCost:"+this.OptCost);
					}

					//stop the loop
					if (row.indexOf("NODE_COORD_SECTION") != -1) { break; }
				}
				
				//Read the points from the file
				double[][] Dpoints = new double[this.dim][2];
				String[] tmp;
				int i=0;
				while((row = br.readLine()) != null)
				{
					if(row.indexOf("EOF") != -1)
					{	break;
					} else{
						//get different points matrix
						if(this.type == 1)
						{
							tmp = row.split(" ");
							//System.out.println(tmp[1]);
							//System.out.println(tmp[2]);
							Dpoints[i][0] = Double.parseDouble(tmp[1]);
							Dpoints[i][1] = Double.parseDouble(tmp[2]);
							i++;
						} else{
							
							tmp = row.split("\\s+");
							//System.out.println(tmp[2]);
							//System.out.println(tmp[3]);
							Dpoints[i][0] = Double.parseDouble(tmp[2]);
							Dpoints[i][1] = Double.parseDouble(tmp[3]);
							i++;
						}
					}
				}


		/*	for(int j=0;j<this.dim;j++)
			{
				System.out.println(Dpoints[j][0]+","+Dpoints[j][1]);
			}*/


			//get the distance matrix


				//get the distance matrix
				this.DisMatrix = new int[this.dim][this.dim]; 




				if(this.type == 1)
				{
					double sum;

					for (int j = 0; j < this.dim; j++)
					{
						for(int k = 0; k < j+1; k++)
						{
							if (j == k){	
								this.DisMatrix[j][k] = 0;
							}else {
								sum = Math.pow((Dpoints[j][0] - Dpoints[k][0]), 2) + Math.pow((Dpoints[j][1] - Dpoints[k][1]), 2);
								//DisMatrix[j][k] = DisMatrix[k][j] = (int)Math.sqrt(sum);
								this.DisMatrix[j][k] = (int)(0.5+Math.sqrt(sum));
								this.DisMatrix[k][j] = this.DisMatrix[j][k];
							}
						}
					}
				}else {

					for(int j = 0; j < this.dim; j++)
					{
						for(int k = 0; k < j+1; k++)
						{
							if (j == k){
								this.DisMatrix[j][k] = 0;
							}else {
								this.DisMatrix[j][k] = (int)GetDistance(Dpoints[j][0], Dpoints[j][1], Dpoints[k][0], Dpoints[k][1]);
								this.DisMatrix[k][j] = this.DisMatrix[j][k];
							}	
						}
					}
				}


	
				//get the ordered matrix
				this.OrderMatrix = new int[this.dim][this.dim];

				int[] TmpArray = new int[this.dim];

				for (int k = 0; k < this.dim; k++)
				{
					HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
					for(int j = 0; j < this.dim ; j++)
					{
						map.put(this.DisMatrix[k][j], j);
						TmpArray[j] = this.DisMatrix[k][j];
					}
					Arrays.sort(TmpArray);
					for(int j = 0; j < this.dim; j++)
					{
						this.OrderMatrix[k][j] = map.get(TmpArray[j]);
					}
				}






















				
				/*for(int j = 0; j < this.dim; j++)
				{
					for(int k = 0; k < this.dim;k++)
						{
						System.out.print(this.OrderMatrix[j][k]+" ");
					}
					System.out.print("\n");
				}

				for(int j = 0; j < this.dim; j++)
				{
					for(int k = 0;k < this.dim;k++)
					{
						System.out.print(this.DisMatrix[j][k]+" ");
					}
					System.out.print("\n");
				} */



			} ////// End of try /////////////////
			catch(IOException e)
			{
				e.printStackTrace();
			}


			return this;
			
		}

//////////////////////////////////////////////////// End of parseTSP ///////////////////////////////////////////////////////



////////////////////////////////////////////////// Beginning of Heuristic Greedy algorithm ////////////////////////////////


		public TSP Heur(){

			this.HeurTour = new int[this.dim];
			this.HeurCost = 0;

			

			int tmpNode = 0;
			int tmpcost = 0;

			int InsertedNum = 0;

			int InsertPosition = 0;
			int[] tmpCost = new int[this.dim];

			boolean[] InsertedNodes = new boolean[this.dim]; 

			//initialization and insert the first vertex: v1
			for(int i=0; i < this.dim; i++) {
				this.HeurTour[i] = 0;
				InsertedNodes[i] = false; 
				tmpCost[i] = 0;}

			InsertedNodes[0] = true;
			InsertedNum++;

			//insert the second vertex which is farthest to v1
			tmpcost = this.DisMatrix[0][1];
			tmpNode = 1;

			for(int i = 2; i < this.dim; i++) {
				if(tmpcost < this.DisMatrix[0][i]){
					tmpcost = this.DisMatrix[0][i];
					tmpNode = i;
				}
			}

			this.HeurTour[InsertedNum] = tmpNode;
			this.HeurCost = this.HeurCost + 2 * tmpcost;
			InsertedNodes[tmpNode] = true;
			InsertedNum++;

			


			while(InsertedNum != this.dim) {

			//Selection step: find node not in the sub-tour farthest from any node in the sub-tour

				int m1 = 0;
				for(int i = 0; i < this.dim; i++) {
					if(InsertedNodes[i] == false) {//the vertex v has not been inserted in the tour

							tmpcost = this.DisMatrix[this.HeurTour[0]][i];

							for(int j = 1; j < InsertedNum; j++) {

								int t = this.DisMatrix[this.HeurTour[j]][i];

								if(tmpcost > t) { tmpcost = t; }/// End of if 2

							}/// End of for 2


							tmpCost[i] = tmpcost;

					}/// End of if 1
				}//// End of for 1

				int m = 0;

				for(int i = 0; i < this.dim; i++) {
					if(InsertedNodes[i] == false) {// the vertex v has not been inserted in the tour
						if(m == 0) {
							tmpNode = i;
							tmpcost = tmpCost[i];
							m = 1;
						}else if(tmpcost < tmpCost[i]){ 
							tmpNode = i;
							tmpcost = tmpCost[i];
						}/// End of if 2


					}/// End of if 1


				}/// End of for 1


			//Insertion step: find the arc (i, j) in the sub-tour which minimizes C(i, r) + c(r, j) - c(i, j)


				



				tmpcost = this.DisMatrix[this.HeurTour[InsertedNum-1]][tmpNode] + this.DisMatrix[this.HeurTour[0]][tmpNode] - this.DisMatrix[this.HeurTour[InsertedNum-1]][this.HeurTour[0]];
				InsertPosition = InsertedNum; // Insert v between v1 and v(InsertedNum) if possible

				for(int j = 0; j < InsertedNum-1; j++) {
					int t = this.DisMatrix[this.HeurTour[j]][tmpNode] + this.DisMatrix[this.HeurTour[j+1]][tmpNode] - this.DisMatrix[this.HeurTour[j]][this.HeurTour[j+1]];

					if(tmpcost > t) {
						tmpcost = t;
						InsertPosition = j+1; // Insert v between vj and v(j+1) if possible
					}/// End of if 1
				}/// End of for 1


				
				
			//Insertion step: insert r between i and j
				

				if(InsertPosition == InsertedNum) {
					this.HeurTour[InsertedNum] = tmpNode;
					this.HeurCost = this.HeurCost + tmpcost;
				}else {
					for(int i = InsertedNum; i > InsertPosition; i--) {
						this.HeurTour[i] = this.HeurTour[i-1];
					}

					this.HeurTour[InsertPosition] = tmpNode;
					this.HeurCost = this.HeurCost + tmpcost;
				}

				InsertedNodes[tmpNode] = true;
				InsertedNum++;


			}/// End of while
			
			/*System.out.println("HeurTour:");
			for(int i =0; i < this.dim; i++) {
				System.out.print("v"+(this.HeurTour[i]+1)+", ");


			}

			System.out.println("v"+(this.HeurTour[0]+1));*/
			//System.out.println("HeurCost:"+Calcost(this.HeurTour, this.dim));


			return this;

		}

////////////////////////////////////////////////// End of Heuristic Greedy algorithm ////////////////////////////////

//////////////////////////////////////////////// Beginning of Hill  Climbing /////////////////////////////////////////////////

/*	public TSP Hill(){

		this.HillTour = new int[this.dim];
		this.HillCost = 0;



			


		int[][] removed = new int[this.dim][3];

		int[][] newDistance;

		int[] nodes;

		

		int size = 0;

		int num = 0;

		int r = 0;

		boolean q = true;


		for(int i = 0; i < this.dim; i++){
			removed[i][0] = 0;
			removed[i][1] = -1;
			removed[i][2] = -1;		
		}
		
		


		while(!(Hami(removed, this.dim))){

			num = numDeg0(removed, this.dim) + numDeg1(removed, this.dim);
			
			newDistance = new int[num][num];
			nodes = new int[num];



			newDistance = updateDistance(this.DisMatrix, nodes, num);


			nodes = newNodes(removed, this.dim);


			size = numDeg1(removed, this.dim)/2;




			if(size > 0) { 
				int[][] M = new int[size][2];
				M = matching(removed, this.dim); 
				newDistance = modifiedDistance(newDistance, M, size, (-1) * this.OptCost); 
			}

			int[][] basic = new int[num][num];

			r = num;

			basic = LocalOpt3(newDistance, num, r);

			if(continueing(basic, num, r)){
				removed = Removal(removed, nodes, basic, num, r);
			}else if(num == 3){
				nodes = newNodes(removed, this.dim);
				int[] Tour = new int[3];
				Tour = nodes;
				basic = allContaining(Tour, num, r);
				removed = Removal(removed, nodes, basic, num, r);
			}else {
				r = 2;
				basic = LocalOpt3(newDistance, num, r);
				q = continueing(basic, num, r);
				if(!q){
					int[] Tour = new int[num];
					Tour = Permutation(num);
					Tour = localOpt3(newDistance, Tour, num, 0, basic);
					basic = allContaining(Tour, num, r);
					removed = Removal(removed, nodes, basic, num, r);
				}
				removed = Removal(removed, nodes, basic, num, r);
			}
		}

		


		this.HillTour = construct(removed, this.dim);

		System.out.println("HillTour:");
			for(int i =0; i < this.dim; i++) {
				System.out.print("v"+(this.HillTour[i]+1)+", ");
			}


			System.out.println("v"+(this.HillTour[0]+1));

		this.HillCost = Calcost(this.HillTour, this.dim);

		System.out.println("HillCost:"+this.HillCost);

		return this;

	}




	public int[][] allContaining(int[] Tour, int N, int r){
		int[][] allIn = new int[N][N];

		for(int i = 0; i < N-1; i++){
			allIn[Tour[i]][Tour[i+1]] = r;
			allIn[Tour[i+1]][Tour[i]] = r;
		}

		allIn[Tour[0]][Tour[N-1]] = r;
		allIn[Tour[N-1]][Tour[0]] = r;

		return allIn;
	}


	public int[] construct(int[][] removed, int N){


		int[] Tour = new int[N];

		int times = 1;
		


		int pair, backNode;

		pair = 0;

		backNode = pair;

		pair = removed[pair][1];

		Tour[times-1] = backNode;

		times++;

		while(times <= N){
			if(backNode == removed[pair][1]){
				backNode = pair;
				pair = removed[pair][2];
			}else {
				backNode = pair;
				pair = removed[pair][1];
			}
			
			Tour[times-1] = backNode;
		

			times++;
		}


		return Tour;


	}



	public boolean Hami(int[][] removed, int N){
		boolean H = true;

		for(int i = 0; i < N; i++){
			if(removed[i][0] != 2) {H = false; break;}
		}

		return H;
	}



	public int[][] matching(int[][] removed, int N){
		int size = numDeg1(removed, N)/2;

		int[] picked = new int[N];

		int[][] matching = new int[size][2];

		int position = 0;

		int pair = 0;

		int backNode;

		for(int i = 0; i < N; i++){
			picked[i] = 0;
		}

		for(int i = 0; i < N; i++){
			if(picked[i] == 0 && removed[i][0] == 1){
				matching[position][0] = i;
				
				pair = removed[i][1];

				picked[i] = 1;

				backNode = i;

				while(removed[pair][0] == 2){
					if(backNode == removed[pair][1]){
						backNode = pair;
						pair = removed[pair][2];
					}else {
						backNode = pair;
						pair = removed[pair][1];
					}

				}





				matching[position][1] = pair;

				picked[pair] = 1;

				position++;
			}

		}


		return matching;
	}

	public int numDeg0(int[][] removed, int N){
		int deg0 = 0;

		for(int i = 0; i < N; i++){
			if(removed[i][0] == 0){
				deg0++;
			}
		}

		return deg0;
	}

	public int numDeg1(int[][] removed, int N){
		int deg1 = 0;

		for(int i = 0; i < N; i++){
			if(removed[i][0] == 1){
				deg1++;
			}
		}

		return deg1;
	}

	public int[] newNodes(int[][] removed, int N){

		int num = numDeg0(removed, N) + numDeg1(removed, N);

		int[] nodes = new int[num];

		int position = 0;

		for(int i = 0; i < N; i++){
			if(removed[i][0] != 2){
				nodes[position] = i;
				position++;
			}
		}

		return nodes;
	}

	public int[][] updateDistance(int[][] distance, int[] nodes, int num){
		int[][] newDistance = new int[num][num];

		for(int i = 0; i < num; i++){
			for(int j = 0; j <= i; j++){
				newDistance[i][j] = distance[nodes[i]][nodes[j]];
				newDistance[j][i] = distance[nodes[j]][nodes[i]];
			}
		}
		
		return newDistance;
	}




	public int[][] modifiedDistance(int[][] newDistance, int[][] matching, int size, int negInfty){

		for(int i = 0; i < size; i++){

			newDistance[matching[i][0]][matching[i][1]] = negInfty;
			newDistance[matching[i][1]][matching[i][0]] = negInfty;
		}
		
		return newDistance;
	}

	public int[][] Removal(int[][] removed, int[] nodes, int[][] basic, int N, int r){

		for(int i = 0; i < N; i++){
			for(int j = 0; j <= i; j++){

				if(basic[i][j] == r){

					removed[nodes[i]][0] = removed[nodes[i]][0] + 1;
					removed[nodes[j]][0] = removed[nodes[j]][0] + 1;
			
					removed[nodes[i]][removed[nodes[i]][0]] = nodes[j];
					removed[nodes[j]][removed[nodes[j]][0]] = nodes[i];
				}
			}
		}










		return removed;
	}

	public boolean continueing(int[][] basic, int N, int r){



		boolean continuity = false;

		for(int i = 0; i < N; i++){
			for(int j = 0; j <= i; j++){
				if(basic[i][j] == r){
					continuity = true;
				}
			}



		}


		return continuity;
	}


	public int[][] LocalOpt3(int[][] distance, int N, int r){

		int[][] basic = new int[N][N];

		int[] Tour = new int[N];

		int q;

		for(int i = 0; i < N; i++){
			for(int j = 0; j <= i; j++){
				basic[i][j] = 0;
				basic[j][i] = 0;
			}



		}

		for(int m = 1; m <= r; m++){
			if(m == 1){
				q = 0;
			}else {
				q =1;
			}






			Tour = Permutation(N);

			Tour = localOpt3(distance, Tour, N, q, basic);

			if(q == 1){
				q = 0; 
				Tour = localOpt3(distance, Tour, N, q, basic);
			}






			for(int i = 0; i < N-1; i++){
				basic[Tour[i]][Tour[i+1]] = basic[Tour[i]][Tour[i+1]] + 1;
				basic[Tour[i+1]][Tour[i]] = basic[Tour[i+1]][Tour[i]] + 1;
			}

	
			basic[Tour[0]][Tour[N-1]] = basic[Tour[0]][Tour[N-1]] + 1;
			basic[Tour[N-1]][Tour[0]] = basic[Tour[N-1]][Tour[0]] + 1;


		}


		return basic;

	}

	

	public int[] localOpt3(int[][] distance, int[] Tour, int N, int q, int[][] basic){

		boolean anti = true;

		int d1, d2, d3, d4, d;

		

		for(int count = 1; count <= N; count++){

			if(q == 0 || basic[Tour[0]][Tour[N-1]] == 0){
				for(int k = 1; k <=N-2; k++){
					for(int j = k+1; j <= N-1; j++){
	
							d1 = distance[Tour[k-1]][Tour[j]] + distance[Tour[0]][Tour[j-1]];
							d2 = distance[Tour[0]][Tour[j]] + distance[Tour[k-1]][Tour[j-1]];
							
							if(d1 <= d2){
								d = d1;
								anti = false;
							}else {
								d = d2; 
								anti = true;
							}
							
							d3 = d + distance[Tour[k]][Tour[N-1]];
							d4 = distance[Tour[0]][Tour[N-1]] + distance[Tour[k-1]][Tour[k]] + distance[Tour[j-1]][Tour[j]];
		
							if(d3 < d4){
								if(anti == false){
									Tour = Clockwise(Tour, N, k, j);
								}else {
									Tour = AntiClockwise(Tour, N, k, j);
								}
								
								count = 1;
								k = 1;
								j = k;
							}
							
					}
				}
			}
			





			Tour = Rotation(Tour, N);
			
		}

		//if(N == this.dim) System.out.println(Calcost(optTour, this.dim));
	
	
		return Tour;
	
	}
	

	

	
	public int[] Clockwise(int[] Tour, int N, int k, int j){ // k < j clockwise: (j+1, ..., N-1, k, ..., j-1, 0, ..., k-1, j)

		int[] PerturbedTour = new int[N];		

		for(int i = 0; i <= N-j-2; i++){
			PerturbedTour[i] = Tour[i+j+1];
		}

		for(int i = N-j-1; i <= N-k-2; i++){
			PerturbedTour[i] = Tour[i-N+j+k+1];
		}

		for(int i = N-k-1; i <= N-2; i++){
			PerturbedTour[i] = Tour[i-N+k+1];
		}

		PerturbedTour[N-1] = Tour[j];

		return PerturbedTour;
	}



	public int[] AntiClockwise(int[] Tour, int N, int k, int j){ // k < j   anticlockwise: (j+1, ..., N-1, k, ..., j-1, k-1, ..., 0, j)

		int[] PerturbedTour = new int[N];		

		for(int i = 0; i <= N-j-2; i++){
			PerturbedTour[i] = Tour[i+j+1];
		}

		for(int i = N-j-1; i <= N-k-2; i++){
			PerturbedTour[i] = Tour[i-N+j+k+1];
		}

		for(int i = N-k-1; i <= N-2; i++){
			PerturbedTour[i] = Tour[N-2-i];
		}

		PerturbedTour[N-1] = Tour[j];

		return PerturbedTour;
	}



	*/

///////////////////////////////////////////////End of Hill Climbing /////////////////////////////////////////////////////////


////////////////////////////////////////////// Beginning of Hill Climbing 2 //////////////////////////////////////////////


	public TSP Hill2(int randseed){

		this.HillTour2 = new int[this.dim];
		this.HillCost2 = 0;	




		TSP heur = this.Heur();

		int opt = 3;
	
		int c1, c2;

		int delta;
		
		int[] currTour = heur.HeurTour;
		int currCost = heur.HeurCost;


		try{
			BufferedWriter outTrace = new BufferedWriter(new FileWriter(tracefile,true));
			long starttime = System.currentTimeMillis();
			outTrace.write("0.00"+" "+currCost+ "\n");
			String tracefile = this.tracefile;
			
			do{

				long endtime = System.currentTimeMillis();
				double limittime = (double)(endtime-starttime) / 1000.0;
				if(limittime > this.cutoff) {break;}
			
				this.HillTour2 = Permutation1(this.dim, randseed);
			
				randseed++;
		
				do{
					if(opt == 3){
						c1 = Calcost(this.HillTour2, this.dim);
						localOpt3_2();
						c2 = Calcost(this.HillTour2, this.dim);
						delta = c2 - c1;
						opt = 2;
					}else {
						c1 = Calcost(this.HillTour2, this.dim);
						localOpt2_2();
						c2 = Calcost(this.HillTour2, this.dim);
						delta = c2 - c1;
						opt = 3;
					}
		
				} while(delta > 0);
				
				if(Calcost(currTour, this.dim) > Calcost(this.HillTour2, this.dim)){
					currTour = this.HillTour2;
					long currtime = System.currentTimeMillis();
					double timestamp = (currtime- starttime) /1000.0 ;
					outTrace.write(String.format("%.2f", timestamp) +" "+Calcost(currTour, this.dim)+ "\n");
					//System.out.println("Hill2 + Tour:"+Calcost(Tour, this.dim));
				}

	
	
	
			} while(Calcost(currTour, this.dim) > this.OptCost);

		
		this.HillTour2 = currTour;

		/*System.out.println("HillTour2:");
			for(int i =0; i < this.dim; i++) {
				System.out.print("v"+(this.HillTour2[i]+1)+", ");
			}

			System.out.println("v"+(this.HillTour2[0]+1));*/

		this.HillCost2 = Calcost(this.HillTour2, this.dim);
			
		outTrace.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}


		return this;		
		

	}



	public TSP localOpt2_2(){
		int N = this.dim;



		int d1, d2;

		

		for(int count = 1; count <= N; count++){

			for(int k = 2; k <=N-2; k++){


					int a = this.HillTour2[0];
					int b = this.HillTour2[k-1];
					int c = this.HillTour2[k];
					int d = this.HillTour2[N-1];

					d1 = this.DisMatrix[a][d] + this.DisMatrix[b][c];
					d2 = this.DisMatrix[a][c] + this.DisMatrix[b][d];	
						
						
	
						if(d2 < d1){
							//System.out.println("Effect of Change Before:"+Calcost(this.HillTour2, N));
							this.HillTour2= Change(this.HillTour2, N, k);
							//System.out.println("Effect of Change After:"+Calcost(this.HillTour2, N));

							count = 1;
							k = 1;
						}		
						
				
			}
			
			this.HillTour2 = Rotation(this.HillTour2, N);
		}


	

		return this;
	
	}

	

	public TSP localOpt3_2(){
		int N = this.dim;

		int changeType;

		int d1, d2, d3, d4, d0;

		

		for(int count = 1; count <= N; count++){
			for(int k = 2; k <=N-4; k++){
				for(int j = k+2; j <= N-2; j++){

					int a = this.HillTour2[0];
					int b = this.HillTour2[k-1];
					int c = this.HillTour2[k];
					int d = this.HillTour2[j-1];
					int e = this.HillTour2[j];
					int f = this.HillTour2[N-1];

					d0 = this.DisMatrix[a][f] + this.DisMatrix[b][c] + this.DisMatrix[d][e];
					d1 = this.DisMatrix[a][c] + this.DisMatrix[d][f] + this.DisMatrix[e][b];
					d2 = this.DisMatrix[a][d] + this.DisMatrix[c][e] + this.DisMatrix[f][b];
					d3 = this.DisMatrix[a][d] + this.DisMatrix[c][f] + this.DisMatrix[e][b];
					d4 = this.DisMatrix[b][d] + this.DisMatrix[c][f] + this.DisMatrix[e][a];

					changeType = 0;
						
						if(d1 < d0){
							d0 = d1;
							changeType = 1;
						}else if(d2 < d0){
							d0 = d2; 
							changeType = 2;
						}else if(d3 < d0){
							d0 = d3; 
							changeType = 3;
						}else if(d4 < d0){
							d0 = d4; 
							changeType = 4;
						}
						
						
	
						if(changeType != 0){
							if(changeType == 1){
								//System.out.println("Effect of Change1 Before:"+Calcost(this.HillTour2, N));
								this.HillTour2= Change1(this.HillTour2, N, k, j);
								//System.out.println("Effect of Change1 After:"+Calcost(this.HillTour2, N));
							}else if(changeType == 2){
								//System.out.println("Effect of Change2 Before:"+Calcost(this.HillTour2, N));
								this.HillTour2= Change2(this.HillTour2, N, k, j);
								//System.out.println("Effect of Change2 After:"+Calcost(this.HillTour2, N));
							}else if(changeType == 3){
								//System.out.println("Effect of Change3 Before:"+Calcost(this.HillTour2, N));
								this.HillTour2= Change3(this.HillTour2, N, k, j);
								//System.out.println("Effect of Change3 After:"+Calcost(this.HillTour2, N));
							}else if(changeType == 4){
								//System.out.println("Effect of Change4 Before:"+Calcost(this.HillTour2, N));
								this.HillTour2= Change4(this.HillTour2, N, k, j);
								//System.out.println("Effect of Change4 After:"+Calcost(this.HillTour2, N));
							}

							count = 1;
							k = 2;
							j = k+1;

						}		
						
						
				}
			}
			
			this.HillTour2 = Rotation(this.HillTour2, N);
		}




	return this;
	
	}

	
	public int[] Change(int[] Tour, int N, int k){ // k, 0 to k to N-1 to k-1, edges: (0, k), (N-1, k-1)

		int[] PerturbedTour = new int[N];	

		int i1, i2;

		i1 = 0;
		i2 = 0;

		while(i1 <= k-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2++;
		}	

		i2 = N-1;

		while(i1 <= N-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		return PerturbedTour;
	}


	
	public int[] Change1(int[] Tour, int N, int k, int j){ // k < j, 0 to k to j-1 to N-1 to j to k-1, edges: (0, k), (j-1, N-1), (j, k-1)

		int[] PerturbedTour = new int[N];	

		int i1, i2;

		i1 = 0;
		i2 = k-1;

		while(i1 <= k-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}	

		i2 = k;

		while(i1 <= j-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2++;
		}

		i2 = N-1;

		while(i1 <= N-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		return PerturbedTour;
	}



	public int[] Change2(int[] Tour, int N, int k, int j){ // k < j,  0 to j-1 to k to j to N-1 to k-1, edges: (0, j-1), (k, j), (N-1, k-1)

		int[] PerturbedTour = new int[N];	

		int i1, i2;

		i1 = 0;
		i2 = k-1;

		while(i1 <= k-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}	

		i2 = j-1;

		while(i1 <= j-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		i2 = j;

		while(i1 <= N-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2++;
		}

		return PerturbedTour;

	}



	public int[] Change3(int[] Tour, int N, int k, int j){ // k < j, 0 to j-1 to k to N-1 to j to k-1, edges: (0, j-1), (k, N-1), (j, k-1)

		int[] PerturbedTour = new int[N];	

		int i1, i2;

		i1 = 0;
		i2 = k-1;

		while(i1 <= k-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}	

		i2 = j-1;

		while(i1 <= j-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		i2 = N-1;

		while(i1 <= N-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		return PerturbedTour;
	}



	public int[] Change4(int[] Tour, int N, int k, int j){ // k < j, k-1 to j-1 to k to N-1 to j to 0, edges: (k-1, j-1), (k, N-1), (j, 0)

		int[] PerturbedTour = new int[N];	

		int i1, i2;

		i1 = 0;
		i2 = 0;

		while(i1 <= k-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2++;
		}	

		i2 = j-1;

		while(i1 <= j-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		i2 = N-1;

		while(i1 <= N-1) {
			PerturbedTour[i1] = Tour[i2];
			i1++;
			i2--;
		}

		return PerturbedTour;
	}

	public int[] Rotation(int[] Tour, int N){


		int[] RotatedTour = new int[N];		

		for(int i = 0; i < N; i++){

			if(i == 0){
				RotatedTour[i] = Tour[N-1];

			}else {
				RotatedTour[i] = Tour[i-1];

			}
		}


		return RotatedTour;

	}


	public int[] Permutation(int N) { 

      		int[] a = new int[N];

			//Random rand = new Random(randseed);
			
      		// insert integers 0...N-1

      		for (int i = 0; i < N; i++) {a[i] = i;}

      		// shuffle
      		for (int i = 0; i < N; i++) {
         		int r = (int) (Math.random() * (i+1));     // int between 0 and i
         		int swap = a[r];
         		a[r] = a[i];
         		a[i] = swap;

      		}

   		//print permutation
   		//for (int i = 0; i < N; i++) {System.out.print(a[i] + " ");}
                //System.out.println(""); 

		//if(N == this.dim) System.out.println(Calcost(a, this.dim));

		return a;

	}

//////////////////////////////////////////// End of Hill Climbing 2 ///////////////////////////////////////////////////////


//////////////////////////////////////Beginning of Local Search SA////////////////////////////////////////////
	public static double acceptanceProbability(int currE, int newE, double temperature){
		if (newE < currE)
			return 1.0;
		else
			return Math.exp((currE-newE) / temperature);
	}
	public TSP SA(int randomseed){
		TSP heur = this.Heur();
		Random rand = new Random(randomseed);
		int N = this.dim;
		int delta = this.OptCost;
		double temp = 15000;
		double coolingRate = 0.0001;
		//int[] currTour = Permutation1(N, randomseed);
		//int currCost = Calcost(currTour,N);
		int[] currTour = heur.HeurTour;
		int currCost = heur.HeurCost;
		this.SATour = currTour;
		this.SACost = currCost;
		try{
			BufferedWriter outTrace = new BufferedWriter(new FileWriter(tracefile,true));//Append mode
			//System.out.println(this.SACost);
			long starttime = System.currentTimeMillis();
			outTrace.write("0.00"+" "+this.SACost+ "\n");
			String tracefile = this.tracefile;
			//int loop = 0;
			do{
				//loop++;
				long endtime = System.currentTimeMillis();
				double limittime = (double)(endtime-starttime) / 1000.0;
				if(limittime > this.cutoff)
					break;
				int t1 = rand.nextInt(N);
				int t2 = rand.nextInt(N);
				int tourPos1;
				int tourPos3;
				while(t1 == t2 || Math.abs(t1-t2) == 1 || Math.abs(t1-t2) == N-1){
					t1 =  rand.nextInt(N);
					t2 =  rand.nextInt(N);
	
	
				}
				if(t1 < t2){
					tourPos1 = t1;
					tourPos3 = t2;
				}
				else{
					tourPos1 = t2;
					tourPos3 = t1;
				}
				int tourPos2 = tourPos1 + 1;
				int tourPos4 = tourPos3 + 1;
				if(tourPos1== N-1)
				{
					tourPos2 = 0;
				}
				if(tourPos3 == N-1)
				{
					tourPos4 = 0;
				}
				//System.out.println("1:"+tourPos1+"  2:"+tourPos2+"  3:"+tourPos3+"  4:"+tourPos4);
				int newCost = currCost-this.DisMatrix[currTour[tourPos1]][currTour[tourPos2]]-this.DisMatrix[currTour[tourPos3]][currTour[tourPos4]]+this.DisMatrix[currTour[tourPos1]][currTour[tourPos3]]+this.DisMatrix[currTour[tourPos2]][currTour[tourPos4]];
				
				//System.out.println("NewCost:"+newCost);
				double pro = rand.nextDouble();
				//System.out.println("pro:"+pro);
				if(this.SACost > newCost)
				{
					
					currCost = newCost;
					currTour = Exchange(currTour, tourPos1, tourPos2, tourPos3, tourPos4);
					this.SACost = currCost;
					this.SATour = currTour;
					long currtime = System.currentTimeMillis();
					double timestamp = (currtime- starttime) /1000.0 ;
					outTrace.write(String.format("%.2f", timestamp) +" "+this.SACost+ "\n");
					//System.out.println("CurrOptCost:"+this.SACost+"  TIME:"+timestamp);
				}
				else if(acceptanceProbability(currCost,newCost,temp) > pro){
						currCost = newCost;
						currTour = Exchange(currTour, tourPos1, tourPos2, tourPos3, tourPos4);				
				}
				temp = (1-coolingRate) * temp;
				delta = this.SACost - this.OptCost;
			}while(delta > 0 && temp > 0.001);
		//System.out.println("loop:"+loop);
			
		outTrace.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}
		//-----------------------This four line is for printing result
		//System.out.println("SAcost:"+this.SACost);
		//double per = (double) (this.SACost - this.OptCost) / (double)this.OptCost;
		//System.out.println("Rate:"+ per);
		//System.out.println("-----------------------------------------");
		return this;			
	}
	public int[] Exchange(int[] currTour, int pos1, int pos2, int pos3, int pos4){
		int N = this.dim;
		int[] newTour = new int[N];
		for(int i = 0; i <= pos1; i++){
			newTour[i] = currTour[i];
		}
		newTour[pos1+1] = currTour[pos3];
		int j = pos1+2;
		for(int i = pos3-1; i >= pos2 ; i--){
			newTour[j] = currTour[i];
			j++;
			//System.out.println("i:"+i+"j:"+j);
		}
		
		//System.out.println("j:"+j);
		if(pos4 != 0){
			newTour[j] = currTour[pos4];
			for(int i = pos4+1; i < N; i++){
				newTour[i] = currTour[i];
			}
		}
		return newTour;
	}	
	
//	public int[] Permutation(int N) { 
//
//  		int[] a = new int[N];
//
//  		// insert integers 0...N-1
//
//  		for (int i = 0; i < N; i++) {a[i] = i;}
//
//  		// shuffle
//  		for (int i = 0; i < N; i++) {
//     		int r = (int) (Math.random() * (i+1));     // int between 0 and i
//     		int swap = a[r];
//     		a[r] = a[i];
//     		a[i] = swap;
//  		}
//
//		//print permutation
//		//for (int i = 0; i < N; i++) {System.out.print(a[i] + " ");}
//            //System.out.println(""); 
//
//	//if(N == this.dim) System.out.println(Calcost(a, this.dim));
//
//  		return a;
//	}
//	public int Calcost(int[] Tour, int N){
//		int cost = DisMatrix[Tour[0]][Tour[N-1]];
//
//
//		for(int i = 0; i < N-1; i++){
//			cost = cost + DisMatrix[Tour[i]][Tour[i+1]];
//		}
//
//		return cost;
//	}

public int[] Permutation1(int N, int randseed) { 
		int[] a = new int[N];
  		Random rand = new Random(randseed);
  		// insert integers 0...N-1


  		for (int i = 0; i < N; i++) {a[i] = i;}


  		// shuffle
  		for (int i = 0; i < N; i++) {
     		int r = (int) (rand.nextDouble() * (i+1));     // int between 0 and i
     		int swap = a[r];
     		a[r] = a[i];
     		a[i] = swap;
  		}

		//print permutation
		//for (int i = 0; i < N; i++) {System.out.print(a[i] + " ");}
            //System.out.println(""); 

	//if(N == this.dim) System.out.println(Calcost(a, this.dim));
	return a;


}


/////////////////////////////////////End of Local Search SA////////////////////////////////////////////


////////////////////////////////////////////// Beginning of GetDistance///////////////////////////////////////////////




	public int GetDistance(double lat1, double lng1, double lat2, double lng2){
		double Earth_Radius =  6378.388;
		double radLat1 = toRadians(lat1);
		double radLat2 = toRadians(lat2);
		double radLng1 = toRadians(lng1);
		double radLng2 = toRadians(lng2);
		
		double latDiff = radLat1 - radLat2;//toRadians(lat1-lat2);
		double lngDiff = radLng1 - radLng2;//toRadians(lng1-lng2);
		double latSum = radLat1 + radLat2;//toRadians(lat1+lat2);
		//double s = Math.sin(radLat1)*Math.sin(radLat2)*Math.cos(b) + Math.cos(radLat1)*Math.cos(radLat2);
		 double q1,q2,q3;
		 q1 = (double) Math.cos(lngDiff);//( radLng1 - radLng2 ); 
		 q2 = (double) Math.cos(latDiff);//( radLat1 - radLat2 ); 
		 q3 = (double) Math.cos(latSum);//( radLat1 + radLat2 ); 
		 int s = (int) ( Earth_Radius * Math.acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0);
		
//		double Earth_Radius = 6378.388;
//		double radLat1 = Math.toRadians(90 - lat1);
//		double radLat2 = Math.toRadians(90 - lat2);
//		double radLng1 = Math.toRadians(lng1);
//		double radLng2 = Math.toRadians(lng2);
//		double a = radLat1 -radLat2;
//		double b = radLng1 - radLng2;
//
//		double s = Math.sin(radLat1)*Math.sin(radLat2)*Math.cos(b) + Math.cos(radLat1)*Math.cos(radLat2);
//
//		s = Math.acos(s);
//		s = s * Earth_Radius;		
		return s;


	}


	public double toRadians(double inp){
		double PI =  3.141592;
		
		int deg = (int) (inp);//(int)
		double min = inp - deg;
		
		double rad =  (PI * (deg + 5.0 * min/3.0) / 180.0); 
		
		return rad;
	}
		


	public int Calcost(int[] Tour, int N){
		int cost = DisMatrix[Tour[0]][Tour[N-1]];


		for(int i = 0; i < N-1; i++){
			cost = cost + DisMatrix[Tour[i]][Tour[i+1]];
		}

		return cost;
	}



///////////////////////////////////////////// End of GetDistance ////////////////////////////////////////////////////////

	
}
		
