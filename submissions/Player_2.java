// updated for 9/16 deliverable - adding the wall
package flip.g4;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Comparator;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.HashMap;
import javafx.util.Pair; 
import java.util.ArrayList;
import java.util.Arrays;
import java.lang.*;

import flip.sim.Point;
import flip.sim.Board;
import flip.sim.Log;
//import flip.g4.HungarianAlgorithm;

public class Player implements flip.sim.Player
{
	private int seed = 42;
	private Random random;
	private boolean isplayer1;
	private Integer n;
	private Double diameter_piece;

    // wall
    private HashMap<Integer, Point> my_pieces;
    private Point[]   ideal_wall_pieces;
    private Integer[] ideal_piece_match;
    private int count;
    private int numWallBuilds = 1;

	public Player()
	{
		random = new Random(seed);
	}

	// Initialization function.
    // pieces: Location of the pieces for the player.
    // n: Number of pieces available. (default 30 in Makefile)
    // t: Total turns available.
	public void init(HashMap<Integer, Point> pieces, int n, double t, boolean isplayer1, double diameter_piece)
	{
		this.n = n;
		this.isplayer1 = isplayer1;
		this.diameter_piece = diameter_piece; // default 2
		
		// wall
		this.count = 0;
		this.my_pieces = pieces;
		this.ideal_piece_match = this.calculateWallStrategy(pieces, isplayer1);
	}

    // ******************
    // wall functions
    private Integer[] calculateWallStrategy(HashMap<Integer, Point> playerPieces, boolean isPlayer1){
        // Generate co-ordinates for ideal wall
        // TODO: Deviate from "ideal" wall to "approx" wall to reduce moves
        int num_wall_pieces = 11;
        this.ideal_wall_pieces = new Point[num_wall_pieces];

        Log.log(String.valueOf(isPlayer1) + playerPieces.toString());

        for(int i=0; i<num_wall_pieces; i++){
            // Calculations done assuming perfect placement for 1cm coins on 40cm board
            ideal_wall_pieces[i] = new Point((isPlayer1)? 21.5: -21.5, 2*i*Math.sqrt(3) + Math.sqrt(3) - 19);
        }

        // Figure out the distance from each coin to the ideal wall piece location
        double[][] cost_matrix = new double[this.n][num_wall_pieces];
        for(int i=0; i<num_wall_pieces; i++){
            for(int j=0; j<this.n; j++){
                double dist = Math.hypot(
                    ideal_wall_pieces[i].x - playerPieces.get(j).x,
                    ideal_wall_pieces[i].y - playerPieces.get(j).y
                );
                
                cost_matrix[j][i] = Board.almostEqual(dist, 2) ? 1: Math.max(2, Math.ceil(dist/2));
            }
        }

        // Solve hungarian algorithm for minimum placement in O(n^3)
        HungarianAlgorithm solver = new HungarianAlgorithm(cost_matrix);
        int[] solution = solver.execute();
        Log.log("CURRENT : Hungarian Algorithm solved " + Arrays.toString(solution));

        // Store calculated solution
        Integer[] ideal_piece_match = new Integer[11];
        for(int pieceID=0; pieceID<this.n; pieceID++){
            if(solution[pieceID]>=0) ideal_piece_match[solution[pieceID]] = pieceID;
        }

        Log.log("CURRENT : Ideal piece match. " + Arrays.toString(ideal_piece_match));
        return ideal_piece_match;
    }
    
    private Integer[] getWallPriority(HashMap<Integer, Point> opponentPieces){
        Integer[] distanceToWall = new Integer[11];
        for(int i=0; i<11; i++)
            distanceToWall[i] = Integer.MAX_VALUE;

        // Iterate through opponent pieces to prioritize wall
        for(int j=0; j<this.n; j++){
            for (int i=0; i<11; i++){
                double dist = Math.hypot(
                    opponentPieces.get(j).x-ideal_wall_pieces[i].x,
                    opponentPieces.get(j).y-ideal_wall_pieces[i].y);
                
                Integer numMoves = (Board.almostEqual(dist, 2)) ? 1: (int) Math.max(2, Math.ceil(dist/2));
                if (numMoves < distanceToWall[i])
                    distanceToWall[i] = numMoves;
            }
        }

        // argsort() function
        // Copied from https://stackoverflow.com/questions/4859261/get-the-indices-of-an-array-after-sorting
        // Will re-implement later
        class ArrayIndexComparator implements Comparator<Integer> {
            private final Integer[] array;
            public ArrayIndexComparator(Integer[] array) { this.array = array;}
            public Integer[] createIndexArray() {
                Integer[] indexes = new Integer[array.length];
                for (int i = 0; i < array.length; i++)
                    indexes[i] = i; // Autoboxing
                return indexes;
            }
            @Override
            public int compare(Integer index1, Integer index2){
                // Autounbox from Integer to int to use as array indexes
                return array[index1].compareTo(array[index2]);
            }
        }

        ArrayIndexComparator comparator = new ArrayIndexComparator(distanceToWall);
        Integer[] indexes = comparator.createIndexArray();
        Arrays.sort(indexes, comparator);

        return indexes;
    }

    public List<Pair<Integer, Point>> getWallMove( Integer numMoves, HashMap<Integer, Point> playerPieces, HashMap<Integer, Point> opponentPieces, boolean isPlayer1) {

        List<Pair<Integer, Point>> moves = new ArrayList<Pair<Integer, Point>>();
        Integer[] bestWallPieces = this.getWallPriority(opponentPieces);
        this.count++;

        for(int i=0; i<11 && moves.size()<2; i++){
            // If piece is in place
            Point wallPiece = ideal_wall_pieces[bestWallPieces[i]];
            int  idealPiece = ideal_piece_match[bestWallPieces[i]];
            Point myPiece   = playerPieces.get(idealPiece);

            double dist = Math.hypot(wallPiece.x-myPiece.x, wallPiece.y-myPiece.y);
            if (dist <= 1e-2) // Check if the piece is in place
                continue;
            else if(Board.almostEqual(dist, 2)) {
                moves.add(new Pair<Integer,Point>(idealPiece, wallPiece));
            }
            else if (dist < 4) {
                double x1 = 0.5 * (wallPiece.x + myPiece.x);
                double y1 = 0.5 * (wallPiece.y + myPiece.y);

                double sqrt_const = Math.sqrt(16/(dist*dist)-1) / 2;
                double x2 = sqrt_const * (wallPiece.y - myPiece.y);
                double y2 = sqrt_const * (myPiece.x - wallPiece.x);

                Pair<Integer, Point> move = new Pair<Integer, Point>( idealPiece, new Point(x1+x2, y1+y2));

                if(check_validity(move, playerPieces, opponentPieces)){
                    moves.add(move);

                    Pair<Integer, Point> move2 = new Pair<Integer, Point>( idealPiece, new Point(wallPiece.x, wallPiece.y));
                    if(check_validity(move2, playerPieces, opponentPieces))
                        moves.add(move2);
                } else{
                    move = new Pair<Integer, Point>( idealPiece, new Point(x1-x2, y1-y2));
                    if(check_validity(move, playerPieces, opponentPieces)){
                        moves.add(move);

                        Pair<Integer, Point> move2 = new Pair<Integer, Point>( idealPiece, new Point(wallPiece.x, wallPiece.y));
                        if(check_validity(move2, playerPieces, opponentPieces))
                            moves.add(move2);
                    }
                }
            } else {
                Pair<Integer, Point> move = new Pair<Integer, Point>(idealPiece, new Point(
                    myPiece.x + 2 * (wallPiece.x - myPiece.x) / dist,
                    myPiece.y + 2 * (wallPiece.y - myPiece.y) / dist
                ));

                if(check_validity(move, playerPieces, opponentPieces)){
                    moves.add(move);
                
                    if (dist > 6 && moves.size()==1){
                        Pair<Integer, Point> move2 = new Pair<Integer, Point>(idealPiece, new Point(
                            myPiece.x + 4 * (wallPiece.x - myPiece.x) / dist,
                            myPiece.y + 4 * (wallPiece.y - myPiece.y) / dist
                        ));
                        
                        if(check_validity(move2, playerPieces, opponentPieces))
                            moves.add(move2);
                    }
                }
            }
        }

        if (moves.size()==0 && this.count < 1000){
            Log.log("WALL COMPLETED IN " + String.valueOf(this.count) + " MOVES");
            this.count = 1000;
        }
        return moves;
    }
    // ******************

    // x coordinate
	public boolean inSoftEndZone( Point piece, boolean isPlayer1 ) {
	    return (isPlayer1? -1 : 1) * piece.x > 20 + 1.75 * diameter_piece + (n / 9) * (diameter_piece / 2);
	}

    // used for getting only pieces for building a wall
    public boolean inSoftNeutralZone( Point piece, boolean isPlayer1 ) {
	    //return (isPlayer1? -1 : 1) * piece.x > 20 + 1.75 * diameter_piece + (n / 9) * (diameter_piece / 2);
	    return (isPlayer1? 1 : -1) * piece.x > 19 * diameter_piece;
	}

    // Returns a list of ALL our player pieces (or opponent's pieces) in order of closest to farthest along x
    private static List<Pair<Integer,Point>> rankXProgress( HashMap<Integer, Point> pieces, boolean arePiecesPlayer1 ) {
    	int sign = arePiecesPlayer1? 1 : -1;
    	List<Integer> indices = new ArrayList<>(pieces.keySet());
    	indices.sort(new Comparator<Integer>() {
            public int compare( Integer i1, Integer i2 ) {
                return (pieces.get(i1).x > pieces.get(i2).x)? sign : -sign;
            }
        });

        // debugging printer: prints (id, x)
        /*
        System.out.print("Pairs:  ");
        for (Integer i : indices) System.out.print(" (" + i + ", " + pieces.get(i).x + ")");
        System.out.println();
        */
        List<Pair<Integer,Point>> orderedPieces = new ArrayList<Pair<Integer,Point>>();
        for (Integer i : indices) {
            orderedPieces.add(new Pair<Integer, Point>(i, pieces.get(i)));
        }
        
        return orderedPieces;
    }
    
    
    // returns our (or our opponent's) sorted pieces list WITHOUT pieces that have already crossed the endzone
    private List<Pair<Integer,Point>> rankXProgressOutsideEndzone( HashMap<Integer,Point> pieces, boolean arePiecesPlayer1 ) {
    	List<Pair<Integer,Point>> orderedPieces = rankXProgress(pieces, arePiecesPlayer1);
        List<Pair<Integer,Point>> outsideEndzonePieces = new ArrayList<Pair<Integer,Point>>();

    	int sign = arePiecesPlayer1? -1 : 1;
    	// use indices.toArray() as a deep copy to avoid changing indices inside for loop
    	for (Pair<Integer,Point> pair : orderedPieces) {
    	    if (!inSoftEndZone(pieces.get(pair.getKey()), arePiecesPlayer1)) {
//    	    if ((pieces.get(pair.getKey()).x *sign) < (20 + diameter_piece/2)) {
    	        outsideEndzonePieces.add(pair);
    	    }
    	}
    	
    	/*
    	// debugging printer
    	System.out.print("Pairs:  ");
        for (Pair<Integer,Point> pair : outsideEndzonePieces) System.out.print(" (" + pair.getKey() + ", " + pair.getValue().x + ")");
        System.out.println();
        */
        
        return outsideEndzonePieces;
    }
    
    // returns our (or our opponent's) sorted pieces list WITHOUT pieces that have already crossed the endzone
    private List<Pair<Integer,Point>> rankXProgressOutsideSoftNeutral( HashMap<Integer,Point> pieces, boolean arePiecesPlayer1 ) {
    	List<Pair<Integer,Point>> orderedPieces = rankXProgress(pieces, arePiecesPlayer1);
        List<Pair<Integer,Point>> outsideNeutralSoftPieces = new ArrayList<Pair<Integer,Point>>();

    	int sign = arePiecesPlayer1? -1 : 1;
    	// use indices.toArray() as a deep copy to avoid changing indices inside for loop
    	for (Pair<Integer,Point> pair : orderedPieces) {
    	    if (!inSoftNeutralZone(pieces.get(pair.getKey()), arePiecesPlayer1)) {
//    	    if ((pieces.get(pair.getKey()).x *sign) < (20 + diameter_piece/2)) {
    	        outsideNeutralSoftPieces.add(pair);
    	    }
    	}
    	
    	/*
    	// debugging printer
    	System.out.print("Pairs:  ");
        for (Pair<Integer,Point> pair : outsideEndzonePieces) System.out.print(" (" + pair.getKey() + ", " + pair.getValue().x + ")");
        System.out.println();
        */
        
        return outsideNeutralSoftPieces;
    }
    

    
    // Unused
    // Return a list of our pieces (in order of closeness to endzone) without us or an opponent's piece blocking straight path (for one move)
    // used to create list of our unblocked pieces in order of closeness to endzone
    private List<Pair<Integer,Point>> findOurUnblockedPieces( HashMap<Integer,Point> ourPieces, HashMap<Integer, Point> opponentPieces, boolean areWePlayer1 ) {
    	List<Pair<Integer,Point>> ourSortedPieces = rankXProgress(ourPieces, areWePlayer1);
    	List<Pair<Integer,Point>> opponentsSortedPieces = rankXProgress(opponentPieces, !areWePlayer1);

        // make a copy of pieces arraylist
        List<Pair<Integer,Point>> freePiecesIgnoringOpp = new ArrayList<Pair<Integer,Point>>();
        for (Pair<Integer,Point> piece : ourSortedPieces) freePiecesIgnoringOpp.add(piece);

        // remove pieces that are blocked by you
    	for (int i = 0; i < ourSortedPieces.size(); i++) {
    	    Point piece_i = ourSortedPieces.get(i).getValue();

    	    for (int j = i +1; j < ourSortedPieces.size(); j++) {
        	    Point piece_j = ourSortedPieces.get(j).getValue();
    	        
    	        // if your y path is blocked AND x path is blocked by yourself
                if (Math.abs(piece_i.y - piece_j.y) < diameter_piece && Math.abs(piece_i.x - piece_j.x) < diameter_piece) {
                    freePiecesIgnoringOpp.remove(ourSortedPieces.get(j));
                }
    	    }
    	}
    	
    	// make a copy of smaller pieces arraylist
        List<Pair<Integer,Point>> freePieces = new ArrayList<Pair<Integer,Point>>();
        for (Pair<Integer,Point> piece : freePiecesIgnoringOpp) freePieces.add(piece);
    	
    	// remove pieces that are blocked by opponent
	    for (int i = 0; i < freePiecesIgnoringOpp.size(); i++) {
            //*
            Point piece_i = freePiecesIgnoringOpp.get(i).getValue();

    	    for (int j = i +1; j < opponentsSortedPieces.size(); j++) {
        	    Point piece_j = opponentsSortedPieces.get(j).getValue();
    	        
    	        // if your y path is blocked AND x path is blocked by opponent
                if (Math.abs(piece_i.y - piece_j.y) < diameter_piece && Math.abs(piece_i.x - piece_j.x) < diameter_piece) {
                    freePieces.remove(opponentsSortedPieces.get(j));
                }
    	    }
    	    //*
	    }
    	
        // Testing
    	//System.out.println("freePieces: " + freePieces);
    	return freePieces;
    }
    

    
    // Unused
    // Return a list of opponent's pieces without an opponent's piece in their "lane" (sharing same y)
    private List<Pair<Integer,Point>> findUnblockedOpponents( HashMap<Integer, Point> opponentPieces, boolean areWePlayer1 ) {
    	boolean arePiecesPlayer1 = areWePlayer1;
    	List<Pair<Integer,Point>> sortedOpponentPieces = rankXProgress(opponentPieces, !arePiecesPlayer1);

        // make a copy of pieces arraylist
        List<Pair<Integer,Point>> freePieces = new ArrayList<Pair<Integer,Point>>();
        for (Pair<Integer,Point> piece : sortedOpponentPieces) freePieces.add(piece);

        // remove pieces that are blocked
    	for (int i = 0; i < sortedOpponentPieces.size(); i++) {
    	    Point piece_i = sortedOpponentPieces.get(i).getValue();

    	    for (int j = i +1; j < sortedOpponentPieces.size(); j++) {
        	    Point piece_j = sortedOpponentPieces.get(j).getValue();
    	        
                if (Math.abs(piece_i.y - piece_j.y) < diameter_piece) {
                    freePieces.remove(sortedOpponentPieces.get(j));
                }
    	    }
    	}
        // Testing
    	System.out.println("freePieces: " + freePieces);
    	return freePieces;
    }
    
    // Unused
    // Returns list of 11 ideal pieces for building wall
    private List<Pair<Integer,Point>> findIdealWallPieces( HashMap<Integer,Point> playerPieces, boolean isplayer1 ) {
    	// Find 11 pieces where x is the most similar (Find location with least absolute offset)
    	// Then, count number of steps to build that ideal wall	    
    	List<Pair<Integer,Point>> orderedXProgress = rankXProgress(playerPieces, isplayer1);
    	int idealNumPiecesForWall = 11;
    	
    	List<Pair<Integer,Point>> idealWallPieces = new ArrayList<Pair<Integer, Point>>();
    	int idealWallPiecesStartingIndex = 0;
    	
    	// by default, set absOffset arbitrarily high
    	Double absOffset = 1000.;
    	for (int i = 0; i < playerPieces.size() -(idealNumPiecesForWall -1); i++) {
    	    Double newAbsOffset = 0.;
    	    for (int j = 0; j < idealNumPiecesForWall -1; j++) {
    	        Point piece_1 = orderedXProgress.get(i+j).getValue();
    	        Point piece_2 = orderedXProgress.get(i+j +1).getValue();
    	        newAbsOffset += Math.abs( piece_1.x - piece_2.x );
    	    }
    	    if (newAbsOffset < absOffset) {
    	        absOffset = newAbsOffset;
    	        idealWallPiecesStartingIndex = i;
    	    }
    	}
        
        for (int i = idealWallPiecesStartingIndex; i < idealNumPiecesForWall; i++) {
            idealWallPieces.add(orderedXProgress.get(i));
        }
        // Testing
    	//System.out.println("idealWallPieces: " + idealWallPieces);
    	return idealWallPieces;
    }
    
    // Unfinished function in progress
    // Returns the minimum number of moves needed to form a viable wall ignoring opponent pieces
    	// Note: Can be used for either yourself or opponent!
    	//Note: Use and assume pieceTowardsPoint works
    /*
    private int stepsToWall( HashMap<Pair<Integer,Point>> playerPieces, boolean isplayer1 ) {
        // in progress: would use findIdealWallPieces()
        int stepsToWall = 0;
        return stepsToWall;
    }
    */

    private Pair<Integer, Point> getForwardMove( HashMap<Integer, Point> playerPieces, HashMap<Integer, Point> opponentPieces, boolean isPlayer1, boolean isFront ) {
        // get all our pieces from closest to farthest
//        List<Pair<Integer,Point>> orderedXProgress = isFront? rankXProgressOutsideEndzone(playerPieces, isPlayer1):
//                                                                            rankXProgress(playerPieces, isPlayer1);
        List<Pair<Integer,Point>> orderedXProgress = rankXProgressOutsideEndzone(playerPieces, isPlayer1);

        // one by one, check validity if we move it forward
        int max_index = orderedXProgress.size() -1;
        for (int i = 0; i <= max_index; i++) {
            Pair<Integer,Point> pair = orderedXProgress.get(isFront? i : (max_index -i));
            Point oldPosition = pair.getValue();
            double dx = (isPlayer1? -1 : 1) * diameter_piece;
            Point newPosition = new Point(oldPosition.x + dx, oldPosition.y);
            Pair<Integer,Point> move = new Pair<Integer,Point>(pair.getKey(), newPosition);

            if (check_validity(move, playerPieces, opponentPieces)) return move;
        }
        return null;        
    }
    
    private Pair<Integer, Point> getRandomMove( HashMap<Integer, Point> playerPieces, HashMap<Integer, Point> opponentPieces, boolean isPlayer1, double spread ) {
        int MAX_TRIALS = 100;

        // get list of all my piece ids
    	List<Integer> ids = new ArrayList<>(playerPieces.keySet());

        for (int trial_num = 0; trial_num < MAX_TRIALS; trial_num++) {
            // select random piece
            int id = ids.get(random.nextInt(ids.size()));
		 	Point oldPosition = playerPieces.get(id);

            // if already in soft endzone, don't use it
            if (inSoftEndZone(oldPosition, isPlayer1)) continue;

		 	// select random angle
		 	double theta = (random.nextDouble() -0.5) *spread *Math.PI/180;
		 	double dx = Math.cos(theta) * (isPlayer1? -1 : 1) * diameter_piece, dy = Math.sin(theta) * diameter_piece;
		 	Point newPosition = new Point(oldPosition.x + dx, oldPosition.y + dy);
            
            Pair<Integer,Point> move = new Pair<Integer,Point>(id, newPosition);
            if (check_validity(move, playerPieces, opponentPieces)) return move;
        }
        return null;  
    }
    
	public List<Pair<Integer, Point>> getMoves( Integer numMoves, HashMap<Integer, Point> playerPieces, HashMap<Integer, Point> opponentPieces, boolean isPlayer1 ) {
		 List<Pair<Integer, Point>> moves = new ArrayList<Pair<Integer, Point>>();
         
         
         List<Pair<Integer, Point>> wallMoves = new ArrayList<Pair<Integer, Point>>();
         wallMoves = getWallMove(numMoves, playerPieces, opponentPieces, isPlayer1);
         for (int j = 0; j < wallMoves.size(); j++) {
              moves.add(wallMoves.get(j));
         }
         
         while (moves.size() < numMoves) {
             Pair<Integer, Point> move = null;


             
             for (int i = 0; i < 20; i++) {
                  move = getForwardMove(playerPieces, opponentPieces, isPlayer1, true);
                   if (move != null) moves.add(move);
              }
                  //wallCount += 1;
             //}
             //move = piece.getValue();
             
             /*
             Pair<Integer,Point> pair = orderedXProgress.get(isFront? i : (max_index -i));
             Point oldPosition = pair.getValue();
             double dx = (isPlayer1? -1 : 1) * diameter_piece;
             Point newPosition = new Point(oldPosition.x + dx, oldPosition.y);
             Pair<Integer,Point> move = new Pair<Integer,Point>(pair.getKey(), newPosition);
             */
             
             // move the farthest away piece that can move forward
             //move = getForwardMove(playerPieces, opponentPieces, isPlayer1, false);
             //if (move != null) moves.add(move);

             // choose best forwardish direction as next option 
//             move = getBestForwardishMove(playerPieces, opponentPieces, isPlayer1, true);
//             if (move != null) moves.add(move);             
//             move = getBestForwardishMove(playerPieces, opponentPieces, isPlayer1, false);
//             if (move != null) moves.add(move);             

             // choose valid random forwardish to less forward directions as next options
             // Can first optimize by looking at only angles you haven't already looked at
             // Ideally, would have an improved function [ getBestForwardishMove(), not yet built ] that finds the best move
             move = getRandomMove(playerPieces, opponentPieces, isPlayer1, 90);
             if (move != null) moves.add(move);             

             move = getRandomMove(playerPieces, opponentPieces, isPlayer1, 180);
             if (move != null) moves.add(move);      
                    
             move = getRandomMove(playerPieces, opponentPieces, isPlayer1, 270);
             if (move != null) moves.add(move);             
         }

         return moves;
     }

    /*   OLD
         //System.out.println("player_pieces: " + player_pieces);
		 int num_trials = 30;
		 int i = 0;

		 while (i<num_trials)
		 {
		 	// create list of your pieces in order of closeness to endzone
		 	List<Pair<Integer,Point>> orderedUnfinishedXProgress = rankXProgressOutsideEndzone(player_pieces, isplayer1);
		 	// create list of your unblocked pieces in order of closeness to endzone
		 	List<Pair<Integer,Point>> ourOrderedUnblockedPieces = findOurUnblockedPieces(player_pieces, opponent_pieces, isplayer1);
		 	
		 	// default integer piece id is random, and will resort to this when better options exhausted
		 	Integer piece_id = random.nextInt(n); 
		 	
		 	// when at least half the number of pieces are unfinished, focus on getting closest unblocked ones to other side
		 	if (orderedUnfinishedXProgress.size() >= player_pieces.size()/2) {
		 	    
		 	    // if there are unblocked unfinished pieces:
		 	    if (ourOrderedUnblockedPieces.size() > 0) {
		 	        // find piece id closest to finish line that is also unblocked
		 	        //Pair pair = ourOrderedUnblockedPieces.get(0);
		 	        piece_id = ourOrderedUnblockedPieces.get(0).getKey();
    		 	    //piece_id = pair.getKey();
    		 	    Point curr_position = player_pieces.get(piece_id);
                    i++;
                    continue;
		 	    }
		 	}
		 	
		 	Point curr_position = player_pieces.get(piece_id);
		 	
		 	//System.out.println("curr_position: " + curr_position);
		 	if(inEndZone(curr_position.x, isplayer1)) {
		 		i++;
		 		continue;
		 	}
		 	
		 	// if not in endzone, and no clear path, set a random piece to original random path settings
		 	Point new_position = new Point(curr_position);
            
            

		 	double theta = -Math.PI/2 + Math.PI * random.nextDouble();
		 	if((isplayer1 && curr_position.x > 20 )|| (!isplayer1 && curr_position.x < -20)) {
		 		theta = 0;
		 	}
		 	double delta_x = diameter_piece * Math.cos(theta);
		 	double delta_y = diameter_piece * Math.sin(theta);

		 	Double val = (Math.pow(delta_x,2) + Math.pow(delta_y, 2));
		 	// System.out.println("delta_x^2 + delta_y^2 = " + val.toString() + " theta values are " +  Math.cos(theta) + " " +  Math.sin(theta) + " diameter is " + diameter_piece);
		 	// Log.record("delta_x^2 + delta_y^2 = " + val.toString() + " theta values are " +  Math.cos(theta) + " " +  Math.sin(theta) + " diameter is " + diameter_piece);

		 	new_position.x = isplayer1 ? new_position.x - delta_x : new_position.x + delta_x;
		 	new_position.y += delta_y;
		 	Pair<Integer, Point> move = new Pair<Integer, Point>(piece_id, new_position);

		 	Double dist = Board.getdist(player_pieces.get(move.getKey()), move.getValue());
		 	// System.out.println("distance from previous position is " + dist.toString());
		 	// Log.record("distance from previous position is " + dist.toString());

		 	if (check_validity(move, player_pieces, opponent_pieces)) {
		 	    if isMovesQuotaReached(num_moves, moves, move) return moves;
		 	}
		 	i++;
		 }
		 // Testing functions:
		 //List<Pair<Integer,Point>> orderedUnfinishedXProgress = rankXProgressOutsideEndzone(player_pieces, isplayer1);
		 //List<Pair<Integer,Point>> orderedUnfinishedXProgress = rankXProgressOutsideEndzone(opponent_pieces, !isplayer1);
		 //List<Pair<Integer,Point>> orderedXProgress = rankXProgress(player_pieces, isplayer1);
		 //List<Pair<Integer,Point>> orderedXProgress = rankXProgress(opponent_pieces, !isplayer1);
		 
		 //List<Pair<Integer,Point>> unblockedOpponents = findUnblockedOpponents(opponent_pieces, isplayer1);
		 //List<Pair<Integer,Point>> idealWallPieces = findIdealWallPieces(player_pieces, isplayer1);
		 return moves;
	}*/

	public boolean check_validity( Pair<Integer, Point> move, HashMap<Integer, Point> player_pieces, HashMap<Integer, Point> opponent_pieces ) {
        boolean valid = true;
       
        // check if move is adjacent to previous position.
        if(!Board.almostEqual(Board.getdist(player_pieces.get(move.getKey()), move.getValue()), diameter_piece))
            {
                return false;
            }
        // check for collisions
        valid = valid && !Board.check_collision(player_pieces, move);
        valid = valid && !Board.check_collision(opponent_pieces, move);

        // check within bounds
        valid = valid && Board.check_within_bounds(move);
        return valid;
    }
}

// Hungarian Algorithm class
class HungarianAlgorithm {
    private final double[][] costMatrix;
    private final int rows, cols, dim;
    private final double[] labelByWorker, labelByJob;
    private final int[] minSlackWorkerByJob;
    private final double[] minSlackValueByJob;
    private final int[] matchJobByWorker, matchWorkerByJob;
    private final int[] parentWorkerByCommittedJob;
    private final boolean[] committedWorkers;
    
    /**
     * Construct an instance of the algorithm.
     * 
     * @param costMatrix
     *          the cost matrix, where matrix[i][j] holds the cost of assigning
     *          worker i to job j, for all i, j. The cost matrix must not be
     *          irregular in the sense that all rows must be the same length; in
     *          addition, all entries must be non-infinite numbers.
     */
     
     HungarianAlgorithm(double[][] costMatrix) {
        this.dim = Math.max(costMatrix.length, costMatrix[0].length);
        this.rows = costMatrix.length;
        this.cols = costMatrix[0].length;
        this.costMatrix = new double[this.dim][this.dim];
        for (int w = 0; w < this.dim; w++) {
            if (w < costMatrix.length) {
                if (costMatrix[w].length != this.cols) {
                    throw new IllegalArgumentException("Irregular cost matrix");
                }
                for (int j = 0; j < this.cols; j++) {
                    if (Double.isInfinite(costMatrix[w][j])) {
                        throw new IllegalArgumentException("Infinite cost");
                    }
                    if (Double.isNaN(costMatrix[w][j])) {
                        throw new IllegalArgumentException("NaN cost");
                    }
                }
                this.costMatrix[w] = Arrays.copyOf(costMatrix[w], this.dim);
            } else {
                    this.costMatrix[w] = new double[this.dim];
            }
        }
        labelByWorker = new double[this.dim];
        labelByJob = new double[this.dim];
        minSlackWorkerByJob = new int[this.dim];
        minSlackValueByJob = new double[this.dim];
        committedWorkers = new boolean[this.dim];
        parentWorkerByCommittedJob = new int[this.dim];
        matchJobByWorker = new int[this.dim];
        Arrays.fill(matchJobByWorker, -1);
        matchWorkerByJob = new int[this.dim];
        Arrays.fill(matchWorkerByJob, -1);
    }

  /**
   * Compute an initial feasible solution by assigning zero labels to the
   * workers and by assigning to each job a label equal to the minimum cost
   * among its incident edges.
   */
    protected void computeInitialFeasibleSolution() {
      for (int j = 0; j < dim; j++) {
        labelByJob[j] = Double.POSITIVE_INFINITY;
      }
      for (int w = 0; w < dim; w++) {
        for (int j = 0; j < dim; j++) {
          if (costMatrix[w][j] < labelByJob[j]) {
            labelByJob[j] = costMatrix[w][j];
          }
        }
      }
    }

  /**
   * Execute the algorithm.
   * 
   * @return the minimum cost matching of workers to jobs based upon the
   *         provided cost matrix. A matching value of -1 indicates that the
   *         corresponding worker is unassigned.
   */
    public int[] execute() {
        /*
         * Heuristics to improve performance: Reduce rows and columns by their
         * smallest element, compute an initial non-zero dual feasible solution and
         * create a greedy matching from workers to jobs of the cost matrix.
         */
        reduce();
        computeInitialFeasibleSolution();
        greedyMatch();
        
        int w = fetchUnmatchedWorker();
        while (w < dim) {
            initializePhase(w);
            executePhase();
            w = fetchUnmatchedWorker();
        }
        int[] result = Arrays.copyOf(matchJobByWorker, rows);
        for (w = 0; w < result.length; w++) {
            if (result[w] >= cols) {
                result[w] = -1;
            }
        }
        return result;
    }

  /**
   * Execute a single phase of the algorithm. A phase of the Hungarian algorithm
   * consists of building a set of committed workers and a set of committed jobs
   * from a root unmatched worker by following alternating unmatched/matched
   * zero-slack edges. If an unmatched job is encountered, then an augmenting
   * path has been found and the matching is grown. If the connected zero-slack
   * edges have been exhausted, the labels of committed workers are increased by
   * the minimum slack among committed workers and non-committed jobs to create
   * more zero-slack edges (the labels of committed jobs are simultaneously
   * decreased by the same amount in order to maintain a feasible labeling).
   * <p>
   * 
   * The runtime of a single phase of the algorithm is O(n^2), where n is the
   * dimension of the internal square cost matrix, since each edge is visited at
   * most once and since increasing the labeling is accomplished in time O(n) by
   * maintaining the minimum slack values among non-committed jobs. When a phase
   * completes, the matching will have increased in size.
   */
    protected void executePhase() {
      while (true) {
        int minSlackWorker = -1, minSlackJob = -1;
        double minSlackValue = Double.POSITIVE_INFINITY;
        for (int j = 0; j < dim; j++) {
          if (parentWorkerByCommittedJob[j] == -1) {
            if (minSlackValueByJob[j] < minSlackValue) {
              minSlackValue = minSlackValueByJob[j];
              minSlackWorker = minSlackWorkerByJob[j];
              minSlackJob = j;
            }
          }
        }
        if (minSlackValue > 0) {
          updateLabeling(minSlackValue);
        }
        parentWorkerByCommittedJob[minSlackJob] = minSlackWorker;
        if (matchWorkerByJob[minSlackJob] == -1) {
          /*
           * An augmenting path has been found.
           */
          int committedJob = minSlackJob;
          int parentWorker = parentWorkerByCommittedJob[committedJob];
          while (true) {
            int temp = matchJobByWorker[parentWorker];
            match(parentWorker, committedJob);
            committedJob = temp;
            if (committedJob == -1) {
              break;
            }
            parentWorker = parentWorkerByCommittedJob[committedJob];
          }
          return;
        } else {
          /*
           * Update slack values since we increased the size of the committed
           * workers set.
           */
          int worker = matchWorkerByJob[minSlackJob];
          committedWorkers[worker] = true;
          for (int j = 0; j < dim; j++) {
            if (parentWorkerByCommittedJob[j] == -1) {
              double slack = costMatrix[worker][j] - labelByWorker[worker]
                  - labelByJob[j];
              if (minSlackValueByJob[j] > slack) {
                minSlackValueByJob[j] = slack;
                minSlackWorkerByJob[j] = worker;
              }
            }
          }
        }
      }
    }
    
    /**
     * 
     * @return the first unmatched worker or {@link #dim} if none.
     */
    protected int fetchUnmatchedWorker() {
      int w;
      for (w = 0; w < dim; w++) {
        if (matchJobByWorker[w] == -1) {
          break;
        }
      }
      return w;
    }
    
    /**
     * Find a valid matching by greedily selecting among zero-cost matchings. This
     * is a heuristic to jump-start the augmentation algorithm.
     */
    protected void greedyMatch() {
      for (int w = 0; w < dim; w++) {
        for (int j = 0; j < dim; j++) {
          if (matchJobByWorker[w] == -1 && matchWorkerByJob[j] == -1
              && costMatrix[w][j] - labelByWorker[w] - labelByJob[j] == 0) {
            match(w, j);
          }
        }
      }
    }
    
    /**
     * Initialize the next phase of the algorithm by clearing the committed
     * workers and jobs sets and by initializing the slack arrays to the values
     * corresponding to the specified root worker.
     * 
     * @param w
     *          the worker at which to root the next phase.
     */
    protected void initializePhase(int w) {
      Arrays.fill(committedWorkers, false);
      Arrays.fill(parentWorkerByCommittedJob, -1);
      committedWorkers[w] = true;
      for (int j = 0; j < dim; j++) {
        minSlackValueByJob[j] = costMatrix[w][j] - labelByWorker[w]
            - labelByJob[j];
        minSlackWorkerByJob[j] = w;
      }
    }
    
    /**
     * Helper method to record a matching between worker w and job j.
     */
    protected void match(int w, int j) {
      matchJobByWorker[w] = j;
      matchWorkerByJob[j] = w;
    }
    
    /**
     * Reduce the cost matrix by subtracting the smallest element of each row from
     * all elements of the row as well as the smallest element of each column from
     * all elements of the column. Note that an optimal assignment for a reduced
     * cost matrix is optimal for the original cost matrix.
     */
    protected void reduce() {
      for (int w = 0; w < dim; w++) {
        double min = Double.POSITIVE_INFINITY;
        for (int j = 0; j < dim; j++) {
          if (costMatrix[w][j] < min) {
            min = costMatrix[w][j];
          }
        }
        for (int j = 0; j < dim; j++) {
          costMatrix[w][j] -= min;
        }
      }
      double[] min = new double[dim];
      for (int j = 0; j < dim; j++) {
        min[j] = Double.POSITIVE_INFINITY;
      }
      for (int w = 0; w < dim; w++) {
        for (int j = 0; j < dim; j++) {
          if (costMatrix[w][j] < min[j]) {
            min[j] = costMatrix[w][j];
          }
        }
      }
      for (int w = 0; w < dim; w++) {
        for (int j = 0; j < dim; j++) {
          costMatrix[w][j] -= min[j];
        }
      }
    }
    
    /**
     * Update labels with the specified slack by adding the slack value for
     * committed workers and by subtracting the slack value for committed jobs. In
     * addition, update the minimum slack values appropriately.
     */
    protected void updateLabeling(double slack) {
    for (int w = 0; w < dim; w++) {
      if (committedWorkers[w]) {
        labelByWorker[w] += slack;
      }
    }
    for (int j = 0; j < dim; j++) {
      if (parentWorkerByCommittedJob[j] != -1) {
        labelByJob[j] -= slack;
      } else {
        minSlackValueByJob[j] -= slack;
      }
    }
  }
}