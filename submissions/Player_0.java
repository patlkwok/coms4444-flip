// With Becky's functions
package flip.players;
import java.util.List;
import java.util.Comparator;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.HashMap;
import javafx.util.Pair; 
import java.util.ArrayList;
import java.lang.*;

import flip.sim.Point;
import flip.sim.Board;
import flip.sim.Log;

public class Player implements flip.sim.Player
{
	private int seed = 42;
	private Random random;
	private boolean isplayer1;
	private Integer n;
	private Double diameter_piece;


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
	}

    // x coordinate
	public boolean inEndZone(double x, boolean isplayer1) {
		if((!isplayer1 && x > 20 + 1.75 * diameter_piece + ((n / 9) * (diameter_piece / 2))) || (isplayer1 && x < -20 - 1.75 * diameter_piece - ((n / 9) * (diameter_piece / 2)))) {

			return true;
		}
		return false;
	}

    // Becky
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
    
    // Becky
    // returns our (or our opponent's) sorted pieces list WITHOUT pieces that have already crossed the endzone
    private List<Pair<Integer,Point>> rankXProgressOutsideEndzone( HashMap<Integer,Point> pieces, boolean arePiecesPlayer1 ) {
    	List<Pair<Integer,Point>> orderedPieces = rankXProgress(pieces, arePiecesPlayer1);
        List<Pair<Integer,Point>> outsideEndzonePieces = new ArrayList<Pair<Integer,Point>>();

    	int sign = arePiecesPlayer1? -1 : 1;
    	// use indices.toArray() as a deep copy to avoid changing indices inside for loop
    	for (Pair<Integer,Point> pair : orderedPieces) {
    	    if ((pieces.get(pair.getKey()).x *sign) < (20 + diameter_piece/2)) {
    	        outsideEndzonePieces.add(pair);
    	    }
    	}
    	// debugging printer
    	System.out.print("Pairs:  ");
        for (Pair<Integer,Point> pair : outsideEndzonePieces) System.out.print(" (" + pair.getKey() + ", " + pair.getValue().x + ")");
        System.out.println();
        
        return outsideEndzonePieces;
    }
    
    // Becky
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
    
    // Becky
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
    
    // Becky: in progress, but Pranav, you may have better intuition about how to complete this
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
    
    
	public List<Pair<Integer, Point>> getMoves(Integer num_moves, HashMap<Integer, Point> player_pieces, HashMap<Integer, Point> opponent_pieces, boolean isplayer1)
	{
		 List<Pair<Integer, Point>> moves = new ArrayList<Pair<Integer, Point>>();

         //System.out.println("player_pieces: " + player_pieces);
		 int num_trials = 30;
		 int i = 0;

		 while(moves.size()!= num_moves && i<num_trials)
		 {
		 	Integer piece_id = random.nextInt(n);
		 	Point curr_position = player_pieces.get(piece_id);
		 	//System.out.println("curr_position: " + curr_position);
		 	if(inEndZone(curr_position.x, isplayer1)) {
		 		i++;
		 		continue;
		 	}
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

		 	if(check_validity(move, player_pieces, opponent_pieces))
		 		moves.add(move);
		 	i++;
		 }
		 // Testing Becky's functions:
		 //List<Pair<Integer,Point>> orderedUnfinishedXProgress = rankXProgressOutsideEndzone(player_pieces, isplayer1);
		 //List<Pair<Integer,Point>> orderedUnfinishedXProgress = rankXProgressOutsideEndzone(opponent_pieces, !isplayer1);
		 //List<Pair<Integer,Point>> orderedXProgress = rankXProgress(player_pieces, isplayer1);
		 //List<Pair<Integer,Point>> orderedXProgress = rankXProgress(opponent_pieces, !isplayer1);
		 
		 //List<Pair<Integer,Point>> unblockedOpponents = findUnblockedOpponents(opponent_pieces, isplayer1);
		 //List<Pair<Integer,Point>> idealWallPieces = findIdealWallPieces(player_pieces, isplayer1);
		 return moves;
	}

	public boolean check_validity(Pair<Integer, Point> move, HashMap<Integer, Point> player_pieces, HashMap<Integer, Point> opponent_pieces)
    {
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