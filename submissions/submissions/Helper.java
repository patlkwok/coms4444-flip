public static double roundEven(double d) {
	return Math.round(d / 2) * 2;
}

public class FlipHelperMethods {

	// Returns a list of x positions and percent completions for possible opponent walls
	public List<Pair<Double,Double>> detectWall(HashMap<Integer,Point> opponentPieces, boolean isplayer1)
	{
		double board_left = Board.x_min;
		double board_right = Board.x_max;
		double diameter_piece = 2.0;  // probably optional
		HashMap<Double,Integer> walls = new HashMap<Double,Integer>;
		for (Point opponentPiece : opponentPieces.values()) {
			double approx_x = roundEven(opponentPiece.x);
			int approx_x_count = walls.containsKey(approx_x) ? walls.get(approx_x) : 0;
			walls.put(approx_x, approx_x_count + 1);
		}
		List<Pair<Double,Double>> possible_walls = new List<Pair<Double,Double>>;
		for (double wall : walls.keySet()) {
			int len = walls.get(wall);
			if (len > 5) {
				double percentage = len / 12.0;
				Pair<Double,Double> possible_wall = new Pair<Double,Double>(wall, percentage);
				possible_walls.add(possible_wall);
			}
		}
		return possible_walls;
	}

	// Compare opponent pieces with opponent pieces from last turn and return list of their moves
	public List<Pair<Integer,Point>> findLastOpponentMoves(HashMap<Integer,Point> currentOpponentPieces, HashMap<Integer,Point> lastOpponentPieces, boolean isplayer1)
	{
		// To make this method work we need to record the board position after each player turn
		List<Pair<Integer,Point>> moves = new List<Pair<Integer,Point>>;
		for (int piece : currentOpponentPieces.keySet()) {
			Point currentLocation = currentOpponentPieces.get(piece);
			Point lastLocation = lastOpponentPieces.get(piece);
			boolean moved = currentLocation.equals(lastLocation);
			if (moved) {
				double delta_x = currentLocation.x - lastLocation.x;
				double delta_y = currentLocation.y - lastLocation.y;
				Point delta = new Point(delta_x, delta_y);
				Pair<Integer,Point> delta_piece = new Pair<Integer,Point>(piece, delta);
				moves.add(delta_piece);
			}
		}
		return moves;
	}
