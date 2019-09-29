package flip.g6;

import java.util.HashMap;
import flip.sim.Point;
import javafx.util.Pair;

public class Aggressive extends Move {

	private HashMap<Integer, Point> player_pieces;
	private HashMap<Integer, Point> opponent_pieces;
	private boolean isplayer1;
	private Integer n;
	private Double diameter_piece;
	public boolean RUNNER_STRATEGY_SET = true;
	private int highPieceID = -1;
	private int lowPieceID = -1;
	int counter = 0;

	public Aggressive(HashMap<Integer, Point> player_pieces, HashMap<Integer, Point> opponent_pieces, boolean isplayer1, Integer n, Double diameter_piece) {
		this.player_pieces = player_pieces;
		this.opponent_pieces = opponent_pieces;
		this.isplayer1 = isplayer1;
		this.n = n;
		this.diameter_piece = diameter_piece;
	}

	@Override
	public void updatePieceInfo(HashMap<Integer, Point> player_pieces, HashMap<Integer, Point> opponent_pieces) {
		this.player_pieces = player_pieces;
		this.opponent_pieces = opponent_pieces;
	}

	@Override
	public Pair<Integer, Point> getMove() {

		Pair<Integer, Point> move = null;
		if(n == 1) {			
			for(Integer piece_id : player_pieces.keySet()) {
				Point currPosition = player_pieces.get(piece_id);
				Point newPosition;

				double theta = 0;
				double numAngles = 180;
				for(int i = 0; i <= numAngles; i++) {
					if(i % 2 == 0)
						theta -= i * Math.PI / numAngles;
					else
						theta += i * Math.PI / numAngles;
					double deltaX = diameter_piece * Math.cos(theta);
					double deltaY = diameter_piece * Math.sin(theta);

					newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
					move = new Pair<Integer, Point>(piece_id, newPosition);

					if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
						player_pieces.put(piece_id, newPosition);
						return move;
					}
				}
			}
			return null;
		}
		else if(n <= 15) {
			HashMap<Integer, Point> unfinished_pieces = getUnfinishedPlayerPieces(player_pieces, isplayer1, Approach.AGGRESSIVE);
			HashMap<Integer, Point> closest_pieces = getClosestPointsToOpponentBoundary(unfinished_pieces.size(), unfinished_pieces, isplayer1);
			double posMarker = isplayer1 ? 15.0 : -15.0;
			for(Integer piece_id : closest_pieces.keySet()) {
				if((isplayer1 && (closest_pieces.get(piece_id).x >= posMarker)) || (!isplayer1 && (closest_pieces.get(piece_id).x <= posMarker))) {
					Point currPosition = player_pieces.get(piece_id);
					Point newPosition;
					double theta = 0;
					double deltaX = diameter_piece * Math.cos(theta);
					double deltaY = diameter_piece * Math.sin(theta);

					newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
					move = new Pair<Integer, Point>(piece_id, newPosition);

					if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
						player_pieces.put(piece_id, newPosition);
						return move;
					}
				}
			}
			for(Integer piece_id : closest_pieces.keySet()) {
				if((isplayer1 && (closest_pieces.get(piece_id).x < posMarker)) || (!isplayer1 && (closest_pieces.get(piece_id).x > posMarker))) {
					Point currPosition = player_pieces.get(piece_id);
					Point newPosition;
					double theta = 0;
					double deltaX = diameter_piece * Math.cos(theta);
					double deltaY = diameter_piece * Math.sin(theta);

					newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
					move = new Pair<Integer, Point>(piece_id, newPosition);

					if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
						player_pieces.put(piece_id, newPosition);
						return move;
					}
				}
			}
			for(Integer piece_id : closest_pieces.keySet()) {
				Point currPosition = player_pieces.get(piece_id);
				Point newPosition;

				double theta = 0;
				double numAngles = 180;
				for(int i = 0; i <= numAngles; i++) {
					if(i % 2 == 0)
						theta -= i * Math.PI / numAngles;
					else
						theta += i * Math.PI / numAngles;
					double deltaX = diameter_piece * Math.cos(theta);
					double deltaY = diameter_piece * Math.sin(theta);

					newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
					move = new Pair<Integer, Point>(piece_id, newPosition);

					if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
						player_pieces.put(piece_id, newPosition);
						return move;
					}
				}
			}
			return null;		}
		else {
			double startingBoundaryForAggressive = isplayer1 ? 19.5 : -19.5;
			HashMap<Integer, Point> unfinished_pieces = getUnfinishedPlayerPieces(player_pieces, isplayer1, Approach.AGGRESSIVE);
			HashMap<Integer, Point> closest_pieces = getClosestPointsToOpponentBoundary(unfinished_pieces.size(), unfinished_pieces, isplayer1);
			for(Integer pieceID : closest_pieces.keySet()) {
				if((isplayer1 && closest_pieces.get(pieceID).x < startingBoundaryForAggressive) || (!isplayer1 && closest_pieces.get(pieceID).x > startingBoundaryForAggressive)) {
					Point currPosition = player_pieces.get(pieceID);
					Point newPosition;

					double theta = 0;
					double numAngles = 180;
					for(int i = 0; i <= numAngles; i++) {
						if(i % 2 == 0)
							theta -= i * Math.PI / numAngles;
						else
							theta += i * Math.PI / numAngles;
						double deltaX = diameter_piece * Math.cos(theta);
						double deltaY = diameter_piece * Math.sin(theta);

						newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
						move = new Pair<Integer, Point>(pieceID, newPosition);

						if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
							player_pieces.put(pieceID, newPosition);
							return move;
						}
					}
				}
			}
		}
		return null;
	}

	public Pair<Integer, Point> run(Integer numMoves, boolean isEvenMove) {
		
		// Top coin move
		if(highPieceID == -1) {
			Point highestPoint = isplayer1 ? new Point(20.0, 20.0) : new Point(-20.0, 20.0);
			highPieceID = 0;
			double distance = Double.MAX_VALUE;
			for(Integer index : player_pieces.keySet()) {
				double newDistance = Math.sqrt(Math.pow(highestPoint.x - player_pieces.get(index).x, 2) + Math.pow(highestPoint.y - player_pieces.get(index).y, 2));
				if(newDistance < distance) {
					distance = newDistance;
					highPieceID = index;
				}
			}
		}
		
		// Lowest coin move
		if(lowPieceID == -1) {
			Point lowestPoint = isplayer1 ? new Point(20.0, -20.0) : new Point(-20.0, -20.0);
			lowPieceID = 0;
			double distance = Double.MAX_VALUE;
			for(Integer index : player_pieces.keySet()) {
				double newDistance = Math.sqrt(Math.pow(lowestPoint.x - player_pieces.get(index).x, 2) + Math.pow(lowestPoint.y - player_pieces.get(index).y, 2));
				if(newDistance < distance) {
					distance = newDistance;
					lowPieceID = index;
				}
			}
			
		}
		
		if(isEvenMove && counter < 5) {
			Point currPosition = player_pieces.get(highPieceID);
			Point newPosition = new Point(currPosition);
			double theta = 0;
			double deltaX = diameter_piece * Math.cos(theta);
			double deltaY = diameter_piece * Math.sin(theta);
			newPosition.x = isplayer1 ? newPosition.x - deltaX: newPosition.x + deltaX;
			newPosition.y += deltaY;
			player_pieces.put(highPieceID, newPosition);
			counter++;
			return new Pair<Integer, Point>(highPieceID, newPosition);
		} 
		else {
			Point currPosition = player_pieces.get(lowPieceID);
			Point newPosition = new Point(currPosition);
			double theta = 0;
			double numAngles = 180;
			for(int i = 0; i <= numAngles; i++) {
				if(i % 2 == 0)
					theta -= i * Math.PI / numAngles;
				else
					theta += i * Math.PI / numAngles;
				double deltaX = diameter_piece * Math.cos(theta);
				double deltaY = diameter_piece * Math.sin(theta);

				newPosition = isplayer1 ? new Point(currPosition.x - deltaX, currPosition.y + deltaY) : new Point(currPosition.x + deltaX, currPosition.y + deltaY);
				Pair<Integer, Point> move = new Pair<Integer, Point>(lowPieceID, newPosition);

				if(checkValidity(move, player_pieces, opponent_pieces, diameter_piece)) {
					player_pieces.put(lowPieceID, newPosition);
					if((isplayer1 && newPosition.x < -22.0) || (!isplayer1 && newPosition.x > 22.0))
						this.RUNNER_STRATEGY_SET = false;
					return move;
				}
			}
			this.RUNNER_STRATEGY_SET = false;
		}
		return null;
	}
	
	public Integer getHighPieceID() {
		return highPieceID;
	}

	public Integer getLowPieceID() {
		return lowPieceID;
	}
	
	public HashMap<Integer, Point> getPlayerPieces() {
		return player_pieces;
	}

	public void setPlayerPieces(HashMap<Integer, Point> player_pieces) {
		this.player_pieces = player_pieces;
	}

	public HashMap<Integer, Point> getOpponentPieces() {
		return opponent_pieces;
	}

	public void setOpponentPieces(HashMap<Integer, Point> opponent_pieces) {
		this.opponent_pieces = opponent_pieces;
	}	
}