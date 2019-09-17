public class FlipHelperMethods {
	// Returns the distance between two points
	public double dist(Point p1, Point p2)
	{
		double d = Math.sqrt(Math.pow(p2.x - p1.x, 2) + Math.pow(p2.y - p1.y, 2));
		return d;
	}

	// Returns target locations given wall pieces
	public List<Point> targetLoc(HashMap<Integer,Point> wallPieces, boolean isplayer1)
	{
		// Step 1: Sort the wall pieces by y-coordinate
		List<Point> wall = new List<Point>(wallPieces.values());  // not required if wall pieces are already provided in a list
		Collections.sort(wall, new Comparator<Point>() {
    			@Override
    			public int compare(Point p1, Point p2) {
        			return p1.y.compareTo(p2.y);
    			}
		});
		// Step 2: Calculate the center locations between the wall pieces
		List<Point> centers = new List<Point>;
		for (int i = 0; i < wall.size() - 1; i++) {
			double xc = (wall.get(i).x + wall.get(i + 1).x) / 2.0;
			double yc = (wall.get(i).y + wall.get(i + 1).y) / 2.0;
			Point c = new Point(xc, yc);  // center between wall pieces i and (i + 1)
			centers.add(c);
		}
		// Step 3: Calculate the target locations behind the wall pieces
		List<Point> targets = new List<Point>;
		for (int i = 0; i < centers.size(); i++) {
			double l = FlipHelperMethods.dist(centers.get(i), wall.get(i));
			double d = 2.0;  // diameter_piece, probably optional
			if (d > l) {
				double v = Math.sqrt(Math.pow(d, 2) - Math.pow(l, 2));  // distance between target and center locations
				double dx = wall.get(i + 1).x - wall.get(i).x;
				double dy = wall.get(i + 1).y - wall.get(i).y;  // dy should be positive
				if (!isplayer1) {
					v = v * -1.0;  // v is positive for player 1 and negative for player 2
				}
				double xt = xc + v * dy / (2.0 * l);
				double yt = yc - v * dx / (2.0 * l);
				Point t = new Point(xt, yt);  // target location behind wall pieces i and (i + 1)
			}
			else {
				Point t = centers.get(i);  // should not happen, just in case
			}
			targets.add(t);
		}
		return targets;
	}

	// Returns the optimal move given locations of non-wall pieces and targets
	public Pair<Integer,Point> optimalMove(HashMap<Integer,Point> ourPieces, HashMap<Integer,Point> wallPieces, List<Point> targets)
	{
		int optimalPiece = 0;  // overall optimal piece to move
		Point optimalTarget = new Point(0, 0);  // overall optimal target location
		double minDist = Double.POSITIVE_INFINITY;  // overall minimum distance to target
		for (int p: ourPieces.keySet()) {
			if (!wallPieces.containsKey(p)) {
				// Only consider pieces that are not wall pieces themselves
				Point optimalTarget_p = new Point(0, 0);  // closest target location from piece p 
				double minDist_p = Double.POSITIVE_INFINITY;  // minimum distance to a target location from piece p
				for (int t = 0; t < targets.size(); t++) {
					double dist_pt = FlipHelperMethods.dist(ourPieces.get(p), targets.get(t));
					if (dist_pt < minDist_p) {
						// Update optimal for piece p
						optimalTarget_p = targets.get(t);
						minDist_p = dist_pt;
					}
				}
				if (minDist_p < minDist) {
					// Update overall optimal
					optimalPiece = p;
					optimalTarget = optimalTarget_p;
					minDist = minDist_p;
				}
			}
		}
		Pair<Integer,Point> optimal = new Pair<Integer,Point>(optimalPiece, optimalTarget);  // optimal piece to move and target
		return optimal;
	}
