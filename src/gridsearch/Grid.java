package gridsearch;

/** Represents a section of R^n space, with a regular grid on it */
class Grid {
    /** The "lower" and "upper" corners */
    ParVec corners[];
    /** For dimension i, the range of x[i] is split into m[i] equal segments */
    int [] m;
    int dim() { return m.length; }
    /** width of one cell of this grid along the i-th dimension */
    double cellWidth(int i) { return (corners[1].x(i) -corners[0].x(i))/m[i];}

    private Grid(ParVec _corners[], int [] _m) {
	corners = _corners;
	m = _m;
    }

    private Grid(double c[][], int [] _m) {
	corners = new ParVec[] { new ParVec(c[0]),new ParVec(c[1])};
	m = _m;
    }
    
    /** Creates a new grid covering a subdomain of the domain covered
	by this grid. Normally, the new grid has its center at the
	specified position, and the width of 2*L cells (L cells from
	the center in each direction) of the original grid in each
	dimension.  The size is reduced in those directions/dimensions
	where it would otherwise go outside the bounds of the original grid. 
	
	@param center The center point of the new grid (if unaffected
	by the old grid's external boundaries)
	@param mfactor Into how many cells should each cell of the original grid
	be divided (in each direction)
	@param L how many cells of the original grid, in each direction,
	should the new grid cover
     */
    Grid vicinityGrid(ParVec center, int mfactor, int L) {
	final int n = dim();
	int mnew[] = new int[n];
	double c[][] = new double[2][];
	for(int k=0; k<2; k++) c[k]=new double[n];
	for(int i=0; i<n; i++) {
	    int md[] = new int[2];
	    double box = cellWidth(i);
	    for(int k=0; k<2; k++) {
		int sign = 2*k - 1;
		double boundary = corners[k].x(i);
		boolean reduce = Math.abs( center.x(i) - boundary) < box*L;
		c[k][i] = reduce?  boundary:  center.x(i) + sign * box*L;
		md[k] = (int)Math.round((Math.abs(center.x(i)-c[k][i])*mfactor)/box);
		// take care of any ill effects of rounding
		if (md[k] * 2 >= (2*L-1)*mfactor) md[k] = L*mfactor;
	    }
	    mnew[i] = md[0] + md[1];
	}
	return new Grid( c, mnew);
    }
}
