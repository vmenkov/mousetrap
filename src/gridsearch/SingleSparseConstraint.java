package gridsearch;

import java.util.*;

/** Represents a linear constraint, of the form 
    sum_j   (a[j]*x[j]) &le; b, 
    with a sparse non-negative a[]
*/
class SingleSparseConstraint implements Constraint{
    int aInd[];
    double aVal[];
    double b;
    private SingleSparseConstraint() {}
    SingleSparseConstraint(	int _aInd[], double _aVal[], double _b) {
	aInd = _aInd;
	aVal = _aVal;
	b = _b;
    }
    
    /** Creates a constraint of the form 
	sum_{i in Z} x_i &le; 1
	for a specified set of index Z
	@param _aInd The set of indexes (Z). The caller is expected to ensure
	uniqueness and order.
    */
    static SingleSparseConstraint simplex(int _aInd[]) {
	SingleSparseConstraint c = new
	    SingleSparseConstraint(_aInd, new double[_aInd.length], 1.0);
	for(int i=0; i<c.aInd.length; i++) {
	    c.aVal[i] = 1.0;
	}
	return c;
    }


    public void validate() {
	if (aInd.length != aVal.length) throw new IllegalArgumentException("SSConstraint: length mismatch");
	for(int i=0; i<aVal.length; i++) {
	    if (aVal[i]<0) throw new IllegalArgumentException("Constraints with negative coefficients not supported");
	}       	
    }

    
    public boolean holds(double x[]) {
	double sum = 0;
	for(int i=0; i<aVal.length; i++) {
	    sum += aVal[i] * x[aInd[i]];
	}
	return sum <= b;
    }
    public boolean holds(int x[]) {
	double sum = 0;
	for(int i=0; i<aVal.length; i++) {
	    sum += aVal[i] * x[aInd[i]];
	}
	return sum <= b;
    }

    public Constraint constraintInt(Grid g) {
	SingleSparseConstraint q = new SingleSparseConstraint();
	q.aVal = new double[aVal.length];
	q.aInd = new int[aInd.length];
	q.b = 1.0;
	for(int i=0; i<aInd.length; i++) {
	    int k = q.aInd[i];
	    if (k>= g.dim()) throw new IllegalArgumentException();
	    q.b -= g.corners[0].x(k);
	    q.aVal[i] = aVal[i] * g.cellWidth(k);
	}
	q.validate();
	return q;
    }
  	
    public boolean equals(Object x) {
	if (!(x instanceof SingleSparseConstraint)) return false;
	SingleSparseConstraint z= (SingleSparseConstraint)x;
	return z.b == b && 
	    Arrays.equals(z.aInd,  aInd) &&
	    Arrays.equals(z.aVal,  aVal);
    }

}

