package gridsearch;

import java.util.*;

/** Represents a linear constraint, of the form 
    sum_j   (a[j]*x[j]) &le; b, 
    with a sparse non-negative a[]
*/
class SingleSparseConstraint extends Constraint{
    final int aInd[];
    final double aVal[];
    final double b;
    /** Createas a constraint with specified indexes and coefficients.
	@param _aInd Must be listed in increasing order, without duplicates
	@param _aVal Must all be non-negative.
    */
    SingleSparseConstraint(int _aInd[], double _aVal[], double _b) {
	aInd = _aInd;
	aVal = _aVal;
	b = _b;
  	validate();
    }
    
    /** Creates a constraint of the form 
	sum_{i in Z} x_i &le; 1
	for a specified set of index Z
	@param _aInd The set of indexes (Z). The caller is expected to ensure
	uniqueness and order; otherwise, an exception is thrown.
    */
    static SingleSparseConstraint simplex(int _aInd[]) {
	double [] _aVal = new double[_aInd.length];
	Arrays.fill(_aVal, 1.0);
	return new SingleSparseConstraint(_aInd, _aVal, 1.0);
    }

    /** 
	@param _aInd The indexes don't have to be ordered or unique. 
	The data will be sorted, and duplicate counts will be
	used as weights
     */
    static SingleSparseConstraint simplexWeighted(int _aInd[]) {
	int z[] = Arrays.copyOfRange(_aInd, 0, _aInd.length);
	Arrays.sort(z);
	int uniqCnt=0;
	for(int i=0; i<z.length; i++) {
	    if (i==0 || z[i] != z[i-1]) uniqCnt++;
	}
	int aInd[] = new int[uniqCnt];
	double aVal[] = new double[uniqCnt];
	int p=0;
	for(int i=0; i<z.length; i++) {
	    if (i==0 || z[i] != z[i-1]) {
		aInd[p] = z[i];
		aVal[p] = 1.0;
		p++;
	    } else {
		aVal[p-1] += 1.0;
	    }
	}
	return new SingleSparseConstraint(aInd, aVal, 1.0);
    }



    /** Checks that the indexes are listed in decreasing order and are unique,
	and that coeffcients are non-negative. */
    public void validate() {
	if (aInd.length != aVal.length) throw new IllegalArgumentException("SSConstraint: length mismatch");
	for(int i=0; i<aVal.length; i++) {
	    if (i>0 && aInd[i]<=aInd[i-1])  throw new IllegalArgumentException("Constraint definition invalid: indexes not in order: " + this);
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
	double[] _aVal =  new double[aVal.length];
	double _b = b;
	for(int i=0; i<aInd.length; i++) {
	    int k = aInd[i];
	    if (k>= g.dim()) throw new IllegalArgumentException();
	    _b -= aVal[i] * g.corners[0].x(k);
	    _aVal[i] = aVal[i] * g.cellWidth(k);
	}
	SingleSparseConstraint q = new SingleSparseConstraint(aInd, _aVal, _b);
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

    /** Returns true if we know that all values in the unit cube satisfy this
	constraint */
    boolean isTrivialInCube() {
	return (aInd.length==1 && aVal[0]==1.0);
    }

    public String toString() {
	StringBuffer s = new StringBuffer("(");
	int cnt=0;
	for(int i=0; i<aInd.length; i++) {
	    if (cnt > 0) s.append(" + ");
	    if (aVal[i] != 1.0) {
		String w = ((double)(int)aVal[i] == aVal[i]) ? "" + (int)aVal[i]: "" + aVal[i];
		s.append(w + "*");
	    }
	    s.append("x["+ aInd[i]+"] ");
	    cnt++;
	}
	s.append(" <= " + b + ")");
	return s.toString();
    }

}

