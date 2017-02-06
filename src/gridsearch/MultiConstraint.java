package gridsearch;

import java.util.Vector;

/** Represents a set of linear constraints. Each one of them may be
    of any form supported by the Constraint interface.
*/
class MultiConstraint implements Constraint {
    /** The underlying constraints */
    Vector<Constraint> v = new Vector<Constraint>();

    void addSimplexConstraint(int[] ind) {
	SingleSparseConstraint c = SingleSparseConstraint.simplex(ind);
	v.add(c);
    }

    /** Checks that the constraint is all-non-negative */
    public void validate() {
	for(Constraint c: v) {
	    c.validate();
	}	
    }

    public boolean holds(double x[]) {
	for(Constraint c: v) {
	    if (!c.holds(x)) return false;
	}
	return true;
    }
    public boolean holds(int x[]) {
	for(Constraint c: v) {
	    if (!c.holds(x)) return false;
	}
	return true;
    }

    public Constraint constraintInt(Grid g) {
	MultiConstraint q = new 	MultiConstraint();
	for(Constraint c: v) {
	    q.v.add( c.constraintInt(g));
	}
	q.validate();
	return q;
    }

}

