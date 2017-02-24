package gridsearch;

import java.util.*;

/** Represents a set of linear constraints. Each one of them may be
    of any form supported by the Constraint interface.
*/
public class MultiConstraint extends Constraint {
    /** The underlying constraints */
    Vector<Constraint> v = new Vector<Constraint>();

    void addSimplexConstraint(int[] ind) {
	SingleSparseConstraint c = SingleSparseConstraint.simplex(ind);
	v.add(c);
    }

    /** Adds a new constraint to this MultiConstraint if the new
	constrain is different from all constraints already in this
	MultiConstraint. This method sorts ind[], and interprets
	non-unique values as having higher weight.
     */
    void addSimplexConstraintIfUnique(int[] ind) {
	SingleSparseConstraint c =  SingleSparseConstraint.simplexWeighted(ind);
	if (c.isTrivialInCube()) {	
	    System.out.println("Not adding constraint " + c + ", because it's trivial");
	    return;
	}
	for(Constraint o: v) {
	    if (o.equals(c)) {
		System.out.println("Not adding constraint " + c + ", because it's already in");
		return;
	    }
	}
	System.out.println("Adding constraint " + c );
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

   public boolean equals(Object x) {
	if (!(x instanceof MultiConstraint)) return false;
	MultiConstraint z= (MultiConstraint)x;
	if (z.v.size() != v.size())  return false;
	for(int i=0; i<v.size(); i++) {
	    if (z.v.elementAt(i).equals(v.elementAt(i))) return false;
	}
	return true;
    }

    public String toString() {
	StringBuffer s = new StringBuffer("(\n");
	for(Constraint c: v) {
	    s.append("" + c + "\n");
	}
	s.append(")");
	return s.toString();  
    }

}

