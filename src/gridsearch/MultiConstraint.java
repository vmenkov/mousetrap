package gridsearch;

import java.util.*;

/** Represents a set of linear constraints. Each one of them may be
    of any form supported by the Constraint interface.
*/
public class MultiConstraint implements Constraint {
    /** The underlying constraints */
    Vector<Constraint> v = new Vector<Constraint>();

    void addSimplexConstraint(int[] ind) {
	SingleSparseConstraint c = SingleSparseConstraint.simplex(ind);
	v.add(c);
    }

    /** Adds a new constraint to this MultiConstraint if the new
	constrain is different from all constraints already in this
	MultiConstraint. This method assumes that all already included
	constraints have their variables ordered in ascending order,
	and it orders the vars in the new constraint in the same way as well.
     */
    void addSimplexConstraintIfUnique(int[] ind) {
	int z[] = Arrays.copyOfRange(ind, 0, ind.length);
	Arrays.sort(z);
	addSimplexConstraintIfUniqueSorted(z);
    }

    void addSimplexConstraintIfUnique(HashSet<Integer> h) {
      int z[] = new int[h.size()];
      int p=0;
      for(Integer a: z) { z[p++] = a.intValue(); }
      addSimplexConstraintIfUnique(z);
    }

    /** Assumes ind[] is sorted */
    private void addSimplexConstraintIfUniqueSorted(int[] ind) {
	SingleSparseConstraint c = SingleSparseConstraint.simplex(ind);
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


}

