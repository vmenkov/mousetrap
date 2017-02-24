package gridsearch;

/** Represents a constraint of the form sum_j (a[j]*x[j]) &le; b,
    with non-negative a[]
*/
public abstract class Constraint {
    /** Checks that the constraint is all-non-negative */
    abstract public void validate();
    abstract public boolean holds(double x[]);
    public boolean holds(ParVec pv) { return holds(pv.getX()); }
    abstract public boolean holds(int x[]);
    /** This constraint is applied to the integer coordinates of grid points (on
	the [0..m] scale, rather than [0..1])  */
    abstract public Constraint constraintInt(Grid g);
    abstract public boolean equals(Object x);
}

