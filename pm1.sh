#!/bin/csh

time java -classpath classes/ -Dgraph.h=3 -Dgraph.cyclic=false -Dgraph.sym=true -Dgrid.mfactor=4 -Dgrid.maxlevel=3 gridsearch.ParametrizedMatrix

