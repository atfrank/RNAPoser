# load needed modules
from sklearn.externals import joblib
import sys, decimal
import numpy as np


# function to read individual decision trees in a forest
def TraverseTreeWeights(tree, node_id=0):
    if (tree.feature[node_id]<0):
        #print " %s %s "%(node_id, tree.value[node_id])
        #s = " "+str(np.argmax(tree.value[node_id]))+" "
        s = " " + str(tree.value[node_id][0][0]) + "-" + str(tree.value[node_id][0][1]) + " "
        sys.stdout.write(s)
        #print " %s %s "%(node_id,tree.value[node_id][0][0])
    else:
        s = " " + str(round(tree.threshold[node_id],10))+":"+str(tree.feature[node_id])+ " "
        sys.stdout.write(s)
        TraverseTreeWeights(tree, tree.children_left[node_id])
        TraverseTreeWeights(tree, tree.children_right[node_id])

# load randomForest predictor
reg = joblib.load("predictor.pkl")

# loop over tree in forest and extract each using TraverseTree
for tree in range(len(reg.estimators_)):
    TraverseTreeWeights(reg.estimators_[tree].tree_)
    sys.stdout.write("\n")
