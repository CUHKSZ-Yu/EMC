from curses.ascii import isdigit
import imp
from operator import le
import numpy as np
from numpy import append
from ete3 import Tree

import random
import six
from ete3 import *
from tree_set import *
from tree_compare import *
from tree_distance import *


# to compute the RF distance matrix for a set of complete trees
def complete_matrix(treesets):
    for i in range(len(treesets)):
        treesets[i] = treesets[i]+";" 
    
    dim  = len(treesets)
    COMPLETE_rf = [[0]*dim for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            result = tree_compare(treesets[i],treesets[j])
            COMPLETE_rf[i][j] = result['rf']

    return COMPLETE_rf


# to approximate the RF distance matrix for a set of incomplete trees
def incomplete_matrix(trees):
    dim = len(trees)
    INCOMPLETE_rf = [[0]*dim for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            result = incomplete_tree_compare(trees[i],trees[j])
            INCOMPLETE_rf[i][j] = result['rf']
    
    return INCOMPLETE_rf