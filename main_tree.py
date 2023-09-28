from curses.ascii import isdigit
from operator import le
import numpy as np
import pandas as pd
from numpy import append
from ete3 import Tree
import random
import six
import sys
sys.path.append("tree_generation/")

from tree_set import *
from tree_compare import *
from tree_matrix import *
from tree_distance import *

tree_shape = "(((((((((2, 9), 7), 6), 3), 0), 4), 8), 1), 5)"
Ntree_list = [100, 200, 500]
ratio_list = [0.4, 0.6, 0.8]
print("Start mission!")

Nleave = len(Tree(tree_shape+";"))
print("Nleave = " + str(Nleave))
for i in range(len(Ntree_list)):
    Ntree = Ntree_list[i]
    print("  Dtrue: Ntree = " + str(Ntree))
    # to generate a set of complete trees for a given tree shape
    T_true = tree_complete_set(tree_shape, Ntree)
    # to calculate the RF distance matrix for a set of complete trees
    D_true = complete_matrix(T_true)

    df_rf = pd.DataFrame(list(D_true))
    df_rf.to_csv('Dtrue_L%d_N%d.csv'%(Nleave,Ntree), index = False, header = False)

    for j in range(len(ratio_list)):
        ratio = ratio_list[j]
        print("  Dmiss: missing ratio = " + str(ratio))
        # to generate a set of incomplete trees with a given missing ratio
        T_miss = tree_incomplete_set(T_true, ratio)
        try:
            # to approximate the RF distance matrix for a set of incomplete trees
            D_miss = incomplete_matrix(T_miss)

            df_rf_miss = pd.DataFrame(list(D_miss))
            df_rf_miss.to_csv('Dmiss_L%d_N%d_R%d.csv'%(Nleave,Ntree,10*ratio), index = False, header = False)
        except:
            pass
    

    