from curses.ascii import isdigit
from operator import le
import numpy as np
from numpy import append
from ete3 import Tree
import random
import six

from tree_compare import *
from tree_matrix import *
from tree_distance import *


# Generate complete trees for given structure

def tree_complete_set(tree_shape, N):
    """
    Generate the set of complete trees for a given structure.
    :param tree_shape:         Given tree structure.
    :param N:                  Number of generated unrooted complete trees.
    :return tree_complete_set: A set of complete trees with same tree structure.
    """
    non_num = ""
    num = ""
    shape = []
    leaf_posi = []
    tree_shape = str(tree_shape)
    for i in range(len(tree_shape)):
        if isdigit(tree_shape[i]):
            num += str(tree_shape[i])
            if non_num!= "":
                shape.append(non_num)
                non_num = ""
        else:
            if num != "":
                leaf_posi.append(len(shape)) #record the leaf position
                shape.append(num)
                num = ""
            non_num += tree_shape[i]
    if non_num!="":
        shape.append(non_num) # the last non empty must be non-numerical string
        non_num = ""     

    #=================generate randomized order===============
    total_leaves = list(range(len(leaf_posi)))
    new_sets = []
    len_of_set = 0
    while len(new_sets) < N :
        random.shuffle(total_leaves)
        new_tree = total_leaves[:]
        if new_tree not in new_sets:
            new_sets.append(new_tree)
            len_of_set += 1

    #================generate complete trees=================
    return print_trees(shape, new_sets, leaf_posi)


# Generate incomplete trees for given missing ratio

def tree_incomplete_set(complete_set, missr):
    """
    Generate the set of incomplete trees for a given complete tree set.
    :param complete_set:         A complete tree set.
    :param missr:                The missing ratio varies from 0% to 100%.
    :return tree_incomplete_set: A set of incomplete trees with some missing leaf labels.
    """
    #======================get shape and leaf position=============
    tree_0 = complete_set[0] 
    non_num = ""
    num = ""
    shape = []
    leaf_posi = []
    tree_0 = str(tree_0)
    for i in range(len(tree_0)):
        if isdigit(tree_0[i]):
            num += str(tree_0[i])
            if non_num!= "":
                shape.append(non_num)
                non_num = ""
        else:
            if num != "":
                leaf_posi.append(len(shape)) #record the leaf position
                shape.append(num)
                num = ""
            non_num += tree_0[i]
    if non_num!="":
        shape.append(non_num) # the last non empty must be non-numerical string
        non_num = ""   
    #=================get each tree's distribution========
    total_leaf_distribution = []
    for tree in complete_set:
        numerical =""
        leaf_distri = []
        for i in range(len(tree)):
            if isdigit(tree[i]):
                numerical+=str(tree[i])
            else:
                if numerical!="":
                    leaf_distri.append(int(numerical))
                    numerical = ""
        total_leaf_distribution.append(leaf_distri)
    #=================generate miss===================
    if 0<missr<1:
        n = round(missr*len(leaf_posi))
        for leaf_distri in total_leaf_distribution:
            a = list(range(len(leaf_posi)))
            random.shuffle(a)
            miss_pos = a[:]
            # print(miss_pos)
            for times in range(n):
                leaf_distri[miss_pos[times]] = "?"
    #=================print incomplete trees==========
    return print_trees(shape,total_leaf_distribution,leaf_posi)


def print_trees(shape, leaf_sets, leaf_posi):
    assembled_total = []
    for tree in leaf_sets:
        leaf_num = 0
        assembled_tree = ""
        for i in range(len(shape)):
            if i in leaf_posi:
                assembled_tree += str(tree[leaf_num])
                leaf_num+=1
            else: assembled_tree += shape[i]
        assembled_total.append(assembled_tree)

    return assembled_total  