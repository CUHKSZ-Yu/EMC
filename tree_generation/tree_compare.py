from curses.ascii import isdigit
from operator import le
import numpy as np
from numpy import append
from ete3 import Tree
import random
import six

from tree_set import *
from tree_distance import *

# to calculate the RF distance between two complete trees
def tree_compare(t1, t2):
    t_1 = Tree(t1)
    t_2 = Tree(t2)
    try:
        result = robinson_foulds(t_1, t_2)
    except:
        result['rf'] = len(t_1) + len(t_2)
    return result


# to approximate the RF distance between two incomplete trees
def incomplete_tree_compare(t1, t2):
    leave_total = Tree(t1+";").get_leaves().__len__()
    try:
        t1_pruned, t2_pruned = prune(t1,t2)
        t1_pruned+=";"
        t2_pruned+=";"
    except:
        result = {}
        result['rf'] = 2*(leave_total)
        return result

    t_1 = Tree(t1_pruned)
    t_2 = Tree(t2_pruned)
    leave_share = t_1.get_leaves().__len__()
    try:
        result = robinson_foulds(t_1, t_2)
        result['rf'] = result['rf'] * (leave_total/leave_share)
    except:
        result = {}
        result['rf'] = (len(t_1) + len(t_2)) * (leave_total/leave_share)
    return result


# to restrict two incomplete trees on common known leaf positions
def prune(t1, t2):
    non_num = ""
    num_1 = ""
    num_2 = ""
    shape = []
    leaf_posi = []
    t1 = str(t1)
    t2 = str(t2)
    missing_pos = []
    t1_leaves = []
    t2_leaves = []
    leave_shape = []
    last_quote = ""
    
    for i in range(len(t1)):
        if isdigit(t1[i]) or t1[i]=="?" :
            num_1 += str(t1[i])
            # num_2 += str(t2[i])
            if t1[i] == "?":
                missing_pos.append(len(leaf_posi))
            if non_num!= "":
                shape.append(non_num)
                non_num = ""      
        else:#not number
            if num_1 != "":
                leaf_posi.append(len(shape))#record the leaf position
                shape.append(num_1);
                t1_leaves.append(num_1);
                # t2_leaves.append(num_2)
                num_1 = ""
                # num_2 = ""
                leave_shape.append(0);
            non_num += t1[i]
            if t1[i] == "(":
                last_quote = "("
            elif t1[i] == ")":
                if last_quote == "(":
                    leave_shape[len(leave_shape)-1]+=2
                    leave_shape[len(leave_shape)-2]+=1
                last_quote = ")"
                
    if non_num!="":
        shape.append(non_num) # the last non empty must be non-numerical string
        non_num = ""  
    
    ct = 0

    for i in range(len(t2)):
        if isdigit(t2[i]) or t2[i]=="?" :
            num_2 += str(t2[i])
            if t2[i] == "?":
                missing_pos.append(ct)
            if non_num!= "":
                non_num = ""
        
        else:# not number
            if num_2 != "":
                t2_leaves.append(num_2)
                num_2 = ""
                ct += 1
            non_num += t2[i]
    missing_pos = set(missing_pos)
    missing_pos = list(missing_pos)  
    missing_pos.sort()
    tt1, tt2, pruned_shape = resemble_trees(t1_leaves, t2_leaves, missing_pos, leave_shape)
    t1_pruned = ""
    t2_pruned = ""
    bitree = False
    bitree_sec = False
    frnum = 0
    flnum = 0
    rnum = 0
    lnum = 0
    pairsnum = 0
    if len(pruned_shape)==1:
        t1_pruned = "("+tt1[0]+")"
        t2_pruned = "("+tt2[0]+")"
        return t1_pruned, t2_pruned
    elif len(pruned_shape) == 0:
        return
    for i in range(len(pruned_shape)):
        if pruned_shape[i] == 1:
            pairsnum+=1
            if bitree == False:
                t1_pruned = t1_pruned + "(" + tt1[i]
                t2_pruned = t2_pruned + "(" + tt2[i]
                bitree = True
            else:
                if pruned_shape[i-1] !=2:
                    t1_pruned = "("+t1_pruned + ", (" + tt1[i]
                    t2_pruned = "("+t2_pruned + ", (" + tt2[i]
                    bitree_sec = True
                    frnum = rnum
                    flnum = lnum
                else:
                    t1_pruned = t1_pruned + ", (" + tt1[i]
                    t2_pruned = t2_pruned + ", (" + tt2[i]
                    bitree_sec = True

        elif pruned_shape[i] == 2 and i != len(pruned_shape)-1:
            t1_pruned = t1_pruned[:t1_pruned.rfind("(")]+"(" + t1_pruned[t1_pruned.rfind("("):]+", " + tt1[i] + ")"
            t2_pruned = t2_pruned[:t2_pruned.rfind("(")]+"(" + t2_pruned[t2_pruned.rfind("("):]+", " + tt2[i] + ")"
        elif pruned_shape[i] == 2 and i == len(pruned_shape)-1:
            t1_pruned = t1_pruned+", " + tt1[i] + ")"
            t2_pruned = t2_pruned+", " + tt2[i] + ")"
            if bitree_sec == True:
                t1_pruned+=")"
                t2_pruned+=")"
        else: 
            if bitree == False:
                t1_pruned = t1_pruned + "(" + tt1[i] + ","
                t2_pruned = t2_pruned + "(" + tt2[i] + ","
                rnum+=1
            else:
                lnum+=1
                t1_pruned =  t1_pruned + ", " + tt1[i] + ")"
                t2_pruned =  t2_pruned + ", " + tt2[i] + ")"
        if i == len(pruned_shape)-1:
            if bitree_sec == False: 
                if lnum >1 and rnum >0:
                    t1_pruned = "("*(lnum-1) + t1_pruned + rnum*")"
                    t2_pruned = "("*(lnum-1) + t2_pruned + rnum*")"
                elif rnum >0:
                    t1_pruned =  t1_pruned + rnum*")"
                    t2_pruned =  t2_pruned + rnum*")"
                elif lnum >1 and rnum ==0:
                        t1_pruned = "("*(lnum-pairsnum) + t1_pruned
                        t2_pruned = "("*(lnum-pairsnum) + t2_pruned

            else:
                if lnum>1 and rnum == 0:
                    if lnum - flnum != 1:
                        t1_pruned = "(" + (flnum-1)*"("+t1_pruned[:t1_pruned.find("), (")+3]+"("*(lnum-flnum-2)+t1_pruned[t1_pruned.find("), (")+3:]+")"
                        t2_pruned = "(" + (flnum-1)*"("+t2_pruned[:t2_pruned.find("), (")+3]+"("*(lnum-flnum-2)+t2_pruned[t2_pruned.find("), (")+3:]+")"
                    else:
                        t1_pruned =  (flnum-1)*"("+t1_pruned[:t1_pruned.find("), (")+3]+"("*(lnum-flnum-2)+t1_pruned[t1_pruned.find("), (")+3:]+")"
                        t2_pruned =  (flnum-1)*"("+t2_pruned[:t2_pruned.find("), (")+3]+"("*(lnum-flnum-2)+t2_pruned[t2_pruned.find("), (")+3:]+")"
                #========================prob to solve====================
                elif lnum >1 and rnum >0:
                    t1_pruned = "("*(lnum-flnum-pairsnum-1) + t1_pruned + rnum*")"
                    t2_pruned = "("*(lnum-flnum-pairsnum-1) + t2_pruned + rnum*")"
                elif rnum > 0 :
                    t1_pruned =  t1_pruned + rnum*")"
                    t2_pruned =  t2_pruned + rnum*")"    
                elif lnum == 1 :
                    try:
                        Tree(t1_pruned+";")
                    except:
                        t1_pruned =  t1_pruned + ")"
                        t2_pruned =  t2_pruned + ")"

            if t1_pruned.find(")") == -1:
                t1_pruned = t1_pruned+")"
                t2_pruned = t2_pruned+")"
            try:
                Tree(t1_pruned+";")
            except:
                if t1_pruned[-1].isdigit():
                    t1_pruned= t1_pruned[:t1_pruned.rfind("(")]+t1_pruned[t1_pruned.rfind("(")+1:]+")"
                    t2_pruned= t2_pruned[:t2_pruned.rfind("(")]+t2_pruned[t2_pruned.rfind("(")+1:]+")"
                else:
                    try:
                        t1_pruned =  t1_pruned[:-1]
                        t2_pruned =  t2_pruned[:-1]
                        Tree(t1_pruned+";")
                    except:
                        t1_pruned =  t1_pruned + ")"
                        t2_pruned =  t2_pruned + ")"
            try:
                Tree(t1_pruned+";")
            except:
                t1_pruned =  t1_pruned + ")"
                t2_pruned =  t2_pruned + ")"
            while True:
                try:
                    Tree(t1_pruned+";")
                    return t1_pruned, t2_pruned  
                except:
                    t1_pruned =  t1_pruned[1:]
                    t2_pruned =  t2_pruned[1:]
    return t1_pruned, t2_pruned      


def resemble_trees(t1, t2, mis_ps, leave_shape):
    j = 0
    skip_flag = False
    for i in mis_ps:
        del t1[i-j]
        del t2[i-j]
        if leave_shape[i-j] == 0:
            del leave_shape[i-j]
            skip_flag = False
        elif leave_shape[i-j] == 1:
            if (i+1) in mis_ps:
                skip_flag = True
            elif i-j > 0:
                if leave_shape[i-j-1] == 2:
                    leave_shape[i-j+1] = 0
                elif leave_shape[i-j-1] == 0 :
                    leave_shape[i-j-1] = 1
            elif i-j==0:
                leave_shape[i-j+1] = 0
                skip_flag = False
            del leave_shape[i-j]
        elif leave_shape[i-j] == 2:
            if skip_flag == True:
                del leave_shape[i-j]
                skip_flag = False
            else:
                if i != 1 and i-j != len(leave_shape)-1 :
                    if (i+1) in mis_ps and skip_flag == False:
                        if leave_shape[i-j+1] == 0:
                           leave_shape[i-j-1] = 0
                        elif leave_shape[i-j+1] == 1:
                            if leave_shape[i-j-1] == 1:
                                leave_shape[i-j-1] = 0
                        del leave_shape[i-j]
                    elif (i+1) not in mis_ps and skip_flag == False:#not consecutive
                        leave_shape[i-j-1] = 0
                        if leave_shape[i-j-1] == 1:
                            leave_shape[i-j-1] = 0
                        del  leave_shape[i-j]
                        skip_flag = False

                elif i-j == len(leave_shape)-1 :
                    leave_shape[i-j-1] = 0
                    del leave_shape[i-j]
                else:
                    if (i+1) in mis_ps and skip_flag == False:
                        if leave_shape[i-j+1] == 0:
                            leave_shape[i-j+2] = 2
                        elif leave_shape[i-j+1] == 1:
                            pass
                        del leave_shape[i-j]
                    elif (i+1) not in mis_ps and skip_flag == False:#not consecutive
                       leave_shape[i-j-1] = 0
                       del  leave_shape[i-j]
                       skip_flag = False
                    elif (i+1) not in mis_ps and skip_flag == True:
                        del  leave_shape[i-j]
                        skip_flag = False
        j += 1 

    if 1 not in leave_shape and 2 not in leave_shape:
        leave_shape[len(leave_shape)-1] = 2
        leave_shape[len(leave_shape)-2] = 1  
    if str(leave_shape).count("1")!= str(leave_shape).count("2"):
        if str(leave_shape).count("1") < str(leave_shape).count("2"):
            two = leave_shape.index(2)
            leave_shape[two-1] = 1
        elif leave_shape[len(leave_shape)-1] == 1:
            leave_shape[len(leave_shape)-1] = 0
        else:
            one = leave_shape.index(1)
            leave_shape[one+1] = 2

    return t1, t2, leave_shape
