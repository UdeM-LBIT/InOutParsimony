from sowing.repr.newick import parse, write
from sowing.traversal import depth
from sowing.zipper import Zipper
from sowing.node import Node
from math import inf
from immutables import Map
import json

def LCA_Content(syntenyTree, x_tilde):
    A = {}
    for nodeZipper in depth(syntenyTree):
        if nodeZipper.is_leaf():
            A[nodeZipper] = x_tilde[nodeZipper.node.data["name"]]
        else: 
            A[nodeZipper] = A[nodeZipper.down(0)] | A[nodeZipper.down(1)]
    B = {}
    B[Zipper(syntenyTree)] = A[Zipper(syntenyTree)]
    x_LCA = {}
    for nodeZipper in depth(syntenyTree, preorder = True):
        if nodeZipper.is_leaf():
            x_LCA[nodeZipper] = B[nodeZipper]
        else:
            B[nodeZipper.down(0)] = B[nodeZipper] - A[nodeZipper.down(1)]
            B[nodeZipper.down(1)] = B[nodeZipper] - A[nodeZipper.down(0)]
            x_LCA[nodeZipper] = (B[nodeZipper] - B[nodeZipper.down(0)]) - B[nodeZipper.down(1)]
    return x_LCA

def Min_Content(syntenyTree, x_tilde, x_LCA):
    x_min = {}
    for nodeZipper in depth(syntenyTree):
        if nodeZipper.is_leaf():
            x_min[nodeZipper] = x_tilde[nodeZipper.node.data["name"]]
        else:
            x_min[nodeZipper] = (x_min[nodeZipper.down(0)] - x_LCA[nodeZipper.down(0)]) | (x_min[nodeZipper.down(1)] - x_LCA[nodeZipper.down(1)])
    return x_min

def loss(parent,child,x_min):
    if x_min[parent] <= x_min[child]:
        return 0
    else:
        return 1

def gain(parent,child,x_min):
    if x_min[child] <= x_min[parent]:
        return 0
    else:
        return 1

def delta_min(parent, child, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra):
    #(Cost,dictChild,(loss_c,gain_c))
    g = gain(parent,child,x_min)
    l = loss(parent,child,x_min)
    case1 = (c_min[child][0] + c_loss*l + c_gain*g, c_min, (l,g))
    case2 = (c_in_extra[child][0] + c_loss*l + c_gain, c_in_extra, (l,1))
    case3 = (c_out_extra[child][0] + c_gain*g + (0 if l else inf), c_out_extra, (0,g))
    case4 = (c_inout_extra[child][0] + c_gain + (0 if l else inf), c_inout_extra, (0,1))
    return min(case1,case2,case3,case4, key = lambda x: x[0])

def delta_out(parent, child, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra):
    #(Cost,dictChild,(loss_c,gain_c))
    g = gain(parent,child,x_min)
    case1 = (c_min[child][0] + c_loss + c_gain*g,c_min,(1,g))
    case2 = (c_in_extra[child][0] + c_loss + c_gain,c_in_extra,(1,1))
    case3 = (c_out_extra[child][0] + c_gain*g,c_out_extra,(0,g))
    case4 = (c_inout_extra[child][0] + c_gain,c_inout_extra,(0,1))
    return min(case1,case2,case3,case4, key = lambda x: x[0])

def delta_in(parent, child, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra):
    #(Cost,dictChild,(loss_c,gain_c))
    g = gain(parent,child,x_min)
    l = loss(parent,child,x_min)
    case1 = (c_min[child][0] + c_loss*l + (0 if g else inf), c_min, (l,0))
    case2 = (c_in_extra[child][0] + c_loss*l, c_in_extra, (l,0))
    case3 = (c_out_extra[child][0] + (0 if (g and l) else inf), c_out_extra, (0,0))
    case4 = (c_inout_extra[child][0] + (0 if l else inf), c_inout_extra, (0,0))
    return min(case1,case2,case3,case4, key = lambda x: x[0])

def delta_inout(parent, child, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra):
    #(Cost,dictChild,(loss_c,gain_c))
    g = gain(parent,child,x_min)
    case1 = (c_min[child][0] + c_loss + (0 if g else inf), c_min, (1,0))
    case2 = (c_in_extra[child][0] + c_loss, c_in_extra, (1,0))
    case3 = (c_out_extra[child][0] + (0 if g else inf), c_out_extra, (0,0))
    case4 = (c_inout_extra[child][0], c_inout_extra, (0,0))
    return min(case1,case2,case3,case4, key = lambda x: x[0])

def pGainLoss(positionGainLoss,nodeZipper,minSol):
    if not nodeZipper.is_leaf():
        positionGainLoss[nodeZipper.down(0)] = minSol[2]
        positionGainLoss[nodeZipper.down(1)] = minSol[4]
        pGainLoss(positionGainLoss, nodeZipper.down(0),minSol[1][nodeZipper.down(0)])
        pGainLoss(positionGainLoss, nodeZipper.down(1),minSol[3][nodeZipper.down(1)])

def x_content(SyntenyTree,x_min,positionGainLoss):
    x = {}
    
    for nodeZipper in depth(SyntenyTree):
        if nodeZipper.is_leaf():
            x[nodeZipper] = x_min[nodeZipper]
        else:
            nodeZipperLeft = nodeZipper.down(0)
            nodeZipperRight = nodeZipper.down(1)
            x_l = x[nodeZipperLeft] if (not positionGainLoss[nodeZipperLeft][1]) else set()
            x_r = x[nodeZipperRight] if (not positionGainLoss[nodeZipperRight][1]) else set()
            x[nodeZipper] = x_min[nodeZipper] | x_l | x_r
    
    for nodeZipper in depth(SyntenyTree, preorder = True):
        if (not positionGainLoss[nodeZipper][0]) and (not nodeZipper.is_root()):
            x[nodeZipper] = x[nodeZipper] | x[nodeZipper.up()] 
                    
    return x

def InOutParsimony(syntenyTree, x_tilde, c_loss, c_gain):
    #Input
    for key in x_tilde:
        x_tilde[key] = set(x_tilde[key])
    x_LCA = LCA_Content(syntenyTree, x_tilde)
    x_min = Min_Content(syntenyTree, x_tilde, x_LCA)
    
    #Cost
    #(Cost,dictLeftChild,(loss_l,gain_l),dictRightChild,(loss_r,gain_r))
    c_min = {}
    c_in_extra = {}
    c_out_extra = {}
    c_inout_extra = {}
    for nodeZipper in depth(syntenyTree):
        if nodeZipper.is_leaf():
            c_min[nodeZipper] = (0,None,None,None,None)
            c_in_extra[nodeZipper] = (inf,None,None,None,None)
            c_out_extra[nodeZipper] = (inf,None,None,None,None)
            c_inout_extra[nodeZipper] = (inf,None,None,None,None)
        else:
            
            leftChild = nodeZipper.down(0)
            
            rightChild = nodeZipper.down(1)
            
            #c_min
            #(Cost,dictChild,(loss_c,gain_c))

            c_min_left = delta_min(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            c_min_right = delta_min(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            
            if x_min[nodeZipper] != set():
                c_min[nodeZipper] = (c_min_left[0] + c_min_right[0], c_min_left[1], c_min_left[2], c_min_right[1], c_min_right[2])
            else:
                c_min[nodeZipper] = inf
            #c_in_extra
            
            #Gains from left only
            c_in_extra_LO_l = delta_in(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_in_extra_LO_r = delta_out(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_in_extra_LO = (c_in_extra_LO_l[0] + c_in_extra_LO_r[0], c_in_extra_LO_l[1], c_in_extra_LO_l[2],
                            c_in_extra_LO_r[1],c_in_extra_LO_r[2])
            
            #Gains from right only
            c_in_extra_RO_r = delta_in(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_in_extra_RO_l = delta_out(nodeZipper, leftChild, x_min, c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_in_extra_RO = (c_in_extra_RO_l[0] + c_in_extra_RO_r[0], c_in_extra_RO_l[1], c_in_extra_RO_l[2],
                            c_in_extra_RO_r[1],c_in_extra_RO_r[2])
            
            #Gains from both sides
            c_in_extra_BS_l = delta_inout(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)

            c_in_extra_BS_r = delta_inout(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_in_extra_BS = (c_in_extra_BS_l[0] + c_in_extra_BS_r[0],c_in_extra_BS_l[1],c_in_extra_BS_l[2],
                             c_in_extra_BS_r[1],c_in_extra_BS_r[2])
            
    
            c_in_extra[nodeZipper] = min(c_in_extra_LO,c_in_extra_RO,c_in_extra_BS, key = lambda x: x[0])
            
            #c_out_extra
            #(Cost,dictChild,(loss_c,gain_c))
            c_out_extra_left = delta_out(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            c_out_extra_right = delta_out(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_out_extra[nodeZipper] = (c_out_extra_left[0] + c_out_extra_right[0],c_out_extra_left[1],
                                       c_out_extra_left[2],c_out_extra_right[1],c_out_extra_right[2])
            #c_inout_extra
            
            #Gains from left only
            c_inout_extra_LO_l = delta_inout(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_LO_r = delta_out(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_LO = (c_inout_extra_LO_l[0] + c_inout_extra_LO_r[0], c_inout_extra_LO_l[1], c_inout_extra_LO_l[2],
                               c_inout_extra_LO_r[1],c_inout_extra_LO_r[2])
            
            #Gains from right only
            c_inout_extra_RO_r = delta_inout(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_RO_l = delta_out(nodeZipper, leftChild, x_min, c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_RO = (c_inout_extra_RO_l[0] + c_inout_extra_RO_r[0], c_inout_extra_RO_l[1], c_inout_extra_RO_l[2],
                               c_inout_extra_RO_r[1],c_inout_extra_RO_r[2])
            
            #Gains from both sides
            c_inout_extra_BS_l = delta_inout(nodeZipper, leftChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_BS_r = delta_inout(nodeZipper, rightChild, x_min,c_loss,c_gain,c_min,c_in_extra,c_out_extra,c_inout_extra)
            
            c_inout_extra_BS = (c_inout_extra_BS_l[0] + c_inout_extra_BS_r[0], c_inout_extra_BS_l[1], c_inout_extra_BS_l[2],
                               c_inout_extra_BS_r[1], c_inout_extra_BS_r[2])
            
            c_inout_extra[nodeZipper] = min(c_inout_extra_LO,c_inout_extra_RO,c_inout_extra_BS, key = lambda x: x[0])
    
    rootZipper = Zipper(syntenyTree) 

    minSol = min(c_min[rootZipper], c_in_extra[rootZipper], key = lambda x: x[0])

    cost = minSol[0] + c_gain

    #Events (loss,gain)
    positionGainLoss = {}
    positionGainLoss[rootZipper] = (0,1)
    pGainLoss(positionGainLoss, rootZipper, minSol)
    
    
    #Contents
    
    x = x_content(syntenyTree, x_min, positionGainLoss)
    
    #Print solution
    numberLosses = 0
    numberGains = 0
    for key in positionGainLoss:
        if positionGainLoss[key][0]:
            numberLosses = numberLosses + 1
        if positionGainLoss[key][1]:
            numberGains = numberGains + 1
    print("Cost: " + str(cost))
    print("Number of losses: " + str(numberLosses))
    print("Number of gains: " + str(numberGains))
    
    
    return (cost, x, positionGainLoss)

def combine_tree(content, left, right):
    return Node(content).add(left).add(right)

def content_To_String(content):
    s = sorted(content)
    s = str(s)
    s = "{" + s[1:-1] + "}"
    return s

#To print solution
def SolutionTreeToPrint(SyntenyTree, x, positionGainLoss):
    solutionTree = {}
    for nodeZipper in depth(SyntenyTree):
        
        if nodeZipper.is_leaf():
            content = nodeZipper.node.data["name"] + " : " + content_To_String(x[nodeZipper])
            solutionTree[nodeZipper] = Node(content)
        else:
            content = content_To_String(x[nodeZipper])
            solutionTree[nodeZipper] = combine_tree(content,solutionTree[nodeZipper.down(0)], solutionTree[nodeZipper.down(1)])
        
        if positionGainLoss[nodeZipper][1]:
            contentUp = set()
            if not nodeZipper.is_root():
                contentUp = x[nodeZipper.up()]
            gain = x[nodeZipper] - contentUp
            content = "Gain : " + content_To_String(gain)
            solutionTree[nodeZipper] = Node(content).add(solutionTree[nodeZipper])                   
        
        if positionGainLoss[nodeZipper][0]:
            loss = x[nodeZipper.up()] - x[nodeZipper]
            content = "Loss : " +  content_To_String(loss)
            solutionTree[nodeZipper] = Node(content).add(solutionTree[nodeZipper])   
    
    return solutionTree[Zipper(SyntenyTree)]

#To write solution
def SolutionTreeToWrite(SyntenyTree, x, positionGainLoss):
    solutionTree = {}
    for nodeZipper in depth(SyntenyTree):
        
        if nodeZipper.is_leaf():
            content = nodeZipper.node.data["name"] +   " : " + content_To_String(x[nodeZipper])
            content = Map({"name":content})
            solutionTree[nodeZipper] = Node(content)
        else:
            content = content_To_String(x[nodeZipper])
            content = Map({"name":content})
            solutionTree[nodeZipper] = combine_tree(content,solutionTree[nodeZipper.down(0)], solutionTree[nodeZipper.down(1)])
        
        if positionGainLoss[nodeZipper][1]:
            contentUp = set()
            if not nodeZipper.is_root():
                contentUp = x[nodeZipper.up()]
            gain = x[nodeZipper] - contentUp
            content = "Gain : " + content_To_String(gain)
            content = Map({"name":content})
            solutionTree[nodeZipper] = Node(content).add(solutionTree[nodeZipper])                   
        
        if positionGainLoss[nodeZipper][0]:
            loss = x[nodeZipper.up()] - x[nodeZipper]
            content = "Loss : " +  content_To_String(loss)
            content = Map({"name":content})
            solutionTree[nodeZipper] = Node(content).add(solutionTree[nodeZipper])   
    
    return write(solutionTree[Zipper(SyntenyTree)])