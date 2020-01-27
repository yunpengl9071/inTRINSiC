import os
import sys
import random
import numpy
import math

wkdir = "/Users/yunpengl/Data/multilayerNetwork/data/"
f = open(wkdir + "9606.protein.actions.v11.0.txt", 'r')
edgelist = f.readlines()
f.close()
edgelist = [line.strip() for line in edgelist]
edgelist = edgelist[1:]

# Tally up edges by mode of action

edgemap = {}
actionTypes = {}

for edge in edgelist:
    edgeinfo = edge.split("\t")
    protA = edgeinfo[0].split(".")[1]
    protB = edgeinfo[1].split(".")[1]
    mode = edgeinfo[2]
    isDirectional = edgeinfo[4]
    AisActing = edgeinfo[5]
    score = float(edgeinfo[6])
    edgetup = (protA, protB)
    if edgetup not in edgemap:
        edgemap[edgetup] = [(mode, isDirectional, AisActing, score)]
        if mode not in actionTypes:
            actionTypes[mode] = 1
        else:
            actionTypes[mode] += 1
    elif (mode, isDirectional, AisActing, score) not in edgemap[edgetup]:
        edgemap[edgetup] += [(mode, isDirectional, AisActing, score)]
        if mode not in actionTypes:
            actionTypes[mode] = 1
        else:
            actionTypes[mode] += 1

for at in actionTypes:
    sys.stdout.write(at + '.\n')

# Make priority order
# {activation, inhibition} > {reaction, catalysis, ptmod} > {binding}
# Expression not considered here since txn regulation network is built separately.
# For binding-only interactions, create bidirectional edges.
# Note that some interaction pairs contain both activation and inhibition
# relationships, here we keep them both in the model.

# Version 2: only consider activation/inhibition relationships.
# For duplicate edges containing both activation & inhibition, look for higher score.
priority = {"activation":3, "inhibition":3, "reaction":2, "catalysis":2, "ptmod":2, "binding":1, "expression":0}
sign = {"activation":1,"inhibition":-1}
edgemap_simp = set()

for pair in edgemap:
    modePriority = [priority[tup[0]] for tup in edgemap[pair]]
    if max(modePriority) == 3:
        maxScore = 0
        best_tup = ()
        for tup in edgemap[pair]:
            if priority[tup[0]] == 3:
                if tup[1]=='t':
                    if tup[2]=='t':
                        edge_simp = (pair[0],pair[1],sign[tup[0]]*tup[3])
                    else:
                        edge_simp = (pair[1],pair[0],sign[tup[0]]*tup[3])
                    if edge_simp not in edgemap_simp:
                        edgemap_simp.add(edge_simp)
                else:
                    edge_simp1 = (pair[0],pair[1],sign[tup[0]]*tup[3])
                    edge_simp2 = (pair[1],pair[0],sign[tup[0]]*tup[3])
                    if edge_simp1 not in edgemap_simp:
                        edgemap_simp.add(edge_simp1)
                    if edge_simp2 not in edgemap_simp:
                        edgemap_simp.add(edge_simp2)
#                if tup[3] > maxScore:
#                    maxScore = tup[3]
#                    best_tup = tup
#        tup = best_tup
#        if tup[1]=='t':
#            if tup[2]=='t':
#                edge_simp = (pair[0],pair[1],sign[tup[0]])
#            else:
#                edge_simp = (pair[1],pair[0],sign[tup[0]])
#            if edge_simp not in edgemap_simp:
#                edgemap_simp.add(edge_simp)
#        else:
#            edge_simp1 = (pair[0],pair[1],sign[tup[0]])
#            edge_simp2 = (pair[1],pair[0],sign[tup[0]])
#            if edge_simp1 not in edgemap_simp:
#                edgemap_simp.add(edge_simp1)
#            if edge_simp2 not in edgemap_simp:
#                edgemap_simp.add(edge_simp2)
##    elif max(modePriority) == 2:
#        ind = modePriority.index(max(modePriority))
#        if edgemap[pair][ind][1]=='t':
#            if edgemap[pair][ind][2]=='t':
#                edge_simp = (pair[0],pair[1],0)
#            else:
#                edge_simp = (pair[1],pair[0],0)
#            edgemap_simp.add(edge_simp)
#        else:
#            edge_simp1 = (pair[0],pair[1],0)
#            edge_simp2 = (pair[1],pair[0],0)
#            edgemap_simp.add(edge_simp1)
#            edgemap_simp.add(edge_simp2)
#    elif max(modePriority) == 1:
#        ind = modePriority.index(max(modePriority))
#        edge_simp1 = (pair[0],pair[1],0)
#        edge_simp2 = (pair[1],pair[0],0)
#        edgemap_simp.add(edge_simp1)
#        edgemap_simp.add(edge_simp2)
#rdcy = [len(key) for key in edgemap_simp]

# The edgelist seems to be symmetric, i.e. the same edge may be represented twice
# where A and B are simply swapped and AisActing negated
edgemap_simp = list(edgemap_simp)

sys.stdout.write('\nTotal number of candidate edges: ' + str(len(edgemap_simp)) + '.\n')

# Save complete edge list and filter elsewhere
f = open(wkdir + "STRING_11.0_human_cleaned_edges.act_inh_directional.txt", 'w')
for edge in edgemap_simp:
    line = '\t'.join([edge[0],edge[1],str(edge[2])]) + '\n'
    f.write(line)

f.close()


