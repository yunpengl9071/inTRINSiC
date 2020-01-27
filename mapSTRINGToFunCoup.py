import os
import sys
import random
import numpy
import math

wkdir = "/Users/yunpengl/Data/multilayerNetwork/data/"
f = open(wkdir + "hs.FunCoupNet", 'r')
alledges_FUNCOUP = f.readlines()
alledges_FUNCOUP = alledges_FUNCOUP[1:len(alledges_FUNCOUP)]
f.close()
f = open(wkdir + "geneSym.mapped.STRING.v11.act_inh.edges.txt", 'r')
alledges_STRING = f.readlines()
f.close()
f = open(wkdir + "ENSP_geneSym.txt", 'r')
allIDs = f.readlines()
allIDs = allIDs[1:len(allIDs)]
f.close()


edgemap_FUNCOUP = {}
for line in alledges_FUNCOUP:
    edgeinfo = line.strip().split('\t')
    if edgeinfo[2] not in edgemap_FUNCOUP:
        edgemap_FUNCOUP[edgeinfo[2]] = {edgeinfo[3]:[float(edgeinfo[1])]}
    elif edgeinfo[3] not in edgemap_FUNCOUP[edgeinfo[2]]:
		edgemap_FUNCOUP[edgeinfo[2]][edgeinfo[3]] = [float(edgeinfo[1])]
    else:
        edgemap_FUNCOUP[edgeinfo[2]][edgeinfo[3]] += [float(edgeinfo[1])]


edgemap_STRING = {}
for line in alledges_STRING:
    edgeinfo = line.strip().split('\t')
    if edgeinfo[0] not in edgemap_STRING:
        edgemap_STRING[edgeinfo[0]] = {edgeinfo[1]:[float(edgeinfo[2])]}
    elif edgeinfo[1] not in edgemap_STRING[edgeinfo[0]]:
        edgemap_STRING[edgeinfo[0]][edgeinfo[1]] = [float(edgeinfo[2])]
    else:
        edgemap_STRING[edgeinfo[0]][edgeinfo[1]] += [float(edgeinfo[2])]


IDmap = {}
for line in allIDs:
    IDinfo = line.strip().split(',')
    if IDinfo[0] not in IDmap and IDinfo[2] != '':
		IDmap[IDinfo[0]] = [IDinfo[2]]
    elif IDinfo[2] != '':
        IDmap[IDinfo[0]] += [IDinfo[2]]

revIDmap = {}
for line in allIDs:
    IDinfo = line.strip().split(',')
    if IDinfo[2] not in revIDmap and IDinfo[2] != '':
		revIDmap[IDinfo[2]] = [IDinfo[0]]
    elif IDinfo[2] != '':
        revIDmap[IDinfo[2]] += [IDinfo[0]]


# Now for each edge in the filtered FUNCOUP network, if
# there is a corresponding edge in the STRING network,
# check for the edge sign with the highest evidence (if multiple)
# and assign that sign to the FUNCOUP edge score.
# For the time being multiple edges between the same pair of nodes
# are not allowed.

# 042118 edit: should use STRING network as backbone and overlay
# scores wherever possible, since the FUNCOUP network is essentially
# undirected and may contain indirect interactions. For now try taking
# the intersection of the two networks.

#for node_s in edgemap_FUNCOUP:
#    for node_t in edgemap_FUNCOUP[node_s]:
#        scores = edgemap_FUNCOUP[node_s][node_t]
#        edgemap_FUNCOUP[node_s][node_t] = [max(scores)]
#        if node_s in IDmap and node_t in IDmap:
#            s_ENSP = IDmap[node_s]
#            t_ENSP = IDmap[node_t]
#            st_scores = []
#            for sp in s_ENSP:
#                for tp in t_ENSP:
#                    if sp in edgemap_STRING:
#                        if tp in edgemap_STRING[sp]:
#                            #sys.stdout.write(".")
#                            st_scores += edgemap_STRING[sp][tp]
#            if len(st_scores):
#                st_scores_abs = [abs(x) for x in st_scores]
#                max_abs_score = max(st_scores_abs)
#                max_scoring_sign = [i for i in range(0,len(st_scores_abs)) if st_scores_abs[i]==max_abs_score]
#                if len(max_scoring_sign) == 1:
#                    positive = (st_scores[max_scoring_sign[0]] > 0)
#                    if not positive:
#                        edgemap_FUNCOUP[node_s][node_t] = [-max(scores)]

edgemap_inters = set()

for node_s in edgemap_STRING:
    for node_t in edgemap_STRING[node_s]:
        s_ENSGs = revIDmap[node_s]
        t_ENSGs = revIDmap[node_t]
        scores_STRING = edgemap_STRING[node_s][node_t]
        abs_scores_STRING = [abs(x) for x in scores_STRING]
        max_abs_score = max(abs_scores_STRING)
        sign = scores_STRING[abs_scores_STRING.index(max_abs_score)]
        for s_ENSG in s_ENSGs:
            for t_ENSG in t_ENSGs:
                if s_ENSG in edgemap_FUNCOUP and t_ENSG in edgemap_FUNCOUP[s_ENSG]:
                    scores = edgemap_FUNCOUP[s_ENSG][t_ENSG]
                    score = max(scores)
                    if sign > 0:
                        edgemap_inters.add((s_ENSG,t_ENSG,score))
                    else:
                        edgemap_inters.add((s_ENSG,t_ENSG,-score))
                elif t_ENSG in edgemap_FUNCOUP and s_ENSG in edgemap_FUNCOUP[t_ENSG]:
                    scores = edgemap_FUNCOUP[t_ENSG][s_ENSG]
                    score = max(scores)
                    if sign > 0:
                        edgemap_inters.add((s_ENSG,t_ENSG,score))
                    else:
                        edgemap_inters.add((s_ENSG,t_ENSG,-score))



outfile = "FUNCOUP_STRING_signed_act_inh_inters.txt"
f = open(wkdir + outfile, 'w')
for edge in edgemap_inters:
    line = '\t'.join([edge[0],edge[1],str(edge[2]),'\n'])
    f.write(line)
#f.write('\t'.join(['0:PFC','1:FBS_max','2:Gene1','3:Gene2']) + '\n')
#for node_s in edgemap_FUNCOUP:
#    for node_t in edgemap_FUNCOUP[node_s]:
#        line = '\t'.join([str(edgemap_FUNCOUP[node_s][node_t][0]),str(edgemap_FUNCOUP[node_s][node_t][0]),node_s,node_t]) + '\n'
#        f.write(line)

f.close()

#sys.stdout.write('\nTotal number of edges mapped: ' + str(nmapped) + '.\n')
#sys.stdout.flush()






