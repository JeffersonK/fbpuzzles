
import random
from pprint import *
import networkx as nx
import matplotlib.pyplot as plt
random.seed()

MAXN = 150 #compounds
MAXM = 150 #Edges

MIN_WEIGHT = 1
MAX_WEIGHT = 1000

#
#
#
def readGraphFromFile(filename):

    f = open(filename)
    lines = f.read().split('\n')
    G = nx.DiGraph()
    for line in lines:
        tok = line.split()
        if len(tok) == 4:
            G.add_edges_from([(tok[1][1:],tok[2][1:], {'label':tok[0][1:], 'weight':int(tok[3])})])
    return G

#
#
#
def writeGraphToFile(filename, G):
    
    solWeight = 0
    solList = []

    f = open(filename+".in", "w+")
    for n,nbrs in G.adjacency_iter():
        for nbr,eattr in nbrs.iteritems():
            weight = eattr['weight']
            label = eattr['label']
            inSol = eattr['solution']
            if inSol:
                solWeight += weight
                solList.append(label)
            print (n, nbr, weight, label, inSol)
            f.write("M%d\tC%d\tC%d\t%d\n" % (label, n, nbr, weight))
    f.close()
            
    solList.sort()
    solStr = "%d\n" % solWeight
    solStr += "%d" % solList[0]
    for m in solList[1:]:
        solStr += " %d" % m
    f = open(filename+".out", "w+")
    f.write(solStr+"\n")
    f.close()  

    #nx.draw(G)
    #plt.savefig(filename+'.png')

    return

def createGraphPng(G, filename):
    nx.draw(G)
    plt.savefig(filename+'.png')
    return

#generate a solution that has c cycles in it
def generateRandomSoln(n, m, c=1, minWeight=MIN_WEIGHT, maxWeight=MAX_WEIGHT):
    solutionE = []

    minpath = []#ordered list of vertices that will be the minimum path through the graph
    minpath = range(0,n)
    random.shuffle(minpath)
    print minpath

    for i in range(1,c):
        r = random.randint(0,len(minpath)-1)
        pos = random.randint(0,len(minpath)-1)
        while not ((r != minpath[pos]) and (r != minpath[pos-1]) and (r != minpath[(pos+1)%len(minpath)])):
            r = random.randint(0,len(minpath)-1)
            pos = random.randint(0,len(minpath)-1)
        #print minpath
        #print r, pos
        minpath = minpath[0:pos] + [r] + minpath[pos:]
        
    #generate edges in the path
    for i in range(0,len(minpath)):
        solutionE.append((minpath[i],minpath[(i+1)%len(minpath)]))
    print solutionE
        
    #generate the edges that won't be part of the solution
    nonSolE = []
    for i in range(0,m-len(solutionE)):
        u = random.randint(0,n-1)
        v = random.randint(0,n-1)
        while (u == v) or ((u,v) in set(solutionE + nonSolE)):
            u = random.randint(0,n-1)
            v = random.randint(0,n-1)
        nonSolE.append((u,v))
    #print nonSolE

    #need to generate edge weights
    machineIDs = range(1,m+1)
    if len(machineIDs) != len(nonSolE) + len(solutionE):
        print "invariant broken."

    #need to assign random machine numbers the first 
    random.shuffle(machineIDs)

    #create at networkx graph object with the solution edgeset        
    G = nx.DiGraph()
    for (u,v) in solutionE:
        #G.add_edges_from([(u,v,{'weight':random.randint(minWeight,maxWeight), 'label':machineIDs[0], 'solution':True})])
        G.add_edge(u,v,random.randint(minWeight,maxWeight))
        #machineIDs = machineIDs[1:]

    for (u,v) in nonSolE:
        l = nx.single_source_dijkstra_path_length(G,u)
        print u,l[v]
        #G.add_edges_from([(u,v,{'label':machineIDs[0],'weight':l[v]+100,'solution':False})])
        G.add_edge(u,v,random.randint(l[v],l[v]+10))
        #machineIDs = machineIDs[1:]

    for n,nbrs in G.adjacency_iter():
        for nbr,eattr in nbrs.iteritems():
            G.remove_edge(n, nbr)
            if (n,nbr) in solutionE:
                G.add_edges_from([(n,nbr,{'label':machineIDs[0],'weight':eattr,'solution':True})])
            else:
                G.add_edges_from([(n,nbr,{'label':machineIDs[0],'weight':eattr,'solution':False})])
            machineIDs = machineIDs[1:]
    #save the current set of edge numbers as the solution set E'
    #fill in with random edges that aren't part of the solution to fill out E
    #every edge we generate after solution must be an edge such that Cost(u,v) >= minpath(u,v)            
    #pprint(edgeset)
    return (G)


def generateDenseGraphSoln(n, m):
    edgeset = []
    return edgeset

#generate a solution where there are
#groups of nodes with cheap edges
#connected by expensive edges
def generateClusteredSolnType1(n, m):
    edgeset = []
    return edgeset

#generate a solution where there are
#groups of nodes with cheap edges
#connected by expensive edges
def generateClusteredSolnType2(n, m):
    edgeset = []
    return edgeset
         

def main():
    G = generateRandomSoln(110, 200, 1)
    writeGraphToFile("graph",G)

    #G2 = readGraphFromFile("graph.in")
    #createGraphPng(G2, 'graph') 
    return


if __name__ == "__main__":
    main()
