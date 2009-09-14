import networkx as nx
import matplotlib.pyplot as plt
from pprint import *
import sys

#Construct list L of all Nodes to visit in Graph
#Construct Node Tree T will be a graph with no cycles
#remove startLoc from L
#u = startLoc
#select v in L such that the path p(u,v) is a minimum path
# 3 Cases to consider
#case 1: p(u,v) is an edge in G
#case 2: p(u,v) goes through nodes already in T (not including v)
#case 3: p(u,v) includes nodes not in T (not including v)
def minpath(u, v, minpaths):
    pprint(minpaths)
    print "minpath(%s,%s)" % (str(u),str(v))
    k = v
    p = [u]
    while k != u:
        k = minpaths[u][k]
        if k == u:
            p.append(v)
            break
        else:
            p.append(k)
    return p
    
def build_path_tree(G, startLoc):
    L = []
    
    T = nx.Graph()
    for n in G:
        L.append(n)
        
    #not optimal because we don't have to generate
    #them all always but conveinet to do it first now
    mindists, minpaths = nx.floyd_warshall(G)
    pprint(mindists)
    
    print "L: %s" % L
    print "startLoc: %s" % str(startLoc)
    T.add_node(startLoc)
    L.remove(startLoc)
    
    u = startLoc
    #while (len(L) > 0):
    if 1:
        print "u: %s" % str(u)
        #neighbors = G[u]
        #print neighbors
        for v in L:#neighbors:
            if v in L:
                #we haven't added this node to T yet
                print "v: %s" % str(v)
                #compute the minpath (u,v)
                minp = minpath(u,v, minpaths)
                print "minpath(u,v): %s" % str(minp)
                mind = mindists[u][v]
                for node in T:
                    if mindists[node][v] < mind:
                        minp = minpath(node,v, minpaths)
                    
                #now check to see if any of the nodes/edges
                #in minp excluding u are already in T
                #need to add the path to the tree
                #starting with the first edge in minp
                #that is not in T
                
                if minp == None:
                    #simple case where we just add the edge (u,v)
                    T.add_node(v)
                    T.add_edge(u,v)
                    #u = v
                else:
                    last = u
                    for node in minp[1:]:
                    #Q: What gaurentees that we don't create a cycle when we add this path?
                    #A: because we are constructing a Tree of minpaths if a node is already in it
                    #   it means there is a shorter path to it, which is a contradiction
                        if node not in G:
                            #add node to G
                            T.add_node(node)
                            #add the edge u,node
                            T.add_edge(u,node)
                        last = node
                    #u = last
                #L.remove(v)
    nx.draw(T)
    plt.savefig("T.png")
    plt.show()
    return
                
def find_optimal_walk(G):
    W = nx.Graph()
    return W

def kaditz_path(G, startLoc):
    T = build_path_tree(G, startLoc)

    W = find_optimal_walk(G)
    return

    
        
def main():

    G = nx.Graph()

    G.add_edges_from([(1,2,{'weight':7}),
                      (2,1,{'weight':7}),
                      (1,3,{'weight':14}),
                      (3,1,{'weight':14}),
                      (1,4,{'weight':9}),
                      (4,1,{'weight':9}),
                      (2,4,{'weight':10}),
                      (4,2,{'weight':10}),
                      (2,5,{'weight':15}),
                      (5,2,{'weight':15}),                  
                      (3,4,{'weight':2}),
                      (4,3,{'weight':2}),
                      (3,6,{'weight':9}),
                      (6,3,{'weight':9}),
                      (4,5,{'weight':11}),
                      (5,4,{'weight':11}),
                      (6,5,{'weight':6}),
                      (5,6,{'weight':6}),])
    
    
    #nx.draw(G)
    #plt.savefig("exp.png")
    #plt.show()
    #dists, paths = nx.floyd_warshall(G)
    #pprint(dists)
    #pprint(paths)
    #print
    
    #path = nx.single_source_dijkstra(G, 1)
    #pprint(path)

    kaditz_path(G, 1)
    return











if __name__ == "__main__":
    main()
