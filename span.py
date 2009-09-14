from pprint import *
import sys
import networkx as nx


def connected_and_nocycles(currentNode, edgeMatrix, nodesVisitedList):
    connected = False
    if len(nodesVisitedList) == len(edgeMatrix.keys()):
        connected = True

    for i in range(len(edgeMatrix[currentNode])):
        if edgeMatrix[currentNode][i] != None:
            if i not in nodesVisitedList:
                isConnected, isCycle = connected_and_nocycles(i, edgeMatrix, nodesVisitedList + [i])
                if isConnected:
                    connected = True

                if isCycle:
                    #exit early if there is a cycle
                    return (connected, True)
            else:
                #there is a cycle
                return (connected, True)
                

    return (connected, False)

 def generate_all_spanning_trees(G):

    #edgelist = G.edges(data=True)
    #print edgelist
              
    C = G.copy()
    items = C.edges(data=True)
    print items
    l = [[x for (pos,x) in zip(range(len(items)), items) if (2**pos) & switches] for switches in range(2**len(items))]

    for edgeset in l:
        print edgeset

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
    
    generate_all_spanning_trees(G)
    return



if __name__ == "__main__":
    main()
