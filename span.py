from pprint import *
import sys
import path as pathalgs
import graph as nx
import perm
import matplotlib.pyplot as plt
import networkx
#
#
#
def safe_toInt(strInt):
	try:
		n = int(strInt)
		return n
	except:
		return None


#
#
#
def safe_toFloat(strFloat):
	try:
		n = float(strFloat)
		return n
	except:
		return None

#
#
#
def safe_readline(fileObj):
	try:
		l = fileObj.readline()
		l = l.strip()
		return (len(l), l)
	except:
		return (-1, None)

#
#
#
def loadfile(filename):
	file = None
	#prLocs = {} #dictionary to hold probabilities
	g = nx.Graph()

	try:
		file = open(filename, 'r')
	except:		
		print "error opening file"
		sys.exit(-1)
	
	#
	#first get probabilities for each location
	(linelen, line) = safe_readline(file)
	num_search_locs = safe_toInt(line)
	if num_search_locs == None:
			print "File Format Error: Integer Expected"
			sys.exit(-1)
			
	startLoc = None
	for i in range(num_search_locs):
		(linelen, line) = safe_readline(file)
		toks = line.split()
		if len(toks) != 2:
			print "File Format Error: 2 tokens expected"
			sys.exit(-1)

		if i == 0:
			startLoc = toks[0]
			
		pr = safe_toFloat(toks[1])
		
		if pr == None:
			print "File Format Error: Float Expected"
			sys.exit(-1)
			
                g.add_node(toks[0], p=pr)

	#	
	#now parse the nodes in the graph
	(linelen, line) = safe_readline(file)
	num_edges = safe_toInt(line)
	if num_edges == None:
			print "File Format Error: Integer Expected"
			sys.exit(-1)
			
	for i in range(num_edges):
		(linelen, line) = safe_readline(file)
		
		toks = line.split()
		if len(toks) != 3:
			print "File Format Error: 3 tokens expected"
			sys.exit(-1)
			
		secs = safe_toInt(toks[2])
		
		if secs == None:
			print "File Format Error: Integer Expected"
			sys.exit(-1)

		g.add_edge(toks[0],toks[1], {'weight':secs})

	file.close()
	return (startLoc, g)

#if there are no cycles in a graph 
#then there is no edge you can delete
#where the graph is still connected
def is_connected(Graph):
    nodes = Graph.nodes()
    #print nodes
    path = pathalgs.single_source_dijkstra_path(Graph, nodes[0])
    if len(nodes) != len(path.keys()):
        return False

    #print "%s == %s" % (nodes, path.keys())
    return True

#
#
#
def ConnectedAndNoCycle(Graph):

    if not is_connected(Graph):
        return False

    edges = Graph.edges(data=True)
    for edge in edges:
        Graph.remove_edge(edge[0],edge[1])
        if is_connected(Graph):
            #must have a cycle
            return False
        else:
            #put the edge back in
            Graph.add_edge(edge[0],edge[1],edge[2])

    #if we get through all the edges
    #and there is no edge we can delete
    #where the graph stays connected if we delete it
    #there is no cycle
    return True

def generate_all_spanning_trees(G):

    #edgelist = G.edges(data=True)
    #print edgelist
    nodelist = G.nodes()
    C = G.copy()
    items = C.edges(data=True)
    #print items
    l = [[x for (pos,x) in zip(range(len(items)), items) if (2**pos) & switches] for switches in range(2**len(items))]
    
    C = None
    SpanningTrees = []
    for edgeset in l:
        
        #print edgeset
        if len(edgeset) == 0:
            #skip the empty set
            continue
        
        C = nx.Graph()
        C.add_nodes_from(nodelist)
        C.add_edges_from(edgeset)
        if ConnectedAndNoCycle(C):
            #print "Spanning Tree: %s" % C.edges()
            SpanningTrees.append(edgeset)
            
        C = None

    return SpanningTrees


#
#
#
def calcMinWalk(startNode, g):#nodeset, edgeset):

    (minExpWalkTime, tTot, pathAccum) = subTreeMinWalk(g, None, startNode, 0)
    return (minExpWalkTime, tTot, pathAccum)


#g is a graph with no cycles becaues it is a
#spanning tree
#return (minExpWalk, totTimeExlapsed)
def subTreeMinWalk(g, previousNode, currentNode, tAccum):
    print "\tPREVIOUS NODE: %s" % previousNode
    print "\tCURRENT NODE:  %s" % currentNode
    
    neighbors = g.neighbors(currentNode)
    if previousNode != None and previousNode in neighbors:
	    neighbors.remove(previousNode)

    print "\tNEIGHBORS: %s" % neighbors

    #we are at a leaf node
    if len(neighbors) == 0:
        dt = g[previousNode][currentNode]['weight']
        min = (tAccum + dt)*g.node[currentNode]['p']
	print "AT LEAF: return(%s, %s)" % (min, tAccum+2*dt)
        return (min, tAccum+2*dt, [currentNode])

    #there is only one choice of path
    if len(neighbors) == 1:
        dtHere = 0.0
	if previousNode == None:
		dtHere = 0.0
	else:
		dtHere = g[previousNode][currentNode]['weight']

	tHere = tAccum + dtHere
	(minExpWalk, tTot, pathAccum) = subTreeMinWalk(g, currentNode, neighbors[0], tHere)
	    
        #send same values back but account for traversing this edge again
        return (minExpWalk + (tHere*g.node[currentNode]['p']), 
		tHere + 2*dtHere, 
		[currentNode] + pathAccum)

    #generate a list of all possible sequences to visit neighbors.
    #should only generate lists of length len(neighbors)
    neighbors_visit_orders = perm.perm(neighbors)
    
    #print "ENUMERATION OF VISIT ORDERS: %s" % neighbors_visit_orders
    minMinExpWalk = 1e9
    minTimeElapsed = None
    minVisitOrder = []
    minPath = None

    dtHere = 0.0
    if previousNode == None:
	    dtHere = 0.0
    else:
	    dtHere = g[previousNode][currentNode]['weight'] 

    tHere = tAccum + dtHere

    for neighbor_visit_order in neighbors_visit_orders:
        if len(neighbor_visit_order) != len(neighbors):
            print "\tSKIPPING VISIT ORDER: %s" % neighbor_visit_order
            continue

        print "\tEXAMINING VISIT ORDER: %s" % neighbor_visit_order
        minExpWalk_order = 0.0
        totTime_order = tHere
	pathAccum = []
        for n in neighbor_visit_order:
	    (minExpWalk, resultT, path) = subTreeMinWalk(g,currentNode, n, totTime_order) 
	    minExpWalk_order += minExpWalk
	    #we should only add this if its not part of the longest
	    #path through the subtree
            totTime_order = resultT + (resultT - totTime_order)
	    pathAccum += path
	print "\t\tVIST ORDER RESULT: (%s, %s)" % (minExpWalk_order, totTime_order)
        #keep track of the best
        if minExpWalk_order < minMinExpWalk:
            minMinExpWalk = minExpWalk_order
            minTimeElapsed = totTime_order
            minVisitOrder = neighbor_visit_order
	    minPath = pathAccum
    print "\t\tMIN VISIT ORDER: %s" % minVisitOrder
    print "\t\tMIN PATH: %s" % minPath
    #now account for the walk from the previous node to this node
    dExpWalk = tHere*g.node[currentNode]['p']
    #dt = minTimeElapsed - tHere (change in time from traversing this subtree)
    return (dExpWalk + minMinExpWalk, minTimeElapsed + 2*dtHere, [currentNode] + minPath)



def testGraph1():
    edges = [(1,2,{'weight':7}),
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
             (5,6,{'weight':6}),]

    G = nx.Graph()
    G.add_edges_from(edges)
    return G
#
#
#
def main():

    #G = testGraph1()
    startNode, G = loadfile(sys.argv[1])

    pprint(G.nodes(data=True))
    pprint(G.edges(data=True))
    networkx.draw(G)
    
    #pos = {}
    #labels = {}
    #i = 1
    #for n in G.nodes():
    #    pos[n] = (i, i)
    #    labels[n] = G.node[n]['p']
    #    i += 1
    #networkx.draw_networkx_labels(G, pos, labels)

    #plt.savefig("C.png")
    #plt.show()
    #plt.close()

    if not is_connected(G):
        #not garaunteed to find her
        print "-1.00"
        return

    nodeset = G.nodes()

    SpanningTrees = generate_all_spanning_trees(G)


    #Global Minimums
    minTree = None
    minWalkTime = 1e9 #INF
    minExpWalkTime = 1e9
    minPath = None
    minGraph = None
    ith_spanningtree = None
    i=0
    for edgeset in SpanningTrees:
        C = nx.Graph()
        C.add_edges_from(edgeset)
        C.add_nodes_from(nodeset)

        for n in nodeset:
            C.node[n]['p'] = G.node[n]['p']
   
        networkx.draw(C)
        plt.savefig("SpanningTree-%d.png"%i)
        #plt.show()
        plt.close()
        #plt.clr()
	print "\nTRYING SPANNING TREE-%d: %s" % (i, C.edges(data=True))
	(expWalkTime, walkTime, path) = calcMinWalk(startNode, C)#nodeset, edgeset)
        print "MIN EXP WALK: %s" % expWalkTime
	print "WALK TIME:    %s" % walkTime
	print "MIN PATH:     %s" % path
	if expWalkTime < minExpWalkTime:
            minExpWalkTime = expWalkTime
            minGraph = C
            minWalkTime = walkTime
	    minPath = path
	    ith_spanningtree = i
	#help out the python garbage collector
	C = None
	i += 1
    print
    print "*** RESULT ***"
    print "SPANNING TREE (%d)" % ith_spanningtree
    print "MIN EXP WALK TIME: %.2f" % minExpWalkTime
    print "MIN WALK TIME:     %.2f" % minWalkTime
    print "MIN PATH:          %s" % minPath
    print "EDGE SET:          %s" % minGraph.edges(data=True)
    print
        
    return #END MAIN



if __name__ == "__main__":
    main()
