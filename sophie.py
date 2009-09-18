#!/usr/bin/python

#DEBUG SWITCHES
WITH_MATPLOT = 0
WITH_NETWORKX = 0
DEBUG = 0

#GLOBALS
INF = 1e9

#from pprint import *
import sys
#import path as pathalgs
import graph as nx

if WITH_MATPLOT and WITH_NETWORKX:
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

#
#
#
def perm(items, n=None):
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[:i] + items[i+1:]
            for p in perm(rest, n-1):
                yield v + p

#
#
#
def comb(items, n=None):
    if n is None:
        n = len(items)
    for i in range(len(items)):
        v = items[i:i+1]
        if n == 1:
            yield v
        else:
            rest = items[i+1:]
            for c in comb(rest, n-1):
                yield v + c

import heapq
def single_source_dijkstra(G,source,target=None,cutoff=None ):
    """Returns shortest paths and lengths in a weighted graph G.

    Uses Dijkstra's algorithm for shortest paths. 
    Returns a tuple of two dictionaries keyed by node.
    The first dicdtionary stores distance from the source.
    The second stores the path from the source to that node.

    Parameters
    ----------
    G : NetworkX graph

    source : node label
       Starting node for path

    target : node label, optional
       Ending node for path 

    cutoff : integer or float, optional
        Depth to stop the search - only
        paths of length <= cutoff are returned.


    Examples
    --------
    >>> G=nx.path_graph(5)
    >>> length,path=nx.single_source_dijkstra(G,0)
    >>> print length[4]
    4
    >>> print length
    {0: 0, 1: 1, 2: 2, 3: 3, 4: 4}
    >>> path[4]
    [0, 1, 2, 3, 4]

    Notes
    ---------
    Distances are calculated as sums of weighted edges traversed.
    Edges must hold numerical values for Graph and DiGraphs.

    Based on the Python cookbook recipe (119466) at 
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/119466

    This algorithm is not guaranteed to work if edge weights
    are negative or are floating point numbers
    (overflows and roundoff errors can cause problems). 
    
    See Also
    --------
    single_source_dijkstra_path()
    single_source_dijkstra_path_length()
    
    """
    if source==target: return (0, [source])
    dist = {}  # dictionary of final distances
    paths = {source:[source]}  # dictionary of paths
    seen = {source:0} 
    fringe=[] # use heapq with (distance,label) tuples 
    heapq.heappush(fringe,(0,source))
    while fringe:
        (d,v)=heapq.heappop(fringe)
        if v in dist: continue # already searched this node.
        dist[v] = d
        if v == target: break
        #for ignore,w,edgedata in G.edges_iter(v,data=True):
        #is about 30% slower than the following
        if G.is_multigraph():
            edata=[]
            for w,keydata in G[v].items():
                edata.append((w,
                             {'weight':min((dd.get('weight',1)
                                            for k,dd in keydata.iteritems()))}))
        else:
            edata=G[v].iteritems()


        for w,edgedata in edata:
            vw_dist = dist[v] + edgedata.get('weight',1)
            if cutoff is not None:
                if vw_dist>cutoff: 
                    continue
            if w in dist:
                if vw_dist < dist[w]:
                    raise ValueError,\
                          "Contradictory paths found: negative weights?"
            elif w not in seen or vw_dist < seen[w]:
                seen[w] = vw_dist
                heapq.heappush(fringe,(vw_dist,w))
                paths[w] = paths[v]+[w]
    return (dist,paths)

#if there are no cycles in a graph 
#then there is no edge you can delete
#where the graph is still connected
#can prove by contradiction trivially
def is_connected(Graph):

    nodes = Graph.nodes()
    #path = pathalgs.single_source_dijkstra_path(Graph, nodes[0])
    length, path = single_source_dijkstra(Graph, nodes[0])

    if len(nodes) != len(path.keys()):
        return False

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
    if DEBUG:
	    print "\tPREVIOUS NODE: %s" % previousNode
	    print "\tCURRENT NODE:  %s" % currentNode
    
    neighbors = g.neighbors(currentNode)
    if previousNode != None and previousNode in neighbors:
	    neighbors.remove(previousNode)

    if DEBUG:
	    print "\tNEIGHBORS: %s" % neighbors

    #we are at a leaf node
    if len(neighbors) == 0:
        dt = g[previousNode][currentNode]['weight']
        min = (tAccum + dt)*g.node[currentNode]['p']
	if DEBUG:
		print "\tAT LEAF: return(%s, %s)" % (min, tAccum+2*dt)
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
		tTot + dtHere, 
		[currentNode] + pathAccum)

    #generate a list of all possible sequences to visit neighbors.
    #should only generate lists of length len(neighbors)
    neighbors_visit_orders = perm(neighbors)
    
    minMinExpWalk = INF
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
            #shouldn't hit this case but in case
            if DEBUG:
		    print "\tSKIPPING VISIT ORDER: %s" % neighbor_visit_order
            continue

	if DEBUG:
		print "\tEXAMINING VISIT ORDER: %s" % neighbor_visit_order
        minExpWalk_order = 0.0
        totTime_order = tHere
	pathAccum = []
	i = 1
        for n in neighbor_visit_order:
	    (minExpWalk, resultT, path) = subTreeMinWalk(g,currentNode, n, totTime_order) 
	    minExpWalk_order += minExpWalk
            totTime_order = resultT# + tHere
	    pathAccum += path
	    
	    #account for the fact that we have to come back through
	    #this subtree if its not the last subtree in the optimal
	    #walk of the spanning tree
	    if i < len(neighbor_visit_order):
		    path.reverse()
		    pathAccum += path[1:] + [currentNode]
	    i += 1

	if DEBUG:
		print "\t\tVIST ORDER RESULT: (%s, %s)" % (minExpWalk_order, totTime_order)
        #keep track of the best
        if minExpWalk_order < minMinExpWalk:
            minMinExpWalk = minExpWalk_order
            minTimeElapsed = totTime_order
            minVisitOrder = neighbor_visit_order
	    minPath = pathAccum
	if DEBUG:
	    print "\t\tMIN VISIT ORDER: %s" % minVisitOrder
	    print "\t\tMIN PATH: %s" % minPath
    #now account for the walk from the previous node to this node
    dExpWalk = tHere*g.node[currentNode]['p']
    #dt = minTimeElapsed - tHere (change in time from traversing this subtree)
    return (dExpWalk + minMinExpWalk, minTimeElapsed + dtHere, [currentNode] + minPath)



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

    #pprint(G.nodes(data=True))
    #pprint(G.edges(data=True))
    
    if WITH_NETWORKX:
	    networkx.draw(G)
    
    if not is_connected(G):
        #not gauranteed to find her
        print "-1.00\n"
        return

    nodeset = G.nodes()

    SpanningTrees = generate_all_spanning_trees(G)

    #Vars to track Global Minimums
    minTree = None
    minWalkTime = INF
    minExpWalkTime = INF
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
   
	if WITH_NETWORKX:
		networkx.draw(C)
	
	if WITH_MATPLOT:
		plt.savefig("SpanningTree-%d.png"%i)
		plt.close()

	if DEBUG:
		print "\nTRYING SPANNING TREE-%d: %s" % (i, C.edges(data=True))
	(expWalkTime, walkTime, path) = calcMinWalk(startNode, C)#nodeset, edgeset)

	if DEBUG:
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

    if DEBUG:
	    print
	    print "*** RESULT ***"
	    print "SPANNING TREE (%d)" % ith_spanningtree
	    print "MIN EXP WALK TIME: %.2f" % minExpWalkTime
	    print "MIN WALK TIME:     %.2f" % minWalkTime
	    print "MIN PATH:          %s" % minPath
	    print "EDGE SET:          %s" % minGraph.edges(data=True)
	    print
        
    print "%.2f\n" % minExpWalkTime
    return #END MAIN



if __name__ == "__main__":
    main()
