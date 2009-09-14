import networkx as NX
import time
import sys

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
	prLocs = {} #dictionary to hold probabilities
	#g = Graph()
        g = NX.Graph()
        loc2int = {}

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
			
		prLocs[toks[0]] = pr
		#g.addVertex(toks[0])
                loc2int[toks[0]] = i
                #g.add_node(i, loc=toks[0], prob=pr)
	#	
	#now parse the nodes in the graph
	#gEdges = {}


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

		#gEdges[(toks[0],toks[1])] = secs
		
		#g.addVertex(toks[0])
		#g.addVertex(toks[1])
		#g.addEdge(toks[0],toks[1], secs)
		#g.addEdge(toks[1],toks[0], secs)
                g.add_edge(loc2int[toks[0]],loc2int[toks[1]], weight=secs)
                #g.add_edge(toks[1],toks[0], weight=secs)
	
	file.close()

	return (startLoc, prLocs, g, prLocs.keys(), loc2int)#gEdges)

#Kaditz's Algorithm
#0.) Use Floyd-Warshall to generate all-pairs shortest paths
#1.) sort paths from shortest to longest
#2.) choose the minimum set of paths such that all nodes in graph are in a path
#3.) determine the longest path without doubling back by appending paths together

#0.)  Use Floyd-Warshall to generate all-pairs shortest paths
#0.5) create a list L of all nodes in the graph
#1.)  Pick the shortest path P between 2 nodes that incorporates the most nodes
#     in (not necessarily the longest path) and add it to the current walk W
#1.5) remove all nodes in P from L
#2.)  choose in L (aka not on the current walk) and find the shortest path 
#     between a node on the current walk and that node
#3.)  add this branch to the current walk (need to account for the walk out and back, should walk the longest branch as part of 'main path' so we don't have to traverse it more than once!)
#     Note that this provably can't create a cycle otherwise this node
#     would have already been on the path (assume no negative weights)
#4.)  goto 2. until all nodes have been used.
#
#5.) figure out what the longest path is in the tree
def kaditz(G, startNode, locs, loc2int):
    print locs
    dists, paths = NX.floyd_warshall(G)    
    print dists
    print paths
    
    L = []
    for n in G:
	    if n != loc2int[startNode]:
		    L.append(n)		    

    locs.remove(startNode)    
    nextLoc = locs[0]
    locs = locs[1:]

    path = NX.dijkstra_path(G, loc2int[startNode], loc2int[nextLoc])
    path_len = NX.dijkstra_path_length(G, loc2int[startNode], loc2int[nextLoc])

    
    x = NX.dijkstra_path(G, loc2int[startNode], loc2int['under_bed'])
    l = NX.dijkstra_path_length(G, loc2int[startNode], loc2int['under_bed'])
    print x,l

#
#
#    
def main():

    (startLoc, prLocs, G, locs, loc2int) = loadfile(sys.argv[1])
    print len(G)
    for n,nbrsdict in G.adjacency_iter(): 
        for nbr,eattr in nbrsdict.iteritems():
            if 'weight' in eattr:
                print (n,nbr,eattr['weight'])


    kaditz(G, startLoc, locs, loc2int)


if __name__ == '__main__':
    main()
