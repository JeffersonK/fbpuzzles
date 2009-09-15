from pprint import *
import sys
import path as pathalgs
import graph as nx
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
	g = Graph()

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
		g.addVertex(toks[0])
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
		g.addEdge(toks[0],toks[1], secs)
		g.addEdge(toks[1],toks[0], secs)

	
	file.close()

	return (startLoc, prLocs, g)#gEdges)

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
def calcMinWalk(nodeset, edgeset):
    return    


#return (minWalkTime, maxWalkTime)
def minWalk(currentNode, totWalkTime):
    return


#
#
#
def main():
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

    if not is_connected(G):
        #not garaunteed to find her
        print "-1.00"
        return

    nodeset = G.nodes()

    SpanningTrees = generate_all_spanning_trees(G)
    pprint(SpanningTrees)
    return

    minTree = None
    minWalkTime = 1e6 #INF
    minWalk = None
    for edgeset in SpanningTrees:
        (walk, walkTime) = calcMinWalk(nodeset, edgeset)
        if walkTime < minWalkTime:
            minWalk = walk
            minTree = tree
    

    
        
    return #END MAIN



if __name__ == "__main__":
    main()
