import sys
import os
from pprint import *
import copy

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

class Graph:
	
	INF = 10e6
	
	#
	#
	#
	def __init__(self, vertexList=[], edgeList=[]):

		#Keep this around for debugging
		self.__vertexList = vertexList
		self.__edgeList = edgeList

		#["vertex id", ....]
		self.__vertices = {}
		for node in vertexList:
			self.__vertices[str(node)] = 1
		
		#[(nodeA,nodeB,weight)
		self.__edges = {}
		for (nodeA,nodeB,weight) in edgeList:
			exists = True
			if nodeA not in self.__vertexList:
				exists = False
				print "vertex %s doesn't exist" % nodeA
			if nodeB not in self.__vertexList:
				exists = False
				print "vertex %s doesn't exist" % nodeB
			
			if type(weight) != type(1):
				exists = False
				print "weight type:%s is not an integer." % type(weight)
			
			if exists:	
				self.__edges[(nodeA,nodeB)] = weight


	#
	#
	#
	def __str__(self):
		#pprint(self.__edges)
		return "{'vertices':%s,'edges':%s}" % (self.__vertices,
						       self.__edges)



	#0. Use Floyd-Warshall to generate all-pairs shortest paths
	#0.5) create a list L of all nodes in the graph
	#1. Pick the shortest path P between 2 nodes that incorporates the most nodes
	#   in (not necessarily the longest path) and add it to the current walk W
	#1.5) remove all nodes in P from L
	#2. choose in L (aka not on the current walk) and find the shortest path 
	#   between a node on the current walk and that node
	#3. add this branch to the current walk (need to account for the walk out and back)
        #   Note that this provably can't create a cycle otherwise this node
	#   would have already been on the path (assume no negative weights)
	#4. goto 2. until all nodes have been used.
	#

	def floyd_warshall(self):
		import networkx as NX
		g = NX.Graph()
		i = 0
		paths = {}
		for v in self.__vertices.keys():
			g.add_node(i,name=v)
			i += 1
			paths[v] = {}
			for y in self.__vertices.keys():
				paths[v][y] = -1
		
		print g



	#
	#
	#
	def numEdges(self):
		return len(self.__edges)

	#
	#
	#
	def numVertices(self):
		return len(self.__vertices)

	#
	#
	#
	def vertexExist(self, vertexID):
		if vertexID in self.__vertices:
			return True
		return False
	#
	#
	#
	def edgeExist(self, a, b):
		if (a,b) in self.__edges:
			return True
		return False
	
	#
	#
	#
	def addVertex(self, vertexID):
		if type(vertexID) != type(""):
			print "vertexID:%s is not a string identifier" % vertexID
			return -1

		if self.vertexExist(vertexID):
			print "vertex: %s already in graph." % vertexID
			return -2

		self.__vertices[vertexID] = 1
		return 0


	#
	#
	#
	def addEdge(self, A, B, weight):
		exists = True
		if A not in self.__vertices:
			exists = False
			print "vertex %s doesn't exist" % A

		if B not in self.__vertices:
			exists = False
			print "vertex %s doesn't exist" % B
			
		if type(weight) != type(1):
			exists = False
			print "weight type:%s is not an integer." % type(weight)
			
		if exists:	
			self.__edges[(A,B)] = weight
			return 0
		
		return -1
	
	#
	#
	#
	def listVertices(self):
		return self.__vertices.keys()

	#
	#
	#
	def listEdges(self):
		return self.__edges.keys()

	#
	#
	#
	def getAdjacent(self, vertexID):
		if not self.vertexExist(vertexID):
			return None
		
		adj = []
		for (A,B) in self.__edges.iterkeys():
			if A == vertexID:
				adj.append((B,self.__edges[(A,B)]))

		return adj
	#
	#
	#
	def minDist(self, A, B, accum, minDist):
		#get L list of adjacent nodes to A
		L = self.getAdjacent(A)
		newMin = accum
		#if B is in L and dist(A,B) < 
		
		#
		return


	#
	#
	#
	def _min(self, d, Q):
		min = self.INF
		minv = None
		for v in d.iterkeys():
			if d[v][0] < min and v in Q:
				min = d[v][0]
				minv = v

		return (minv, min)

	def dijkstra(self, source):
		g = {}
		for v in self.__vertices.iterkeys():
			#(distance, previous)
			if source == v:
				g[v] = (0, None)
			else:
				g[v] = (self.INF, None)
		
		L = copy.deepcopy(g)
		Q = self.listVertices()
		#print Q
		while len(Q) > 0:
			u, dist = self._min(g,Q)
			#print u, dist
			if dist == self.INF:
				break
			#print "Q(%d) = %s" % (len(L),L)
			Q.remove(u)
			neighbors = self.getAdjacent(u)
			#print neighbors
			for v,w in neighbors:
				if v in Q: #hasn't been removed yet
					alt = dist + self.__edges[(u,v)]
					if alt < g[v][0]:
						g[v] = (alt, u)
		
		#print g
		return g

	#
	#
	#
	def bellman_ford(self, source):
		g = {}
		L = self.listVertices
		for v in self.__vertices.iterkeys():
			distance = None
			predecessor = None
			if v == source:
				distance = 0
			else:
				distance = self.INF
			g[v] = (distance, predecessor)

		for i in range(len(self.__vertices)):
			for uv in self.__edges.iterkeys():
				#print uv
				u = uv[0]
				v = uv[1]
				weight = self.__edges[(u,v)]
				if g[u][0] + weight < g[v][0]:
					g[v] = (g[u][0] + weight, g[v][1])
					g[v] = (g[v][0], u)


		#pprint(g)
		return g



	#
	#
	#
	def isconnected(self):
		#just pick forst vertex in list as
		#starting poit if its connected it doesn't
		#matter
		source = self.__vertices.keys()[0]
		g = self.dijkstra(source)
		for (dst,(l,src)) in g.iteritems():
			if l == self.INF:
				#print "%s not reachable" % dst
				return False

		return True

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
###########
#
#
#
###########
def main():
	startLoc, prLocs, graph = loadfile(sys.argv[1])
	
	
	#DEBUG PRINT
	pprint(startLoc)
	pprint(prLocs)
	#pprint(gEdges)
	print graph
	#print graph.numVertices()
	#print graph.numEdges()
	#print graph.getAdjacent("behind_blinds")

	print "BELLMAN-FORD"
	print graph.bellman_ford(startLoc)

	print "DIJKSTRA"
	print graph.dijkstra(startLoc)

	#check to make sure all vertices are reachable
	#aka the graph is connected
	if not graph.isconnected():
		print "-1.00"

	#sys.setrecursionlimit(10000)
	#print graph.allwalks([startLoc], 0)
	graph.floyd_warshall()

if __name__ == "__main__":
	main()





#NP HARD SO JUST STOP
	#currentwalk is a list of nodes we have visited previously to node
        #walklen is the length of the walk thus far 
	def allwalks(self, currentwalk, walklen):
		print "currentwalk: %s" % currentwalk
	#if all nodes are in currentwalk
	#we are done with this walk
	#print currentwalk, walklen
	#return [(currentwalk, walklen)]
		c = set(currentwalk)
		t = set(self.__vertices.keys())
		#print t
		#print c
		x = len(t) - len(c)
		if x == 0:
			print "%s => %d" % (currentwalk, walklen)
			return (currentwalk, walklen)

		else:
	#for each edge attached to this node call
			edges = []
			for (a,b) in self.__edges.iterkeys():
				if currentwalk[-1] == a:
				#store the edge and the weight
					edges += [(a,b,self.__edges[(a,b)])]

			minwalk = []
			minlen = self.INF
			print "edges: %s" % edges
			for (a,b,l) in edges:
			
	#CHECK TO AVOID PING-PONG-ING	
        #if we are at node x and can go to node n and the 
	#last three nodes visited look like x,n,x skip n
				if len(currentwalk) > 2 and \
					    currentwalk[-3] == a and \
					    currentwalk[-2] == b and \
					    currentwalk[-1] == a:
					pass
				else:			
	#else: make the recursive call
	#allwalks(currentwalk+nextvertex, walklen+nextlen,
					(walk, length) = self.allwalks(currentwalk+[b], walklen+l)
					if length < minlen:
						minlen = length
						minwalk = walk

	#make a list of the returned results and return it
			print "minwalk: %s minlen: %d" % (minwalk, minlen)
			return (minwalk, minlen)
	
