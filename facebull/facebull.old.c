
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define DEBUG 0
#define DEBUG_ASSERT 0
#define DEBUG_LOADFILE 0
#define DEBUG_ADJMATRIX 0
#define DEBUG_DIJKSTRA 0
#define DEBUG_MACHINELIST 0
#define DEBUG_PROGRESS 0

#define INFINITY (HUGE_VAL)

#define MAX_NODES 256
#define MAX_EDGES 2048
#define MAX_NAME_LEN 8

#if DEBUG_ASSERT
#include <assert.h>
#endif

typedef struct {
  short machIndex;
  float edgecost;
  float minpathcost;
  short crossed;
  short predecessor;
} machine_t;

short machineIndex[MAX_EDGES];
short compoundIndex[MAX_NODES];
machine_t adjMatrix[MAX_NODES][MAX_NODES];
float best; //lowest cost
char visited[MAX_NODES];
short startNode; //nodes we started from
short nNodes = 0, nEdges = 0;
short bestMachineList[MAX_EDGES];
short curMachineList[MAX_EDGES];

void printAdjMatrix(void) {
  
  int i, j;
  for (i=0; i < nNodes; i++) {
    for (j=0; j < nNodes; j++) {
      printf("%0.2lf(%0.2f) ", 
	     adjMatrix[i][j].edgecost,
	     adjMatrix[i][j].minpathcost);
    }
    printf("\n");
  }
  
}

float getMarginalPathCost(int u, int v)
{
  if(u==v) return 0.0;
 
  if (adjMatrix[adjMatrix[u][v].predecessor][v].crossed == 0)
    return adjMatrix[adjMatrix[u][v].predecessor][v].edgecost + getMarginalPathCost(u, adjMatrix[u][v].predecessor);
  
  return getMarginalPathCost(u, adjMatrix[u][v].predecessor);  
}

void incCrossedEdges(int u, int v){

  if (u == v) return;

#if DEBUG_ASSERT
  assert(u>=0 && u<nNodes);
  assert(v>=0 && v<nNodes);
  assert(adjMatrix[u][v].predecessor>=0 && adjMatrix[u][v].predecessor<nNodes);
#endif

  adjMatrix[adjMatrix[u][v].predecessor][v].crossed++;
  curMachineList[adjMatrix[adjMatrix[u][v].predecessor][v].machIndex]++;
  //printf("+%d ", adjMatrix[adjMatrix[u][v].predecessor][v].name);
  incCrossedEdges(u, adjMatrix[u][v].predecessor);
}

void decCrossedEdges(int u, int v){

  if (u == v) return;

#if DEBUG_ASSERT
  assert(u>=0 && u<nNodes);
  assert(v>=0 && v<nNodes);
  assert(adjMatrix[u][v].predecessor>=0 && adjMatrix[u][v].predecessor<nNodes);
#endif

  adjMatrix[adjMatrix[u][v].predecessor][v].crossed--;
  curMachineList[adjMatrix[adjMatrix[u][v].predecessor][v].machIndex]--;
  //printf("-%d ", adjMatrix[adjMatrix[u][v].predecessor][v].name);
  decCrossedEdges(u, adjMatrix[u][v].predecessor);
}

void printpath(int x,int i)
{
  printf("\n");
  if(i==x)
    {
      printf("%d",x);
    }
  else if(adjMatrix[x][i].predecessor==-1)
    printf("no path from %d to %d",x,i);
  else
    {
      printpath(x,adjMatrix[x][i].predecessor);
      printf("..%d",i);
    }
}

int minimum(machine_t *m,int k)
{
  float mi=INFINITY;
  int i,t=-1;
  for(i=0;i<k;i++)
    {
      if(m[i].crossed!=1)
	{
	  if(mi>=m[i].minpathcost)
	    {
	      mi=m[i].minpathcost;
	      t=i;
	    }
	}
    }

#if DEBUG_ASSERT
  assert(t >= 0);
#endif

  return t;
}

void dijkstra_single_source_shortest_paths(int source){
  int i, count, u;
  for(i=0;i<nNodes;i++)
    {
      adjMatrix[source][i].crossed=0;
      //adjMatrix[source][i].minpathcost should already be initialized
      adjMatrix[source][i].predecessor=-1;
    }
  adjMatrix[source][source].minpathcost = 0;
  
  count = 0;
  while(count<nNodes)
    {
      u=minimum((machine_t *)&adjMatrix[source],nNodes);
      count++;
      adjMatrix[source][u].crossed=1;
#if DEBUG_ASSERT
      assert(source >= 0 && source < nNodes);
      assert(u>=0 && u < nNodes);
#endif
      for(i=0;i<nNodes;i++)
	{
	  //printf("u=%d\ti=%d\tgraph[u][i]=%d\tmark[i]=%d\n", u, i, graph[u][i], mark[i]);
	  if(adjMatrix[u][i].edgecost > 0)
	    {
	      if(adjMatrix[source][i].crossed!=1)
		{
		  //printf("pathestimate[i]=%d\tpathestime[u]=%d\n", pathestimate[i], pathestimate[u]);
		  if(adjMatrix[source][i].minpathcost > (adjMatrix[source][u].minpathcost + adjMatrix[u][i].edgecost))
		    {
		   
		      //printf(".");
		      adjMatrix[source][i].minpathcost = adjMatrix[source][u].minpathcost + adjMatrix[u][i].edgecost;
		      adjMatrix[source][i].predecessor=u;
		    }
		}
	    }
	}
    }

  return;
}

/*
short name2index(short name, short * list, int len){

  int i;

  for (i = 0; i < len; i++)
    if (name == list[i])
      return i;

  return -1;
  }*/

#define MEMORY_OVERFLOW -666
short insert2index(short name, short * list, short * listlen, int maxlen){
  int i;

  for(i=0; i < *listlen; i++){
    if(list[i] == name)
      return i; //already exists
  }

  if(*listlen < maxlen){
    list[*listlen] = name;
    *listlen += 1;
    return *listlen-1;
  } 
  return MEMORY_OVERFLOW;
}

int loadFile(char * filename) {

  FILE * f = NULL;
  char compoundA[MAX_NAME_LEN];
  char compoundB[MAX_NAME_LEN];
  char machineName[MAX_NAME_LEN];
  int i, j, u, uindex, v, vindex, machIntName, machIndex;
  int cost;
  
  best = INFINITY;

  memset(machineIndex, 0xff, sizeof(machineIndex));
  memset(compoundIndex, 0xff, sizeof(compoundIndex));

  /* initialize the matrix */
  for (i=0; i<MAX_NODES; i++)
    for (j=0; j<MAX_NODES; j++){
      if (i == j) {
	adjMatrix[i][j].edgecost = adjMatrix[i][j].minpathcost = 0;
	
      } else {
	adjMatrix[i][j].edgecost = INFINITY;
	adjMatrix[i][j].minpathcost = INFINITY;
      }
      adjMatrix[i][j].crossed = 0;
      adjMatrix[i][j].machIndex = -1;//the index in the machineList where the real name is stored
  }
  
  f = fopen(filename, "r");
  if (f == NULL)
    return -1;
  
  nEdges = 0;
  while( fscanf(f, "%s %s %s %d\n",
		machineName,
		compoundA,
		compoundB,
		&cost) == 4){

    //this assumes that if the first machine has a name of 0
    //there will be compounds with name 0
    //in the specification of the problem we don't have to handle
    //this but it is useful for testing
    machIntName = atoi(machineName+1);
    //if(nEdges == 0)
    //  machIdShift = machId;
    
    if((machIndex=insert2index(machIntName, machineIndex, &nEdges, MAX_EDGES)) < 0){
      return MEMORY_OVERFLOW;//we ran out of memory allocated
    }
    if (machIndex+1 != nEdges)
      return -5;//this means we had a duplicate machine name

    u = atoi(compoundA+1);
    if( (uindex=insert2index(u, compoundIndex, &nNodes, MAX_NODES)) == MEMORY_OVERFLOW){
      return MEMORY_OVERFLOW;//ran out of memory allocated
    }
    v = atoi(compoundB+1);
    if((vindex=insert2index(v, compoundIndex, &nNodes, MAX_NODES)) == MEMORY_OVERFLOW){
      return MEMORY_OVERFLOW;//ran out of memory allocated
    }

#if DEBUG_LOADFILE
    printf("add (M(%d=>%d) C(%d=>%d) C(%d=>%d) %d\n", machIntName, machIndex, u, uindex, v, vindex, cost);
#endif

    //if (u > max) max = u;
    //if (v > max) max = v;
    adjMatrix[u][v].machIndex = machIndex;
    adjMatrix[u][v].edgecost = (float) cost;
    //nEdges ++;
  }
  //nNodes = max+1;//the highest compound we see is the max

  fclose(f);
  return 0;
}

void printMachineList(void){
  int i, j, firstPrint=0;
  int curMinIndex = -1;
  //for now we do a bubble sort because O(n^2) nothing compared to heart of algorithm
  //could also make the outer loop smaller if we kept track of how many bits were set in the bestMachineList
  for(i=0;i<nEdges;i++){
    curMinIndex = -1;
    for(j=0;j<nEdges;j++){
      if (bestMachineList[j]){
	if(curMinIndex == -1){
	  curMinIndex = j;
	} 
	else if(machineIndex[curMinIndex] > machineIndex[j])
	  curMinIndex = j;	
      }
    }
    if(curMinIndex == -1)
      break;
    if(firstPrint++ == 0)
      printf("%d", machineIndex[curMinIndex]);
    else
      printf(" %d", machineIndex[curMinIndex]);
    bestMachineList[curMinIndex] = 0; //unset this bit so we don't show this machine again
  }
  printf("\n");
}

void printResults(void){
  printf("%d\n", (int) best);
  printMachineList();
}
  

//end condition is all nodes have been visited at least once and we
//are back to the startNode  
void walk(int node, float totCost) {
  int i, end = 1;
  float dCost = 0;
  visited[node] = 1; 
  for(i=0;i<nNodes;i++){
    if(!visited[i]) {
      end = 0;
#if DEBUG_ASSERT
      //shouldnt be possible otherwise graph is not strongly connected to begin with
      assert(adjMatrix[node][i].minpathcost != INFINITY);
#endif
      dCost = getMarginalPathCost(node, i);
      if ((totCost + dCost) < best){
	incCrossedEdges(node, i);
	walk(i, totCost + dCost);
	decCrossedEdges(node, i);
      }
    }
  }
  visited[node] = 0;
  if (end) {
    //now we need to append the path back to the start node
    //remember not to count edges we have already crossed
    totCost += getMarginalPathCost(node, startNode);
    if (((totCost < best) && (best != -1)) || (best < 0)){
      best = totCost;
      incCrossedEdges(node, startNode);
      memmove(bestMachineList, curMachineList, sizeof(short)*nEdges);
      //bzero((void *)curMachineList, sizeof(int)*nEdges);
#if DEBUG_MACHINELIST
      for(i=0;i<nNodes;i++)
	for(end=0;end<nNodes;end++){
	  if(adjMatrix[i][end].machIndex != -1 &&
	     adjMatrix[i][end].crossed > 0){
	    //bestMachineList[adjMatrix[i][end].machIndex] = 1;
	    printf("%d ", adjMatrix[i][end].machIndex);
	  } 
	}
      printf("\n");
#endif
      decCrossedEdges(node, startNode);

    }//update best 

  }//(if(end)
  return;
}//walk()


int main(int argc, char ** argv) {
  int ret, i, j; 
  ret = loadFile(argv[1]);
  
  bzero(curMachineList, sizeof(short)*nEdges);

#if DEBUG_LOADFILE
  printf("loadFile: %d\n", ret);
  printf("nNodes: %d\n", nNodes);
  printf("nEdges: %d\n", nEdges);
#endif

#if DEBUG_ADJMATRIX
  printAdjMatrix();
#endif

  for(i=0;i<nNodes;i++) 
    dijkstra_single_source_shortest_paths(i);

#if DEBUG_DIJKSTRA
  //print the paths
  for(i=0;i<nNodes;i++) 
    for(j=0;j<nNodes;j++) {
      //source, i, predecessor
      if (i != j) {
	printpath(i,j);
	if(adjMatrix[i][j].minpathcost != INFINITY)
	  printf("->(%lf)\n",adjMatrix[i][j].minpathcost);
      }
    }
#endif

  /* Check for any unreachable locations */
  for (i = 0; i < nNodes; i++) {
    for (j = 0; j < nNodes; j++) {
      adjMatrix[i][j].crossed = 0;//re-inialize this field
    }
  }

#if DEBUG_ADJMATRIX
  printf("\n-----\n");
  printAdjMatrix();
#endif

  for(i=0;i<nNodes;i++){    
    //reset visited list
    for(j = 0;j<nNodes; j++)     
      visited[j] = 0;

#if DEBUG_PROGRESS
    printf("Beginning exploration with node %d of %d (%0.2f%%)\n", i+1, nNodes, 100*((float)i/nNodes));
#endif
    startNode = i;
    walk(i, 0);
    
#if DEBUG_ASSERT
    int k,z;
    for(k=0;k<nNodes;k++)
      for(z=0;z<nNodes;z++)
	assert(adjMatrix[k][z].crossed == 0);
#endif
  }

  printResults();
  
  /*
  FILE * file = fopen("out","w+");
  printf("%0.0f\n", best);
  fprintf(file, "%d\n", (int) best);
  for(i=0; i < nNodes;i++){
    if(bestMachineList[i])
      fprintf(file, "%d ", i);
  }
  fprintf(file,"\n");
  
  fclose(file);
  */
  return 0;
}
