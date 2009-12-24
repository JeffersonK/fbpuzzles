#include <stdlib.h>
#include <stdio.h>
#include <string.h> //for bzero()
#include <math.h> //for INFINITY

#define DEBUG_LOADFILE 0
//#define DEBUG_MEGBB 1
//#define WITH_PRUNING 1
//#define DEBUG_PRUNING 1

#define INFINITY (HUGE_VAL)
//#define INF ((2^31)-1)
#define MAX_NODES (512)
#define MAX_EDGES (MAX_NODES*MAX_NODES) //number of elements in adjacency matrix (we can do better than this with a linked lists
#define MAX_NAME_LEN (8)

typedef struct MACHINE_T {
  float edgecost; 
  short edgeRealName;
  short crossed;
  short inAbar; //whether edge is part of optimal solution or not 
  short ci;
  short cj;
  struct MACHINE_T * adjNext;
  struct MACHINE_T * adjPrev;	
  struct MACHINE_T * kthNext;//pointer to next highest edge weight
  struct MACHINE_T * nameNext;//pointer to next highest edge weight
} machine_t;

typedef struct COMPOUND_T {
  short mark;
  short predecessor;
  float minpathcost;  
  machine_t * m;
  //struct COMPOUND_T * predecessor;
} compound_t;
compound_t A[MAX_NODES][MAX_NODES];

#define CALCK(i,j) (i*nNodes + j)
#define DONE (99)
#define NOT_DONE (0)
#define MEMORY_OVERFLOW (-666)
float best;
int maxOutDegree;
float graphWeight;
short startNode;
short uniqueVisited;
short visited[MAX_NODES];
machine_t * EdgesAdjHead = NULL;//replaces the adjaceny matrix
machine_t * EdgesNameHead = NULL;//edges sorted by real name (machine name)
machine_t * EdgesCostHead = NULL;//edges sorted in decreasing weight
short compoundIndex[MAX_NODES];
int nNodes = 0;
int nEdges = 0;


/*************************
 *
 * DEBUG FUNCTIONS
 *
 ************************/

/****
 *
 ****/   
void printMachine(machine_t * m){
  fprintf(stderr, "Machine:%hi (@%p)\n", m->edgeRealName, m);
  fprintf(stderr, "\tedgecost:%.0f\n", m->edgecost);
  fprintf(stderr, "\trow:%hi\n", m->ci);
  fprintf(stderr, "\tcol: %hi\n", m->cj);
  fprintf(stderr, "\tadjPrev: @%p\n", m->adjPrev);
  fprintf(stderr, "\tadjNext: @%p\n", m->adjNext);
  fprintf(stderr, "\tkthNext: @%p\n", m->kthNext);
  fprintf(stderr, "\tnameNext: @%p\n", m->nameNext);
}

/****
 *
 ****/   
void printEdgeNames(void){
  machine_t * m = EdgesNameHead;
  fprintf(stderr, " *** Edge Names (Asc Order) ***\n");
  while(m != NULL){
    printMachine(m);
    m = m->nameNext;
  }
}

/****
 *
 ****/   
void printEdgeWeights(void){
  machine_t * m = EdgesCostHead;
  fprintf(stderr, " *** Edge Weights (Desc Order) ***\n");
  while(m != NULL){
    printMachine(m);
    m = m->kthNext;
  }
}

/****
 *
 ****/   
void printEdgeAdj(void){
  machine_t * ptr;
  ptr = EdgesAdjHead;
  fprintf(stderr, " *** Edge Adj List ***\n");
  fprintf(stderr, "\t*** FORWARDS ***\n");
  while(ptr->adjNext != NULL){
    printMachine(ptr);
    ptr = ptr->adjNext;
  }
  printMachine(ptr);

  fprintf(stderr, "\t*** REVERSE ***\n");
  while(ptr->adjPrev != NULL){
    printMachine(ptr);
    ptr = ptr->adjPrev;
  }
  printMachine(ptr);
  fprintf(stderr, "ptr:%p EdgesAdjHead:%p\n", ptr, EdgesAdjHead);
}

/*************************
 *
 * HELPER FUNCTIONS
 *
 ************************/

/****
 *
 ****/   
void freeEdges(machine_t * m){
  if (m == NULL)
    return;  
  //free everyone after you
  freeEdges(m->adjNext);
  //free yourself
  free(m);
  return;
}

/****
 *
 ****/   
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

/****
 *
 * insert from highest to lowest
 *
 ****/   
int insertEdgeCost(machine_t * m){
  machine_t * ptr, * last;
  last = NULL;
  m->kthNext = NULL;  
  if (EdgesCostHead == NULL){
    EdgesCostHead = m;
    return (0);
  }
  ptr = EdgesCostHead;
  while (ptr != NULL) 
    {
      if (ptr->edgecost < m->edgecost){
	//insert to left
	m->kthNext = ptr;
	if (ptr == EdgesCostHead)
	  EdgesCostHead = m;
	else
	  last->kthNext = m;
	return (0);
      } 
      last = ptr;
      ptr = ptr->kthNext;
    }
  last->kthNext = m;
  return (0);
}


/****
 *
 * insert from lowest to highest
 * TODO: combine with insertEdgeCost()
 ****/   
int insertEdgeName(machine_t * m){
  machine_t * ptr, * last;
  last = NULL;
  m->nameNext = NULL;  
  if (EdgesNameHead == NULL){
    EdgesNameHead = m;
    return (0);
  }
  ptr = EdgesNameHead;
  while (ptr != NULL) 
    {
      if (ptr->edgeRealName > m->edgeRealName){
	//insert to left
	m->nameNext = ptr;
	if (ptr == EdgesNameHead)
	  EdgesNameHead = m;
	else
	  last->nameNext = m;
	return (0);
      } 
      last = ptr;
      ptr = ptr->nameNext;
    }
  last->nameNext = m;
  return (0);
}

/***
 *
 ***/  
int insertEdgeAdj(machine_t * m) {
  machine_t * ptr;
  m->adjNext = m->adjPrev = NULL;	
  if (EdgesAdjHead == NULL)
    {
      EdgesAdjHead = m;
      return (0);
    }
  ptr = EdgesAdjHead;	
  while( (ptr->adjNext != NULL) && (CALCK(ptr->ci,ptr->cj) < CALCK(m->ci,m->cj)) )
    ptr = ptr->adjNext;

  if(CALCK(ptr->ci, ptr->cj) ==  CALCK(m->ci,m->cj)){
    fprintf(stderr, "Duplicate Edge.\n");
    return -1;
  }

  if(CALCK(ptr->ci, ptr->cj) > CALCK(m->ci,m->cj)){
    //insert before ptr
    m->adjNext = ptr;
    m->adjPrev = ptr->adjPrev;//NULL
    ptr->adjPrev = m;
    if(ptr == EdgesAdjHead && ptr->adjPrev == NULL){
      EdgesAdjHead = m;
    }  else
      m->adjPrev->adjNext = m;
    if(ptr == EdgesAdjHead && ptr->adjPrev != NULL){
      fprintf(stderr, "invariant broken %d\n", __LINE__);
      return -1;
    }	
    return (0);
  }
  //insert after ptr
  m->adjPrev = ptr;
  m->adjNext = ptr->adjNext;
  if (ptr->adjNext != NULL)
    ptr->adjNext->adjPrev = m;
  ptr->adjNext = m;
  return (0);
}

/****
 *
 *
 ****/
void printAnswer(void){
  machine_t * p = EdgesNameHead;
  int first = 1;
  //printf("%d\n", graphWeight - sumEbarprime);  
  printf("%.0f\n", best);  
  while (p != NULL){
    if (p->inAbar){
      if (first){
	first = 0;		
	printf("%hi", p->edgeRealName);		
      } else
	printf(" %hi", p->edgeRealName);		
    }
    p = p->nameNext;
  }
  printf("\n");
}

/****
 *
 *
 ****/
int loadFile(char * filename) {

  FILE * f = NULL;
  char compoundA[MAX_NAME_LEN];
  char compoundB[MAX_NAME_LEN];
  char machineName[MAX_NAME_LEN];
  int u, uindex, v, vindex;
  short machIntName;
  int cost;
  machine_t * ptr; 

  for(u=0; u<MAX_NODES; u++)
    for(v=0; v<MAX_NODES; v++)
      {
	/*
	if (u == v)
	  A[u][v].edgecost = 0;
	else
	  A[u][v].edgecost = INFINITY;
	*/
	if (u == v)
	  A[u][v].minpathcost = INFINITY;
	else
	  A[u][v].minpathcost = 0;
	A[u][v].predecessor = -1;
	A[u][v].mark = 0;
	A[u][v].m = NULL;
	//A[u][v].edgeRealNmae = -1;
	//A[u][v].inAbar = 0;

      }
  graphWeight = 0;
  memset(compoundIndex, 0xff, sizeof(compoundIndex));

  f = fopen(filename, "r");
  if (f == NULL)
    return -1;
  
  nEdges = 0;
  while( fscanf(f, "%s %s %s %d\n",
		machineName,
		compoundA,
		compoundB,
		&cost) == 4){

    machIntName = (short) atoi(machineName+1);
    u = atoi(compoundA+1);
    if( (uindex=insert2index(u, compoundIndex, &nNodes, MAX_NODES)) == MEMORY_OVERFLOW){
      return MEMORY_OVERFLOW;//ran out of memory allocated
    }

    v = atoi(compoundB+1);
    if((vindex=insert2index(v, compoundIndex, &nNodes, MAX_NODES)) == MEMORY_OVERFLOW){
      return MEMORY_OVERFLOW;//ran out of memory allocated
    }

#if DEBUG_LOADFILE
    fprintf(stderr, "add (M(%d) C(%d=>%d) C(%d=>%d) %d\n", machIntName, u, uindex, v, vindex, cost);
#endif

    if (nNodes >= MAX_NODES)
      return MEMORY_OVERFLOW;

    if (nEdges >= MAX_EDGES)
      return MEMORY_OVERFLOW;

    ptr = (machine_t*) malloc(sizeof(machine_t));  
    if (ptr == NULL)
      return MEMORY_OVERFLOW;
    memset((void *) ptr, 0x0, sizeof(machine_t));
    ptr->crossed = 0;
    ptr->inAbar = 1;
    ptr->edgeRealName = machIntName;
    ptr->ci = uindex;
    ptr->cj = vindex;
    ptr->edgecost = (float)cost;
    graphWeight += (int) cost;
    A[uindex][vindex].m = ptr;

    insertEdgeAdj(ptr);
    insertEdgeName(ptr);
    insertEdgeCost(ptr);
    nEdges ++;
  }	
  best = graphWeight;
  fclose(f);
  return 0;
}

/*****
 *
 *
 *****/
void updateOptimumSoln(float ebar){
  machine_t * p = EdgesAdjHead;
  best = ebar;
  while (p != NULL){
    if (p->crossed)
      p->inAbar = 1;
    else
      p->inAbar = 0;    
    p = p->adjNext;
  }
}

/************
 *
 *
 *
 ************/
short visitedLastCross[MAX_EDGES];
void walk2(int node, float totCost){
  int j;
  visited[node] ++;
  if (visited[node] == 1)
    uniqueVisited ++;

  //base case
  if ( (node == startNode) && (uniqueVisited == nNodes)){
    //this is a possible solution
    if (totCost < best)
      updateOptimumSoln(totCost);
    return;
  }

  //first try to go to nodes you haven't visisted yet
  for(j=0;j<nNodes;j++){
    if( (A[node][j].m != NULL) && !visited[j]){
      if( (A[node][j].m->crossed == 0) && ((totCost + A[node][j].m->edgecost) < best)){
	A[node][j].m->crossed ++;
	visitedLastCross[CALCK(node, j)] = uniqueVisited;
	walk2(j, totCost + A[node][j].m->edgecost);
	A[node][j].m->crossed --;
      }
      else if ( (A[node][j].m->crossed > 0) && (A[node][j].m->crossed < maxOutDegree) ){
	if(visitedLastCross[CALCK(node,j)] != uniqueVisited){
	  A[node][j].m->crossed ++;
	  visitedLastCross[CALCK(node, j)] = uniqueVisited;
	  walk2(j, totCost);
	  A[node][j].m->crossed --;
	}
      }
    }
  }

  //next try revisiting nodes that have already been visited
  for(j=0; j<nNodes; j++){
    if ( (A[node][j].m != NULL) && visited[j]){      
      if( (A[node][j].m->crossed == 0) && ((totCost + A[node][j].m->edgecost) < best)){
	A[node][j].m->crossed ++;
	visitedLastCross[CALCK(node, j)] = uniqueVisited;
	walk2(j, totCost + A[node][j].m->edgecost);
	A[node][j].m->crossed --;
      }
      else if ( (A[node][j].m->crossed > 0) && (A[node][j].m->crossed < maxOutDegree) ){
	if(visitedLastCross[CALCK(node,j)] != uniqueVisited){
	  A[node][j].m->crossed ++;
	  visitedLastCross[CALCK(node, j)] = uniqueVisited;
	  walk2(j, totCost);
	  A[node][j].m->crossed --;
	}
      }
    }
  }
  visited[node] --;
  if(visited[node] == 0)
    uniqueVisited --;
  return;
}


/*************************
 *
 *
 *
 *************************/

/*********************************
 *
 * Main()
 *
 *********************************/
int main(int argc, char ** argv){
  int ret, i, j;
  
  if( (ret=loadFile(argv[1])) < 0){
    fprintf(stderr, "input file format error '%s' (ret=%d).\n", argv[1], ret);
    return -1;
  }
#if DEBUG_LOADFILE
  fprintf(stderr, "loadFile: %d\n", ret);
  fprintf(stderr, "nNodes: %d\n", nNodes);
  fprintf(stderr, "nEdges: %d\n", nEdges);
  fprintf(stderr, "graphWeight: %d\n", graphWeight);
  printEdgeAdj();
  printEdgeNames();
  printEdgeWeights();//desc order
#endif

  if (nEdges == nNodes){
    //must be hamiltonian cycle if it is strongly connected
    //and problem statement said that we could assume it was
    printAnswer();
    goto done;
  }

  maxOutDegree = 0;
  ret = 0;
  for (i=0; i<nNodes; i++)
    {
      for(j=0; j<nNodes; j++)
	if (A[i][j].m != NULL)
	  ret ++;
      
      if (ret > maxOutDegree)
	maxOutDegree = ret;
    }

  //if the ratio of edges to nodes is < 2n use megbb
  //for(i=0; i<nNodes; i++)
  //  dijkstra_single_source_shortest_paths(i);
  //printMinPaths();

  uniqueVisited = 0;
  //reset visited list
  for(j = 0;j<nNodes; j++)     
    visited[j] = 0;
  
  for(j = 0;j<nEdges; j++)     
    visitedLastCross[j] = 0;
  
  startNode = 0;
  walk2(startNode, 0);
  //walk(startNode, 0);
  
  printAnswer();

 done:
  //cleanup for good form
  freeEdges(EdgesAdjHead);
  EdgesAdjHead = EdgesNameHead = EdgesCostHead = NULL;
  
  return (0);
}






//other algorithm

//try to insert nodes in between edges (u,v) in current cycle, look for path (u,k,v), k can't exist, 
//because, that would imply there is a short path between u and v
