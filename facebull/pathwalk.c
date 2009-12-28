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
#define CALCK(i,j) (i*nNodes + j)
#define DONE (99)
#define NOT_DONE (0)
#define MEMORY_OVERFLOW (-666)

typedef struct MACHINE_T {
  float edgecost; 
  int edgeRealName;
  int crossed;
  int inAbar; //whether edge is part of optimal solution or not 
  int ci;
  int cj;
  struct MACHINE_T * adjNext;
  struct MACHINE_T * adjPrev;	
  struct MACHINE_T * kthNext;//pointer to next highest edge weight
  struct MACHINE_T * nameNext;//pointer to next highest edge weight
} machine_t;

typedef struct COMPOUND_T {
  int mark;
  int predecessor;
  float minpathcost;  
  machine_t * m;
} compound_t;
compound_t A[MAX_NODES][MAX_NODES];

float best;
float graphWeight;
float theoryMin;

int startNode;
int uniqueVisited;
int visited[MAX_NODES];
machine_t * EdgesAdjHead = NULL;//replaces the adjaceny matrix
machine_t * EdgesNameHead = NULL;//edges sorted by real name (machine name)
machine_t * EdgesCostHead = NULL;//edges sorted in decreasing weight
int compoundIndex[MAX_NODES];
int nNodes = 0;
int nEdges = 0;
machine_t * rowMin[MAX_NODES];

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

void printpath(int x,int i)
{
  printf("\n");
  if(i==x)
    {
      printf("%d",x);
    }
  else if(A[x][i].predecessor==-1)
    printf("no path from %d to %d",x,i);
  else
    {
      printpath(x,A[x][i].predecessor);
      printf("..%d",i);
    }
}


void printMinPaths(void){
  int i, j;
  fprintf(stderr, "*** MIN PATHS ***\n");
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++)
      fprintf(stderr, "%04.0f(%d)\t", A[i][j].minpathcost, A[i][j].predecessor);
    fprintf(stderr, "\n");
  }
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
  //clear him out for good measure, we might get this reallocated to us later
  //bzero(m, sizeof(machine_t));
  //free yourself
  free(m);
  return;
}

/****
 *
 ****/   
int insert2index(int name, int * list, int * listlen, int maxlen){
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
  int machIntName;
  int cost;
  machine_t * ptr; 

  for(u=0; u<MAX_NODES; u++){
    rowMin[u] = NULL;
    for(v=0; v<MAX_NODES; v++)
      {
	if (u == v)
	  A[u][v].minpathcost = INFINITY;
	else
	  A[u][v].minpathcost = 0;
	A[u][v].predecessor = -1;
	A[u][v].mark = 0;
	A[u][v].m = NULL;
      }
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

    machIntName = (int) atoi(machineName+1);
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
    //ptr = &A[uindex][vindex];    
    memset((void *) ptr, 0x0, sizeof(machine_t));
    //ptr->mark = 0;
    //ptr->predecessor = -1;
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

    if(rowMin[uindex] == NULL)
      rowMin[uindex] = ptr;
    else if(rowMin[uindex]->edgecost > ptr->edgecost) 
      rowMin[uindex] = ptr;

  }

  theoryMin = 0;
  for(u=0;u<nNodes;u++)
    theoryMin += rowMin[u]->edgecost;

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
  //sumEbarprime = sumEbar;
  //nEbarprime = nEbar;
  best = ebar;
  while (p != NULL){
    if (p->crossed)//p->edgeset == E_EDGE)
      p->inAbar = 1;
    else
      p->inAbar = 0;    
    p = p->adjNext;
  }
}

/*****
 *
 *
 *****/
int minimum(compound_t *m,int k)
{
  float mi=INFINITY;
  int i,t=-1;
  for(i=0;i<k;i++)
    {
      if(m[i].mark!=1)
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


void dijkstra_single_source_shortest_paths(int source){//, machine_t **A, int nNodes){
  int i, count, u;
  float edgecost = INFINITY;
  for(i=0;i<nNodes;i++)
    {
      A[source][i].mark=0;
      A[source][i].minpathcost = INFINITY;
      A[source][i].predecessor = -1;
    }
  A[source][source].minpathcost = 0;

  count = 0;
  while(count<nNodes)
    {
      u=minimum((compound_t *)&A[source][0],nNodes);
      count++;
      A[source][u].mark=1;
#if DEBUG_ASSERT
      assert(source >= 0 && source < nNodes);
      assert(u>=0 && u < nNodes);
#endif
      for(i=0;i<nNodes;i++)
        {
          //printf("u=%d\ti=%d\tgraph[u][i]=%d\tmark[i]=%d\n", u, i, graph[u][i], mark[i]);
	  if(A[u][i].m == NULL)
	    edgecost = INFINITY;
	  else if(u == i)
	    edgecost = 0;
	  else
	    edgecost = A[u][i].m->edgecost;
	  
          if(/*A[u][i].*/edgecost > 0)
            {
              if(A[source][i].mark!=1)
                {
		  
                  //printf("pathestimate[i]=%d\tpathestime[u]=%d\n", pathestimate[i], pathestimate[u]);
                  if(A[source][i].minpathcost > (A[source][u].minpathcost + edgecost /*A[u][i].edgecost*/))
                    {
                      //printf(".");
                      A[source][i].minpathcost = A[source][u].minpathcost + /*A[u][i].*/edgecost;
                      A[source][i].predecessor=u;
                    }
                }
            }
        }
    }

  return;
}


float getMarginalPathCost(int u, int v, float totCost){
  float accum = 0.0;
  int vtmp = v;
  while(u != vtmp){
     if (A[A[u][vtmp].predecessor][vtmp].m->crossed == 0)
       accum += A[A[u][vtmp].predecessor][vtmp].m->edgecost;      
     vtmp = A[u][vtmp].predecessor;
  }
  return accum;
}
/*
float getMarginalPathCost_r(int u, int v, float accum)
{
  if(u==v) return 0.0;//accum;//0.0;
  
  if (A[A[u][v].predecessor][v].m->crossed == 0){
    return A[A[u][v].predecessor][v].m->edgecost + getMarginalPathCost_r(u, A[u][v].predecessor, accum +  A[A[u][v].predecessor][v].m->edgecost);
  }
  return getMarginalPathCost_r(u, A[u][v].predecessor, accum);  
  }*/
/*
 void incCrossedEdges_r(int u, int v, int dst){

  if (u == v) {
    return;
  }

  A[A[u][v].predecessor][v].m->crossed++;
  //don't increment for last node on path, or we'll end up
  //double counting the visit
  if (v != dst) {
    if (!visited[v])
      uniqueVisited ++;
    visited[v] ++;
  }
  incCrossedEdges_r(u, A[u][v].predecessor, dst);
  }*/
 
void incCrossedEdges(int u, int v, int dst){
  int vtmp = v;
  while(u != vtmp){
    A[A[u][vtmp].predecessor][vtmp].m->crossed++;
    //printf("xing - %d (%d,%d)\n", A[A[u][vtmp].predecessor][vtmp].m->edgeRealName, A[u][vtmp].predecessor, vtmp);
    //don't increment for last node on path, or we'll end up
     //double counting the visit
    if (vtmp != dst) {
      if (!visited[vtmp])
	uniqueVisited ++;
      visited[vtmp] ++;
     
    }
    vtmp = A[u][vtmp].predecessor;
  }
}
/*
void decCrossedEdges_r(int u, int v, int dst){
  if (u == v) {
    return;
  }
  A[A[u][v].predecessor][v].m->crossed--;
  if (v != dst){
  visited[v]--;
  if(!visited[v])
    uniqueVisited --;
  }
  decCrossedEdges_r(u, A[u][v].predecessor, dst);
  }*/

void decCrossedEdges(int u, int v, int dst){
  int vtmp = v;
  while(u != vtmp){
    A[A[u][vtmp].predecessor][vtmp].m->crossed--;
    //printf("unxing - %d (%d,%d)\n", A[A[u][vtmp].predecessor][vtmp].m->edgeRealName, A[u][vtmp].predecessor, vtmp);

     //don't increment for last node on path, or we'll end up
     //double counting the visit
    if (vtmp != dst) {
      visited[vtmp]--;
      if (!visited[vtmp])
	uniqueVisited--;
    }
    vtmp = A[u][vtmp].predecessor;
  }
}

void walk(int node, float totCost) {
  int i,end = 1; 
  float dCost = 0;

  if(!visited[node])
    uniqueVisited++;
  visited[node]++; 
  
  if ( (node == startNode) 
	 && (uniqueVisited == nNodes) ){
    if(totCost < best)
      updateOptimumSoln(totCost);
    return;
  }

  for(i=0;i<nNodes;i++){
    if(node!=i && (uniqueVisited < nNodes) && !visited[i]) {
      end = 0;
      dCost = getMarginalPathCost(node, i, totCost);
      if ((totCost + dCost) < best){
	incCrossedEdges(node, i, i);
	walk(i, totCost + dCost);
	decCrossedEdges(node, i, i);
      }
      //if(best == theoryMin)
      //return;

    }
  }

  if (0 && end) {
    //now we need to append the path back to the start node
    //remember not to count edges we have already crossed
    dCost = getMarginalPathCost(node, startNode, totCost);
    if ((totCost + dCost) < best){
      incCrossedEdges(node, startNode, startNode);
      updateOptimumSoln(totCost + dCost);
      decCrossedEdges(node, startNode, startNode);
    }//update best 
  }//(if(end)

  visited[node]--;
  if(!visited[node])
    uniqueVisited--;

  return;
}//walk()


/*********************************
 *
 * Main()
 *
 *********************************/
int main(int argc, char ** argv){
  int ret, i, j;
  
  if( (ret=loadFile(argv[1])) < 0){
    fprintf(stderr, "input file format error (ret=%d).\n", ret);
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

  for(i=0; i<nNodes; i++)
    dijkstra_single_source_shortest_paths(i);
  //printMinPaths();

  //for(i=0; i<nNodes; i++){
  uniqueVisited = 0;
  //reset visited list
  for(j = 0;j<nNodes; j++)     
    visited[j] = 0;
  
  best = graphWeight;
  startNode = 0;
    walk(startNode, 0);
    printAnswer();
    //}
 done:
  //cleanup for good form
  freeEdges(EdgesAdjHead);
  EdgesAdjHead = EdgesNameHead = EdgesCostHead = NULL;
  
  return (0);
}






//other algorithm

//try to insert nodes in between edges (u,v) in current cycle, look for path (u,k,v), k can't exist, 
//because, that would imply there is a short path between u and v
