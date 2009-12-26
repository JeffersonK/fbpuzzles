#include <stdlib.h>
#include <stdio.h>
#include <string.h> //for bzero()
#include <math.h>

#define DEBUG_LOADFILE 0
#define DEBUG_MEGBB 0
#define WITH_PRUNING 1
#define DEBUG_PRUNING 0

#define INFINITY (HUGE_VAL)
#define INF ((2^31)-1)
#define MAX_NODES (512)
#define MAX_EDGES (MAX_NODES*MAX_NODES) //number of elements in adjacency matrix (we can do better than this with a linked lists
#define MAX_NAME_LEN (8)

typedef struct MACHINE_T {
  short edgeset;//which edge set its part of {NO_EDGE, E_EDGE, EBAR_EDGE}
  short edgeRealName;
  //short edgelabel;//label k: for k in [0 ... n-1]
  short inAbar; //whether edge is part of optimal solution or not
  int edgecost;
  short ci;
  short cj;
  struct MACHINE_T * adjNext;
  struct MACHINE_T * adjPrev;	
  struct MACHINE_T * kthNext;//pointer to next highest edge weight
  struct MACHINE_T * nameNext;//pointer to next highest edge weight
  struct MACHINE_T * nextColMin;
  struct MACHINE_T * nextRowMin;
} machine_t;

#define DONE (99)
#define NOT_DONE (0)
#define MEMORY_OVERFLOW (-666)

#define NO_EDGE ( (short) 0)
#define E_EDGE ( (short) 1)
#define EBAR_EDGE ( (short)-1)

#define INSERT_EBAR(m) m->edgeset=EBAR_EDGE;sumEbar+=m->edgecost;nEbar++;Vout[m->ci]--;Vin[m->cj]--
#define INSERT_E(m) m->edgeset=E_EDGE;sumEbar-=m->edgecost;nEbar--;Vout[m->ci]++;Vin[m->cj]++
#define INSERT_S(m) nS++;sumS+=m->edgecost
#define REMOVE_S(m) nS--;sumS-=m->edgecost
#define EDGE_IN_E(m) (m->edgeset == E_EDGE)
#define EDGE_IN_EBAR(m) (m->edgeset == EBAR_EDGE)

#define CALCK(i,j) (i*nNodes + j)

typedef struct COMPOUND_T {
  short mark;
  short predecessor;
  float minpathcost;  
  machine_t * m;
} compound_t;

compound_t A[MAX_NODES][MAX_NODES];

int nS;//cardinality of edges in S, where S is set of edges for i >= k in E
int sumS;//sum of edge weights for edges in S

//Ebar is the set of edges that have been deleted from E
int nEbar;//cardinality of Ebar
int sumEbar;//sum of edge weights in Ebar, we want to maximize this in our solution

//current bests
int nEbarprime;//current number of edges in Ebarprime
int sumEbarprime;//current max sum of edge weights we can delete

machine_t * EdgesAdjHead = NULL;//replaces the adjaceny matrix
machine_t * EdgesNameHead = NULL;//edges sorted by real name (machine name)
machine_t * EdgesCostHead = NULL;//edges sorted in decreasing weight

short compoundIndex[MAX_NODES];

//the ith element contains the sum of the min edges in the last i rows
int minSumRows[MAX_NODES];
//the ith element contains the sum of the min edges in the last i rows
int maxSumRows[MAX_NODES];

machine_t * rowMin[MAX_NODES];
machine_t * rowMax[MAX_NODES];
machine_t * colMin[MAX_NODES];
machine_t * colMax[MAX_NODES];

//sum of all edges in graph
int graphWeight;

//useful book keeping datastructures to know quickly if we can delete an edge or not
//out degree of node i in current graph 
short Vout[MAX_NODES];
//out degree of node i in current graph 
short Vin[MAX_NODES];

int nNodes = 0;
int nEdges = 0;

/*
#define ROWI(k) (k / nNodes)
#define COLJ(k) (k % nNodes)
#define KTH_EDGE(k) (A[ROWI(k)][COLJ(k)])
*/

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
  fprintf(stderr, "\tedgeset:%d\n", m->edgeset);
  fprintf(stderr, "\tcost:%d\n", m->edgecost);
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

void printMinPaths(void){
  int i, j;
  fprintf(stderr, "*** MIN PATHS ***\n");
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++)
      fprintf(stderr, "%04.0f(%d)\t", A[i][j].minpathcost, A[i][j].predecessor);
    fprintf(stderr, "\n");
  }
}

/****
 *
 ****/   
void printState(){
  int i,j;
  machine_t * p = EdgesAdjHead;
  fprintf(stderr, "--- STATE ---\n");
  //printColsFilled();
  fprintf(stderr, "Best(%d) sumEbarprime(%d) sumEbar(%d) nEbar(%d) nS(%d) sumS(%d)\n", 
	  graphWeight-sumEbarprime, sumEbarprime, sumEbar, nEbar, nS, sumS);
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++){
      if (p!=NULL && p->ci==i && p->cj==j)
	{
	  if(p->edgeset < 0)
	    fprintf(stderr, "%hi(%hi) ", p->edgeset*p->edgecost, p->inAbar);
	  else
	    fprintf(stderr, " %hi(%hi) ", p->edgeset*p->edgecost, p->inAbar);
	  p = p->adjNext;
	}
      else {//no edge
	fprintf(stderr, " 0000(0) ");
      }
    }//j
    fprintf(stderr, "\n");
  }//i
  fprintf(stderr, "\n");
  
}

/*************************
 *
 * HELPER FUNCTIONS
 *
 ************************/

inline int getColMin(machine_t * m){
  machine_t * p = colMin[m->cj];
  while(p != NULL){
    if(EDGE_IN_E(p) && (p->ci >= m->ci) ){
      return p->edgecost;
    }
    p = p->nextColMin;
  }
  return 0;
}

inline int getRowMin(int u){
  if (u < 0 || u >= nNodes) fprintf(stderr, "broken invariant\n");
  machine_t * p = rowMin[u];
  while(p != NULL){
    if(EDGE_IN_E(p)){
      return p->edgecost;
    }
    p = p->nextRowMin;
  }
  return 0;
}

/****
 *
 ****/   
void insertColMin(machine_t * m){
  machine_t * last, * p = colMin[m->cj];  
  if(p == NULL){
    m->nextColMin = NULL;
    colMin[m->cj] = m;
    return;
  }
  
  if(p->edgecost > m->edgecost){
    m->nextColMin = p;
    colMin[m->cj] = m;
    return;
  }

  last = p;
  while(p!=NULL){
    if (p->edgecost > m->edgecost){
      m->nextColMin = p;
      last->nextColMin = m;
      return;
    }
    last = p;
    p = p->nextColMin;
  }
    
  m->nextColMin = NULL;
  last->nextColMin = m;
  return;
}

/****
 *
 ****/   
void insertRowMin(machine_t * m){
  machine_t * last, * p = rowMin[m->ci];  
  if(p == NULL){
    m->nextRowMin = NULL;
    rowMin[m->ci] = m;
    return;
  }
  
  if(p->edgecost >= m->edgecost){
    m->nextRowMin = p;
    rowMin[m->ci] = m;
    return;
  }

  last = p;
  p = p->nextRowMin;
  while(p!=NULL){
    if (p->edgecost >= m->edgecost){
      m->nextRowMin = p;
      last->nextRowMin = m;
      return;
    }
    last = p;
    p = p->nextRowMin;
  }
    
  m->nextRowMin = NULL;
  last->nextRowMin = m;
  return;
}

/****
 *
 ****/   
void freeEdges(machine_t * m){
  if (m == NULL)
    return;  
  //free everyone after you
  freeEdges(m->adjNext);
  //clear him out for good measure, we might get this reallocated to us later
  bzero(m, sizeof(machine_t));
  //free yourself
  free(m);
  return;
}

/****
 *
 ****/   
short insert2index(short name, short * list, int * listlen, int maxlen){
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
  printf("%d\n", graphWeight - sumEbarprime);  
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

  graphWeight = 0;
  memset(compoundIndex, 0xff, sizeof(compoundIndex));
  bzero(colMin, sizeof(machine_t *)*MAX_NODES);
  bzero(colMax, sizeof(machine_t *)*MAX_NODES);
  bzero(rowMin, sizeof(machine_t *)*MAX_NODES);
  bzero(rowMin, sizeof(machine_t *)*MAX_NODES);

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
    if (nEdges >= MAX_NODES)
      return MEMORY_OVERFLOW;

    ptr = (machine_t*) malloc(sizeof(machine_t));  
    if (ptr == NULL)
      return MEMORY_OVERFLOW;
    memset((void *) ptr, 0x0, sizeof(machine_t));

    A[uindex][vindex].m = ptr;
    ptr->edgeset = E_EDGE;
    ptr->edgeRealName = machIntName;
    ptr->ci = uindex;
    ptr->cj = vindex;
    ptr->edgecost = cost;

    graphWeight += cost;
    
    insertEdgeAdj(ptr);
    insertEdgeName(ptr);
    insertEdgeCost(ptr);
    
    insertRowMin(ptr);
    insertColMin(ptr);

    nEdges ++;
  }	
  fclose(f);
  return 0;
}

/*****
 *
 *
 *****/
void updateOptimumSoln(int ebar){
  machine_t * p = EdgesAdjHead;
  sumEbarprime = sumEbar;
  nEbarprime = nEbar;
  while (p != NULL){
    if (p->edgeset == E_EDGE)
      p->inAbar = 1;
    else
      p->inAbar = 0;    
    p = p->adjNext;
  }
#if DEBUG_MEGBB
  printState();
#endif
}

/*****
 *
 * breadth-first search path existance algorithm
 *
 *****/
char mark[MAX_NODES];
int pathExists(short src, short dst){
  int new_marks = 0;
  machine_t * p = EdgesAdjHead;
  bzero(mark, nNodes);
  mark[src] = 1;  
  while(1){
    if(EDGE_IN_E(p) && mark[p->ci] && !mark[p->cj]){
      if (p->cj == dst)
	return 1;
      mark[p->cj]++;
      new_marks++;
    }
    if(p->adjNext == NULL){
      if(new_marks == 0)
	return mark[dst];
      new_marks = 0;
      p = EdgesAdjHead;
    } else 
      p = p->adjNext;
  }
  return mark[dst];
}

/***
 *
 * Minimum Equaivalent Branch and Bound Algorithm
 *
 ***/
void megbbinit(void){
  int i, min, max;
  machine_t * p;
  
  for (i=0; i<nNodes; i++){
    maxSumRows[i] = minSumRows[i] = Vout[i] = Vin[i] = 0;
  }
  //STEP 1
  p = EdgesAdjHead;
  while(p != NULL)
    {
      if (p->edgeset == E_EDGE){
	Vout[p->ci]++;
	Vin[p->cj]++;
      }
      p = p->adjNext;
    }
  p = EdgesAdjHead;
  i = p->ci;
  max = min = p->edgecost;
  while(p != NULL){
    if(p->ci > i){//row changed

#if DEBUG_MEGBB
      fprintf(stderr, "min in row %d is %d\n", i, min);
      fprintf(stderr, "max in row %d is %d\n", i, max);
#endif
      minSumRows[i] = min;
      maxSumRows[i] = max;
      i = p->ci;
      min = p->edgecost;
      max = p->edgecost;

    } else if (p->edgecost < min){
      min = p->edgecost;

    } else if (p->edgecost > max){
      max = p->edgecost;
    }
    p = p->adjNext;
  }
  minSumRows[i] = min;
  maxSumRows[i] = max;

#if DEBUG_MEGBB
  fprintf(stderr, "min in row %d is %d\n", i, min);
  fprintf(stderr, "max in row %d is %d\n", i, max);
#endif
  max = maxSumRows[nNodes-1];
  min = minSumRows[nNodes-1];
  maxSumRows[nNodes] = 0;
  minSumRows[nNodes] = 0;
  for(i=nNodes-1;i>=0;i--){
    min = minSumRows[i];
    max = maxSumRows[i];
    if (i == nNodes - 1){
      minSumRows[i] = min;
      maxSumRows[i] = max;
    } else {
      minSumRows[i] = min + minSumRows[i+1];
      maxSumRows[i] = max + maxSumRows[i+1];
    }
#if DEBUG_MEGBB
    fprintf(stderr, "minSumRows[%d] = %d\n", i, minSumRows[i]);
    fprintf(stderr, "maxSumRows[%d] = %d\n", i, maxSumRows[i]);
#endif
  }
#if DEBUG_MEGBB
  fprintf(stderr, "theoretical minimum solution is %d\n", minSumRows[0]);
#endif
  nS = 0;
  sumS = 0;
  sumEbar = 0;
  nEbar = 0;
  sumEbarprime = 0;

  p = EdgesCostHead;
  while(p != NULL){
#if DEBUG_MEGBB
    fprintf(stderr, "(%hi,%hi)=%hi Vout[%hi]=%hi Vin[%hi]=%hi\n",
	    p->ci,p->cj,p->edgeset,p->ci,Vout[p->ci],p->cj,Vin[p->cj]);
#endif
    if( (EDGE_IN_E(p)) && (Vout[p->ci] != 1) && (Vin[p->cj] != 1) ){      
      p->edgeset = NO_EDGE;//set hypothesis
      if(pathExists(p->ci, p->cj)){
	INSERT_EBAR(p);
	
      } 
      else {
	p->edgeset = E_EDGE;//unset hyptohesis
      }
    }
    p = p->kthNext;
  }
  updateOptimumSoln(sumEbar);

  //reset the edges
  p = EdgesCostHead;
  while (p != NULL){
    if (EDGE_IN_EBAR(p)){
      INSERT_E(p);
    }
    p = p->kthNext;
  }

  p = EdgesAdjHead;
  //STEP 2: do a first pass to get the first solution
  while(p != NULL){
#if DEBUG_MEGBB
    fprintf(stderr, "(%hi,%hi)=%hi Vout[%hi]=%hi Vin[%hi]=%hi\n",
	    p->ci,p->cj,p->edgeset,p->ci,Vout[p->ci],p->cj,Vin[p->cj]);
#endif
    if( (EDGE_IN_E(p)) && (Vout[p->ci] != 1) && (Vin[p->cj] != 1) ){
      
      p->edgeset = NO_EDGE;//set hypothesis
      if(pathExists(p->ci, p->cj)){
	INSERT_EBAR(p);
	//colsFilled[p->cj]--;
      } 
      else {
	p->edgeset = E_EDGE;//unset hyptohesis
      }
    }
    p = p->adjNext;
  }
  //save current solution
  if (sumEbar > sumEbarprime)
    updateOptimumSoln(sumEbar);
  return;
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


void dijkstra_single_source_shortest_paths(int source){
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

/********
 *
 * 
 *
 *******/
void megbb(machine_t * tail){
  machine_t * k = tail;
  int path_existed;

  if (nEbar == 0)//we couldn't delete any edges 
    return;

  nS = 0;
  sumS = 0;
  while(k!=NULL){
    if(EDGE_IN_E(k)){
      INSERT_S(k);
      k = k->adjPrev;
    }
    else if (k->edgecost > A[k->ci][k->cj].minpathcost){
      //no point in adding it back if we know its not part of soln
      INSERT_S(k);
      k = k->adjPrev;

      } else 
      {//it's an edge we previously deleted, try forward move
	    INSERT_E(k);

#if WITH_PRUNING	
	    if((sumEbar + sumS - minSumRows[k->ci+1]) <= sumEbarprime)
	  //if ( (nEbar + nS - (nNodes - (k->ci+1)) ) <= nEbarprime)
	      {
		INSERT_S(k);	  
	      } 
	    else 
#endif 
	  {
	    //examine k+1 ... n edges after restoring edge k
	    k = k->adjNext;	    
	    while(k != NULL)//(k < max_edges)
	      {
		if(!EDGE_IN_E(k))
		 goto next_edge;

		REMOVE_S(k);
		path_existed = 0;
		if ( (Vout[k->ci] != 1) && (Vin[k->cj] != 1) )
		  {
		    k->edgeset = NO_EDGE;//set hypothesis
		    if(pathExists(k->ci,k->cj))
		    {
		      INSERT_EBAR(k);
		      //colsFilled[k->cj] --;
		      path_existed = 1;
		    } 
		  else 
		    {
		      k->edgeset = E_EDGE;//unset hypothesis
		    }
		}

#if WITH_PRUNING		
		if(!path_existed && ((sumEbar + sumS - minSumRows[k->ci+1]) <= sumEbarprime))
		  //(!path_existed && ((nEbar + nS - (nNodes - (k->ci+1)) ) <= nEbarprime) )
		  {
#if DEBUG_PRUNING
		    fprintf(stderr, "!!! PRUNING BRANCH (machine:%hi (%hi,%hi,%d) !!!\n", 
			    k->edgeRealName, k->ci, k->cj, k->edgecost);
		    printState();
#endif
		    INSERT_S(k);
		    break;
		  }		
		//we currently have the minimum number of edges
		//the graph can contain and be strongly connected
		if ( (nS - (nNodes-(k->ci+1))) == 0){
		  //doens't mean is optimal solution however
		  if (sumEbar > sumEbarprime)
		    updateOptimumSoln(sumEbar);
		  break;
		}

		//we have reached the theoretical 
		//minimum solution for this graph
		//i.e. sum of all row minimums
		if((graphWeight - sumEbar) == minSumRows[0]){
		  updateOptimumSoln(sumEbar);
		  return;
		}
#endif
		if (sumEbar > sumEbarprime){
		  updateOptimumSoln(sumEbar);
		}
	      next_edge:
	      k = k->adjNext;
	      }//while (kk != NULL)

	    if (k==NULL)
	      k = tail;
	    	    
	    if(sumEbar > sumEbarprime){
	      updateOptimumSoln(sumEbar);
	    }
	 }
      }
  }//while (k > 0)

#if DEBUG_MEGGB
  if (sumS != graphWeight)
    fprintf(stderr, "sumS: %d != graphWeight: %d\n", sumS, graphWeight);
#endif

  return;
}

int main(int argc, char ** argv){
  int ret, i;
  machine_t * tail;
  
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
    updateOptimumSoln(graphWeight);
    printAnswer();
    return DONE;
  }

  for(i=0; i<nNodes; i++)
    dijkstra_single_source_shortest_paths(i);
  //printMinPaths();

  //find last element in adjacency list
  tail = EdgesAdjHead;
  while (tail->adjNext != NULL)
    tail = tail->adjNext;

  //initialize stuff and do first pass solution
  megbbinit();

  //if we can't delete any edges or we have already found the theoretically 
  //minimum solution in the first pass (sum of min weight edges from each row)
  //we're done
  if( (nEbar == 0) || (minSumRows[0] == (graphWeight - sumEbarprime)) ){
    printAnswer();
    goto done;
  } 
  megbb(tail);
  
  printAnswer();
  done:
  //cleanup for good form
  freeEdges(EdgesAdjHead);
  EdgesAdjHead = EdgesNameHead = EdgesCostHead = NULL;
  return (0);
}
