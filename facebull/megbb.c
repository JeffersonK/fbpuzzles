#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h> //for memmove

#define DEBUG_LOADFILE 1

#define INF (HUGE_VAL)
#define MAX_NODES (128)
#define MAX_EDGES (MAX_NODES*MAX_NODES) //number of elements in adjacency matrix (we can do better than this with a linked lists
#define MAX_NAME_LEN (8)

typedef struct MACHINE_T {
  short edgeset;//which edge set its part of {NO_EDGE, E_EDGE, EBAR_EDGE}
  short edgeRealName;
  short edgelabel;//label k: for k in [0 ... n-1]
  float edgecost;
  short ci;
  short cj;
  struct MACHINE_T * adjNext;
  struct MACHINE_T * adjPrev;	
  struct MACHINE_T * kthNext;//pointer to next highest edge weight
} machine_t;

machine_t * maxEdgePerRow = NULL;//linked list of the heaviest edge per row
machine_t * minEdgePerCol = NULL;//linked list of the lightest edge per row
machine_t * EdgesAdjHead = NULL;//replaces the adjaceny matrix
machine_t * EdgesCostHead = NULL;//edges sorted by cost from low to high

short machineIndex[MAX_EDGES];
short compoundIndex[MAX_NODES];

float theoryMin = 0;

//working copy
machine_t A[MAX_NODES][MAX_NODES];
//short r; //the number of edges in Ebar
short t; //the number of edges in S; edges in E, greater than the kth labled edge

//current optimal solution
machine_t Abar[MAX_NODES][MAX_NODES];
float rbar;

//number of edges in/out per node in A
short Vout[MAX_NODES];
short Vin[MAX_NODES];

short nNodes = 0;
short nEdges = 0;

int max_edges; 

#define ROWI(k) (k / nNodes)
#define COLJ(k) (k % nNodes)
#define KTH_EDGE(k) (A[ROWI(k)][COLJ(k)])
#define NO_EDGE ( (short) 0)
#define E_EDGE ( (short) 1)
#define EBAR_EDGE ( (short)-1)
   
void printMachine(machine_t * m){
  printf("Machine:%hi (@%p)\n", m->edgeRealName, m);
  printf("\tlabel:%hi\n", m->edgelabel);
  printf("\tcost:%.0f\n", m->edgecost);
  printf("\trow:%hi\n", m->ci);
  printf("\tcol: %hi\n", m->cj);
  printf("\tadjPrev: @%p\n", m->adjPrev);
  printf("\tadjNext: @%p\n", m->adjNext);
}

void printEdgeCosts(void){
  machine_t * m = EdgesCostHead;
  while(m != NULL){
    printMachine(m);
    m = m->kthNext;
  }
}

void printEdgeAdj(void){
  machine_t * ptr;
  ptr = EdgesAdjHead;

  printf(" *** FORWARDS ***\n");
  while(ptr->adjNext != NULL){
    printMachine(ptr);
    ptr = ptr->adjNext;
  }
  printMachine(ptr);

  printf(" *** REVERSE ***\n");
  while(ptr->adjPrev != NULL){
    printMachine(ptr);
    ptr = ptr->adjPrev;
  }
  printMachine(ptr);
  printf("ptr:%p EdgesAdjHead:%p\n", ptr, EdgesAdjHead);
}

void printState(){
  int i,j;
  printf("\n");
  printf("rbar = %.0f\n", rbar);
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++)
      if(Abar[i][j].edgeset < 0)
	printf("%hi ", Abar[i][j].edgeset);
      else
	printf(" %hi ", Abar[i][j].edgeset);
    printf("\n");
  }
  printf("\n");
  
  /*
  printf("r = %hi\n", r);
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++)
      if(A[i][j].edgeset < 0)
	printf("%hi ", A[i][j].edgeset);
      else
	printf(" %hi ", A[i][j].edgeset);
    printf("\n");
    }
*/
}

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

//insert from lowest to highest
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
      if (ptr->edgecost > m->edgecost){
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
  

#define CALCK(i,j) (i*nNodes + j)

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
      printf("invariant broken %d\n", __LINE__);
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

int loadFile(char * filename) {

  FILE * f = NULL;
  char compoundA[MAX_NAME_LEN];
  char compoundB[MAX_NAME_LEN];
  char machineName[MAX_NAME_LEN];
  int i, j, u, uindex, v, vindex, machIntName, machIndex;
  int cost;
  machine_t * ptr; 
  rbar = 0;
  memset(machineIndex, 0xff, sizeof(machineIndex));
  memset(compoundIndex, 0xff, sizeof(compoundIndex));
  
  /* initialize the matrix */
  for (i=0; i<MAX_NODES; i++)
    for (j=0; j<MAX_NODES; j++){
      if (i == j) {
	A[i][j].edgecost = 0;
      } else {
	A[i][j].edgecost = INFINITY;
      }
      A[i][j].edgelabel = -1;//the index in the machineList where the real name is stored
      A[i][j].edgeset = NO_EDGE;
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

    machIntName = atoi(machineName+1);

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
    A[uindex][vindex].edgeset = (short)E_EDGE;
    A[uindex][vindex].edgelabel = machIndex;
    A[uindex][vindex].edgecost = (float) cost;
    rbar += cost; //this first solution is all edges

    /*
      ptr = (machine_t*) malloc(sizeof(machine_t));  
      if (ptr == NULL)
      return MEMORY_OVERFLOW;
      memset((void *) ptr, 0x0, sizeof(machine_t));
    */
    ptr = &(A[uindex][vindex]);
    ptr->edgeset = E_EDGE;
    ptr->edgelabel = machIndex;
    ptr->edgeRealName = machIntName;
    ptr->ci = uindex;
    ptr->cj = vindex;
    ptr->edgecost = (float) cost;
    insertEdgeAdj(ptr);
    insertEdgeCost(ptr);
  }	
  //compute the sum minimum edge weights from each row is lower bound on solution
  max_edges = nNodes * nNodes;//TODO remove later
  fclose(f);
  return 0;
}

void printAnswer(void){
  int u, v;
  int first = 1;
  printf("%.0f\n", rbar);
  for(u=0;u<nNodes;u++)
    for(v=0;v<nNodes;v++)
      {
	if(Abar[u][v].edgeset == E_EDGE){
	  if (first){
	    first = 0;		
	    printf("%hi", machineIndex[A[u][v].edgelabel]);		
	  } else
	    printf(" %hi", machineIndex[A[u][v].edgelabel]);		
	}
      }
  
  printf("\n");
}

#define DONE (99)
#define NOT_DONE (0)
int updateOptimumSoln(int rr){
  int i, j;

  if (rbar < rr)
    return NOT_DONE;
	
  rbar = rr;
  //alphabar = alpha;
  //not sure if memmove will be faster than loop or not here
  //for (i=0; i<nNodes; i++)
  //memmove((void *) &Abar[i], (void *)&A[i], sizeof(short)*nNodes); 
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++){
      Abar[i][j] = A[i][j];
    }
  }

  //this is probably not a good check if are looking for
  //a minimum over edge weights because there could be
  //more than one solution with |nNodes| edges
  /*if( (nEdges - rbar) == nNodes) 
    return DONE;
  else
  return NOT_DONE;*/
  return NOT_DONE;
}

//A is the adjacency Matrix
//n is the number of nodes in G
//ii is the source node
//jj is the destination node
//returns:
//     1 if alternate path exists
//     0 otherwise
float F[MAX_NODES];
float H[MAX_NODES];
short pathExists(short src, short dst){
  short i, j, l, in, is, js, k;
  for(j=0; j<nNodes; j++){
    if(j != dst)
      F[j] = INF;
    if(j != src)
      H[j] = j;
  } 
  F[dst] = 0;
  H[src] = nNodes-1;//we count from node 0
  k = nNodes-2;// k starts as 1 less than last node
  l = src;//source node
  for(i=0; i<=k; i++){
    j = H[i];
    if (A[l][j].edgeset == E_EDGE)
      F[j] = 0;

    if (F[j] == 0){
      is = i + 1;
      if (is <= k){
	for(in=is; in<=k; in++){
	  js = H[in];
	  if (A[l][js].edgeset == E_EDGE)
	    F[js] = 0;
	}
      }
      if (F[dst] == 0)
	return 1;
      H[i] = H[k];
      l = j;
      k--;
      if ( !(k > 0) ){
	if (A[l][dst].edgeset != E_EDGE)
	  return 0;
	else
	  return 1;
      }
    }//if F[j] == 0
  }//for i in [0, ... ,k-1]
  return 0;
}



/***
 *
 *
 ***/
int megbbinit(void){
  int i, k;
  float min;
  machine_t * p;
  int r = rbar;

  for (i=0; i<nNodes; i++)
    Vout[i] = Vin[i] = 0;
  
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
  min = p->edgecost;
  theoryMin = 0;
  while(p != NULL){
    if(p->ci > i){//row changed
      printf("min in row %d is %.0f\n", i, min);
      i = p->ci;
      theoryMin += min;
      min = p->edgecost;
    } else if (p->edgecost < min){
      min = p->edgecost;
    }
    p = p->adjNext;
  }
  printf("min in row %d is %.0f\n", i, min);
  theoryMin += min;
  printf("theoretical min is %.0f\n", theoryMin);

  p = EdgesAdjHead;
  //STEP 2: do a first pass to get the first solution
  for(k=0; k<max_edges; k++){
    //TODO: while(p != NULL){
    //printf("A[%hi][%hi]=%hi Vout[%hi]=%hi Vin[%hi]=%hi\n",
    //   ROWI(k),COLJ(k),A[ROWI(k)][COLJ(k)].edgeset,ROWI(k),Vout[ROWI(k)],COLJ(k),Vin[COLJ(k)]);
    if( (A[ROWI(k)][COLJ(k)].edgeset == E_EDGE) &&
	(Vout[ROWI(k)] != 1) &&
	(Vin[COLJ(k)] != 1) )
      {
	
	A[ROWI(k)][COLJ(k)].edgeset = NO_EDGE;//set hypothesis
	if(pathExists(ROWI(k), COLJ(k))){
	  A[ROWI(k)][COLJ(k)].edgeset = EBAR_EDGE;
	  Vout[ROWI(k)]--;
	  Vin[COLJ(k)]--;
	  r -= A[ROWI(k)][COLJ(k)].edgecost;
	} 
	else {
	  A[ROWI(k)][COLJ(k)].edgeset = E_EDGE;//unset hyptohesis
	}
      }
    //TODO: p = p->adjNext;
  }
  //save current solution
  updateOptimumSoln(r);
  return NOT_DONE;
}

#define BACKWARD 0
#define FORWARD 1
void megbb_r(machine_t * k, int t, int rr, int direction){
  
  if (k == NULL)
    {
      if(rr < rbar)
	updateOptimumSoln(rr);
      return;
    }

  //TODO: rr here should be sum of min edge of each row
  if (rr == theoryMin)//nNodes == (nEdges - rr)){//theoretically optimal solution
    updateOptimumSoln(rr);
    return;
  }

  if(direction == FORWARD)/*k->edgeset == E_EDGE*/
    {
      if(k->edgeset == E_EDGE && Vout[k->ci] != 1 && Vin[k->cj] != 1)
      {  
	/*if by deleting this edge we can't do better than rbar*/   
	k->edgeset = NO_EDGE;
	if(pathExists(k->ci, k->cj))
	  {
	    k->edgeset = EBAR_EDGE;
	    Vout[k->ci]--;
	    Vin[k->cj]--;
	    megbb_r(k->adjNext, t-1, rr - k->edgecost, FORWARD);
	    k->edgeset = E_EDGE;//put edge back
	    Vout[k->ci]++;
	    Vin[k->cj]++;
	  } 
	else 
	  {
	    //no path existed, so we can't delete this edge, but continue forward move
	    k->edgeset = E_EDGE;//restore the edge
	  }
     
	/*else k->edgeset == EBAR_EDGE */
      }
    megbb_r(k->adjNext, t-1, rr, FORWARD);      
    return;
  }

  if (direction == BACKWARD)
    {
    if(k->edgeset == EBAR_EDGE)
      {
	//restore the edge then start a forward move
	k->edgeset = E_EDGE;
	Vout[k->ci]++;
	Vin[k->cj]++;
	megbb_r(k->adjNext, t - 1, rr + k->edgecost, FORWARD);
	k->edgeset = EBAR_EDGE;
	Vout[k->ci]--;
	Vin[k->cj]--;
	//don't have to do this we just came from there megbb_r(k->adjNext, t - 1, rr, FORWARD);
      } 
    else if(Vout[k->ci] != 1 && Vin[k->cj] != 1)
      {  
	//see if an alternate path exists if you delete  this edge
	k->edgeset = NO_EDGE;
	if(pathExists(k->ci, k->cj))
	  {
	    k->edgeset = EBAR_EDGE;
	    Vout[k->ci]--;
	    Vin[k->cj]--;
	    megbb_r(k->adjNext, t, rr - k->edgecost, FORWARD);
	    k->edgeset = E_EDGE;//put edge back
	    Vout[k->ci]++;
	    Vin[k->cj]++;
	  } 
	else 
	  {
	    k->edgeset = E_EDGE;//put edge back
	  }
      }
    //continue backtrack
    megbb_r(k->adjPrev, t + 1, rr, BACKWARD);
    }
  
  return;
}


int main(int argc, char ** argv){
  int ret;
  machine_t * tail;

  if( (ret=loadFile(argv[1])) < 0)
    return -1;

#if DEBUG_LOADFILE
  printf("loadFile: %d\n", ret);
  printf("nNodes: %d\n", nNodes);
  printf("nEdges: %d\n", nEdges);
  printEdgeAdj();
  //printEdgeCosts();
#endif
  
  if (nEdges == nNodes)//must be hamiltonian cycle if it is strongly connected
    return DONE;

  tail = EdgesAdjHead;
  while (tail->adjNext != NULL)
    tail = tail->adjNext;

  bzero(Abar, sizeof(Abar));
  if(megbbinit() != DONE)
    megbb_r(tail, 0, rbar, BACKWARD);
  else
    printf("early done!");
  printState();
  printAnswer();
  return (0);
}




/*
void megbb(void){

  //given the number of nodes, this give number of elements in the complete graph, 
  //or the adjacency matrix, we can optimize this later and better represent the,
  //rows as linked lists which will save memory and time
  int i;
  int k, kprime, t;//t is cardinality of set S 
  int path_existed;

  if (r == 0)//we couldn't delete any edges 
    return;
  
  t = 0;  
  max_edges = nNodes*nNodes;//don't need this when we get rid of matrix
  k = kprime = max_edges - 1;
  while(k >= 0){
    if (A[ROWI(k)][COLJ(k)].edgeset == NO_EDGE){
      k--;
    }
    else if(A[ROWI(k)][COLJ(k)].edgeset == E_EDGE){
      t++;
      k--;
    }
    else {//it's an edge we previously deleted, try forward move
      A[ROWI(k)][COLJ(k)].edgeset = E_EDGE;
      Vout[ROWI(k)]++;
      Vin[COLJ(k)]++;
      r--;	  
      if ( (r + t - (nNodes - (ROWI(k)+1)) ) <= rbar)
	{
	  //we can't delete this edge
	  t++;
	} 
      else 
	{  
	  k ++; //examine k+1 ... n edges after restoring edge k
	  while(k < max_edges)
	    {
	      if (A[ROWI(k)][COLJ(k)].edgeset != NO_EDGE)//E_BAR
		{
		  t--;
		  path_existed = 0;
		  if ( (Vout[ROWI(k)] != 1) && (Vin[COLJ(k)] != 1) )
		    {
		      A[ROWI(k)][COLJ(k)].edgeset = NO_EDGE;//hypothesis
		      if(pathExists(ROWI(k),COLJ(k)))
			{
			  A[ROWI(k)][COLJ(k)].edgeset = EBAR_EDGE;//hypothesis
			  //TODO: subtract edge weight from alpha only if wasn't deleted already
			  //alpha -= A[ROWI(k)][COLJ(k)].edgecost;
			  Vout[ROWI(k)]--;
			  Vin[COLJ(k)]--;
			  r++;//increment count of edges in Ebar
			  path_existed = 1;
			} 
		      else 
			{
			  A[ROWI(k)][COLJ(k)].edgeset = E_EDGE;//undo hypothesis
			  //add edge weight back to alpha
			  //alpha += A[ROWI(k)][COLJ(k)].edgecost;
	
			}
		    }
		  
		  if (!path_existed && ((r + t - (nNodes - (ROWI(k)+1)) ) <= rbar) )
		    {
		      t ++;
		      break; //we can't remove this edge
		    }
		  
		  if ( (t - (nNodes-ROWI(k))) == 0)
		    {
		      if(updateOptimumSoln(r) == DONE)
			return;
		      break;
		    }   		
		}
	      	      
	      k++;
	    }//while (k < max_edges)
	  
	  if(updateOptimumSoln(r) == DONE)
	    return;
	}
      
      
    }

  }//while (k > 0)

  return;
}
*/
