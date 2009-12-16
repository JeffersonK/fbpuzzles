#include <stdlib.h>
#include <stdio.h>
//#include <math.h>
#include <string.h> //for memmove

#define DEBUG_LOADFILE 1
#define DEBUG_MEGBB_R 1

#define INF ((2^31)-1)
#define MAX_NODES (128)
#define MAX_EDGES (MAX_NODES*MAX_NODES) //number of elements in adjacency matrix (we can do better than this with a linked lists
#define MAX_NAME_LEN (8)

typedef struct MACHINE_T {
  short edgeset;//which edge set its part of {NO_EDGE, E_EDGE, EBAR_EDGE}
  short edgeRealName;
  short edgelabel;//label k: for k in [0 ... n-1]
  short inAbar; //whether edge is part of optimal solution or not
  int edgecost;
  short ci;
  short cj;
  struct MACHINE_T * adjNext;
  struct MACHINE_T * adjPrev;	
  struct MACHINE_T * kthNext;//pointer to next highest edge weight
} machine_t;

#define BACKWARD 0
#define FORWARD 1

#define DONE (99)
#define NOT_DONE (0)

int nS;//cardinality of edges in S, where S is set of edges for i > k in E
int sumS;//sum of edge weights for edges in S
//Ebar is the set of edges that have been deleted from E
int nEbar;//cardinality of Ebar
int sumEbar;//sum of edge weights in Ebar, we want to maximize this in our solution
int sumEbarprime;//current max sum of edges we can delete

#define INSERT_EBAR(m) m->edgeset=EBAR_EDGE;sumEbar+=m->edgecost;nEbar++;Vout[m->ci]--;Vin[m->cj]--
#define INSERT_E(m) m->edgeset=E_EDGE;sumEbar-=m->edgecost;nEbar--;Vout[m->ci]++;Vin[m->cj]++
#define INSERT_S(m) nS++;sumS+=m->edgecost
#define REMOVE_S(m) nS--;sumS-=m->edgecost
#define EDGE_IN_E(m) (m->edgeset == E_EDGE)
#define EDGE_IN_EBAR(m) (m->edgeset == EBAR_EDGE)



//machine_t * maxEdgePerRow = NULL;//linked list of the heaviest edge per row
//machine_t * minEdgePerCol = NULL;//linked list of the lightest edge per row
machine_t * EdgesAdjHead = NULL;//replaces the adjaceny matrix
machine_t * EdgesCostHead = NULL;//edges sorted by cost from low to high

short machineIndex[MAX_EDGES];
short compoundIndex[MAX_NODES];

//the ith element contains the sum of the min edges in the last i rows
int minSumRows[MAX_NODES];
int theoryMin;//TODO: this is redundant remove
int graphWeight;

//working copy
machine_t A[MAX_NODES][MAX_NODES];

//int rbar;

//number of edges in/out per node in A
short Vout[MAX_NODES];
short Vin[MAX_NODES];

short nNodes = 0;
short nEdges = 0;


#define ROWI(k) (k / nNodes)
#define COLJ(k) (k % nNodes)
#define KTH_EDGE(k) (A[ROWI(k)][COLJ(k)])
#define NO_EDGE ( (short) 0)
#define E_EDGE ( (short) 1)
#define EBAR_EDGE ( (short)-1)
   
void printMachine(machine_t * m){
  printf("Machine:%hi (@%p)\n", m->edgeRealName, m);
  printf("\tlabel:%hi\n", m->edgelabel);
  printf("\tcost:%d\n", m->edgecost);
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
  printf("--- STATE ---\n");
  printf("Best(%d) sumEbarprime(%d) sumEbar(%d) nEbar(%d) nS(%d) sumS(%d)\n", 
	 graphWeight-sumEbarprime, sumEbarprime, sumEbar, nEbar, nS, sumS);
  for(i=0; i<nNodes; i++){
    for(j=0; j<nNodes; j++)
      if(A[i][j].edgeset < 0)
	printf("%hi(%hi) ", A[i][j].edgeset, A[i][j].inAbar);
      else
	printf(" %hi(%hi) ", A[i][j].edgeset, A[i][j].inAbar);
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
  int i, j, u, uindex, v, vindex;
  short machIntName, machIndex;
  int cost;
  machine_t * ptr; 
  //rbar = 0;
  graphWeight = 0;
  memset(machineIndex, 0xff, sizeof(machineIndex));
  memset(compoundIndex, 0xff, sizeof(compoundIndex));
  
  /* initialize the matrix */
  for (i=0; i<MAX_NODES; i++)
    for (j=0; j<MAX_NODES; j++){
      if (i == j) {
	A[i][j].edgecost = 0;
      } else {
	A[i][j].edgecost = INF;//max signed int
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
    A[uindex][vindex].edgecost = cost;
    graphWeight += cost;
    //    rbar += cost; //this first solution is all edges

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
    ptr->edgecost = cost;
    
    insertEdgeAdj(ptr);
    insertEdgeCost(ptr);
  }	
  //compute the sum minimum edge weights from each row is lower bound on solution
  fclose(f);
  return 0;
}

void printAnswer(void){
  machine_t * p = EdgesAdjHead;
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
    p = p->adjNext;
  }
  printf("\n");
}

void updateOptimumSoln(int ebar){
  
  machine_t * p = EdgesAdjHead;
  sumEbarprime = ebar;
  while (p != NULL){
    if (p->edgeset == E_EDGE)
      p->inAbar = 1;
    else
      p->inAbar = 0;    
    p = p->adjNext;
  }
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
  int min;
  machine_t * p;
  //int r = rbar;//rbar is sum of all edges at this point
  
  for (i=0; i<nNodes; i++)
    minSumRows[i] = Vout[i] = Vin[i] = 0;
  
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
      printf("min in row %d is %d\n", i, min);
      minSumRows[i] = min;
      i = p->ci;
      theoryMin += min;
      min = p->edgecost;
    } else if (p->edgecost < min){
      min = p->edgecost;
    }
    p = p->adjNext;
  }
  minSumRows[i] = min;
  printf("min in row %d is %d\n", i, min);
  theoryMin += min;
  printf("theoretical minimum solution is %d\n", theoryMin);

  min = minSumRows[nNodes-1];
  for(i=nNodes-1;i>=0;i--){
    min = minSumRows[i];
    if (i == nNodes - 1)
      minSumRows[i] = min;
    else
      minSumRows[i] = min + minSumRows[i+1];
    printf("minSumRows[%d] = %d\n", i, minSumRows[i]);
  }
  //minSumRows[0] = theoryMin;

  nS = 0;
  sumS = 0;
  sumEbar = 0;
  nEbar = 0;
  sumEbarprime = 0;
  p = EdgesAdjHead;
  //STEP 2: do a first pass to get the first solution
  while(p != NULL){
    //printf("A[%hi][%hi]=%hi Vout[%hi]=%hi Vin[%hi]=%hi\n",
    //   ROWI(k),COLJ(k),A[ROWI(k)][COLJ(k)].edgeset,ROWI(k),Vout[ROWI(k)],COLJ(k),Vin[COLJ(k)]);
    if( (EDGE_IN_E(p)) && (Vout[p->ci] != 1) && (Vin[p->cj] != 1) ){
      
      p->edgeset = NO_EDGE;//set hypothesis
      if(pathExists(p->ci, p->cj)){
	INSERT_EBAR(p);
      } 
      else 
	p->edgeset = E_EDGE;//unset hyptohesis
    }
    p = p->adjNext;
  }
  //save current solution
  updateOptimumSoln(sumEbar);
  return NOT_DONE;
}

void megbb_r(machine_t * k, int direction){

  if (k == NULL && direction == FORWARD)
    {//we have exhausted the forward move
#if DEBUG_MEGBB_R
      if ( (nS != 0) && (sumS != 0) ) 
	printf("Forward Move Finished: invariant broken ns=%d sumS=%d\n", nS, sumS);
#endif
      if(sumEbar > sumEbarprime)
	updateOptimumSoln(sumEbar);
      return;
    }
  else if (k == NULL){
#if DEBUG_MEGBB_R
    if (sumEbar != 0)
      printf("Backtrack Finished: invariant broken sumEbar = %d\n", sumEbar);
#endif
    //megbb_r(EdgesAdjHead, FORWARD);
    return;//done back tracking
  }

  if ( (graphWeight - sumEbar) == minSumRows[0]){
    //we have found the theoretical optimum solution for this graph
    //i.e. - the current solution has the minimum edge from each row
    updateOptimumSoln(sumEbar);
    return;
  }

  printf("@%d(%d) ", k->edgeRealName, direction);
  if(direction == FORWARD){/*k->edgeset == E_EDGE*/
    //if ( (rr + t - minSumRows[k->ci]) <= rbar)
    //  return;
    if(EDGE_IN_E(k)){
      if( (Vout[k->ci] != 1) && (Vin[k->cj] != 1) )
      {  
	/*if by deleting this edge we can't do better than rbar*/   
	k->edgeset = NO_EDGE;
	if(pathExists(k->ci, k->cj)){
	  INSERT_EBAR(k);
	  REMOVE_S(k);
	  megbb_r(k->adjNext, FORWARD);    
	  INSERT_S(k);
	  INSERT_E(k);
	} 
	else 
	  { 
	    //no path existed, so we can't delete this edge, but continue forward move
	    k->edgeset = E_EDGE;//restore the edge
	    REMOVE_S(k);
	    megbb_r(k->adjNext, FORWARD);
	  }
      }
    } 
    else {//EDGE_IN_EBAR(k) == true
      megbb_r(k->adjNext, FORWARD);    
    }
    
    if (k->adjNext == NULL){
      //start back tracking move when we have exhausted forward move
      	INSERT_S(k);
	megbb_r(k->adjPrev, BACKWARD);
	REMOVE_S(k);
	//don't have to update optimal because we are at last edge and just increased the current solution
    }
    return;
  }//FORWARD

  //DIRECTION IS BACKWARD
  if(EDGE_IN_EBAR(k)){
    //restore the edge then start a forward move
    INSERT_E(k);
    //! REMOVE_S(k) because k can't be in S
    megbb_r(k->adjNext, FORWARD);
    INSERT_S(k);
    megbb_r(k->adjPrev, BACKWARD);
    REMOVE_S(k);
  }
  else {//if (EDGE_IN_E(k))
    //REMOVE_S(k);
    //megbb_r(k->adjNext, FORWARD);
  }
 
  megbb_r(k->adjPrev, BACKWARD);  
  return;
}

void forward(machine_t *k){


}

void megbb_r2(machine_t * k){

  machine_t * p = k;

  if (k == NULL){
    if (sumEbar > sumEbarprime)
      updateOptimumSoln(sumEbar);
    return;
  }

  while (p != NULL){
    if (EDGE_IN_E(k)){
      //INSERT_S(k);
      megbb_r2(k->adjPrev);
      //try to remove it
      if( (Vout[k->ci] != 1) && (Vin[k->cj] != 1) )
	{  
	  /*if by deleting this edge we can't do better than rbar*/   
	  k->edgeset = NO_EDGE;
	  if(pathExists(k->ci, k->cj)){
	    INSERT_EBAR(k);
	    //REMOVE_S(k);
	    megbb_r2(k->adjNext);
	    //INSERT_S(k);
	    //INSERT_E(k);
	  }
	else 
	  { 
	    //no path existed, so we can't delete this edge, but continue forward move
	    k->edgeset = E_EDGE;//restore the edge
	    //REMOVE_S(k);
	    megbb_r2(k->adjNext);
	  }
      
	}
    }
 
    if (EDGE_IN_EBAR(k)){
      INSERT_E(k);
      megbb_r2(k->adjNext);
    }

    //INSERT_S(p);
    p = p->adjPrev;
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
  printf("graphWeight: %d\n", graphWeight);
  printEdgeAdj();
  //printEdgeCosts();
#endif
  
  if (nEdges == nNodes){
    //must be hamiltonian cycle if it is strongly connected
    updateOptimumSoln(graphWeight);
    printAnswer();
    return DONE;
  }

  //find last element in adjacency list
  tail = EdgesAdjHead;
  while (tail->adjNext != NULL)
    tail = tail->adjNext;

  if(megbbinit() != DONE){
    printState();
    megbb_r2(tail);//, BACKWARD);

  }
  else
    printf("early done!");
  
  printState();
  printAnswer();
  return (0);
}

/*
int r = rbar;
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
	  //kprime = k;
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
      //k = kprime - 1;
      
    }

  }//while (k > 0)

  return;
}

*/
