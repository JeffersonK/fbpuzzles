#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define DEBUG 0

//#define INFINITY 0xffffffffffffffff
#define INFINITY (HUGE_VAL)

#define MAX_NODES 100
#define MAX_NAME_LEN 32

int n;

double Pr[MAX_NODES];
double adjMatrix[MAX_NODES][MAX_NODES];
char names[MAX_NODES][MAX_NAME_LEN];
double best;
int visited[MAX_NODES];
int nNodes, nEdges;

double min(double a, double b) {

  if ((a == INFINITY) &&
      (b == INFINITY))
    return INFINITY;
		       
  if (a == INFINITY) return b;
  if (b == INFINITY) return a;

  if (a < b) 
    return a;
  return b;
}

int name2index(char * name){

  int i;

  for (i = 0; i < nNodes; i++)
    if (memcmp(name, names[i], strlen(name)) == 0)
      return i;

  return -1;
}

void printAdjMatrix(void) {
  
  int i, j;

  for (i=0; i < nNodes; i++) {
    for (j=0; j < nNodes; j++) {
      if (adjMatrix[i][j] == INFINITY)
	printf("INF ");
      else
	printf("%0.2lf ", adjMatrix[i][j]);
    }
    printf("\n");
  }
  
}


void floydWarshall(){
  int i, j, k;
  
  for (k = 0; k < nNodes; k++) {
    for (i = 0; i < nNodes; i++) {
      for (j = 0; j < nNodes; j++) {
	if (i != j)
	  adjMatrix[i][j] = min(adjMatrix[i][j], 
				adjMatrix[i][k] + adjMatrix[k][j]);
      } 		
    } 
  }
}

int loadFile(char * filename) {

  FILE * f = NULL;
  char locA[MAX_NAME_LEN];
  char locB[MAX_NAME_LEN];
  int i, j, u, v;
  unsigned int weight;
  
  best = INFINITY;
  
  f = fopen(filename, "r");
  if (f == NULL)
    return -1;
  
  if (fscanf(f, "%d\n", &nNodes) < 1){
    return -2;
  }
  
  /* initialize the matrix */
  for (i=0; i<nNodes; i++)
    for (j=0; j<nNodes; j++)
      if (i == j)
	adjMatrix[i][j] = 0;
      else
	adjMatrix[i][j] = INFINITY;
  
  /* read location probability */
  for (i = 0; i < nNodes; i++){
    if (fscanf(f, "%s %lf", names[i], &Pr[i]) < 2)
      return -3;

#if DEBUG
    printf("add node: %s(%i) %lf\n", names[i], i, Pr[i]);
#endif
  }

  if (fscanf(f, "%d\n", &nEdges) < 1)
    return -4;

  /* read locA locB weight */
  for (i = 0; i < nEdges; i++){
    if (fscanf(f, "%s %s %ui", locA, locB, &weight) < 3)
      return -3;
    
    u = name2index(locA);
    v = name2index(locB);
    adjMatrix[u][v] = adjMatrix[v][u] = (double) weight;
    
#if DEBUG
    printf("add edge: %s %s %i\n", locA, locB, weight);
#endif

  }

  fclose(f);
  
  return 0;
}


void walk(int node, double evDt, double probability) {
  int i, end;
  double dDt;
  visited[node] = 1;
  probability -= Pr[node];
  
  end = 1;
  for(i=0;i<nNodes;i++){
    if(!visited[i]) {
      end = 0;
      dDt = probability*adjMatrix[node][i];
      if ((evDt + dDt) < best)
	walk(i, evDt + dDt, probability);
    }      
  }
    visited[node] = 0;
    if (end) {
      //printf("%0.2lf %0.2lf\n", best, evDt);
      if (((evDt < best) && (best != -1)) || (best < 0)){
	best = evDt;
      } 
    }
}

int main(int argc, char ** argv) {

  int ret, i, j; 

  ret = loadFile(argv[1]);

#if DEBUG
  printf("loadFile: %d\n", ret);
  printf("nNodes: %d\n", nNodes);
  printf("nEdges: %d\n", nEdges);
  printAdjMatrix();
#endif

  floydWarshall();

  /* Check for any unreachable locations */
  for (i = 0; i < nNodes; i++) {
    for (j = i; j < nNodes; j++) {
      if(adjMatrix[i][j] == INFINITY) {
	printf("-1.00\n");
	return (0);
      }
    }
  }

  
  for(i = 0;i< nNodes; i++)
    visited[i] = 0;
  
#if DEBUG
  printf("\n-----\n");
  printAdjMatrix();
#endif

  walk(0, 0, 1);
  printf("%0.2lf\n", best);
  return 0;
}
