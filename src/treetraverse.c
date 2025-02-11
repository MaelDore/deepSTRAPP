// Source: R package BAMMtools
// Authors: Dan Rabosky, Mike Grundler
// References: http://bamm-project.org/

#include <R.h>
#include <Rinternals.h>

void setrecursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit);
void recursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit);

static int zkzkz;

void setrecursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit) {
	zkzkz = 0;
	recursivesequence(anc,desc,node,ne,downseq,lastvisit);
}

void recursivesequence(int * anc, int * desc, int * node, int * ne, int * downseq, int * lastvisit) {
	downseq[zkzkz] = *node; zkzkz++; 
	int i, d = 0;
	int * children;
	children = R_Calloc(2, int);
	for (i = 0; i < *ne; i++) {
		if (anc[i] == *node) {
			children[d] = desc[i];
			d++;
		}
		if (d == 2) break;
	}
	if (children[0] != 0 && children[1] != 0) {
		int * child;
		child = R_Calloc(1, int);
		for (i = 0; i < 2; i++)
		{
			*child = children[i];
			recursivesequence(anc,desc,child,ne,downseq,lastvisit);
		}
		R_Free(child);
	}
	for (i = 0; i < *ne+1; i++) if (downseq[i] == 0) break;
	lastvisit[*node-1] = downseq[i-1];  
	R_Free(children);	
}

