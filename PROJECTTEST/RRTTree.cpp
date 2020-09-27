#include "StdAfx.h"
#include "RRTTree.h"
#include <windows.h>

RRTTree::RRTTree(void)
{
	vertex = NULL;
	edge = NULL;
	nbrVertex = nbrEdge = 0;
}

RRTTree::~RRTTree(void)
{
/*
	Vertexes *p;
	p=vertex;
	while(p->indiceVert!=0)
	{
		struct Vertexes* l;
		p=p->suivant;
		l=p;
		//if(l->vertex.ConsultVertex()!=NULL)l->vertex.~RRTVertexe();
		delete [] l;      		
	}
//	if(vertex->indiceVert!=0)delete vertex;
    Edges *l;
	l=edge;
	while(l!=NULL)
	{
		struct Edges *p;
		l=l->suivant;
		p=l;
		//if(p->edge.ConsultDist()!=NULL)p->edge.~RRTEdge();
		delete [] p;
	}
	delete [] edge;

	*/}
void RRTTree::ParcoruTree(Vertexes *l)
{  
	if(l)vertex=l;
}
	
//consultations
int RRTTree::ConsultNbrVertex()
{ 
	return nbrVertex;
}

int RRTTree::ConsultNbrEdge()
{
	return this->nbrEdge;
}
double* RRTTree:: ConsultVerto(int i)
{
	struct Vertexes *l;
	l=vertex;
	int k=0;
	while(k<i)
	{
		l=l->suivant;
		k++;
	}
	return l->vertex.ConsultVertex();
}
	
RRTTree::Vertexes * RRTTree::ConsultVertexes()
{
	return vertex;
}

RRTTree::Edges * RRTTree::ConsultEdges()
{
	return edge;
}
	
//ajouter des vertexes et des edges
void RRTTree::ajoutVertexes(RRTVertexe &vert)
{
	RRTTree::Vertexes *p;
	RRTTree::Vertexes *l; 
	l=this->vertex;
	p = new struct Vertexes;
	p->vertex=vert; 
	p->indiceVert = ++nbrVertex;
	p->box=false;
	p->suivant=l;
	this->vertex = p;
}
void RRTTree::ajoutEdges(RRTEdge &edg)
{
	Edges *p ,*l;
	l= this->edge;
	p = new struct Edges;
      
	p->edge = edg;
	p->indiceEdge = ++nbrEdge;
	p->suivant=l;
    edge = p;
}
