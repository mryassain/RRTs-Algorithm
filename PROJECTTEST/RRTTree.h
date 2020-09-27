#pragma once
#include "RRTVertexe.h"
#include "RRTEdge.h"

class RRTTree
{
public:

	//structures
	struct Vertexes {
		RRTVertexe vertex;
		int indiceVert;
		bool box;
		Vertexes *suivant;
	};
	
	struct Edges {
		RRTEdge edge;
		int indiceEdge;
		Edges *suivant;
	};

	
	//constructeurs
	RRTTree(void);
	~RRTTree(void);//modifier1
	
	//consultations
	int ConsultNbrVertex();
	int ConsultNbrEdge();
	double* ConsultVerto(int);

	
	Vertexes *ConsultVertexes();
	Edges *ConsultEdges();
	
    //parcori
	void ParcoruTree(Vertexes * l);

	//ajouter des vertexes et des edges
	void ajoutVertexes(RRTVertexe &);
	void ajoutEdges(RRTEdge &);

protected:

	struct Vertexes *vertex;
	struct Edges *edge;
	int nbrVertex;
	int nbrEdge;
};
