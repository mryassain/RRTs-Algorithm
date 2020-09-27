#pragma once


class RRTVertexe
{
public:
	RRTVertexe(void); //  ructeur par default
	RRTVertexe(double* ,int ); //  ructeur de copie
	RRTVertexe(RRTVertexe & );

	virtual ~RRTVertexe(void); // destructeur

	// les méthodes permettant de consulter les attributs
	int ConsultDim(void);
    double *ConsultVertex(void);
	double ConsultVertex(int);


	// la méthode permettant de modifier les attributs
	void ModiftVertex(double*,int);
	 
	//la modification au sein de tableau
    void ModifVertex(double ,int);

	//la modification de vertex 
	void ModifVertex(RRTVertexe & );

	//operateur pour tester l'egalite des deux vertex
    bool operator==(RRTVertexe &);
	bool operator!=(RRTVertexe &);

	//operateur d'affectation
	RRTVertexe & operator=(  RRTVertexe &) ;

	double Metrique(  RRTVertexe &);

protected:
	
	double *Vertex; // coordoné des points  
 	int dim ;      // dimension de lespace

};
