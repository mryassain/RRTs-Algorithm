#pragma once
#include "StdAfx.h"
#include "RRTTree.h"
#include <windows.h>

class RRTAlgorithmes
{
public:

	//constructeurs
	RRTAlgorithmes(void);
	RRTAlgorithmes(double*,double,int,int,int,int,double,double,double,double,double);
	~RRTAlgorithmes(void);

	//consultation
	RRTTree ConsultTrees();

    //consultation test
	double *ConsultstardiscrepencyforGRRT();
	RRTVertexe & consultexinit();
	RRTVertexe & consultexgoal();
	RRTVertexe * consultexgoals();
	RRTVertexe & consultminV();
	RRTVertexe & consultmaxV();
	RRTVertexe & consultexnew();
	RRTVertexe *consulteXsolution();
	RRTVertexe & consultexnear();
    RRTVertexe *ConsultUcontrole();
	RRTVertexe *ConsultXgoals();
	int consultNobmreOfItération();


	//modification
	void ModiftRRTreesOfAlgo(RRTVertexe a);
	void ModiftRRTAlgorithmes(RRTAlgorithmes A);
	//void InPutSystem(double tmps[],int dim,int ki,int NbrOfU,int dimu,double pas,double min1,double max1,double minu,double maxu);
		

	//fcts random
	void Random(double ,double,int );
	void Random(RRTVertexe aLower ,RRTVertexe bUpper,int i);

	//recherche l'existence
	bool IsNotHere(RRTVertexe );

	//near fction return lindice de xnear trouver dans l'arbre
	int nearfct();
   
	//fction 
    double fct(double t1,RRTVertexe x,RRTVertexe u,int i);

	//fonction pour generer le controles U
	void genereU(double min2u,double max2u,int nbrOfu,int dimu);
    //une fonction qui teste l'apprtenance du point solution dans l'ensemble d'etat 
	bool AppartienEspaceDetat(RRTVertexe a);
	double ifct(double r);

	//jai pas encore definit le K de runge kutta
	void runge4(double t2, RRTVertexe xnewrs,RRTVertexe u);
	void krunge4(double t3, RRTVertexe xnear,RRTVertexe u);
	void kPrimrunge4(double t3, RRTVertexe xnear,RRTVertexe u);
	void searchSolutionU(double t, RRTVertexe xnear);
	void searchSolutionPrimU(double t, RRTVertexe xnear);
	int CONSTRUCTREE(RRTVertexe xgl,RRTVertexe xnr,int indiceXnear);
	void CONSTRUCTREENewS(RRTVertexe xgl,RRTVertexe xnr,int indiceXnear);
	void ConstructionOfTrees(int eta);
	void ConstructionTreforFirst();
	void RemplirTri(int i);
	void GeneralRRT();
	void GeneralRRTT();
	void GeneralRRTNews();
	void GeneralGRRT();
	void GeneralGRRTT();
	void GeneralMGRRT();
	void GeneralKRRT();
	//////////////////////////////////////////////
	void GenereRandomC(int k);
	void GenereRandomFaure(int k);
	void GenereRandomHalton(int k);
	////////////////////////////////////////////mehtode 1
	void GeneralGRRTTLimite();
	double stardiscrpancy(int k);
	void localise();
	/////////////////////////////////////////////methode 2
	void Partit(RRTVertexe a,RRTVertexe b ,int k);
	void localise(int k);
	void RRTMethode2();
	double CalculVolume(RRTVertexe a ,RRTVertexe b);
	void cherchentre(RRTVertexe a,RRTVertexe b,double *tab1,double *tab2);

protected:
	RRTVertexe Xsol;
	RRTVertexe xnew;
	RRTVertexe xgoal;
	RRTVertexe xnear;                      //fabrication interne
	RRTVertexe *U;                         //input of interface
	RRTVertexe xint;                       //input of interface
    RRTVertexe *xgoals;
	RRTVertexe *xsolution;                 //fabrication interne
	RRTTree T;                            //output résult
	double t;                             //condition initiale 
	int dimension;                         //dimension of probleme 
	int K,nbrOfU,dimU;                          //k : nbre iteration and nombre of U controle
	double h,min,max;                      //h pas espace free Xfree min and max 
	double minU,maxU;
	RRTVertexe minV;
	RRTVertexe maxV;
	RRTVertexe *Partiti;
	double *stardiscrepency;
	RRTTree Tre1;
	RRTTree Tre2;
	int *PG;


};
