// PROJECTTEST.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "math.h"
#include "RRTTree.h"
#include "RRTEdge.h"
#include "RRTVertexe.h"
#include "RRTAlgorithmes.h"
#include "RRTBoxS.h"
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
    cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
    cout << " \t=============================================================" <<endl;
	cout << " \t==========       Impelmentation du :               ==========" <<endl;
    cout << " \t==========              ALGORITHME RRT             ==========" <<endl;
    cout << " \t==========              PARTITION OF BOXS          ==========" <<endl;
	cout << " \t==========              CALCUL STAR DISCREPENCY    ==========" <<endl;
	cout << " \t==========              MESURE OF COVERGE          ==========" <<endl;
	cout << " \t==========              ALGORITHME GRRT            ==========" <<endl;
	cout << " \t==========              CALCUL DISPARITY           ==========" <<endl;
	cout << " \t==========              ERROR ESTIMATION           ==========" <<endl;
	cout << " \t==========                                         ==========" <<endl;
	cout << " \t==========                                         ==========" <<endl;
	cout << " \t==========    Encadrant:        Etudiant:          ==========" <<endl;
    cout << " \t==========  MR.TARIQ NAHAL  &  MR.KODAD YASSINE    ==========" <<endl;
    cout << " \t=============================================================\n\n\n\n" <<endl;
	cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
    cout << " \t      WAIT GENERATION OF RRT FOR VERTEXS...                  " <<endl;  
	cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
	double *xinit;
	cout<<"donner la dimension du Sytéme :";
	int dim;
	cin>>dim;
	cout<<"donner le nombre d'iteration :";
	int k;
	cin>>k;
	xinit=new double[dim];
	for(int j=0;j<dim;j++)xinit[j]=0;
	////////////////////////////////////////////////////////////////////////
	RRTAlgorithmes  A;
	A=RRTAlgorithmes(xinit,0.002,dim,k,100,2,0.01,-3,3,-0.2,0.2);
	//A.GeneralGRRTTLimite();
	//A.GeneralGRRTT();
	A.GeneralRRTT();
	//A.RRTMethode2();
	//cout<<A.stardiscrpancy(k)<<endl;

	FILE* output;
	FILE* outputs;
	FILE* outputst;
	output=fopen("T.dat", "w");
	outputs=fopen("G.dat", "w");
	outputst=fopen("S.dat", "w");
	int kg=0;
	for(int i=0;i<A.ConsultTrees().ConsultNbrVertex();i++)
	{
		//star discrepacy
		fprintf(outputst, " %d ",i+1);
		fprintf(outputst, " %f \n",A.ConsultstardiscrepencyforGRRT()[i]);	

		for(int g=0;g<A.consultexinit().ConsultDim();g++)
		{					
			if(kg!=dim-1)
			{
				//  vertex tree et xgols
				fprintf(output, " %f ",A.ConsultTrees().ConsultVerto(i)[g]);
				if(i<A.ConsultTrees().ConsultNbrVertex()-1)
					fprintf(outputs, " %f ",A.consultexgoals()[i].ConsultVertex()[g]);
				kg++;
			}
			else
			{
				fprintf(output, " %f \n",A.ConsultTrees().ConsultVerto(i)[g]);
				if(i<A.ConsultTrees().ConsultNbrVertex()-1)
				fprintf(outputs, " %f \n",A.consultexgoals()[i].ConsultVertex()[g]);
				kg=0;
			}
		}
	}
/*	///////////////////////////////////
	RRTBoxS b;
	b=RRTBoxS(dim,0.1,0,5);
	b.PartitionOfBoxs();
	b.CalculVolumeS();
	RRTBoxS::box *l;
	l=b.ConsultBoxs();
	FILE* BM;
	FILE* BP;
	BM=fopen("BM.dat", "w");
	BP=fopen("BP.dat", "w");
	int i=0;
	while(l!=NULL)
	{
		fprintf(BM, " %f",l->LowerVertex.ConsultVertex()[0]);
		fprintf(BM, " %f\n",l->LowerVertex.ConsultVertex()[1]);
		fprintf(BP, " %f",l->UpperVertex.ConsultVertex()[0]);
		fprintf(BP, " %f\n",l->UpperVertex.ConsultVertex()[1]);
		cout<<"volume de boite plus "<<i<<" "<<l->volumeBplus<<endl;
		cout<<"volume de boite moins "<<i++<<" "<<l->VolumeBmoins<<endl;
		l=l->suivant;
	}*/
	///////////////for rrt et discrpency////////////////////////////////////////////////////
	//A.GeneralMGRRT();
	//A.GeneralRRTNews();
   /* FILE* output;
	FILE* outputst;
	output=fopen("ATTRE.dat", "w");
	//outputst=fopen("MSTTRE.dat", "w");
	int kg=0;
	for(int i=0;i<A.ConsultTrees().ConsultNbrVertex();i++)
		{
			for(int g=0;g<A.consultexinit().ConsultDim();g++)
			{
				if(kg!=dim-1)
				{
					fprintf(output, " %f ",A.ConsultTrees().ConsultVerto(i)[g]);
					kg++;
				}
				else
				{
					fprintf(output, " %f \n",A.ConsultTrees().ConsultVerto(i)[g]);
					kg=0;
				}
			}
			//fprintf(outputst, " %d ",i+1);
		   // fprintf(outputst, " %f \n",A.ConsultstardiscrepencyforGRRT()[i]);
		}*/
	/*
	for(int i=0;i<1500;i++)
	{
		b.GeneralBoxsForRRT(A,i+1);
		fprintf(output, " %d ",i+1);
		fprintf(output, " %f \n",b.ConsultstarDiscrepancyMax());
		//cout<<b.ConsultstarDiscrepancyMax()<<endl;
	}*/
	
/////////////////////////////////////////////////////////////////////////////////////////////////


	/*RRTBoxS::box *l;
	l=B.ConsultBoxs();
	while(l!=NULL)
	{
		cout<<"("<<l->LowerVertex.ConsultVertex()[0]<<";"<<l->LowerVertex.ConsultVertex()[1]<<")";
		cout<<"-"<<"("<<l->UpperVertex.ConsultVertex()[0]<<";"<<l->UpperVertex.ConsultVertex()[1] <<")"<<endl;
		l=l->suivant;
	}*/
	
	/*
	cout<<"donner la dimension de U :";
	int dimu;
	cin>>dimu; 
	cout<<"donner le nombre de U a generer pseudo-aleatoirement: ";
	int nbreOfU;
	cin>>nbreOfU;
	cout<<"donnez le nombre de points a generer dans l'espace ou bien le nbre d'iteration :";
    int nbrOfiteration;
	cin>>nbrOfiteration;
	cout<<"voulez vous generer le Rapidly Random Exploring Trees :";
	cout<<"si oui taper 1 sinon taper 0"<<endl;
	int rrt;
	cin>>rrt;
	if(rrt==1)
	{
	RRTAlgorithmes  B;
	B=RRTAlgorithmes(xinit,0.02,dim,nbrOfiteration,nbreOfU,dimu,0.01,-3,3,-0.2,0.2);
	B.GeneralRRT();
	
	
	cout<<"voulez vous voir les resultats ? pour oui taper 1 sinon 0"<<endl;
	int l;
	cin>>l;
	if(l==1)
	{
		cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
	cout << " \t                                                             " <<endl;
    cout<<  " \t      WAIT GENERATION OF RRT FOR VERTEXS...                  " <<endl;  
	cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;






		int kg=0;
		for(int i=0;i<B.ConsultTrees().ConsultNbrVertex();i++)
		{ 
			for(int g=0;g<B.consultexinit().ConsultDim();g++)
			{
				if(kg!=dim-1)
				{
					cout<<" \t"<<B.ConsultTrees().ConsultVerto(i)[g]<<" ";
					kg++;
				}
				else 
				{
					cout<<" \t"<<B.ConsultTrees().ConsultVerto(i)[g]<<endl;
					kg=0;
				}
			}
		}
	}

	////////pour ecrire dans un fichier RRT.data
	cout<<"vouler vous Enregistrer les donnes dans le fichier? taper svp 1 si oui,sinon 0 "<<endl;
	int fich;
	cin>>fich;
	if(fich==1)
	{
		FILE* output;
		int kg=0;
		output=fopen("rrt.dat", "w");
		for(int i=0;i<B.ConsultTrees().ConsultNbrVertex();i++)
		{
			for(int g=0;g<B.consultexinit().ConsultDim();g++)
			{
				if(kg!=dim-1)
				{
					fprintf(output, " %f ",B.ConsultTrees().ConsultVerto(i)[g]);
					kg++;
				}
				else
				{
					fprintf(output, " %f \n",B.ConsultTrees().ConsultVerto(i)[g]);
					kg=0;
				}
			}
		}
  //output.close(); 
	} 
}
	else
		cout<<"voulez vous voir les partitions?"<<endl;
	cout<<"dsl en cours.. "; 

	cout << " \t                                                             " <<endl;
    cout << " \t                                                             " <<endl;
	cout<<  " \t                                                     THE END. "<<endl;*/
	/*  random C
	cout<<"la star discrepancy de random c "<<A.stardiscrpancy(100);


	for(int i=0;i<100;i++)
	{
		fprintf(outputs, " %f",A.ConsultXgoals()[i].ConsultVertex()[0]);
		fprintf(outputs, " %f\n",A.ConsultXgoals()[i].ConsultVertex()[1]);
	}*/
/*//pour la sequence de halton 
	FILE* outputs;
	outputs=fopen("H.dat", "w");
	cout<<A.stardiscrpancy(k);
	for(int j=0;j<k;j++)
	{
		fprintf(outputs, " %f",A.consultexgoals()[j].ConsultVertex()[0]);
		fprintf(outputs, " %f\n",A.consultexgoals()[j].ConsultVertex()[1]);
	}*/
system("PAUSE");
return EXIT_SUCCESS;
}