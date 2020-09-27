#pragma once
#include "RRTVertexe.h"
#include "RRTTree.h"
#include "RRTAlgorithmes.h"
#include <windows.h>

class RRTBoxS
{
public:

	//constructeurs 
	RRTBoxS(void);
	RRTBoxS(int dm,double PasU,double MinEsacedeta,double MaxEpacedeta);//,struct Vertexes *p);
	~RRTBoxS(void);

	//structure boxs

	struct box{
		//////////////////////////13
		//int NbrOfVertex;//fai just pour calculer le nombre de vertexs
		int NbrOfVertexAplus;//fai pour discrepency plus
		int NbrOfVertexAmoins;//fai pour discrepency moins
		long double volumeBplus;//fai pour discr
		long double VolumeBmoins;//fai pour discr
		long double MCVmax2;//fait pour discr
		long double MCVmin1;//fait pour discr
		long double MCVmin2;//fait pour discr
		long double MCVmax1;//fait pour discr
		long double DiscrAplus;//fai pour discr
		long double DiscrAMoins;//fai pour discr 
		long double VolumeBox;//fai pour discr
		long double PoidsPartition;//fai pour discr
		////////////////////////////8///////////////////////////////////////////////////////////
		int AstarPlus;// la A star plus
		int AstarMoins;// la A star moins
		int DeltaPlus;//delta plus
		int DeltaMoins;//delta moins
		long double Epsilone;//fonction potentiel pour la borne inferieur
		long double betaGlobal;//faite pour la guidage
		long double betaLocal;//faite pour la guidage
		long double PotentielGenral;//faite pour la guidage
		/////////////////////////////19////////////////////////////////////////////////////////
		long double MCVmax1Plus;//fai
		long double MCVmax2Plus;//fai
		long double MCVmin1Moin;//fai
		long double MCVmin2Moin;//fai
		int NbrOfSequencePlusP;//fai 
		int NbrOfSequenceMoinsP;//fai 
		int NbrOfSequencePlusQ;//fai  
		int NbrOfSequenceMoinsQ;//fai 
		long double MesureOfdisparityPlus;//fai
		long double mesureOfdisparityMoins;//fai
		long double DeltaDispartyPlus;
		long double DeltaDispartyMoins;
		long double EpsiloneForDispartit;
		long double Muc;//fai
		long double Mu0;//fai
		long double Mum;//fai
		long double PotentielGenralForKRRT;
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		RRTVertexe LowerVertex;//fai
		RRTVertexe UpperVertex;//fai
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool boiteLimite;// parametre limité la boite b. 
		bool frontieres;
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		box *suivant;//fai
	};
	
	//consultation 
	int ConsultDim();
	double ConsultPasU();
	double Consultmaxetat(); 
	double ConsultMinetat();
	double ConsultstarDiscrepancyMax();
	double ConsultStarDiscrepancyMin();
	double ConsultVolumTotal();
	RRTVertexe ConsultationMinGRRT();
	RRTVertexe ConsultationMaxGRRT();
	RRTVertexe ConsultationMinKRRT();
	RRTVertexe ConsultationMaxKRRT();
	RRTAlgorithmes ConsultAlgorithmes();
	
	int ConsultNbrOfVertex();
	int ConsultNbrOfVertexBplus();
	int ConsultNbrOfVertexBmoins();
	double ConsultDiscrAplus();
	double ConsultDiscrAMoins();
	double ConsultPoidsPartition();
	double ConsultVolumeBox();
	double ConsultDisparityMax();
	double ConsultDisparityMin();
	double ConsultPoit();
	void LocaliseSystem(double minvertX,double maxvertX,double minvertY,double maxverY);
	

	//Consultation des contenu des boites
	struct box *ConsultBoxs();

	RRTVertexe & ConsultUpperVertexs();
	RRTVertexe & ConsultLowerVertexs();

	//Modification les attributs
	void ModifDim(int d);
	void ModifPoidsEpsilone(double eps);
	void ModifPoidsBetaLocal(double betl);
	void ModifPoidsBetaGolabal(double betg);
	void ModifPasU(double Pasu);
	void ModificationEspaceMax(double maxx);
	void ModificationEspaceMin(double minn);
	void ModificationNbrOfIteration(double minn);
	//affectation des resultats la relation avec la classe  Algorithme
	void AjoutObjetAlgorithm(RRTAlgorithmes A);

	//methodes pour les partions des boites 
	void PartitionOfBoxs();
	void poidOfPartition();
	void ajouterBox(RRTVertexe LowerVertex,RRTVertexe UpperVertex);
	void IsInBox(RRTVertexe a);
	int chercheBoite();
	void BoiteLimites();

	
	//Les methodes pour calculer des volumes de boxs:
	double CalculVolume(RRTVertexe a ,RRTVertexe b);//fai
	void CalculVolumeS();//fai
	void CalculAstarPlusMoins(int kNbreOfvertex);

	//methode pour calculer le nombre de vertaxe pour une boite IL NE FAU PAS OUBLIER QUE je doi coonstruire une structure apartir de la classe RRTtress
	void IsInBox(RRTVertexe a,box &b);
	void CalculMesureOfDiscrpency(int kNbreOfvertex);//
	void CalculMesureOfDiscrpence(RRTTree::Vertexes *l,int kNbreOfvertex);
	void CalculFonctionPotentielBeta();

	void GeneralBoxsForRRT(RRTVertexe a,int kNbreOfvertex);
	void CalculMesureOfDiscrpency(RRTVertexe a,int kNbreOfvertex);
	//calcule de la discreancy pour chaque boite
	//void CalculBoundedStarDiscripancy(box &b);
	void StarDiscrepancancy(::RRTVertexe * xgoals,int k);
	double Star();
	
	// calcule de la star discrepancy bounded
	void StarDiscrepancy();
	void InitialiserBoxs();
	void InitialiserBoxsForRRT();
	void InitialiserBoxsForGRRT();

	//chercher le maximum du fonction potentiel
	void SearchfctPotentiel();
	void SearchfctPotentielforKRRT();
	
	//fonction General pour RRT
	void GeneralBoxsForRRT(RRTAlgorithmes A,int kNbreOfvertex);
	void GeneralBoxsForGRRT(RRTTree::Vertexes *l,int kNbreOfvertex);
	void GeneralBoxsForKRRT(RRTTree::Vertexes *l,int kNbreOfvertex,RRTTree::Vertexes *Tre1,RRTTree::Vertexes *Tre2);
	void CalculMesureOfDiscrpency(int kNbreOfvertex,RRTAlgorithmes A);
	void CalculMesureOfDiscrpence(RRTVertexe a,int kNbreOfvertex);
	void GeneralBoxsForGRRTT(RRTVertexe a,int kNbreOfvertex);
	void GeneralBoxsForRRTTLimit(RRTVertexe a,int kNbreOfvertex);
	void NouvelMesureCoverge();
	void SearchfctPotentielLimte();
	
	//for diparity
	void RandomForSequence(int TailleOfsequenceQ,int TailleOfsequenceP);
	void CalculMesureOfDiparity();
	void CalculDisparityMinMax();
	void IsInBoxDisparityP(RRTVertexe a,box &b);
	void IsInBoxDisparityQ(RRTVertexe a,box &b);
	void CalculMesureOfDiparity(RRTVertexe a,bool v);
	void IsInBoxDisparity(RRTVertexe a,bool v);
	void GeneralDisparity(int TailleOfsequenceQ,int TailleOfsequenceP);
	void CalculMesureOfDiparity(RRTTree::Vertexes *Tre1,RRTTree::Vertexes *Tre2);
	void CalculDiscr(int kNbreOfvertex);
	void CalculMesureOfDiscrpenceGRRT(RRTTree::Vertexes *l,int kNbreOfvertex);
	double arrpos_i(double X, int N);
	void CalculNouvelleMesure();
	/////////////////////////////////////////////////////////////////////////////////
	double Tranc(double d);

protected :
struct box *Boxs;//fai
int dim;//fai
int NbrOfIteration;//le nombre de vertexes dans la structures vertexes
int TailleOfsequenceP;//Fai
int TailleOfsequenceQ;//fai
double PasU;//fai
double VolumeTotal;//fai
double MinEspacedeta;//fai
double MaxEspacedeta;//fai
double StarDiscrepancyMax;//fai
double StarDiscrepancyMin;//fai
double PoiOfpartitionMax;//fai
double AvergStarDiscrepancyMax;//fait
double MesureOfcoverge;//fait
double DisparityMax;//fai
double DisparityMin;//fai
double AverageDIsparity;//fai
double poieps;//fai
double poibetg;//fai
double poibetal;//fai
double PotG;
double *stardiscrepency;
RRTVertexe *SequenceP;//fai avec random c
RRTVertexe *IdealSequenceOfHalton;
RRTVertexe Origine;//fai
RRTVertexe minGRRT;//fai
RRTVertexe maxGRRT;//fai
RRTVertexe minKRRT;//fai
RRTVertexe maxKRRT;//fai
RRTAlgorithmes TreAlgo;//apartir de la classe algorithme on construi un lobjet mm pour le consulter apre
RRTTree::Vertexes Tree1;
RRTTree::Vertexes Tree2;
};
