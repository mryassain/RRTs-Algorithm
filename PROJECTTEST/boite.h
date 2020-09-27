#pragma once
#include "RRTVertexe.h"
#include "RRTTree.h"
#include "RRTAlgorithmes.h"
#include <windows.h>
#include <vector>
using namespace std;
class boite
{
public:
//constructeurs 
	boite();
	boite(int dm,double PasU,double MinEsacedeta,double MaxEpacedeta);
	~boite(void);
	struct box{
		int NbrOfVertex;//fai just pour calculer le nombre de vertexs
		int NbrOfVertexAplus;//fai pour discrepency plus
		int NbrOfVertexAmoins;//fai pour discrepency moins
		double volumeBplus;//fai pour discr
		double VolumeBmoins;//fai pour discr
		double MCVmax2;//fait pour discr
		double MCVmin1;//fait pour discr
		double MCVmin2;//fait pour discr
		double MCVmax1;//fait pour discr
		double DiscrAplus;//fai pour discr
		double DiscrAMoins;//fai pour discr 
		double VolumeBox;//fai pour discr
		double PoidsPartition;//fai pour discr
		RRTVertexe LowerVertex;//fai
		RRTVertexe UpperVertex;//fai
	};
	void PartitionOfBoxs();

protected :
	std::vector <struct box> Boxs;
	int dim;//fai
    int NbrOfIteration;//fai
	double PasU;//fai
	double VolumeTotal;
	double MinEspacedeta;//fai
	double MaxEspacedeta;//fai
	double StarDiscrepancyMax;
	double StarDiscrepancyMin;
	double PoiOfpartitionMax;
	double AvergStarDiscrepancyMax;
	double MesureOfcoverge;
	double *stardiscrepency;
	RRTVertexe Origine;//fai
};
