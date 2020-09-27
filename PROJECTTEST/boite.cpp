#include "StdAfx.h"
#include "boite.h"
#include "StdAfx.h"
#include "RRTBoxS.h"
#include "RRTTree.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>

boite::boite(int dm,double PasU,double MinEsacedeta,double MaxEpacedeta)
{
	double *tab;
	tab=new double[dim];
	double *tabUpp;
	tabUpp=new double[dim];
	dim=dm;
	PasU=PasU;
	MinEspacedeta=MinEsacedeta;
	MaxEspacedeta=MaxEpacedeta;
	Origine=RRTVertexe();
	VolumeTotal=1;
	for(int i=0;i<dim;i++)
	{
		tabUpp[i]=MaxEspacedeta;
		tab[i]=MinEspacedeta;
		VolumeTotal*=abs(MaxEspacedeta-MinEspacedeta);
	}
	Origine.ModiftVertex(tab,dim);
	this->Boxs[0].LowerVertex.ModiftVertex(tab,dim);
	this->Boxs[0].UpperVertex.ModiftVertex(tabUpp,dim);

}

boite::~boite(void)
{
}


void boite::PartitionOfBoxs()
{
	RRTVertexe uper2,lower2,uper1,lower1;
	uper2=RRTVertexe();
	lower2=RRTVertexe();
	uper1=RRTVertexe();
	lower1=RRTVertexe();

}