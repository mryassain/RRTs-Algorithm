#include "StdAfx.h"
#include "RRTAlgorithmes.h"
#include "RRTVertexe.h"
#include "RRTBoxS.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <cstdlib>
using namespace std;

RRTAlgorithmes::RRTAlgorithmes(double *xinit,double tmps,int dim,int ki,int NbrOfU,int dimu,double pas,double min1,double max1,double minu,double maxu)
{
	this->t=tmps;
	this->dimU=dimu;
	this->dimension=dim;
	this->K=ki;
	this->nbrOfU=NbrOfU;
	this->h=pas;
	this->min=min1;
	this->max=max1;
	U=NULL;
	xgoals=NULL;
	//new RRTVertexe[nbrOfU];
    xsolution=NULL;
		//new RRTVertexe[nbrOfU];
	this->minU=minu;
	this->maxU=maxu;
	T=RRTTree();
	Xsol=RRTVertexe();
	xnew=RRTVertexe();
	xnear=RRTVertexe();
	xgoal=RRTVertexe();
	xint=RRTVertexe();
//construction 
	double *tab1,*tab2;
	tab1=new double[dim];
	tab2=new double[dim];
	for(int i=0;i<dim;i++)
	{
		tab1[i]=min1;
		tab2[i]=max1;
	}
////////////////////////////////////////////////////////////::
	maxV=RRTVertexe(tab2,dim);
	minV=RRTVertexe(tab1,dim);
	xint.ModiftVertex(xinit,dimension);
	T.ajoutVertexes(xint);
}
RRTVertexe & RRTAlgorithmes::consultexinit()
{
	RRTVertexe &l=xint;
	return l;
}
RRTAlgorithmes::RRTAlgorithmes(void)
{   U = NULL;  
	xsolution = NULL;
	T = RRTTree();
    Xsol = RRTVertexe();
	xgoal = RRTVertexe();
	xnew = RRTVertexe();
    xnear = RRTVertexe();
    xint = RRTVertexe();                             
	minU = 0;
	maxU = h = min = 0;
	max=K=0;
	nbrOfU=0;
	dimU=0;
	t=dimension=0; 
 }
RRTTree RRTAlgorithmes::ConsultTrees()
{	
	return T;
}

RRTAlgorithmes::~RRTAlgorithmes()
{ /*
	if(xsolution!=NULL)
	{
		delete []xsolution;
	}
	if(U!=NULL)
	{
		delete []U;
	}*/
}
RRTVertexe *RRTAlgorithmes::ConsultUcontrole()
{
	return U;
}
RRTVertexe *RRTAlgorithmes::consultexgoals()
{
	return this->xgoals;
}
RRTVertexe &RRTAlgorithmes::consultmaxV()
{
	return this->maxV;
}
RRTVertexe &RRTAlgorithmes::consultminV()
{
	return this->minV;
}
RRTVertexe &RRTAlgorithmes::consultexgoal()
{ RRTVertexe &l=xgoal;
	return l;
}
RRTVertexe &RRTAlgorithmes::consultexnew()
{
	RRTVertexe &l=xnew;
	return l;
}
RRTVertexe &RRTAlgorithmes::consultexnear()
{
	RRTVertexe &l=xnear;
	return l;
}
RRTVertexe *RRTAlgorithmes::consulteXsolution()
{
	return xsolution;
}

int RRTAlgorithmes::consultNobmreOfItération()
{
	return this->K;
}
double *RRTAlgorithmes::ConsultstardiscrepencyforGRRT()
{
	return stardiscrepency;
}
//modification
void RRTAlgorithmes::ModiftRRTAlgorithmes(RRTAlgorithmes A)
{
	this->T=A.ConsultTrees();
}
void RRTAlgorithmes::ModiftRRTreesOfAlgo(RRTVertexe a)
{
	this->T.ajoutVertexes(a);

}
bool RRTAlgorithmes::AppartienEspaceDetat(RRTVertexe a)
{ 
   bool l=true;
   int  i=0;
   while(i<a.ConsultDim())
    {
		if( (a.ConsultVertex()[i] > this->max ) || (a.ConsultVertex()[i] < this->min) ) 
		{
			l=false;
		}
		i++;
   }
   return l;
}
void RRTAlgorithmes::genereU(double min2u,double max2u,int nbrOfu,int dimu)
{
	double *tab;
	U=new RRTVertexe[nbrOfU];
	tab=new double[dimu];
	for(int j=0;j<nbrOfu;j++)
	{
		for (int i=0;i<dimu;i++)
		{
			tab[i]= min2u + ((max2u - min2u) * (rand () / (double) RAND_MAX));
		}
		U[j].ModiftVertex(tab,this->dimension);
	}

}
void RRTAlgorithmes::Random(double min1 ,double max1,int k)
{
	double *tab;
	tab=new double[dimension];
	for (int i=0;i<this->dimension;i++)
	{
          //tab[i]= min1 + (int) ((double) rand() / (RAND_MAX+1.0) * (max1-min1+1));
		if(i==0)tab[i]= min1 + ((max1 - min1) * (rand () / (double) RAND_MAX));
		else tab[i]= 0 + ((3) * (rand () / (double) RAND_MAX));//derniere simulation
	}
	xgoal.ModiftVertex(tab,this->dimension); 
	this->xgoals[k]=xgoal;
}
void RRTAlgorithmes::Random(RRTVertexe aLower ,RRTVertexe bUpper,int k)
{
	double *tab;
	tab=new double[dimension];
	for (int i=0;i<bUpper.ConsultDim();i++)
	{
		//tab[i]= min1 + (int) ((double) rand() / (RAND_MAX+1.0) * (max1-min1+1));
		tab[i]=aLower.ConsultVertex()[i]+((bUpper.ConsultVertex()[i]-aLower.ConsultVertex()[i])*(rand () / (double) RAND_MAX));
	}
	xgoal.ModiftVertex(tab,this->dimension); 
	this->xgoals[k].ModifVertex(xgoal);
}
// fction pour tester si la vertexe est dans l'arbre ou non
bool RRTAlgorithmes::IsNotHere(RRTVertexe a){
	RRTTree::Vertexes *l;
	l=T.ConsultVertexes();
	while( (l!=NULL)&&((l->vertex)!=a) )
	{
		l=l->suivant;
	}
	if(l==NULL)return 1;
	else return 0;
}
int RRTAlgorithmes::nearfct()
{
	double d=0;
	int indicexnear=0;
	double dist=0;

	RRTTree::Vertexes *p;

	if( (T.ConsultVertexes()->suivant==NULL))
	{
		xnear.ModifVertex(xint);//on peu faire une methode qui cherche l'indice de vertex donnée
		indicexnear=T.ConsultVertexes()->indiceVert;
	}
	else 
	{
		d=xgoal.Metrique(xint/*T.ConsultVertexes()->vertex*/);
		p=T.ConsultVertexes();  
		indicexnear=1;//p=p->suivant;

	   while(p!=NULL)
	   {
		   dist=xgoal.Metrique(p->vertex);
 		   if(d>dist)
		   {
			   d=dist;
			   xnear.ModifVertex(p->vertex);
			   indicexnear=p->indiceVert;
		    }
		   p=p->suivant;
        }
    } 	
  return indicexnear;	 
}
double RRTAlgorithmes::ifct(double r)
{
	double res=((17.76)*r-(103.79)*pow(r,2)+(229.62)*pow(r,3)-(226.3)*pow(r,4)+(83.72)*pow(r,5));
	return(res);
}
double RRTAlgorithmes::fct(double t1,RRTVertexe x,RRTVertexe u,int i)
{	
	double lc=-1/10;
	double rc=-1.5/5;

////////////////////////////////////////////////////////////////////////////////////////////////////
switch (i)
	{
	case 0:
		{
			return (-x.ConsultVertex()[1]+u.ConsultVertex()[0]);
		}
		break;
	case 1:
		{
			return (lc*x.ConsultVertex()[0]+rc*x.ConsultVertex()[1]+(1.5*lc)+u.ConsultVertex()[1]);
		}
		break;
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////


	/*switch (i) {
case 0:
	{
		double C=2;
		return(   ( - ifct(x.ConsultVertex()[0])+ x.ConsultVertex()[1] ) /(C) + u.ConsultVertex()[0] );
	}
	  break;
case 1: 
     {
		 double R=1.5;
		 double E=1.2;
		 double L=5;
		 return ( (E-(R*x.ConsultVertex()[1])-x.ConsultVertex()[0])/(L)+ u.ConsultVertex()[1]);
	 }
	 break;
	}*//*
	switch (i) {
case 0:
	{
		return (u.ConsultVertex()[0]);
	}
	break;
case 1:
	{
		return (u.ConsultVertex()[1]);
	}
	break;
case 2:
	{
		return (u.ConsultVertex()[2]);
	}
	break;
case 3:
	{
		return (u.ConsultVertex()[3]);
	}
	break;
}*/
	/*
	// a définir pour chaque system  //van der pol
    switch (i) {
case 0:
	{
		return (x.ConsultVertex()[1]+u.ConsultVertex()[0]);
	}
	break;
case 1:
	{
		return (-1*(x.ConsultVertex()[0])+(1-(x.ConsultVertex()[0])*(x.ConsultVertex()[0]))*(x.ConsultVertex()[1]) + u.ConsultVertex()[1]);
	}
	break;
	}*/
	/* 
switch (i){	

 case 0:					
      {
    return(2*(x.ConsultVertex()[0])*(1-(x.ConsultVertex()[0])/2)-(x.ConsultVertex()[0])*(x.ConsultVertex()[1]) + u.ConsultVertex()[0]);
      }
       break;
 case 1: 
     {
		 return (3*x.ConsultVertex()[1]*(1-x.ConsultVertex()[1]/3) - 2*x.ConsultVertex()[0]*x.ConsultVertex()[1] + u.ConsultVertex()[1]);
   } 
	 break;
 }
*/

//Tuneel diode

	/*
case 0:
      {
       return  (-ApproxTD(x[0])+x[1])/(C) +u[0];
      }
break;
case 1: 
     {
      return (E-R*x[1]-x[0])/(L)+ u[1];
     }
break;	
*/
// from d/dt 
	/*
 switch (i) {

 case 0:					
       {
		   return -x.ConsultVertex()[0]-4*x.ConsultVertex()[1]+ u.ConsultVertex()[0]  ;
       }
 break;
 case 1: 
      {
		  return 4*x.ConsultVertex()[0]-x.ConsultVertex()[1]+   u.ConsultVertex()[1]  ;
      }
	break;
 
 case 2: 
      {
		  return 0.5*x.ConsultVertex()[2] +  u.ConsultVertex()[2]  ;
      }
 break;	
 }*/
// linear system dim 2
/*  
 case 0:					
      {
      return    u[0]  ;
      }
break;
case 1: 
     {
      return    u[1]  ;
     }
break;

	
/*
case 2: 
     {
      return    u[2]  ;
     }
break;	

case 3: 
     {
      return    u[3]  ;
     }
break;	
case 4: 
     {
      return    u[4]  ;
     }
break;	
*/
  //sys 10
/*
case 0: 
     {
      return  -0.8*x[0]+ x[1]+ u[0]  ;
     }
break;

case 1: 
     {
      return  -0.8*x[1] + x[2]+ u[1]  ;
     }
break;

case 2: 
     {
      return  -0.8*x[2]+x[3]+ u[2]  ;
     }
break;

case 3: 
     {
      return  -0.8*x[3]+ x[4]+ u[3]  ;
     }
break;

case 4: 
     {
      return  -0.8*x[4]+x[5]+  u[4]  ;
     }
break;
case 5: 
     {
      return  -0.8*x[5]+x[6]+  u[5]  ;
     }
break;

case 6: 
     {
      return  -0.8*x[6]+x[7]+  u[6]  ;
     }
break;
case 7: 
     {
      return  -0.8*x[7]+ x[8]+ u[7]  ;
     }
break;
case 8: 
     {
      return  -0.8*x[8]+ x[9]+ u[8]  ;
     }
break;
case 9: 
     {
      return  -0.8*x[9]+ u[9]  ;
     }
break;
*/

//system dim 2
/*
   case 0:					
      {
      return    -1*x[0]-4*x[1] + u[0]  ;
      }
break;
case 1: 
     {
      return   4*x[0] - x[1] + u[1]  ;
     }
break;	

case 2: 
     {
      return   0.5*x[2] + u[2]  ;
     }
break;	
*/
// lorenz system
/*  
   case 0:					
      {
      return   s*( x[1]-x[0]) + u[0]  ;
      }
break;
case 1: 
     {
      return   r*x[0]-x[1] - x[0]*x[2] + u[1]  ;
     }
break;	
case 2: 
     {
      return   x[0]*x[1] - beta*x[2] + u[2]  ;
     }
break;	
*/

// linear system dim 2 totalement contolable
/*
   case 0:					
      {
      return    u[0]  ;
      }
break;
case 1: 
     {
      return   u[1]  ;
     }
break;	
 
*/
  //
/*

      case 0:					
      {
      return    x[1]*x[2] + u[0]  ;
      }
break;
case 1: 
     {
      return   -1*x[0] +x[1] +u[1]  ;
     }
break;	

case 2: 
     {
      return   x[2]+x[0] +u[2]  ;
     }
break;	

case 3: 
     {
      return   x[0]*x[3] +u[2]  ;
     }
break;	
*/
// 4d avec U	

  /*  
 case 0:					
      {
      return    -1*x[1]- x[2] + u[0]  ;
      }
break;
case 1: 
     {
      return   x[0]+a*x[1] +u[1]  ;
     }
break;

case 2: 
     {
      return       b+x[2]*(x[0]-c) + u[2]  ;
     }
break;

case 3: 
     {
      return      -4*x[2] + u[3]  ;
     }
break;

*/
/*
case 0:					
      {
      return    x[1] + u[0]  ;
      }
break;
case 1: 
     {
      return   -8*x[0] +u[1]  ;
     }
break;

case 2: 
     {
      return       x[3] + u[2]  ;
     }
break;

case 3: 
     {
      return      -4*x[2] + u[3]  ;
     }
break;
default:
        {
               
        }
        
 }
	*/
/*switch (i)
	{
	case 0:
		{
			return (-x.ConsultVertex()[0]-(1.9)*x.ConsultVertex()[1]+u.ConsultVertex()[0]);
		}
		break;
	case 1:
		{
			return ((1.9)*x.ConsultVertex()[0]-x.ConsultVertex()[1]+u.ConsultVertex()[1]);
		}
		break;
	}*/
 }
void RRTAlgorithmes::runge4(double t2,RRTVertexe xnewrs,RRTVertexe u)
{
	double h1=h/2.0, k1, k2, k3, k4; 
    double sauv,e;  

	xnew.ModifVertex(xnewrs);
	int dimx = xnew.ConsultDim();

	for (int i=0;i<u.ConsultDim(); i++)
	{
		//step1
		sauv=xnew.ConsultVertex()[i];
		k1=h*(fct(t2,xnew,u,i));
		//step 2
		e=sauv+k1/2.0;
		xnew.ModifVertex(e,i);
		k2=h*(fct(t2+h1,xnew,u,i));
		//step 3
		e=sauv+k2/2.0;
		xnew.ModifVertex(e,i);
		k3=h*(fct(t2+h1,xnew,u,i));
		//step 4
		e=sauv+k3;
		xnew.ModifVertex(e,i);
		k4=h*(fct(t2+h,xnew,u,i));
		//step 5
		e=sauv+(k1+2*k2+2*k3+k4)/6.0;
		xnew.ModifVertex(e,i);
	}
}
void RRTAlgorithmes::krunge4(double t3, RRTVertexe xnear,RRTVertexe u)
{
	/*
	double d=xgoal.Metrique(xnear);//?
    double dist;//essentiel le nouveau point se cree ici
	xnew.ModifVertex(xnear);//2
	do
	{
		runge4(t3,xnew,u);
	 	dist=xnew.Metrique(xgoal);
		t3=t3+h;
		if(d>dist)
		{
			Xsol=xnew;
			d=dist;
		}
		else break;
	}
	while(AppartienEspaceDetat(xnew));*/
	xnew.ModifVertex(xnear);//2
	for(int i=0;i<3;i++)
	{
		runge4(t3,xnew,u);
		t3=t3+h;
		Xsol=xnew;
	}
}
void RRTAlgorithmes::searchSolutionU(double t, RRTVertexe xnear)
{
	xsolution=new RRTVertexe[nbrOfU];
	 for(int i=0;i<nbrOfU;i++)
	 {
		krunge4(t,xnear,U[i]);          	
        xsolution[i]=Xsol;
	 }
}
int RRTAlgorithmes::CONSTRUCTREE(RRTVertexe xgl,RRTVertexe xnr,int indiceXnear)
{
	double d,dist;
	int indisol=0;
    d=xgl.Metrique(xsolution[0]);
	for(int i=1;i<nbrOfU;i++)
	{
		dist=(xgl.Metrique(xsolution[i]));
		if (dist<d)
		{ 
			d=dist;
			indisol=i;
		}
	}
	T.ajoutVertexes(xsolution[indisol]);
	xnew.ModifVertex(xsolution[indisol]);
	/*RRTEdge eg;
    eg=RRTEdge(); 
    eg.ModifDist(d);
	eg.ModifIndiceOrg(indiceXnear);
	eg.ModifIndiceDest(T.ConsultNbrVertex());
	eg.ModifIndofu(indisol); //car indice solution = indice de U solution car le nbre de solution=nbre de U
	T.ajoutEdges(eg);*/
return indisol;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::CONSTRUCTREENewS(RRTVertexe xgl,RRTVertexe xnr,int indiceXnear)
{
	double d,dist;
	int indisol=0;
    d=xgl.Metrique(xsolution[0]);
	for(int i=1;i<nbrOfU;i++)
	{
		dist=(xgl.Metrique(xsolution[i]));
		if (dist<d)
		{ 
			d=dist;
			indisol=i;
		}
	}
	T.ajoutVertexes(xsolution[indisol]);
	searchSolutionU(t,xsolution[indisol]);
	CONSTRUCTREE(xgl,xsolution[indisol],indisol);
	T.ajoutVertexes(xsolution[indisol]);
	/*RRTEdge eg;
    eg=RRTEdge(); 
    eg.ModifDist(d);
	eg.ModifIndiceOrg(indiceXnear);
	eg.ModifIndiceDest(T.ConsultNbrVertex());
	eg.ModifIndofu(indisol); //car indice solution = indice de U solution car le nbre de solution=nbre de U
	T.ajoutEdges(eg);*/
}
void RRTAlgorithmes::GeneralRRTNews()
{
	int indiceXnear=0;
	int i=0;
	genereU(minU,maxU,nbrOfU,dimU);
	for(int i=0 ; i<K ; i++)
	{
		Random(min,max,i);
		indiceXnear=nearfct();
		searchSolutionU(t,xnear);
		CONSTRUCTREENewS(xgoal,xnear,indiceXnear);	   
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void RRTAlgorithmes::GeneralRRT()
{
	int indiceXnear=0;
	int i=0;
	genereU(minU,maxU,nbrOfU,dimU);
	for(int i=0 ; i<K ; i++)
	{
		Random(min,max,i);
		indiceXnear=nearfct();
		searchSolutionU(t,xnear);
		CONSTRUCTREE(xgoal,xnear,indiceXnear);	   
	}
}*/
///////////////////////////////////////////////////////////////////////////////////////////
//////////////////RRRT avec star discrepancy//////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::GeneralRRTT()
{
	xgoals=new RRTVertexe[this->K];
	stardiscrepency=new double[K];
	RRTBoxS B;
	int indiceXnear=0;
	int indsol;
	B=RRTBoxS(this->dimension,0.1,min,max);
	genereU(minU,maxU,nbrOfU,dimU);
	
	for(int i=0 ; i<K ; i++)
	{
		Random(min,max,i);
		indiceXnear=nearfct();//cherchons voisinages
		searchSolutionU(t,xnear);//cherchons la solution
		indsol=CONSTRUCTREE(xgoal,xnear,indiceXnear);//construir l'arbre
		B.GeneralBoxsForRRT(xnew,i+1);//calcul mesure de couverture
		stardiscrepency[i]=B.ConsultstarDiscrepancyMax();//
		//if(i>0 && stardiscrepency[i]>stardiscrepency[i-1])stardiscrepency[i]=stardiscrepency[i-1];
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////
 //ca marche
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////GRRT//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void RRTAlgorithmes::GeneralGRRT()
{
	stardiscrepency=new double[K];
	RRTBoxS B;
	B=RRTBoxS(this->dimension,0.1,min,max);
	int indiceXnear=0;
	int i=0;
	RRTTree::Vertexes *l;
	genereU(minU,maxU,nbrOfU,dimU);//F5
	for(int i=0 ; i<K ; i++)
	{ 
		l=this->ConsultTrees().ConsultVertexes();
		B.GeneralBoxsForGRRT(l,i+1);
		stardiscrepency[i]=B.ConsultstarDiscrepancyMax();
		minV=B.ConsultationMinGRRT();
		maxV=B.ConsultationMaxGRRT();
		Random(minV,maxV,i);//F4
		indiceXnear=nearfct();//F3
		searchSolutionU(t,xnear);//f2
		CONSTRUCTREE(xgoal,xnear,indiceXnear);//f1	   
	}
}*/
/////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::GeneralGRRTT()
{
	////////////////////////////////////////////////
	FILE* b;////////////////////////////////
	b=fopen("B.dat", "w");////////////////////
	////////////////////////////////////////////////
	RRTBoxS B;
	int indsol=0;
	int indiceXnear=0;
	xgoals=new RRTVertexe[this->K];
	stardiscrepency=new double[K];
	double stardisc;
	/////////////////////////////////////////////
	B=RRTBoxS(this->dimension,0.1,min,max);//good
	genereU(minU,maxU,nbrOfU,dimU);//good
	/////////////////////////////////////////////
	for(int i=0 ; i<K ; i++)
	{
		if(i>0)
		{
			B.GeneralBoxsForGRRTT(xnew,i+1);
			stardiscrepency[i]=B.ConsultstarDiscrepancyMax();
			minV=B.ConsultationMinGRRT();
			maxV=B.ConsultationMaxGRRT();
			Random(minV,maxV,i);
			////////////////////////////////////////////////////////////////////////////
			fprintf(b," %f ",minV.ConsultVertex()[0]);
			fprintf(b, " %f\n ",minV.ConsultVertex()[1]);
			fprintf(b," %f ",maxV.ConsultVertex()[0]);
			fprintf(b, " %f\n ",maxV.ConsultVertex()[1]);
			/////////////////////////////////////////////////////////////////////////////
		}
		else Random(this->min,this->max,i);
		//if(i>0 && stardiscrepency[i]>stardiscrepency[i-1])stardiscrepency[i]=stardiscrepency[i-1];	
		indiceXnear=nearfct();//F3
		searchSolutionU(t,xnear);//f2
		indsol=CONSTRUCTREE(xgoal,xnear,indiceXnear);//f1	
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::GeneralMGRRT()
{
	stardiscrepency=new double[K];
	RRTBoxS B;
	double starMGRTT;
	B=RRTBoxS(this->dimension,0.5,min,max);
	int indiceXnear=0;
	int i=0;
	RRTTree::Vertexes *l;
	genereU(minU,maxU,nbrOfU,dimU);//F5
	for(int i=0 ; i<K ; i++)
	{ 
		l=this->ConsultTrees().ConsultVertexes();
		B.GeneralBoxsForGRRT(l,i+1);
		if(i==0)
		{
			starMGRTT=B.ConsultstarDiscrepancyMax();
		}
		else
		{
			if(starMGRTT>B.ConsultstarDiscrepancyMax())
			{
				stardiscrepency[i]=B.ConsultstarDiscrepancyMax();
			}
			else
			{
				stardiscrepency[i]=starMGRTT;
			}
		}
		minV=B.ConsultationMinGRRT();
		maxV=B.ConsultationMaxGRRT();
		Random(minV,maxV,i);//F4
		indiceXnear=nearfct();//F3
		searchSolutionU(t,xnear);//f2
		CONSTRUCTREE(xgoal,xnear,indiceXnear);//f1	   
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::ConstructionTreforFirst()
{
	this->Tre1=RRTTree();
	this->Tre2=RRTTree();
	RRTTree::Vertexes *l;
	RRTTree l1,l2;
	l=T.ConsultVertexes();
	l1=Tre1;
	l2=Tre2;
	int j=0;
	while(l!=NULL)
	{
		if(j%2==0)
		{
			int k=0;
			while(l!=NULL && k<4)
			{
				if(l->suivant==NULL)
				{
					l1.ajoutVertexes(l->vertex);
					break;
				}
				l1.ajoutVertexes(l->vertex);
				l=l->suivant;
				k++;
			}
			j++;
		}
		else
		{
			l2.ajoutVertexes(l->vertex);
			l=l->suivant;
			j++;
		}
		if((l!=NULL)&&(l->suivant==NULL))l1.ajoutVertexes(l->vertex);

	}
}
void RRTAlgorithmes::ConstructionOfTrees(int eta)
{
	this->Tre1=RRTTree();
	this->Tre2=RRTTree();
	RRTTree::Vertexes *l;
	l=T.ConsultVertexes();
	int j=0;
	while(j<eta)
	{
		l=l->suivant;
		j++;
	}
	int i=0;
	while (l!=NULL)
	{
		if(i%2==0)
		{
			Tre1.ajoutVertexes(l->vertex);
		}
		else
		{
			Tre2.ajoutVertexes(l->vertex);
		} 
		i++;
		l=l->suivant;
	}

}
void RRTAlgorithmes::RemplirTri(int i)
{
	if(i%2==0)
	{
		this->Tre2.ajoutVertexes(this->T.ConsultVertexes()->vertex);
	}
	else
	{
		this->Tre1.ajoutVertexes(this->T.ConsultVertexes()->vertex);
	}

}
/*void RRTAlgorithmes::GeneralKRRT()
{
	//////////////////////////////////////////////////////////////////////////////////////
	stardiscrepency=new double[K];////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	int indiceXnear=0;////////////////////////////////////////////////////////////////////
	int i=0;//////////////////////////////////////////////////////////////////////////////
	int h=0;//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	RRTTree::Vertexes *l;/////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	genereU(minU,maxU,nbrOfU,dimU);///////////////////////////////////génére 100 controles
	///////////////////////////////////////////////////////////////////////////////////////
	RRTBoxS B;////////////////////////////////////////////////////////////////////////////
	B=RRTBoxS(this->dimension,0.1,min,max);///////////////////////////////////////////////
	B.PartitionOfBoxs();//////////////////////////////////////////////////////////////////
	B.CalculVolumeS();////////////////////////////////////////////////////////////////////
	B.InitialiserBoxs();///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	l=this->ConsultTrees().ConsultVertexes();//////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for(int i=0 ; i<K ; i++)
	{ 
		////////////////////////////////////////////////////////////////////////////
		B.CalculMesureOfDiscrpence(a,i+1);//////////////
		B.StarDiscrepancy();/////////////////////////////////////    pour la fction potentiel
		B.CalculAstarPlusMoins(i+2);/////////////////////////////          de la star discrépency
		B.CalculFonctionPotentielBeta();/////////////////////////
		//////////////////////////////////////////////////////////////////////////////
		B.CalculMesureOfDiparity(a,i%2);//////////
		B.CalculDisparityMinMax();/////////////////////////////
		B.PotentielDisparity();
		B.CalculNouvelleMesure();//////////////////////////////
		B.SearchfctPotentielforKRRT();//////////////////////////
		////////////////////////////////////////////////////////////////////////////////

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/*if(i<=100)
		{
			B.CalculAstarPlusMoins(i+1);
			B.CalculMesureOfDiscrpence(l,i+1);
			B.StarDiscrepancy();
			B.CalculFonctionPotentielBeta();
			B.SearchfctPotentiel();
			minV=B.ConsultationMinGRRT();
			maxV=B.ConsultationMaxGRRT();
		}
		else 
		{
			if(i==101){ConstructionTreforFirst();//construction des tre1  et tre2 pour la premiere fois
			}
			B.InitialiserBoxs();
			B.CalculAstarPlusMoins(i+1);//pour la mesure de couverture
			B.CalculMesureOfDiscrpence(l,i);//pour la mesure de couverture
			B.StarDiscrepancy();//pour la mesure de couverture
			B.CalculFonctionPotentielBeta();//pour la mesure de couverture
			B.CalculMesureOfDiparity(Tre1.ConsultVertexes(),Tre2.ConsultVertexes());//tre 1 et tre 2
			B.CalculNouvelleMesure();//calcule la nouvelle mesure
			B.SearchfctPotentielforKRRT();//chercher la fonction potentiel qui maximise
			minV=B.ConsultationMaxKRRT();
			maxV=B.ConsultationMaxKRRT();
		}
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		stardiscrepency[i]=B.ConsultstarDiscrepancyMax();//////////////////////////////////////////////////////
		Random(minV,maxV,i);///////////////////////////////////////////////////////////////////////////////////
		indiceXnear=nearfct();/////////////////////////////////////////////////////////////////////////////////
		searchSolutionU(t,xnear);//////////////////////////////////////////////////////////////////////////////
		CONSTRUCTREE(xgoal,xnear,indiceXnear);/////////////////////////////////////////////////////////////////
		RemplirTri(i);//fonction pour ajouter les vertexes au tres aussi//////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////
	}

}*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
RRTVertexe *RRTAlgorithmes::ConsultXgoals()
{
	return this->xgoals;
}
void RRTAlgorithmes::GenereRandomC(int k)
{
	xgoals=new RRTVertexe[k];
	for(int j=0;j<k;j++)
	{
		double *tab;
		tab=new double[dimension];
		for (int i=0;i<this->dimension;i++)
		{
			tab[i]=(rand()/(double)RAND_MAX);
		}
		this->xgoal.ModiftVertex(tab,dimension);
		this->xgoals[j]=xgoal;
	}
}
void RRTAlgorithmes::GenereRandomHalton(int k)
{
	double H,HP;
	double half,digit;
	int I;
	double *halton;
	double *halton2;
	halton=new double[k];
	halton2=new double[k];
	//////////////////////////////////////////////////
	for(int i=0;i<k;i++)
	{
		H=0;
		half=1/2.0;
		I=i+1;
		while(I!=0)
		{
			 digit=I%2;
			 H=H+digit * half;
			 I=(I-digit)/2.0;
			 half=half/2.0;

		}
		halton[i]=H;
	}
	/////////////////////////////////////////////////////
	for(int i=0;i<k;i++)
	{
		HP=0;
		half=1/3.0;
		I=i+1;
		while(I!=0)
		{
			 digit=I%3;
			 HP=HP+digit*half;
			 I=(I-digit)/3.0;
			 half=half/3.0;

		}
		halton2[i]=HP;
	}
	/////////////////////////////////////////////////////////
	double *tab;
	tab=new double[2];
	this->xgoals=new RRTVertexe[k];
	for(int i=0;i<k;i++)
	{
		tab[0]=halton[i];
		tab[1]=halton2[i];
		this->xgoal.ModiftVertex(tab,2);
		this->xgoals[i]=xgoal;
	}
}

double RRTAlgorithmes::stardiscrpancy(int k)
{
	//GenereRandomC(k);
	GenereRandomHalton(k);
	RRTBoxS B;
	B=RRTBoxS(2,0.03,0,1);
	B.PartitionOfBoxs();
	B.CalculVolumeS();
	B.InitialiserBoxs();
	B.StarDiscrepancancy(xgoals,k);
	return B.Star();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::kPrimrunge4(double t3, RRTVertexe xnear,RRTVertexe u)
{

	double d=xgoal.Metrique(xnear);//?
    double dist;//essentiel le nouveau point se cree ici
	xnew.ModifVertex(xnear);//2
	do
	{
		runge4(t3,xnew,u);
	 	dist=xnew.Metrique(xgoal);
		t3=t3+h;
		if(d>dist)
		{
			Xsol=xnew;
			d=dist;
		}
		else break;
	}
	while(AppartienEspaceDetat(xnew));
}
void RRTAlgorithmes::searchSolutionPrimU(double t, RRTVertexe xnear)
{
	xsolution=new RRTVertexe[nbrOfU];
	
	 for(int i=0;i<nbrOfU;i++)
	 {
		kPrimrunge4(t,xnear,U[i]);          	
        xsolution[i]=Xsol;
	 }
}
void RRTAlgorithmes::GeneralGRRTTLimite()
{
	RRTBoxS B;
	int indsol=0;
	int indiceXnear=0;
	xgoals=new RRTVertexe[this->K];
	stardiscrepency=new double[K];
	double stardisc;
	/////////////////////////////////////////////
	B=RRTBoxS(this->dimension,0.1,min,max);//good
	genereU(minU,maxU,nbrOfU,dimU);//good
	/////////////////////////////////////////////
	for(int i=0 ; i<K ; i++)
	{
		if(i>0)
		{
			B.GeneralBoxsForRRTTLimit(xnew,i);
			//stardiscrepency[i]=B.ConsultstarDiscrepancyMax();
			if(i<3000)
			{
				minV=B.ConsultationMinGRRT();
				maxV=B.ConsultationMaxGRRT();
			}
			else
			this->localise();
			Random(minV,maxV,i);
		}
		else Random(this->min,this->max,i);
		//if(i>0 && stardiscrepency[i]>stardiscrepency[i-1])stardiscrepency[i]=stardiscrepency[i-1];	
		indiceXnear=nearfct();//F3
		this->searchSolutionU(t,xnear);
		//searchSolutionPrimU(t,xnear);//f2
		indsol=CONSTRUCTREE(xgoal,xnear,indiceXnear);//f1	
	}
}
void RRTAlgorithmes::localise()
{
	///////////////////////////////////////////// dimension 2////////////////////////////////////////////////
	RRTTree ::Vertexes *t;//////////////////////////////////////////////////////////////////////
	t=T.ConsultVertexes();//////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	double minx=this->T.ConsultVertexes()->vertex.ConsultVertex()[0];
	double maxx=this->T.ConsultVertexes()->vertex.ConsultVertex()[0];
	double miny=this->T.ConsultVertexes()->vertex.ConsultVertex()[1];
	double maxy=this->T.ConsultVertexes()->vertex.ConsultVertex()[1];
	while(t!=NULL)
	{
		if(minx>t->vertex.ConsultVertex()[0])minx=t->vertex.ConsultVertex()[0];
		if(maxx<t->vertex.ConsultVertex()[0])maxx=t->vertex.ConsultVertex()[0];
		if(miny>t->vertex.ConsultVertex()[1])miny=t->vertex.ConsultVertex()[1];
		if(maxy<t->vertex.ConsultVertex()[1])maxy=t->vertex.ConsultVertex()[1];
		///////////////////////////////////////////////////////////////////////
		t=t->suivant;
	}
//////////////////////////////////////////////////////////////////////////////////
/* pour n dimension
	double *min;//////////////////////////////////////////////////////////////////////
	double *max;////////////////////////////////////////////////////////////////////////////
	max=new double[this->dimension];/////////////////////////////////////////////////////////////////:
	min=new double[this->dimension];/////////////////////////////////////////////////////////////
	for(int i=0;i<this->dimension;i++)////////////////////////////////////////////////////////
	{/////////////////////////////////////////////////////////////////////////////////////////////////
		min[i]=this->T.ConsultVertexes()->vertex.ConsultVertex()[i];//////////////////////////////:
		max[i]=this->T.ConsultVertexes()->vertex.ConsultVertex()[i];///////////////////////////////////
	}//////////////////////////////////////////////////////////////////////////////////////////////////
	while(t!=NULL)///////////////////////////////////////////////////////////////////////////////:
	{/////////////////////////////////////////////////////////////////////////////////////////
		for(int i=0;i<this->dimension;i++)///////////////////////////////////////////////////
		{//////////////////////////////////////////////////////////////////////////////////////////////////////:
			if(min[i]>t->vertex.ConsultVertex()[i])min[i]=t->vertex.ConsultVertex()[i];///////////////////////////:
			if(max[i]<t->vertex.ConsultVertex()[i])max[i]=t->vertex.ConsultVertex()[i];///////////////////////////////:
		}/////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////
		t=t->suivant;///////////////////////////////////////////////////////////////////////////////////////////////
	}/////////////////////////////////////////////////////////////////////////////////////////////:
*/
/////////////////////////////////////////////////////////////////////////////
	double *min;
	double *max;
	max=new double[this->dimension];
	min=new double[this->dimension];
	min[0]=minx;
	min[1]=miny;
	max[0]=maxx;
	max[1]=maxy;
/////////////////////////////////////////////////////////////////////////////
	this->minV=RRTVertexe(min,this->dimension);
	this->maxV=RRTVertexe(max,this->dimension);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RRTAlgorithmes::cherchentre(RRTVertexe a,RRTVertexe b,double *tab1,double *tab2)
{
	RRTTree ::Vertexes *t;
	t=T.ConsultVertexes();
	while(t!=NULL)
	{
		if((t->vertex.ConsultVertex()[0]>=a.ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]<=b.ConsultVertex()[0]))
		{
			for(int i=0;i<this->dimension;i++)
			{
				tab1[i]=t->vertex.ConsultVertex()[i];
				tab2[i]=t->vertex.ConsultVertex()[i];
			}
			break;
		}
		else t=t->suivant;

	}
}
void RRTAlgorithmes::localise(int k)
{
	double *min;
	double *max;
	min=new double[this->dimension];
	max=new double[this->dimension];
	int h=0;
	RRTTree ::Vertexes *t;
	for(int i=0;i<k;i++)
	{
		cherchentre(Partiti[h],Partiti[h+1],min,max);
		t=T.ConsultVertexes();
		///////////////////////////////////dimension 2/////////////////////////////////////////////////
		while(t!=NULL)
		{
			if((min[0]>t->vertex.ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]>=this->Partiti[h].ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]<=this->Partiti[h+1].ConsultVertex()[0]))min[0]=t->vertex.ConsultVertex()[0];
			if((max[0]<t->vertex.ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]>=this->Partiti[h].ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]<=this->Partiti[h+1].ConsultVertex()[0]))max[0]=t->vertex.ConsultVertex()[0];
			if((min[1]>t->vertex.ConsultVertex()[1])&&(t->vertex.ConsultVertex()[0]>=this->Partiti[h].ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]<=this->Partiti[h+1].ConsultVertex()[0]))min[1]=t->vertex.ConsultVertex()[1];
			if((max[1]<t->vertex.ConsultVertex()[1])&&(t->vertex.ConsultVertex()[0]>=this->Partiti[h].ConsultVertex()[0])&&(t->vertex.ConsultVertex()[0]<=this->Partiti[h+1].ConsultVertex()[0]))max[1]=t->vertex.ConsultVertex()[1];
			///////////////////////////////////////////////////////////////////////
			t=t->suivant;
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		for(int j=0;j<this->dimension;j++)
		{
			/////////////////////////////////////////
			this->Partiti[h].ModifVertex(min[j],j);
			/////////////////////////////////////////
			this->Partiti[h+1].ModifVertex(max[j],j);
			///////////////////////////////////////////
		}
		//////////////////////////////////////////////////////////////////////////////////////////
		h+=2;
	}

}
void RRTAlgorithmes::Partit(RRTVertexe a,RRTVertexe b,int k)
{
	///////////////aprés localisation
	this->Partiti=new RRTVertexe[2*k];
	//////////////////////////////////
	for(int i=0;i<2*k;i++)
	{
		if(i%2==0){Partiti[i]=a;}
		else{Partiti[i]=b;}
	}
	double hp=(maxV.ConsultVertex()[0]-minV.ConsultVertex()[0])/4.0;
	double s=minV.ConsultVertex()[0];
	for(int i=0;i<2*k;i++)
	{
		this->Partiti[i].ModifVertex(s,0);
		if(i%2==0)s+=hp;
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	this->localise(k);//////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	double *vol,somvol=0;
	vol=new double[k];
	int h=0;
	for(int i=0;i<k;i++)
	{
		vol[i]=CalculVolume(Partiti[h],Partiti[h+1]);//fonction friend 
		somvol+=vol[i];
		h=h+2;
	}
	/////////////////////////////////////nombre des etats  goals utilisé dans chaque boite ///////////////
	PG=new int[k];
	for(int i=0;i<k;i++)
	{
		this->PG[i]=((double)vol[i]/(double)somvol)*(K-3000);
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////

}
void RRTAlgorithmes::RRTMethode2()
{
	RRTBoxS B;
	int indsol=0;
	int indiceXnear=0;
	//////////////////////////////////////////
	int h=0;//////////////////////////////////
	int d=0;//////////////////////////////////
	//////////////////////////////////////////
	xgoals=new RRTVertexe[this->K];
	stardiscrepency=new double[K];
	double stardisc;
	/////////////////////////////////////////////
	B=RRTBoxS(this->dimension,0.1,min,max);//////good
	genereU(minU,maxU,nbrOfU,dimU);//////////////good
	/////////////////////////////////////////////
	for(int i=0;i<K;i++)
	{
		if(i==0){Random(this->min,this->max,i);}
		else
		{
			B.GeneralBoxsForRRTTLimit(xnew,i);
			if(i<3000)
			{
				//stardiscrepency[i]=B.ConsultstarDiscrepancyMax();
				minV=B.ConsultationMinGRRT();
				maxV=B.ConsultationMaxGRRT();
				Random(minV,maxV,i);
			}
			else
			{
				if(i==3000)
				{
					this->localise();
					Partit(minV,maxV,4);
				}
				if(this->PG[d]!=0)
				{
					Random(this->Partiti[h],this->Partiti[h+1],i);
					this->PG[d]--;
					if(this->PG[d]==0)
					{
						if(d<3)
						{
							d++;
							h+=2;
						}
						else 
						{
							break;
						}
						
					}

				}
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//if(i>0 && stardiscrepency[i]>stardiscrepency[i-1])stardiscrepency[i]=stardiscrepency[i-1];
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		indiceXnear=nearfct();///////////////////////////////////////////////////////////////////////////////////////F3
		this->searchSolutionU(t,xnear);////////////////////////////////////////////////////////////////////////////////
		indsol=CONSTRUCTREE(xgoal,xnear,indiceXnear);////////////////////////////////////////////////////////////////f1
	}	
}
double RRTAlgorithmes::CalculVolume(RRTVertexe a ,RRTVertexe b)
{
	////////////////////////////////////////////////////////////////////////////
	double V=1;/////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	for(int i=0;i<a.ConsultDim();i++)
	{
		///////////////////////////////////////////////////////////////////////
		V*=abs(b.ConsultVertex()[i]-a.ConsultVertex()[i]);/////////////////////
		///////////////////////////////////////////////////////////////////////
	}
	////////////////////////////////////////////////////////////////////////////
	return V;///////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////
} 