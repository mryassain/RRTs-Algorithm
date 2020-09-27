#include "StdAfx.h"
#include "RRTVertexe.h"
#include <math.h>
#include <windows.h>

RRTVertexe::RRTVertexe(void)
{
	Vertex = NULL;
	dim = 0;
}

RRTVertexe::RRTVertexe(RRTVertexe &a)
 {
	 this->dim=a.ConsultDim();
	 this->Vertex= new double[a.ConsultDim()];
	 for(int i=0;i<this->dim;i++)
	 {
		 this->Vertex[i]=a.ConsultVertex()[i];
	 }
}
RRTVertexe::RRTVertexe(double *ver,int n)
{
	this->dim=n;
    this->Vertex= new double[dim];
    for(int i=0;i<n;i++)
		this->Vertex[i]=ver[i];
}

RRTVertexe::~RRTVertexe(void)
{
	if(Vertex!=NULL) delete Vertex;
}

int RRTVertexe::ConsultDim(void)  
{
	return dim;
}

double *RRTVertexe::ConsultVertex(void)  
{
	return Vertex;
}
double RRTVertexe::ConsultVertex(int i)
{ 
	if(i<this->dim)return this->Vertex[i];
	else return -1;
}

void RRTVertexe::ModiftVertex(double *ver,int n)
{
	/*if ((dim==0) || (dim!=n))
	{
		if (Vertex) delete[] Vertex;
		Vertex = new double[n];
	}*/
	
	this->dim = n;
	if (Vertex)  
	{
		delete[] Vertex;
	}
	Vertex = new double[dim];
	for (int i=0;i<dim;i++)
	{
		this->Vertex[i] = ver[i];
	}
}

void RRTVertexe::ModifVertex(double element,int ind)
{
	if ((ind>=0) && (ind<dim))
		Vertex[ind] = element;
}

void RRTVertexe::ModifVertex(RRTVertexe & ver)
{
	int i;
	if (ver.ConsultDim()==dim)
	{
		for(i = 0;i<dim;i++)
			Vertex[i] = ver.ConsultVertex()[i];
	}
	else {
		delete[] Vertex;
		Vertex = new double[dim = ver.ConsultDim()];
		for(i = 0;i<dim;i++)
			Vertex[i] = ver.ConsultVertex()[i];
	}
}

bool RRTVertexe::operator==(  RRTVertexe & ver)
{
	if (dim != ver.ConsultDim())
		return false;
	else
		for (int i = 0;i<dim;i++)
			if (Vertex[i]!=ver.ConsultVertex()[i]) return false;

	return true;
}

bool RRTVertexe::operator!=(  RRTVertexe & ver)
{
	return !(this->operator ==(ver));
}

RRTVertexe & RRTVertexe::operator=(  RRTVertexe & ver) 
{
	if (this->operator!=(ver))
		ModifVertex(ver);

	return ver;
}

double RRTVertexe::Metrique(  RRTVertexe & ver)
{
	double s=0,d=0;
	
	if(dim == ver.ConsultDim())
	{	
	    for(int i=0; i<dim;i++)
		{
               d = Vertex[i] - ver.ConsultVertex()[i]; 
	           s += pow(d,2);
		}
		s=sqrt(s);
		return s;
	}

	return -1;
}