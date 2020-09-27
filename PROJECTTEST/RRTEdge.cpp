#include "StdAfx.h"
#include "RRTEdge.h"

RRTEdge::RRTEdge(void)
{
	Indofu = 0;
    Dist = 0; 
	IndiceDest = 0;
	IndiceOrg = 0;
}

RRTEdge::RRTEdge(int indu,double distance,int indd,int indo)
{
	Indofu = indu;
    Dist = distance; 
	IndiceDest = indd;
	IndiceOrg = indo;
}

RRTEdge::~RRTEdge(void)
{
}

int RRTEdge::ConsultIndiceDest()
{
	return IndiceDest;
}

int RRTEdge::ConsultIndiceOrg()
{
	return IndiceOrg;
}

double RRTEdge::ConsultDist()
{
	return Dist;
}

int RRTEdge::ConsultIndofu()
{
	return Indofu;
}

void RRTEdge::ModifDist(double distance)
{
	Dist = distance;
}

void RRTEdge::ModifIndiceDest(int indd)
{
	IndiceDest = indd;
}

void RRTEdge::ModifIndiceOrg(int indo)
{
	IndiceOrg = indo;
}

void RRTEdge::ModifIndofu(int indu)
{
	Indofu = indu;
}

RRTEdge & RRTEdge::operator=(RRTEdge & edge)
{
	Indofu = edge.ConsultIndofu();
    Dist = edge.ConsultDist(); 
	IndiceDest = edge.ConsultIndiceDest();
	IndiceOrg = edge.ConsultIndiceOrg();
	
	return edge;
}