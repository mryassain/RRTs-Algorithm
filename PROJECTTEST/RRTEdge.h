#pragma once

class RRTEdge
{
public:
	// ructeurs destructeur
	RRTEdge(void);
	RRTEdge(int ,double ,int ,int );

	~RRTEdge(void);

	//consultation
	
	int ConsultIndiceDest() ;
	int ConsultIndiceOrg() ;
	double ConsultDist() ;
    int ConsultIndofu() ;

	//modification
	void ModifDist(double);
	void ModifIndiceDest(int);
	void ModifIndiceOrg(int);
    void ModifIndofu(int);

	RRTEdge & operator=(  RRTEdge &);

protected:
	int Indofu;
    double Dist; 
	int IndiceDest,IndiceOrg;
};
