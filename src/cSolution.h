#ifndef _CSOLUTION_H_
#define _CSOLUTION_H_

#include <vector>
#include "cGrid.h"

enum eSolutionType { POISSON, EIKONAL, BRUTE_FORCE};	//Method of solution
enum eSolutionData { DISTANCE, PHI/*,VARIABLE*/};			//Data type to save


struct sSolution{
	std::vector<double>		d;				//Distance
	std::vector<double>		phi;			//Phi
	//std::vector<double>		variable;		//
};


class cSolution : public cGrid{
	sSolution	mSolution;
	bool		mCaseSolved;		
public:		
	cSolution();
	void		addValue(double inValue, eSolutionData valueType= DISTANCE);
	double		getValue(int k, eSolutionData valueType = DISTANCE);
	double		getValue(int i, int j, eSolutionData valueType = DISTANCE);
	void		SolutionDone();	
};

#endif