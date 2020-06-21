#include <iostream>
#include "cGrid.h"


void cGrid::addNode(const cPoint & inPoint, eGridMode mode){
	if(mode == PROFILE)
		++mPointsProfile;

	mPoints.push_back(inPoint);	
}

//Returns number of nodes 
//mode:		eGridMode::PROFILE	- profile
//			eGridMode::GRID		- entire grid
int cGrid::nodeCount(eGridMode mode) const {
	if(mode == PROFILE)
		return mPointsProfile;
	else
		return mPoints.size();
}

//Returns number of rows
int cGrid::rowCount() const {
	return mRows;
}

//double cGrid::distanceBetween(int N1_i, int N1_j, int N2_i, int N2_j){
//	double dx, dy;
//
//	cGrid & grid = *this;	//Shortens  (*this).operator()(N1_i, N1_j)
//
//	dx = grid(N1_i, N1_j).getX() - grid(N2_i, N2_j).getX();
//	dy = grid(N1_i, N1_j).getY() - grid(N2_i, N2_j).getY();
//
//	return sqrt(dx*dx + dy*dy);
//}

void cGrid::setRows(int inRows){
	mRows = inRows;
}


cPoint & cGrid::operator()(const int k){
	return mPoints[k];
}

cPoint cGrid::operator()(const int k) const{
	return mPoints[k];
}

cPoint & cGrid::operator()(const int i, const int j){
	if(i < 0)
		return mPoints[i + (j+1)*nodeCount(PROFILE)];	//Modulo used to make grid continuous e.g. i=N => i=0
	else if(i >= nodeCount(PROFILE))
		return mPoints[i + (j-1)*nodeCount(PROFILE)];	//Modulo used to make grid continuous
	else
		return mPoints[i + j*nodeCount(PROFILE)];
}

cPoint cGrid::operator()(const int i, const int j) const{
	return mPoints[i + j*mPointsProfile];
}

std::ostream & operator<<(std::ostream & os, cGrid & inGrid){
	std::vector<cPoint>::iterator i;
	int j;

	for(i = inGrid.mPoints.begin(), j = 0; i != inGrid.mPoints.end(); ++i, ++j)
		os << "[" << j << "] " << *i << std::endl;

	return os;
}