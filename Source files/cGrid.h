#ifndef _CGRID_H_
#define _CGRID_H_

#include <iostream>
#include <vector>
#include "cPoint.h"

enum eGridMode {PROFILE, GRID};	//PROFILE	- operations on profile's nodes
								//GRID		- operations on grid's nodes

class cGrid{	
	std::vector<cPoint> mPoints;			//Container for nodes
	int					mPointsProfile;		//Number of profile's points
	int					mRows;				//Number of rows in the grid
public:
	cGrid() : mPointsProfile(0), mRows(0) {};

	void				addNode(const cPoint & inPoint, eGridMode mode = GRID);		//Adds node to the grid
	int					nodeCount(eGridMode mode = GRID) const;						//Returns number of nodes in the grid
	int					rowCount() const;											//Returns number of rows
	//double			distanceBetween(int N1_i, int N1_j, int N2_i, int N2_j);	//Returns distance between nodes
	void				setRows(int inRows);

	cPoint &				operator()(const int k);
	cPoint					operator()(const int k) const;
	cPoint &				operator()(const int i, const int j); 
	cPoint					operator()(const int i, const int j) const; 	
	friend std::ostream &	operator<<(std::ostream & os, cGrid & inGrid);
};

#endif