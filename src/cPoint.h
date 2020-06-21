#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>

class cPoint{
	double	mX;		//x-coordinate
	double	mY;		//y-coordinate
public:
	cPoint(double inX = 0.0, double inY = 0.0) : mX(inX), mY(inY) {};

	cPoint &				setX(double inX);
	cPoint &				setY(double inY);
	cPoint &				setXY(double inX, double inY);

	double					getX(void) const { return mX; }
	double					getY(void) const { return mY; }

	cPoint					operator+(cPoint & p);								//Adds coordinates of two points (vector behaviour)
	cPoint					operator*(double c);								//Multiplies coordinates by c (vector behaviour)
	friend std::ostream &	operator<<(std::ostream & os, const cPoint & inP);
};


#endif