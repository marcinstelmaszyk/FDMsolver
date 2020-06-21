#include <iostream>
#include <iomanip>
#include "cPoint.h"

cPoint & cPoint::setX(double inX){
	mX = inX;

	return *this;
}

cPoint & cPoint::setY(double inY){
	mY = inY;

	return *this;
}

cPoint & cPoint::setXY(double inX, double inY){
	mX = inX;
	mY = inY;

	return *this;
}

cPoint cPoint::operator+(cPoint & p){
	cPoint point;

	point.setXY(mX + p.getX(), mY + p.getY());
	
	return point;
}

cPoint cPoint::operator*(double c){
	cPoint point;

	point.setXY(c*mX, c*mY);

	return point;
}

std::ostream & operator<<(std::ostream & os, const cPoint & inP){
	os << std::showpoint << std::setprecision(6);
	os << "(" << inP.mX << "; " << inP.mY << ")";

	return os;
}