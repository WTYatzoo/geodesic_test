#ifndef _FACE_
#define _FACE_

#include <math.h>
#include "myvector.h"

class face 
{
public:

	int index_HE;
	double area;
	myvector normal;
	int is_non_obtuse; //是不是钝三角形

	myvector gradientOfu; //热函数的梯度
	myvector gradientOfdis; //距离函数的梯度 

	face(){}
	~face(){}
	
	face(int index_HE,double area,int is_non_obtuse,myvector normal):index_HE(index_HE),area(area),is_non_obtuse(is_non_obtuse),normal(normal){}
};
#endif