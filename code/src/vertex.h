#ifndef _VERTEX_
#define _VERTEX_

#include <math.h>
#include "myvector.h"
class vertex
{
public:
	myvector location;
	myvector normal;
	int index_HE_towards; //指向的一条halfedge的索引
	double area_mixed; //one ring 后的混合型有限面积域

	double u0; //初始热函数值
	double ut; //t时间后的热函数值
	double divergence; //散度值
	double valueOfdis; //距离函数值
	vertex();
	vertex(myvector location);
	~vertex();
};

#endif