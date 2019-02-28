#ifndef _ISOLINE_
#define _ISOLINE_

#include "myvector.h"
class  isoline
{
public:
	myvector begin_loc;
	myvector end_loc;
	isoline(){}
	~isoline() {}
	isoline(myvector begin_loc,myvector end_loc):begin_loc(begin_loc),end_loc(end_loc) {} 
};
#endif