#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "vertex.h"
using namespace std;

vertex::vertex()
{
	;
}

vertex::vertex(myvector location)
{
	this->location=location;
	index_HE_towards=-1; //-1表示并没有找到一条边由其出发
}

vertex::~vertex()
{
	;
}
