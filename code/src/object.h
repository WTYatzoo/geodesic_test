#ifndef _OBJECT_
#define _OBJECT_

#include <vector>
#include "halfedge.h"
#include "vertex.h"
#include "face.h"
#include "testface.h"
#include "isoline.h"
using namespace std;

class object
{
public:
	vector<vertex > myvertexs;
	vector<face > myfaces;
	vector<halfedge > myhalfedges;
	vector<testface > mytestfaces;
        vector<isoline> myisolines;

	int num_vertex;
	int num_face;
	int num_halfedge;

	double disMax,disMin;

	double t; //¿©…¢¡À ±º‰t
	object();
	~object();
	int calData(int index[3],double data[3][3],double &area,myvector &normal);
	void getObjData();
	void calNormal();
	void calGradientOfDis();
        void getIsoline();

	void testdraw();
	void drawwithNormal();
        void drawIsoline();
	
};
#endif
