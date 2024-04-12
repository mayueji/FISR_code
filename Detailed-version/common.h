// Author: Josiah Manson

#pragma once

#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <cstdio>

#include "vect.h"

using namespace std;

struct OctNode;
struct OctTree;

struct Triangle
{
	vect3f verts[3];
};

struct Mesh
{
	vector<Triangle> triangles;

	void clear()
	{
		triangles.clear();
	}
};

struct OrientedPoint
{
	vect3f norm, pos;
	float ds;
	float depth;
};

struct Index
{
	Index()
	{
		val = 0;
	}
	Index(const int v)
	{
		val = v;
	}
	Index(const unsigned char x_, const unsigned char y_, const unsigned char z_)
	{
		val = x_ | (y_ << 1) | (z_ << 2);
	}

	union
	{
		struct
		{
			unsigned int x : 1;
			unsigned int y : 1;
			unsigned int z : 1;
		};
		unsigned char val;
	};
};

struct OctNode
{
	unsigned char *data;

	static int data_size;
	static int node_count;
	

	static int o_children;
	static int o_coeffs;
	static int o_pt_num; 
	static int o_pts_node; 
	static int o_func_vals;
	static int o_blur_vals;
	static int o_neighbor_cells;
	static int o_neighbors;
	static int o_above_below_zero;
	static void calc_offsets();

	OctNode(){data = NULL;}
	~OctNode(){}

	void allocate();
	void clear();


	template <class T> T* gp(int offset){assert(data!=0); return (T*)(data+offset);} 
	template <class T> T& ga(int offset, Index i){assert(data != NULL); return ((T*)(data+offset))[i.val];} 
	template <class T> T& gv(int offset){assert(data!=0); return *(T*)(data+offset);} 
	

	OctNode& neighbors(int d, int v){return ((OctNode*)(data+o_neighbors))[(d<<1)+v];}
	unsigned char& above_below_zero(){return gv<unsigned char>(o_above_below_zero);}
	int& pts_node(){return gv<int>(o_pts_node);}
	float& blur_vals(int x, int y, int z){return ga<float>(o_blur_vals, Index(x,y,z));}
	float& blur_vals_linear(int i){return gp<float>(o_blur_vals)[i];}
	int& pt_num(int x, int y, int z){return ga<int>(o_pt_num, Index(x,y,z));}
	int& pt_num_linear(int i){return gp<int>(o_pt_num)[i];}
	float& func_vals(int x, int y, int z){return ga<float>(o_func_vals, Index(x,y,z));}
	float& func_vals_linear(int i){return gp<float>(o_func_vals)[i];}
	float func_vals_linear2(int i){return gp<float>(o_func_vals)[i];}
	unsigned char& neighbor_cells(int x, int y, int z){return ga<unsigned char>(o_neighbor_cells, Index(x,y,z));}
	unsigned char& neighbor_cells_linear(int i){return gp<unsigned char>(o_neighbor_cells)[i];}
	float& coeffs(int x, int y, int z){return ga<float>(o_coeffs, Index(x,y,z));}
	float& coeffs_linear(int i){return gp<float>(o_coeffs)[i];}
	OctNode& children(int x, int y, int z){return ga<OctNode>(o_children, Index(x,y,z));}
	OctNode& children_linear(int i){return gp<OctNode>(o_children)[i];}
	OctNode getNeighbor(int x, int y, int z); 
	OctNode getNeighborUnsafe(int x, int y, int z); 

	void removeChildren();
	void detachFromNeighbors();
	
	float getValueAtPoint(vect3f &pos, vect3f mine, vect3f ext); 
	float getValueAtPointhaar(vect3f &pos, vect3f mine, vect3f ext,float val);

	void countPoints(vect3f &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int num);

	void coeffs_haar(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth);
	void coeffs_haar2(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext);
	void changecoeffs_haar(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth);
	void coeffs_haar_restart(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth);
	void normal_haar(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth);
	int getdepth(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth);
	float getstep(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth,float &s1,float &s2);
	void coeffs_daub4(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void coeffs_haar_smooth(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void normal_daub4(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void coeffs_daub4_restart(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void changecoeffs_daub4(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void coeffs_haar_smooth_restart(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void changecoeffs_haar_smooth(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	void coeffs_daub4ds(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth);
	bool prune(int from, int to, int minz, int extz);
	void getValueAtPointdaub4(vect3f &pos,int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth,float &val);
	
};


struct OctTree
{
	void countPoints(vect3f &pt, int num = 1);

	void coeffs_haar(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	void changecoeffs_haar(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	void changecoeffs_daub4(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	void changecoeffs_haar_smooth(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	void coeffs_haar_restart(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	void coeffs_daub4_restart(OrientedPoint &pt, int minDepth = 1, int maxDepthRecur = 100);
	void coeffs_haar_smooth_restart(OrientedPoint &pt, int minDepth = 1, int maxDepthRecur = 100);
	void normal_haar(OrientedPoint &pt, int minDepth = -1, int maxDepthRecur = 100);
	int normal_haar2(OrientedPoint &pt, int N ,int minDepth = -1, int maxDepthRecur = 100);
	int normal_haar3(OrientedPoint &pt, int N ,float av,int minDepth = -1, int maxDepthRecur = 100);
	int normal_daub4_3(OrientedPoint &pt, int N, float av,int minDepth=-1, int maxDepthRecur2=100); 
	int normal_daub4_4(OrientedPoint &pt, int N, float av,int minDepth=-1, int maxDepthRecur2=100); 
	void normal_daub4(OrientedPoint &pt, int minDepth = 1, int maxDepthRecur = 100);
	int normal_daub4_2(OrientedPoint &pt, int N, float av,int minDepth = 1, int maxDepthRecur = 100);
	void coeffs_daub4(OrientedPoint &pt, int minDepth = 1, int maxDepthRecur = 100);
	void coeffs_haar_smooth(OrientedPoint &pt, int minDepth = 1, int maxDepthRecur = 100);
	float getvalue(OrientedPoint &pt);
	void outputvalue();
	void removeHighRes(int zvalue);
	void removeHighRes(OctNode node, int size, int z, int zvalue, int depth);
	void removeUpTo(int zvalue);
	void removeUpTo(OctNode node, const int size, const int z, const int zvalue);
	OctNode root;
	vect3f mine, maxe; 
};
struct TraversalData
{
	OctNode node;
	float value;
	int pt_num;

	unsigned char *neighborCells; 
	float *blur_val; 

	short x,y,z; 
	short cellsize; 

	unsigned char depth;
	bool isCopy;
};

struct Visitor;

struct Traversal
{
	Traversal()
	{
		do_evaluation = false;
	}
	bool do_evaluation;

	void traverse(OctTree &tree, Visitor *v);
	
protected:
	void genNodeChild(TraversalData &ctd, TraversalData &td, int x, int y, int z);

	void traverseNode(Visitor *v, TraversalData &td);
	void traverseFace(Visitor *v, TraversalData td[2], int orientation);
	void traverseEdge(Visitor *v, TraversalData td[2][2], int orientation);
	void traverseVertex(Visitor *v, TraversalData td[2][2][2]);
	
	void recFace0(Visitor *v, TraversalData ch[2][2][2]);
	void recFace1(Visitor *v, TraversalData ch[2][2][2]);
	void recFace2(Visitor *v, TraversalData ch[2][2][2]);

	void recEdge0(Visitor *v, TraversalData ch[2][2][2]);
	void recEdge1(Visitor *v, TraversalData ch[2][2][2]);
	void recEdge2(Visitor *v, TraversalData ch[2][2][2]);
};

struct Visitor
{
	virtual bool onNode(TraversalData &n) {return true;}
	virtual bool onFace(TraversalData td[2], int orientation) {return true;}
	virtual bool onEdge(TraversalData td[2][2], int orientation) {return true;}
	virtual bool onVertex(TraversalData td[2][2][2]) {return true;}
};

struct StreamingVisitor : public Visitor
{
	int z_low, z_high;

	StreamingVisitor(int low, int high);
	StreamingVisitor(int z_border);
	TraversalData *minSize(TraversalData *a, TraversalData *b);
	bool insideOrCrossing(int zl, int zm, int zh);
	virtual bool onNode(TraversalData &td);
	virtual bool onFace(TraversalData td[2], int orientation);
	virtual bool onEdge(TraversalData td[2][2], int orientation);
	virtual bool onVertex(TraversalData td[2][2][2]);
};

struct Process
{
	Process(void (*f)(int from, int to))
	{
		completed = 0;
		finished = false;
		func = f;
	}
	Process()
	{
		completed = 0;
		finished = false;
		func = 0;
	}
	int completed; 
	bool finished; 
	void (*func)(int from, int to); 
};
struct Process2
{
	Process2(void (*f)(int from, int to,vector<float>  &normalnew))
	{
		completed = 0;
		finished = false;
		func = f;
		vector <float>  normalnew={};
	}
	Process2()
	{
		completed = 0;
		finished = false;
		func = 0;
		vector <float>  normalnew={};
	}
	int completed;
	bool finished; 
	void (*func)(int from, int to,vector<float>  &normalnew); 
	vector <float>  normalnew;
};

struct StreamInfo
{
	int cells_total; 
	int cells_mem; 
	int stride; 
};

void reconstruct();
void readtable();
void initMCTable();
void readBBox(vect3f &mine, vect3f &maxe, ifstream &file);
bool readPointUpTononormal(OrientedPoint &pt_read, ifstream &file, const vect3f& center, const float scale_factor, OctTree &tree, int &pts_index, int pts_num, int which_stream, int upto);
bool readPointUpTotrue(OrientedPoint &pt_read, ifstream &file, const vect3f& center, const float scale_factor, OctTree &tree, int &pts_index, int pts_num, int which_stream, int upto);
bool readPointUpTo(OrientedPoint &pt_read, ifstream &file, const vect3f& center, const float scale_factor, OctTree &tree, int &pts_index, int pts_num, int which_stream, int upto);
bool readdsUpTo(OrientedPoint &pt_read, ifstream &file, const vect3f& center, const float scale_factor, OctTree &tree, int &pts_index, int pts_num, int which_stream, int upto);
void readPointsForDisplay(string filename);

void connect_neighbors_first(int from, int to);
void connect_neighbors_second(int from, int to);
void blur_second(int from, int to);
void blur_in_memory();
void eval_daub4(int start, int end);
void eval_daub4change(int start, int end,vector <float> &normalnew);
void eval_haar(int start, int end);
void eval_clear(int start, int end);
void blur_clear(int start, int end);
void eval_haarlimit(int start, int end);
double calc_average();
void create_surface(int from, int to);
void prune(int from, int to);
void colorPlane();
void average_nonstreaming(int from, int to);
double average(int from, int to);
void removeHighResAll(OctNode node, int depth);
void output(int j,vector <float> ::iterator it,int pts_num,int to);
void output3(int j,vector <float> ::iterator it,int pts_num,int to);
void counttrue(int j,vector <float> ::iterator it,int pts_num,int to);
void PCAtrue(int j,vector <float> ::iterator it,int pts_num,int to);
void Large_area_reverse(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float av,int i=10);
void overturn(int j,vector <float> ::iterator it,int pts_num,int to);
int Telescopic_flipping(int j,vector <float> ::iterator it,int pts_num,int to,int i);
float getloss(int j,vector <float> ::iterator it,int pts_num,int to,float average);
void getnewcoeffs(int j,vector <float> ::iterator it,int pts_num,int to);
void overturn_daub4(int j,vector <float> ::iterator it,int pts_num,int to);
int Telescopic_flipping_daub4(int j,vector <float> ::iterator it,int pts_num,int to,float avalue,int i);
void Large_area_reverse_daub4(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float av,int i=10);
void Large_area_reverse_daub4_big(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float av,int i=10);
void overturn_haar_smooth(int j,vector <float> ::iterator it,int pts_num,int to);
int Telescopic_flipping_haar_smooth(int j,vector <float> ::iterator it,int pts_num,int to,float avalue,int i);
void Large_area_reverse_haar_smooth(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float av,int i=10);
void changestep();
void Large_area_reverse_haar_smooth_big(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float av,int i=10);
