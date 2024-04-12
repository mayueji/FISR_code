// Author: Josiah Manson

#include "global.h"
#include "matrix.h"
#include "parse.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <float.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Point_set_3/IO/XYZ.h>
// #include <GL/glut.h>

using namespace std;

extern double average_val;


volatile __int32 pts_num;
vect3f ext; 
float scale_factor;
vect3f center;

ifstream file;
ifstream file2;
ifstream file3;
ifstream file4;
int getnodenumber(OctNode od)
{
	if(od.data==0)
	{
		return 0;
	}
	else
	{
		int number=1;
		for(int i=0;i<8;i++)
		{
			number=number+getnodenumber(od.children_linear(i));
		}

		return number;
	}

}

void count_points_second(int from, int to)
{
	static int pts_index = 0;
	OrientedPoint pt_read;

	// report progress
	printf("read slab %d/%d\n", from, streaminfo.cells_total);
	srand(1010101);
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index, pts_num, 0, to))
	{
		g.tree.countPoints(pt_read.pos);
		pts_index++;
	}
	int number=0;
	prune(from, to);  
}
void count_points_second_change(int from, int to, vector <float> &normalnew)
{
	static int pts_index = 0;
	OrientedPoint pt_read;

	printf("read slab %d/%d\n", from, streaminfo.cells_total);

	// count the points
	srand(1010101);
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index, pts_num, 0, to))
	{
		g.tree.countPoints(pt_read.pos);
		pts_index++;
	}
}


void count_points_first(int from, int to)
{
	static int pts_index = 0;
	OrientedPoint pt_read;

	// report progress
	printf("read slab %d/%d\n", from, streaminfo.cells_total);

	// count the points
	srand(1010101);
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index, pts_num, 0, to))
	{
		g.tree.countPoints(pt_read.pos);
		pts_index++;
	}
	int number=0;
	number=getnodenumber(g.tree.root);
	cout<<"number:"<<number<<endl;
}

void calc_coeffs_daub4_second(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;

	connect_neighbors_second(from, to);

	// calc coeffs
	srand(1010101);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		pt_read.ds = 0;
		g.tree.coeffs_daub4(pt_read, g.maxDepthMem);
		pts_index2++;
	}
}
void calc_coeffs_haar_smooth_second(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;

	connect_neighbors_second(from, to);

	// calc coeffs
	srand(1010101);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		pt_read.ds = 0;
		g.tree.coeffs_haar_smooth(pt_read, g.maxDepthMem);
		pts_index2++;
	}
}
void calc_coeffs_daub4_secondchange(int from, int to,vector <float> &normalnew)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;

	connect_neighbors_second(from, to);

	// calc coeffs
	srand(1010101);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		pt_read.ds = 0;
		g.tree.coeffs_daub4(pt_read, g.maxDepthMem);
		pts_index2++;
	}
}

void calc_coeffs_daub4_first(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;

	connect_neighbors_first(from, to);

	srand(1010101);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		pt_read.ds = 0;
		g.tree.coeffs_daub4(pt_read, 1, g.maxDepthMem);
		pts_index2++;
	}
	g.tree.removeHighRes(to);
}
void normal_haar_third(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;
	OrientedPoint pt_read2;
	int size = 1 << (g.maxDepth + 1);
	int z = 0;
	float avalue=0.0;
	int start=1;
	OctNode n = g.tree.root;
	static int num=to-from;
	 vector <float>  normalnew;
	 vector <float> ::iterator it,it2;
	 it=normalnew.begin();
	 it2=normalnew.begin();
	float s11,s21=0;
	float normalcal2[3]={0,0,0};
	srand(1010101);
	int j;
	j=file3.tellg();
	blur_clear(from,to);
	eval_clear(from,to);
	eval_haar(from,to);
	blur_second(from,to);
	avalue=average(from,to);

	float loss=0.0;


	pts_index2=0;
	file3.seekg(j,ios::beg);
	while (readPointUpTo(pt_read, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		normalnew.push_back(pt_read.norm[0]);
		normalnew.push_back(pt_read.norm[1]);
		normalnew.push_back(pt_read.norm[2]);
		pts_index2++;
	}

	if(avalue<0)
	{
		it=normalnew.begin();
		overturn(j,it,pts_num,to);
		blur_clear(from,to);
		eval_clear(from,to);
		eval_haar(from,to);
		blur_second(from,to);
		avalue=average(from,to);
	}
	
	it=normalnew.begin();
	counttrue(j,it,pts_num,to);

	static int i=0;
	int flag=0;
	int N=0;
	static int sigmanum=1;
	file3.seekg(j,ios::beg);
	pts_index2=0;
	it=normalnew.begin();
	it2=normalnew.begin();
	float temp[3];
	
	int sigma0=1;
	sigmanum=1;		
	while(sigmanum!=0 && i<N)
	{
		if(i<N-1)
		{
		blur_clear(from,to);
		eval_clear(from,to);
		blur_second(from,to);
		avalue=average(from,to);
		it=normalnew.begin();
		it2=normalnew.begin();
		Large_area_reverse(j,it,it2,pts_num,from,to,avalue,i);
		}
	
		file3.seekg(j,ios::beg);
		pts_index2=0;
		it=normalnew.begin();
		it2=normalnew.begin();
		int sigma=1;
		sigmanum=Telescopic_flipping(j,it,pts_num,to,i);
		i++;

	}
	blur_clear(from,to);
	eval_clear(from,to);
	eval_haar(from,to);
	blur_second(from,to);
	avalue=average(from, to);
	it=normalnew.begin();
	output(j,it,pts_num,to);

	file4.seekg(j,ios::beg);
	it=normalnew.begin();
	counttrue(j,it,pts_num,to);

		if(0)
		{
		it=normalnew.begin();
		getnewcoeffs(j,it,pts_num,to);
		}
		blur_clear(from,to);
		eval_clear(from,to);
		eval_haar(from,to);
		blur_second(from,to);
	cout<<"num:"<<pts_num<<endl;
				
}

void normal_daub4_third(int from, int to)
{
	connect_neighbors_second(from, to);
	static int pts_index2 = 0;
	OrientedPoint pt_read;
	OrientedPoint pt_read2;
	int size = 1 << (g.maxDepth + 1);
	int z = 0;
	OctNode n = g.tree.root;
	
	static vector <float>  normalnew;
	vector <float> ::iterator it,it2;
	it=normalnew.begin();
	it2=normalnew.begin();

	int j;
	j=file3.tellg();
	while (readPointUpTo(pt_read, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		normalnew.push_back(pt_read.norm[0]);
		normalnew.push_back(pt_read.norm[1]);
		normalnew.push_back(pt_read.norm[2]);
		pts_index2++;
	}

	if(average_val<0)
	{
		it=normalnew.begin();
		overturn_daub4(j,it,pts_num,to);
	}
	int i=0;
	blur_clear(from,to);
	eval_clear(from,to);
	eval_daub4(from,to);
	blur_second(from,to);
	it=normalnew.begin();
	it2=normalnew.begin();
	Large_area_reverse_daub4(j,it,it2,pts_num,from,to,average_val,i);
	
	static int sigmanum=1;
    int N=0;
	N=g.N;

	while(i<g.iternum*N)
	{
		if((i%N)<N-7)
		{
		blur_clear(from,to);
		eval_clear(from,to);
		eval_daub4(from,to);
		blur_second(from,to);
		it=normalnew.begin();
		it2=normalnew.begin();
		Large_area_reverse_daub4(j,it,it2,pts_num,from,to,average_val,i);
		}
		file3.seekg(j,ios::beg);
		pts_index2=0;
		it=normalnew.begin();
		it2=normalnew.begin();
		sigmanum=Telescopic_flipping_daub4(j,it,pts_num,to,average_val,i);
		i++;
		if(sigmanum==0)
		{
			i=(int(i/N)+1)*N;
		}
		if(i%N==0)
		{
			changestep();
		}
		if(average_val<0)
		{
		it=normalnew.begin();
		overturn_daub4(j,it,pts_num,to);
		}
	}

	blur_clear(from,to);
	eval_clear(from,to);
	eval_daub4(from,to);
	blur_second(from,to);

	it=normalnew.begin();
	output(j,it,pts_num,to);

}

void normal_haar_smooth_third(int from, int to)
{
	connect_neighbors_second(from, to);
	static int pts_index2 = 0;
	OrientedPoint pt_read;
	OrientedPoint pt_read2;
	int size = 1 << (g.maxDepth + 1);
	int z = 0;
	OctNode n = g.tree.root;
	
	static vector <float>  normalnew;
	vector <float> ::iterator it,it2;
	it=normalnew.begin();
	it2=normalnew.begin();
	int j;
	j=file3.tellg();
	float avalue;
	blur_clear(from,to);
	eval_clear(from,to);
	eval_daub4(from,to);
	blur_second(from,to);
	avalue=average(from,to);

	static int flag=0;
	while (readPointUpTo(pt_read, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		normalnew.push_back(pt_read.norm[0]);
		normalnew.push_back(pt_read.norm[1]);
		normalnew.push_back(pt_read.norm[2]);
		pts_index2++;
	}

	if(avalue<0)
	{
		it=normalnew.begin();
		overturn_haar_smooth(j,it,pts_num,to);
		blur_clear(from,to);
		eval_clear(from,to);
		eval_daub4(from,to);
		blur_second(from,to);
		avalue=average(from,to);
		cout<<"the avalue is:"<<avalue<<endl;
	}
	

	it=normalnew.begin();
	counttrue(j,it,pts_num,to);
	int i=0;
	static int sigmanum=1;
	it=normalnew.begin();
	it2=normalnew.begin();
	Large_area_reverse_haar_smooth_big(j,it,it2,pts_num,from,to,avalue,i);


    int N=20;
	while(sigmanum!=0 && i<N)
	{

		if(i<N-8)
		{
		blur_clear(from,to);
		eval_clear(from,to);
		eval_daub4(from,to);
		blur_second(from,to);
		avalue=average(from,to);
		cout<<avalue<<endl;
		it=normalnew.begin();
		it2=normalnew.begin();
		Large_area_reverse_haar_smooth(j,it,it2,pts_num,from,to,avalue,i);
		it=normalnew.begin();
		it2=normalnew.begin();
		Large_area_reverse_haar_smooth_big(j,it,it2,pts_num,from,to,avalue,i);
		}
		file3.seekg(j,ios::beg);
		pts_index2=0;
		it=normalnew.begin();
		it2=normalnew.begin();
		int sigma=1;
		sigmanum=Telescopic_flipping_haar_smooth(j,it,pts_num,to,avalue,i);
		i++;
		if(avalue<0)
		{
		it=normalnew.begin();
		overturn_haar_smooth(j,it,pts_num,to);
		blur_clear(from,to);
		eval_clear(from,to);
		eval_daub4(from,to);
		blur_second(from,to);
		avalue=average(from,to);
		cout<<"the avalue is:"<<avalue<<endl;
		}
	
	}
	blur_clear(from,to);
	eval_clear(from,to);
	eval_daub4(from,to);
	blur_second(from,to);
	avalue=average(from, to);
	cout<<avalue<<endl;
	it=normalnew.begin();
	output(j,it,pts_num,to);
	g.tree.outputvalue();

	file4.seekg(j,ios::beg);
	it=normalnew.begin();
	counttrue(j,it,pts_num,to);
	cout<<"num:"<<pts_num<<endl;

}

void read_normal_daub4(int from, int to, vector <float>  &normalnew)
{
	OrientedPoint pt_read;
	OrientedPoint pt_read2;
	int size = 1 << (g.maxDepth + 1);
	int z = 0;
	OctNode n = g.tree.root;
	vector <float> ::iterator it,it2;
	it=normalnew.begin();
	it2=normalnew.begin();
    int pts_index2=0;
	srand(1010101);
	int j;
	j=file3.tellg();
	static int flag=0;
	while (readPointUpTo(pt_read, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		normalnew.push_back(pt_read.norm[0]);
		normalnew.push_back(pt_read.norm[1]);
		normalnew.push_back(pt_read.norm[2]);
		pts_index2++;
	}
}
void normal_daub4_change(int from, int to,vector <float>  &normalnew)
{
	connect_neighbors_second(from, to);
	static int pts_index2 = 0;
	OrientedPoint pt_read;
	OrientedPoint pt_read2;
	int size = 1 << (g.maxDepth + 1);
	int z = 0;
	OctNode n = g.tree.root;
	int j;
	j=file3.tellg();
	
	vector <float> ::iterator it,it2;
	it=normalnew.begin();
	it2=normalnew.begin();
	srand(1010101);
	int i;
	j=file3.tellg();
	if(1)
	{
		file3.seekg(j,ios::beg);
		pts_index2=0;
		it=normalnew.begin();
		it2=normalnew.begin();
		int sigma=1;
		float temp[3];

		while (readPointUpTo(pt_read2, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.norm[0]=*it;
			temp[0]=*it;
			it++;
			pt_read2.norm[1]=*it;
			temp[1]=*it;
			it++;
			pt_read2.norm[2]=*it;
			temp[2]=*it;
			it++;
			g.tree.normal_daub4(pt_read2, g.maxDepthMem);
			pts_index2++;
			if((temp[0]*pt_read2.norm[0]+temp[1]*pt_read2.norm[1]+temp[2]*pt_read2.norm[2])<0)
			{
				sigma=-1;
			}else
			{
				sigma=1;
			}
			if(i%5<0)
			{
			*it2=pt_read2.norm[0];
			it2++;
			*it2=pt_read2.norm[1];
			it2++;
			*it2=pt_read2.norm[2];
			it2++;
			}
			else
			{
			*it2=temp[0]*sigma;
			it2++;
			*it2=temp[1]*sigma;
			it2++;
			*it2=temp[2]*sigma;
			it2++;
			}
		}		
	}

}




void coeffs_daub4_change(int from, int to,vector <float> &normalnew)
{

			connect_neighbors_second(from, to);
			OrientedPoint pt_read;
			OrientedPoint pt_read2;
			int pts_index2=0;	
			int j;
			j=file3.tellg();
			vector <float> ::iterator it,it2;
			file3.seekg(j,ios::beg);
			it=normalnew.begin();
			it2=normalnew.begin();
			while (readPointUpTo(pt_read2, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
			{
				pt_read2.norm[0]=*it;
				it++;
				pt_read2.norm[1]=*it;
				it++;
				pt_read2.norm[2]=*it;
				it++;
				pt_read2.ds = 0;
				g.tree.coeffs_daub4(pt_read2, 1,g.maxDepthMem);
				pts_index2++;
			}

}


void coeffs_daub4_restart(int from, int to,vector <float>  &normalnew)
{
		connect_neighbors_second(from, to);
		OrientedPoint pt_read;
		OrientedPoint pt_read2;
		int j;
		j=file3.tellg();
		int pts_index2=0;	
		vector <float> ::iterator it,it2;
		file3.seekg(j,ios::beg);
		it=normalnew.begin();
		it2=normalnew.begin();
		while (readPointUpTo(pt_read2, file3, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			g.tree.coeffs_daub4_restart(pt_read2, 1,g.maxDepthMem);
			pts_index2++;
		}

}

void calc_coeffs_haar_second(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;
	ifstream fileds;
	fileds.open("D:\\Surface_Reconstruction\\meshtoxyz\\meshpoint\\results\\spherearea2.txt");
	if(!fileds.is_open())
    {
        cerr << "open fileds file failed" << endl;
        return;
    }
	vector <float>  dsnew;
	vector <float> ::iterator it,it2;
	float ds;
	int k=0;
	while (fileds >> ds)
		{
			dsnew.push_back(ds);
			k++;
		}
	it=dsnew.begin();
	it2=dsnew.begin();
	cout<<"k:"<<k<<endl;


	// calc coeffs
	srand(1010101);
	int j;
	j=file3.tellg();
	file2.seekg(j,ios::beg);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		g.tree.coeffs_haar(pt_read, g.maxDepthMem);
		pts_index2++;
	}
	fileds.close();
}

void calc_coeffs_haar_first(int from, int to)
{
	static int pts_index2 = 0;
	OrientedPoint pt_read;

	// calc coeffs
	srand(1010101);
	while (readPointUpTo(pt_read, file2, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		pt_read.ds = 0;
		g.tree.coeffs_haar(pt_read, -1, g.maxDepthMem);
		pts_index2++;
	}

	g.tree.removeHighRes(to);
}

void count_and_coeffs_haar_first(int from, int to)
{
	count_points_first(from, to);
	prune(from,to); 
	calc_coeffs_haar_first(from, to);
}

void count_and_coeffs_haar_second(int from, int to)
{
	count_points_second(from, to);
	calc_coeffs_haar_second(from, to);
}

void execPipeline(vector<Process> &pipeline)
{
	srand(1010101);
	
	printf("init marching cubes table\n");
	initMCTable();
	file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
		return;
	file2.open(g.pointsfile.c_str(), ios::binary);
	if (file2.fail())
		return;
	file3.open(g.pointsfile.c_str(), ios::binary);
	if (file3.fail())
	{
		return;
	}
	file4.open(g.pointsfile2.c_str(), ios::binary);
	if (file4.fail())
	{
		return;
	}
		

	readBBox(g.tree.mine, g.tree.maxe, file);
	readBBox(g.tree.mine, g.tree.maxe, file2);
	readBBox(g.tree.mine, g.tree.maxe, file3);
	readBBox(g.tree.mine, g.tree.maxe, file4);
	

	ext = g.tree.maxe - g.tree.mine; 
	scale_factor = 1.0 / max(max(ext[0], ext[1]), ext[2]);
	center = (g.tree.mine + g.tree.maxe) * 0.5;
	if (g.wavelet == daub4)
	{
		g.tree.maxe(3,3,3);
		g.tree.mine(-1,-1,-1);
	}
	else if (g.wavelet == haar)
	{
		g.tree.maxe(1,1,1);
		g.tree.mine(0,0,0);
	}
	ext = g.tree.maxe - g.tree.mine;

	readInline(file, pts_num);
	readInline(file2, pts_num);
	readInline(file3, pts_num);
	readInline(file4, pts_num);


	streaminfo.stride=streaminfo.stride*100;
	while (!pipeline[pipeline.size()-1].finished)
	{
		for (int i = 0; i < pipeline.size(); i++)
		{
			Process &proc = pipeline[i];

			if (!proc.finished)
			{
				if (i == 0 || pipeline[i-1].completed > proc.completed+streaminfo.stride || pipeline[i-1].finished)
				{
					proc.func(proc.completed, proc.completed+streaminfo.stride);
					proc.completed += streaminfo.stride;
					if (proc.completed >= streaminfo.cells_total)
						proc.finished = true;
				}
			}
		}
	}
	file.close();
	file2.close();
	file3.close();
	file4.close();
	
}




void execPipelinechange(vector<Process2> &pipeline)
{
	srand(1010101);
	
	printf("[init] init marching cubes table\n");
	initMCTable();

	file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
		return;
	file2.open(g.pointsfile.c_str(), ios::binary);
	if (file2.fail())
		return;
	file3.open(g.pointsfile.c_str(), ios::binary);
	if (file3.fail())
	{
		return;
	}
	file4.open(g.pointsfile2.c_str(), ios::binary);
	if (file4.fail())
	{
		return;
	}
		

	readBBox(g.tree.mine, g.tree.maxe, file);
	readBBox(g.tree.mine, g.tree.maxe, file2);
	readBBox(g.tree.mine, g.tree.maxe, file3);
	readBBox(g.tree.mine, g.tree.maxe, file4);
	

	ext = g.tree.maxe - g.tree.mine; 
	scale_factor = 1.0 / max(max(ext[0], ext[1]), ext[2]);
	center = (g.tree.mine + g.tree.maxe) * 0.5;
	if (g.wavelet == daub4)
	{
		g.tree.maxe(3,3,3);
		g.tree.mine(-1,-1,-1);
	}
	else if (g.wavelet == haar)
	{
		g.tree.maxe(1,1,1);
		g.tree.mine(0,0,0);
	}
	ext = g.tree.maxe - g.tree.mine;

	readInline(file, pts_num);
	readInline(file2, pts_num);
	readInline(file3, pts_num);
	readInline(file4, pts_num);

	// run the pipeline
	while (!pipeline[pipeline.size()-1].finished)
	{
		for (int i = 0; i < pipeline.size(); i++)
		{
			Process2 &proc = pipeline[i];

			if (!proc.finished)
			{
				if (i == 0 || pipeline[i-1].completed > proc.completed+streaminfo.stride || pipeline[i-1].finished)
				{
					proc.func(proc.completed, proc.completed+streaminfo.stride,proc.normalnew);
					proc.completed += streaminfo.stride;
					if (proc.completed >= streaminfo.cells_total)
						proc.finished = true;
				}
			}
		}
	}
	file.close();
	file2.close();
	file3.close();
	file4.close();
	
}

void fill_in_mem(OctNode n, int d = 0)
{
	if (d >= g.maxDepthMem)
		return;

	for (int i = 0; i < 8; i++)
	{
		if (n.children_linear(i).data == 0)
			n.children_linear(i).allocate();
			
		fill_in_mem(n.children_linear(i), d+1);
	}
}

void average_nonstreaming(int start, int end)
{
	 average_val= calc_average();
}
double average(int start, int end)
{
	double average_val = calc_average();
	return average_val;

}

extern int pruned;

void reconstruct()
{
	g.do_eval_pass=true;


	OctNode::calc_offsets();

	// create tree
	OctTree &tree = g.tree;
	tree.root.allocate();

	if (g.wavelet == daub4)
	{
		g.maxDepth = g.depth + 1;
		if (g.do_streaming)
			g.maxDepthMem = (g.maxDepth-2)*.6 + 2.5;
		else
			g.maxDepthMem = 1;
		g.res = 1 << (g.maxDepth + 1); // 2 cells per 'leaf' node
		g.cellscale = 4.0 / (g.res);
		g.celloffset = 1 + .5*g.cellscale;
	}
	else if (g.wavelet == haar)
	{
		g.maxDepth = g.depth - 1;
		if (g.do_streaming)
			g.maxDepthMem = g.maxDepth*.6 + 0.5;
		else
			g.maxDepthMem = -1;
		g.res = 1 << (g.maxDepth + 1); // 2 cells per 'leaf' node
		g.cellscale = 1.0 / (g.res);
		g.celloffset = 0;
	}

	streaminfo.cells_total = g.res; 
	streaminfo.cells_mem = 1 << max(0, g.maxDepthMem);
	streaminfo.stride = streaminfo.cells_total / streaminfo.cells_mem;

	clock_t t_start = clock();
	vector<Process> second_pass;
	readtable();
	if (g.wavelet == haar && !g.do_prune)
	{
		second_pass.push_back(Process(count_and_coeffs_haar_second));
	}
	else
	{
		second_pass.push_back(Process(count_points_second));

		if (g.wavelet == daub4)
		{
			if (g.smooth)
			{
				second_pass.push_back(Process(calc_coeffs_haar_smooth_second));
			}
			else
			{
				second_pass.push_back(Process(calc_coeffs_daub4_second));
			}

		}
			
		else if (g.wavelet == haar)
		{
			second_pass.push_back(Process(calc_coeffs_haar_second));
			
		}
			
	}
	
	if (g.do_eval_pass)
	{
		if (g.wavelet == daub4)
		{
			if (g.smooth)
			{
			second_pass.push_back(Process(eval_daub4));
			second_pass.push_back(Process(normal_haar_smooth_third));
			}else
			{
			second_pass.push_back(Process(normal_daub4_third));
			}
			



		}

		else if (g.wavelet == haar)
		{
			second_pass.push_back(Process(normal_haar_third));
		}

	}
	second_pass.push_back(Process(create_surface));
	execPipeline(second_pass);
	clock_t t_end = clock();
	printf("average function value: %f\n", average_val);
	printf("\033[94m [iso-value] Final Average Function Value: %f \n \033[0m", average_val);
	printf("\033[94m [Timer] Total: %f \n \033[0m", float(t_end - t_start) / CLOCKS_PER_SEC);

}


