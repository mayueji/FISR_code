#include "global.h"
#include "global.h"
#include "matrix.h"
#include "parse.h"
#include <iostream>
#include <fstream>
void overturn_daub4(int j,vector <float> ::iterator it,int pts_num,int to)
{
    ifstream file;
    file.open(g.pointsfile2.c_str(), ios::binary);
	vector <float> ::iterator it2=it;
    file.seekg(j,ios::beg);
    int pts_index2=0;
	int Num=0;
    OrientedPoint pt_read2;
		pts_index2=0;
		it=it2;
		file.seekg(j,ios::beg);
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			
			*it=-*it;
			pt_read2.norm[0]=2*(*it);
			it++;
			*it=-*it;
			pt_read2.norm[1]=2*(*it);
			it++;
			*it=-*it;
			pt_read2.norm[2]=2*(*it);
			it++;
			pt_read2.ds=0;
			g.tree.coeffs_daub4(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
		average_val=-average_val;


}

int  Telescopic_flipping_daub4(int j,vector <float> ::iterator it,int pts_num,int to,float avalue,int i)
{
    int pts_index2=0;
	vector <float> ::iterator beginit=it;
	vector <float> ::iterator it2=it;
    ifstream file;
	float temp[3];
	int sigma=1;
	int sigmanum=0;
    file.open(g.pointsfile2.c_str(), ios::binary);
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;
    while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
            sigma=g.tree.normal_daub4_2(pt_read2,i,avalue);
			if(sigma<0)
			{
				pt_read2.ds=0;
				sigmanum=sigmanum+1;
				pt_read2.norm[0]=-2*pt_read2.norm[0];
				pt_read2.norm[1]=-2*pt_read2.norm[1];
				pt_read2.norm[2]=-2*pt_read2.norm[2];
				g.tree.coeffs_daub4(pt_read2, g.maxDepthMem);
			}
            *it2=temp[0]*sigma;
			it2++;
			*it2=temp[1]*sigma;
			it2++;
			*it2=temp[2]*sigma;
			it2++;
			pts_index2++;

		}
		i++;
		cout<< "Telescopic_Flipping:"<< sigmanum <<endl;
   return sigmanum;
} 

float getloss_daub4(int j,vector <float> ::iterator it,int pts_num,int to,float average)
{
	float loss=0.0;
	int pts_index2=0;
    ifstream file;
	OrientedPoint pt_read;
	file.open(g.pointsfile2.c_str(), ios::binary);
	file.seekg(j,ios::beg);
	float pointvalue=0.0;
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		vect3f ext = g.tree.maxe-g.tree.mine;
		g.tree.root.getValueAtPointdaub4(pt_read.pos, 100, 0, g.tree.mine, ext,-1,pointvalue);
		loss=loss+pow(pointvalue*scale_factor*scale_factor-average,2);
		pointvalue=0.0;
		pts_index2++;	
	}
	return pow(loss/(float)pts_num,0.5);
}


void Large_area_reverse_daub4(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float avalue,int i)
{
	int sigma0=1;
	int pts_index2=0;
    ifstream file;
	vector <float> ::iterator beginit=it;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;



	float temp[3]={0,0,0};
	int sigmanum=0;
	
	while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
			sigma0=g.tree.normal_daub4_3(pt_read2, i,avalue);
			if(sigma0<0)
			{
				sigmanum=sigmanum+1;
				pt_read2.ds=0;
				pt_read2.norm[0]=-2*pt_read2.norm[0];
				pt_read2.norm[1]=-2*pt_read2.norm[1];
				pt_read2.norm[2]=-2*pt_read2.norm[2];
				g.tree.coeffs_daub4(pt_read2, g.maxDepthMem);
				*it2=temp[0]*sigma0;
				it2++;
				*it2=temp[1]*sigma0;
				it2++;
				*it2=temp[2]*sigma0;
				it2++;
			}else
			{
				it2++;
				it2++;
				it2++;
			}
			
			pts_index2++;
		}

	cout<<"Large_area_reverse:"<<sigmanum<<endl;
}
void Large_area_reverse_daub4_big(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float avalue,int i)
{
	int sigma0=1;
	int pts_index2=0;
    ifstream file;
	vector <float> ::iterator beginit=it;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;



	float temp[3]={0,0,0};
	int sigmanum=0;
	
	while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
			sigma0=g.tree.normal_daub4_4(pt_read2, i,avalue);
			if(sigma0<0)
			{
				sigmanum=sigmanum+1;
				pt_read2.ds=0;
				pt_read2.norm[0]=-2*pt_read2.norm[0];
				pt_read2.norm[1]=-2*pt_read2.norm[1];
				pt_read2.norm[2]=-2*pt_read2.norm[2];
				g.tree.coeffs_daub4(pt_read2, g.maxDepthMem);
				*it2=temp[0]*sigma0;
				it2++;
				*it2=temp[1]*sigma0;
				it2++;
				*it2=temp[2]*sigma0;
				it2++;
			}else
			{
				it2++;
				it2++;
				it2++;
			}
			
			pts_index2++;
		}

	cout<<"Large_area_reverse_big:"<<sigmanum<<endl;

}
void Large_area_reverse_haar_smooth(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float avalue,int i)
{
	int sigma0=1;
	int pts_index2=0;
    ifstream file;
	vector <float> ::iterator beginit=it;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;



	float temp[3]={0,0,0};
	int sigmanum=0;
	
	while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
			sigma0=g.tree.normal_daub4_3(pt_read2, i,avalue);
			if(sigma0<0)
			{
				sigmanum=sigmanum+1;
				pt_read2.ds=0;
				g.tree.changecoeffs_haar_smooth(pt_read2, g.maxDepthMem);
			}
			*it2=temp[0]*sigma0;
			it2++;
			*it2=temp[1]*sigma0;
			it2++;
			*it2=temp[2]*sigma0;
			it2++;
			pts_index2++;
		}

	cout<<"Large_area_reverse:"<<sigmanum<<endl;

		file.seekg(j,ios::beg);
		pts_index2=0;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.ds = 0;
			g.tree.coeffs_haar_smooth_restart(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
		file.seekg(j,ios::beg);
		pts_index2=0;
		it=beginit;
		it2=beginit;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			
			pt_read2.norm[0]=*it;
			it++;
			pt_read2.norm[1]=*it;
			it++;
			pt_read2.norm[2]=*it;
			it++;
			pt_read2.ds = 0;
			g.tree.coeffs_haar_smooth(pt_read2, g.maxDepthMem);
			pts_index2++;
		}

}

void Large_area_reverse_haar_smooth_big(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float avalue,int i)
{
	int sigma0=1;
	int pts_index2=0;
    ifstream file;
	vector <float> ::iterator beginit=it;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;



	float temp[3]={0,0,0};
	int sigmanum=0;
	
	while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
			sigma0=g.tree.normal_daub4_4(pt_read2, i,avalue);
			if(sigma0<0)
			{
				sigmanum=sigmanum+1;
				pt_read2.ds=0;
				g.tree.changecoeffs_haar_smooth(pt_read2, g.maxDepthMem);
			}
			*it2=temp[0]*sigma0;
			it2++;
			*it2=temp[1]*sigma0;
			it2++;
			*it2=temp[2]*sigma0;
			it2++;
			pts_index2++;
		}
	cout<<"Large_area_reverse_big:"<<sigmanum<<endl;
		file.seekg(j,ios::beg);
		pts_index2=0;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.ds = 0;
			g.tree.coeffs_haar_smooth_restart(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
		file.seekg(j,ios::beg);
		pts_index2=0;
		it=beginit;
		it2=beginit;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			
			pt_read2.norm[0]=*it;
			it++;
			pt_read2.norm[1]=*it;
			it++;
			pt_read2.norm[2]=*it;
			it++;
			pt_read2.ds = 0;
			g.tree.coeffs_haar_smooth(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
}

int  Telescopic_flipping_haar_smooth(int j,vector <float> ::iterator it,int pts_num,int to,float avalue,int i)
{
    int pts_index2=0;
	vector <float> ::iterator beginit=it;
	vector <float> ::iterator it2=it;
    ifstream file;
	float temp[3];
	int sigma=1;
	int sigmanum=0;
    file.open(g.pointsfile2.c_str(), ios::binary);
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;
    while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
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
            sigma=g.tree.normal_daub4_2(pt_read2,i,avalue);
			if(sigma<0)
			{
				pt_read2.ds=0;
				sigmanum=sigmanum+1;
				g.tree.changecoeffs_haar_smooth(pt_read2, g.maxDepthMem);
			}
            *it2=temp[0]*sigma;
			it2++;
			*it2=temp[1]*sigma;
			it2++;
			*it2=temp[2]*sigma;
			it2++;
			pts_index2++;

		}
		i++;
		cout<<sigmanum<<endl;
		file.seekg(j,ios::beg);
		pts_index2=0;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.ds=0;
			g.tree.coeffs_haar_smooth_restart(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
		file.seekg(j,ios::beg);
		pts_index2=0;
		it2=beginit;
		it=beginit;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.norm[0]=*it;
			it++;
			pt_read2.norm[1]=*it;
			it++;
			pt_read2.norm[2]=*it;
			it++;
			pt_read2.ds=0;
			g.tree.coeffs_haar_smooth(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
   return sigmanum;
} 
void overturn_haar_smooth(int j,vector <float> ::iterator it,int pts_num,int to)
{
    ifstream file;
    file.open(g.pointsfile2.c_str(), ios::binary);
	vector <float> ::iterator it2=it;
    file.seekg(j,ios::beg);
    int pts_index2=0;
	int Num=0;
    OrientedPoint pt_read2;
	while (readPointUpTotrue(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			*it=-*it;
			it++;
			*it=-*it;
			it++;
			*it=-*it;
			it++;
			pts_index2++;
		}
		file.seekg(j,ios::beg);
		pts_index2=0;
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			g.tree.coeffs_haar_smooth_restart(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
	
		pts_index2=0;
		it=it2;
		file.seekg(j,ios::beg);
		while (readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read2.norm[0]=*it;
			it++;
			pt_read2.norm[1]=*it;
			it++;
			pt_read2.norm[2]=*it;
			it++;
			pt_read2.ds=0;
			g.tree.coeffs_haar_smooth(pt_read2, g.maxDepthMem);
			pts_index2++;
		}

}
void changestep()
{
	if(g.stepA==0.005&&g.stepB==0.005&&g.stepC==0.005)
	{
		g.stepA=0.003;
		g.stepB=0.003;
		g.stepC=0.003;
	}else if (g.stepA>=0.005&&g.stepB>=0.005&&g.stepC>=0.005)
	{
		g.stepA=g.stepA/2>0.005?g.stepA/2:0.005;
		g.stepB=g.stepB/2>0.005?g.stepB/2:0.005;
		g.stepC=g.stepC/2>0.005?g.stepC/2:0.005;

	}else if(g.stepA==0.003&&g.stepB==0.003&&g.stepC==0.003)
	{
		g.stepA=0.002;
		g.stepB=0.002;
		g.stepC=0.002;

	}else
	{
		
	}
		g.stepA=g.stepA/2;
		g.stepB=g.stepB/2;
		g.stepC=g.stepC/2;
	
}

void OctTree::outputvalue()
{
	vect3f ext = maxe-mine;
	float value=0.0;
	char*  const fileNameOut2 = "outvalue.txt";
	ofstream osm(fileNameOut2,ios::out);
	OrientedPoint ptt;
	int i,j;
	for(i=0;i<400;i++)
	{
		for(j=0;j<400;j++)
		{
			ptt.pos[0]=i*0.0025;
			ptt.pos[1]=j*0.0025;
			ptt.pos[2]=0.5;
			value=0;
			root.getValueAtPointdaub4(ptt.pos, 400, 0, mine, ext, 1,value);
			osm<<value<<endl;
		}
	}
}