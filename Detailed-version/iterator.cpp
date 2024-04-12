#include "global.h"
#include "global.h"
#include "matrix.h"
#include "parse.h"
#include <iostream>
#include <fstream>
void output(int j,vector <float> ::iterator it,int pts_num,int to)
{
    char*  const fileNameOut = "output2.xyz";
	ofstream osm(fileNameOut,ios::out);
    ifstream file;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;
	int pts_index2=0;   
	while(readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{		
		pt_read2.norm[0]=*it;
		it++;
		pt_read2.norm[1]=*it;
		it++;
		pt_read2.norm[2]=*it;
		it++;
		osm << pt_read2.pos[0] <<" "<< pt_read2.pos[1] <<" "<< pt_read2.pos[2] <<" "<< pt_read2.norm[0] <<" "<< pt_read2.norm[1] <<" "<< pt_read2.norm[2] << endl;
		pts_index2++;	
	}
}
void output3(int j,vector <float> ::iterator it,int pts_num,int to)
{
    char*  const fileNameOut = "output3.xyz";
	ofstream osm(fileNameOut,ios::out);
    ifstream file;
    file.open(g.pointsfile.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;
	int pts_index2=0;   
	while(readPointUpTo(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{		
		pt_read2.norm[0]=*it;
		it++;
		pt_read2.norm[1]=*it;
		it++;
		pt_read2.norm[2]=*it;
		it++;
		osm << pt_read2.pos[0] <<" "<< pt_read2.pos[1] <<" "<< pt_read2.pos[2] <<" "<< pt_read2.norm[0] <<" "<< pt_read2.norm[1] <<" "<< pt_read2.norm[2] << endl;
		pts_index2++;	
	}
}

void counttrue(int j,vector <float> ::iterator it,int pts_num,int to)
{
    ifstream file;
    file.open(g.pointsfile2.c_str(), ios::binary);
    file.seekg(j,ios::beg);
    int pts_index2=0;
	int Num=0;
	float s1,s2=0;
	float normalcal[3]={0,0,0};
    OrientedPoint pt_read2;
	while (readPointUpTotrue(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			normalcal[0]=*it;
			it++;
			normalcal[1]=*it;
			it++;
			normalcal[2]=*it;
			it++;
			s1=sqrt(pt_read2.norm[0]*pt_read2.norm[0]+pt_read2.norm[1]*pt_read2.norm[1]+pt_read2.norm[2]*pt_read2.norm[2]);
			s2=sqrt(normalcal[0]*normalcal[0]+normalcal[1]*normalcal[1]+normalcal[2]*normalcal[2]);
			if(((normalcal[0]*pt_read2.norm[0]+normalcal[1]*pt_read2.norm[1]+normalcal[2]*pt_read2.norm[2])/s1/s2)>0)
			{
				Num++;
			}
			pts_index2++;
        }      
		cout<<"Normal accuracy is :"<<double(Num/float(pts_index2))<<endl;
		cout<<"Normal accuracy number is :"<<Num<<endl;
}

void PCAtrue(int j,vector <float> ::iterator it,int pts_num,int to)
{
    int pts_index2=0;
	
	vector <float> ::iterator it2=it;
    ifstream file;
    file.open(g.pointsfile2.c_str(), ios::binary);
	if (file.fail())
	{
		return;
	}
    file.seekg(j,ios::beg);
    OrientedPoint pt_read2;
    float normalcal2[3]={0,0,0};
    float s11,s21=0;
	while (readPointUpTotrue(pt_read2, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			normalcal2[0]=*it;
			it++;
			normalcal2[1]=*it;
			it++;
			normalcal2[2]=*it;
			it++;
			s11=sqrt(pt_read2.norm[0]*pt_read2.norm[0]+pt_read2.norm[1]*pt_read2.norm[1]+pt_read2.norm[2]*pt_read2.norm[2]);
			s21=sqrt(normalcal2[0]*normalcal2[0]+normalcal2[1]*normalcal2[1]+normalcal2[2]*normalcal2[2]);
			if(((normalcal2[0]*pt_read2.norm[0]+normalcal2[1]*pt_read2.norm[1]+normalcal2[2]*pt_read2.norm[2])/s11/s21)<0)
			{
				*it2=-normalcal2[0];
				it2++;
				*it2=-normalcal2[1];
				it2++;
				*it2=-normalcal2[2];
				it2++;
			}else
			{
				it2++;
				it2++;
				it2++;

			}
			pts_index2++;
		}

        cout<<"This is the best condition of PCA "<<endl;
} 
void Large_area_reverse(int j,vector <float> ::iterator it,vector <float> ::iterator it2,int pts_num,int from,int to,float avalue,int i)
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
			sigma0=g.tree.normal_haar3(pt_read2, i,avalue);
			if(sigma0<0)
			{
				sigmanum=sigmanum+1;
				g.tree.changecoeffs_haar(pt_read2, g.maxDepthMem);
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
			g.tree.coeffs_haar_restart(pt_read2, g.maxDepthMem);
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
			g.tree.coeffs_haar(pt_read2, g.maxDepthMem);
			pts_index2++;
		}

}


void overturn(int j,vector <float> ::iterator it,int pts_num,int to)
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
			g.tree.coeffs_haar_restart(pt_read2, g.maxDepthMem);
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
			g.tree.coeffs_haar(pt_read2, g.maxDepthMem);
			pts_index2++;
		}

}


int  Telescopic_flipping(int j,vector <float> ::iterator it,int pts_num,int to,int i)
{
    int pts_index2=0;
	vector <float> ::iterator beginit=it;
	vector <float> ::iterator it2=it;
    ifstream file;
	float temp[3];
	int sigma=0;
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
			sigma=g.tree.normal_haar2(pt_read2, i);
			if(sigma<0)
			{
				sigmanum=sigmanum+1;
				g.tree.changecoeffs_haar(pt_read2, g.maxDepthMem);
	
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
			g.tree.coeffs_haar_restart(pt_read2, g.maxDepthMem);
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
			g.tree.coeffs_haar(pt_read2, g.maxDepthMem);
			pts_index2++;
		}
   return sigmanum;

} 

float getloss(int j,vector <float> ::iterator it,int pts_num,int to,float average)
{
	float loss=0.0;
	int pts_index2=0;
    ifstream file;
	OrientedPoint pt_read;
	file.open(g.pointsfile2.c_str(), ios::binary);
	file.seekg(j,ios::beg);
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
		loss=loss+pow(g.tree.getvalue(pt_read)*scale_factor*scale_factor-average,2);
		pts_index2++;	
	}
	return pow(loss/(float)pts_num,0.5);
}



void getnewcoeffs(int j,vector <float> ::iterator it,int pts_num,int to)
{
	float loss=0.0;
	int pts_index2=0;
    ifstream file;
	vector <float> ::iterator beginit=it;
	vector <float> ::iterator it2=it;
	OrientedPoint pt_read;
	file.open(g.pointsfile2.c_str(), ios::binary);
	file.seekg(j,ios::beg);
	file.seekg(j,ios::beg);
	pts_index2=0;
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
	{
			g.tree.coeffs_haar_restart(pt_read, g.maxDepthMem);
			pts_index2++;
	}
	file.seekg(j,ios::beg);
	pts_index2=0;
	it=beginit;
	it2=beginit;
	while (readPointUpTo(pt_read, file, center, scale_factor, g.tree, pts_index2, pts_num, 1, to))
		{
			pt_read.norm[0]=*it;
			it++;
			pt_read.norm[1]=*it;
			it++;
			pt_read.norm[2]=*it;
			it++;
			g.tree.coeffs_haar(pt_read, g.maxDepthMem);
			pts_index2++;

		}

}
