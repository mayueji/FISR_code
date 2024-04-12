#include "global.h"

void OctTree::countPoints(vect3f &pt, int num) 
{
	vect3f ext = maxe-mine;
	root.countPoints(pt, g.maxDepth, 0, mine, ext, num);
}
void OctTree::coeffs_daub4(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_daub4(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::coeffs_haar_smooth(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_haar_smooth(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::coeffs_haar(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_haar(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::changecoeffs_haar(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.changecoeffs_haar(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::changecoeffs_daub4(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.changecoeffs_daub4(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::changecoeffs_haar_smooth(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.changecoeffs_haar_smooth(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::coeffs_haar_restart(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_haar_restart(pt, maxDepthRecur, 0, mine, ext, minDepth);
	}
void OctTree::coeffs_daub4_restart(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_daub4_restart(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::coeffs_haar_smooth_restart(OrientedPoint &pt, int minDepth, int maxDepthRecur) 
{
	vect3f ext = maxe-mine;
	root.coeffs_haar_smooth_restart(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
void OctTree::normal_haar(OrientedPoint &pt, int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int maxDepthRecur= ((int(maxDepthRecur2/8))%6)+6;
	root.normal_haar(pt, maxDepthRecur, 0, mine, ext, minDepth);
}
int OctTree::normal_haar2(OrientedPoint &pt, int N, int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int d=-1000;
	int sigma=0;
	minDepth=-1;
	d=root.getdepth(pt, maxDepthRecur2, 0, mine, ext, minDepth);
	float st1=0,st2=0;
	root.getstep(pt, maxDepthRecur2, 0, mine, ext, minDepth,st1,st2);

	d=pt.ds;
	OrientedPoint pt1,pt2;
	int i=0;
	if(d==-1000)
	{
		cout<<"error"<<endl;
	}

	float step=pow(0.5,d)*(1.0+(N%120)/20.0);
	for (i=0;i<3;i++)
	{
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pt1.norm[i]=pt.norm[i];
		pt2.norm[i]=pt.norm[i];
	}
	float pointvalue[3];
	pointvalue[0]=root.getValueAtPointhaar(pt.pos, mine, ext,0);
	pointvalue[1]=root.getValueAtPointhaar(pt1.pos, mine, ext,0);
	pointvalue[2]=root.getValueAtPointhaar(pt2.pos, mine, ext,0);
	while((pointvalue[1]-pointvalue[0])*(pointvalue[2]-pointvalue[0])>=0 && st1 <pow(0.5,d)*10)
	{
		st1=st1+pow(0.5,d);
		st2=st2+pow(0.5,d);
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pointvalue[1]=root.getValueAtPointhaar(pt1.pos, mine, ext,0);
		pointvalue[2]=root.getValueAtPointhaar(pt2.pos, mine, ext,0);

	}

	if(pointvalue[1]>pointvalue[2])
	{
		sigma=-1;
	}
	else
	{
		sigma=1;
	}
 return(sigma);	
}
int OctTree::normal_haar3(OrientedPoint &pt, int N, float av,int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int d=-1000;
	int sigma=0;
	minDepth=-1;
	d=root.getdepth(pt, maxDepthRecur2, 0, mine, ext, minDepth);
	float st1=0,st2=0;
	root.getstep(pt, maxDepthRecur2, 0, mine, ext, minDepth,st1,st2);
	OrientedPoint pt1,pt2;
	int i=0;
	if(d==-1000)
	{
		cout<<"error"<<endl;
	}
	float step=pow(0.5,d);
	for (i=0;i<3;i++)
	{
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pt1.norm[i]=pt.norm[i];
		pt2.norm[i]=pt.norm[i];
	}
	float pointvalue[3];
	pointvalue[0]=root.getValueAtPointhaar(pt.pos, mine, ext,0);
	pointvalue[1]=root.getValueAtPointhaar(pt1.pos, mine, ext,0);
	pointvalue[2]=root.getValueAtPointhaar(pt2.pos, mine, ext,0);
	int flag=0;
	while((st1< pow(0.5,d)*10 &&(pointvalue[1]-av)<0 && (pointvalue[2]-av)<0 && pointvalue[0]<av))
	{
		st1=st1+pow(0.5,d);
		st2=st2+pow(0.5,d);
		pt1.pos[i]=pt.pos[i]+2*step*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-2*step*pt.norm[i];
		pointvalue[1]=root.getValueAtPointhaar(pt1.pos, mine, ext,0);
		pointvalue[2]=root.getValueAtPointhaar(pt2.pos, mine, ext,0);
		flag=flag+1;

	}

	if(flag>1)
	{
		sigma=-1;
	}
	else
	{
		sigma=1;
	}
 return(sigma);

	
}



void OctTree::normal_daub4(OrientedPoint &pt, int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int maxDepthRecur= ((int(maxDepthRecur2/8)+4)%6)+4;
	root.normal_daub4(pt, maxDepthRecur, 0, mine, ext, minDepth);
}

int OctTree::normal_daub4_2(OrientedPoint &pt,int N,float avalue,int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int sigma=0;

	float st1=0,st2=0;
	OrientedPoint pt1,pt2;
	int i=0;

    float step=g.stepC;
	st1=step*2;
	st2=step*2;
	
	for (i=0;i<3;i++)
	{
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pt1.norm[i]=pt.norm[i];
		pt2.norm[i]=pt.norm[i];
	}
	float pointvalue[3]={0.0,0.0,0.0};
	root.getValueAtPointdaub4(pt.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[0]);
	root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
	root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);
	while((pointvalue[1]-pointvalue[0])*(pointvalue[2]-pointvalue[0])>0 && st1 <step*10)
	{
		st1=st1+step;
		st2=st2+step;
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pointvalue[1]=0.0;
		pointvalue[2]=0.0;
		root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
		root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);

	}
	
	if(pointvalue[1]>pointvalue[2])
	{	
		sigma=-1;
	}
	else
	{
		sigma=1;
	}
 return(sigma);	
}
int OctTree::normal_daub4_3(OrientedPoint &pt, int N, float av,int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int sigma=0;
	
	float st1=0,st2=0;

	OrientedPoint pt1,pt2;
	int i=0;
	float step=g.stepA;
	st1=step*2;
	st2=step*2;
	for (i=0;i<3;i++)
	{
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pt1.norm[i]=pt.norm[i];
		pt2.norm[i]=pt.norm[i];
	}
	float pointvalue[3]={0,0,0};
	
	root.getValueAtPointdaub4(pt.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[0]);
	root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
	root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);
	int flag=0;
	while(((st1< step*10 &&(pointvalue[1]-av)<0 && (pointvalue[2]-av)<0))||(st1< step*10 &&(pointvalue[1]-av)>0 && (pointvalue[2]-av)>0))
	{
		st1=st1+1*step;
		st2=st2+1*step;
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pointvalue[1]=0.0;
		pointvalue[2]=0.0;
		root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
		root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);
		flag=flag+1;

	}

	if(flag>2)
	{
		
		sigma=-1;
	}
	else
	{
		sigma=1;
	}
 return(sigma);

	
}

int OctTree::normal_daub4_4(OrientedPoint &pt, int N, float av,int minDepth, int maxDepthRecur2) 
{
	vect3f ext = maxe-mine;
	int d=-1000;
	int sigma=0;
	
	d=pt.ds+3;
	float st1=0,st2=0;

	OrientedPoint pt1,pt2;
	int i=0;
	if(d==-1000)
	{
		cout<<"error"<<endl;
	}
	int level=0;
	level=1;
	float step=g.stepB/level;
	st1=step*2;
	st2=step*2;
	for (i=0;i<3;i++)
	{
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pt1.norm[i]=pt.norm[i];
		pt2.norm[i]=pt.norm[i];
	}
	float pointvalue[3]={0,0,0};
	
	root.getValueAtPointdaub4(pt.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[0]);

	root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
	root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);
	int flag=0;
	while((st1< step*10 &&(pointvalue[1]-av)>0 && (pointvalue[2]-av)>0 && pointvalue[0]>av))
	{
		st1=st1+1*step;
		st2=st2+1*step;
		pt1.pos[i]=pt.pos[i]+st1*pt.norm[i];
		pt2.pos[i]=pt.pos[i]-st2*pt.norm[i];
		pointvalue[1]=0.0;
		pointvalue[2]=0.0;
		root.getValueAtPointdaub4(pt1.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[1]);
		root.getValueAtPointdaub4(pt2.pos, maxDepthRecur2, 0, mine, ext, minDepth=1,pointvalue[2]);
		flag=flag+1;

	}

	if(flag>3)
	{
		
		sigma=-1;
	}
	else
	{
		sigma=1;
	}
 return(sigma);	
}
void OctTree::removeUpTo(int zvalue)
{
	removeUpTo(root, 1 << (g.maxDepth+1), 0, zvalue);
}

void OctTree::removeUpTo(OctNode node, const int size, const int z, const int zvalue)
{
	if (!node.data)
		return;

	const int size2 = size >> 1;

	if (z < zvalue)
	{
		if (z + size <= zvalue)
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (node.children(i,j,0).data)
					{
						node.children(i,j,0).detachFromNeighbors();
						node.children(i,j,0).clear();
					}
				}
			}
		}
		else
		{
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
					removeUpTo(node.children(i,j,0), size2, z, zvalue);
			}
		}
	}

	if (z + size2 < zvalue)
	{
		if (z + size + size2 <= zvalue)
		{
			// remove high children
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (node.children(i,j,1).data)
					{
						node.children(i,j,1).detachFromNeighbors();
						node.children(i,j,1).clear();
						//node.children(i,j,1).data = 0;
					}
				}
			}
		}
		else
		{
			// recurse high
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
					removeUpTo(node.children(i,j,1), size2, z+size2, zvalue);
			}
		}
	}
}


// removal of high resolution nodes
void OctTree::removeHighRes(int zvalue)
{
	removeHighRes(root, 1 << (g.maxDepth+1), 0, zvalue, 0);
}

void OctTree::removeHighRes(OctNode node, int size, int z, int zvalue, int depth)
{
	if (!node.data)
		return;

	if (depth >= g.maxDepthMem)
	{
		if (z + size <= zvalue)
		{
			// remove children
			node.removeChildren();
		}

		return;
	}

	// recurse
	int size2 = size / 2;

	if (z < zvalue)
	{
		// recurse low
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				removeHighRes(node.children(i,j,0), size2, z, zvalue, depth+1);
		}
	}
	
	if (z + size2 < zvalue)
	{
		// recurse high
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				removeHighRes(node.children(i,j,1), size2, z+size2, zvalue, depth+1);
		}
	}
}

void removeHighResAll(OctNode node, int depth)
{
	if (depth >= g.maxDepthMem)
		node.removeChildren();

	// recur
	for (int i = 0; i < 8; i++)
		if (node.children_linear(i).data)
			removeHighResAll(node.children_linear(i), depth+1);
}

float OctTree::getvalue(OrientedPoint &pt) 
{
	vect3f ext = maxe-mine;
	float value=root.getValueAtPoint(pt.pos, mine, ext);
	return value;
}


