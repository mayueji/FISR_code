#include "global.h"
#include <math.h>
#include <fstream>


void OctNode::normal_haar(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth)
{
	vect3f ext2 = ext*.5;
	vect3f center = mine + ext2;
	float val[8];
	int ind [ 3 ] = {0, 0, 0};
	vect3f lp; // local position

	for (int i = 0; i < 3; i++)
	{
		lp[i] = (pt.pos[i] - mine[i]) / ext[i];
		if ( lp [ i ] > 0.5 )
		{
			ind [ i ] = 1;
		}
	}
	if (children(ind[0],ind[1],ind[2]).data == 0 || depth >= maxDepth)
	{
		float n[3];
		for (int child = 0; child < 8; child++)
		{
			val[child]=func_vals_linear(child);
		}
		float s2[4];
		float s[4];
	if(lp[1]<0.5)
	{
			if(lp[2]<0.5)
			{
				n[0]=val[1]-val[0];
			}
			else
			{
				n[0]=val[5]-val[4];
			}
	}
	else
	{
			if(lp[2]<0.5)
			{
				n[0]=val[3]-val[2];
			}
			else
			{
				n[0]=val[7]-val[6];
			}

	}
		if(lp[0]<0.5)
	{
			if(lp[2]<0.5)
			{
				n[1]=val[2]-val[0];
			}
			else
			{
				n[1]=val[6]-val[4];
			}
	}
	else
	{
			if(lp[2]<0.5)
			{
				n[1]=val[3]-val[1];
			}
			else
			{
				n[1]=val[7]-val[5];
			}

	}
		if(lp[0]<0.5)
	{
			if(lp[1]<0.5)
			{
				n[2]=val[4]-val[0];
			}
			else
			{
				n[2]=val[6]-val[2];
			}
	}
	else
	{
			if(lp[1]<0.5)
			{
				n[2]=val[5]-val[1];
			}
			else
			{
				n[2]=val[7]-val[3];
			}

	}

	


		n[0]=(val[7]-val[6])*lp[1]*lp[2]+(val[1]-val[0])*(1-lp[1])*(1-lp[2])+(val[5]-val[4])*(1-lp[1])*lp[2]+(val[3]-val[2])*lp[1]*(1-lp[2]);
		n[1]=(val[7]-val[5])*lp[0]*lp[2]+(val[2]-val[0])*(1-lp[0])*(1-lp[2])+(val[6]-val[4])*(1-lp[0])*lp[2]+(val[3]-val[1])*lp[0]*(1-lp[2]);
		n[2]=(val[7]-val[3])*lp[0]*lp[1]+(val[4]-val[0])*(1-lp[0])*(1-lp[1])+(val[6]-val[2])*(1-lp[0])*lp[1]+(val[5]-val[1])*lp[0]*(1-lp[1]);
		float l2;
		l2=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		if (l2 !=0)
		{
			n[0]=-n[0]/l2;
			n[1]=-n[1]/l2;
			n[2]=-n[2]/l2;
			pt.ds=sqrt(pt.norm[0]*pt.norm[0]+pt.norm[1]*pt.norm[1]+pt.norm[2]*pt.norm[2]);
			int m1,m2;
			m1=1;
			m2=0;
			pt.norm[0]=(n[0]*m1+pt.norm[0]/pt.ds*m2)/(m1+m2);
			pt.norm[1]=(n[1]*m1+pt.norm[1]/pt.ds*m2)/(m1+m2);
			pt.norm[2]=(n[2]*m1+pt.norm[2]/pt.ds*m2)/(m1+m2);

			pt.ds=sqrt(pt.norm[0]*pt.norm[0]+pt.norm[1]*pt.norm[1]+pt.norm[2]*pt.norm[2]);
			pt.norm[0]=pt.norm[0]/pt.ds;
			pt.norm[1]=pt.norm[1]/pt.ds;
			pt.norm[2]=pt.norm[2]/pt.ds;

		
		}
		else
		{
			cout<<"000"<<endl;
		}
	}
	else
	{
		vect3f off(ext2[0]*ind[0], ext2[1]*ind[1], ext2[2]*ind[2]);
		children(ind[0],ind[1],ind[2]).normal_haar(pt, maxDepth, depth+1, mine+off, ext2, minDepth);
	}		


}


void OctNode::normal_daub4(OrientedPoint &pt, int maxDepth, int depth, vect3f &mine, vect3f &ext, int minDepth)
{
	if (g.do_pts_node)
	{
		if (pts_node() == 0)
			return;
	}
	else
	{
		if (pt_num(0,0,0) + pt_num(0,0,1) + pt_num(0,1,0) + pt_num(0,1,1) +
			pt_num(1,0,0) + pt_num(1,0,1) + pt_num(1,1,0) + pt_num(1,1,1) == 0)
			return;
	}
	vect3f ext2 = ext*.5;
	vect3f center = mine + ext2;
	int flag=0;
	int ix = 0, iy = 0, iz = 0;

	if (pt.pos[0] > center[0])
		ix = 1;
	if (pt.pos[1] > center[1])
		iy = 1;
	if (pt.pos[2] > center[2])
		iz = 1;
	float val[8]={0,0,0,0,0,0,0,0};
	// push down
	if (children(ix,iy,iz).data)
	{
		vect3f corner2(mine[0] + ext2[0]*ix, mine[1] + ext2[1]*iy, mine[2] + ext2[2]*iz);
		children(ix,iy,iz).normal_daub4(pt, maxDepth, depth+1, corner2, ext2, minDepth);
		flag=1;
	}
	if (flag==0||depth == maxDepth)
	{
		vect3f lp; // local position
		for (int i = 0; i < 3; i++)
			lp[i] = (pt.pos[i] - mine[i]) / ext[i];

		// add to coefficients
		float factor = (1 << (depth << 1))*pt.ds;
		bool do_zero = (depth == 2);
		for (int child = 0; child < 8; child++)
		{
			int x = child&1;
			int y = (child>>1)&1;
			int z = child>>2;
			val[child]=func_vals_linear(child);;
	

		}
	
		float n[3];

		n[0]=(val[7]-val[6])*lp[1]*lp[2]+(val[1]-val[0])*(1-lp[1])*(1-lp[2])+(val[5]-val[4])*(1-lp[1])*lp[2]+(val[3]-val[2])*lp[1]*(1-lp[2]);
		n[1]=(val[7]-val[5])*lp[0]*lp[2]+(val[2]-val[0])*(1-lp[0])*(1-lp[2])+(val[6]-val[4])*(1-lp[0])*lp[2]+(val[3]-val[1])*lp[0]*(1-lp[2]);
		n[2]=(val[7]-val[3])*lp[0]*lp[1]+(val[4]-val[0])*(1-lp[0])*(1-lp[1])+(val[6]-val[2])*(1-lp[0])*lp[1]+(val[5]-val[1])*lp[0]*(1-lp[1]);


		float l2;
		l2=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
		
		if (l2 !=0)
		{
			n[0]=-n[0]/l2;
			n[1]=-n[1]/l2;
			n[2]=-n[2]/l2;
			pt.ds=sqrt(pt.norm[0]*pt.norm[0]+pt.norm[1]*pt.norm[1]+pt.norm[2]*pt.norm[2]);

			int m1,m2;
			m1=1;
			m2=0;
			pt.norm[0]=(n[0]*m1+pt.norm[0]/pt.ds*m2)/(m1+m2);
			pt.norm[1]=(n[1]*m1+pt.norm[1]/pt.ds*m2)/(m1+m2);
			pt.norm[2]=(n[2]*m1+pt.norm[2]/pt.ds*m2)/(m1+m2);

			pt.ds=sqrt(pt.norm[0]*pt.norm[0]+pt.norm[1]*pt.norm[1]+pt.norm[2]*pt.norm[2]);
			pt.norm[0]=pt.norm[0]/pt.ds;
			pt.norm[1]=pt.norm[1]/pt.ds;
			pt.norm[2]=pt.norm[2]/pt.ds;
		}
	}	
}


int OctNode::getdepth(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth)
{
	vect3f ext2 = ext*.5;
	vect3f center = mine + ext2;
	int d=depth;

	if ( depth > minDepth )
	{
		int ind [ 3 ] = {0, 0, 0};
		vect3f lp; // local position

		for (int i = 0; i < 3; i++)
		{
			lp[i] = (pt.pos[i] - mine[i]) / ext[i];
			if ( lp [ i ] > 0.5 )
			{
				ind [ i ] = 1;
			}
		}
		if ((children(ind[0],ind[1],ind[2]).data == 0)||(pt_num(ind[0],ind[1],ind[2])<=2 && pt_num(ind[0],ind[1],ind[2])>0))
		{
			pt.depth=d;

		}
		else
		{
			vect3f off(ext2[0]*ind[0], ext2[1]*ind[1], ext2[2]*ind[2]);
			d=children(ind[0],ind[1],ind[2]).getdepth(pt, maxDepth, depth+1, mine+off, ext2, minDepth);
		}
		
	}
	return(d);
	
}
float OctNode::getstep(OrientedPoint &pt, int maxDepth, int depth, vect3f mine, vect3f &ext, int minDepth,float &s1,float &s2)
{
	vect3f ext2 = ext*.5;
	vect3f center = mine + ext2;
	int d=depth;
	float dis=0.25;

	if ( depth > minDepth )
	{
		int ind [ 3 ] = {0, 0, 0};
		vect3f lp; 

		for (int i = 0; i < 3; i++)
		{
			lp[i] = (pt.pos[i] - mine[i]) / ext[i];
			if ( lp [ i ] > 0.5 )
			{
				ind [ i ] = 1;
			}
		}
		if ((children(ind[0],ind[1],ind[2]).data == 0)||(pt_num(ind[0],ind[1],ind[2])<2 && pt_num(ind[0],ind[1],ind[2])>0))
		{
			float forward[3]={0,0,0};
			float back[3]={0,0,0};
			if(g.do_pts_node)
			{
				for(int i=0;i<3;i++)
				{
					forward[i]=fabs((pt.norm[i]>0)? (1+dis-lp[i])/pt.norm[i]:(lp[i]+dis)/(-pt.norm[i]));
					back[i]=fabs((pt.norm[i]>0)? (lp[i]+dis)/(pt.norm[i]):(1+dis-lp[i])/(-pt.norm[i]));
				}
				if(forward[0]<forward[1])
				{
					if(forward[0]<forward[2])
					{
						s1=forward[0]*pow(0.5,d);
					}else
					{
						s1=forward[2]*pow(0.5,d);
					}
				}else
				{
					if(forward[1]<forward[2])
					{
						s1=forward[1]*pow(0.5,d);
					}else
					{
						s1=forward[2]*pow(0.5,d);
					}
				}
				if(back[0]<back[1])
				{
					if(back[0]<back[2])
					{
						s2=back[0]*pow(0.5,d);
					}else
					{
						s2=back[2]*pow(0.5,d);
					}
				}else
				{
					if(back[1]<back[2])
					{
						s2=back[1]*pow(0.5,d);
					}else
					{
						s2=back[2]*pow(0.5,d);
					}
				}


			}else
			{
				for(int i=0;i<3;i++)
				{
					float upbound=0.0;
					float lowbound=0.0;
					if(ind[i]==1)
					{
						upbound=1+dis;
						lowbound=0.5-dis;
					}
					else
					{
						upbound=0.5+dis;
						lowbound=-dis;
					}
					forward[i]=fabs((pt.norm[i]>0)? (upbound-lp[i])/pt.norm[i]:(lp[i]-lowbound)/(-pt.norm[i]));
					back[i]=fabs((pt.norm[i]>0)? (lp[i]-lowbound)/(pt.norm[i]):(upbound-lp[i])/(-pt.norm[i]));
				}
				if(forward[0]<forward[1])
				{
					if(forward[0]<forward[2])
					{
						s1=forward[0]*pow(0.5,d);
					}else
					{
						s1=forward[2]*pow(0.5,d);
					}
				}else
				{
					if(forward[1]<forward[2])
					{
						s1=forward[1]*pow(0.5,d);
					}else
					{
						s1=forward[2]*pow(0.5,d);
					}
				}
				if(back[0]<back[1])
				{
					if(back[0]<back[2])
					{
						s2=back[0]*pow(0.5,d);
					}else
					{
						s2=back[2]*pow(0.5,d);
					}
				}else
				{
					if(back[1]<back[2])
					{
						s2=back[1]*pow(0.5,d);
					}else
					{
						s2=back[2]*pow(0.5,d);
					}
				}


			}		

		}
		else
		{
			vect3f off(ext2[0]*ind[0], ext2[1]*ind[1], ext2[2]*ind[2]);
			d=children(ind[0],ind[1],ind[2]).getstep(pt, maxDepth, depth+1, mine+off, ext2, minDepth,s1,s2);
		}
		
	}
	return(s1);
	
}

