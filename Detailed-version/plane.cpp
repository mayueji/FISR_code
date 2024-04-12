#include "global.h"
// #include <GL/glut.h>

void colorPlane()
{
	// find x,y directions and the corner of the texture
	vect3f corner(0,0,0);
	corner[g.plane_dir] = float(g.plane_offset) / float(g.plane_img.width);
	vect3f X,Y;

	if (g.plane_dir == 0)
	{
		X(0,1,0);
		Y(0,0,1);
	}
	else if (g.plane_dir == 1)
	{
		X(1,0,0);
		Y(0,0,1);
	}
	else if (g.plane_dir == 2)
	{
		X(1,0,0);
		Y(0,1,0);
	}

	X /= g.plane_img.width;
	Y /= g.plane_img.height;

	corner += X*.5 + Y*.5;

	// fill the image
	for (int i = 0; i < g.plane_img.width; i++)
	{
		for (int j = 0; j < g.plane_img.height; j++)
		{
			vect3f pos = corner + X*i + Y*j;
			float val = g.tree.root.getValueAtPoint(pos, g.tree.mine, g.tree.maxe - g.tree.mine);

			if (val >= average_val)
				g.plane_img(i,j)(val/(average_val) - 1, 0, 0);
			else
				g.plane_img(i,j)(0, 0, 1 - val/(average_val));
		}
	}

}
