// Author: Josiah Manson

// #define OPENGL_DISPLAY // when undefined, the glut library will not need to be linked to

#ifdef OPENGL_DISPLAY
//  #include <GL/glut.h>
#endif
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <fstream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "vect.h"
#include "global.h"

using namespace std;

#ifdef OPENGL_DISPLAY
void reshape(int width, int height)
{
	glViewport(0, 0, width, height);
}

void drawTris()
{
	glBegin(GL_TRIANGLES);
	glColor3f(0, 0, .7);
	for (int i = 0; i < g.reconstructed.triangles.size(); i++)
	{
		vect3f &v0 = g.reconstructed.triangles[i].verts[0];
		vect3f &v1 = g.reconstructed.triangles[i].verts[1];
		vect3f &v2 = g.reconstructed.triangles[i].verts[2];

		vect3f norm = (v0 - v1) % (v2 - v1);
		norm.unit();
		glNormal3fv(norm.v);

		glVertex3fv(v2.v);
		glVertex3fv(v1.v);
		glVertex3fv(v0.v);
	}
	glEnd();
}

void display() 
{
	float color_pts[4] = {1,.3,.3, 1};
	float color_surf[4] = {.8,.8,.8, 1};
	float color_lines[4] = {.2,.2,.8, 1};

	// Set color and clear the buffers
	glClearColor(g.bg_color[0], g.bg_color[1], g.bg_color[2], 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHT0);

	// Set up matrices
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(35, float(g.width)/g.height, .0001, 10);
	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glTranslatef(0, 0, -g.zoom);
	glMultMatrixf(g.rotation);
	glTranslatef(-g.focus[0], -g.focus[1], -g.focus[2]);

	// draw reconstructed isosurface
	if (g.showSurf)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_CULL_FACE);
		glDisable(GL_TEXTURE_2D);

		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_surf);

		drawTris();

//		glDisable(GL_LIGHTING);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_lines);

		// draw mesh lines of isosurface (offset forward)
		if (g.showLines)
		{
			glEnable(GL_POLYGON_OFFSET_LINE);
			glPolygonOffset(-1, -1);
			glPolygonMode(GL_FRONT, GL_LINE);

			drawTris();
			
			glPolygonMode(GL_FRONT, GL_FILL);
			glDisable(GL_POLYGON_OFFSET_LINE);
		}
	}
	else if (g.showLines)
	{
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_surf);
		glPolygonMode(GL_FRONT, GL_LINE);

		drawTris();
		
		glPolygonMode(GL_FRONT, GL_FILL);
	}

	// draw points
	if (g.showPoints && g.do_load_points)
	{
		const int ngon = 10;
		const float ptsize = .003;
		glEnable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_pts);

		for (int i = 0; i < g.points.size(); i++)
		{
			vect3f &norm = g.points[i].norm;

			vect3f dir = vect3f(1,1,0);
			vect3f x = ~(dir - norm*(dir*norm));
			vect3f y = norm % x;

			x *= ptsize;
			y *= ptsize;

			if (sizeof(float) == sizeof(float))
				glNormal3fv((GLfloat*)norm.v);
			else
				glNormal3dv((GLdouble*)norm.v);

			glBegin(GL_POLYGON);
			for (int n = 0; n < ngon; n++)
			{
				vect3f pos = g.points[i].pos + x*cos(n*2*3.14159/ngon) + y*sin(n*2*3.14159/ngon);
				
				if (sizeof(float) == sizeof(float))
					glVertex3fv((GLfloat*)pos.v);
				else
					glVertex3dv((GLdouble*)pos.v);
			}
			glEnd();
		}
		
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color_surf);

		glColor3f(1,1,1);
		glDisable(GL_LIGHTING);
	}

	
	// draw plane slice
	if (g.showPlane && g.do_load_plane)
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

		// draw textured quad
		glDisable(GL_LIGHTING);
		glEnable(GL_TEXTURE_2D);
		glDisable(GL_CULL_FACE);

		glBindTexture(GL_TEXTURE_2D, g.plane_tex);

		glColor3f(1,1,1);
		glBegin(GL_QUADS);

		glTexCoord2f(0,0);
		glVertex3fv((corner).v);

		glTexCoord2f(1,0);
		glVertex3fv((corner + X).v);

		glTexCoord2f(1,1);
		glVertex3fv((corner + X + Y).v);

		glTexCoord2f(0,1);
		glVertex3fv((corner + Y).v);

		glEnd();
	}

	// Signal that the scene is complete
	glFlush();
	glutSwapBuffers();
}

void keyboard(unsigned char key, int x, int y)
{
	if (key == 'l')
		g.showLines = !g.showLines;
	else if (key == 's')
		g.showSurf = !g.showSurf;
	else if (key == 'p')
		g.showPoints = !g.showPoints;
	else if (key == 'o')
		g.showPlane = !g.showPlane;
	else if (key == '[')
	{
		g.plane_offset--;
		colorPlane();
	}
	else if (key == ']')
	{
		g.plane_offset++;
		colorPlane();
	}
	else if (key == '-')
	{
		g.plane_offset -= 10;
		colorPlane();
	}
	else if (key == '=')
	{
		g.plane_offset += 10;
		colorPlane();
	}
	else if (key == 'x')
	{
		g.plane_dir = 0;
		colorPlane();
	}
	else if (key == 'y')
	{
		g.plane_dir = 1;
		colorPlane();
	}
	else if (key == 'z')
	{
		g.plane_dir = 2;
		colorPlane();
	}
	else if (key == 'r')
		g.readView();
	else if (key == 'w')
		g.writeView();
	else if (key >= '0' && key <= '9')
		g.viewFile = key - '0';
	
	display();
}

void mouse(int button, int state, int x, int y)
{
	y = g.height - y;

	// Mouse state that should always be stored on pressing
	if (state == GLUT_DOWN)
	{
		g.mouse_home_x = x;
		g.mouse_home_y = y;
		g.mouse_prev_x = x;
		g.mouse_prev_y = y;
	}

	// Keep track of the mouse state (GLUT seems not to have a polling function. Boo!)
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		g.mouse_left = true;
	}
	if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		g.mouse_left = false;
	}

	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		g.mouse_right = true;
	}
	if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
	{
		g.mouse_right = false;
	}
	
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
	{
		g.mouse_middle = true;
	}
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_UP)
	{
		g.mouse_middle = false;
	}
}

void motion(int x, int y)
{
	// flip the mouse position to be a bit nicer
	y = g.height - y;

	// rotate the scene
	if (g.mouse_left)
	{
		float dx = (x - g.mouse_prev_x);
		float dy = (y - g.mouse_prev_y);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(dx, 0, 1, 0);
		glRotatef(dy, -1, 0, 0);
		glMultMatrixf(g.rotation);
		glGetFloatv(GL_MODELVIEW_MATRIX, g.rotation);
	}
	else if (g.mouse_middle)
	{
		float dx = (x - g.mouse_prev_x);
		float dy = (y - g.mouse_prev_y);

		vect3f X(g.rotation[0], g.rotation[4], g.rotation[8]);
		vect3f Y(g.rotation[1], g.rotation[5], g.rotation[9]);

		g.focus -= X*(2e-3*dx) + Y*(2e-3*dy);
	}
	else if (g.mouse_right)
	{
		float dy = (y - g.mouse_prev_y);
		g.zoom -= 1e-2*dy;
	}
	
	// Store previous mouse positions
	g.mouse_prev_x = x;
	g.mouse_prev_y = y;
	
	display();
}
#endif

bool atob(string s)
{
	if (s == "t" || s == "T" || s == "true" || s == "1")
		return true;
	return false;
}

string get_tok(ifstream &f)
{
	string s;
	char c;
	char buf[1025];

	while (f.get(c) && !f.eof())
	{
		if (c == ' ' || c == '\n' || c == '\r' || c == '\t')
			break;
		else if (c == '#')
		{
			f.getline(buf, 1024, '\n');
			break;
		}
		else
			s += c;
	}

	return s;
}

vector<string> parse_file(string fn)
{
	vector<string> v;
	ifstream f(fn.c_str());

	while (!f.eof())
	{
		string s = get_tok(f);
		if (!s.empty())
			v.push_back(s);
	}

	f.close();
	return v;
}

void init(vector<string> args)
{
	if (args.empty())
	{
		printf("Parameters can be passed to the command line by first\n");
		printf("giving the name of the parameter to modify followed\n");
		printf("by its value.\n");
		printf("\n");
		printf("pts1         [.pts]\n");
		printf("pts2         [.pts]\n");
		printf("depth       [int]\n");
		printf("N           [int]\n");
		printf("iterstep    [int]\n");
		printf("stepA       [float]\n");
		printf("stepB       [float]\n");
		printf("stepC       [float]\n");
		printf("wavelet     (haar|daub4)\n");
		printf("to_screen   (t|f)\n");
		printf("to_file     (t|f)\n");
		printf("prune       (t|f)\n");
		printf("blur        (t|f)\n");
		printf("pts_node    (t|f)\n");
		printf("above_below (t|f)\n");
		printf("eval_pass   (t|f)\n");
		printf("stream      (t|f)\n");
		printf("surf_at_pts (t|f)\n");
		printf("fill_in_mem (t|f)\n");
		printf("stay_open   (t|f)\n");
		printf("load_plane  (t|f)\n");
		printf("load_points (t|f)\n");
		printf("cfg         [config_file]\n");

		exit(0);
	}

	for (int i = 0; i < args.size(); i++)
	{
		if (args[i] == "pts")
			g.pointsfile = args[++i];
		else if (args[i] == "depth")
			g.depth = atoi(args[++i].c_str());
		else if (args[i] == "iternum")
			g.iternum = atoi(args[++i].c_str());
		else if (args[i] == "N")
			g.N = atoi(args[++i].c_str());
		else if (args[i] == "stepA")
			g.stepA = atof(args[++i].c_str());
		else if (args[i] == "stepB")
			g.stepB = atof(args[++i].c_str());
		else if (args[i] == "stepC")
			g.stepC = atof(args[++i].c_str());
		else if (args[i] == "wavelet")
		{
			i++;
			if (args[i] == "haar")
				g.wavelet = haar;
			else if (args[i] == "daub4")
				g.wavelet = daub4;
		}
		else if (args[i] == "cfg")
			init(parse_file(args[++i]));
		else if (args[i] == "to_screen")
			g.output_screen = atob(args[++i]);
		else if (args[i] == "to_file")
			g.output_file = atob(args[++i]);
		else if (args[i] == "prune")
			g.do_prune = atob(args[++i]);
		else if (args[i] == "blur")
			g.do_blur = atob(args[++i]);
		else if (args[i] == "pts_node")
			g.do_pts_node = atob(args[++i]);
		else if (args[i] == "above_below")
			g.do_above_below = atob(args[++i]);
		else if (args[i] == "eval_pass")
			g.do_eval_pass = atob(args[++i]);
		else if (args[i] == "stream")
			g.do_streaming = atob(args[++i]);
		else if (args[i] == "surf_at_pts")
			g.do_surf_only_at_points = atob(args[++i]);
		else if (args[i] == "fill_in_mem")
			g.do_fill_in_mem = atob(args[++i]);
		else if (args[i] == "stay_open")
			g.stay_open = atob(args[++i]);
		else if (args[i] == "load_plane")
			g.do_load_plane = atob(args[++i]);
		else if (args[i] == "load_points")
			g.do_load_points = atob(args[++i]);
	}
}

vector<string> arglist(int argc, char **argv)
{
	vector<string> args;

	for (int i = 1; i < argc; i++)
		args.push_back(argv[i]);

	return args;
}




int main(int argc, char **argv)
{
	init(arglist(argc, argv));
	g.pointsfile2 = g.pointsfile;
	g.do_prune=true;
	g.do_above_below = false;
	g.do_eval_pass = true;
	g.do_blur=true;
	g.output_file=true;
	g.wavelet=daub4;
	g.do_streaming=false;
	g.do_load_plane = false;
	g.depth=9;
	g.do_pts_node=false;
	g.smooth=false;
	reconstruct();

#ifdef OPENGL_DISPLAY
	if (g.output_screen)
	{
		// Initialize GLUT
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
		glutInitWindowPosition(0, 0);
		glutInitWindowSize(g.width, g.height);
		glutCreateWindow("Streaming Wavelet Surface Reconstruction");

		// read points
		if (g.do_load_points)
			readPointsForDisplay(g.pointsfile);
		
		// create plane slice through function
		int plane_size = g.res >> (g.wavelet == haar ? 0 : 2);
		g.plane_img.resize(plane_size, plane_size);
		g.plane_offset = plane_size / 2;
		g.plane_dir = 0;
		if (g.do_load_plane)
		{
			glGenTextures(1, &g.plane_tex);
			colorPlane();
		}

		// display
		display();

		// Register callback functions
		glutDisplayFunc(display);
		glutReshapeFunc(reshape);
		glutKeyboardFunc(keyboard);
		glutMouseFunc(mouse);
		glutMotionFunc(motion);

		// Enter event processing loop
		glutMainLoop();
	}
#endif

	if (g.stay_open)
		system("pause");
}
