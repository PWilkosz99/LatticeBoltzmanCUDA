#include "glut.h"
#include <math.h>
#include <stdlib.h>
#include "Visualization.h"
#include "LBM0.h"

void KeyboardEsc(unsigned char key, int x, int y);

int WinId, win1, win2;
bool FirstCycle = true;
bool LastCycle = false;
int icycle = 0;

static void TimeEvent(int te) /////////////////////////////////////////////////

{
	glutSetWindow(win1);
	glutPostRedisplay();
	glutSetWindow(win2);
	glutPostRedisplay();
	glutTimerFunc(10, TimeEvent, 1); // Reset our timer.
}

int main(int argc, char** argv) ///////////////////////////////////////////////
{
	// Initialization

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(768, 768);
	glutInitWindowPosition(0, 20);
	glutTimerFunc(10, TimeEvent, 1);
	win1 = glutCreateWindow("LBM-CUDA [Vy]");

	// Registration

	glutDisplayFunc(DrawLBM1);
	Initialize();
	glutInitWindowPosition(768, 20);
	win2 = glutCreateWindow("LBM-CUDA [Vx]");
	glutDisplayFunc(DrawLBM2);
	Initialize();
	glutKeyboardFunc(KeyboardEsc);
	glutMainLoop();
	return 0;
}

void KeyboardEsc(unsigned char key, int x, int y) //////////////////////
{
	switch (key)
	{
	case 27: //Escape key
		glutDestroyWindow(WinId);
		exit(0);
		break;
	}
	glutPostRedisplay();
}

void DrawLBM1() /////////////////////////////////////////////////////////
{
	auto rgb = new char[Npic * Ny][Npic * Nx][3];
	if (mainLBM(FirstCycle) == 1) LastCycle = true;
	DrawAtmosVer(rgb);
	const auto bmpfilelen = 11;
	char bmpfile[bmpfilelen];
	snprintf(bmpfile, bmpfilelen, "Uy%04d.bmp", icycle);
	BMPout(rgb, bmpfile);
	delete[] rgb;
	FirstCycle = false;
}

void DrawLBM2() ///////////////////////////////////////////////////////
{
	auto rgb = new char[Npic * Ny][Npic * Nx][3];
	DrawAtmosHor(rgb);
	const auto bmpfilelen = 11;
	char bmpfile[bmpfilelen];
	snprintf(bmpfile, bmpfilelen, "Ux%04d.bmp", icycle);
	BMPout(rgb, bmpfile);
	delete[] rgb;
	FirstCycle = false;
	icycle++;
	if (LastCycle) { exit(0); }
}