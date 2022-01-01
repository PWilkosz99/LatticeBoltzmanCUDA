#ifndef Visualization_h
#define Visualization_h

#include "glut.h"
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include <iostream>
#include <fstream>
#include "LBM0.h"

void Initialize();

void DrawLBM1();

void DrawLBM2();

void DrawAtmosHor(char[Npic * Ny][Npic * Nx][3]);

void DrawAtmosVer(char[Npic * Ny][Npic * Nx][3]);

void Timer(int value);

void BMPout(char[Npic * Ny][Npic * Nx][3], char[10]);

void DrawLine(int x0, int y0, int x1, int y1, char rgbp[Npic * Ny][Npic * Nx][3]);

void DrawLine(int x0, int y0, int x1, int y1, char rgbp[Npic * Ny][Npic * Nx][3], int[3]);


#endif // !Visualization_h