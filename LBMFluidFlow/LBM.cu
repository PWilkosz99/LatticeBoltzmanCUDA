#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "glut.h"
#include <string.h>

#include "LBM0.h"
#include "Visualization.h"
#include "LBM_Atmos.h"

std::fstream infile; // file with input data
float timev; // real time from cycle beginning
int npasses, ipass; // number of cycles and cycle number
float timecycle[40]; //duration of cycle

int mainLBM(bool FirstCycle)

{
	using namespace std;
	dim3 block(8, 8, 1); // it could be different
	dim3 grid(Nx / block.x, Ny / block.y, 1);
	dim3 gridX(1, Nx / block.y, 1);
	dim3 gridY(1, Ny / block.y, 1);
	if (FirstCycle)
	{
		FirstCycle = false;
		if (!NewSim) InRes();
		InitialAtmos << < grid, block >> > ();
		cudaDeviceSynchronize();
		timem = 0.0f; timev = 0.0f;

		infile.open("Schedule.txt", ios::in);
		infile >> npasses; ipass = 0; // read number of cycles
		for (int i = 0; i < npasses; i++) { infile >> timecycle[i]; }
		infile.close();
		return 0;
	}
	else
	{
		while ((timev < timecycle[ipass]) && (ipass < npasses))
		{
			stept = 0.01; // time step
			EquiRelaxAtmos << < grid, block >> > ();
			cudaDeviceSynchronize();

			StreamingAtmos << < grid, block >> > ();
			cudaDeviceSynchronize();

			BoundaryEast << < gridX, block >> > ();
			BoundaryWest << < gridX, block >> > ();
			BoundarySouth << < gridY, block >> > ();
			BoundaryNord << < gridY, block >> > ();
			cudaDeviceSynchronize();
			timem += stept; timev += stept;
		}
		{ timev = 0.0f; ipass++; }
	}
	if (ipass >= npasses) { OutRes(); return 1; }; // exit if it was the last cycle
	return 0;
}