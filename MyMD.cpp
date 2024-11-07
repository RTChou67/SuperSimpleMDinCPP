#include <iostream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <random>
#include <iomanip>
#include <cstdlib>
#include <string>
// Writen by Chengze Zhou, 2024.
// My first MD programme, writen in C++, Last update:2024.10.31
/*
Some NOTICE:
Compiling option -O3 and -Ofast is enabled.
The current single core performance ~ 0.15s/100steps for a 1000-atom-system, so neighboring list is not implemented.
/*
After loading the .gro file to VMD, run the following 2 consoles:
mol modstyle 0 top VDW
pbc box -on
The first one is to make the atoms look more clearly, and the second one toggles the visualization of the box.
*/
using namespace std;
struct Mol
{
	// This structure contains all the molecular information needed.
	double x[1000], y[1000], z[1000];	 // position
	double vx[1000], vy[1000], vz[1000]; // velocity
	double fx[1000], fy[1000], fz[1000]; // force
	string type[1000];					 // type
};
static double Input(string VarName, string Type, bool Judge)
{
	// This function determines the legitimacy of the input.
	// All three variables should be positive, but only StepNum should be an integer, so isInteger is introduced.
	// For StepNum, Input should be positive, and isInteger should be 0(False)
	// For StepSize and Cutoff, Judge is false, so whatever isInteger is, the input is legal.
	double Input;
	cout << "Input the " << VarName << endl;
	cin >> Input;
	bool isInteger = (Input - floor(Input));
	if (Input <= 0 or (Judge && isInteger) == true)
	{
		cout << "ERROR: The " << VarName << " should be a " << Type << endl;
		exit(1);
	}
	return Input;
}
static bool ifTraj()
{
	// This function determines whether to output trajectory file.
	string x;
	cout << "Record trajectory?[YES/NO]" << endl;
	cin >> x;
	if (x == "YES")
	{
		return true;
	}
	else if (x == "NO")
	{
		return false;
	}
	else
	{
		cout << "Please enter \"YES\" or \"NO\"";
		exit(1);
	}
}
static void InitialCoord(Mol &mol)
{
	// This function generates initial coordinates.
	for (int i = 0; i < 500; i++)
	{
		mol.x[i] = i / 100 + 0.5;
		mol.y[i] = (i / 10) % 10 + 0.5;
		mol.z[i] = (i % 10) + 0.5;
		mol.type[i] = "HE";
	}
	for (int i = 500; i < 1000; i++)
	{
		mol.x[i] = i / 100 + 0.5;
		mol.y[i] = (i / 10) % 10 + 0.5;
		mol.z[i] = (i % 10) + 0.5;
		mol.type[i] = "NE";
	}
}
static void InitialVelo(Mol &mol)
{
	// This function generates initial velocities expressed in spherical coordinates.
	double PI = 3.14159265358979323846264;
	double PI2 = PI * 2;
	double Vtot = sqrt(3);
	for (int i = 0; i < 1000; i++)
	{
		int seed = rand();
		mt19937 engine(seed); // Mersenne Twister random number generator
		uniform_real_distribution<double> a(0.0, PI);
		uniform_real_distribution<double> b(0.0, PI2);
		double theta = a(engine);
		double phi = b(engine);
		mol.vx[i] = sin(theta) * cos(phi) * Vtot;
		mol.vy[i] = sin(theta) * sin(phi) * Vtot;
		mol.vz[i] = cos(theta) * Vtot;
	}
}
static void VeloVerlet(Mol &mol, double HalfStepSize, double StepSize, int Status)
{
	// This function executes the first and the last step of the Velocity-Verlet process.
	for (int n = 0; n < 1000; n++)
	{
		mol.vx[n] += mol.fx[n] * HalfStepSize;
		mol.vy[n] += mol.fy[n] * HalfStepSize;
		mol.vz[n] += mol.fz[n] * HalfStepSize;
		if (Status == 1)
		{
			mol.x[n] += mol.vx[n] * StepSize;
			mol.y[n] += mol.vy[n] * StepSize;
			mol.z[n] += mol.vz[n] * StepSize;
		}
	}
}
static void Force(Mol &mol, double Cutoff2)
{
	// This function calculates the force according to the current atomic coordinates.
	for (int n = 0; n < 1000; n++)
	{ // initializing force
		mol.fx[n] = mol.fy[n] = mol.fz[n] = 0.0;
	}
	for (int i = 0; i < 1000; i++)
	{
		for (int j = i + 1; j < 1000; j++)
		{ // calculating force according to vdW potential
			double xij = mol.x[j] - mol.x[i];
			double yij = mol.y[j] - mol.y[i];
			double zij = mol.z[j] - mol.z[i];
			double r2 = xij * xij + yij * yij + zij * zij;
			if (r2 < Cutoff2)
			{ // determines whether this pair should be skipped according to the cutoff distance
				continue;
			}
			double r2inv = 1 / r2;
			double r8inv = pow(r2inv, 4);
			double r14inv = pow(r2inv, 7);
			double f_ij = 6 * r8inv - 12 * r14inv;
			mol.fx[i] += f_ij * xij;
			mol.fx[j] -= f_ij * xij;
			mol.fy[i] += f_ij * yij;
			mol.fy[j] -= f_ij * yij;
			mol.fz[i] += f_ij * zij;
			mol.fz[j] -= f_ij * zij;
		}
	}
}
static double Periodic1(double x, double length)
{
	while (x < 0)
	{
		x += length;
	}
	while (x > length)
	{
		x -= length;
	}
	return x;
}

static void Periodic(Mol &mol, double length)
{
	for (int i = 0; i < 1000; i++)
	{
		mol.x[i] = fmod(fmod(mol.x[i], length) + length, length);
		mol.y[i] = fmod(fmod(mol.y[i], length) + length, length);
		mol.z[i] = fmod(fmod(mol.z[i], length) + length, length);
	}
}
static void ScaleVelo(Mol &mol)
{
	// This function scales the velocities of atoms acoording to the NVT ensemble.
	double Ek = 0;
	for (int i = 0; i < 1000; i++)
	{ // calculating kinetic energy
		Ek += mol.vx[i] * mol.vx[i] + mol.vy[i] * mol.vy[i] + mol.vz[i] * mol.vz[i];
	}
	double scaleV = sqrt(3000 / Ek); // calculating scale
	for (int i = 0; i < 1000; i++)
	{ // scaling velocity;
		mol.vx[i] *= scaleV;
		mol.vy[i] *= scaleV;
		mol.vz[i] *= scaleV;
	}
}
static void WriteCoord(Mol &mol, string filename, int empty, int status)
{
	// This function is responsible for outputting the atomic coordinates.
	// Different code for different conditions, as listed below:
	// 1. initial or result coordinates, then directly covers the previous file if it exists.
	// 2. trajectory file:
	// 2.1 if it is the first step, then clear the previous file
	// 2.2 if it is not the first step, then append to the previous file, and write the current Step number as the title
	fstream coord;
	coord.precision(3);
	string filename1 = filename + ".gro";
	if (empty == 1)
	{
		coord.open(filename1.c_str(), ios::out);
		coord.close();
	}
	coord.open(filename1.c_str(), ios::app);
	if (status == -100)
	{
		coord << filename << " Coordinates" << endl
			  << "  1000" << endl;
	}
	else
	{
		coord << "Step: " << status << endl
			  << "  1000" << endl;
	}
	for (int i = 0; i < 1000; i++)
	{
		coord << "    1MOL     " << mol.type[i] << right << setw(5) << fixed << i + 1 << "   " << mol.x[i] << "   " << mol.y[i] << "   " << mol.z[i] << endl;
	}
	coord << "   10.000   10.000   10.000" << endl;
	coord.close();
}

int main()
{
	Mol mol;
	int StepNum = (int)Input("number of steps", "POSITIVE INTEGER!!!", true);  // inputing the number of steps
	double StepSize = Input("size of each step", "POSITIVE NUMBER!!!", false); // inputing the size of each step
	double Cutoff = Input("cutoff distance", "POSITIVE NUMBER!!!", false);	   // inputing cutoff distance
	bool iftraj = ifTraj();
	double length = 10.0; // defining box size
	InitialCoord(mol);	  // generating initial coordinates
	WriteCoord(mol, "initial", 1, -100);
	double HalfStepSize = 1 / 2 * StepSize;
	double Cutoff2 = Cutoff * Cutoff;
	double tElapsed1 = 0;
	double tElapsed2 = 0;
	cout << "Start simulation:" << endl;
	auto TimeStart = chrono::steady_clock::now();
	auto TimeNow = chrono::steady_clock::now();
	InitialVelo(mol); // generating initial random
	for (int step = 1; step <= StepNum; step++)
	{												// start simulation
		VeloVerlet(mol, HalfStepSize, StepSize, 1); // velocity verlet step 1
		Force(mol, Cutoff2);						// calculating force
		VeloVerlet(mol, HalfStepSize, StepSize, 2); // velocity verlet step 2
		Periodic(mol, length);						// periodic boundary conditions
		ScaleVelo(mol);								// velocity scaling
		if (iftraj == true && step % 50 == 0)
		{ // output trajectory file if wanted
			if (step == 50)
			{
				WriteCoord(mol, "traj", 1, step); // if it is the first step, clear traj file first
			}
			else
			{
				WriteCoord(mol, "traj", 0, step); // if it is not the first step,append to the previous file
			}
		}
		if (step % 100 == 0)
		{ // timing module, output time every 100 frames
			/*
			Here I used a little trick to produce the time used for the current 100 steps.
			By using different refreshing time, the time duration of 2 adjacent 100 steps are successfully produced.
			This trick was neccesary because an array cannot be declared with the "auto" keyword.
			*/
			tElapsed1 = chrono::duration<double>(TimeNow - TimeStart).count();
			TimeNow = chrono::steady_clock::now();
			tElapsed2 = chrono::duration<double>(TimeNow - TimeStart).count();
			double tElapsed100 = tElapsed2 - tElapsed1;
			double tExpect = tElapsed2 * ((double)StepNum / (double)step - 1);
			int num = (int)floor(log10(StepNum) + 2);
			cout << left << setw(num) << step << "steps completed using: " << left << setw(10) << tElapsed2 << "s. Time per 100 steps: " << left << setw(10) << tElapsed100 << "s. Estimated remaining time: " << tExpect << "s." << endl;
		}
	}
	WriteCoord(mol, "output", 1, -100); // generating output coordinates
	TimeNow = chrono::steady_clock::now();
	tElapsed2 = chrono::duration<double>(TimeNow - TimeStart).count();
	cout << "Simulation successfully ended using " << tElapsed2 << "s." << endl;
	system("pause");
	return 0;
}