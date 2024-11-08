#include <chrono>    //timing module
#include <cmath>     //mathematical functions
#include <cstdlib>   //prevent from automatic shutdown
#include <fstream>   //output files
#include <iomanip>   //output format
#include <iostream>  //input and output
#include <random>    //mt19937 engine
#include <string>    //string
// Writen by Chengze Zhou, 2024.
// My first MD program, writen in C++, Last update:2024.10.31
/*
Some NOTICE:
This program supports customized number of steps, size of each step and cutoff distance.
Compiling option -O3 and -Ofast is enabled.
The current single core performance ~ 0.15s/100steps for a 1000-atom-system, so neighboring list is not implemented.
After loading the .gro file to VMD, run the following 2 commands:
    mol modstyle 0 top VDW
    pbc box -on
The first one is to adjust the visualization of  the atoms, and the second one toggles the visualization of the box.
*/
using namespace std;

// Structure is used to record all the molecular information needed.
struct Mol {
    double x[1000], y[1000], z[1000];        // position
    double vx[1000], vy[1000], vz[1000];     // velocity
    double fx[1000], fy[1000], fz[1000];     // forcew
    string type[1000];   // type
};

// Input function judges the legibility of and record the parameters needed in the simulation. StepNum should be a positive integer, while StepNum and Cutoff should be a positive number.
// In this case, isInteger is introduced to judge whether the input is an integer when variable Judge is true.
static double Input(string VarName, string Type, bool Judge) {
    double Input;
    cout << "Input the " << VarName << endl;
    cin >> Input;
    bool isInteger = (Input - floor(Input));
    if (Input <= 0 or (Judge && isInteger) == true) {
        cout << "ERROR: The " << VarName << " should be a " << Type << endl;
        exit(1);
    }
    return Input;
}

// This function determines whether to output trajectory file.
static bool ifTraj() {
    string x;
    cout << "Record trajectory?[YES/NO]" << endl;
    cin >> x;
    if (x == "YES") {
        return true;
    } else if (x == "NO") {
        return false;
    } else {
        cout << "Please enter \"YES\" or \"NO\"";
        exit(1);
    }
}

// This function generates initial coordinates.
static void InitialCoord(Mol& mol) {
    for (int i = 0; i < 500; i++) {
        mol.x[i] = i / 100 + 0.5;
        mol.y[i] = (i / 10) % 10 + 0.5;
        mol.z[i] = (i % 10) + 0.5;
        mol.type[i] = "HE";
    }
    for (int i = 500; i < 1000; i++) {
        mol.x[i] = i / 100 + 0.5;
        mol.y[i] = (i / 10) % 10 + 0.5;
        mol.z[i] = (i % 10) + 0.5;
        mol.type[i] = "NE";
    }
}

// This function generates the initial velocities of atoms under spherical coordinates using a MT19937 engine within the <random> library.
static void InitialVelo(Mol& mol) {
    double PI = 3.14159265358979323846264;
    double PI2 = PI * 2;
    double Vtot = sqrt(3);
    for (int i = 0; i < 1000; i++) {
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

// This function executes the first and the last step of the Velocity-Verlet process.
// If Status == 1, then it's the first step; If Status == 2, then it's the last step.
static void VeloVerlet(Mol& mol, double HalfStepSize, double StepSize, int Status) {
    for (int n = 0; n < 1000; n++) {
        mol.vx[n] += mol.fx[n] * HalfStepSize;
        mol.vy[n] += mol.fy[n] * HalfStepSize;
        mol.vz[n] += mol.fz[n] * HalfStepSize;
        if (Status == 1) {
            mol.x[n] += mol.vx[n] * StepSize;
            mol.y[n] += mol.vy[n] * StepSize;
            mol.z[n] += mol.vz[n] * StepSize;
        }
    }
}

// This function calculates the force according to the current atomic coordinates.
static void Force(Mol& mol, double Cutoff2) {
    for (int n = 0; n < 1000; n++) { // initializing force
        mol.fx[n] = mol.fy[n] = mol.fz[n] = 0.0;
    }
    for (int i = 0; i < 1000; i++) {
        for (int j = i + 1; j < 1000; j++) { // calculating force according to vdW potential
            double xij = mol.x[j] - mol.x[i];
            double yij = mol.y[j] - mol.y[i];
            double zij = mol.z[j] - mol.z[i];
            double r2 = xij * xij + yij * yij + zij * zij;
            if (r2 < Cutoff2) { // determines whether this pair should be skipped according to the cutoff distance
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

// A double-layer Modulo function is used because of the definition. For example, fmod(-7.0,3.0) = -1 instead of the 2 we wanted.
static void Periodic(Mol& mol, double length) {
    for (int i = 0; i < 1000; i++) {
        mol.x[i] = fmod(fmod(mol.x[i], length) + length, length);
        mol.y[i] = fmod(fmod(mol.y[i], length) + length, length);
        mol.z[i] = fmod(fmod(mol.z[i], length) + length, length);
    }
}

// This function scales the velocities of atoms acoording to the NVT ensemble.
static void ScaleVelo(Mol& mol) {
    double Ek = 0;
    for (int i = 0; i < 1000; i++) { // calculating kinetic energy
        Ek += mol.vx[i] * mol.vx[i] + mol.vy[i] * mol.vy[i] + mol.vz[i] * mol.vz[i];
    }
    double scaleV = sqrt(3000 / Ek); // calculating scale
    for (int i = 0; i < 1000; i++) { // scaling velocity;
        mol.vx[i] *= scaleV;
        mol.vy[i] *= scaleV;
        mol.vz[i] *= scaleV;
    }
}

// If it is the initial or result coordinates, or the first step of the trajectory file, then empty the file first.
// Then, append to the original file.
// The variable filename and status determines the title and file name of the .gro file.
static void WriteCoord(Mol& mol, string filename, int empty, int status) {
    fstream coord;
    coord.precision(3);
    string filename1 = filename + ".gro";
    if (empty == 1) {
        coord.open(filename1.c_str(), ios::out);
        coord.close();
    }
    coord.open(filename1.c_str(), ios::app);
    if (status == -100) {
        coord << filename << " Coordinates" << endl
              << "  1000" << endl;
    } else {
        coord << "Step: " << status << endl
              << "  1000" << endl;
    }
    for (int i = 0; i < 1000; i++) {
        coord << "    1MOL     " << mol.type[i] << right << setw(5) << fixed << i + 1 << "   " << mol.x[i] << "   " << mol.y[i] << "   " << mol.z[i] << endl;
    }
    coord << "   10.000   10.000   10.000" << endl;
    coord.close();
}

// Since the auto type cannot be used as a function's reference, the timing module is written inside the main function timing every 100 steps.
// Since the auto type can not be stored in an array, two asynchronous timing is used, and the difference between these two are the time used in this 100 steps.
// Based on the current time the expected elapsing time is calculated and outputted.
int main() {
    Mol mol;
    int StepNum = (int)Input("number of steps", "POSITIVE INTEGER!!!", true); // inputing the number of steps
    double StepSize = Input("size of each step", "POSITIVE NUMBER!!!", false); // inputing the size of each step
    double Cutoff = Input("cutoff distance", "POSITIVE NUMBER!!!", false); // inputing cutoff distance
    bool iftraj = ifTraj();
    double length = 10.0; // defining box size
    InitialCoord(mol); // generating initial coordinates
    WriteCoord(mol, "initial", 1, -100);
    double HalfStepSize = 1 / 2 * StepSize;
    double Cutoff2 = Cutoff * Cutoff;
    double tElapsed1 = 0;
    double tElapsed2 = 0;
    cout << "Start simulation:" << endl;
    auto TimeStart = chrono::steady_clock::now();
    auto TimeNow = chrono::steady_clock::now();
    InitialVelo(mol); // generating initial random
    for (int step = 1; step <= StepNum; step++) { // start simulation
        VeloVerlet(mol, HalfStepSize, StepSize, 1); // velocity verlet step 1
        Force(mol, Cutoff2); // calculating force
        VeloVerlet(mol, HalfStepSize, StepSize, 2); // velocity verlet step 2
        Periodic(mol, length); // periodic boundary conditions
        ScaleVelo(mol); // velocity scaling
        if (iftraj == true && step % 50 == 0) { // output trajectory file if wanted
            if (step == 50) {
                WriteCoord(mol, "traj", 1, step); // if it is the first step, clear traj file first
            } else {
                WriteCoord(mol, "traj", 0, step); // if it is not the first step,append to the previous file
            }
        }
        if (step % 100 == 0) { // timing module, output time every 100 frames
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