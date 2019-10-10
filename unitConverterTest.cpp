#include <iostream>
#include <stdlib.h>
#include <cmath>
#include "unitConverter.h"

using namespace std;

typedef double T;

#define converterFromPhysVelocityAndReynoldsNumber
#define converterFromLatticeVelocityAndReynoldsNumber
#define converterFromMachNumberAndReynoldsNumber
#define converterFromPhysVelocityAndRelaxationTime
#define converterFromLatticeVelocityAndRelaxationTime
#define converterFromMachNumberAndRelaxationTime
#define converterFromPhysVelocityAndPhysViscosity

int main()
{

    #if defined(converterFromPhysVelocityAndReynoldsNumber)
    UnitConverterFromPhysVelocityAndReynoldsNumber<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     100,        // charPhysVelocity: maximal/highest expected velocity during simulation in [m/sec]
        (T)     100         // reynoldsNumber: specified Reynolds number
    );
    #elif defined(converterFromLatticeVelocityAndReynoldsNumber)
    UnitConverterFromLatticeVelocityAndReynoldsNumber<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     0.166293,   // charLatticeVelocity: maximal/highest expected velocity during simulation in lattice units
        (T)     100         // reynoldsNumber: specified Reynolds number
    );
    #elif defined(converterFromMachNumberAndReynoldsNumber)
    UnitConverterFromMachNumberAndReynoldsNumber<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     0.288028,   // machNumber: specified Mach number
        (T)     100         // reynoldsNumber: specified Reynolds number
    );
    #elif defined(converterFromPhysVelocityAndRelaxationTime)
    UnitConverterFromPhysVelocityAndRelaxationTime<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     100,        // charPhysVelocity: maximal/highest expected velocity during simulation in [m/sec]
        (T)     0.998879    // latticeRelaxationTime: specified lattice relaxation time
    );
    #elif defined(converterFromLatticeVelocityAndRelaxationTime)
    UnitConverterFromLatticeVelocityAndRelaxationTime<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     0.166293,   // charLatticeVelocity: maximal/highest expected velocity during simulation in lattice units
        (T)     0.998879    // latticeRelaxationTime: specified lattice relaxation time
    );
    #elif defined(converterFromMachNumberAndRelaxationTime)
    UnitConverterFromMachNumberAndRelaxationTime<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     0.288028,   // machNumber: specified Mach number
        (T)     0.998879    // latticeRelaxationTime: specified lattice relaxation time
    );
    #elif defined(converterFromPhysVelocityAndPhysViscosity)
    UnitConverterFromPhysVelocityAndPhysViscosity<T> converter(
        (int)   100,        // resolution: number of voxels per charPhysLength
        (T)     1.0,        // charPhysLength: reference length of simulation geometry in [m]
        (T)     100,        // charPhysVelocity: maximal/highest expected velocity during simulation in [m/sec]
        (T)     1.0         // physViscosity: specified physical kinematic viscosity in [m^2/sec] 
    );

    #endif

    // Prints the converter log as console output
    converter.print();

    cout << endl;
    cout << "The following results are intended to test the outputs functions!" << endl;
    cout << endl;
    const double N = converter.getResolution();
    cout << "input resolution: " << N << endl;

    const double latticeDensity = converter.getLatticeDensity();
    cout << "lattice density: " << latticeDensity << endl;

    const double latticeVelocity = converter.getLatticeVelocity(2.0*converter.getCharPhysVelocity());
    cout << "calculated lattice velocity doubling: " << latticeVelocity << endl;

    return 0;
}
