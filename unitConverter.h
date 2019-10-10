/*
*      DATE:    OCTOBER 09, 2019
*        BY:    LatticeX
*   VERSION:    <1.0>
* CASE NAME:    LBM UNIT CONVERTER IN C++ (ABSTRACTED)
*  FEATURES:    MOST ACCURATE AND EASY LBM UNIT CONVERTER IN PUBLIC
*   REMARKS:    ABSTRACTED UNIT CONVERTER
*               SUITABLE FOR LOW MACH NUMBER
*               CAN BE INTEGRATED INTO LBM CODE
*               REFER unitConverterTest.cpp FOR USAGE
*/

/*
*Conversion between physical and lattice units, as well as discretization.
* We distingish between physical (dimensioned) and lattice (dimensionless) values.
* A specific conversion factor maps the two different scopes,
* e.g. physLength = conversionLength * latticeLength
*/

#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H

#include <iostream>
#include <stdlib.h>
#include <cmath>


//====================  DO NOT EDIT THIS SECTION  ====================//
// LATTICE BOLTZMANN MODEL PARAMETERS WITH FIXED VALUES
extern const double DELTA_L_LATTICE         =   1.0;                                                // UNIT LATTICE SPACE
extern const double DELTA_T_LATTICE         =   1.0;                                                // UNIT LATTICE TIME
extern const double DELTA_C_LATTICE         =   DELTA_L_LATTICE/DELTA_T_LATTICE;                    // LATTICE SPEED
extern const double CS_LATTICE              =   DELTA_C_LATTICE*sqrt(1.0/3.0);                      // SOUND SPEED IN LATTICE UNIT, FROM CS_LATTICE = √(γRT)

extern const double GAMMA_LATTICE           =   1.0;                                                // ADIABATIC CONSTANT IN LATTICE UNIT OR NO UNIT
extern const double R_GAS_LATTICE           =   1.0;                                                // SPECIFIC GAS CONSTANT IN LATTICE UNIT
extern const double DENSITY_LATTICE         =   0.22;                                               // DENSITY IN LATTICE UNIT, HARD CODED, FOLLOWING POWERFLOW (ONLY AFFECT CONVERSION FACTOR OF DENSITY&PRESSURE)

extern const double PRESSURE_LATTICE        =   CS_LATTICE*CS_LATTICE*DENSITY_LATTICE/GAMMA_LATTICE;// PRESSURE IN LATTICE UNIT (CS^2 = GAMMA*P/RHO)
//extern const double PRESSURE_LATTICE      =   R_GAS_LATTICE*DENSITY_LATTICE*TEMPERATURE_LATTICE;  // PRESSURE IN LATTICE UNIT (P = RHO*R_GAS*T)
extern const double TEMPERATURE_LATTICE     =   CS_LATTICE*CS_LATTICE/GAMMA_LATTICE/R_GAS_LATTICE;  // TEMPERATURE IN LATTICE UNIT (CS^2 = GAMMA*R_GAS*T)
//====================================================================//

//===================  PHYSICAL CONSTANTS DEFAULT  ===================//
// PHYSICAL PARAMETERS WITH DEFAULT VALUES
extern const double GAMMA_PHYS              =   1.4;                                                // ADIABATIC CONSTANT IN PHYSICAL UNIT OR NO UNIT
extern const double R_GAS_PHYS              =   287;                                                // SPECIFIC GAS CONSTANT IN PHYSICAL UNIT [J/(KG*DEG_K)]
extern const double PRESSURE_PHYS           =   101325;                                             // PRESSURE IN PHYSICAL UNIT [PA]
extern const double TEMPERATURE_PHYS        =   300;                                                // TEMPERATURE IN PHYSICAL UNIT [DEG_K]
//====================================================================//

using namespace std;

template <typename T>
class UnitConverter {
public:
  /** Documentation of constructor:
    *  \param physDeltaX                spacing between two lattice cells in [m]
    *  \param physDeltaT                time step in [sec]
    *  \param charPhysLength            reference/characteristic length of simulation geometry in [m]
    *  \param charPhysVelocity          maximal or highest expected velocity during simulation in [m/sec]
    *  \param physDensity               physical density in [kg/m^3]
    *  \param physViscosity             physical kinematic viscosity in [m^2/sec] 
    *  \param charPhysPressure          reference/characteristic physical pressure in [Pa] or [N/m^2] or [kg/(m*sec^2)]
    *  \param charPhysTemperature       reference/characteristic physical temperature in [degK]; degC = (degK - 273.15)
    *  \param gamma                     adiabatic constant
    *  \param r_gas                     specific gas constant in [J/(kg*degK)] or [Pa/(kg*degK/(m^3))]
    *                                   universal gas constant r_gas_u = 8.3143 [J/(mol*degK)]; gas molecular weight M_gas = 28.97 [g/mol] = 28.97 [kg/kmol]
    *                                   specific gas constant r_gas = r_gas_u/M_gas
    */

    constexpr UnitConverter(T physDeltaX, T physDeltaT, T charPhysLength, T charPhysVelocity,
                            T physDensity, T physViscosity,
                            T charPhysPressure = PRESSURE_PHYS, T charPhysTemperature = TEMPERATURE_PHYS,
                            T gamma = GAMMA_PHYS, T r_gas = R_GAS_PHYS)
        :   _conversionLength(physDeltaX),
            _conversionTime(physDeltaT),
            _conversionVelocity(_conversionLength/_conversionTime),
            _conversionDensity(physDensity/DENSITY_LATTICE),
            _conversionMass(_conversionDensity*pow(_conversionLength,3)),
            _conversionViscosity(_conversionLength*_conversionLength/_conversionTime),
            _conversionForce(_conversionMass*_conversionLength/(_conversionTime*_conversionTime)),
            _conversionPressure(_conversionForce/pow(_conversionLength,2)/gamma),
            _conversionForceDensity(_conversionForce/pow(_conversionLength,3)),
            _charPhysLength(charPhysLength),
            _charPhysVelocity(charPhysVelocity),
            _physDensity(physDensity),
            _physViscosity(physViscosity),
            _charPhysPressure(charPhysPressure),
            _charPhysTemperature(charPhysTemperature),
            _gamma(gamma),
            _r_gas(r_gas),
            _resolution((int)(_charPhysLength/_conversionLength+0.5)),
            _latticeRelaxationTime((_physViscosity/_conversionViscosity)/CS_LATTICE/CS_LATTICE/DELTA_T_LATTICE+0.5),
            _charLatticeVelocity(_charPhysVelocity/_conversionVelocity)
    {
    }

    virtual ~UnitConverter() = default;
    
    // nice terminal output for conversion factors, characteristical and physical data
    virtual void print() const;

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // length related
    // returns grid spacing (voxel length) in [m]
    constexpr T getPhysDeltaX() const
    {
        return _conversionLength;
    }
    // access (read-only) to private member variable
    constexpr T getConversionFactorLength() const
    {
        return _conversionLength;
    }
    // return characteristic length in physical units
    constexpr T getCharPhysLength() const
    {
        return _charPhysLength;
    }
    // conversion from lattice to  physical length
    constexpr T getPhysLength(int latticeLength) const
    {
        return _conversionLength*latticeLength;
    }
    // return resolution
    constexpr int getResolution() const
    {
        return _resolution;
    }
    // conversion from physical to lattice length, returns number of voxels for given physical length
    constexpr int getLatticeLength(T physLength) const
    {
        return int(physLength/_conversionLength + 0.5);
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // time related
    // returns time spacing (timestep length) in [sec]
    constexpr T getPhysDeltaT() const
    {
        return _conversionTime;
    }
    // access (read-only) to private member variable
    constexpr T getConversionFactorTime() const
    {
        return _conversionTime;
    }
    // conversion from lattice to  physical time
    constexpr T getPhysTime(int latticeTime) const
    {
        return _conversionTime*latticeTime;
    }
    // conversion from physical to lattice time
    constexpr int getLatticeTime(T physTime) const
    {
        return int(physTime/_conversionTime+0.5);
    }
    // return relaxation time in lattice units
    constexpr T getLatticeRelaxationTime() const
    {
        return _latticeRelaxationTime;
    }
    // return relaxation frequency in lattice units
    constexpr T getLatticeRelaxationFrequency() const
    {
        return 1.0/_latticeRelaxationTime;
    }
    // return relaxation frequency in lattice units computed from given physical diffusivity in [m^2/sec]
    template <typename DESCRIPTOR_>
    constexpr T getLatticeRelaxationFrequencyFromDiffusivity(const T physDiffusivity) const
    {
        T latticeDiffusivity = physDiffusivity/_conversionViscosity;
        return 1.0/(latticeDiffusivity/CS_LATTICE/CS_LATTICE/DELTA_T_LATTICE+0.5);
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // velocity related
    // access (read-only) to private member variable
    constexpr T getConversionFactorVelocity() const
    {
        return _conversionVelocity;
    }
    // return characteristic velocity in physical units
    constexpr T getCharPhysVelocity() const
    {
        return _charPhysVelocity;
    }
    // return sound speed in physical units
    constexpr T getPhysSoundSpeed() const
    {
        return sqrt(_gamma*_r_gas*_charPhysTemperature);
    }
    // conversion from lattice to  physical velocity
    constexpr T getPhysVelocity(T latticeVelocity) const
    {
        return _conversionVelocity*latticeVelocity;
    }
    // return characteristic velocity in lattice units
    constexpr T getCharLatticeVelocity() const
    {
        return _charLatticeVelocity;
    }
    // return sound speed in lattice units
    constexpr T getLatticeSoundSpeed() const
    {
        return CS_LATTICE;
    }
    // conversion from physical to lattice velocity
    constexpr T getLatticeVelocity(T physVelocity) const
    {
        return physVelocity/_conversionVelocity;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // density related
    // access (read-only) to private member variable
    constexpr T getConversionFactorDensity() const
    {
        return _conversionDensity;
    }
    // return density in physical units
    constexpr T getPhysDensity() const
    {
        return _physDensity;
    }
    // conversion from lattice to  physical density
    constexpr T getPhysDensity(T latticeDensity) const
    {
        return _conversionDensity*latticeDensity;
    }
    // return density in lattice units
    constexpr T getLatticeDensity() const
    {
        return _physDensity/_conversionDensity;
    }
    // conversion from physical to lattice density
    constexpr T getLatticeDensity(T physDensity) const
    {
        return physDensity/_conversionDensity;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // mass related
    // access (read-only) to private member variable
    constexpr T getConversionFactorMass() const
    {
        return _conversionMass;
    }
    // conversion from lattice to  physical mass
    constexpr T getPhysMass(T latticeMass) const
    {
        return _conversionMass*latticeMass;
    }
    // conversion from physical to lattice mass
    constexpr T getLatticeMass(T physMass) const
    {
        return physMass/_conversionMass;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // viscosity related
    // access (read-only) to private member variable
    constexpr T getConversionFactorViscosity() const
    {
        return _conversionViscosity;
    }
    // return viscosity in physical units
    constexpr T getPhysViscosity() const
    {
        return _physViscosity;
    }
    // conversion from lattice to  physical viscosity
    constexpr T getPhysViscosity(T latticeViscosity) const
    {
        return _conversionViscosity*latticeViscosity;
    }
    // return viscosity in lattice units
    constexpr T getLatticeViscosity() const
    {
        return _physViscosity/_conversionViscosity;
    }
    // conversion from physical to lattice viscosity
    constexpr T getLatticeViscosity(T physViscosity) const
    {
        return physViscosity/_conversionViscosity;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // force related
    // access (read-only) to private member variable
    constexpr T getConversionFactorForce() const
    {
        return _conversionForce;
    }
    // conversion from lattice to  physical force
    constexpr T getPhysForce(T latticeForce) const
    {
        return _conversionForce*latticeForce;
    }
    // conversion from physical to lattice force
    constexpr T getLatticeForce(T physForce) const
    {
        return physForce/_conversionForce;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // pressure ralated
    // access (read-only) to private member variable
    constexpr T getConversionFactorPressure() const
    {
        return _conversionPressure;
    }
    // return characteristic pressure in physical units
    constexpr T getCharPhysPressure() const
    {
        return _charPhysPressure;
    }
    // conversion from lattice to  physical pressure
    constexpr T getPhysPressure(T latticePressure) const
    {
        return _conversionPressure*latticePressure;
    }
    // conversion from physical to lattice pressure
    constexpr T getLatticePressure( T physPressure ) const
    {
        return (physPressure)/_conversionPressure;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // force-density related
    // access (read-only) to private member variable
    constexpr T getConversionFactorForceDensity() const
    {
        return _conversionForceDensity;
    }
    // conversion from lattice to  physical force
    constexpr T getPhysForceDensity(T latticeForceDensity) const
    {
        return _conversionForceDensity*latticeForceDensity;
    }
    // conversion from physical to lattice force
    constexpr T getLatticeForceDensity(T physForceDensity) const
    {
        return physForceDensity/_conversionForce;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // temperature related
    // return characteristic temperature in physical units
    constexpr T getCharPhysTemperature() const
    {
        return _charPhysTemperature;
    }

    // ----------~~~~~~~~~~----------~~~~~~~~~~----------~~~~~~~~~~----------
    // numbers related
    // return Reynolds number
    constexpr T getReynoldsNumber() const
    {
        return _charPhysVelocity*_charPhysLength/_physViscosity;
    }
    // return Mach number
    constexpr T getMachNumber() const
    {
        return _charPhysVelocity/getPhysSoundSpeed();
    }
    // return Knudsen number
    constexpr T getKnudsenNumber() const
    {
        // This calculates the lattice Knudsen number.
        // See e.g. (7.22) in "The Lattice Boltzmann Method: Principles and Practice" [kruger2017lattice].
        return getMachNumber()/getReynoldsNumber();
    }

protected:
    // conversion factors
    const T _conversionLength;          // [m]
    const T _conversionTime;            // [sec]
    const T _conversionVelocity;        // [m/s]
    const T _conversionDensity;         // [kg/m^3]
    const T _conversionMass;            // [kg]
    const T _conversionViscosity;       // [m^2/sec]
    const T _conversionForce;           // [N] or [kg*m/sec^2]
    const T _conversionPressure;        // [Pa] or [N/m^2] or [kg/(m*sec^2)]
    const T _conversionForceDensity;    // [N/m^3] or [Pa/m]; force density, force per volume, pressure-gradient

    // physical units, e.g characteristic or reference values
    const T _charPhysLength;            // [m]
    const T _charPhysVelocity;          // [m/sec]
    const T _physDensity;               // [kg/m^3]
    const T _physViscosity;             // [m^2/sec]
    const T _charPhysPressure;          // [N/m^2] or [kg/(m*sec^2)]
    const T _charPhysTemperature;       // [degK], char temperature, degC = (degK - 273.15)

    // physicall parameters
    const T _gamma;                     // adiabatic constant
    const T _r_gas;                     // [Pa/(kg*degK/(m^3))], [J/(kg*degK)], specific gas constant for air

    // lattice units, discretization parameters
    const int _resolution;
    const T _latticeRelaxationTime;
    const T _charLatticeVelocity; 

}; // End of UnitConverter Class Defination

template <typename T>
void UnitConverter<T>::print() const
{

    cout << "==================== UnitConverter Information ====================" << endl;
    cout << endl;

    cout << "-------------------------------------------------------------------" << endl;
    cout << "PHYSICAL PARAMETERS ..." << endl;
    cout << "characteristic pressure [Pa]:          charPressure    =   " << getCharPhysPressure() << endl;
    cout << "characteristic temperature [degK]:     charTemperature =   " << getCharPhysTemperature() << endl;
    cout << "charactreistic length [m]:             charL           =   " << getCharPhysLength() << endl;
    cout << "charactreistic velocity [m/sec]:       charU           =   " << getCharPhysVelocity() << endl;
    cout << "physical density [kg/m^3]:             charRho         =   " << getPhysDensity() << endl;
    cout << "physical kinematic viscoity [m^2/sec]: charNu          =   " << getPhysViscosity() << endl;
    cout << "voxel length [m]:                      physDeltaX      =   " << getConversionFactorLength() << endl;
    cout << "time step [sec]:                       physDeltaT      =   " << getConversionFactorTime() << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << endl;

    cout << "-------------------------------------------------------------------" << endl;
    cout << "LATTICE PARAMETERS ..." << endl;
    cout << "resolution:                            N               =   " << getResolution() << endl;
    cout << "lattice velocity:                      latticeU        =   " << getCharLatticeVelocity() << endl;
    cout << "lattice relaxation time:               tau             =   " << getLatticeRelaxationTime() << endl;
    cout << "lattice relaxation frequency:          omega           =   " << getLatticeRelaxationFrequency() << endl;
    cout << "lattice kinematic viscoity:            latticeNu       =   " << getLatticeViscosity() << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << endl;

    cout << "-------------------------------------------------------------------" << endl;
    cout << "DIMLESS NUMBERS ..." << endl;
    cout << "mach number:                           machNumber      =   " << getMachNumber() << endl;
    cout << "reynolds number:                       reynoldsNumber  =   " << getReynoldsNumber() << endl;
    cout << "knudsen number:                        knudsenNumber   =   " << getKnudsenNumber() << endl;
    cout << "-------------------------------------------------------------------" << endl;
    cout << endl;

    cout << "==================== +++++++++++++++++++++++++ ====================" << endl;

}

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: physViscosity = charPhysLength*charPhysVelocity/reynoldsNumber
template <typename T>
class UnitConverterFromPhysVelocityAndReynoldsNumber : public UnitConverter<T> {
public:
    constexpr UnitConverterFromPhysVelocityAndReynoldsNumber(
    int resolution,
    T charPhysLength,
    T charPhysVelocity,
    T reynoldsNumber,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (charPhysVelocity),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        (charPhysLength*charPhysVelocity/reynoldsNumber),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: charPhysVelocity = charLatticeVelocity*conversionVelocity = charLatticeVelocity*(physSoundSpeed/latticeSoundSpeed)
// Remark: physViscosity = charPhysLength*charPhysVelocity/reynoldsNumber
template <typename T>
class UnitConverterFromLatticeVelocityAndReynoldsNumber : public UnitConverter<T> {
public:
    constexpr UnitConverterFromLatticeVelocityAndReynoldsNumber(
    int resolution,
    T charPhysLength,
    T charLatticeVelocity,
    T reynoldsNumber,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (charLatticeVelocity*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        (charPhysLength*(charLatticeVelocity*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE))/reynoldsNumber),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: charPhysVelocity = machNumber*physSoundSpeed
// Remark: physViscosity = charPhysLength*(machNumber*physSoundSpeed)/reynoldsNumber
template <typename T>
class UnitConverterFromMachNumberAndReynoldsNumber : public UnitConverter<T> {
public:
    constexpr UnitConverterFromMachNumberAndReynoldsNumber(
    int resolution,
    T charPhysLength,
    T machNumber,
    T reynoldsNumber,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (machNumber*sqrt(gamma*r_gas*charPhysTemperature)),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        (charPhysLength*(machNumber*sqrt(gamma*r_gas*charPhysTemperature))/reynoldsNumber),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: physViscosity = conversionLength*conversionVelocity*latticeViscosity = (charPhysLength/resolution)*(physSoundSpeed/latticeSoundSpeed)*latticeViscosity
// Remark: latticeViscosity = (latticeRelaxationTime-0.5)*latticeSoundSpeed*latticeSoundSpeed*latticeDeltaT
// Remark: If specify lattice relaxation time, the Reynolds number cannot be guaranteed since the physViscosity will be changed (if want guanratee Reynolds number, charPhysVelocity should change as well as physViscosity)
template <typename T>
class UnitConverterFromPhysVelocityAndRelaxationTime : public UnitConverter<T> {
public:
    constexpr UnitConverterFromPhysVelocityAndRelaxationTime(
    int resolution,
    T charPhysLength,
    T charPhysVelocity,
    T latticeRelaxationTime,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (charPhysVelocity),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        ((charPhysLength/resolution)*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)*((latticeRelaxationTime-0.5)*CS_LATTICE*CS_LATTICE*DELTA_T_LATTICE)),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: charPhysVelocity = charLatticeVelocity*conversionVelocity = charLatticeVelocity*(physSoundSpeed/latticeSoundSpeed)
// Remark: physViscosity = conversionLength*conversionVelocity*latticeViscosity = (charPhysLength/resolution)*(physSoundSpeed/latticeSoundSpeed)*latticeViscosity
// Remark: latticeViscosity = (latticeRelaxationTime-0.5)*latticeSoundSpeed*latticeSoundSpeed*latticeDeltaT
// Remark: If specify lattice relaxation time, the Reynolds number cannot be guaranteed since the physViscosity will be changed (if want guanratee Reynolds number, charPhysVelocity should change as well as physViscosity)
template <typename T>
class UnitConverterFromLatticeVelocityAndRelaxationTime : public UnitConverter<T> {
public:
    constexpr UnitConverterFromLatticeVelocityAndRelaxationTime(
    int resolution,
    T charPhysLength,
    T charLatticeVelocity,
    T latticeRelaxationTime,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (charLatticeVelocity*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        ((charPhysLength/resolution)*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)*((latticeRelaxationTime-0.5)*CS_LATTICE*CS_LATTICE*DELTA_T_LATTICE)),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
// Remark: charPhysVelocity = machNumber*physSoundSpeed
// Remark: physViscosity = conversionLength*conversionVelocity*latticeViscosity = (charPhysLength/resolution)*(physSoundSpeed/latticeSoundSpeed)*latticeViscosity
// Remark: latticeViscosity = (latticeRelaxationTime-0.5)*latticeSoundSpeed*latticeSoundSpeed*latticeDeltaT
// Remark: If specify lattice relaxation time, the Reynolds number cannot be guaranteed since the physViscosity will be changed (if want guanratee Reynolds number, charPhysVelocity should change as well as physViscosity)
template <typename T>
class UnitConverterFromMachNumberAndRelaxationTime : public UnitConverter<T> {
public:
    constexpr UnitConverterFromMachNumberAndRelaxationTime(
    int resolution,
    T charPhysLength,
    T machNumber,
    T latticeRelaxationTime,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (machNumber*sqrt(gamma*r_gas*charPhysTemperature)),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        ((charPhysLength/resolution)*(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)*((latticeRelaxationTime-0.5)*CS_LATTICE*CS_LATTICE*DELTA_T_LATTICE)),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

// Remark: UnitConverter(physDeltaX, physDeltaT, charPhysLength, charPhysVelocity, physDensity, physViscosity, charPhysPressure, charPhysTemperature, gamma, r_gas)
// Remark: physDeltaT = conversionLength/conversionVelocity = (charPhysLength/resolution)/(physSoundSpeed/latticeSoundSpeed)
// Remark: physSoundSpeed = sqrt(gamma*r_gas*charPhysTemperature)
// Remark: physDensity = charPhysPressure/(r_gas*charPhysTemperature)
// Remark: Usually Only Update <charPhysVelocity>, <physViscosity>
// Remark: ----------~~~~~~~~~~----------~~~~~~~~~~----------
template <typename T>
class UnitConverterFromPhysVelocityAndPhysViscosity : public UnitConverter<T> {
public:
    constexpr UnitConverterFromPhysVelocityAndPhysViscosity(
    int resolution,
    T charPhysLength,
    T charPhysVelocity,
    T physViscosity,
    T charPhysPressure = PRESSURE_PHYS,
    T charPhysTemperature = TEMPERATURE_PHYS,
    T gamma = GAMMA_PHYS,
    T r_gas = R_GAS_PHYS) : UnitConverter<T>(
        (charPhysLength/resolution),
        (charPhysLength/resolution/(sqrt(gamma*r_gas*charPhysTemperature)/CS_LATTICE)),
        (charPhysLength),
        (charPhysVelocity),
        (charPhysPressure/(r_gas*charPhysTemperature)),
        (physViscosity),
        (charPhysPressure),
        (charPhysTemperature),
        (gamma),
        (r_gas))
    {
    }

};

#endif