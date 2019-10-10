/*
*      DATE:    OCTOBER 09, 2019
*        BY:    LatticeX
*   VERSION:    <1.0>
* CASE NAME:    LBM UNIT CONVERTER IN C++ (STANDALONE)
*  FEATURES:    MOST ACCURATE AND EASY LBM UNIT CONVERTER IN PUBLIC
*   REMARKS:    STANDALONE UNIT CONVERTER
*               SUITABLE FOR LOW MACH NUMBER
*               CAN BE USED FOR CALCULATION
*/

#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;


//====================  DO NOT EDIT THIS SECTION  ====================//

const double DELTA_L_LATTICE        =   1.0;                                                // UNIT LATTICE SPACE
const double DELTA_T_LATTICE        =   1.0;                                                // UNIT LATTICE TIME
const double DELTA_C_LATTICE        =   DELTA_L_LATTICE/DELTA_T_LATTICE;                    // LATTICE SPEED
const double CS_LATTICE             =   DELTA_C_LATTICE*sqrt(1.0/3.0);                      // SOUND SPEED IN LATTICE UNIT, FROM CS_LATTICE = √(γRT)

const double GAMMA_LATTICE          =   1.0;                                                // ADIABATIC CONSTANT IN LATTICE UNIT, REDUCED VALUE TO REDUCE MACH NUMBER
const double R_GAS_LATTICE          =   1.0;                                                // SPECIFIC GAS CONSTANT IN LATTICE UNIT
const double DENSITY_LATTICE        =   1.0;                                                // DENSITY IN LATTICE UNIT, HARD CODED (ONLY AFFECT CONVERSION FACTOR OF DENSITY&PRESSURE)

const double TEMPERATURE_LATTICE    =   CS_LATTICE*CS_LATTICE/GAMMA_LATTICE/R_GAS_LATTICE;  // TEMPERATURE IN LATTICE UNIT (CS^2 = GAMMA*R_GAS*T)
const double PRESSURE_LATTICE       =   CS_LATTICE*CS_LATTICE*DENSITY_LATTICE/GAMMA_LATTICE;// PRESSURE IN LATTICE UNIT (CS^2 = GAMMA*P/RHO)
//const double PRESSURE_LATTICE     =   R_GAS_LATTICE*DENSITY_LATTICE*TEMPERATURE_LATTICE;  // PRESSURE IN LATTICE UNIT (P = RHO*R_GAS*T)


//====================  START EDIT FROM HERE  ====================//

// STEP 0: NON-CRITICAL PARAMETERS, BETTER NOT CHANGE
const double gamma                  =   1.4;                                                // adiabatic constant
const double r_gas                  =   287;                                                // [Pa/(kg*degK/(m^3))], [J/(kg*degK)], specific gas constant for air
const double char_temperature       =   300;                                                // [degK], char temperature, = (300 - 273.15) = 26.85 degC

const double sound_speed            =   sqrt(gamma*r_gas*char_temperature);                 // [m/sec], sound speed
const double C_U                    =   sound_speed/CS_LATTICE;                             // [M/SEC], CONVERSION FACTOR OF VELOCITY

// STEP 1: CRITICAL PARAMETERS, CHANGE FREELY
const unsigned int resolution       =   100;                                                // resolution (cells along char length)
const double char_length            =   1.0;                                                // [m], char length
const double char_pressure          =   101325;                                             // [Pa], char pressure

// STEP 2: SPECIFY PHYSICAL VELOCITY OR LATTICE VELOCITY
const double char_velocity          =   1.0;                                                // [m/sec], char velocity
const double char_velocity_lattice  =   char_velocity/C_U;                                  // velocity in lattice unit
//const double char_velocity_lattice  =   0.1;                                              // velocity in lattice unit
//const double char_velocity          =   char_velocity_lattice*C_U;                        // [m/sec], char velocity

// STEP 3: SELECT ONE FROM THE FOLLOWING TWO METHODS AND QUOTE ANOTHER
#define reynolds_control                                                                    // specify reynolds number
//#define relaxation_time_control                                                           // specify relaxation time

// STEP 4: SPECIFY REYNOLDS NUMBERS OR RELAXATION TIME ACCORDING TO STEP 3
#if defined(reynolds_control) && !defined(relaxation_time_control)
const double reynolds               =   100;                                                // reynolds number
#elif defined(relaxation_time_control) && !defined(reynolds_control)
const double tau                    =   1.0;                                                // relaxation time
#endif


//====================  STOP ANY EDIT FROM HERE  ====================//

// DERIVED PARAMETERS
const double mach                   =   char_velocity/sound_speed;                          // mach number
const double density                =   char_pressure/(r_gas*char_temperature);             // [kg/m^3], dentsity

// DERIVED CONVERTERS
const double C_L                    =   char_length/resolution;                             // [M], CONVERSION FACTOR OF LENGTH
const double C_T                    =   C_L/C_U;                                            // [SEC], CONVERSION FACTOR OF TIME
const double C_R                    =   density/DENSITY_LATTICE;                            // [KG/M^3], CONVERSION FACTOR OF DENSITY
const double C_F                    =   C_R*C_L/(C_T*C_T);                                  // [N/M^3] or [PA/M], CONVERSION FACTOR OF FORCE PER VOLUME (FORCE DENSITY, PRESSURE-GRADIENT)
const double C_A                    =   C_L/(C_T*C_T);                                      // [M/S^2], CONVERSION FACTOR OF ACCELERATION
const double C_P                    =   C_F*C_L/gamma;                                      // [N/M^2], CONVERSION FACTOR OF PRESSURE

// OUTPUT AND CALCULATION OF MORE DERIVED PARAMETERS
int main()
{
    // OUTPUT
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "INPUT PARAMETERS ..." << endl;
    cout << "characteristic pressure:           " << char_pressure << " [Pa]" << endl;
    cout << "characteristic temperature:        " << char_temperature << " [degK]" << endl;
    cout << "characteristic physical velocity:  " << char_velocity << " [m/sec]" << endl;
    cout << "characteristic lattice velocity:   " << char_velocity_lattice << ";" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << endl;

    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "DERIVED PARAMETERS ..." << endl;
    cout << "calculated sound speed:            " << sound_speed << " [m/sec];" << endl;
    cout << "calculated mach number:            " << mach << ";" << endl;
    cout << "calculated density:                " << density << " [kg/m^3];" << endl;
    cout << "calculated physical time scaling:  " << C_T << " [sec]" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << endl;

    /* ********** CALCULATION - BASED ON REYNOLDS NUMBER ********** */
    #if defined(reynolds_control) && !defined(relaxation_time_control)
    const double viscosity = char_length*char_velocity/reynolds; // [m^2/sec], kinematic viscosity
    const double viscosity_lattice = char_velocity_lattice*resolution/reynolds; // kinematic viscosity in lattice unit
    //const double viscosity_lattice = viscosity/C_U/C_L; // kinematic viscosity in lattice unit
    const double tau = viscosity_lattice/CS_LATTICE/CS_LATTICE/DELTA_T_LATTICE+0.5; // relaxation time in lattice unit (should be lattice sec)

    // OUTPUT
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "REYNOLDS NUMBER SPECIFIED ..." << endl;
    cout << "specified reynolds number:         " << reynolds << ";" << endl;
    cout << "calculated relaxation time:        " << tau << ";" << endl;
    cout << "calculated physical viscosity:     " << viscosity << " [m^2/sec];" << endl;
    cout << "calculated lattice viscosity:      " << viscosity_lattice << ";" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << endl;

    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "CONVERSION FACTORS ..." << endl;
    cout << "CONVERSION FACTOR OF LENGTH:         " << C_L << " [M]" << endl;
    cout << "CONVERSION FACTOR OF DENSITY:        " << C_R << " [KG/M^3]" << endl;
    cout << "CONVERSION FACTOR OF PRESSURE:     " << C_P << " [N/M^2];" << endl;
    cout << "calculated physical pressure:      " << C_P*PRESSURE_LATTICE << " [N/M^2]" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << endl;

    /* ********** CALCULATION - BASED ON RELAXATION TIME ********** */
    #elif defined(relaxation_time_control) && !defined(reynolds_control)
    const double viscosity_lattice = (tau-0.5)*CS_LATTICE*CS_LATTICE*DELTA_T_LATTICE; // kinematic viscosity in lattice unit
    const double reynolds = char_velocity_lattice*resolution/viscosity_lattice; // reynolds number
    const double viscosity = char_velocity*char_length/reynolds; // [m^2/sec], kinematic viscosity
    //const double viscosity = viscosity_lattice*C_U*C_L; // [m^2/sec], kinematic viscosity

    // OUTPUT
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "RELAXATION TIME SPECIFIED ..." << endl;
    cout << "specified relaxation time:         " << tau << ";" << endl;
    cout << "calculated reynolds number:        " << reynolds << ";" << endl;
    cout << "calculated physical viscosity:     " << viscosity << " [m^2/sec];" << endl;
    cout << "calculated lattice viscosity:      " << viscosity_lattice << ";" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << endl;

    #else
    cout << "==========++~~~~~~~~~~++==========" << endl;
    cout << "WRONG METHOD! CALCULATION STOPPED!" << endl;
    cout << "==========++~~~~~~~~~~++==========" << endl;
    exit(0); // wrong calculation method defined, program terminated
    #endif

    return 0;
}
