/*
 * main.cpp
 *
 *  Created on: 04.01.2012
 *      Author: Stefan Becker <stefan.becker@mv.uni-kl.de>
 *
 *  The method main() is called via the command line allowing for multiple flags.
 *  Operation: 	- processing the the input information: prefix of output file, what flags are set,
 *  			  adoption of the values of the flags (if set), otherwise default values
 *  		   	- checking the input arguments for completeness and consistency
 *  			- generating class instances and calling the appropriate methods: model parameters, writing methods
 *
 */

#include "ConfigWriter.h"
#include "PhaseSpaceWriter.h"
#include "Component.h"


extern const std::string FLUID_AR;
extern const std::string FLUID_CH4;
extern const std::string WALL_CU_LJ; 		// copper by Lennard-Jones
extern const std::string WALL_TERSOFF;
extern const std::string WALL_MAT_AKA_FIT;
extern const unsigned short THERMOSTAT_VELSCALE;
extern const unsigned short THERMOSTAT_ANDERSEN;



//@todo: Konstanten sauber als "extern const float DEFAULT_ALPHA = 3" definieren und in einer header-Datei deklarieren
const double DEFAULT_ALPHA  = 3;
const double DEFAULT_BETA = 1.8;
const double DEFAULT_GAMMA = 3;
const double DEFAULT_EDGE = 1;
const double PI = 3.1415927;
const unsigned DEFAULT_PROFILE_PHI = 1;
const unsigned DEFAULT_PROFILE_R = 200;
//const unsigned DEFAULT_PROFILE_H = 200;



int main(int argc, char** argv){
	const char* usage = "usage: \nmkSD -p <prefix> [-a <alpha>] [-b <beta>] [-g <gamma>] [-c initCanonical] [-d <vapor density factor>] [-e <eta12 Lorentz>] [-f <fluid>] [-h <edge_proportion>] [-m] \n"
			"-N <fluid_number> [-o <output time steps>] [-P <n_phi n_r n_h>] [-s <sigma wall-wall>] [-S <number of stripes>] -T <temperatur>  [-u <LJ units>] "
			"-w <wall thickness> \n[-W <wall interaction model>] -x12 <xi_12 Berthelot> [-x13 <xi_13 Berthelot>]\n\n"
			"-a \t alpha: box length in x- direction as a mulitple of the liquid cuboid length in x-direction, default: alpha == 3\n"
			"-b \t beta: effective width of the fluid space in y-direction as a multiple of the liquid cuboid length in y-direction, default: beta ==1.8 \n"
			"-g \t gamma: box length in z-direction as a mulitple of the liquid cuboid length in z-direction, default: gamma == 3\n"
			"-c \t Init Canonical: number of time steps\n"
			"-d \t vapor density factor: correction factor for the vapor density (initial thought: metastable states?) \n"
			"-e \t eta_12 according to the extended Lorentz combining rule. \n"
			"-f \t fluid: choosing different types of fluids, Ar is the default, other options: CH4, C2H6, N2, CO2, C6H14. So far not implemented! \n"
			"-h \t liquid cuboid: edge lengths in x,z-direction as a multiple of the edge length in y-direction.\n"
			"-m \t movie option: For a movie the timestep lenth is increased and the number of output time steps is set to 500. \n"
			"-N \t number of fluid particles. \n"
			"-o \t number of output timesteps.\n"
			"-P \t number of profile units employed for recording the density profile.\n"
			"-s \t sigma of the wall-wall interaction => obsolete, not to be used! \n"
			"-S \t choosing the wall to be structured in a stripes shaped manner, the number of stripes is bounded by the ratio of stripe width to the lattice constant. \n"
			"-t \t Thermostat: -t factor == Andersen thermostat with factor * nu, defalt: Velocity Scaling \n"
			"-T \t temperature in atomic units \n"
			"-u \t using the Lennard-Jones units set. \n"
			"-w \t wall thickness as (w-0.5) times the lattice constant of the solid matter. In fact: -w = 1 => h_Wall = 0.5 * lattice constant;  -w 2 => h_Wall = 1.5*lattice constant,\n"
			"-W \t choosing a wall interaction model => So far the LJTS is the only one implemented. \n"
			"-x12 \t the Berthelot combining rule interaction parameter xi between the components 1 and 2. \n"
			"-x13 \t the Berthelot combining rule interaction parameter xi between the components 1 and 3. \n"
			"\n*********************************************************************************************************\n\n";
	if((argc<10)||(argc>35)){
		std::cout << "\n\n*********************************************************************************************************\n"
		<< "There are "<< floor(0.5*(argc-1)) << " complete arguments where 5 to 19 should be given.\n\n";
		std::cout << usage ;
		return 1;
	}

// initializing the c++ pseudo-random number generator
srand(time(NULL));

//flags set? =>default: not set	// flags:
bool in_prefix = false;			// -p => prefix of the output file
bool in_alpha = false;			// -a => width of the simulation box in x-direction as a multiple of the liquid cuboid width (in x-direction), allowed stepwidth: 0.1
bool in_beta = false;			// -b => height of the simulation box in y-direction as a multiple of the liquid cuboid height (in y-direction), stepwidth: 0.1
bool in_gamma = false;			// -g => width of the simulation box in z-direction as a multiple of the liquid cuboid width (in z-direction), allowed stepwidth: 0.1
bool in_initCanon = false;		// -c => number of time steps used for init canonical => simulated annealing
bool in_density = false;		// -d => factor changing the vapor density in the start configuration, allowing for start density different from Kedia et al.
bool in_eta12 = false;			// -e => eta12_fluid_wall of the Lorentz combinig rule
bool in_fluid = false;			// -f => type of fluid: Argon by default
bool in_edgeProp  = false;		// -h => edge proportions of the liquid cuboid, default: 1
bool movie = false;			// -m => movie of the system's time evolution intended? => output of VisittWriter
bool in_N = false;			// -N => number of fluid particles
bool in_outputTime = false;		// -o => number of time steps after that the density profile is written out
bool in_numProfileUnits = false;// -P => numer of profile units (bins) per direction, in the cylindrical coordinate system
bool in_sigWall = false;		// -s => sigma_wall_wall: (i) LJ-wall changes sigma of the wall, (ii) otherwise (e.g. in case of harmonic oszillators) the equilibrium bond length
bool stripes = false;			// -S => strpies shaped plane wall, the stripes exhibit different values of xi
bool in_xi12 = false;			// -xi12 => the standard argument to specify the xi of the Berthelot combining rule for two components
bool in_xi13 = false;			// -xi13 => to be set only if the flag -S (stripes shaped wall) is set. Then, (kind of)three components are present!
bool in_temperature = false;	// -T => Temperature in atomic units
bool LJunits = true; 			// -u => input in Lennard-Jones units
bool in_wallThick = false;		// -w => wall thickness as a multiple of 0.5 times cristal lattice constant (fcc cristal presumed!) @todo: realization/significance for different kinds of cristal lattices
bool in_wallModel = false;		// -W => choice of models describing the solid wall: first LJ-wall (fcc), representing copper, @todo: other models to be implemented!!!
//bool in_xi = false;			// -x => xi_fluid_wall of the Berthelot combining rule
bool LJShifted = true;			// by default a LJTS potential is used


// declaration of global variables
unsigned N;
unsigned initCanon;
unsigned profilePhi, profileR, profileH; // number of discrete (density-) profile elements in the corresponding direction
unsigned numberOfStripes;
unsigned short thermostat = THERMOSTAT_VELSCALE;	// ==1 --> applying velocity scaling, ==2 --> applying the Andersen Thermostat
int outTimeSteps, wallThick;
double alpha, beta, gamm, edgeProp, densFac, eta12, sigWall, temp, xi12, xi13, nuFactor;
char* prefix = (char*)0;		// name of the output file as a C-string, initialized by NULL pointer
std::string fluid, wall, prefixStr;

// processing the arguments put in the command line
for(int i = 1; i < argc; i++){
	if(*argv[i] != '-'){
		std::cout <<"\nFlag expected where '"<< argv[i] << "' was given. \n\n";
		std::cout << usage;
		return 2;
	}
	for(unsigned j=1; argv[i][j]; j++)
	{
		if(argv[i][j] == 'p')
		{
			in_prefix = true;
			i++;
			prefix = argv[i];
			break;
		}
		if(argv[i][j] == 'a')
		{
			in_alpha = true;
			i++;
			alpha = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'b')
		{
			in_beta= true;
			i++;
			beta = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'g'){
			in_gamma = true;
			i++;
			gamm = atof(argv[i]);
			if (gamm < 1.1){
				std::cout << "Value chosen for gamma too small! Choose value bigger than 1.1!\n\n";
				return 14;
			}
			break;
		}
		else if(argv[i][j] == 'c')
		{
			in_initCanon = true;
			i++;
			initCanon = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'd')
		{
			in_density = true;
			i++;
			densFac = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'e')
		{
			in_eta12 = true;
			i++;
			eta12 = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'f')
		{
			in_fluid = true;
			i++;
			if(!strcmp(argv[i], "Ar")) fluid = FLUID_AR;
			//else if(!strcmp(argv[i], "CH4")) fluid = FLUID_CH4;
			/*else if(!strcmp(argv[i], "C2H6")) fluid = FLUID_C2H6;
			else if(!strcmp(argv[i], "N2")) fluid = FLUID_N2;
			//else if(!strcmp(argv[i], "CO2")) fluid = FLUID_CO2;
			else if(!strcmp(argv[i], "C6H14")) fluid = FLUID_C6H14;*/
			else
			{
				std::cout << "Fluid " << argv[i] << "is not available. \n\n" << usage;
				return 3;
			}
			break;
		}
		else if (argv[i][j] == 'h'){
		 in_edgeProp = true;
		 i++;
		 edgeProp = atof(argv[i]);
		 break;
		}
		else if (argv[i][j] == 'm') movie = true;
		else if(argv[i][j] == 'N')
		{
			in_N = true;
			i++;
			N = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'o')
		{
			in_outputTime = true;
			i++;
			outTimeSteps = atoi(argv[i]);
			break;
		}
		else if(argv[i][j] == 'P'){
			in_numProfileUnits = true;
			i++;
			profilePhi = atoi(argv[i]);
			i++;
			profileR = atoi(argv[i]);
			i++;
			profileH = atoi(argv[i]);
			break;
		}
		else if(argv[i][j] == 's')
		{
			in_sigWall = true;
			i++;
			sigWall = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'S'){
			stripes = true;
			i++;
			numberOfStripes = atof(argv[i]);
			break;
		}
		else if (argv[i][j] == 't'){
			thermostat = THERMOSTAT_ANDERSEN;
			i++;
			nuFactor = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'T')
		{
			in_temperature = true;
			i++;
			temp = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'u') LJunits = true;
		else if(argv[i][j] == 'w')
		{
			in_wallThick = true;
			i++;
			wallThick = atoi(argv[i]);
			break;
		}
		else if(argv[i][j] == 'W') //@todo: Tersoff model to be implemented!
		{
			in_wallModel = true;
			i++;
			if(!strcmp(argv[i], "LJ")) wall = WALL_CU_LJ;
			else if(!strcmp(argv[i],"Tersoff")) wall = WALL_TERSOFF;
			else if(!strcmp(argv[i], "MatAkaFit")) wall = WALL_MAT_AKA_FIT;
			else
			{
				std::cout << "Wall model" << argv[i] << "is not available. \n\n" << usage;
				return 4;
			}
			break;
		}
		else if(argv[i][j] == 'x'&& argv[i][j+1] == '1'&& argv[i][j+2] == '2')
		{
			in_xi12 = true;
			i++;
			xi12 = atof(argv[i]);
			break;
		}
		else if(argv[i][j] == 'x'&& argv[i][j+1] == '1'&& argv[i][j+2] == '3')
		{
			in_xi13 = true;
			i++;
			xi13 = atof(argv[i]);
			break;
		}
	} // end for(j...)
} // end for(i...)

// checking the input arguments for completeness and consistency
// (i) completeness => mandatory arguments: file name (i.e. prefix), number of fluid particles, temperature, wall thickness, xi_fluid_wall
if(in_prefix == false)
{
	std::cout << "No output prefix specified. \n\n" << usage;
	return 5;
}
if(in_N == false)
{
	std::cout << "Number of fluid particles not specified. \n\n" << usage;
	return 6;
}
if(in_temperature == false)
{
	std::cout << "No temeprature specified. \n\n" << usage;
	return 7;
}
if(in_wallThick == false)
{
	std::cout << "Wall thickness not specified. \n\n" << usage;
	return 8;
}
if(in_xi12 == false)
{
	std::cout << "Berthelot combining rule: xi_12 not specified. \n\n" << usage;
	return 9;
}
// (ii) consistency
if(in_wallModel == true)
{
	if(wall == WALL_TERSOFF)
	{
		std::cout << "The Tersoff potential is not implemented yet.\n\n" << usage;
		return 10;
	}
	else if(wall == WALL_MAT_AKA_FIT)
	{
		std::cout << "The model rendering TiO2 by fitting the Matsui+Akaogi potential is not implemented yet. \n\n" << usage;
		return 11;
	}
}
if(!stripes && in_xi13){
	std::cout << "\n\n*********************************************************************************************************\n";
	std::cout << "Surplus specification of the Berthelot xi_13:\nHas to be specified for the mixing rule of components 1 and 3 only \n"
			"if there is a stripes shaped wall! (or any other kind of three component system)\n\n" << usage;
	return 12;
}
else if(stripes && ! in_xi13){
	std::cout << "\n\n*********************************************************************************************************\n";
	std::cout << "No xi_13 specified: \nIf a stripes shaped wall is employed the Berthelot xi_13 must be specified.\n\n" << usage ;
	return 12;
}

/*

 //@todo: stripes shaped wall
if(stripes == true)
{
	std::cout << "A wall exhibiting different values of xi in a stripes shaped manner is not implemented yet. \n\n" << usage;
	return 12;
}
*/

//@todo: input in Lennard-Jones units
/*
if(LJunits == true)
{
	std::cout<<"Input in Lennard-Jones units not enabled yet. Atomic units to be used instead. \n\n" << usage;
	return 13;
}*/


// setting the default values
if(!in_alpha) alpha = DEFAULT_ALPHA;
if(!in_beta) beta = DEFAULT_BETA;
if(!in_gamma) gamm = DEFAULT_GAMMA;
if(!in_edgeProp) edgeProp = DEFAULT_EDGE;
if(!in_density) densFac = 1.0;
if(!in_eta12) eta12 = 1.0;
if(!in_xi13) xi13 = 0.0;
if(!in_fluid) fluid = FLUID_AR;
if(!in_outputTime) outTimeSteps = 500*1000;
if(!in_wallModel) wall = WALL_CU_LJ;
if(!in_initCanon) initCanon = 0;









//@brief: defining the cutoff radii of the different interactions/models
double  chargeCutoffRadius;
chargeCutoffRadius = 30;



prefixStr = prefix;
//cout << "prefixStr: " << prefixStr<< "\n";
PhaseSpaceWriter PSW(prefixStr, temp, densFac, N, fluid, wall, wallThick, xi12, xi13, eta12,  alpha, beta, gamm, edgeProp, stripes, numberOfStripes, LJShifted, LJunits);
PSW.write();
//@todo: the call of Configwriter is not done in a wrong way: the wall LJ cut off radius does not exist!!! the fluid r_c is used instead!
// the same holds for the overall cutoff radius!!!
Component fluidComp(fluid, LJunits);
double refTime;
refTime = fluidComp.gRefTime(LJunits);
double boxLengthY = PSW.gBoxLengthY();
if(!in_numProfileUnits){
	profilePhi = DEFAULT_PROFILE_PHI;
	profileR = DEFAULT_PROFILE_R;
	profileH = (int)10.0*boxLengthY / fluidComp.gSigma(0); // default profile spacing in y-direction: 0.5*sigma
}

// generating an instance of ConfigWriter and calling the write method
if(thermostat == THERMOSTAT_VELSCALE){
  std::cout << "Velocity scaling applied \n";
  ConfigWriter CfgWriter(prefix, wall, wallThick, refTime, profilePhi, profileR, profileH, outTimeSteps, initCanon, movie, fluidComp);
  CfgWriter.write();
}
else if(thermostat == THERMOSTAT_ANDERSEN){
  ConfigWriter CfgWriter(prefix, wall, wallThick, refTime, profilePhi, profileR, profileH, outTimeSteps, initCanon, movie, PSW, fluidComp, nuFactor);
  CfgWriter.write();
}

//cout << "\n**********************************\nConfig file written\n**********************************\n";
} // end main
