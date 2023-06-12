/*
 * ConfigWriter.cpp
 *
 *  Created on: 05.01.2012
 *      Author: becker
 */


#include"ConfigWriter.h"

const double DT = 0.030620; // corresponds to 1 fs


extern const std::string WALL_TERSOFF;
extern const std::string WALL_CU_LJ;
extern double LATTICE_CONST_WALL_LJTS;



// @brief: implementing the constructor and destructor, respectively
ConfigWriter::ConfigWriter(
		char* in_prefix, std::string in_wall, int in_wallLays, double in_refTime,
		unsigned in_profilePhi, unsigned in_profileR, unsigned in_profile_H,
		unsigned in_profileOutputTimesteps, unsigned initCanon, bool in_movie, Component& fluidComp
		)
{
	std::cout << "\n**********************************\nConfigwriter opened\n**********************************\n";
	thermostat = THERMOSTAT_VELSCALE;
	sPrefix(in_prefix);
	sTimestepLength(in_refTime);
	wallLays = in_wallLays;
	sigFluid = fluidComp.gSigma(0);
	//unknown: what kind of cutoff radius???
	sCutoffRadius(fluidComp.gRCutLJ());
	sLjCutoffRadius(fluidComp.gRCutLJ());
	// the only information needed in the config file:
	// what wall model applied => writing out the corresponding cutoff radius
	/*if(in_wall == WALL_TERSOFF)
	{
		wall = WALL_TERSOFF;
		sTersoffCutoffRadius(fluidComp.gRCutLJ());	// in case the wall is modeled by the Tersoff potential
	}
	//@brief: obsolete!
	else if (in_wall == WALL_CU_LJ)
	{
		wall = WALL_CU_LJ;
		sLjWallCutoffRadius(in_wallCutoffRadius);
	}
	// the error message following should not be necessary currently,
	//implemented for error detection and further models implemented
	//@todo: if no error reported during test runs => removing the else-branch
	else
	{
		std::cout << "wall model: only Lennard-Jones TS.\n";
		//return 51;
	}*/
	sProfile (in_profilePhi, in_profileR, in_profile_H);
	sProfileOutputTimesteps(in_profileOutputTimesteps);
	// values set as default, if required they can be passed by main.cpp, i.e. as a user input => additional constructor
	sInitCanonical(initCanon);					// temperature raise during equilibration applied? (ask M.Horsch) => short equilibration time needed otherwise vaporisation of the drop
	sInitStatistics(500000);
	sProfileRecordingTimesteps(2);
	sOutputResWriter(40);
	sOutputXyzWriter(500000);
	// if no movie is to be made, there's no need for the visitt-writer, by default: movie = false; then also no VisItt output is generated: if(movie) wOutputVisittWriter();
	_movie = in_movie;
	if(_movie){
		sOutputVisittWriter(500);
		sOutputMmspdWriter(500);
	}

	//	std::cout << "\n**********************************\nConstructor of Configwriter finished\n**********************************\n";
}

ConfigWriter:: ConfigWriter(	char* in_prefix, std::string in_wall, int in_wallLays, double in_refTime,
			unsigned in_profilePhi, unsigned in_profileR, unsigned in_profile_H,
			unsigned in_profileOutputTimesteps, unsigned initCanon, bool in_movie, PhaseSpaceWriter& psw, Component& fluidComp,
			double nuAndFac	){
  std::cout << "\n**********************************\nConfigwriter opened\n**********************************\n";
	thermostat = THERMOSTAT_ANDERSEN;
	sPrefix(in_prefix);
	sTimestepLength(in_refTime);
	wallLays = in_wallLays;
	sigFluid = fluidComp.gSigma(0);
	//unknown: what kind of cutoff radius???
	sCutoffRadius(fluidComp.gRCutLJ());
	sLjCutoffRadius(fluidComp.gRCutLJ());
	// the only information needed in the config file:
	// what wall model applied => writing out the corresponding cutoff radius
	/*if(in_wall == WALL_TERSOFF)
	{
		wall = WALL_TERSOFF;
		sTersoffCutoffRadius(fluidComp.gRCutLJ());	// in case the wall is modeled by the Tersoff potential
	}
	//@brief: obsolete!
	else if (in_wall == WALL_CU_LJ)
	{
		wall = WALL_CU_LJ;
		sLjWallCutoffRadius(in_wallCutoffRadius);
	}
	// the error message following should not be necessary currently,
	//implemented for error detection and further models implemented
	//@todo: if no error reported during test runs => removing the else-branch
	else
	{
		std::cout << "wall model: only Lennard-Jones TS.\n";
		//return 51;
	}*/
	sProfile (in_profilePhi, in_profileR, in_profile_H);
	sProfileOutputTimesteps(in_profileOutputTimesteps);
	// values set as default, if required they can be passed by main.cpp, i.e. as a user input => additional constructor
	sInitCanonical(initCanon);					// temperature raise during equilibration applied? (ask M.Horsch) => short equilibration time needed otherwise vaporisation of the drop
	sInitStatistics(500000);
	sProfileRecordingTimesteps(2);
	sOutputResWriter(40);
	sOutputXyzWriter(500000);
	// if no movie is to be made, there's no need for the visitt-writer, by default: movie = false; then also no VisItt output is generated: if(movie) wOutputVisittWriter();
	_movie = in_movie;
	if(_movie){
		sOutputVisittWriter(500);
		sOutputMmspdWriter(500);
	}

	double averageMassPerParticle = 1;// psw.gAverageMassPerParticle();
	double diffCoeffLJ = 0.05; // estimate of the self-diffusion coefficient of the LJ-Fluid
	nuAndersenSingle = psw.gTemperature()*timestepLength/ averageMassPerParticle/diffCoeffLJ;
	nuAndersen = nuAndersenSingle*pow(psw.gNTotal(), -2.0/3.0)*nuAndFac;
	std::cout << "nu Andersen: " << nuAndersen << "\n";

	//	std::cout << "\n**********************************\nConstructor of Configwriter finished\n**********************************\n";
}

ConfigWriter::~ConfigWriter()
{
	confStrm.close();
	//cout << "Config file completed.\n";
}

//@brief: superior writing method: handling the streams, calling the single write-methods
void ConfigWriter::write(){

//	std::cout << "\n**********************************\nwrite() method of Configwriter started\n**********************************\n";
	confFile << prefix << "_1R.cfg";					// building the file prefix.cfg
	confStrm.open(confFile.str().c_str(), std::ios::trunc);	// linking the file prefix.cfg with confStrm => confStrm writes in prefix.cfg

	std::cout << "\n**********************************\nWriting the config file \n**********************************\n\n";

	confStrm << "MDProjectConfig\n";
	wTimestepLength();
	wCutoffRadius();
	wLjCutoffRadius();
	confStrm <<"\n";
	// set the appropriate cutoff radius, the appropriate wall model, respectively
	//if(wall == WALL_CU_LJ)	wLjWallCutoffRadius();
	//else if (wall == WALL_TERSOFF) wTersoffCutoffRadius();
	wInitCanonical();
	wInitStatistics();
	wPhaseSpaceFileName();
        //confStrm << "parallelization DomainDecomposition \n";
	confStrm << "parallelization KDDecomposition 2000 3 \n";
	confStrm << "# for LinkedCells, the cellsInCutoffRadius has to be provided\n"
			 << "datastructure\tLinkedCells \t1\n";
	wOutputResWriter();
	wOutputXyzWriter();
	// only if a movie is to be made
	if(_movie){
	  wOutputMmspdWriter();
	}
	wProfile();
	confStrm << "yOffset\t" << (wallLays -0.5 + 0.05)*LATTICE_CONST_WALL_LJTS << "\n"; // +0.05*latticeConstant => little offset to avoid particles beeing placed outside the simulation box
	confStrm << "profileVirial\n";
	wProfileRecordingTimesteps();
	wProfileOutputTimesteps();
	confStrm <<"profiledComponent\t1\n";
	confStrm << "SessileDrop\n";
	wProfileOutputPrefix();
	confStrm << "AlignCentre\t25\t1\n" ;
	confStrm << "ComponentForYShift\t2 3\n";
	//confStrm << "nomomentum\t200\n";
	confStrm << "thermostat 2 0.05\n";
	if(thermostat == THERMOSTAT_ANDERSEN){
	confStrm << "thermostat 2 " << nuAndersen << "\n";
	}
	confStrm << "NumberOfFluidComponents\t1\n";
	confStrm << "WallFun_LJ_9_3 2 1.0741 1 100 0 100 0.05 1 1 1\n";
	confStrm << "Mirror 50\n";
}


//@brief: implemeting the set-method s...(arg)
	void	ConfigWriter::sPrefix(char* input)
	{
		prefix = input;
	}
	void 	ConfigWriter::sTimestepLength(double input)
	{
		timestepLength = DT/input;
	}
	void	ConfigWriter::sCutoffRadius(double input)
	{
		cutoffRadius = input;
	}
	void	ConfigWriter::sLjCutoffRadius(double input)
	{
		ljCutoffRadius = input;
	}
	void	ConfigWriter::sTersoffCutoffRadius(double input)
	{
		tersoffCutoffRadius = input;
	}
	void	ConfigWriter::sLjWallCutoffRadius(double input)
		{
			ljWallCutoffRadius = input;
		}
	void	ConfigWriter::sInitCanonical(unsigned input)
	{
		initCanonical = input;
	}
	void	ConfigWriter::sInitStatistics(unsigned input)
	{
		initStatistics = input;
	}
	void	ConfigWriter::sProfileOutputTimesteps(unsigned input)
	{
		profileOutputTimesteps = input;
	}
	void	ConfigWriter::sProfileRecordingTimesteps(unsigned input)
	{
		profileRecordingTimesteps = input;
	}
	void	ConfigWriter::sProfile (unsigned profilePhi, unsigned profileR, unsigned profileH)
	{
		profile[0] = profilePhi;
		profile[1] = profileR;
		profile[2] = profileH;
	}
	void 	ConfigWriter::sHProfile(unsigned profileH){
		profile[2] = profileH;
	}
	void 	ConfigWriter::sOutputResWriter(unsigned input)
			{
				outputResWriter = input;
			}
	void	ConfigWriter::sOutputXyzWriter(unsigned input)
	{
		outputXyzWriter = input;
	}
	void	ConfigWriter::sOutputVisittWriter(unsigned input)
	{
		outputVisittWriter = input;
	}
	void 	ConfigWriter::sOutputMmspdWriter(unsigned input)
	{
		outputMmspdWriter = input;
	}



//@brief: implementing the write-methods w...()
	char*	ConfigWriter::wPrefix()
		{
		return prefix;
		}
	void 	ConfigWriter::wTimestepLength()
		{
		confStrm << "timestepLength\t"<<timestepLength<<"\n";
		}
	void 	ConfigWriter::wCutoffRadius()
		{
		confStrm << "cutoffRadius\t"<< cutoffRadius<<"\n";
		}
	void 	ConfigWriter::wLjCutoffRadius()
		{
		confStrm<<"LJCutoffRadius\t"<< ljCutoffRadius<<"\n";
		}
	void	ConfigWriter::wTersoffCutoffRadius()
		{
		confStrm<<"tersoffCutoffRadius\t"<< tersoffCutoffRadius<<"\n";
		}
	void	ConfigWriter::wLjWallCutoffRadius()
			{
			confStrm<<"ljWallCutoffRadius\t"<< ljWallCutoffRadius<<"\n";
			}
	void 	ConfigWriter::wInitCanonical()
		{
		confStrm<<"initCanonical\t"<< initCanonical<<"\n";
		}
	void 	ConfigWriter::wInitStatistics()
		{
		confStrm<<"initStatistics\t"<< initStatistics<<"\n";
		}
	void 	ConfigWriter::wProfileOutputTimesteps()
		{
		confStrm<<"profileOutputTimesteps\t"<< profileOutputTimesteps<<"\n";
		}
	void 	ConfigWriter::wProfileRecordingTimesteps()
		{
		confStrm<<"profileRecordingTimesteps\t"<< profileRecordingTimesteps<<"\n";
		}
	void 	ConfigWriter::wProfile ()
		{
		confStrm<<"profile\t"<< profile[0]<<"\t"<<profile[1]<<"\t"<<profile[2]<<"\n";
		}
	void	ConfigWriter::wPhaseSpaceFileName()
		{
		confStrm<<"phaseSpaceFile\t"<< "OldStyle\t"<< prefix<<".inp\n";
		}
	void	ConfigWriter::wOutputResWriter()
		{
		confStrm<<"output\tResultWriter\t"<< outputResWriter<<"\t"<< wPrefix()<<"\n";
		}
	void	ConfigWriter::wOutputXyzWriter()
		{
		confStrm<<"output\tXyzWriter\t"<< outputXyzWriter<< "\t"<< wPrefix()<<".buxyz\n";
		}
	void	ConfigWriter::wOutputVisittWriter()
		{
		confStrm<<"output\tVisittWriter\t"<< outputVisittWriter<<"\t"<<wPrefix()<<"\n";
		}
	void	ConfigWriter::wOutputMmspdWriter()
		{
		confStrm<<"output\tMmspdWriter\t"<< outputMmspdWriter<<"\t"<<wPrefix()<<"\n";
		}
	void	ConfigWriter::wProfileOutputPrefix()
		{
		confStrm<<"profileOutputPrefix\t"<< wPrefix()<<"\n";
		}



/*
// implementing the merge(args)-functions
	char* merge(char* s1, char* s2, char* s3)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		return merged;
	}
	char* merge(char* s1, char* s2, char* s3, char*s4)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		merged = strcat(merged, s4);
	}
	char* merge(char* s1, char* s2, char* s3, char*s4, char* s5)
	{
		char* merged;
		merged = strcat(s1, s2);
		merged = strcat(merged, s3);
		merged = strcat(merged, s4);
		merged = strcat(merged, s5);
	}
*/

