#include "io/XyzWriter.h"

#include <fstream>
#include <sstream>
#include <map>
#include <sys/stat.h>

#include "Common.h"
#include "Domain.h"
#include "molecules/Molecule.h"
#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "utils/Logger.h"
#include "Simulation.h"

#include "io/VisitWriter.h"

using Log::global_log;
using namespace std;

XyzWriter::XyzWriter(unsigned long writeFrequency, string outputPrefix, bool incremental) {
	_outputPrefix= outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	}
	else {
		_appendTimestamp = false;
	}
}

XyzWriter::~XyzWriter(){}

void XyzWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	global_log->info() << "Incremental numbers: " << _incremental << endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if(appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	global_log->info() << "Append timestamp: " << _appendTimestamp << endl;
}

void XyzWriter::initOutput(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain) {
	// Check if file path already exists; important for independent use of xyz-writer and profile-writer
	FILE *f;
	FILE *g;
	FILE *h;
	f=fopen("./Results/Profile","r");
	g=fopen("./Results/ParticleData","r");
	h=fopen("./Results/Cavity","r");
	if(f==0)
	  mkdir("./Results/Profile", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(g==0)
	  mkdir("./Results/ParticleData", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(h==0)
	  mkdir("./Results/Cavity", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void XyzWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain, unsigned long simstep, list<ChemicalPotential>* lmu, std::map<unsigned, CavityEnsemble>* mcav   ) {
	if( simstep % _writeFrequency == 0) {
                map<unsigned, CavityEnsemble>::iterator ceit;
		stringstream filenamestream;
                map<unsigned, stringstream*> cav_filenamestream;
		stringstream vtk;
		filenamestream << "./Results/Profile/";
		filenamestream << _outputPrefix;
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                {
                   cav_filenamestream[ceit->first] = new stringstream;
		   *cav_filenamestream[ceit->first] << "./Results/Cavity/";
                   *cav_filenamestream[ceit->first] << _outputPrefix << "-c" << ceit->first;
                }
		vtk << "./Results/ParticleData/";
		vtk << _outputPrefix;
		
		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			filenamestream << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
                        for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                           *cav_filenamestream[ceit->first] << "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
			if(_simulation.isRecordingSlabProfile())
			  vtk <<  "-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			filenamestream << "-" << gettimestring();
			if(_simulation.isRecordingSlabProfile())
			  vtk << "-" << gettimestring();
		}
		filenamestream << ".xyz";
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                   *cav_filenamestream[ceit->first] << ".cav.xyz";

		vtk << ".vtk";
		
		
		int ownRank = domainDecomp->getRank();
		if( ownRank == 0 ) {
			ofstream xyzfilestream( filenamestream.str(). c_str() );
			xyzfilestream << domain->getglobalNumMolecules() << endl;
			xyzfilestream << "comment line" << endl;
			xyzfilestream.close();
                        
                        for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                        {
                           ofstream cavfilestream( cav_filenamestream[ceit->first]->str().c_str() );
                           cavfilestream << ceit->second.numCavities() << endl;
                           cavfilestream << "comment line" << endl;
                           cavfilestream.close();
			   
                        }

		}
		
		//NEW: Molecule-ID is written in the XYZ-File; afterwards Rank == 0 sorts the Molecules in each Timestep according to their ID
		for( int process = 0; process < domainDecomp->getNumProcs(); process++ ){
			domainDecomp->barrier();
			if( ownRank == process ){
			  
				ofstream xyzfilestream( filenamestream.str().c_str(), ios::app );
				ofstream particleDataFilestream( vtk.str().c_str(), ios::app );
				Molecule* tempMol;
				double temp = 0.0;
				double v[3];
				for(int d = 0; d < 3; d++)
				  v[d] = 0.0;
				for( tempMol = particleContainer->begin(); tempMol != particleContainer->end(); tempMol = particleContainer->next()){
					if( tempMol->componentid() == 0) { xyzfilestream <<  "Ar ";}
					else if( tempMol->componentid() == 1 ) { xyzfilestream <<  "Xe ";}
					else if( tempMol->componentid() == 2 ) { xyzfilestream <<  "C ";}
					else if( tempMol->componentid() == 3 ) { xyzfilestream <<  "O ";}
					else { xyzfilestream << "H ";}
					xyzfilestream << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << "\t" << tempMol->id() << endl;
					 temp = tempMol->getAveragedTemperature()/tempMol->getAverageCount();
					 if(_simulation.isRecordingSlabProfile()){
					  if (temp != temp)
					    temp = 0.0;
					  for(int d = 0; d < 3; d++){
					    v[d] = tempMol->getAveragedVelocity(d)/tempMol->getAverageCount();
					    if(v[d] != v[d])
					      v[d] = 0.0;
					  }
							
					  particleDataFilestream << tempMol->r(0) << "\t" << tempMol->r(1) << "\t" << tempMol->r(2) << "\t"
					  << temp << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << tempMol->componentid() << "\t" << tempMol->id() << endl;

					  tempMol->setAveragedVelocity(0,0);
					  tempMol->setAveragedVelocity(1,0);
					  tempMol->setAveragedVelocity(2,0);
					  tempMol->setAverageCount(0);
					  tempMol->setAveragedTemperature(0);
					}
				}
				xyzfilestream.close();
				particleDataFilestream.close();
				
                                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                                {
                                   ofstream cavfilestream( cav_filenamestream[ceit->first]->str().c_str(), ios::app );
                                   
                                   map<unsigned long, Molecule*> tcav = ceit->second.activeParticleContainer();
                                   map<unsigned long, Molecule*>::iterator tcit;
                                   for(tcit = tcav.begin(); tcit != tcav.end(); tcit++)
                                   {
                                      if( ceit->first == 0 ) { cavfilestream << "C ";}
                                      else if( ceit->first == 1 ) { cavfilestream << "N ";}
                                      else if( ceit->first == 2 ) { cavfilestream << "O ";}
                                      else if( ceit->first == 3 ) { cavfilestream << "F ";}
                                      else { cavfilestream << "Ne "; }
                                      cavfilestream << tcit->second->r(0) << "\t" << tcit->second->r(1) << "\t" << tcit->second->r(2) << "\n";
                                   }
                                
                                   cavfilestream.close();
				}
			}
		}
		//NEW: Sorting of Molecules in the same Order as in Timestep t=0 for a continously Post-Processing in VMD 
		stringstream arrangedFileName;
		arrangedFileName << "./Results/Profile/";
		arrangedFileName << _outputPrefix;
		
		//NEW: Initilisation of a new File
		if(_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil( log( double( numTimesteps / _writeFrequency ) ) / log(10.) );
			arrangedFileName << "_arranged-" << aligned_number( simstep / _writeFrequency, num_digits, '0' );
		}
		if(_appendTimestamp) {
			arrangedFileName << "-" << gettimestring();
		}
		arrangedFileName << ".xyz";
		// NEW: old XYZ-File is recorded to String-Array fileOutput[]; subsequently, Content of fileOutput[] is transfered to the new XYZ-File (arranged by the Molecule-ID)
		if( ownRank == 0 ) {
			std::fstream oldXyzFilestream;
			std::fstream oldVTK;
			oldXyzFilestream.open(filenamestream.str(). c_str() );
			oldVTK.open(vtk.str().c_str());
			string token;
			unsigned long count = 0;
			unsigned long numberOfMolecules;
			string comment;
			string element;
			double xPos, yPos, zPos;
			unsigned long molID;
			double tempVTK;
			double vTotal[3];
			int compID;
			
			std::map<int,std::string> fileOutput;
			unsigned long maxMolID = 0;
			while(!oldXyzFilestream.eof( )) {
			  
			  token.clear();
			  if(count == 0){
			    oldXyzFilestream >> numberOfMolecules;
			    count++;
			  }
			  if(count == 1){
			    oldXyzFilestream >> comment >> token;
			    count++;
			  }
			  token.clear();
			  oldXyzFilestream >> token;
			  //NEW: Cache for XYZ-Data of each Molecule
			  if(token == "Ar" || token == "Xe" || token == "C" || token == "O" || token == "H"){
			      stringstream cache;
			      string cache2;
			      oldXyzFilestream >> xPos >> yPos >> zPos >> molID;
			      cache << token << "\t" << xPos << "\t" << yPos << "\t" << zPos << "\n";
			      cache2 = cache.str();
			      fileOutput[molID] = cache2;
			      if(molID > maxMolID)
				maxMolID = molID;
			  }
			}
			
			
			oldXyzFilestream.close();
			
			//NEW: Output from fileOutput[] to new XYZ-File
			// The number of Molecules may differ from the highest IDs from the molecules
			for(unsigned long i = 0; i <= maxMolID; i++ ){
			 	ofstream newXyzFilestream( arrangedFileName.str().c_str(), ios::app );
				if(i == 0){
				  newXyzFilestream << numberOfMolecules << "\n" << comment << "\n";
				}
				// Demands if there is a content in map<i,string> at position "i"
				if(fileOutput.count(i)>0){
				  newXyzFilestream << fileOutput[i];
				  newXyzFilestream.close();
				}
			}
			//NEW: Deletes the old XYZ-File
			remove(filenamestream.str().c_str());
			fileOutput.clear();
			
			// VTK Molecule Data for VMD
			if(_simulation.isRecordingSlabProfile()){
			  float coordinates[3*domain->getglobalNumMolecules()];
			  int nvars = 2;
			  int vardims[] = {1, 3/*, 1*/};
			  const char *varnames[] = {"temperature", "velocity"/*, "component"*/};
			  float temperature[domain->getglobalNumMolecules()];
// 			  int componentID[domain->getglobalNumMolecules()];
			  float velocity[domain->getglobalNumMolecules()][3];
			  float *vars[] = {(float *)temperature, (float *)velocity/*, (int *)componentID*/};
			  
			
			  while(!oldVTK.eof()) {
			    stringstream cache;
			    string cache2;
			    oldVTK >> xPos >> yPos >> zPos >> tempVTK >> vTotal[0] >> vTotal[1] >> vTotal[2] >> compID >> molID;
			    cache << xPos << "\t" << yPos << "\t" << zPos << "\t" << tempVTK << "\t" << vTotal[0] << "\t" << vTotal[1]  << "\t" << vTotal[2] << "\t" << compID << "\t" << molID << "\n";
			    cache2 = cache.str();
			    fileOutput[molID] = cache2;
			  }
			  unsigned long countVTK = 1;
			  for (std::map<int,std::string>::iterator it=fileOutput.begin(); it!=fileOutput.end(); ++it){
				  stringstream cache3;
				  
				  cache3 << it->second;
				  cache3 >> xPos >> yPos >> zPos >> tempVTK >> vTotal[0] >> vTotal[1] >> vTotal[2] >> compID >> molID;
				  
				  coordinates[-3+3*countVTK] = xPos;
				  coordinates[-3+3*countVTK+1] = yPos;
				  coordinates[-3+3*countVTK+2] = zPos;
// 			    	  componentID[molID-1] = compID;
				  temperature[countVTK-1] = tempVTK;
				  velocity[countVTK-1][0] = vTotal[0];
				  velocity[countVTK-1][1] = vTotal[1];
				  velocity[countVTK-1][2] = vTotal[2];
				  countVTK++;
			  }
			
			  //NEW: Deletes the old XYZ-File
			  fileOutput.clear();
			
			  oldVTK.close();
			  remove(vtk.str().c_str());
			
			  // Output to VTK
			  write_point_mesh(vtk.str().c_str(), 0, numberOfMolecules, coordinates, nvars, vardims, varnames, vars);
			}
		}
	}
}

void XyzWriter::finishOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain ) {}
