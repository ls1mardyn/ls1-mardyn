/*
 * MPI_IOCheckpointWriter.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: andal
 */

#include "io/MPI_IOCheckpointWriter.h"

#include <sstream>
#include <string>
#include <iostream>

#include "Common.h"
#include "Domain.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "utils/mardyn_assert.h"

#include "molecules/Molecule.h"
#include "particleContainer/ParticleCell.h"
#ifdef ENABLE_MPI
#include "parallel/ParticleData.h"
#include "parallel/DomainDecompBase.h"
#include <mpi.h>
#endif
//#include <time.h>

#include "particleContainer/LinkedCells.h"
#include "particleContainer/Cell.h"

#include <cassert>

MPI_IOCheckpointWriter::MPI_IOCheckpointWriter(unsigned long writeFrequency,
		std::string outputPrefix, bool incremental) {
	_outputPrefix = outputPrefix;
	_writeFrequency = writeFrequency;
	_incremental = incremental;

	if (outputPrefix == "default") {
		_appendTimestamp = true;
	} else {
		_appendTimestamp = false;
	}
}

MPI_IOCheckpointWriter::~MPI_IOCheckpointWriter() {
	// TODO Auto-generated destructor stub
}

void MPI_IOCheckpointWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	Log::global_log->info() << "Write frequency: " << _writeFrequency << std::endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	Log::global_log->info() << "Output prefix: " << _outputPrefix << std::endl;

	int incremental = 1;
	xmlconfig.getNodeValue("incremental", incremental);
	_incremental = (incremental != 0);
	Log::global_log->info() << "Incremental numbers: " << _incremental << std::endl;

	int appendTimestamp = 0;
	xmlconfig.getNodeValue("appendTimestamp", appendTimestamp);
	if (appendTimestamp > 0) {
		_appendTimestamp = true;
	}
	Log::global_log->info() << "Append timestamp: " << _appendTimestamp << std::endl;
}

void MPI_IOCheckpointWriter::init(ParticleContainer *particleContainer,
                                  DomainDecompBase *domainDecomp, Domain *domain) {

}

void MPI_IOCheckpointWriter::endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp,
                                     Domain *domain,
                                     unsigned long simstep) {
#ifdef ENABLE_MPI
	if (simstep % _writeFrequency == 0) {
		//get the file name
		std::stringstream filenamestream;
		filenamestream << _outputPrefix;

		if (_incremental) {
			/* align file numbers with preceding '0's in the required range from 0 to _numberOfTimesteps. */
			unsigned long numTimesteps = _simulation.getNumTimesteps();
			int num_digits = (int) ceil(
					log(double(numTimesteps / _writeFrequency)) / log(10.));
			filenamestream << "-" << aligned_number(simstep / _writeFrequency,
					num_digits, '0');
		}
		if (_appendTimestamp) {
			filenamestream << "-" << gettimestring();
		}
		filenamestream << ".restart";
		domain->writeCheckpointHeader(filenamestream.str(), particleContainer,
				domainDecomp, _simulation.getSimulationTime());
		filenamestream << ".mpi";
		std::string filename = filenamestream.str();


		//some debug stuff to gather cell information from the LinkedCells Class
		/*
		int *boxCellDimension = lcContainer->boxWidthInNumCells();
		long realLocalNumCells = boxCellDimension[0]*boxCellDimension[1]*boxCellDimension[2];
		long realGlobalNumCells = 0;
		MPI_Allreduce(&realLocalNumCells, &realGlobalNumCells, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
		std::cout << "Rank: " << domainDecomp->getRank() << " realGlobalNumCells: " << realGlobalNumCells << std::endl;
		std::cout << "Rank: " << domainDecomp->getRank() << " BoundingBoxMin: " << domainDecomp->getBoundingBoxMin(0, domain) << "," << domainDecomp->getBoundingBoxMin(1, domain) << "," << domainDecomp->getBoundingBoxMin(2, domain) << " BoundingBoxMax: " << domainDecomp->getBoundingBoxMax(0, domain) << "," << domainDecomp->getBoundingBoxMax(1, domain) << "," << domainDecomp->getBoundingBoxMax(2, domain) << std::endl;
		domainDecomp->getBoundingBoxMin(0, domain);
		*/


		//cell length for the cell structure in the output file
		//here each cell has the same radius in x-,y- and z-direction
		double *cellLength = particleContainer->getCellLength();

		//compute lengths and sizes of the domain
		int lengthInCells[3];
		for (unsigned short i = 0; i < 3; i++) {
			lengthInCells[i] = domain->getGlobalLength(i) / cellLength[i];
		}

		long globalNumCells = 1;
		size_t localNumCells = 0;

		for (unsigned short i = 0; i < 3; i++) {
			mardyn_assert(lengthInCells[i] > 0);
//			mardyn_assert((lengthInCells[i] - (domain->getGlobalLength(i) / cellLength[i])) == 0);

			if (lengthInCells[i] != 0) {
				globalNumCells *= lengthInCells[i];
			}
		}

		//compute localNumParticles, isCellOfProcess and localNumCells
		std::vector<int> globalNumParticlesPerCell(globalNumCells);
		std::vector<int> localNumParticlesPerCell(globalNumCells);
		std::vector<bool> isCellOfProcess(globalNumCells); //true if the cell is located in this process

		for (int i = 0; i < globalNumCells; i++) {
			globalNumParticlesPerCell[i] = 0;
			localNumParticlesPerCell[i] = 0;
			isCellOfProcess[i] = false;
		}


		for (auto tempMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMolecule.isValid(); ++tempMolecule) {
			int cellIndex[3];
			for (unsigned short i = 0; i < 3; i++) {
				cellIndex[i]
						= (int) floor(tempMolecule->r(i) / cellLength[i]);
			}
			unsigned long index = (cellIndex[2] * lengthInCells[1]
					+ cellIndex[1]) * lengthInCells[0] + cellIndex[0];
			if (localNumParticlesPerCell[index] == 0) {
				localNumCells++;
				isCellOfProcess[index] = true;
			}
			localNumParticlesPerCell[index]++;

		}

		/*
		//Timer:
		timeval timer1, timer2;
		double timeDiff;

		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer1, NULL);
		}

		*/
		//every process has to know how many particles each cell in the linked cell has,
		//as every process has to know the same header data
		MPI_Allreduce(localNumParticlesPerCell.data(), globalNumParticlesPerCell.data(),
				globalNumCells, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		for (int i = 0; i < globalNumCells; i++) {
			if (localNumParticlesPerCell[i] > 0) {
				if (globalNumParticlesPerCell[i] != localNumParticlesPerCell[i]) {
					Log::global_log->info() << "Allreduce ist fehlgeschlagen"
							<< std::endl;
				}
			}
		}

		/*
		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer2, NULL);
			timeDiff = timer2.tv_sec - timer1.tv_sec + (timer2.tv_usec
					- timer1.tv_usec) / 1.E6;
			Log::global_log->info() << "Das MPI-IO Allreduce hat " << timeDiff
					<< " Sekunden benötigt" << std::endl;
		}
		*/

		//Parallel IO
		int ret, size;
		MPI_Offset offset = 0;
		MPI_Status status;

		const char* fileName = filename.c_str();

		MPI_File fh;
		MPI_Info info;

		MPI_Info_create(&info);

		/*
		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer1, NULL);
		}
		*/

		ret = MPI_File_open(MPI_COMM_WORLD, fileName,
				MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		//setGlobalNumCells, cellLength, NumParticles
		long numCellsAndMolecules[2];
		numCellsAndMolecules[0] = globalNumCells;
		numCellsAndMolecules[1] = domain->getglobalNumMolecules(true, particleContainer, domainDecomp);

		if (domainDecomp->getRank() == 0) {
			ret = MPI_File_seek(fh, 0, MPI_SEEK_SET);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}

			ret
					= MPI_File_write(fh, numCellsAndMolecules, 2, MPI_LONG,
							&status);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}

			ret = MPI_Type_size(MPI_LONG, &size);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}
			offset = 2 * size;
			ret = MPI_File_seek(fh, offset, MPI_SEEK_SET);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}

			ret = MPI_File_write(fh, cellLength, 3, MPI_DOUBLE, &status);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}

			ret = MPI_Type_size(MPI_DOUBLE, &size);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}
			offset += (3 * size);
			ret = MPI_File_seek(fh, offset, MPI_SEEK_SET);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}
			ret = MPI_File_write(fh, globalNumParticlesPerCell.data(), globalNumCells,
					MPI_INT, &status);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}

			ret = MPI_Type_size(MPI_INT, &size);
			if (ret != MPI_SUCCESS) {
				handle_error(ret);
			}
			offset += (globalNumCells * size);
		}

		ret = MPI_Bcast(&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		/*
		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer2, NULL);
			timeDiff = timer2.tv_sec - timer1.tv_sec + (timer2.tv_usec
					- timer1.tv_usec) / 1.E6;
			Log::global_log->info() << "Das Schreiben des MPI-IO Headers hat "
					<< timeDiff << " Sekunden benötigt" << std::endl;
		}

		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer1, NULL);
		}

		*/

		//Writing:
		MPI_Datatype mpiParticleData;
		ParticleData::getMPIType(mpiParticleData);

		//First dimension represents the cell
		//Second dimension contains the molecules of this cell
		std::vector<ParticleData *> writeArray;
		writeArray.resize(localNumCells);

		//define the sizes of the second dimension and a mapping from a global cell index
		//to a local cell index(cell contained in the local process)
		//
		//entry of the array is -1 if the cell is not located at this process
		std::vector<int> globalToLocalCell(globalNumCells);
		int j = 0;
		for (size_t i = 0; i < localNumCells; i++) {
			while (isCellOfProcess[j] == false) {
				globalToLocalCell[j] = -1;
				j++;
			}
			//check, if the I/O cells are each only on one process
			mardyn_assert(localNumParticlesPerCell[j] == globalNumParticlesPerCell[j]);

			globalToLocalCell[j] = i;
			if (globalNumParticlesPerCell[j] > 0) {
				writeArray[i] = new ParticleData[localNumParticlesPerCell[j]];
			}
			j++;
		}

		//counter for each cell
		std::vector<int> cellCounter;
		cellCounter.resize(localNumCells);
		for (size_t i = 0; i < localNumCells; i++) {
			cellCounter[i] = 0;
		}

		//filling of the writeArray
		for (auto tempMolecule = particleContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); tempMolecule.isValid(); ++tempMolecule) {

			int cellIndex[3];
			for (unsigned short i = 0; i < 3; i++) {
				cellIndex[i] = floor(tempMolecule->r(i) / cellLength[i]);
			}
			int index = (cellIndex[2] * lengthInCells[1] + cellIndex[1])
					* lengthInCells[0] + cellIndex[0];
			ParticleData::MoleculeToParticleData(
					writeArray[globalToLocalCell[index]][cellCounter[globalToLocalCell[index]]],
					*tempMolecule);
			cellCounter[globalToLocalCell[index]]++;
		}

		for (int i = 0; i < globalNumCells; i++) {
			if (isCellOfProcess[i] == true) {
				mardyn_assert(cellCounter[globalToLocalCell[i]] == globalNumParticlesPerCell[i]);
			}
		}

		/*
		if (domainDecomp->getRank() == 0) {
			gettimeofday(&timer2, NULL);
			timeDiff = timer2.tv_sec - timer1.tv_sec + (timer2.tv_usec
					- timer1.tv_usec) / 1.E6;
			Log::global_log->info() << "Das Belegen des MPI-IO Schreibearrays hat "
					<< timeDiff << " Sekunden benötigt" << std::endl;
		}

		//so we know that the size of each local cell is in cellCounter


		gettimeofday(&timer1, NULL);

		*/


		ret = MPI_Type_size(mpiParticleData, &size);
		if (ret != MPI_SUCCESS) {
			handle_error(ret);
		}

		for (int i = 0; i < globalNumCells; i++) {
			if (globalNumParticlesPerCell[i] == 0) {
				continue;
			}
			if (isCellOfProcess[i] == true) {

				ret = MPI_File_seek(fh, offset, MPI_SEEK_SET);

				ret = MPI_File_write(fh, writeArray[globalToLocalCell[i]],
						globalNumParticlesPerCell[i], mpiParticleData, &status);
			}

			offset += (globalNumParticlesPerCell[i] * size);
		}

		/*
		gettimeofday(&timer2, NULL);
		timeDiff = timer2.tv_sec - timer1.tv_sec
				+ (timer2.tv_usec - timer1.tv_usec) / 1.E6;
		double timeDiffGlobal = 0;
		MPI_Reduce(&timeDiff, &timeDiffGlobal, 1, MPI_DOUBLE, MPI_MAX, 0,
					MPI_COMM_WORLD);

		*/

		/*
		if (domainDecomp->getRank() == 0) {

			Log::global_log->info() << "Das Lesen der Zellen hat " << timeDiffGlobal
					<< " Sekunden benötigt" << std::endl;
		}

		*/

		//delete tempMolecule;
		for (size_t i = 0; i < localNumCells; i++) {
			delete[] writeArray[i];
		}
		MPI_Type_free(&mpiParticleData);

		MPI_File_close(&fh);
	}
#endif
}

void MPI_IOCheckpointWriter::finish(ParticleContainer *particleContainer,
									DomainDecompBase *domainDecomp, Domain *domain) {

}

void MPI_IOCheckpointWriter::handle_error(int i) {
#ifdef ENABLE_MPI
	char error_string[BUFSIZ];
	int length_of_error_string;

	MPI_Error_string(i, error_string, &length_of_error_string);

	Log::global_log->error() << "Writing of file was not successfull " << " , " << i
			<< " , " << error_string << std::endl;
	MARDYN_EXIT(1);
#endif
}
