/*
 * benchmark.cpp
 *
 *  Created on: Aug 17, 2011
 *      Author: andreas
 */
#include <stdio.h>

#include "benchmark.h"

void SimulationStats::writeFrameStats( const std::string &frameFile ) {
	FILE *file;

	file = fopen( frameFile.c_str(), "wt" );

	fprintf( file, "potential,virial" );
#ifdef COMPARE_TO_CPU
	fprintf( file, ",forceRMSError,torqueRMSError" );
#endif
	fprintf( file, "\n" );

	for( int i = 0 ; i < potentials.getCount() ; i++ ) {
		fprintf( file, "%.18e,%.18e", potentials[i], virials[i] );
#if defined(COMPARE_TO_CPU) && !defined(NO_CUDA)
		fprintf( file, ",%.18e,%.18e", forceRMSError[i], torqueRMSError[i] );
#endif
		fprintf( file, "\n" );
	}

	fclose( file );
}

void SimulationStats::writeRunStats( const std::string &buildFile ) {
	// write header?
	FILE *file;
	if( (file = fopen( buildFile.c_str(), "rt" )) != NULL ) {
		fclose(file);
		file = fopen( buildFile.c_str(), "at" );
	}
	else {
		file = fopen( buildFile.c_str(), "wt" );
		fprintf( file, "timeSteps,moleculeCount,cutOffRadius,totalTime" );
#ifndef NO_CUDA
		fprintf( file, ",numWarps,maxRegisterCount,"
				"warpBlockCellProcessor,doubleMode,"
				"maxNumComponents,maxNumLJCenters,maxNumCharges,maxNumDipoles,"
				"CUDA_totalTime,CUDA_preTime,CUDA_postTime,CUDA_singleTime,CUDA_pairTime,CUDA_processingTime" );
#endif
		fprintf( file, ",name\n" );
	}


	fprintf( file, "%i,%i,%e,%e", timeSteps, moleculeCount, cutOffRadius, (double) totalTime );
#ifndef NO_CUDA
	fprintf( file, ",%i,%i,"
			"%i,%i,%i,%i,%i,"
			"%i,%i,%i,%i,"
			"%e,%e,%e,%e,%e,%e",
			NUM_WARPS, MAX_REGISTER_COUNT,
#ifdef CUDA_WARP_BLOCK_CELL_PROCESSOR
			1,
#else
			0,
#endif
#ifdef CUDA_DOUBLE_MODE
			1,
#else
			0,
#endif
			MAX_NUM_COMPONENTS, MAX_NUM_LJCENTERS, MAX_NUM_CHARGES, MAX_NUM_DIPOLES,
			(double) CUDA_frameTime, (double) CUDA_preTime, (double) CUDA_postTime, (double) CUDA_singleTime, (double) CUDA_pairTime, (double) CUDA_processingTime
			);
#endif
	fprintf( file, ",%s\n", name.c_str() );

	fclose( file );
}

void SimulationStats::writeCellStats(const std::string &domainInfoFile, const std::vector<Cell> &cells, const std::string &domainInfo ) {
	FILE *file = fopen( domainInfoFile.c_str(), "at" );

	fprintf( file, "\"%s\"", domainInfo.c_str() );

	for( int i = 0 ; i < cells.size() ; i++ ) {
		fprintf( file, ",%i", cells[i].getMoleculeCount() );
	}

	fprintf( file, "\n" );
	fclose( file );
}

SimulationStats simulationStats;
