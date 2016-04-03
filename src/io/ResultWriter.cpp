/***************************************************************************
 *   Copyright (C) 2010 by Martin Bernreuther <bernreuther@hlrs.de> et al. *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "io/ResultWriter.h"

#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "Simulation.h"
#include "utils/Logger.h"
#include "ensemble/CavityEnsemble.h"

using Log::global_log;
using namespace std;


ResultWriter::ResultWriter(unsigned long writeFrequency, string outputPrefix)
: _writeFrequency(writeFrequency),
	_outputPrefix(outputPrefix)
{
	size_t ACC_STEPS = 1000;
	_U_pot_acc = new Accumulator<double>(ACC_STEPS);
	_p_acc = new Accumulator<double>(ACC_STEPS);
        
        _v1x_acc = new Accumulator<double>(ACC_STEPS);
        _v1y_acc = new Accumulator<double>(ACC_STEPS);
        _v1z_acc = new Accumulator<double>(ACC_STEPS);
        _b_v2ll_acc = new Accumulator<double>(ACC_STEPS);
        _b_v2lm_acc = new Accumulator<double>(ACC_STEPS);
        _aa_vv11dia_acc = new Accumulator<double>(ACC_STEPS);
        _aa_vv11off_acc = new Accumulator<double>(ACC_STEPS);
}

ResultWriter::~ResultWriter(){}

void ResultWriter::readXML(XMLfileUnits& xmlconfig) {
	_writeFrequency = 1;
	xmlconfig.getNodeValue("writefrequency", _writeFrequency);
	global_log->info() << "Write frequency: " << _writeFrequency << endl;

	_outputPrefix = "mardyn";
	xmlconfig.getNodeValue("outputprefix", _outputPrefix);
	global_log->info() << "Output prefix: " << _outputPrefix << endl;

	size_t acc_steps = 1000;
	xmlconfig.getNodeValue("accumulation_steps", acc_steps);
	_U_pot_acc = new Accumulator<double>(acc_steps);
	_p_acc = new Accumulator<double>(acc_steps);
	global_log->info() << "Accumulation steps: " << acc_steps << endl;
}

void ResultWriter::initOutput(ParticleContainer* particleContainer,
			      DomainDecompBase* domainDecomp, Domain* domain){
	 
	// initialize result file
	string resultfile(_outputPrefix+".res");
	time_t now;
	time(&now);
	if(domainDecomp->getRank()==0){
		_resultStream.open(resultfile.c_str());
		_resultStream << "# ls1 MarDyn simulation started at " << ctime(&now) << endl;
		_resultStream << "#step t\t\tU_pot <U_pot>\tp <p>\t\tT\tbeta_trans beta_rot\tc_v\t\tN\t(N_cav*)";
                
                double thrfactor = pow((double)MAX_THRESHOLD/MIN_THRESHOLD, 1.0/(NUM_THRESHOLD - 1.0));
                double tthreshold = MIN_THRESHOLD;
                for(unsigned i=0; i < NUM_THRESHOLD; i++)
                {
                   _resultStream << "\t(N_cav >= " << (unsigned)round(tthreshold) << ")*";
                   tthreshold *= thrfactor;
                }
                _resultStream << "\tmax(i_cav)\t\t<v1x> <v1y> <v1z>\t<b_v1> <b_v2ll> <b_v2lm> <b>\t<a^2_v1^2_dia> <a^2_v1^2_off> <-a^2/2T>\t<b> - <a^2/2T>\n";
	}
}

void ResultWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
	unsigned long simstep, list<ChemicalPotential>* lmu, map<unsigned, CavityEnsemble>* mcav )
{
	_U_pot_acc->addEntry(domain->getAverageGlobalUpot());
	_p_acc->addEntry(domain->getGlobalPressure());
        
        double T = domain->getGlobalCurrentTemperature();
        double v1x = domain->getGlobalVirial(0);
        double v1y = domain->getGlobalVirial(1);
        double v1z = domain->getGlobalVirial(2);
        
        _v1x_acc->addEntry(v1x);
        _v1y_acc->addEntry(v1y);
        _v1z_acc->addEntry(v1z);
        _b_v2ll_acc->addEntry(0.5 * domain->getGlobalVirialIILL());
        _b_v2lm_acc->addEntry(-0.5 * domain->getGlobalVirialIILM());
        _aa_vv11dia_acc->addEntry(-0.25 * (v1x*v1x + v1y*v1y + v1z*v1z) / T);
        _aa_vv11off_acc->addEntry(0.25 * (v1x*v1y + v1x*v1z + v1y*v1z) / T);
        
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_resultStream << simstep << " " << _simulation.getSimulationTime()
		              << "\t\t" << domain->getAverageGlobalUpot() << " " << _U_pot_acc->getAverage()
					  << "\t" << domain->getGlobalPressure() << " " << _p_acc->getAverage()
		              << "\t\t" << T << "\t" << domain->getGlobalBetaTrans() << " " << domain->getGlobalBetaRot()
		              << "\t" << domain->cv() << "\t\t" << domain->getglobalNumMolecules();
                 
                map<unsigned, CavityEnsemble>::iterator ceit;
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                {
                   _resultStream << "\t" << ceit->second.numCavities();
                }
                
                for(unsigned i=0; i < NUM_THRESHOLD; i++)
                {
                   for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                   {
                      _resultStream << "\t" << ceit->second.getThresholdPopulation(i);
                   }
                }
                unsigned largest_cavity = 0;
                for(ceit = mcav->begin(); ceit != mcav->end(); ceit++)
                {
                   unsigned tmax = ceit->second.getLargestCavity();
                   if(tmax > largest_cavity) largest_cavity = tmax;
                }
                _resultStream << "\t" << largest_cavity << "\t\t";
                
                double v1x_avg = _v1x_acc->getAverage();
                double v1y_avg = _v1y_acc->getAverage();
                double v1z_avg = _v1z_acc->getAverage();
                double b_v1_avg = -0.5 * (v1x_avg + v1y_avg + v1z_avg);
                double b_v2ll_avg = _b_v2ll_acc->getAverage();
                double b_v2lm_avg = _b_v2lm_acc->getAverage();
                double b_avg = b_v1_avg + b_v2ll_avg + b_v2lm_avg;
                double aa_vv11dia_avg = _aa_vv11dia_acc->getAverage();
                double aa_vv11off_avg = _aa_vv11off_acc->getAverage();
                double minaa_inv2T_avg = aa_vv11dia_avg + aa_vv11off_avg;
                double b_minaa_inv2T_avg = b_avg + minaa_inv2T_avg;
                
                _resultStream << v1x_avg << " " << v1y_avg << " "
                              << v1z_avg << "\t" << b_v1_avg << " "
                              << b_v2ll_avg << " " << b_v2lm_avg
                              << " " << b_avg << "\t"
                              << aa_vv11dia_avg << " "
                              << aa_vv11off_avg << " "
                              << minaa_inv2T_avg << "\t" 
                             << b_minaa_inv2T_avg << "\n";
	}
}

void ResultWriter::finishOutput(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain){
	time_t now;
	time(&now);
	_resultStream << "# ls1 MarDyn simulation finished at " << ctime(&now) << endl;
	_resultStream.close();
}
