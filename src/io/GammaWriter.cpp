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

#include "io/GammaWriter.h"

#include "particleContainer/ParticleContainer.h"
#include "parallel/DomainDecompBase.h"
#include "Domain.h"

#include <ctime>
using std::ios_base;

using namespace std;

GammaWriter::GammaWriter(unsigned long writeFrequency, string outputPrefix)
: _writeFrequency(writeFrequency),
	_outputPrefix(outputPrefix)
{ }

GammaWriter::~GammaWriter(){}

void GammaWriter::initOutput(ParticleContainer* particleContainer,
			      DomainDecompBase* domainDecomp, Domain* domain){
	 
	// initialize result file
	string resultfile(_outputPrefix+".gamma");
	_gammaStream.precision(6);
	time_t now;
	time(&now);
	if(domainDecomp->getRank()==0){
		_gammaStream.open(resultfile.c_str());
		_gammaStream << "# mardyn MD simulation starting at " << ctime(&now) << endl;
		_gammaStream << "#\tgamma" << endl;
	}
}

void GammaWriter::doOutput( ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain,
			     unsigned long simstep, std::list<ChemicalPotential>* lmu ) 
{
	domain->calculateGamma(particleContainer,domainDecomp);
	if((domainDecomp->getRank() == 0) && (simstep % _writeFrequency == 0)){
		_gammaStream << simstep << "\t"; 
		for (unsigned i=0; i<domain->getNumberOfComponents(); i++){
			_gammaStream << domain->getGamma(i)/_writeFrequency << "\t";
		}
		_gammaStream << endl;
		domain->resetGamma();
	}
}

void GammaWriter::finishOutput(ParticleContainer* particleContainer,
				DomainDecompBase* domainDecomp, Domain* domain){
	_gammaStream.close();
}
