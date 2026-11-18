/**
 * @file DomainDecompAnalyser.h
 * @author Georg von Bismarck
 * @date 18.11.2026
 */

#pragma once
#include <iosfwd>
#include <string>
#include <vector>

#ifdef ENABLE_MPI
	#include "parallel/DomainBox.h"
#endif
#include "plugins/PluginBase.h"


/**
 * Writes per-process domain decomposition data to a *.domainDecompAnalysis file.
 *
 * The output files can be read into pandas DataFrames with tools/DomainDecompAnalyser2Pandas.py.
 * Changes to the output format must be reflected there (and formatVersion increased).
 *
 * Output File structure:
 * @code
 * formatVersion: 1
 * [General information]
 * [Data]
 * simstep: [...]
 * [Process 0]
 * traversalTime: [...]
 * ...
 * [Process 1]
 * ...
 * @endcode
 * The file is line based:
 * - A line of the form "[name]" starts a new section.
 * - Every other non-empty line has the form "key: value".
 * - Lists are written as "[v0,v1,...]", an empty list as "[]".
 * - Every line, including the last one, ends with '\n'.
 */
class DomainDecompAnalyser : public PluginBase  {
public:

    DomainDecompAnalyser();
	~DomainDecompAnalyser() {}


    void readXML(XMLfileUnits& xmlconfig);


	void init(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain);

	void endStep(ParticleContainer *particleContainer, DomainDecompBase *domainDecomp, Domain *domain, unsigned long simstep);

	void finish(ParticleContainer *particleContainer,
				DomainDecompBase *domainDecomp, Domain *domain);

	std::string getPluginName() { return std::string("DomainDecompAnalyser");}
	static PluginBase* createInstance() { return new DomainDecompAnalyser();}

private:
	bool _appendTimestamp;
	std::string _outputPrefix;

	std::vector<unsigned long> _simstep;
	std::vector<double> _traversalTime;
	std::vector<double> _averageTraversalTime;
	std::vector<unsigned long> _numberParticles;
	
	
	std::vector<unsigned long> _rebuildStep;
	#ifdef ENABLE_MPI
		std::vector<DomainBox> _rebuildDomainBox;
	#endif

	/**
	 * Writes a vector as "name: [v0,v1,...]" followed by a newline.
	 * @param os The stream to write to
	 * @param name The name written in front of the values
	 * @param values The values to write
	 */
	template <typename T>
	static void writeVector(std::ostream& os, const std::string& name, const std::vector<T>& values);
};
