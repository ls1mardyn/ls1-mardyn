/*
 * Created on Wed May 08 2024
 * Author: Jose A. Pinzon Escobar
 * Email to: jose.escobar@hsu-hh.de
 * Copyright (c) 2024 Helmut-Schmidt University, Hamburg
 */

#pragma once

#include "../features/FTH_Grid.h"
#include "../Interpolation.h"
#include "plugins/AdResS/util/Averager.h"

class ParticleContainer;
class DomainDecompBase;
class Domain;

#if defined(__GNUC__)
	#define s_inline inline __attribute__((always_inline))
#else
	#define s_inline inline
#endif

class SamplerBase {
public:
	virtual ~SamplerBase() = default;

	/**
	 * Initializes all necessary fields
	 * */
	virtual void init(Domain *domain) = 0;

	/**
	 * Writes the currently stored sample data to the specified file
	 * @param filename filename
	 * @param simstep simstep
	 * */
	virtual void writeSample(const std::string &filename, int simstep) = 0;

	/**
	 * Returns a reference to the currently stored sample data
	 * */
	std::vector<double>& getSampledData() { return _sampled_data; }

	/**
	 * Samples data and stores in buffer
	 * */
	virtual void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) = 0;

	/**
	 * Returns true iff the implementation uses the grid for density computation
	 * */
	virtual bool usesGrid() { return false; }

	/**
	 * Returns true iff the implementation computes a 3d density, else a 1d projection
	 * */
	virtual bool is3D() { return false; }

	/**
	 * If the implementation is NOT 3d and NOT using the grid, then this method provides the output
	 * @param fun function destination
	 * */
	virtual void getSampleFunction(Interpolation::Function& fun) { };

	/**
	 * If the implementation is 3d and NOT using the grid, then this method provides the output
	 * @param fun function destination
	 * */
	virtual void getSampleFunction(Interpolation::Function3D& fun) { };

	/**
	 * Returns number of bins if implementation is projected
	 * */
	virtual int getNumSamplePoints() { return 0; };

	/**
	 * Returns true if sampleData was called at least once before
	 * */
	[[nodiscard]] bool wasSampled() const { return _was_sampled; }

protected:
	/// sample buffer
	std::vector<double> _sampled_data;
	/// flag whether sampleData was already called at least once
	bool _was_sampled = false;
};

class GridSampler : public SamplerBase {
public:
	/**
	 * Construct a GridSampler
	 * @param grid ptr to used grid
	 * @param rad radius of sphere around each node for measurement
	 * */
    GridSampler(FTH::grid_t *grid, double rad);

    /**
     * Do resizing and other stuff
    */
	void init(Domain *domain) override;

	/**
	 * Populates sample buffer
	 * */
    void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

    /**
     * A grid sampler would print all data on every node
    */
    void writeSample(const std::string &filename, int simstep) override;

	bool usesGrid() override { return true; }

	bool is3D() override { return true; }

protected:
	/**
     * Samples densities at all nodes using only the local particles
     * @param pc active particle container
     * */
	void sampleAtNodes(ParticleContainer* pc);
	/**
	 * \brief Check if particle within distance using Euclidean norm
	*/
	[[nodiscard]] s_inline
	bool ParticleInsideMeasuringSpace(const std::array<double,3>& n_pos, const std::array<double,3>& l_pos) const;

	/// pointer to used grid, not owned by this class
	FTH::grid_t *_grid;
	/// radius of sphere around each node for measurement
	double _measure_radius;
};

class AveragedGridSampler : public GridSampler {
public:
	/**
	 * Construct an AveragedGridSampler
	 * @param grid ptr to used grid
	 * @param rad radius of sphere around each node for measurement
	 * */
    AveragedGridSampler(FTH::grid_t *grid, double rad);

	/**
	 * Returns a mutable reference to the averager
	 * */
    Averager<std::vector<double>>& getAverager();

	/**
     * Do resizing and other stuff
    */
	void init(Domain *domain) override;

	/**
	 * Populates sample buffer. Calls super class sampleData and uses results for averaging
	 * */
    void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
     * A grid sampler would print all data on every node
    */
    void writeSample(const std::string &filename, int simstep) override;
protected:
	/// averages data across all sampleData calls
	Averager<std::vector<double>> _averager;
};

/**
 * Projects all molecular positions down to the selected dimension and performs density sampling on that
 * */
class ProjectedSampler : public SamplerBase {
public:
	/**
	 * Creates a projection sampler, that projects to the selected dimension
	 * @param dim dimension
	 * */
	explicit ProjectedSampler(int dim) : _dim(dim) { }

protected:
	/**
	 * Populates _global_mol_pos with the molecular positions of all ranks
	 * */
	void loadGlobalMolPos(ParticleContainer* particleContainer, DomainDecompBase* domainDecomp, Domain* domain);

	/// projection dimension
	int _dim;
	/// MPI buffer - mol pos of all ranks
	std::vector<double> _global_mol_pos;
private:
	/// MPI buffer - mol pos of only local rank
	std::vector<double> _local_mol_pos;
};

/**
 * Splits selected 1D axis into bins of selected width and counts molecules per bin
 * */
class DirectProjectedSampler : public ProjectedSampler {
public:
	/**
	 * Creates direct projected sampler
	 * @param dim projection dimension
	 * @param bin_width width of each bin
	 * */
	DirectProjectedSampler(int dim, double bin_width);

	// init some fields that need the domain
	void init(Domain *domain) override;

	// "uses grid" by having bins, will output data in _sampled_data
	bool usesGrid() override { return true; }

	/**
	 * Samples all molecules into their bins and stores result in _sampled_data
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
	 * Writes _sampled_data to file
	 * */
	void writeSample(const std::string &filename, int simstep) override;

	int getNumSamplePoints() override { return static_cast<int>(_bins); }

protected:
	/// number of bins
	unsigned long _bins;
	/// width of bin in domain
	double _binWidth;
	/// volume of each bin
	double _binVolume;
	/// MPI buffer - densities of only local rank
	std::vector<double> _localDensities;
};

/**
 * Smoothes the output of DirectProjectedSampler using a Gaussian kernel
 * */
class SmoothedProjectedSampler : public DirectProjectedSampler {
public:
	/**
	 * Creates smoothed projected sampler
	 * @param dim projection dimension
	 * @param bin_width width of each bin
	 * @param smoothing_strength strength of smoothing kernel
	 * */
	SmoothedProjectedSampler(int dim, double bin_width, double smoothing_strength);

	// init some fields that need the domain
	void init(Domain *domain) override;

	/**
	 * Calls DirectProjectedSampler::sampleData and smoothes its results
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

private:
	/// smoothing kernel
	Interpolation::Matrix _smoothingFilter;
	/// smoothing kernel strength
	double _filterStrength;
};

/**
 * Assumes the domain is a space of dirac delta functions. Everywhere, where a molecule is, a dirac is 1.
 * Uses a real FT to filter out high frequencies.
 * */
class FTProjectedSampler : public ProjectedSampler {
public:
	/**
	 * Creates a ft projected sampler
	 * @param dim projections dimension
	 * @param frequencies number of considered frequencies
	 * @param samples number of resampling points
	 * */
	FTProjectedSampler(int dim, int frequencies, int samples);

	// init some fields that need the domain
	void init(Domain *domain) override;

	/**
	 * Samples density by using real FT filtering
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
	 * Writes _fun to file
	 * */
	void writeSample(const std::string &filename, int simstep) override;

	/**
	 * Writes into fun the current sampled function
	 * */
	void getSampleFunction(Interpolation::Function& fun) override;

	int getNumSamplePoints() override { return _samples; }

private:
	/// number of discrete frequencies
	int _frequencies;
	/// number of samples to resample the ift output
	int _samples;
	/// buffer for real FT
	std::vector<std::complex<double>> _realFT;
	/// current output function
	Interpolation::Function _fun;
	/// volume of each bin
	double _binVolume;
};

/**
 * Creates GMM around all molecular positions
 * */
class GMMProjectedSampler : public ProjectedSampler {
public:
	/**
	 * Creates a gmm projected sampler
	 * @param dim projections dimension
	 * @param samples number of resampling points
	 * @param smoothing_strength strength of filtering
	 * */
	GMMProjectedSampler(int dim, int samples, double smoothing_strength);

	/// Does nothing
	void init(Domain *domain) override { }

	/**
	 * Samples density by using GMM filtering
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
	 * Writes _fun to file
	 * */
	void writeSample(const std::string &filename, int simstep) override;

	/**
	 * Writes into fun the current sampled function
	 * */
	void getSampleFunction(Interpolation::Function& fun) override;

	int getNumSamplePoints() override { return _samples; }

private:
	/// number of samples for resampling
	int _samples;
	/// strength of filtering
	double _filterStrength;
	/// current output function
	Interpolation::Function _fun;
};
