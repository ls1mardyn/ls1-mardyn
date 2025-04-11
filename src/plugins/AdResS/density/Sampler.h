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
    /**
     * @param window average window size
	 * @param useReset should avg be reset
     * */
    explicit SamplerBase(int window, bool useReset) : _window(window), _sampled_data(), _was_sampled(false), _averager(), _useReset(useReset) {}

	virtual ~SamplerBase() = default;

	/**
	 * Initializes all necessary fields
	 * */
	virtual void init(Domain *domain) = 0;

    /**
	 * Writes the currently stored averaged data to the specified file
	 * @param filename filename
	 * @param simstep simstep
	 * */
    void writeAverage(const std::string &filename, int simstep);

	/**
	 * Samples data and stores in buffer
	 * */
	virtual void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) = 0;

	/**
	 * Returns true if sampleData was called at least once before
	 * */
	[[nodiscard]] bool wasSampled() const { return _was_sampled; }

    /**
     * Returns the average data
     * */
    virtual std::vector<double> getAverage() = 0;

    /**
     * Performs the reset check
     * */
    void checkReset(unsigned long simstep);

    /**
     * Contains the actual reset logic
     * */
    bool shouldReset(unsigned long simstep);

protected:
    /// averaging window size
    int _window;
	/// sample buffer (containing last instantaneous measurement)
	std::vector<double> _sampled_data;
	/// flag whether sampleData was already called at least once
	bool _was_sampled;
    /// averages data across all sampleData calls
    Averager<std::vector<double>> _averager;
    /// should averages be reset after _window samples?
    bool _useReset;
};

class GridSampler : public SamplerBase {
public:
	/**
	 * Construct a GridSampler
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param grid ptr to used grid
	 * @param rad radius of sphere around each node for measurement
	 * */
    GridSampler(int window, bool useReset, FTH::grid_t *grid, double rad);

    /**
     * Do resizing and other stuff
    */
	void init(Domain *domain) override;

	/**
	 * Populates sample buffer
	 * */
    void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

    /**
     * Returns average of grid
     * */
    std::vector<double> getAverage() override;

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

/**
 * Projects all molecular positions down to the selected dimension and performs density sampling on that
 * */
class ProjectedSampler : public SamplerBase {
public:
	/**
	 * Creates a projection sampler, that projects to the selected dimension
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param dim dimension
	 * */
	explicit ProjectedSampler(int window, bool useReset, int dim) : SamplerBase(window, useReset), _dim(dim) { }

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
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param dim projection dimension
	 * @param bin_width width of each bin
	 * */
	DirectProjectedSampler(int window, bool useReset, int dim, int bins);

	// init some fields that need the domain
	void init(Domain *domain) override;

	/**
	 * Samples all molecules into their bins and stores result in _sampled_data
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

    /**
     * Returns average of direct samples
     * */
    std::vector<double> getAverage() override;

protected:
	/// number of bins
	unsigned long _bins;
	/// width of bin in domain
	double _binWidth;
	/// volume of each bin
	double _binVolume;
};

/**
 * Smoothes the output of DirectProjectedSampler using a Gaussian kernel
 * */
class SmoothedProjectedSampler : public DirectProjectedSampler {
public:
	/**
	 * Creates smoothed projected sampler
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param dim projection dimension
	 * @param bin_width width of each bin
	 * @param smoothing_strength strength of smoothing kernel
	 * */
	SmoothedProjectedSampler(int window, bool useReset, int dim, int bins, double smoothing_strength);

	// init some fields that need the domain
	void init(Domain *domain) override;

    /**
     * Returns smoothed version of direct sampling
     * */
    std::vector<double> getAverage() override;

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
class FTProjectedSampler : public DirectProjectedSampler {
public:
	/**
	 * Creates a ft projected sampler
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param dim projections dimension
	 * @param bins number of resampling points
	 * @param frequencies number of considered frequencies
	 * */
	FTProjectedSampler(int window, bool useReset, int dim, int bins, int frequencies);

	// init some fields that need the domain
	void init(Domain *domain) override;

    /**
     * Returns FT filtered average
     * */
    std::vector<double> getAverage() override;

private:
	/// number of discrete frequencies
	int _frequencies;
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
	 * @param window average window size
	 * @param useReset should avg be reset
	 * @param dim projections dimension
	 * @param bins number of resampling points
	 * @param smoothing_strength strength of filtering
	 * */
	GMMProjectedSampler(int window, bool useReset, int dim, int bins, double smoothing_strength);

	void init(Domain *domain) override;

	/**
	 * Samples density by using GMM filtering
	 * */
	void sampleData(ParticleContainer *pc, DomainDecompBase *domainDecomp, Domain *domain) override;

	/**
	 * Returns GMM filtered average
	 * */
    std::vector<double> getAverage() override;

private:
	/// number of samples for resampling
	int _bins;
	/// strength of filtering
	double _filterStrength;
	/// current output function
	Interpolation::Function _fun;
};
