#ifndef ENSEMBLE_BASE_H
#define ENSEMBLE_BASE_H


//! list of updateable values
enum GlobalVariable {
	NUM_PARTICLES      = 1<<0,
	ENERGY             = 1<<1,
	VOLUME             = 1<<2,
	CHEMICAL_POTENTIAL = 1<<3,
	TEMPERATURE        = 1<<4,
	PRESSURE           = 1<<5
};


//! @brief Base class for ensembles
//! @author Christoph Niethammer <niethammer@hlrs.de>
//! 
//! Each ensemble should provide access to extensive (NVE) and intensive 
//! (\mu p t) variables as well as a function to update global variables.
class Ensemble {
public:
	Ensemble() {}
	virtual ~Ensemble() {}

	//! @brief Returns the global number of Molecules of the ensemble.
	virtual unsigned long N() = 0;
	//! @brief Returns the global volume of the ensemble
	virtual double V() = 0;
	//! @brief Returns the global energy of the ensemble
	virtual double E() = 0;
	//! @brief Returns the global chemical potential of the ensemble
	virtual double mu() = 0;
	//! @brief Returns the global presure of the ensemble.
	virtual double p() = 0;
	//! @brief Returns the global Temperature of the ensemble.
	virtual double T() = 0;

	//! @brief Calculate global variables
	//! @param variable Variable to be updated.
	virtual void updateGlobalVariable(GlobalVariable variable) = 0;
};

#endif
