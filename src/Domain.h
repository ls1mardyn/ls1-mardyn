#ifndef DOMAIN_H_
#define DOMAIN_H_

#include <string>
#include <fstream>
#include <climits>

#include "molecules/Comp2Param.h"
#include "molecules/Component.h"
#include "utils/Log.h"

class Molecule;

namespace datastructures{
  template<class ParticleType>
  class ParticleContainer;
}

namespace parallel{
  class DomainDecompBase; 
}


using namespace std;

//! @brief This class is used to read in the phasespace and to handle macroscopic values
//!
//! The main purpose of this class is responsible for all macroscopic value (Potential,...).
//! It is important to differentiate between local and global version of those values
//! Macroscopic values are values that aggregate some property of a set of molecules.
//! As this program is designed to run on a parallel computer, there are typically 
//! several processes. Each process has an instance of this class, but with a different 
//! subset of molecules. Therefore, also the macroscopic values are only representative 
//! for the "local" domain of that process. Thus, local member variables that represent 
//! "local" macroscopic values, begin with _local...
//! 
//! At some points of the simulation, macroscopic values for the whole set of molecules 
//! have to be calculated. Those values are stored in member variables beginning with _global. 
//!
//!  
//! The phasespace file consists of two parts, a header (including component description)
//! and the molecule data. Header and data have to be read in with two different methods. 
//! Thats' because the molecules from the data part have to be stored in a molecule 
//! container, which can only be created with information from the header part.
//! 
//! Example for the Phasespace file:
//! <TABLE><TR><TD><CODE>
//! MDProject 20070123\n
//! currenttime "t (double)" \n
//! Temperature "T (double)" \n
//! Length "x (double)" "y (double)" "z (double)" \n
//! NumberOfComponents "C (int)" \n
//!    [#LJ]_C1 [#DP]_C1 [#QP]_C1
//!    [xlj_1]_C1 [ylj_1]_C1 [zlj_1]_C1 [m_1]_C1 [eps_1]_C1 [sigma_1]_C1 \n
//!    ... \n
//!    [xlj_#LJ]_C1 [ylj_#LJ]_C1 [zlj_#LJ]_C1 [m_#LJ]_C1 [eps_#LJ]_C1 [sigma_#LJ]_C1 \n
//!    [xdp_1]_C1 [ydp_1]_C1 [zdp_1]_C1 [eMyx_1]_C1 [eMyy_1]_C1 [eMyz_1]_C1 [absMy_1]_C1 \n
//!    ... \n
//!    [xdp1_#DP]_C1 [ydp1_#DP]_C1 [zdp1_#DP]_C1 [eMyx_#DP]_C1 [eMyy_#DP]_C1 [eMyz_#DP]_C1 [absMy_#DP]_C1 \n
//!    [xqp_1]_C1 [yqp_1]_C1 [zqp_1]_C1 [eMyx_1]_C1 [eMyy_1]_C1 [eMyz_1]_C1 [absMy_1]_C1 \n
//!    ... \n
//!    [xqp1_#QP]_C1 [yqp1_#QP]_C1 [zqp1_#QP]_C1 [eQx_#QP]_C1 [eQy_#QP]_C1 [eQz_#QP]_C1 [absQ_#QP]_C1 \n
//!    [I11]_C1 [I22]_C1 [I33]_C1 \n
//!    
//!    ... \n
//!  
//!    [#LJ]_Cn [#DP]_Cn [#QP]_Cn
//!    [xlj_1]_Cn [ylj_1]_Cn [zlj_1]_Cn [m_1]_Cn [eps_1]_Cn [sigma_1]_Cn \n
//!    ... \n
//!    [xlj_#LJ]_Cn [ylj_#LJ]_Cn [zlj_#LJ]_Cn [m_#LJ]_Cn [eps_#LJ]_Cn [sigma_#LJ]_Cn \n
//!    [xdp_1]_Cn [ydp_1]_Cn [zdp_1]_Cn [eMyx_1]_Cn [eMyy_1]_Cn [eMyz_1]_Cn [absMy_1]_Cn \n
//!    ... \n
//!    [xdp1_#DP]_Cn [ydp1_#DP]_Cn [zdp1_#DP]_Cn [eMyx_#DP]_Cn [eMyy_#DP]_Cn [eMyz_#DP]_Cn [absMy_#DP]_Cn \n
//!    [xqp_1]_Cn [yqp_1]_Cn [zqp_1]_Cn [eMyx_1]_Cn [eMyy_1]_Cn [eMyz_1]_Cn [absMy_1]_Cn \n
//!    ... \n
//!    [xqp1_#QP]_Cn [yqp1_#QP]_Cn [zqp1_#QP]_Cn [eQx_#QP]_Cn [eQy_#QP]_Cn [eQz_#QP]_Cn [absQ_#QP]_Cn \n
//!    [I11]_Cn [I22]_Cn [I33]_Cn \n
//!
//! NumberOfMolecules "N" \n
//! [id_1] [type_1] [x_1] [y_1] [z_1] [vx_1] [vy_1] [vz_1] [q0_1] [q1_1] [q2_1] [q3_1] [Dx_1] [Dy_1] [Dz_1]\n
//! [id_2] [type_2] [x_2] [y_2] [z_2] [vx_2] [vy_2] [vz_2] [q0_2] [q1_2] [q2_2] [q3_2] [Dx_2] [Dy_2] [Dz_2]\n
//! ... \n
//! [id_N] [type_N] [x_N] [y_N] [z_N] [vx_N] [vy_N] [vz_N] [q0_N] [q1_N] [q2_N] [q3_N] [Dx_N] [Dy_N] [Dz_N]\n 
//! </CODE></TD></TR></TABLE>
class Domain{
  public:
    //! The constructor sets _localRank to rank and initializes all member variables 
    Domain(int rank);

    //! @brief gets a filename and opens an ifstream associated with the given file
    //! 
    //! As the reading of the phasespace file is separated into two parts,
    //! but each line of the file should only be parsed once, not the filename
    //! itself is stored, but a stream (_phaseSpaceFileStream) which is associated with
    //! the file
    //! @param filename full path to the input file
    void setPhaseSpaceFile(string filename);  
        
    //! @brief reads in header of the input file (including component description)
    //!
    //! The Header in the input file consists of several elements. An element starts
    //! in a new line with the element name followed by some whitespace (e.g. "Temperature ").
    //! After that, there can be any number of tokens belonging to that element. 
    //! A description of the different elements follows below. But first
    //! some rules dealing with the order of the elements:
    //! \li The first header element must start with the string "MDProject" 
    //!     followed by a timestamp.
    //! \li The last header element must start with the string NumberOfMolecules.
    //!     After that, the header is over and the molecule data follows
    //! \li The order of the remaining header lines is not important.
    //! \li Header elements beginning with "#" are ignored. 
    //!
    //! Now the description of the different elements:
    //! \li MDProject: One token follows with a version number (see description of _inpversion);
    //! \li currenttime: One token follows with the start time
    //! \li Temperature: One token follows with the temperature
    //! \li Length: Three tokens follow with the length of the simulation box
    //!     in the three dimensions (x, y and z)
    //! \li NumberOfComponents: Here follow several tokens, the first one is the
    //!     actual Number of Components.
    //!     Then the values describing the components have to follow, seperated by 
    //!     whitespace. For each component, the following values have to be provided:
    //!     - Number of Lennard-Jones-Centers, Number of Dipoles, Number of Quadrupoles (all int)
    //!     - For each LJ-Center: x-coord., y-coord., z-coord., mass, epsilon, sigma (all double)
    //!     - For each Dipole: x-coord., y-coord., z-coord., eMyx, eMyy, eMyz, absMy (all double)
    //!     - For each Quadrupole: x-coord., y-coord., z-coord., eQx, eQy, eQz, absQ (all double)
    //!     - moments of inertia for principal axes: I11, I22, I33 (all double)
    //!     - For each pair of different components: xi, eta (both double)
    //!     - epsilonRF (double)
    //! \li NumberOfMolecules: One token follows with the number of molecules
    //!
    //! An example can be seen in the documentation of this class
    void readPhaseSpaceHeader();

    //! @brief reads in the data of all molecules
    //! 
    //! The Molecule Data starts in a new line with the string "MoleculeFormat"
    //! followed by whitespace and a string representing the format.
    //! In the standard case (format ICRVQD), the following values are provided for each molecule:
    //! \li id of the molecule (int)
    //! \li id of the component of the molecule (int)
    //! \li Coordinates: x, y, z (all double)
    //! \li velocities: vx, vy, vz (all double)
    //! \li Orientation (quaternion): q0, q1, q2, q3 (all double)
    //! \li Angular Momentum: Dx, Dy, Dz (all double)
    //! 
    //! An example can be seen in the documentation of this class
    //! @param particleContainer Here the Molecules from the input file are stored 
    void readPhaseSpaceData(datastructures::ParticleContainer<Molecule>* particleContainer);

    //! @brief writes a checkpoint file that can be used to continue the simulation
    //!
    //! The format of the checkpointfile written by this method is the same as the format
    //! of the input file. 
    //! @param filename Name of the checkpointfile (including path)
    //! @param particleContainer The molecules that have to be written to the file are stored here
    //! @param domainDecomp In the parallel version, the file has to be written by more than one process.
    //!                     Methods to achieve this are available in domainDecomp 
    void writeCheckpoint(string filename, datastructures::ParticleContainer<Molecule>* particleContainer,
                         parallel::DomainDecompBase* domainDecomp) const;

    //! @brief initialize far field correction parameters
    //! 
    //! By limiting the calculation to pairs of particles which have
    //! less distance than the given cutoff radius, an error is made
    //! By calculating approximations for the neglected pairs, the
    //! error can be reduced
    //! @param cutoffRadius cutoff radius
    //! @todo How does the Correction work? Give reference to some paper,
    //!       documentation in the implementation
    void initFarFieldCorr(double cutoffRadius); // CHECKED

    //! @brief initialize parameter streams
    //!
    //! This method should only be called, after the the component information
    //! and all molecule data have been read in
    //! @param cutoffRadius cutoff radius
    void initParameterStreams(double cutoffRadius);

    //! @brief set the potential of the local process
    void setLocalUpot(double Upot); 

    //! @brief get the potential of the local process
    double getLocalUpot() const;

    //! @brief set the virial of the local process
    void setLocalVirial(double Virial); 

    //! @brief get the virial of the local process
    double getLocalVirial() const;   
    
    //! @brief get thermostat scaling for translations 
    double getGlobalBetaTrans() const;

    //! @brief get thermostat scaling for rotations 
    double getGlobalBetaRot() const;

    //! @brief return the length of the domain
    //!
    //! @param index dimension for which the length should be returned
    double getGlobalLength(int index) const;
    
    //! @brief set the length of the domain
    //!
    //! @param index dimension for which the length should be set
    //! @param index value which should be set
    void setGlobalLength(int index, double length);

    //! @brief get the global temperature
    double getGlobalTemperature() const;

    //! @brief set the global temperature
    void setGlobalTemperature(double temp);

    //! @brief get the mixcoeff
    vector<double> & getmixcoeff();

    //! @brief get the epsilonRF
    double getepsilonRF() const;

    //! @brief set the epsilonRF
    void setepsilonRF(double erf);

    //! @brief get globalNumMolecules
    unsigned long getglobalNumMolecules() const;

    //! @brief set globalNumMolecules
    void setglobalNumMolecules(unsigned long glnummol);

    //! @brief get the global pressure
    //!
    //! @todo provide justification for the formula
    double getGlobalPressure() const;

    //! @brief get the global average potential per particle
    //!
    //! Before this method is called, it has to be sure that the
    //! global potential has been calculated (method calculateGlobalValues)
    double getAverageGlobalUpot() const;

    //! @brief get the global average virial per particle
    //!
    //! Before this method is called, it has to be sure that the
    //! global virial has been calculated (method calculateGlobalValues)
    double getAverageGlobalVirial() const;

    //! @brief sets _localSummv2 to the given value
    void setLocalSummv2(double summv2);
    
    //! @brief sets _localSumIw2 to the given value
    void setLocalSumIw2(double sumIw2);
	
	//! @brief get local rank 
    int getlocalRank();
    
    //! @brief get inpversion
    unsigned long getinpversion();

    //! @brief set inpversion
    void setinpversion(unsigned long inpv);

    //! @brief get globalRho
    double getglobalRho();

    //! @brief set globalRho
    void setglobalRho(double grho);

    //! @brief get globalRotDOF
    unsigned long getglobalRotDOF();

    //! @brief set globalRotDOF
    void setglobalRotDOF(unsigned long grotdof);

    //! @brief get logging interface
    utils::Log getlog();
    
    //! @brief get the current time
    double getCurrentTime();
    
    //! @brief get the current time
    void setCurrentTime(double curtime);
    
    //! @brief advance the current time by timestep
    void advanceTime(double timestep);

    //! @brief get a reference to the vector of Components
    vector<Component>& getComponents();
    
    //! @brief get the parameter streams
    Comp2Param& getComp2Params();
       
    //! @brief calculate the global macroscopic values
    //!
    //! @param domainDecomp domain decomposition
    //! @param particleContainer particle Container
    //! @todo more detailled description
    void calculateGlobalValues(parallel::DomainDecompBase* domainDecomp,
                              datastructures::ParticleContainer<Molecule>* particleContainer);
    
    //! @brief calculate _localSummv2 and _localSumIw2
    //! @todo more detailled description
    void calculateVelocitySums(datastructures::ParticleContainer<Molecule>* partCont);
    
  private:
    //! Logging interface
    static utils::Log _log;  
  
    //! rank of the local process
    int _localRank;

    //! file stream associatted to the input file
    ifstream _phaseSpaceFileStream;

    //! @brief Version of the input file
    //!
    //! even though it is desirable, that the format of the input file
    //! doesn't change, is sometimes does change. When that happens,
    //! the code which reads in the input file (parser) has to be changed as well.
    //! old versions of the input file then can't be read any more.
    //! So whenever the parser is changed, _inpversion is set to the
    //! date of the change (YYYYMMDD) (hard-coded). Only input files
    //! with the same version are sure to be processed correctly
    unsigned long _inpversion;
    
    //! Potential of the local process
    double _localUpot;
    //! Virial of the local process
    double _localVirial;
    //! global Potential
    double _globalUpot;
    //! global virial
    double _globalVirial;
    //! @todo WHAT IS THIS? 
    double _globalBetaTrans;
    //! @todo WHAT IS THIS?
    double _globalBetaRot;
    //! global density
    double _globalRho;
    //! global Temperature
    double _globalTemperature;
    //! global Number of rotational degrees of freedom
    unsigned long _globalRotDOF;
    //! global Number of Molecules
    //! @todo redundancy?
    unsigned long _globalNumMolecules;
    //! side length of the cubic simulation box
    double _globalLength[3];

    //! local sum (over all molecules) of the mass multiplied with the squared velocity
    //! @todo rename this 
    double _localSummv2;
    //! local sum (over all molecules) of the moment of inertia 
    //! multiplied with the squared  rotational velocity
    //! @todo rename this
    double _localSumIw2; 
    
    //! reaction field
    //! @todo WHAT IS THIS?
    //! @todo local or global?
    double _epsilonRF;

    //! Potential correction for the error made by the cutoff
    //! @todo local or global?
    double _UpotCorr;
    //! Virial correction for the error made by the cutoff
    //! @todo local or global?
    double _VirialCorr;     

    //! @todo comment
    double _currentTime;

    //! Components resp. molecule types
    vector<Component> _components;
    //! parameter streams for each possible pair of molecule-types
    Comp2Param _comp2params;
    //! modified Lorentz-Berthelot mixing rule parameters
    //! @todo more explanation
    vector<double> _mixcoeff;

};


#endif /*DOMAIN_H_*/
