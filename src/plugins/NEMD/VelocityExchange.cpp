
#include "plugins/NEMD/VelocityExchange.h"
#include "particleContainer/ParticleContainer.h"
#include "Domain.h"
#include "parallel/DomainDecompBase.h"
#include "molecules/Molecule.h"
#include "utils/Logger.h"
#include "utils/CommVar.h"

#include <vector>
#include <cmath>
#include <algorithm>

#include "utils/xmlfileUnits.h"
#include "utils/FileUtils.h"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

using namespace std;
using Log::global_log;

VelocityExchange::VelocityExchange(){
}

VelocityExchange::~VelocityExchange() {
}

void VelocityExchange::init(ParticleContainer* particleContainer,
		  DomainDecompBase* domainDecomp, Domain* domain)
{
	global_log->info() << "VelocityExchange enabled" << endl;

  for(unsigned d = 0; d < 3; d++){
    _boxLength[d] = domain->getGlobalLength(d);
  }

}

void VelocityExchange::readXML(XMLfileUnits& xmlconfig)
{
	// Timestep control
	_control.start = 0;
	_control.freq = 10000;
	_control.stop = 1000000000;
	xmlconfig.getNodeValue("control/start", _control.start);
	xmlconfig.getNodeValue("control/frequency", _control.freq);
	xmlconfig.getNodeValue("control/stop", _control.stop);
	global_log->info() << "[VelocityExchange] takes place start:freq:stop = "
			<< _control.start << ":" << _control.freq << ":" << _control.stop << endl;

  	// range
	Domain* domain = global_simulation->getDomain();
	_cold_range.xmin = 0.;
	_cold_range.xmax = domain->getGlobalLength(0);
	_cold_range.ymin = 0.;
	_cold_range.ymax = domain->getGlobalLength(1);
	_cold_range.zmin = 0.;
	_cold_range.zmax = domain->getGlobalLength(2);
  _warm_range.xmin = 0.;
	_warm_range.xmax = domain->getGlobalLength(0);
	_warm_range.ymin = 0.;
	_warm_range.ymax = domain->getGlobalLength(1);
	_warm_range.zmin = 0.;
	_warm_range.zmax = domain->getGlobalLength(2);
  symmetry = false;

	xmlconfig.getNodeValue("coldrange/xmin", _cold_range.xmin);
	xmlconfig.getNodeValue("coldrange/xmax", _cold_range.xmax);
	xmlconfig.getNodeValue("coldrange/ymin", _cold_range.ymin);
	xmlconfig.getNodeValue("coldrange/ymax", _cold_range.ymax);
	xmlconfig.getNodeValue("coldrange/zmin", _cold_range.zmin);
	xmlconfig.getNodeValue("coldrange/zmax", _cold_range.zmax);

  xmlconfig.getNodeValue("warmrange/symmetric", symmetry);
  xmlconfig.getNodeValue("warmrange/xmin", _warm_range.xmin);
	xmlconfig.getNodeValue("warmrange/xmax", _warm_range.xmax);
	xmlconfig.getNodeValue("warmrange/ymin", _warm_range.ymin);
	xmlconfig.getNodeValue("warmrange/ymax", _warm_range.ymax);
	xmlconfig.getNodeValue("warmrange/zmin", _warm_range.zmin);
	xmlconfig.getNodeValue("warmrange/zmax", _warm_range.zmax);

  if ( symmetry ) {
    global_log->info() << "[VelocityExchange] symmetric warm domain wrt. x-z-plane "
        << _control.start << ":" << _control.freq << ":" << _control.stop << endl;
  }
  
  uint32_t numComp = global_simulation->getEnsemble()->getComponents()->size();

  resizeExactly(v_c_abs_local, numComp);
  resizeExactly(v_w_abs_local, numComp);
  resizeExactly(v_c_abs_global, numComp);
  resizeExactly(v_w_abs_global, numComp);

  resizeExactly(cold_mol_local, numComp);
  resizeExactly(warm_mol_local, numComp);
  resizeExactly(cold_mol_global, numComp);
  resizeExactly(warm_mol_global, numComp);
  
  resizeExactly(v_cold_local, numComp*3);
  resizeExactly(v_warm_local, numComp*3);
  resizeExactly(D_cold_local, numComp*3);
  resizeExactly(D_warm_local, numComp*3);

  resizeExactly(v_cold_global, numComp*3);
  resizeExactly(v_warm_global, numComp*3);
  resizeExactly(D_cold_global, numComp*3);
  resizeExactly(D_warm_global, numComp*3);

  std::fill(v_c_abs_local.begin(), v_c_abs_local.end(), 1000.0 );
  std::fill(v_w_abs_local.begin(), v_w_abs_local.end(), 0.0 );
  std::fill(v_c_abs_global.begin(), v_c_abs_global.end(), 0.0);
  std::fill(v_w_abs_global.begin(), v_w_abs_global.end(), 0.0);
  
  std::fill(cold_mol_local.begin(), cold_mol_local.end(), 0.0); 
  std::fill(cold_mol_global.begin(), cold_mol_global.end(), 0.0);
  std::fill(warm_mol_local.begin(), warm_mol_local.end(), 0.0); 
  std::fill(warm_mol_global.begin(), warm_mol_global.end(), 0.0);

}

void VelocityExchange::beforeForces(
		ParticleContainer* particleContainer, DomainDecompBase* domainDecomp,
		unsigned long simstep
	)
{
	if (simstep < _control.start || simstep > _control.stop
			|| simstep % _control.freq != 0)
		return;
	this->exchangeVelocities(particleContainer, domainDecomp, simstep);
}

void VelocityExchange::exchangeVelocities(ParticleContainer* particleContainer,
		DomainDecompBase* domainDecomp, unsigned long simstep)
{

  uint32_t numComp = global_simulation->getEnsemble()->getComponents()->size();

  std::fill(v_c_abs_local.begin(), v_c_abs_local.end(), 1000.0 );
  std::fill(v_w_abs_local.begin(), v_w_abs_local.end(), 0.0 );
  std::fill(v_c_abs_global.begin(), v_c_abs_global.end(), 0.0);
  std::fill(v_w_abs_global.begin(), v_w_abs_global.end(), 0.0);
  std::fill(cold_mol_global.begin(), cold_mol_global.end(), 0.0);
  std::fill(warm_mol_global.begin(), warm_mol_global.end(), 0.0);

	double c_regionLowCorner[3], c_regionHighCorner[3];

	// if linked cell in the region of interest
	for (unsigned d = 0; d < 3; d++) {
		c_regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		c_regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}

	// ensure that we do not iterate over things outside of the container.
	c_regionLowCorner[0] = std::max(_cold_range.xmin, c_regionLowCorner[0]);
	c_regionLowCorner[1] = std::max(_cold_range.ymin, c_regionLowCorner[1]);
	c_regionLowCorner[2] = std::max(_cold_range.zmin, c_regionLowCorner[2]);

	c_regionHighCorner[0] = std::min(_cold_range.xmax, c_regionHighCorner[0]);
	c_regionHighCorner[1] = std::min(_cold_range.ymax, c_regionHighCorner[1]);
	c_regionHighCorner[2] = std::min(_cold_range.zmax, c_regionHighCorner[2]);

  double w_regionLowCorner[3], w_regionHighCorner[3];

	// if linked cell in the region of interest
	for (unsigned d = 0; d < 3; d++) {
		w_regionLowCorner[d] = particleContainer->getBoundingBoxMin(d);
		w_regionHighCorner[d] = particleContainer->getBoundingBoxMax(d);
	}

	// ensure that we do not iterate over things outside of the container.
	w_regionLowCorner[0] = std::max(_warm_range.xmin, w_regionLowCorner[0]);
	w_regionLowCorner[1] = std::max(_warm_range.ymin, w_regionLowCorner[1]);
	w_regionLowCorner[2] = std::max(_warm_range.zmin, w_regionLowCorner[2]);

	w_regionHighCorner[0] = std::min(_warm_range.xmax, w_regionHighCorner[0]);
	w_regionHighCorner[1] = std::min(_warm_range.ymax, w_regionHighCorner[1]);
	w_regionHighCorner[2] = std::min(_warm_range.zmax, w_regionHighCorner[2]);

  double v_dummy[3] = {0.0, 0.0, 0.0};


  Molecule* warm_mol_ptr[numComp];
  Molecule* cold_mol_ptr[numComp];

  // find warmest mol in cold region 
  auto begin_c = particleContainer->regionIterator(c_regionLowCorner, c_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
	for(auto it = begin_c; it.isValid(); ++it) {
    uint32_t cid = it->componentid();   
    v_dummy[0] = it->v(0);
    v_dummy[1] = it->v(1);
    v_dummy[2] = it->v(2);
    double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
    if ( v_dummy_abs > v_w_abs_local[cid] ){
      v_w_abs_local[cid] = v_dummy_abs;
      warm_mol_ptr[cid] = &(*it);
    }
  }

  // find coldest mol in warm region
  auto begin_w = particleContainer->regionIterator(w_regionLowCorner, w_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
	for(auto it = begin_w; it.isValid(); ++it) {
    uint32_t cid = it->componentid();   
    v_dummy[0] = it->v(0);
    v_dummy[1] = it->v(1);
    v_dummy[2] = it->v(2);
    double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
    if ( v_dummy_abs < v_c_abs_local[cid] ){
      v_c_abs_local[cid] = v_dummy_abs;
      cold_mol_ptr[cid] = &(*it);
    }
  }
  
  if ( symmetry ) { // warm region is symmetric

    double ymin_symm = _boxLength[1] - _warm_range.ymax;
    double ymax_symm = _boxLength[1] - _warm_range.ymin;
    w_regionLowCorner[1] = particleContainer->getBoundingBoxMin(1);
    w_regionHighCorner[1] = particleContainer->getBoundingBoxMax(1);

	  w_regionLowCorner[1] = std::max(ymin_symm, w_regionLowCorner[1]);
	  w_regionHighCorner[1] = std::min(ymax_symm, w_regionHighCorner[1]);

    auto begin_w = particleContainer->regionIterator(w_regionLowCorner, w_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
    for(auto it = begin_w; it.isValid(); ++it) {
      uint32_t cid = it->componentid();   
      v_dummy[0] = it->v(0);
      v_dummy[1] = it->v(1);
      v_dummy[2] = it->v(2);
      double v_dummy_abs = v_dummy[0]*v_dummy[0]+v_dummy[1]*v_dummy[1]+v_dummy[2]*v_dummy[2];
      if ( v_dummy_abs < v_c_abs_local[cid] ){
        v_c_abs_local[cid] = v_dummy_abs;
        cold_mol_ptr[cid] = &(*it);
      }
    }
  }

// find rank with coldest/warmest molecule and read velocities

// std::cout << domainDecomp->getRank() << " cl " << v_c_abs_local[0] << " " <<  v_c_abs_local[1] << std::endl;
// std::cout << domainDecomp->getRank() << " wl " << v_w_abs_local[0] << " " <<  v_w_abs_local[1] << std::endl;
// std::cout << domainDecomp->getRank() << " cg " << v_c_abs_global[0] << " " <<  v_c_abs_global[1] << std::endl;
// std::cout << domainDecomp->getRank() << " wg " << v_w_abs_global[0] << " " <<  v_w_abs_global[1] << std::endl;


#ifdef ENABLE_MPI

  MPI_Barrier(MPI_COMM_WORLD);  

  MPI_Allreduce( v_c_abs_local.data(), v_c_abs_global.data(), numComp, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce( v_w_abs_local.data(), v_w_abs_global.data(), numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else

  v_c_abs_global = v_c_abs_local;
  v_w_abs_global = v_w_abs_local;

#endif

// std::cout << domainDecomp->getRank() << " cg " << v_c_abs_global[0] << " " <<  v_c_abs_global[1] << std::endl;
// std::cout << domainDecomp->getRank() << " wg " << v_w_abs_global[0] << " " <<  v_w_abs_global[1] << std::endl;

// std::cout << domainDecomp->getRank() << " p3 " << std::endl;
std::fill(v_cold_local.begin(), v_cold_local.end(), -100000.0 );
std::fill(v_warm_local.begin(), v_warm_local.end(), -100000.0 );
std::fill(D_cold_local.begin(), D_cold_local.end(), -100000.0 );
std::fill(D_warm_local.begin(), D_warm_local.end(), -100000.0 );
std::fill(v_cold_global.begin(), v_cold_global.end(), 0.0 );
std::fill(v_warm_global.begin(), v_warm_global.end(), 0.0 );
std::fill(D_cold_global.begin(), D_cold_global.end(), 0.0 );
std::fill(D_warm_global.begin(), D_warm_global.end(), 0.0 );
std::fill(cold_mol_local.begin(), cold_mol_local.end(), 0.0); 
std::fill(cold_mol_global.begin(), cold_mol_global.end(), 0.0);
std::fill(warm_mol_local.begin(), warm_mol_local.end(), 0.0); 
std::fill(warm_mol_global.begin(), warm_mol_global.end(), 0.0);

  for(unsigned k = 0; k < numComp; k++){
    if (abs(v_c_abs_local[k]-v_c_abs_global[k])<1e-9 && v_c_abs_local[k]<1000.0 ) {
      cold_mol_local[k] = cold_mol_ptr[k]->getID();
      for(unsigned d = 0; d < 3; d++){
        v_cold_local[d+3*k] = cold_mol_ptr[k]->v(d);
        D_cold_local[d+3*k] = cold_mol_ptr[k]->D(d);  
      }
    }

    if (abs(v_w_abs_local[k]-v_w_abs_global[k])<1e-9 && v_w_abs_local[k]>0.0 ) {
      warm_mol_local[k] = warm_mol_ptr[k]->getID();
      for(unsigned d = 0; d < 3; d++){
        v_warm_local[d+3*k] = warm_mol_ptr[k]->v(d);
        D_warm_local[d+3*k] = warm_mol_ptr[k]->D(d); 
      }
    }

  }


// // distribute pID and velocities to all ranks
#ifdef ENABLE_MPI

// std::cout << domainDecomp->getRank() << " mpi " << std::endl;

  MPI_Allreduce( warm_mol_local.data(), warm_mol_global.data(), numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce( cold_mol_local.data(), cold_mol_global.data(), numComp, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  MPI_Allreduce( v_cold_local.data(), v_cold_global.data(), 3*numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce( v_warm_local.data(), v_warm_global.data(), 3*numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce( D_cold_local.data(), D_cold_global.data(), 3*numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce( D_warm_local.data(), D_warm_global.data(), 3*numComp, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else
// std::cout << domainDecomp->getRank() << " hm " << std::endl;

  warm_mol_global = warm_mol_local;
  cold_mol_global = cold_mol_local;
  v_cold_local = v_cold_global;
  v_warm_local = v_warm_global;
  D_cold_local = D_cold_global;
  D_warm_local = D_warm_global;

#endif

// for(unsigned d = 0; d < 2; d++){
//   std::cout << "cold " <<  cold_mol_global[d] << std::endl;
// }
// for(unsigned d = 0; d < 2; d++){
//   std::cout << "warm " << warm_mol_global[d] << std::endl;
// }
// for(unsigned d = 0; d < 3; d++){
//   std::cout << "D_cold " << D_cold_global[d] << std::endl;
// }
// for(unsigned d = 0; d < 3; d++){
//   std::cout << "D_warm " << D_warm_global[d] << std::endl;
// }

	for(auto it = begin_c; it.isValid(); ++it) {
    uint32_t cid = it->componentid();   
    if ( it->getID() == warm_mol_global[cid] ) {
        // std::cout << domainDecomp->getRank() << " p5 " << cid << std::endl;
        // std::cout << it->getID() << " vw1 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
        // std::cout << it->getID() << " Dw1 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
      for(unsigned d = 0; d < 3; d++){
        it->setv(d,v_cold_global[d+3*cid]);
        it->setD(d,D_cold_global[d+3*cid]);
      }
        // std::cout << it->getID() << " vw2 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
        // std::cout << it->getID() << " Dw2 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
    }
  }

  for(auto it = begin_w; it.isValid(); ++it) {
    uint32_t cid = it->componentid();   
    if ( it->getID() == cold_mol_global[cid] ) {
      // std::cout << domainDecomp->getRank() << " p5 " << cid << std::endl;
      // std::cout << it->getID() << " vc1 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
      // std::cout << it->getID() << " Dc1 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
      for(unsigned d = 0; d < 3; d++){
        it->setv(d,v_warm_global[d+3*cid]);
        it->setD(d,D_warm_global[d+3*cid]);
      }
      // std::cout << it->getID() << " vc2 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
      // std::cout << it->getID() << " Dc2 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
    }
  }

  if ( symmetry ) {
    auto begin_w_s = particleContainer->regionIterator(w_regionLowCorner, w_regionHighCorner, ParticleIterator::ONLY_INNER_AND_BOUNDARY);  // over all cell types
    for(auto it = begin_w_s; it.isValid(); ++it) {
      uint32_t cid = it->componentid();   
      if ( it->getID() == cold_mol_global[cid] ) {
        // std::cout << domainDecomp->getRank() << " p5 " << cid << std::endl;
        // std::cout << it->getID() << " vs1 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
        // std::cout << it->getID() << " Ds1 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
        for(unsigned d = 0; d < 3; d++){
          it->setv(d,v_warm_global[d+3*cid]);
          it->setD(d,D_warm_global[d+3*cid]);
        }
        // std::cout << it->getID() << " vs2 " << it->v(0) << " " << it->v(1) << " "<< it->v(2) <<   std::endl;
        // std::cout << it->getID() << " Ds2 " << it->D(0) << " " << it->D(1) << " "<< it->D(2) <<   std::endl;
      }
    }
  }

  for(unsigned d = 0; d < numComp; d++){
    global_log->info() << "[VelocityExchange] flip molecule velocities cold:warm " << d << " "
          << cold_mol_global[d] << ":" << warm_mol_global[d] << endl;
  }
}



