#include "ChemicalPotential.h"

#include "parallel/DomainDecompBase.h"
#include "particleContainer/ParticleContainer.h"
#include "particleContainer/adapter/ParticlePairs2PotForceAdapter.h"
#include "utils/Logger.h"


ChemicalPotential::ChemicalPotential()
{
	_ownrank = -1;
	_muTilde = 0.0;
	_T = 1.0;
	_interval = (unsigned) ((int) -1);
	_instances = 0;
	for (int d = 0; d < 3; d++) {
		_system[d] = 1.0;
		_minredco[d] = 0.0;
		_maxredco[d] = 1.0;
		_control_bottom[d] = 0.0;
		_control_top[d] = 1.0;
	}
	_nextid = 10000000;
	_globalN = 1;
	_globalV = 1.0;
	_restrictedControlVolume = false;

	_remainingDeletions = std::list<unsigned>();
	for (int d = 0; d < 3; d++)
		_remainingInsertions[d] = std::list<double>();
	_remainingInsertionIDs = std::list<unsigned long>();
	_remainingDecisions = std::list<float>();
	_reservoir = NULL;
	_id_increment = 1;
	_lambda = 1.0;

	_widom = false;

	_localInsertionsMinusDeletions = 0;
}

void ChemicalPotential::setSubdomain(int rank, double x0, double x1, double y0,
		double y1, double z0, double z1)
{
	_ownrank = rank;
	_rnd.init(8624);
	_rndmomenta.init(8623);
	if (!_restrictedControlVolume) {
		_globalV = _system[0] * _system[1] * _system[2];

		for (int d = 0; d < 3; d++) {
			_control_bottom[d] = 0.0;
			_control_top[d] = _system[d];
		}
	}

	_minredco[0] = (x0 - _control_bottom[0])
			/ (_control_top[0] - _control_bottom[0]);
	_minredco[1] = (y0 - _control_bottom[1])
			/ (_control_top[1] - _control_bottom[1]);
	_minredco[2] = (z0 - _control_bottom[2])
			/ (_control_top[2] - _control_bottom[2]);
	_maxredco[0] = (x1 - _control_bottom[0])
			/ (_control_top[0] - _control_bottom[0]);
	_maxredco[1] = (y1 - _control_bottom[1])
			/ (_control_top[1] - _control_bottom[1]);
	_maxredco[2] = (z1 - _control_bottom[2])
			/ (_control_top[2] - _control_bottom[2]);
}

void ChemicalPotential::setSystem(double x, double y, double z, double m)
{
	_system[0] = x;
	_system[1] = y;
	_system[2] = z;
	_molecularMass = m;
	if (!_restrictedControlVolume) {
		_globalV = x * y * z;

		for (int d = 0; d < 3; d++) {
			_control_bottom[d] = 0.0;
			_control_top[d] = _system[d];
		}
	}
}

// note that *C must not contain the halo
// but when the decisions are evaluated, the halo must be taken into account!
//
void ChemicalPotential::prepareTimestep(ParticleContainer* moleculeContainer,
		DomainDecompBase* comm)
{
	_remainingDeletions.clear();
#ifndef NDEBUG
	for (int d = 0; d < 3; d++)
		mardyn_assert(_remainingInsertions[d].empty());
#endif
	_remainingDecisions.clear();

	// get information on the system decomposition
	//
	unsigned localN;
	if ((_maxredco[0] < 0.0) || (_maxredco[1] < 0.0) || (_maxredco[2] < 0.0)
			|| (_minredco[0] > 1.0) || (_minredco[1] > 1.0)
			|| (_minredco[2] > 1.0))
		localN = 0;
	else if ((_minredco[0] < 0.0) || (_minredco[1] < 0.0) || (_minredco[2] < 0.0)
			|| (_maxredco[0] > 1.0) || (_maxredco[1] > 1.0)
			|| (_maxredco[2] > 1.0))
		localN = this->countParticles(moleculeContainer, _componentid, _control_bottom,
				_control_top);
	else
		localN = this->countParticles(moleculeContainer, _componentid);
	float minrnd = 0.0;
	float maxrnd = 1.0;
	_globalN = comm->Ndistribution(localN, &minrnd, &maxrnd);
#ifndef NDEBUG
	Log::global_log->debug() << " believes N(" << _componentid << ")=" << _globalN
			<< ", rho=" << _globalN / _globalV
			<< ", the decisive density quotient equals "
			<< (float) _globalN / _globalReducedVolume << "\n";
#endif

	// construct deletions (disabled for Widom test particle method)
	//
	float sel, dec;
	unsigned localIndex;
	if (!_widom) {
		for (unsigned i = 0; i < _instances; i++) {
			sel = _rnd.rnd();
			dec = _rnd.rnd();
#ifndef NDEBUG
			// if(!ownrank) cout << "global index " << sel << " chosen for deletion.\n";
#endif
			if ((sel >= minrnd) && (sel < maxrnd)) {
				localIndex = (unsigned) floor(
						localN * (sel - minrnd) / (maxrnd - minrnd));
#ifndef NDEBUG
				// cout << "rank " << ownrank << " will try to delete index " << localIndex << ".\n";  // \\ //
#endif
				_remainingDeletions.push_back(localIndex);
				_remainingDecisions.push_back(dec);
			}
		}
	}

	int insertions = _instances;
#ifndef NDEBUG
	Log::global_log->debug() << "Number of insertions: " << insertions << ".\n";
#endif

	// construct insertions
	//
	float redc[3];
	double tc[3];
	for (int i = 0; i < insertions; i++) {
		for (int d = 0; d < 3; d++)
			redc[d] = _rnd.rnd();
		dec = _rnd.rnd();
		if ((redc[0] >= _minredco[0]) && (redc[1] >= _minredco[1])
				&& (redc[2] >= _minredco[2]) && (redc[0] < _maxredco[0])
				&& (redc[1] < _maxredco[1]) && (redc[2] < _maxredco[2])) {
			for (int d = 0; d < 3; d++) {
				tc[d] = _control_bottom[d]
						+ redc[d] * (_control_top[d] - _control_bottom[d]);
			}
#ifndef NDEBUG
			// cout << "rank " << ownrank << " will try to insert ID "
			//      << nextid << " (" << tc[0] << "/" << tc[1]
			//      << "/" << tc[2] << ").\n";  // \\ //
#endif
			if (moleculeContainer->isInBoundingBox(tc)) { // necessary because of rounding errors
				for (int d = 0; d < 3; d++) {
					_remainingInsertions[d].push_back(tc[d]);
				}
				_remainingDecisions.push_back(dec);
				_remainingInsertionIDs.push_back(_nextid);
			}
		}
		_nextid += _id_increment;
	}
}

bool ChemicalPotential::moleculeStrictlyNotInBox(const Molecule& m,
		const double l[3], const double u[3]) const
{
	bool out = false;
	for (int d = 0; d < 3; ++d) {
		out |= m.r(d) < l[d] or m.r(d) > u[d];
	}
	return out;
}

ParticleIterator ChemicalPotential::getDeletion(ParticleContainer* moleculeContainer, double* minco, double* maxco)
{
	if (_remainingDeletions.empty()) {
		ParticleIterator m;
		// a default constructed iterator is always invalid.
		return m; // DELETION_FALSE (always occurring for Widom)
	}
	unsigned idx = *_remainingDeletions.begin();
	_remainingDeletions.erase(_remainingDeletions.begin());
	double tminco[3];
	double tmaxco[3];
	if (_restrictedControlVolume)
		for (int d = 0; d < 3; d++) {
			tminco[d] = (minco[d] > _control_bottom[d]) ? minco[d] : _control_bottom[d];
			tmaxco[d] = (maxco[d] < _control_top[d])    ? maxco[d] : _control_top[d];
		}
	else
		for (int d = 0; d < 3; d++) {
			tminco[d] = minco[d];
			tmaxco[d] = maxco[d];
		}

	auto m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
	if (not m.isValid()){
		return m; // DELETION_INVALID
	}

	int j = 0;
	for (unsigned i = 0; (i < idx); i++) {
		while ((moleculeStrictlyNotInBox(*m, tminco, tmaxco) or (m->componentid() != _componentid))
				and m.isValid())
		{
			++m;
			if (not m.isValid()) {
				if (j == 0)
					return m; // DELETION_FALSE
				m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
				j = 0;
			}
		} /*end while*/

		++m;
		j++;
		if (not m.isValid()) {
			m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			j = 0;
		}
	}

	while (moleculeStrictlyNotInBox(*m, tminco, tmaxco) or (m->componentid() != _componentid))
	{
		++m;
		if (not m.isValid()) {
			if (j == 0)
				return m; // DELETION_FALSE
			m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
		}
	}

#ifndef NDEBUG
	Log::global_log->debug() << "ID " << m->getID() << " selected for deletion (index " << idx << ")." << std::endl;
#endif

	mardyn_assert(m->getID() < _nextid);
	return m; // DELETION_TRUE
}

// returns 0 if no insertion remains for this subdomain
unsigned long ChemicalPotential::getInsertion(double* ins)
{
	if (_remainingInsertionIDs.empty())
		return 0;

	for (int d = 0; d < 3; d++) {
		ins[d] = *_remainingInsertions[d].begin();
		_remainingInsertions[d].erase(
				_remainingInsertions[d].begin());
	}
	unsigned long nextid = *_remainingInsertionIDs.begin();
	_remainingInsertionIDs.erase(_remainingInsertionIDs.begin());
	return nextid;
}

bool ChemicalPotential::decideDeletion(double deltaUTilde)
{
	mardyn_assert(!_widom); // the Widom test particle method should never call decideDeletion ...

	if (_remainingDecisions.empty()) {
		if (_widom) {
			Log::global_log->error()
					<< "SEVERE WARNING: The Widom method is (erroneously) trying to carry out test deletions.\n";
			return false;
		}
		std::ostringstream error_message;
		error_message << "No decision is possible." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	float dec = *_remainingDecisions.begin();
	_remainingDecisions.erase(_remainingDecisions.begin());
	float acc = ((float) (_globalN)) * exp(_muTilde - deltaUTilde)
			/ _globalReducedVolume;
	bool ans;
	if (dec < 0.000001)
		ans = true;
	else if (dec > 0.999999)
		ans = false;
	else
		ans = (acc > dec);
	// ans = (acc > 1.0e-05)? (acc > dec): false;
#ifndef NDEBUG
	// cout << "rank " << ownrank << (ans? " accepted ": " rejected ")
	//      << "deletion with deltaUtilde = " << deltaUTilde << " (P = "
	//      << ((acc > 1.0)? 1.0: acc) << ").\n"; // \\ //
#endif
	if (ans)
		_globalN -= (2 * _ownrank + 1); // estimate, the precise value is communicated later
	return ans;
}

bool ChemicalPotential::decideInsertion(double deltaUTilde)
{
	if (_remainingDecisions.empty()) {
		if (_widom) {
			Log::global_log->error() << "!!! SEVERE WARNING on rank " << _ownrank
					<< ": no decision is possible !!!\n";
			return false;
		}
		std::ostringstream error_message;
		error_message << "No decision is possible." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	double acc = _globalReducedVolume * exp(_muTilde - deltaUTilde)
			/ (1.0 + (double) (_globalN));

	bool ans;
	if (_widom)
		ans = false; // the Widom method does not actually insert any particles ...
	else {
		float dec = *_remainingDecisions.begin();
		if (dec < 0.000001)
			ans = true;
		else if (dec > 0.999999)
			ans = false;
		else
			ans = (acc > dec);
		// ans = (acc > 1.0e-05)? (acc > dec): false;
#ifndef NDEBUG
		// cout << "rank " << ownrank << (ans? " accepted ": " rejected ")
		//      << "insertion with deltaUtilde = " << deltaUTilde
		//      << " (P = " << ((acc > 1.0)? 1.0: acc) << ").\n"; // \\ //
#endif
		if (ans)
			_globalN += (2 * _ownrank + 1); // estimate, the precise value is communicated later
	}
	_remainingDecisions.erase(_remainingDecisions.begin());
	return ans;
}

void ChemicalPotential::submitTemperature(double T_in)
{
	_T = T_in;
	_muTilde = _mu / _T;
	_lambda = 0.39894228 * _h / sqrt(_molecularMass * _T);
	_globalReducedVolume = _globalV / (_lambda * _lambda * _lambda);
	_decisive_density = (float) _globalN / _globalReducedVolume;
	double doOutput = _rnd.rnd();
#ifdef NDEBUG
	if(_ownrank) return;
#endif
	if (doOutput >= 0.01)
		return;
	std::cout << "rank " << _ownrank << " sets mu~ <- " << _muTilde;
	std::cout << ", T <- " << _T << ", lambda <- " << _lambda;
	std::cout << ", and Vred <- " << _globalReducedVolume << "\n";
}

void ChemicalPotential::assertSynchronization(DomainDecompBase* comm) {
	comm->assertIntIdentity(_rnd.getIX());
}

void ChemicalPotential::setControlVolume(double x0, double y0, double z0,
		double x1, double y1, double z1)
{
	if ((x0 >= x1) || (y0 >= y1) || (z0 >= z1)) {
		std::ostringstream error_message;
		error_message << "\nInvalid control volume (" << x0 << " / " << y0
				<< " / " << z0 << ") to (" << x1 << " / " << y1 << " / " << z1
				<< ")." << std::endl;
		MARDYN_EXIT(error_message.str());
	}
	_restrictedControlVolume = true;
	_globalV = (x1 - x0) * (y1 - y0) * (z1 - z0);
	_control_bottom[0] = x0;
	_control_top[0] = x1;
	_control_bottom[1] = y0;
	_control_top[1] = y1;
	_control_bottom[2] = z0;
	_control_top[2] = z1;
}

Molecule ChemicalPotential::loadMolecule()
{
	mardyn_assert(_reservoir != NULL);
	Molecule tmp = *_reservoir;
	unsigned rotdof = tmp.component()->getRotationalDegreesOfFreedom();
	mardyn_assert(tmp.componentid() == _componentid);
#ifndef NDEBUG
	tmp.check(tmp.getID());
#endif
	if (!_widom) {
		double v[3];
		double vv = 0.0;
		for (int d = 0; d < 3; d++) {
			v[d] = -0.5 + _rndmomenta.rnd();
			vv += v[d] * v[d];
		}
		double vnorm = sqrt(3.0 * _T / (vv * tmp.mass()));
		for (int d = 0; d < 3; d++)
			tmp.setv(d, v[d] * vnorm);

		if (rotdof > 0) {
			double qtr[4];
			double qqtr = 0.0;
			for (int d = 0; d < 4; d++) {
				qtr[d] = -0.5 + _rndmomenta.rnd();
				qqtr += qtr[d] * qtr[d];
			}
			double qtrnorm = sqrt(1.0 / qqtr);
			Quaternion tqtr = Quaternion(qtr[0] * qtrnorm, qtr[1] * qtrnorm,
					qtr[2] * qtrnorm, qtr[3] * qtrnorm);
			tmp.setq(tqtr);

			std::array<double,3> D;
			double Dnorm = 0.0;
			for (int d = 0; d < 3; d++)
				D[d] = -0.5 + _rndmomenta.rnd();
			std::array<double, 3> w = tqtr.rotateinv(D);
			double Iw2 = w[0] * w[0] * tmp.component()->I11()
					+ w[1] * w[1] * tmp.component()->I22()
					+ w[2] * w[2] * tmp.component()->I33();
			Dnorm = sqrt(_T * rotdof / Iw2);
			for (int d = 0; d < 3; d++)
				tmp.setD(d, D[d] * Dnorm);
		}
	}

	return tmp;
}

int ChemicalPotential::grandcanonicalBalance(DomainDecompBase* comm) {
	auto collComm = makeCollCommObjAllreduceAdd(comm->getCommunicator(), _localInsertionsMinusDeletions);
	collComm.communicate();
	auto [universalInsertionsMinusDeletions] = collComm.get();

	return universalInsertionsMinusDeletions;
}

void ChemicalPotential::grandcanonicalStep(
		ParticleContainer* moleculeContainer, double T, Domain* domain,
		CellProcessor* cellProcessor)
{
	bool accept = true;
	double DeltaUpot;

	ParticlePairs2PotForceAdapter particlePairsHandler(*domain);

	_localInsertionsMinusDeletions = 0;

	this->submitTemperature(T);
	double minco[3];
	double maxco[3];
	for (int d = 0; d < 3; d++) {
		minco[d] = moleculeContainer->getBoundingBoxMin(d);
		maxco[d] = moleculeContainer->getBoundingBoxMax(d);
	}

	bool hasDeletion = true;
	bool hasInsertion = true;
	double ins[3];
	unsigned nextid = 0;
	while (hasDeletion || hasInsertion) {
		if (hasDeletion) {
			auto m = this->getDeletion(moleculeContainer, minco, maxco);
			if(m.isValid()) {
				DeltaUpot = -1.0 * moleculeContainer->getEnergy(&particlePairsHandler, &(*m), *cellProcessor);

				accept = this->decideDeletion(DeltaUpot / T);
#ifndef NDEBUG
				if (accept) {
					std::cout << "r" << this->rank() << "d" << m->getID() << " with energy " << DeltaUpot << std::endl;
					std::cout.flush();
				}
				/*
				 else
				 cout << "   (r" << this->rank() << "-d" << m->getID()
				 << ")" << std::endl;
				 */
#endif
				if (accept) {
					// m->upd_cache(); TODO what to do here? somebody deleted the method "upd_cache"!!! why???
					// reset forces and momenta to zero
					{
						double zeroVec[3] = {0.0, 0.0, 0.0};
						m->setF(zeroVec);
						m->setM(zeroVec);
						m->setVi(zeroVec);
					}

					this->storeMolecule(*m);

					moleculeContainer->deleteMolecule(m, true/*rebuildCaches*/);
					m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);
					_localInsertionsMinusDeletions--;
				}
			} else{
				hasDeletion = false;
			}
		} /* end of second hasDeletion */

		if (!this->hasSample()) {
			for (auto mit = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); mit.isValid(); ++mit) {
				if (mit->componentid() == this->getComponentID()) {

					this->storeMolecule(*mit);
					break;
				}
			}
		}
		if (hasInsertion) {
			nextid = this->getInsertion(ins);
			hasInsertion = (nextid > 0);
		}
		if (hasInsertion) {
			Molecule tmp = this->loadMolecule();
			for (int d = 0; d < 3; d++)
				tmp.setr(d, ins[d]);
			tmp.setid(nextid);
			// TODO: I (tchipevn) am changing the former code
			// (add particle, compute energy and if rejected, delete particle)
			// to: add particle only if accepted
//			_particles.push_back(tmp);

//			std::list<Molecule>::iterator mit = _particles.end();
			Molecule* mit = &tmp;
			// m->upd_cache(); TODO: what to do here? somebody deleted the method "upd_cache"!!! why??? -- update caches do no longer exist ;)
			// reset forces and torques to zero
			if (!this->isWidom()) {
				double zeroVec[3] = { 0.0, 0.0, 0.0 };
				mit->setF(zeroVec);
				mit->setM(zeroVec);
				mit->setVi(zeroVec);
			}
			mit->check(nextid);
#ifndef NDEBUG
			/*
			 cout << "rank " << this->rank() << ": insert "
			 << m->getID() << " at the reduced position (" << ins[0] << "/"
			 << ins[1] << "/" << ins[2] << ")? " << std::endl;
			 */
#endif

//			unsigned long cellid = moleculeContainer->getCellIndexOfMolecule(m);
//			moleculeContainer->_cells[cellid].addParticle(m);
			DeltaUpot = moleculeContainer->getEnergy(&particlePairsHandler, mit,
					*cellProcessor);
			domain->submitDU(this->getComponentID(), DeltaUpot, ins);
			accept = this->decideInsertion(DeltaUpot / T);

#ifndef NDEBUG
			if (accept) {
				std::cout << "r" << this->rank() << "i" << mit->getID()
						<< " with energy " << DeltaUpot << std::endl;
				std::cout.flush();
			}
			/*
			 else
			 cout << "   (r" << this->rank() << "-i"
			 << mit->getID() << ")" << std::endl;
			 */
#endif
			if (accept) {
				this->_localInsertionsMinusDeletions++;
				double zeroVec[3] = { 0.0, 0.0, 0.0 };
				tmp.setVi(zeroVec);
				bool inBoxCheckedAlready = false, checkWhetherDuplicate = false, rebuildCaches = true;
				moleculeContainer->addParticle(tmp, inBoxCheckedAlready, checkWhetherDuplicate, rebuildCaches);
			} else {
				// moleculeContainer->deleteMolecule(m->getID(), m->r(0), m->r(1), m->r(2));
//				moleculeContainer->_cells[cellid].deleteMolecule(m->getID());
				double zeroVec[3] = { 0.0, 0.0, 0.0 };
				mit->setF(zeroVec);
				mit->setM(zeroVec);
				mit->setVi(zeroVec);
				mit->check(mit->getID());
//				moleculeContainer->_particles.erase(mit);
			}
		}
	}
#ifndef NDEBUG
	for (auto m = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY); m.isValid(); ++m) {
		// cout << *m << "\n";
		// cout.flush();
		m->check(m->getID());
	}
#endif
}

unsigned ChemicalPotential::countParticles(
		ParticleContainer* moleculeContainer, unsigned int cid) const
{
	// ParticleContainer::countParticles functionality moved here as it was:
	// i.e. halo-particles are NOT counted

	unsigned N = 0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:N)
	#endif
	{

		auto begin = moleculeContainer->iterator(ParticleIterator::ONLY_INNER_AND_BOUNDARY);

		for (auto m = begin; m.isValid(); ++m) {
			if (m->componentid() == cid)
				++N;
		}

	}

	return N;
}

unsigned ChemicalPotential::countParticles(
		ParticleContainer* moleculeContainer, unsigned int cid,
		double* cbottom, double* ctop) const
{
	// ParticleContainer::countParticles functionality moved here as it was:
	// i.e. halo-particles NOT counted

	unsigned N = 0;

	#if defined(_OPENMP)
	#pragma omp parallel reduction(+:N)
	#endif
	{
		for (auto m = moleculeContainer->regionIterator(cbottom, ctop, ParticleIterator::ONLY_INNER_AND_BOUNDARY);
			 m.isValid(); ++m) {
			if (m->componentid() == cid) {
				++N;
			}
		}
	}

	return N;
}
