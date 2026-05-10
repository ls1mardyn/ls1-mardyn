/**
 * @file NeighborAcquirer.cpp
 * @author seckler
 * @date 06.05.19
 */

#include <tuple>
#include "NeighborAcquirer.h"
#include "Domain.h"
#include "HaloRegion.h"

/*
 * 1. Initial Exchange of all desired regions.
 * 2. Each process checks whether he owns parts of the desired regions and will save those regions in partners02.
 * 3. Each process will notify the other processes whether they own parts of their desired domains (i.e. whether they
 * are a partner for them)
 * 4. The processes talk with each other to specify the exact domains they will communicate. Received parts will be
 * saved in partners01.
 */
std::tuple<std::vector<CommunicationPartner>, std::vector<CommunicationPartner>> NeighborAcquirer::acquireNeighbors(
	const std::array<double, 3> &globalDomainLength, HaloRegion *ownRegion, const std::vector<HaloRegion> &desiredRegions,
	const MPI_Comm &comm, bool excludeOwnRank) {

	int my_rank{};  // my rank
	MPI_Comm_rank(comm, &my_rank);
	int num_processes{};  // the number of processes in comm
	MPI_Comm_size(comm, &num_processes);

	const auto num_regions = desiredRegions.size();  // the number of regions I would like to acquire from other processes

	// tell the other processes how much you are going to send
	// how many bytes am I going to send to all the other processes
	const int num_bytes_send =
		sizeof(int) * 2 + (sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) * 1) * num_regions;

	// create byte send buffer
	std::vector<unsigned char> outgoingDesiredRegionsVector(num_bytes_send);  // outgoing byte buffer

	// msg format: rank | number_of_regions | region_01 | region_02 | ...
	// fill the buffer
	int bufferPosition = 0;
	memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, &my_rank, sizeof(int));
	bufferPosition += sizeof(int);
	memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, &num_regions, sizeof(int));
	bufferPosition += sizeof(int);

	for (auto &region : desiredRegions) {  // filling up the outgoing byte buffer
		memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, region.rmin, sizeof(double) * 3);
		bufferPosition += sizeof(double) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, region.rmax, sizeof(double) * 3);
		bufferPosition += sizeof(double) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, region.offset, sizeof(int) * 3);
		bufferPosition += sizeof(int) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + bufferPosition, &region.width, sizeof(double));
		bufferPosition += sizeof(double);
	}

	// set up structure information data for the Allgatherv operation
	// vector of number of bytes I am going to receive
	std::vector<int> num_bytes_receive_vec(num_processes, 0);
	MPI_Allgather(&num_bytes_send, 1, MPI_INT, num_bytes_receive_vec.data(), 1, MPI_INT, comm);
	// vector of offsets (=displacement in MPI) in the receive buffer
	std::vector<int> num_bytes_displacements(num_processes, 0);
	int num_bytes_receive = 0;
	for (int j = 0; j < num_processes; j++) {
		num_bytes_displacements[j] = num_bytes_receive;
		num_bytes_receive += num_bytes_receive_vec[j];
	}

	std::vector<unsigned char> incomingDesiredRegionsVector(num_bytes_receive);  // the incoming byte buffer

	// send your regions
	MPI_Allgatherv(outgoingDesiredRegionsVector.data(), num_bytes_send, MPI_BYTE, incomingDesiredRegionsVector.data(),
				   num_bytes_receive_vec.data(), num_bytes_displacements.data(), MPI_BYTE, comm);

	std::vector<int> numberOfRegionsToSendToRank(num_processes, 0);       // outgoing row

	// parse / deserialize received data
	constexpr int bytesOneRegion =
		sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) + sizeof(double) * 3;
	// the regions I own and want to send: ranks<regions<regionData>>
	std::vector<std::vector<std::vector<unsigned char>>> sendingList(num_processes);

	std::vector<CommunicationPartner> comm_partners02{};
	bufferPosition = 0;
	while (bufferPosition < num_bytes_receive /*== buffer length*/) {

		int rank{};
		memcpy(&rank, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(int));
		bufferPosition += sizeof(int);  // 4
		int regions{};
		memcpy(&regions, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(int));
		bufferPosition += sizeof(int);  // 4

		for (int regionId = 0; regionId < regions; ++regionId) {
			HaloRegion unshiftedRegion{};
			memcpy(unshiftedRegion.rmin, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(double) * 3);
			bufferPosition += sizeof(double) * 3;  // 24
			memcpy(unshiftedRegion.rmax, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(double) * 3);
			bufferPosition += sizeof(double) * 3;  // 24
			memcpy(unshiftedRegion.offset, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(int) * 3);
			bufferPosition += sizeof(int) * 3;  // 12
			memcpy(&unshiftedRegion.width, incomingDesiredRegionsVector.data() + bufferPosition, sizeof(double));
			bufferPosition += sizeof(double);  // 4

			// msg format one region: rmin | rmax | offset | width | shift
			const auto [regionsToTest, shifts] = getPotentiallyShiftedRegions(globalDomainLength, unshiftedRegion);
			// Before every set of push_backs make sure there is enough space for this set + all remaining.
			// This guarantees that there is enough space for the current set of push_backs, and, if subsequent sets
			// are smaller, further reallocations can be avoided. This potentially leads to an overestimate but comes
			// with the advantage of fewer resizes.
			sendingList[rank].reserve(sendingList[rank].size() + ((regions - regionId) * regionsToTest.size()));
			comm_partners02.reserve(comm_partners02.size() + ((regions - regionId) * regionsToTest.size()));
			for(size_t regionIndex = 0; regionIndex < regionsToTest.size(); ++regionIndex){
				auto regionToTest = regionsToTest[regionIndex];
				if ((not excludeOwnRank or rank != my_rank) and isIncluded(ownRegion, &regionToTest)) {
					auto currentShift = shifts[regionIndex];

					numberOfRegionsToSendToRank[rank]++;  // this is a region I will send to rank

					const auto overlappedRegion = overlap(*ownRegion, regionToTest);  // different shift for the overlap?

					// make a note in partners02 - don't forget to squeeze partners02
					constexpr bool enlarged[3][2] = {{false}};
					for (int k = 0; k < 3; k++) {
						currentShift[k] *= -1;
					}

					comm_partners02.emplace_back(rank, overlappedRegion.rmin, overlappedRegion.rmax, overlappedRegion.rmin,
												 overlappedRegion.rmax, currentShift.data(), overlappedRegion.offset, enlarged);

					for (int k = 0; k < 3; k++) {
						currentShift[k] *= -1;
					}

					// Undo the shift. So it is again in the perspective of the rank we got this region from.
					// We cannot use unshiftedRegion, as it is not overlapped and thus potentially too big.
					HaloRegion unshiftedOverlappedRegion{overlappedRegion};
					for (int dimI = 0; dimI < 3; ++dimI) {
						unshiftedOverlappedRegion.rmax[dimI] -= currentShift[dimI];
						if (fabs(unshiftedOverlappedRegion.rmax[dimI]) < 1e-10 and currentShift[dimI] != 0.) {
							// we have to ensure that if we shifted, then the rmax, etc. are correct!
							unshiftedOverlappedRegion.rmax[dimI] = 0.;
						}
						unshiftedOverlappedRegion.rmin[dimI] -= currentShift[dimI];
						if (fabs(unshiftedOverlappedRegion.rmin[dimI] - globalDomainLength[dimI]) < 1e-10 and
							currentShift[dimI] != 0.) {
							// we have to ensure that if we shifted, then the rmax, etc. are correct!
							unshiftedOverlappedRegion.rmin[dimI] = globalDomainLength[dimI];
						}
					}

					std::vector<unsigned char> singleRegion(bytesOneRegion);

					int p = 0;
					memcpy(&singleRegion[p], unshiftedOverlappedRegion.rmin, sizeof(double) * 3);
					p += sizeof(double) * 3;
					memcpy(&singleRegion[p], unshiftedOverlappedRegion.rmax, sizeof(double) * 3);
					p += sizeof(double) * 3;
					memcpy(&singleRegion[p], unshiftedOverlappedRegion.offset, sizeof(int) * 3);
					p += sizeof(int) * 3;
					memcpy(&singleRegion[p], &unshiftedOverlappedRegion.width, sizeof(double));
					p += sizeof(double);
					memcpy(&singleRegion[p], currentShift.data(), sizeof(double) * 3);
					//p += sizeof(double) * 3;

					sendingList[rank].emplace_back(std::move(singleRegion));
				}
			}
		}
	}

	std::vector<std::vector<unsigned char>> merged(num_processes);  // Merge each list of char arrays into one char array
	for (int j = 0; j < num_processes; j++) {
		if (numberOfRegionsToSendToRank[j] > 0) {
			std::vector<unsigned char> mergedRegions(numberOfRegionsToSendToRank[j] * bytesOneRegion);

			for (int k = 0; k < numberOfRegionsToSendToRank[j]; k++) {
				memcpy(&mergedRegions[k * bytesOneRegion], sendingList[j][k].data(), bytesOneRegion);
			}

			merged[j] = std::move(mergedRegions);
		}
	}

	// We cannot know how many regions we are going to receive from each process a-priori.
	// So we need to figure this out.
	/* e.g.: 4x4
	 *
	 * row := sender
	 * column := receiver
	 *
	 *    0 1 2 3 4
	 *   -----------
	 * 0 | | |3| | |
	 *   -----------
	 * 1 | | | |2| |
	 *   -----------
	 * 2 | | | | | |
	 *   -----------
	 * 3 | | | | | |
	 *   -----------
	 * reduce
	 *   | | |2|2| |
	 *
	 * Each process has a horizontal vector, where it marks how many regions it is going to send to another process.
	 * In this case, process 0 will send 3 regions to process 2 and process 1 will send 2 regions to process 3.
	 * After the Allreduce step every process has the information how many regions it will receive.
	 */
	std::vector<int> numberOfRegionsToReceive(num_processes, 0);  // how many bytes does each process expect?
	MPI_Allreduce(numberOfRegionsToSendToRank.data(), numberOfRegionsToReceive.data(), num_processes, MPI_INT, MPI_SUM, comm);

	// all the information for the final information exchange has been collected -> final exchange

	std::vector<MPI_Request> requests(num_processes, MPI_REQUEST_NULL);
	MPI_Status probe_status;
	MPI_Status rec_status;

	// sending (non blocking)
	for (int j = 0; j < num_processes; j++) {
		if (numberOfRegionsToSendToRank[j] > 0) {
			MPI_Isend(merged[j].data(), numberOfRegionsToSendToRank[j] * bytesOneRegion, MPI_BYTE, j, 1, comm,
					  &requests[j]);  // tag is one
		}
	}

	std::vector<CommunicationPartner> comm_partners01;  // the communication partners

	// receive data (blocking)
	/**
	 * We now receive as many regions as we previously determined that we will receive.
	 * For that we keep track, how many regions we received and increase this according to the number of regions
	 * received per MPI operation.
	 */
	for (int byte_counter = 0; byte_counter < numberOfRegionsToReceive[my_rank] * bytesOneRegion; ) {
		// MPI_PROBE
		MPI_Probe(MPI_ANY_SOURCE, 1, comm, &probe_status);
		// interpret probe
		const auto source = probe_status.MPI_SOURCE;
		int bytes{};
		MPI_Get_count(&probe_status, MPI_BYTE, &bytes);
		// we have receive `bytes` bytes. So we increase the byte_counter.
		byte_counter += bytes;
		// create buffer
		std::vector<unsigned char> raw_neighbours(bytes);
		MPI_Recv(raw_neighbours.data(), bytes, MPI_BYTE, source, 1, comm, &rec_status);
		// Interpret Buffer and add neighbours
		const auto numRegionsToReceive = bytes / bytesOneRegion;
		comm_partners01.reserve(std::max(comm_partners01.size(), static_cast<size_t>(numberOfRegionsToReceive[my_rank] * numRegionsToReceive)));
		for (int regionId = 0; regionId < numRegionsToReceive; ++regionId) {  // number of regions from this process
			HaloRegion region{};
			bufferPosition = regionId * bytesOneRegion;

			memcpy(region.rmin, raw_neighbours.data() + bufferPosition, sizeof(double) * 3);
			bufferPosition += sizeof(double) * 3;
			memcpy(region.rmax, raw_neighbours.data() + bufferPosition, sizeof(double) * 3);
			bufferPosition += sizeof(double) * 3;
			memcpy(region.offset, raw_neighbours.data() + bufferPosition, sizeof(int) * 3);
			bufferPosition += sizeof(int) * 3;
			memcpy(&region.width, raw_neighbours.data() + bufferPosition, sizeof(double));
			bufferPosition += sizeof(double);

			double shift[3];
			memcpy(shift, raw_neighbours.data() + bufferPosition, sizeof(double) * 3);
			// bufferPosition += sizeof(double) * 3;

			constexpr bool enlarged[3][2] = {{false}};

			comm_partners01.emplace_back(source, region.rmin, region.rmax, region.rmin, region.rmax, shift,
										 region.offset, enlarged);
		}
	}

	// ensure that all sends have been finished.
	for (int j = 0; j < num_processes; j++) {
		if (numberOfRegionsToSendToRank[j] > 0) MPI_Wait(&requests[j], MPI_STATUS_IGNORE);
	}

	// barrier for safety.
	MPI_Barrier(comm);

	return std::make_tuple(squeezePartners(comm_partners01), squeezePartners(comm_partners02));
}

std::vector<CommunicationPartner> NeighborAcquirer::squeezePartners(const std::vector<CommunicationPartner> &partners) {
	std::vector<CommunicationPartner> squeezedPartners;
	std::vector<bool> used(partners.size(),
						   false);  // flag table, that describes, whether a certain comm-partner has already been added
	for (unsigned int i = 0; i < partners.size(); i++) {
		if (used[i]) continue;  // if we already added the neighbour, don't add it again!
		int rank = partners[i].getRank();
		CommunicationPartner tmp = partners[i];
		for (unsigned int j = i + 1; j < partners.size(); j++) {
			if (partners[j].getRank() != rank) continue;  // only add those with same rank
			tmp.add(partners[j]);
			used[j] = true;
		}
		squeezedPartners.push_back(tmp);
	}
	return squeezedPartners;
}

bool NeighborAcquirer::isIncluded(HaloRegion *myRegion, HaloRegion *inQuestion) {
	return myRegion->rmax[0] > inQuestion->rmin[0] && myRegion->rmin[0] < inQuestion->rmax[0] &&
		   myRegion->rmax[1] > inQuestion->rmin[1] && myRegion->rmin[1] < inQuestion->rmax[1] &&
		   myRegion->rmax[2] > inQuestion->rmin[2] && myRegion->rmin[2] < inQuestion->rmax[2];
	// myRegion->rmax > inQuestion->rmin
	// && myRegion->rmin < inQuestion->rmax
}

HaloRegion NeighborAcquirer::overlap(const HaloRegion& myRegion, const HaloRegion& inQuestion) {
	/*
	 * Choose the overlap of myRegion and inQuestion.
	 */
	HaloRegion overlap{inQuestion};

	for (int i = 0; i < 3; i++) {
		overlap.rmax[i] = std::min(myRegion.rmax[i], inQuestion.rmax[i]);
		overlap.rmin[i] = std::max(myRegion.rmin[i], inQuestion.rmin[i]);
	}

	return overlap;
}

std::pair<std::vector<HaloRegion>, std::vector<std::array<double, 3>>> NeighborAcquirer::getPotentiallyShiftedRegions(
				  const std::array<double, 3> &domainLength, const HaloRegion &region) {
	std::vector<HaloRegion> haloRegions;
	std::vector<std::array<double, 3>> shifts;

	std::array<std::vector<int>,3> doShiftsVector;
	for (unsigned dim = 0; dim < 3; ++dim) {
		// if rmin is small enough, include wrapping over bottom of domain -> positive shift
		if (region.rmin[dim] < 0) {
			doShiftsVector[dim].emplace_back(1);
		}
		// if rmax is big enough, include wrapping over top of domain -> negative shift
		if (region.rmax[dim] > domainLength[dim]) {
			doShiftsVector[dim].emplace_back(-1);
		}
		// if halo region is not completely outside, include non-wrapped halo -> zero shift
		if (region.rmax[dim] > 0 and region.rmin[dim] < domainLength[dim]) {
			doShiftsVector[dim].emplace_back(0);
		}
		// the shift vector should never be empty!
		mardyn_assert(not doShiftsVector[dim].empty());
	}

	auto num_regions = doShiftsVector[0].size() * doShiftsVector[1].size() * doShiftsVector[2].size();
	haloRegions.reserve(num_regions);
	shifts.reserve(num_regions);

	// Calculate and apply the shifts.
	for (int x_index_shift : doShiftsVector[0]) {
		for (int y_index_shift : doShiftsVector[1]) {
			for (int z_index_shift : doShiftsVector[2]) {
				std::array<int, 3> indexShifts{x_index_shift, y_index_shift, z_index_shift};
				std::array<double, 3> shift{};
				auto shiftedRegion = region;
				for(unsigned dim = 0; dim < 3; ++dim){
					shift[dim] = domainLength[dim] * indexShifts[dim];
					shiftedRegion.rmin[dim] += shift[dim];
					shiftedRegion.rmax[dim] += shift[dim];
				}
				haloRegions.emplace_back(shiftedRegion);
				shifts.emplace_back(shift);
			}
		}
	}
	return std::make_pair(haloRegions, shifts);
}
