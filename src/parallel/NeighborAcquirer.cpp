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
	const std::array<double,3>& globalDomainLength, HaloRegion *ownRegion, std::vector<HaloRegion> &desiredRegions, double skin, const MPI_Comm& comm) {

	HaloRegion ownRegionEnlargedBySkin = *ownRegion;
	for(unsigned int dim = 0; dim < 3; ++dim){
		ownRegionEnlargedBySkin.rmin[dim] -= skin;
		ownRegionEnlargedBySkin.rmax[dim] += skin;
	}

	int my_rank;  // my rank
	MPI_Comm_rank(comm, &my_rank);
	int num_processes;  // the number of processes in comm
	MPI_Comm_size(comm, &num_processes);

	int num_regions = desiredRegions.size();  // the number of regions I would like to acquire from other processes

	// tell the other processes how much you are going to send
	int num_bytes_send =
		sizeof(int) * 2 + (sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) * 1) *
							  num_regions;  // how many bytes am I going to send to all the other processes?
	std::vector<int> num_bytes_receive_vec(num_processes, 0);  // vector of number of bytes I am going to receive
	// MPI_Allreduce(&num_bytes_send, &num_bytes_receive, 1, MPI_INT, MPI_SUM, comm);
	MPI_Allgather(&num_bytes_send, 1, MPI_INT, num_bytes_receive_vec.data(), 1, MPI_INT, comm);

	// create byte buffer
	std::vector<unsigned char> outgoingDesiredRegionsVector(num_bytes_send);  // outgoing byte buffer
	int i = 0;
	int p = 0;

	// msg format: rank | number_of_regions | region_01 | region_02 | ...

	memcpy(outgoingDesiredRegionsVector.data() + i, &my_rank, sizeof(int));
	i += sizeof(int);
	memcpy(outgoingDesiredRegionsVector.data() + i, &num_regions, sizeof(int));
	i += sizeof(int);

	for (auto &region : desiredRegions) {  // filling up the outgoing byte buffer
		memcpy(outgoingDesiredRegionsVector.data() + i, region.rmin, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + i, region.rmax, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + i, region.offset, sizeof(int) * 3);
		i += sizeof(int) * 3;
		memcpy(outgoingDesiredRegionsVector.data() + i, &region.width, sizeof(double));
		i += sizeof(double);
	}

	int num_bytes_receive = 0;
	std::vector<int> num_bytes_displacements(num_processes, 0);  // vector of number of bytes I am going to receive
	for (int j = 0; j < num_processes; j++) {
		num_bytes_displacements[j] = num_bytes_receive;
		num_bytes_receive += num_bytes_receive_vec[j];
	}

	std::vector<unsigned char> incomingDesiredRegionsVector(num_bytes_receive);  // the incoming byte buffer

	// send your regions
	MPI_Allgatherv(outgoingDesiredRegionsVector.data(), num_bytes_send, MPI_BYTE, incomingDesiredRegionsVector.data(),
				   num_bytes_receive_vec.data(), num_bytes_displacements.data(), MPI_BYTE, comm);

	std::vector<int> numberOfRegionsToSendToRank(num_processes, 0);       // outgoing row

	int bytesOneRegion =
		sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) + sizeof(double) * 3;
	std::vector<std::vector<std::vector<unsigned char>>> sendingList(num_processes);  // the regions I own and want to send
	std::vector<CommunicationPartner> comm_partners02;

	i = 0;
	while (i != num_bytes_receive) {
		int rank;
		int regions;

		memcpy(&rank, incomingDesiredRegionsVector.data() + i, sizeof(int));
		i += sizeof(int);  // 4
		memcpy(&regions, incomingDesiredRegionsVector.data() + i, sizeof(int));
		i += sizeof(int);  // 4

		for (int j = 0; j < regions; j++) {
			HaloRegion unshiftedRegion{};
			memcpy(unshiftedRegion.rmin, incomingDesiredRegionsVector.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;  // 24
			memcpy(unshiftedRegion.rmax, incomingDesiredRegionsVector.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;  // 24
			memcpy(unshiftedRegion.offset, incomingDesiredRegionsVector.data() + i, sizeof(int) * 3);
			i += sizeof(int) * 3;  // 12
			memcpy(&unshiftedRegion.width, incomingDesiredRegionsVector.data() + i, sizeof(double));
			i += sizeof(double);  // 4

			// msg format one region: rmin | rmax | offset | width | shift
			std::array<double, 3> shift{};
			auto shiftedRegion = getPotentiallyShiftedRegion(globalDomainLength, unshiftedRegion, shift.data(), skin);

			std::vector<HaloRegion> regionsToTest;
			std::vector<std::array<double, 3>> shifts;
			if(skin == 0.){
				regionsToTest = {shiftedRegion};
				shifts = {shift};
			} else {
				// if the skin is not zero, the neighbor which owns a particle is no longer uniquely defined.
				// We thus have to check multiple different possible shifts.
				std::tie(regionsToTest, shifts) =
					getAllShiftedAndNonShiftedRegionsAndShifts(unshiftedRegion, shiftedRegion, shift);
			}

			for(size_t regionIndex = 0; regionIndex < regionsToTest.size(); ++regionIndex){
				auto regionToTest = regionsToTest[regionIndex];
				if (rank != my_rank && isIncluded(&ownRegionEnlargedBySkin, &regionToTest)) {
					auto currentShift = shifts[regionIndex];

					numberOfRegionsToSendToRank[rank]++;  // this is a region I will send to rank

					auto overlappedRegion = overlap(ownRegionEnlargedBySkin, regionToTest);  // different shift for the overlap?

					// make a note in partners02 - don't forget to squeeze partners02
					bool enlarged[3][2] = {{false}};
					for (int k = 0; k < 3; k++) currentShift[k] *= -1;

					if (skin != 0.) {
						for (size_t dim = 0; dim < 3; ++dim) {
							if (overlappedRegion.offset[dim] == -1 and overlappedRegion.rmax[dim] == ownRegion->rmax[dim]) {
								overlappedRegion.rmax[dim] = ownRegionEnlargedBySkin.rmax[dim];
							} else if (overlappedRegion.offset[dim] == 1 and overlappedRegion.rmin[dim] == ownRegion->rmin[dim]) {
								overlappedRegion.rmin[dim] = ownRegionEnlargedBySkin.rmin[dim];
							}
						}
					}

					comm_partners02.emplace_back(rank, overlappedRegion.rmin, overlappedRegion.rmax, overlappedRegion.rmin,
												 overlappedRegion.rmax, currentShift.data(), overlappedRegion.offset, enlarged);

					for (int k = 0; k < 3; k++) currentShift[k] *= -1;

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

					p = 0;
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

					sendingList[rank].push_back(std::move(singleRegion));
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
	int byte_counter = 0;

	/**
	 * We now receive as many regions as we previously determined that we will receive.
	 * For that we keep track, how many regions we received and increase this according to the number of regions
	 * received per MPI operation.
	 */
	while (byte_counter < numberOfRegionsToReceive[my_rank] * bytesOneRegion) {
		// MPI_PROBE
		MPI_Probe(MPI_ANY_SOURCE, 1, comm, &probe_status);
		// interpret probe
		int source = probe_status.MPI_SOURCE;
		int bytes;
		MPI_Get_count(&probe_status, MPI_BYTE, &bytes);
		// we have receive `bytes` bytes. So we increase the byte_counter.
		byte_counter += bytes;
		// create buffer
		std::vector<unsigned char> raw_neighbours(bytes);
		MPI_Recv(raw_neighbours.data(), bytes, MPI_BYTE, source, 1, comm, &rec_status);
		// Interpret Buffer and add neighbours
		for (int k = 0; k < (bytes / bytesOneRegion); k++) {  // number of regions from this process
			HaloRegion region{};
			double shift[3];
			i = k * bytesOneRegion;

			memcpy(region.rmin, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;
			memcpy(region.rmax, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;
			memcpy(region.offset, raw_neighbours.data() + i, sizeof(int) * 3);
			i += sizeof(int) * 3;
			memcpy(&region.width, raw_neighbours.data() + i, sizeof(double));
			i += sizeof(double);

			memcpy(shift, raw_neighbours.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3;

			bool enlarged[3][2] = {{false}};

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

HaloRegion NeighborAcquirer::getPotentiallyShiftedRegion(const std::array<double,3>& domainLength, const HaloRegion &region,
														 double *shiftArray, double skin) {
	for (int i = 0; i < 3; i++) {  // calculating shift
		if (region.rmin[i] >= domainLength[i] - skin) {
			shiftArray[i] = -domainLength[i];
		}
	}

	for (int i = 0; i < 3; i++) {  // calculating shift
		if (region.rmax[i] <= skin) {
			shiftArray[i] = domainLength[i];
		}
	}

	auto shiftedRegion = region;
	for (int i = 0; i < 3; i++) {  // applying shift
		shiftedRegion.rmax[i] += shiftArray[i];
		shiftedRegion.rmin[i] += shiftArray[i];
	}
	return shiftedRegion;
}

std::tuple<std::vector<HaloRegion>, std::vector<std::array<double, 3>>>
NeighborAcquirer::getAllShiftedAndNonShiftedRegionsAndShifts(HaloRegion nonShiftedRegion, HaloRegion shiftedRegion,
															 std::array<double, 3> defaultShift) {
	std::vector<HaloRegion> haloRegionVector;
	std::vector<std::array<double, 3>> shiftVector;
	std::array<std::vector<bool>, 3> doUnshiftVector;
	for (int dim = 0; dim < 3; ++dim) {
		bool isShifted = defaultShift[dim] != 0.;
		if (isShifted) {
			doUnshiftVector[dim] = {false, true};
		} else {
			doUnshiftVector[dim] = {false};
		}
	}

	for (auto doUnShiftX : doUnshiftVector[0]) {
		for (auto doUnShiftY : doUnshiftVector[1]) {
			for (auto doUnShiftZ : doUnshiftVector[2]) {
				std::array<bool, 3> doUnShift{doUnShiftX, doUnShiftY, doUnShiftZ};
				HaloRegion region{shiftedRegion};
				std::array<double, 3> shift{defaultShift};

				for (int dim = 0; dim < 3; ++dim) {
					if (doUnShift[dim]) {
						region.rmin[dim] = nonShiftedRegion.rmin[dim];
						region.rmax[dim] = nonShiftedRegion.rmax[dim];
						shift[dim] = 0;
					}
				}
				haloRegionVector.push_back(region);
				shiftVector.push_back(shift);
			}
		}
	}
	return std::make_tuple(haloRegionVector, shiftVector);
}
