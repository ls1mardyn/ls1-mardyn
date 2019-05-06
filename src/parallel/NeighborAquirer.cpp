/**
 * @file NeighborAquirer.cpp
 * @author seckler
 * @date 06.05.19
 */

#include "NeighborAquirer.h"
#include "HaloRegion.h"
#include "Domain.h"


/*
 * 1. Initial Exchange of all desired regions.
 * 2. Feedback from processes which own part of the region.
 *
 */
void NeighborAquirer::aquireNeighbours(Domain *domain, HaloRegion *myRegion, std::vector<HaloRegion> &desiredRegions,
									   std::vector<CommunicationPartner> &partners01,
									   std::vector<CommunicationPartner> &partners02) {
	int my_rank; // my rank
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	int num_incoming; // the number of processes in MPI_COMM_WORLD
	MPI_Comm_size(MPI_COMM_WORLD, &num_incoming);

	int num_regions = desiredRegions.size(); // the number of regions I would like to aquire from other processes


	// tell the other processes how much you are going to send
	int num_bytes_send =  sizeof(int) * 2 + (sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) * 1) * num_regions; // how many bytes am I going to send to all the other processes?
	std::vector<int> num_bytes_receive_vec(num_incoming, 0); // vector of number of bytes I am going to receive
	//MPI_Allreduce(&num_bytes_send, &num_bytes_receive, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allgather(&num_bytes_send, 1, MPI_INT, num_bytes_receive_vec.data(), 1, MPI_INT, MPI_COMM_WORLD);


	// create byte buffer
	std::vector<unsigned char> outgoing(num_bytes_send); // outgoing byte buffer
	int i = 0;
	int p = 0;

	// msg format: rank | number_of_regions | region_01 | region_02 | ...

	memcpy(outgoing.data() + i, &my_rank, sizeof(int));
	i += sizeof(int);
	memcpy(outgoing.data() + i, &num_regions, sizeof(int));
	i += sizeof(int);

	for(auto &region : desiredRegions) { // filling up the outgoing byte buffer
		memcpy(outgoing.data() + i, region.rmin, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoing.data() + i, region.rmax, sizeof(double) * 3);
		i += sizeof(double) * 3;
		memcpy(outgoing.data() + i, region.offset, sizeof(int) * 3);
		i += sizeof(int) * 3;
		memcpy(outgoing.data() + i, &region.width, sizeof(double));
		i += sizeof(double);
	}

	int num_bytes_receive = 0;
	std::vector<int> num_bytes_displacements(num_incoming, 0); // vector of number of bytes I am going to receive
	for (int j = 0; j < num_incoming; j++) {
		num_bytes_displacements[j] = num_bytes_receive;
		num_bytes_receive += num_bytes_receive_vec[j];
	}

	std::vector<unsigned char> incoming(num_bytes_receive); // the incoming byte buffer

	// send your regions
	//MPI_Allgather(&outgoing, num_bytes_send, MPI_BYTE, &incoming, num_bytes_receive, MPI_BYTE, MPI_COMM_WORLD);
	MPI_Allgatherv(outgoing.data(), num_bytes_send, MPI_BYTE, incoming.data(), num_bytes_receive_vec.data(), num_bytes_displacements.data(), MPI_BYTE, MPI_COMM_WORLD);


	std::vector<int> candidates(num_incoming, 0); // outgoing row
	std::vector<int> rec_information(num_incoming, 0); // how many bytes does each process expect?
	int bytesOneRegion = sizeof(double) * 3 + sizeof(double) * 3 + sizeof(int) * 3 + sizeof(double) + sizeof(double) * 3;
	std::vector<std::vector <unsigned char*>> sendingList (num_incoming); // the regions I own and want to send
	std::vector<CommunicationPartner> comm_partners02;

	i = 0;
	while(i != num_bytes_receive) {

		int rank;
		int regions;



		memcpy(&rank, incoming.data() + i, sizeof(int));
		i += sizeof(int); // 4
		memcpy(&regions, incoming.data() + i, sizeof(int));
		i += sizeof(int); // 4


		for(int j = 0; j < regions; j++) {
			HaloRegion region{};
			memcpy(region.rmin, incoming.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3; // 24
			memcpy(region.rmax, incoming.data() + i, sizeof(double) * 3);
			i += sizeof(double) * 3; // 24
			memcpy(region.offset, incoming.data() + i, sizeof(int) * 3);
			i += sizeof(int) * 3; // 12
			memcpy(&region.width, incoming.data() + i, sizeof(double));
			i += sizeof(double); // 4


			// msg format one region: rmin | rmax | offset | width | shift
			std::vector<double> shift(3, 0);
			double domainLength[3] = { domain->getGlobalLength(0), domain->getGlobalLength(1), domain->getGlobalLength(2) }; // better for testing
			shiftIfNeccessary(domainLength, &region, shift.data());

			if(rank != my_rank && isIncluded(myRegion, &region)) {
				candidates[rank]++; // this is a region I will send to rank

				overlap(myRegion, &region); // different shift for the overlap?

				// make a note in partners02 - don't forget to squeeze partners02
				bool enlarged[3][2] = {{ false }};
				for(int k = 0; k < 3; k++) shift[k] *= -1;

				CommunicationPartner myNewNeighbour(rank, region.rmin, region.rmax, region.rmin, region.rmax, shift.data(), region.offset, enlarged);
				comm_partners02.push_back(myNewNeighbour);

				for(int k = 0; k < 3; k++) shift[k] *= -1;

				for(int k = 0; k < 3; k++) { // shift back
					region.rmax[k] -= shift[k];
					region.rmin[k] -= shift[k];
				}

				auto* singleRegion = new unsigned char[bytesOneRegion];

				p = 0;
				memcpy(singleRegion + p, region.rmin, sizeof(double) * 3);
				p += sizeof(double) * 3;
				memcpy(singleRegion + p, region.rmax, sizeof(double) * 3);
				p += sizeof(double) * 3;
				memcpy(singleRegion + p, region.offset, sizeof(int) * 3);
				p += sizeof(int) * 3;
				memcpy(singleRegion + p, &region.width, sizeof(double));
				p += sizeof(double);
				memcpy(singleRegion + p, shift.data(), sizeof(double) * 3);
				p += sizeof(double) * 3;


				sendingList[rank].push_back(singleRegion);
			}
		}
	}


	// squeeze here
	if(not comm_partners02.empty()) {
		std::vector<CommunicationPartner> squeezed = squeezePartners(comm_partners02);
		partners02.insert(partners02.end(), squeezed.begin(), squeezed.end());
	}

	std::vector<unsigned char *> merged (num_incoming); // Merge each list of char arrays into one char array
	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j]  > 0) {
			auto* mergedRegions = new unsigned char[candidates[j] * bytesOneRegion];

			for(int k = 0; k < candidates[j]; k++) {
				memcpy(mergedRegions + k * bytesOneRegion, sendingList[j][k], bytesOneRegion);
			}

			merged[j] = mergedRegions;
		}
	}

	// delete sendingList

	for(const auto& one : sendingList){
		for (auto two : one){
			delete[] two;
		}
	}


	// The problem is, that I do not know how many regions I am going to receive from each process
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
	 *
	 * Each process has a horizontal vector, where it marks how many regions it is going to send to another process.
	 * As allToAll is too slow or apparently too unpredictable regarding runtime,
	 * AllReduce is used in combination with probe.
	 *
	 */

	MPI_Allreduce(candidates.data(), rec_information.data(), num_incoming,MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	// all the information for the final information exchange has been collected -> final exchange

	std::vector<MPI_Request> requests(num_incoming, MPI_REQUEST_NULL);
	MPI_Status probe_status;
	MPI_Status rec_status;

	// sending (non blocking)
	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j] > 0) {
			MPI_Isend(merged[j], candidates[j] * bytesOneRegion, MPI_BYTE, j, 1, MPI_COMM_WORLD, &requests[j]); // tag is one
		}
	}

	std::vector<CommunicationPartner> comm_partners01; // the communication partners

	// receive data (blocking)
	int byte_counter = 0;

	while(byte_counter < rec_information[my_rank] * bytesOneRegion) { // TODO: I do not know yet, how many regions I received

		// MPI_PROBE
		MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &probe_status);
		// interpret probe
		int source = probe_status.MPI_SOURCE;
		int bytes;
		MPI_Get_count(&probe_status, MPI_BYTE, &bytes);


		byte_counter += bytes;
		// create buffer
		std::vector<unsigned char> raw_neighbours(bytes);
		MPI_Recv(raw_neighbours.data(), bytes, MPI_BYTE, source, 1, MPI_COMM_WORLD, &rec_status);
		// Interpret Buffer and add neighbours
		for(int k = 0; k < (bytes / bytesOneRegion); k++) { // number of regions from this process
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

			bool enlarged[3][2] = {{ false }};

			// CommunicationPartner(const int r, const double hLo[3], const double hHi[3], const double bLo[3], const double bHi[3], const double sh[3], const int offset[3], const bool enlarged[3][2]) {
			CommunicationPartner myNewNeighbour(source, region.rmin, region.rmax, region.rmin, region.rmax, shift, region.offset, enlarged); // DO NOT KNOW ABOUT THE 0s
			comm_partners01.push_back(myNewNeighbour);

		}
	}


	for(int j = 0; j < num_incoming; j++) {
		if(candidates[j] > 0)
			MPI_Wait(&requests[j], MPI_STATUS_IGNORE);
	}

	for (auto two : merged) {
		delete[] two;
	}

	if(not comm_partners01.empty()) {
		std::vector<CommunicationPartner> squeezed = squeezePartners(comm_partners01);
		partners01.insert(partners01.end(), squeezed.begin(), squeezed.end());
	}


	MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<CommunicationPartner> NeighborAquirer::squeezePartners(const std::vector<CommunicationPartner>& partners) {
	std::vector<CommunicationPartner> squeezedPartners;
	std::vector<bool> used(partners.size(), false); // flag table, that describes, whether a certain comm-partner has already been added
	for (unsigned int i = 0; i < partners.size(); i++) {
		if (used[i])
			continue;  // if we already added the neighbour, don't add it again!
		int rank = partners[i].getRank();
		CommunicationPartner tmp = partners[i];
		for (unsigned int j = i + 1; j < partners.size(); j++) {
			if (partners[j].getRank() != rank)
				continue;  // only add those with same rank
			tmp.add(partners[j]);
			used[j] = true;
		}
		squeezedPartners.push_back(tmp);
	}
	return squeezedPartners;
}

bool NeighborAquirer::isIncluded(HaloRegion *myRegion, HaloRegion *inQuestion) {
	return myRegion->rmax[0] > inQuestion->rmin[0] && myRegion->rmin[0] < inQuestion->rmax[0]
	       && myRegion->rmax[1] > inQuestion->rmin[1] && myRegion->rmin[1] < inQuestion->rmax[1]
	       && myRegion->rmax[2] > inQuestion->rmin[2] && myRegion->rmin[2] < inQuestion->rmax[2];
	// myRegion->rmax > inQuestion->rmin
	// && myRegion->rmin < inQuestion->rmax
}

void NeighborAquirer::overlap(HaloRegion *myRegion, HaloRegion *inQuestion) {
	/*
	 * m = myRegion, q = inQuestion, o = overlap
	 * i)  m.max < q.max ?
	 * ii) m.min < q.min ?
	 *
	 * i) | ii) | Operation
	 * -------------------------------------------
	 *  0 |  0  | o.max = q.max and o.min = m.min
	 *  0 |  1  | o.max = q.max and o.min = q.min
	 *  1 |  0  | o.max = m.max and o.min = m.min
	 *  1 |  1  | o.max = m.max and o.min = q.min
	 *
	 */
	HaloRegion overlap{};


	for(int i = 0; i < 3; i++) {
		if(myRegion->rmax[i] < inQuestion->rmax[i]) { // 1
			if(myRegion->rmin[i] < inQuestion->rmin[i]) { // 1 1
				overlap.rmax[i] = myRegion->rmax[i];
				overlap.rmin[i] = inQuestion->rmin[i];
			} else { // 1 0
				overlap.rmax[i] = myRegion->rmax[i];
				overlap.rmin[i] = myRegion->rmin[i];
			}
		} else { // 0
			if(myRegion->rmin[i] < inQuestion->rmin[i]) { // 0 1
				overlap.rmax[i] = inQuestion->rmax[i];
				overlap.rmin[i] = inQuestion->rmin[i];
			} else { // 0 0
				overlap.rmax[i] = inQuestion->rmax[i];
				overlap.rmin[i] = myRegion->rmin[i];
			}
		}
	}

	// adjust width and offset?
	memcpy(inQuestion->rmax, overlap.rmax, sizeof(double) * 3);
	memcpy(inQuestion->rmin, overlap.rmin, sizeof(double) * 3);
}

void NeighborAquirer::shiftIfNeccessary(const double *domainLength, HaloRegion *region, double *shiftArray) {
	for(int i = 0; i < 3; i++) // calculating shift
		if(region->rmin[i] >= domainLength[i])
			shiftArray[i] = -domainLength[i];

	for(int i = 0; i < 3; i++) // calculating shift
		if(region->rmax[i] <= 0)
			shiftArray[i] = domainLength[i];

	for(int i = 0; i < 3; i++) { // applying shift
		region->rmax[i] += shiftArray[i];
		region->rmin[i] += shiftArray[i];
	}
}