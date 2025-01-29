#ifndef COLLECTIVECOMMUNICATION_PERSISTENT_H_
#define COLLECTIVECOMMUNICATION_PERSISTENT_H_

#include "CollectiveCommunicationPersistent_helper.h"
#include "CollectiveCommBase.h"
#include "CollectiveCommunicationInterface.h"

/** @brief This class can be used as a replacement to the default implementation of the collective communication.
* @author Mike Söhner
* 
* The CollCommObj class is a more optimized version of the default implementation (found in CollectiveCommunication.h) 
* and allows for persistent collective communication (for this feature enable ENABLE_PERSISTENT).
* If all information that is required to create a persistent is given during the object construction, a persistent
* request will be generated. This is achieved by requiring all necessary information for creating a MPI_Type and MPI_Op
* at compile time. If ENABLE_PERSISTENT is disabled the implementation falls back to the default MPI collectives 
* (MPI_Allreduce, MPI_Bcast, MPI_Scan), but this will still be more optimized than the default implementation. 
* This improvement is due to a static storage of the generated MPI_Type, MPI_Op and MPI_Request which allows for a
* reuse of those constructs. Given the same configuration must be communicated again. The default implementation 
* created and deleted the MPI constructs every time a collective communication had to be performed. This implementation
* eliminates that overhead. 
*
* By putting the values we want to communicate and the operation as template parameters we
* create a unique class each time we want to comunicate a new configuration. This unique class stores the data
* necessary for the communication in static members. These members will remain, even if the current instance of the
* class goes out of scope. Therefore once the same configuration is called again the data can be reused and there is
* no need to reconstruct it.
*
* A problem occurs when it is time the destruct these MPI constructs. Static variables are destructed at the end of
* the program and since we use the RAII programming idiom, the destruction of the MPI construct would also happen at
* that time. The problem is that this happens after MPI_Finalize is called and beyond that point MPI_x_free functions
* are no longer allowed. Hence we also wrap the MPI_Init and MPI_Finalize functions into a static object.
* This object is created before the MPI constructs and therefore will be destructed after the other static objects, 
* restoring the correct order creating and freeing MPI constructs. The wrapper can be found in 
* CollectiveCommunicationPersistent_helper.h.
*
* Besides the points mentioned above, does CollCommObj the same MPI functionalities as the already existing 
* CollectiveCommunication object.
* Currently supported commands are:
* - reduce using add, max or min as reduce operation
* - broadcast
* - scan using add as operation
*
* Currently supported datatypes are:
* - most basic integral and floating-point MPI datatypes 
*           (see transform functions in CollectiveCommunicationPersistent_helper.h)
* 
* Further mpi operations could be added very easily.
* A typical usage of this class could look like this:
* @code
*   // create object, select reduce operation and store values to be sent
*   auto collComm = makeCollCommObjAllreduceAdd(MPI_COMM_WORLD, 5, 8, 1.3, 9l);
*
*   // execute collective communication
*   collComm.communicate();
*
*   // create variables and get values (IMPORTANT: same order as for storing)
*   auto [i0, i1, d, l] = collComm.get(); // or use std::tie() if the variables were already declared.
* @endcode
*
* One last thing that has to be noted.
* We differentiate the different classes by their template parameters (which includes the types of the message).
* If we try to do communication and want to use the same communication configuration on two communication object
* and let them run simultaneously, the same buffer will be used for both communication operations. 
* This will most likely lead to an issue with data loss. For this reason the tag argument was introduced.
* This int argument is also part of the signture and if used will then create separate classes.
* (This is only necessary in the ENABLE_PERSISTENT=On case, as otherwise nonstatic buffers will be used.)
* It can be used as follows:
* @code
*   // create objects, use tags to differentiate between the two classes
*   auto collComm0 = makeCollCommObjAllreduceMax<0> (MPI_COMM_WORLD, 12ul, 3);
*   auto collComm1 = makeCollCommObjAllreduceMax<1> (MPI_COMM_WORLD, 42ul, 2);
* @endcode
*
*
* @tparam tag Used to change the template parameters, so a different class will be generated.
* @tparam Fn Determines which collective communication function is used.
* @tparam Op Determines the data processing of the collective communication (eg. a sum for allreduce).
* @tparam Root Determines the root rank.
* @tparam Ts... A list of all elements that will take part in the collective communication.
*
* @param comm MPI communicator used for communication.
*/

template<int tag, int Fn, typename Op, typename Root, typename... Ts>
class CollCommObj {
public:

    /**
     * Construcor of the class
     *
     * @param op Contains the operation that is performed during the collective communication (eg. sum for allreduce).
     * @param root The root rank of the collective communication.
     * @param comm The MPI communicator used for the collective communication.
     * @param args... The elements that will take part in the collective communication.
     */

    CollCommObj(Op op, Root root, MPI_Comm comm, Ts... args) {
        // fill the buffer, this always has to hapen
        free_fill_buffer(_buffer.data(), args...);
        // check if this combination of values has already been constructed (this is for MaMico compatibility)
        if ( _mpi_members == nullptr ) {
            // set communicator
            _mpi_comm = comm;
            // create the static members
            static MPI_Members members;
            _mpi_members = &members;
            // create variables necessary for MPI Datatype
            std::array<int, sizeof...(Ts)> array_of_blocklengths;
            array_of_blocklengths.fill(1);
            std::array<MPI_Aint, sizeof...(Ts)> array_of_displacements;
            std::array<MPI_Datatype, sizeof...(Ts)> array_of_types;
            // fill arrays necessary for custom mpi type        
            fill_displs_array(array_of_displacements.data(), args...);
            fill_type_array(array_of_types.data(), args...);
            // create custom mpi type
            MPI_Type_create_struct(sizeof...(Ts), array_of_blocklengths.data(), array_of_displacements.data(),
                                    array_of_types.data(), &_mpi_members->get_type());
            MPI_Type_commit(&_mpi_members->get_type());
            // check if we were given a custom mpi operation
            if constexpr( is_op_v<Op> ) {
                // Create MPI Operation
                MPI_Op_create(Op::template func<Ts...>, 1, &_mpi_members->get_op());
            }
            // check if we were given a root rank
            if constexpr( std::is_same_v<Root, int> )
                _root = root;
#ifdef ENABLE_PERSISTENT
            // check if we should create a persistent communication request
            _init_persistent_request();
#endif
	    }
    }

    /**
     * Starts the communication and it is only available if enough information has been given to the constructor.
     */

    auto communicate()
        // Checks if all necessary data for either an allreduce, bcast or scan has been given in the constructor.
        -> std::enable_if_t<
            (is_allreduce_v<Fn> && is_op_v<Op>) ||
            (is_bcast_v<Fn> && std::is_same_v<Root, int>) ||
            (is_scan_v<Fn> && is_op_v<Op>)
        > {
#ifdef ENABLE_PERSISTENT
        MPI_Start( &_mpi_members->get_request() );
        MPI_Wait( &_mpi_members->get_request(), MPI_STATUS_IGNORE );
#else
        if constexpr ( Fn == static_cast<int>(MPI_CollFunctions::Allreduce) )
            MPI_Allreduce( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), 
                            _mpi_members->get_op(), _mpi_comm );

        else if constexpr ( Fn == static_cast<int>(MPI_CollFunctions::Bcast) )
            MPI_Bcast( _buffer.data(), 1, _mpi_members->get_type(), _root, _mpi_comm );
        
        else if constexpr ( Fn == static_cast<int>(MPI_CollFunctions::Scan) )
            MPI_Scan( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), 
                            _mpi_members->get_op(), _mpi_comm );
#endif
    }

    /**
     * Get the communicated values out of the buffer.
     */
    std::tuple<Ts...> get() {
        return helper_get<Ts...> (_buffer.data());
    }


private:
    /**
     * Creates an persistent collective communication request, if the correct parameters were given.
     * 
     * We hide this function inside an ENABLE_PERSISTENT ifdef, because the persistent collective communication
     * requests require >= MPI 4.0. Disabling ENABLE_PERSISTENT allows lower MPI versions to compile the code.
     */
#ifdef ENABLE_PERSISTENT
    auto _init_persistent_request()
        // Checks if all necessary data for either an allreduce, bcast or scan has been given in the constructor.
        -> std::enable_if_t<
            (is_allreduce_v<Fn> && is_op_v<Op>) ||
            (is_bcast_v<Fn> && std::is_same_v<Root, int>) ||
            (is_scan_v<Fn> && is_op_v<Op>)
        > {
        if constexpr (Fn == static_cast<int>(MPI_CollFunctions::Allreduce)) {
            MPI_Allreduce_init( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), _mpi_members->get_op(),
                            _mpi_comm, MPI_INFO_NULL, &_mpi_members->get_request() );
        }
        else if constexpr (Fn == static_cast<int>(MPI_CollFunctions::Bcast)) {
            MPI_Bcast_init( _buffer.data(), 1, _mpi_members->get_type(), _root, 
                            _mpi_comm, MPI_INFO_NULL, &_mpi_members->get_request() );
        }
        else if constexpr (Fn == static_cast<int>(MPI_CollFunctions::Scan)) {
            MPI_Scan_init( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), _mpi_members->get_op(), 
                            _mpi_comm, MPI_INFO_NULL, &_mpi_members->get_request() );
        }
    }
#endif

    // class to hide MPI types as private members
    class MPI_Members {
    public:
        MPI_Members() = default;

        MPI_Datatype& get_type() { return _type; }
        MPI_Op& get_op() { return _op; }
        MPI_Request& get_request() { return _request; }

        ~MPI_Members() {
            if (_type != MPI_DATATYPE_NULL)
                MPI_Type_free(&_type);
            if (_op != MPI_OP_NULL)
                MPI_Op_free(&_op);
            if (_request != MPI_REQUEST_NULL)
                MPI_Request_free(&_request);
        }

    private:
        MPI_Datatype _type = MPI_DATATYPE_NULL;
        MPI_Op _op = MPI_OP_NULL;
        MPI_Request _request = MPI_REQUEST_NULL;
    };


    // private member variables

    // it is more portable to make the communcation buffer non-static, 
    // but then persistent collective would not be possible in the current implementation
#ifdef ENABLE_PERSISTENT
    alignas(16) static inline std::array<std::byte, sizeof...(Ts) * get_max_size(Ts{}...)> _buffer{};
#else
    alignas(16) std::array<std::byte, sizeof...(Ts) * get_max_size(Ts{}...)> _buffer{};
#endif
    static inline MPI_Members* _mpi_members {};
    static inline MPI_Comm _mpi_comm {};
    static inline int _root {};
};

// this struct is used to leave some template arguments "empty"
// this will trigger the SFINAE to reject inappropriate functions
struct empty_template_t{};


/**
 * Function to explicitly create a CollCommObj that does an allreduce with the sum operation.
 *
 * @tparam tag Used to change the template parameters, so a different class will be generated.
 * 
 * @param comm The MPI communicator used for the collective communication.
 * @param args... The elements that will take part in the collective communication.
 */
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceAdd(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created. 0 means none.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, add_struct, empty_template_t, Ts...> 
                        (add_struct{}, empty_template_t{}, comm, args...);
}

/**
 * Function to explicitly create a CollCommObj that does an allreduce with the max operation.
 *
 * @tparam tag Used to change the template parameters, so a different class will be generated.
 * 
 * @param comm The MPI communicator used for the collective communication.
 * @param args... The elements that will take part in the collective communication.
 */
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceMax(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, max_struct, empty_template_t, Ts...>
                        (max_struct{}, empty_template_t{}, comm, args...);
}

/**
 * Function to explicitly create a CollCommObj that does an allreduce with the min operation.
 *
 * @tparam tag Used to change the template parameters, so a different class will be generated.
 * 
 * @param comm The MPI communicator used for the collective communication.
 * @param args... The elements that will take part in the collective communication.
 */
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceMin(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, min_struct, empty_template_t, Ts...>
                        (min_struct{}, empty_template_t{}, comm, args...);
}

/**
 * Function to explicitly create a CollCommObj that does a bcast.
 *
 * @tparam tag Used to change the template parameters, so a different class will be generated.
 * 
 * @param root The root rank of the MPI communication.
 * @param comm The MPI communicator used for the collective communication.
 * @param args... The elements that will take part in the collective communication.
 */
template<int tag = 0, typename... Ts>
auto makeCollCommObjBcast(int root, MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Bcast);
    return CollCommObj<tag, fn, empty_template_t, int, Ts...>
                        (empty_template_t{}, root, comm, args...);
}

/**
 * Function to explicitly create a CollCommObj that does a scan with the sum operation.
 *
 * @tparam tag Used to change the template parameters, so a different class will be generated.
 * 
 * @param comm The MPI communicator used for the collective communication.
 * @param args... The elements that will take part in the collective communication.
 */
template<int tag = 0, typename... Ts>
auto makeCollCommObjScanAdd(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Scan);
    return CollCommObj<tag, fn, add_struct, empty_template_t, Ts...>
                        (add_struct{}, empty_template_t{}, comm, args...);
}

#endif
