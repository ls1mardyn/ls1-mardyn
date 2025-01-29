#ifndef COLLECTIVECOMMUNICATION_PERSISTENT_H_
#define COLLECTIVECOMMUNICATION_PERSISTENT_H_

#include "CollectiveCommunicationPersistent_helper.h"
#include "CollectiveCommBase.h"
#include "CollectiveCommunicationInterface.h"


//! @brief The Coll_Comm_Obj class can be used as a replacement to the default implementation of the collective communication in ls1
//! @author Mike Söhner
//! 
//! To enable this feature the ENABLE_PERSISTENT variable was to be set to true (for example via the ccmake). @par
//! The Coll_Comm_Obj class is a more optimized version of the default implementation (found in CollectiveCommunication.h) and allows for persistent collective communication.
//! This is achieved by constructing all necessary information for creating a MPI_Type and MPI_Op at compile time.
//! If all information that is required to create a persistent is given during the object construction, a persistent request will be generated.
//! It is, however, not necessary to give the constructor all information required for a persistent request. It is possible to omit the communicator,
//! the operation or the root rank. It should however be noted that the values that will be send always have to be given. It is therefore not possible to use this class
//! if the number of values communicated is decided at runtime. @par
//! In case not all information is given the usual MPI collective functions can still be used.
//! The second improvement that was implemented is a static storage of the generated MPI_Types, MPI_Ops and MPI_Requests. In the previous implementation
//! the MPI constructs were created and destroyed each time they were used, which creates overhead. By putting the values we want to communicate, the communicator
//! and the operation as template parameters we create a unique class each time we want to comunicate a new configuration. This unique class stores the data
//! necessary for the communication in static members. These members will remain, even if the current instance of the class goes out of scope. Therefore once
//! the same configuration is called again the data can be reused and there is no need to reconstruct it. @par
//! A problem occurs when it is time the destruct these MPI constructs. Static variables are destructed at the end of the program and since
//! we use the RAII programming idiom the destruction of the MPI construct would also happen then. The problem is that this happens after MPI_Finalize
//! is called and beyond that point MPI_x_free functions are no longer allowed. This is a general problem of letting modern idioms like RAII clash with older 
//! concepts such as the manual initialization and termination of constructs in C.
//! There are 2 ways to solve this problem and both try to destruct the constructs before
//! MPI_Finalize is called. The first would be to also wrap the MPI_Init and MPI_Finalize into static members of a class and make sure that class is initialized before
//! the other static members. Therefore this class would be destructed last and MPI_Finalize would be called after the other free functions. This requires an almost trivial, but
//! nontheless invasive change of the original source code. This is why the default is the second option. Nevertheless is such a wrapper provided in the 
//! CollectiveCommunicationPersistent_helper.h file.
//! It can be used in the following way:
//! @code
//!   // Replace MPI_Init(&argc, &argv); with the following code
//!   MPI_Env_Wrapper::init_environment(&argc, &argv);
//!   // Also remove calls to MPI_Finalize();
//! @endcode
//!
//! If this first option is not chosen the second option is to manually free the static members before MPI_Finalize. 
//! For this the class Coll_Comm_Deallocator in CollectiveCommunicationPersistent_helper.h
//! is provided. Coll_Comm_Obj does automatically track all create MPI constructs and registers them in the deallocator class.
//! It is used in the following way:
//! @code
//!   Coll_Comm_Deallocator::deallocate();
//!   // The above call has to happen before MPI_Finalize
//!   MPI_Finalize();
//! @endcode
//!
//! This second option is a bit more tedious and does not really fit into the modern C++ style, but is less invasive. If the first option is chosen and MPI_Init
//! and MPI_Finalize is replaced, the registration process to Coll_Comm_Deallocator can also be removed. If the first option is adopted at a larger scale, Coll_Comm_Deallocator
//! could also be removed permanently.
//!
//! Besides the points mentioned above, does Coll_Comm_Obj the same MPI functionalities as the already existing CollectiveCommunication object.
//! Currently supported commands are:
//! - reduce using add, max or min as reduce operation
//! - broadcast
//! - scan using add as operation
//!
//! Currently supported datatypes are:
//! - most basic integral and floating-point MPI datatypes (see transform functions in CollectiveCommunicationPersistent_helper.h)
//! 
//! Further mpi operations could be added very easily.
//! A typical usage of this class could look like this:
//! @code
//!   // variables to store values
//!   int i0, i1;
//!   double d;
//!   unsigned long l;
//!   // create object, select reduce operation and store values to be sent
//!   auto collComm = make_CollCommObj_AllreduceAdd(MPI_COMM_WORLD, 5, 8, 1.3, 9);
//!
//!   // execute persistent collective communication
//!   collComm.persistent();
//!
//!   // read values (IMPORTANT: same order as for storing)
//!   collComm.get(i0, i1, d, l);
//! @endcode
//! 
//! Another use case would be the following:
//! @code
//!   // variables to store values
//!   auto d0 = 1.3;
//!   auto d1 = 9.2;
//!   // create object, omit the communicator and store the values
//!   auto collComm = make_CollCommObj_AllreduceAdd(d0, d1);
//!   
//!   // execute collective communication (IMPORTANT: since we omitted the comunicator, no persistent request has been created)
//!   // we also have to give the communicator now (IMPORTANT: if we want to use the allreduce() here, we have to at least specify the operation in the constructor)
//!   collComm.allreduce(MPI_COMM_WORLD);
//!  
//!  // read values
//!  collComm.get(d0, d1);
//! @endcode
//! 
//! A final use case would be the following:
//! @code
//!   // variables to store values
//!   auto i = 2;
//!   // create object and give nothing but the values
//!   auto collComm = make_CollCommObj(i);
//!   
//!   // execute collective communication (IMPORTANT: since we omitted the comunicator, no persistent request has been created)
//!   // we also have to specify the communicator and the root rank (default is also 0) now 
//!   collComm.bcast(MPI_COMM_WORLD, 0);
//!  
//!  // read values
//!  collComm.get(i);
//! @endcode
//!
//! One last thing that has to be noted.
//! We differentiate the different classes by their signature (which includes the types of the message).
//! If we try to do nonblocking communication and want to use the same communication configuration on two 
//! different parts of the code, the same buffer will be used for both communication operations. This 
//! will most likely lead an issue with data loss. For this reason the tag argument was introduced.
//! This int argument is also part of the signture and if used will then create separate classes.
//! It can be used as follows:
//! @code
//!   // variables to store values
//!   unsigned long long ull0 = 12;
//!   short s0 = 3;
//!   unsigned long long ull1 = 46;
//!   short s1 = 2;
//!  
//!   // create objects, use tags to differentiate between the two classes
//!   auto collComm0 = make_CollCommObj_AllreduceMax<0> (MPI_COMM_WORLD, ull0, s0);
//!   auto collComm1 = make_CollCommObj_AllreduceMax<1> (MPI_COMM_WORLD, ull1, s1);
//!   
//!   // start collective communication 
//!   collComm0.persistent_start(); 
//!   collComm1.persistent_start();
//!   // wait for collective communication 
//!   collComm0.persistent_end(); 
//!   collComm1.persistent_end();
//!  
//!  // read values
//!  collComm0.get(ull0, s0);
//!  collComm1.get(ull1, s1);
//! @endcode

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
            MPI_Allreduce( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), _mpi_members->get_op(), _mpi_comm );
        else if constexpr ( Fn == static_cast<int>(MPI_CollFunctions::Bcast) )
            MPI_Bcast( _buffer.data(), 1, _mpi_members->get_type(), _root, _mpi_comm );
        else if constexpr ( Fn == static_cast<int>(MPI_CollFunctions::Scan) )
            MPI_Scan( MPI_IN_PLACE, _buffer.data(), 1, _mpi_members->get_type(), _mpi_members->get_op(), _mpi_comm );
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
     * We hide this function inside an ENABLE_PERSISTENT ifdef, because the persistent collective communication requests
     * require >= MPI 4.0. Disabling ENABLE_PERSISTENT allows lower MPI versions to compile the code.
     */
#ifdef ENABLE_PERSISTENT
    auto _init_persistent_request()
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


// Functions to explicitly create a persistent collective communication request.
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceAdd(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created. 0 means none.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, add_struct, empty_template_t, Ts...> 
                        (add_struct{}, empty_template_t{}, comm, args...);
}


// Function to explicitly create a persistent collective communication request.
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceMax(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, max_struct, empty_template_t, Ts...>
                        (max_struct{}, empty_template_t{}, comm, args...);
}

// Function to explicitly create a persistent collective communication request.
template<int tag = 0, typename... Ts>
auto makeCollCommObjAllreduceMin(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Allreduce);
    return CollCommObj<tag, fn, min_struct, empty_template_t, Ts...>
                        (min_struct{}, empty_template_t{}, comm, args...);
}


// Function to explicitly create a persistent collective communication request.
template<int tag = 0, typename... Ts>
auto makeCollCommObjBcast(int root, MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Bcast);
    return CollCommObj<tag, fn, empty_template_t, int, Ts...>
                        (empty_template_t{}, root, comm, args...);
}

// Function to explicitly create a persistent collective communication request.
template<int tag = 0, typename... Ts>
auto makeCollCommObjScanAdd(MPI_Comm comm, Ts... args) {
    // This variable determines if and which persistent request should be created.
    constexpr auto fn = static_cast<int>(MPI_CollFunctions::Scan);
    return CollCommObj<tag, fn, add_struct, empty_template_t, Ts...>
                        (add_struct{}, empty_template_t{}, comm, args...);
}

#endif
