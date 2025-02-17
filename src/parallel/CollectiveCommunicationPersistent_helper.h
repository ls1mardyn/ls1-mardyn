#ifndef COLLECTIVECOMMUNICATION_PERSISTENT_HELPER_H_
#define COLLECTIVECOMMUNICATION_PERSISTENT_HELPER_H_

#include <iostream>
#include <vector>
#include <array>
#include <cstddef>
#include <array>
#include <cstring>
#include <bit>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <mpi.h>

/** @brief This file contains auxiliary classes and functions for the CollCommObj in CollectiveCommunicationPersistent.h
*   @author Mike Söhner
*/


/** This is a wrapper class that wraps MPI_Init and MPI_finalize statically.
* 
* This class wraps the MPI_Init and MPI_Finalize functions into static members of a class. Since we call this wrapper
* before all the other MPI functions we make sure that the static members of this class are constructed before the
* other static class members and destructed after the other static class members. Therefore
* this class would be destructed last and MPI_Finalize would be called after the other free functions.
* It can be used in the following way:
* @code
*   // Replace MPI_Init(&argc, &argv); with the following code
*   MPI_Env_Wrapper::init_environment(&argc, &argv);
*   // Also remove calls to MPI_Finalize();
* @endcode
 */
class MPI_Env_Wrapper {
public:
    static auto init_environment(int* argc, char*** argv) {
        static MPI_Environment mpi_env(argc, argv);
    }
    static auto init_thread_environment(int* argc, char*** argv, int required) {
        static MPI_Environment mpi_env(argc, argv, required);
    }
private:
    // class to manage MPI_Init and MPI_Finalize
    class MPI_Environment {
    public:
        MPI_Environment(int* argc, char*** argv) {
            MPI_Init(argc, argv);
        }
        MPI_Environment(int* argc, char*** argv, int required) {
            int provided;
            MPI_Init_thread(argc, argv, required, &provided);
            if (provided < required) {
                std::cerr << "Cannot provide requested level of MPI thread support." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        ~MPI_Environment() {
            MPI_Finalize();
        }
    };
};


// enums to determine which MPI function should be used
enum class MPI_CollFunctions {
    None, Allgather, Allreduce, Alltoall, Bcast, Gather, Reduce, Reduce_scatter, Scatter, Scan
};


// use this function to get the combined size of all the types inside a variadic template
// overload
template<typename U>
constexpr size_t get_size(U head) {
    return sizeof(head);
}
// base case: used when pack is non-empty
template<typename U, typename... Us>
constexpr size_t get_size(U head, Us... tail) {
    return sizeof(head) + get_size(tail...);
}

// use this function to get the maximum size of all the types inside a variadic template
// overload
template<typename U>
constexpr size_t get_max_size(U head) {
    return sizeof(head);
}
// base case: used when pack is non-empty
template<typename U, typename... Us>
constexpr size_t get_max_size(U head, Us... tail) {
    return std::max(sizeof(head), get_max_size(tail...));
}

class Type2MPI {
public:
    static MPI_Datatype transform(char) { return MPI_CHAR; }
    static MPI_Datatype transform(short) { return MPI_SHORT; }
    static MPI_Datatype transform(int) { return MPI_INT; }
    static MPI_Datatype transform(long) { return MPI_LONG; }
    static MPI_Datatype transform(long long) { return MPI_LONG_LONG; }
    static MPI_Datatype transform(unsigned char) { return MPI_UNSIGNED_CHAR; }
    static MPI_Datatype transform(unsigned short) { return MPI_UNSIGNED_SHORT; }
    static MPI_Datatype transform(unsigned int) { return MPI_UNSIGNED; }
    static MPI_Datatype transform(unsigned long) { return MPI_UNSIGNED_LONG; }
    static MPI_Datatype transform(unsigned long long) { return MPI_UNSIGNED_LONG_LONG; }
    static MPI_Datatype transform(float) { return MPI_FLOAT; }
    static MPI_Datatype transform(double) { return MPI_DOUBLE; }
    static MPI_Datatype transform(long double) { return MPI_LONG_DOUBLE; }
};


// fill arrays to create a MPI_Datatype with them
// overload
template<typename T>
constexpr void free_fill_type_array(MPI_Datatype* t, T head) {
    *t = Type2MPI::transform(head);
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
constexpr void free_fill_type_array(MPI_Datatype* t, T head, Ts... tail) {
    *t = Type2MPI::transform(head);
    free_fill_type_array(t + 1, tail...);
}
// helper function called by the constructor
template<typename... Ts>
constexpr void fill_type_array(MPI_Datatype* t, Ts... args) {
    free_fill_type_array( t, args... );
}

// overload
template<typename T>
constexpr void free_fill_displs_array(MPI_Aint* d, size_t size, size_t offset, T head) {
    *d = offset;
}
// base case
template<typename T, typename... Ts>
constexpr void free_fill_displs_array(MPI_Aint* d, size_t size, size_t offset, T head, Ts... tail) {
    *d = offset;
    offset += size;
    free_fill_displs_array( d + 1, size, offset, tail... );
}
// helper function used to fill displacement array required to create custom MPI type
template<typename... Ts>
constexpr void fill_displs_array(MPI_Aint* d, Ts... args) {
    size_t size = get_max_size(args...);
    size_t offset = 0;
    free_fill_displs_array( d, size, offset, args... );
}

// put values of variadic template into array of bytes
// overload
template<typename T>
constexpr void fill_buffer(std::byte* buffer, size_t size, size_t offset, T head) {
    std::memcpy(buffer + offset, &head, sizeof(head));
}
// base case
template<typename T, typename... Ts>
constexpr void fill_buffer(std::byte* buffer, size_t size, size_t offset, T head, Ts... tail) {
    std::memcpy(buffer + offset, &head, sizeof(head));
    offset += size;
    fill_buffer(buffer, size, offset, tail...);
}
// helper function
template<typename... Ts>
constexpr void free_fill_buffer(std::byte* buffer, Ts... args) {
    size_t size = get_max_size(args...);
    size_t offset = 0;
    fill_buffer(buffer, size, offset, args...);
}

// get values out of byte vector into variadic template values
// overload
template<int id, int size, typename T, typename... TTs>
constexpr void free_get(std::byte* buffer, std::tuple<TTs...>& tuple, size_t offset, T head) {
    std::memcpy(&head, buffer + offset, sizeof(head));
    std::get<id>(tuple) = head;
}
// base case: used when pack is non-empty
template<int id, int size, typename T, typename... Ts, typename... TTs>
constexpr void free_get(std::byte* buffer, std::tuple<TTs...>& tuple, size_t offset, T head, Ts... tail) {
    std::memcpy(&head, buffer + offset, sizeof(head));
    std::get<id>(tuple) = head;

    offset += size;

    free_get<id+1, size>(buffer, tuple, offset, tail...);
}
// helper function
template<typename... Ts>
constexpr std::tuple<Ts...> helper_get(std::byte* buffer) {
    std::tuple<Ts...> tuple;

    constexpr size_t size = get_max_size(Ts{}...);
    size_t offset = 0;
    free_get<0, size>(buffer, tuple, offset, Ts{}...);

    return tuple;
}

// print variadic template values
// overload
template<typename T>
void print(std::vector<std::byte>& buffer, size_t offset, T head) {
    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    std::cout << value << std::endl;
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
void print(std::vector<std::byte>& buffer, size_t& offset, T head, Ts... tail) {
    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    offset += sizeof(T);
    std::cout << value << std::endl;
    print(buffer, offset, tail...);
}

// struct for adding
// overload
template<typename T>
void free_add( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec += *pinvec;
}
// base case
template<typename T, typename... Ts>
void free_add( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec += *pinvec;
    
    offset += size;
    free_add(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this adds the values of different types together
struct add_struct {
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype ) {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_add( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// struct for maxing
// overload
template<typename T>
void free_max( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec = std::max(*pinoutvec, *pinvec);
    
    offset += size;
}
// base case
template<typename T, typename... Ts>
void free_max( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec = std::max(*pinoutvec, *pinvec);
    
    offset += size;
    free_max(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this maxs the values of different types together
struct max_struct {
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype ) {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_max( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// struct for mining
// overload
template<typename T>
void free_min( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec = std::min(*pinoutvec, *pinvec);
    
    offset += size;
}
// base case
template<typename T, typename... Ts>
void free_min( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail ) {
    T *pinvec = reinterpret_cast<T*> (invec + offset);
    T *pinoutvec = reinterpret_cast<T*> (inoutvec + offset);

    *pinoutvec = std::min(*pinoutvec, *pinvec);
    
    offset += size;
    free_min(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this mins the values of different types together
struct min_struct {
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype ) {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_min( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// Pre-C++20 constraint to formulate the requirement that a template arg is an operation type
template < typename T >
struct is_op
: public std::false_type {};

template <>
struct is_op < add_struct >
: public std::true_type {};

template <>
struct is_op < max_struct >
: public std::true_type {};

template <>
struct is_op < min_struct >
: public std::true_type {};

template< typename T >
constexpr bool is_op_v = is_op<T>::value;


// Pre-C++20 constraint to formulate the requirement that a template arg is a communicator
template < typename T >
struct is_comm
: public std::false_type {};

template <>
struct is_comm < MPI_Comm >
: public std::true_type {};

template< typename T >
constexpr bool is_comm_v = is_comm<T>::value;

// Pre-C++20 constraint to formulate the requirement that a function is an allreduce
template < int I >
struct is_allreduce
: public std::false_type {};

template <>
struct is_allreduce < static_cast<int>(MPI_CollFunctions::Allreduce) >
: public std::true_type {};

template< int T >
constexpr bool is_allreduce_v = is_allreduce<T>::value;

// Pre-C++20 constraint to formulate the requirement that a function is a bcast
template < int I >
struct is_bcast
: public std::false_type {};

template <>
struct is_bcast < static_cast<int>(MPI_CollFunctions::Bcast) >
: public std::true_type {};

template< int T >
constexpr bool is_bcast_v = is_bcast<T>::value;

// Pre-C++20 constraint to formulate the requirement that a function is a bcast
template < int I >
struct is_scan
: public std::false_type {};

template <>
struct is_scan < static_cast<int>(MPI_CollFunctions::Scan) >
: public std::true_type {};

template< int T >
constexpr bool is_scan_v = is_scan<T>::value;


#endif
