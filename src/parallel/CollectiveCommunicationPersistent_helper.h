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
#include <mpi.h>

//! @brief This file contains auxiliary classes and functions for the Coll_Comm_Obj in mpi_perscomm_obj.h
//! @author Mike Söhner

// class to manage MPI_Init and MPI_Finalize
class MPI_Environment
{
public:
    MPI_Environment(int* argc, char*** argv)
    {
        MPI_Init(argc, argv);
    }
    MPI_Environment(int* argc, char*** argv, int required)
    {
        int provided;
        MPI_Init_thread(argc, argv, required, &provided);
        if (provided < required)
        {
            std::cerr << "Cannot provide requested level of MPI thread support." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    ~MPI_Environment()
    {
        MPI_Finalize();
    }
};


// enums to determine which MPI function should be used
enum class MPI_CollFunctions
{
    None, Allgather, Allreduce, Alltoall, Bcast, Gather, Reduce, Reduce_scatter, Scatter, Scan
};


// use this function to get the combined size of all the types inside a variadic template
// overload
template<typename U>
constexpr size_t get_size(U head)
{
    return sizeof(head);
}
// base case: used when pack is non-empty
template<typename U, typename... Us>
constexpr size_t get_size(U head, Us... tail)
{
    return sizeof(head) + get_size(tail...);
}

// use this function to get the maximum size of all the types inside a variadic template
// overload
template<typename U>
constexpr size_t get_max_size(U head)
{
    return sizeof(head);
}
// base case: used when pack is non-empty
template<typename U, typename... Us>
constexpr size_t get_max_size(U head, Us... tail)
{
    return std::max(sizeof(head), get_max_size(tail...));
}

class Type2MPI
{
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
constexpr void free_fill_type_array(MPI_Datatype* t, T head)
{
    *t = Type2MPI::transform(head);
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
constexpr void free_fill_type_array(MPI_Datatype* t, T head, Ts... tail)
{
    *t = Type2MPI::transform(head);
    free_fill_type_array(t + 1, tail...);
}
// helper function called by the constructor
template<typename... Ts>
constexpr void fill_type_array(MPI_Datatype* t, Ts... args)
{
    free_fill_type_array( t, args... );
}

// overload
template<typename T>
constexpr void free_fill_displs_array(MPI_Aint* d, size_t size, size_t offset, T head)
{
    *d = offset;
}
// base case
template<typename T, typename... Ts>
constexpr void free_fill_displs_array(MPI_Aint* d, size_t size, size_t offset, T head, Ts... tail)
{
    *d = offset;
    offset += size;
    free_fill_displs_array( d + 1, size, offset, tail... );
}
// helper function used to fill displacement array required to create custom MPI type
template<typename... Ts>
constexpr void fill_displs_array(MPI_Aint* d, Ts... args)
{
    size_t size = get_max_size(args...);
    size_t offset = 0;
    free_fill_displs_array( d, size, offset, args... );
}

// put values of variadic template into array of bytes
// overload
template<typename T>
constexpr void fill_buffer(std::byte* buffer, size_t size, size_t offset, T head)
{
    std::memcpy(buffer + offset, &head, sizeof(head));
}
// base case
template<typename T, typename... Ts>
constexpr void fill_buffer(std::byte* buffer, size_t size, size_t offset, T head, Ts... tail)
{
    std::memcpy(buffer + offset, &head, sizeof(head));
    offset += size;
    fill_buffer(buffer, size, offset, tail...);
}
// helper function
template<typename... Ts>
constexpr void free_fill_buffer(std::byte* buffer, Ts... args)
{
    size_t size = get_max_size(args...);
    size_t offset = 0;
    fill_buffer(buffer, size, offset, args...);
}

// get values out of byte vector into variadic template values
// overload
template<typename T>
void free_get(std::byte* buffer, size_t size, size_t offset, T& head)
{
    std::memcpy(&head, buffer + offset, sizeof(head));
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
void free_get(std::byte* buffer, size_t size, size_t offset, T& head, Ts&... tail)
{
    std::memcpy(&head, buffer + offset, sizeof(head));
    offset += size;
    free_get(buffer, size, offset, tail...);
}
// helper function
template<typename... Ts>
void helper_get(std::byte* buffer, Ts&... args)
{
    size_t size = get_max_size(args...);
    size_t offset = 0;
    free_get(buffer, size, offset, args...);
}

// print variadic template values
// overload
template<typename T>
void print(std::vector<std::byte>& buffer, size_t offset, T head)
{
    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    std::cout << value << std::endl;
}
// base case: used when pack is non-empty
template<typename T, typename... Ts>
void print(std::vector<std::byte>& buffer, size_t& offset, T head, Ts... tail)
{
    T value{};
    std::memcpy(&value, buffer.data() + offset, sizeof(T));
    offset += sizeof(T);
    std::cout << value << std::endl;
    print(buffer, offset, tail...);
}

// struct for adding
// overload
template<typename T>
void free_add( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec += *pinvec;
}
// base case
template<typename T, typename... Ts>
void free_add( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec += *pinvec;
    
    offset += size;
    free_add(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this adds the values of different types together
struct add_struct
{
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype )
    {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_add( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// struct for maxing
// overload
template<typename T>
void free_max( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec = std::max(*pinoutvec, *pinvec);
    
    offset += size;
}
// base case
template<typename T, typename... Ts>
void free_max( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec = std::max(*pinoutvec, *pinvec);
    
    offset += size;
    free_max(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this maxs the values of different types together
struct max_struct
{
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype )
    {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_max( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// struct for mining
// overload
template<typename T>
void free_min( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec = std::min(*pinoutvec, *pinvec);
    
    offset += size;
}
// base case
template<typename T, typename... Ts>
void free_min( std::byte *invec, std::byte *inoutvec, size_t size, size_t offset, T head, Ts... tail )
{
    T *pinvec = (T *) (invec + offset);
    T *pinoutvec = (T *) (inoutvec + offset);

    *pinoutvec = std::min(*pinoutvec, *pinvec);
    
    offset += size;
    free_min(invec, inoutvec, size, offset, tail...);
}
// class that contains function for MPI_Op
// this mins the values of different types together
struct min_struct
{
    template<typename... Ts>
    static void func( void *invec, void *inoutvec, int *len, MPI_Datatype *datatype )
    {
        size_t size = get_max_size(Ts{}...);
        size_t offset = 0;
        free_min( (std::byte *) invec, (std::byte *) inoutvec, size, offset, Ts{}... );
    }
};


// class to handle deallocation of MPI allocations
class Coll_Comm_Deallocator
{
public:
    static auto emplace_back(MPI_Request* request)
    {
        _allocated_requests.emplace_back(request);
    }

    static auto emplace_back(MPI_Op* op)
    {
        _allocated_ops.emplace_back(op);
    }

    static auto emplace_back(MPI_Datatype* type)
    {
        _allocated_types.emplace_back(type);
    }

    static void deallocate()
    {
        for (auto ptr : _allocated_requests)
            MPI_Request_free(ptr);
        for (auto ptr : _allocated_ops)
            MPI_Op_free(ptr);
        for (auto ptr : _allocated_types)
            MPI_Type_free(ptr);
    }

private:
    static inline std::vector<MPI_Request*> _allocated_requests;
    static inline std::vector<MPI_Op*> _allocated_ops;
    static inline std::vector<MPI_Datatype*> _allocated_types;
};

#endif