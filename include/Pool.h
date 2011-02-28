#include <stdlib.h>
#include <stddef.h>
 
/** 
 *
 * @class   Pool
 * @author Buote Xu (Buote.Xu@stud.uni-heidelberg.de)
 * @date   January, 2010
 *
 * @brief  This memory pool can be used for as a container for a limited number of "small" uniformingly sized objects.
 * A large boost is gained by saving a lot of malloc/free calls, as the objects will all be placed in this container.
 * The cache-hit percentage will also benefit from it because the container (and therefore its elements) will be continuous in the memory instead of scattered.
 * 
 * This class has a very specific use for ETS, a global object of class Pool is created to save temporary tuples created by pepsplice.
 * During its initialization, a sufficiently large block of memory is allocated via malloc.
 * To allow a class to use this container, <tt>void * Classname::operator new(size_t)</tt> has to be overloaded.
 * If more objects are concurrently placed into the container than slots are available, an exception is thrown, free'd slots will be re-used of course.
 * Memory Usage will be the size of the container itself + POOLSIZE * sizeof(int) to memorize the list of free slots.
 * Example Usage:
 * 
 * @code
 * Pool<Pepsplice::Tuple, POOLSIZE> __POOL__;
 *   void * Tuple::operator new(size_t)
 *  {
 *      return __POOL__.getNewObject();
 *  }
 *  void Tuple::operator delete(void * oldpointer)
 *  {
 *      __POOL__.deleteObject(oldpointer);
 *  }
 *
 * @endcode
 *
 */

template<typename T, size_t N>
class Pool
{

    size_t sType;
    void * memory;
    size_t * freeslot; //list
    size_t current_list_start;
public:

    /**
        Standard constructor, will allocate memory according to the template parameters, initialize the free slots.
    */
    Pool(){

        sType = sizeof(T);

        memory = malloc(N*sType);

        freeslot = (size_t*) malloc(N * sizeof(size_t) );
        if (freeslot == 0 || memory == 0)
        {
            throw "not enough contiguous free memory blocks";
        }
        reinit();
    }

    ~Pool()
    {
        free(memory);
        free(freeslot);
    }



    /**
        You get a pointer to the slot where you can put in your object, this should be the returnValue of your overloaded new operator.
        After claiming that memory chunk, it will get tagged as "used" via freeslot.
        According to specification, a void pointer has to be returned in the operator new(size_t) call, be careful as only char* can be used for calculation.
        @retval Address to memory chunk / object
    */
    void * getNewObject()
    {
        if (current_list_start == N) //No slot left
        {
            throw "Pool is already full, can't allocate";
        }
        size_t index = current_list_start; //current free index
        current_list_start = freeslot[index]; //next time, the next free entry in the list is returned.
        freeslot[index] = index;
        return (char*) memory + index * sType; //returns free place in the memory.
    }


    /**
        The equivalent to operator delete, given an address, it will tag the according chunk as "free".
        @param[in] mempointer *this
    */
    void deleteObject(void * mempointer)
    {
#ifdef DEBUG
        if (mempointer > memory + sType * N || mempointer < memory) throw "Pointer is bad";
#endif

        size_t index = ((char*)mempointer - (char*) memory) / sType;
        if (freeslot[index] != index) throw "Memory already deallocated, check your delete";
        freeslot[index] = current_list_start;

        current_list_start = index;
    }

    /**
        A reinitialization will render all objects in the container invalid, as all slots are tagged as free. Destructors will not be called.
    */
    void reinit(){
        current_list_start = 0;

        for (size_t i = 0; i < N; ++i)
        {
            freeslot[i] = i + 1;
        }
    }
    /**
        After a lot of new/delete operations, the free slots will be scattered, calling new repeatedly will not get you continuous chunks of memory.
        defrag will not move used slots, but link free slots in order according to their index / memory placement.
    */

    void defrag()
    {
        size_t oldindex = N;

        for(size_t i = N-1;i>=0;i--)

        {
            if (freeslot [i] == i) {continue;}

            freeslot[i] = oldindex;
            oldindex = i;
        }
        current_list_start = oldindex;
    }



};

