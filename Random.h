/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Interface to different RANDOM NUMBER GENERATORS                             //
// Make generators interchangeable                                             //                 
//                                                                             //
// Burkhard Militzer                                    Urbana 1998            //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

#ifndef _RANDOM_
#define _RANDOM_

#include "Standard.h"

#ifdef USE_SPRNG

// remember, the MPI version is SPRNG must be linked if mpi.h is included here.
// otherwise strange link error
#ifdef USE_MPI
#include <mpi.h>
#endif

// #define SIMPLE_SPRNG // use the simpler SPRNG interface (one stream only)
#include "sprng.h"          

#ifndef SIMPLE_SPRNG
extern int* sprngStream;
#endif 

inline double Random() { 
#ifdef SIMPLE_SPRNG
  return sprng();
#else
  return sprng(sprngStream);
#endif
}

#ifdef SIMPLE_SPRNG
inline void InitRandom(const bool randomSeedFlag) {
  if (randomSeedFlag) {
    int seed = make_sprng_seed();   
    //    init_sprng(seed,SPRNG_DEFAULT);                      // SPRNG version 1.0
    init_sprng(DEFAULT_RNG_TYPE,seed,SPRNG_DEFAULT);  // SPRNG version 2.0
  } else {
    warning("SIMPLE_SPRNG: Check if running without init_sprng() is ok!!!");
  }
}

#else

inline void InitRandom(const bool randomSeedFlag, const int stream=0) {
  int seed;
  if (randomSeedFlag) {
    seed = make_sprng_seed(); // make new seed depending of date and time
  } else {
    //    seed = 985456376;  // Ashok's default seed
    seed = 0;                // Still will reproduce the "SIMPLE_SPRING" results without explicite initialization
  }
  int gType = 0;
  sprngStream = init_sprng(gType,stream,stream+1,seed,SPRNG_DEFAULT); // initialize stream 
}
#endif

#else // USE_SPRNG 

#ifdef USE_CRANF32
#include "cranf32.h"
inline double Random() { 
  return getranf();
}
inline void InitRandom(const bool randomSeedFlag, const int stream=0){
  if (randomSeedFlag || stream>0) error("No implemented for CRANF32");
  int seed[]={12345,54321,-1};
  setseed(seed);
};

#else // USE_CRANF32

#include <stdlib.h>
#include <unistd.h> 
#include <sys/times.h>
#include <sys/time.h>
inline void InitRandom(const bool randomSeedFlag, const int stream=0) {
  if (stream>0) error("No implemented for DRAND48");

  // added to run properly on G5 using gcc 3.3
  unsigned short seed16v[3]={0,0,0};
  //  unsigned short *j;

  if (randomSeedFlag) {
    struct timeval val;
    struct timezone zone;
    gettimeofday(&val,&zone);
    long t = val.tv_sec*1000000 + val.tv_usec;
    unsigned short* tt = (unsigned short*) &t;
    seed16v[0] = tt[0];
    seed16v[1] = tt[1];
    seed16v[2] = tt[2];
    //    srand (time (NULL));
  }

  //  j = seed48(seed16v);
  seed48(seed16v);
}

inline double Random() { 
  return drand48();
}

#endif // USE_CRANF32
#endif // USE_SPRNG 
#endif // _RANDOM_
