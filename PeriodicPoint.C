//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Define a Point in n dimensions WITH periodic boundary conditions     //
//                                                                      //
// Burkhard Militzer                                    Paris 4-23-99   //
// Modified class structure                         Livermore 8-08-02   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "PeriodicPoint.h"

template <int n> double PeriodicPoint<n>::boxSize[n];
template <int n> double PeriodicPoint<n>::halfBoxSize[n];
template <int n> double PeriodicPoint<n>::oneOverBoxSize[n];
template <int n> double PeriodicPoint<n>::maxDistance=0.0;
template <int n> double PeriodicPoint<n>::maxBoxSize=0.0;
template <int n> double PeriodicPoint<n>::minBoxSize=0.0;
template <int n> double PeriodicPoint<n>::volume=0.0;

template double PeriodicPoint<NDIM>::boxSize[NDIM];
template double PeriodicPoint<NDIM>::halfBoxSize[NDIM];
template double PeriodicPoint<NDIM>::oneOverBoxSize[NDIM];
template double PeriodicPoint<NDIM>::maxDistance;
template double PeriodicPoint<NDIM>::maxBoxSize;
template double PeriodicPoint<NDIM>::minBoxSize;
template double PeriodicPoint<NDIM>::volume;
