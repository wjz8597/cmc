//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Define a Point in n dimensions WITH periodic boundary conditions     //
//                                                                      //
// Burkhard Militzer                                    Paris 4-23-99   //
// Modified class structure                         Livermore 8-08-02   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "PeriodicPointGen.h"

template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::a;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::b;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::c;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::aMin;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::bMin;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::cMin;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::aa;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::bb;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::cc;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::aaMin;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::bbMin;
template<int DIM> PeriodicPointGen<DIM> PeriodicPointGen<DIM>::ccMin;
template<int DIM> double PeriodicPointGen<DIM>::volume;
template<int DIM> bool PeriodicPointGen<DIM>::isentropic;
template<int DIM> bool PeriodicPointGen<DIM>::orthogonal;

template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::a;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::b;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::c;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::aMin;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::bMin;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::cMin;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::aa;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::bb;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::cc;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::aaMin;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::bbMin;
template PeriodicPointGen<NDIM> PeriodicPointGen<NDIM>::ccMin;
template double PeriodicPointGen<NDIM>::volume;
template bool PeriodicPointGen<NDIM>::isentropic;
template bool PeriodicPointGen<NDIM>::orthogonal;
