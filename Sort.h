//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
// Interface to apply the STL function sort for Array1                          //
//                                                                              //
// B. Militzer                                                 Hamburg, 7-17-06 //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#ifndef _SORT_
#define _SORT_

#include <algorithm>
#include "Vector.h"
#include "Array.h"

template <class T>
inline void SortArray(Array1 <T> & c) {
  int n = c.Size();
  Vector <T> v(n);
  for(int i=0; i<n; i++) v[i]=c[i];
  sort(v.begin(),v.end());
  for(int i=0; i<n; i++) c[i]=v[i];
}

template <class T>
inline void SortArrayReverse(Array1 <T> & c) {
  int n = c.Size();
  Vector <T> v(n);
  for(int i=0; i<n; i++) v[i]=c[i];
  sort(v.begin(),v.end());
  for(int i=0; i<n; i++) c[i]=v[n-1-i];
}

template <class T>
class Index {
public:
  Index():a(0),j(0){}
  const Array1 <T> * a;
  int j;
  bool operator< (const Index & i) const {
    if (a==0) error("Index sorting: Pointer not set");
    if (i.a==0) error("Index sorting: Pointer not set");
    return (*a)[j] < (*i.a)[i.j];
  }
};

template <class T>
inline Array1 <int> SortArrayIntoIndex(const Array1 <T> & ta) {
  int n = ta.Size();
  Vector < Index<T> > v(n);
  for(int i=0; i<n; i++) {
    v[i].j=i;
    v[i].a=&ta;
  }
  sort(v.begin(),v.end());
  Array1 <int> index(n);
  for(int i=0; i<n; i++) index[i]=v[i].j;
  return index;
}

template <class T>
inline Array1 <int> SortArrayIntoIndexReverse(const Array1 <T> & ta) {
  int n = ta.Size();
  Vector < Index<T> > v(n);
  for(int i=0; i<n; i++) {
    v[i].j=i;
    v[i].a=&ta;
  }
  sort(v.begin(),v.end());
  Array1 <int> index(n);
  for(int i=0; i<n; i++) index[n-1-i]=v[i].j;
  return index;
}
///////////////////////////////////////////////////////////////////////////////////

template <class T, class C>
inline void SortArrayGen(Array1 <T> & c) {
  int n = c.Size();
  Vector <T> v(n);
  for(int i=0; i<n; i++) v[i]=c[i];
  sort(v.begin(),v.end(),C());
  for(int i=0; i<n; i++) c[i]=v[i];
}

template <class T, class C>
class IndexGen {
public:
  IndexGen():a(0),j(0){}
  const Array1 <T> * a;
  int j;
  bool operator< (const IndexGen & i) const {
    if (a==0) error("Index sorting: Pointer not set");
    if (i.a==0) error("Index sorting: Pointer not set");
    return C((*a)[j],(*i.a)[i.j]);
  }
};

template <class T, class C>
inline Array1 <int> SortArrayIntoIndexGen(const Array1 <T> & ta) {
  int n = ta.Size();
  Vector < IndexGen<T,C> > v(n);
  for(int i=0; i<n; i++) {
    v[i].j=i;
    v[i].a=&ta;
  }
  sort(v.begin(),v.end(),C());
  Array1 <int> index(n);
  for(int i=0; i<n; i++) index[i]=v[i].j;
  return index;
}

#endif // _SORT_
