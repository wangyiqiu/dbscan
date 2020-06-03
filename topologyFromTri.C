// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include "parallel.h"
#include "gettime.h"
#include "geometry.h"
#include "topology.h"
#include "deterministicHash.h"
using namespace std;

typedef pair<intT, intT> pairInt;
typedef pair<pairInt, tri*> edge;

// Hash table to store skinny triangles
struct hashEdges {
  typedef pairInt kType;
  typedef edge* eType;
  eType empty() {return NULL;}
  kType getKey(eType v) { return v->first;}
  uintT hash(kType s) { return utils::hash(s.first)+3*(utils::hash(s.second)); }
  int cmp(kType s1, kType s2) {
    return ((s1.first > s2.first) ? 1 : 
	    (s1.first < s2.first) ? -1 : 
	    (s1.second > s2.second) ? 1 :
	    (s1.second < s2.second) ? -1 : 0);
  }
  bool replaceQ(eType s, eType s2) {return 0;}
};

typedef Table<hashEdges,intT> EdgeTable;
EdgeTable makeEdgeTable(intT m) {return EdgeTable(m,hashEdges());}

void topologyFromTriangles(triangles<point2d> Tri, vertex** vr, tri** tr) {
  intT n = Tri.numPoints;
  point2d* P = Tri.P;

  intT m = Tri.numTriangles;
  triangle* T = Tri.T;

  if (*vr == NULL) *vr = newA(vertex,n);
  vertex* v = *vr;
  parallel_for (intT i=0; i < n; i++)
    v[i] = vertex(P[i],i);

  if (*tr == NULL) *tr = newA(tri,m);
  tri* Triangs = *tr;
  edge* E = newA(edge, m*3);
  EdgeTable ET = makeEdgeTable(m*6);
  parallel_for (intT i=0; i < m; i++)
    for (int j=0; j<3; j++) {
      E[i*3 + j] = edge(pairInt(T[i].C[j], T[i].C[(j+1)%3]), &Triangs[i]);
      ET.insert(&E[i*3+j]);
      Triangs[i].vtx[(j+2)%3] = &v[T[i].C[j]];
    }

  parallel_for (intT i=0; i < m; i++) {
    Triangs[i].id = i;
    Triangs[i].initialized = 1;
    Triangs[i].bad = 0;
    for (int j=0; j<3; j++) {
      pairInt key = pairInt(T[i].C[(j+1)%3], T[i].C[j]);
      edge *Ed = ET.find(key);
      if (Ed != NULL) Triangs[i].ngh[j] = Ed->second;
      else { Triangs[i].ngh[j] = NULL;
	//Triangs[i].vtx[j]->boundary = 1;
	//Triangs[i].vtx[(j+2)%3]->boundary = 1;
      }
    }
  }
  
  ET.del();
  free(E);
}

// Note that this is not currently a complete test of correctness
// For example it would allow a set of disconnected triangles, or even no
// triangles
bool checkDelaunay(tri *triangs, intT n, intT boundarySize) {
  intT *bcount = newA(intT,n);
  parallel_for (intT j=0; j<n; j++) bcount[j] = 0;
  intT insideOutError = n;
  intT inCircleError = n;
  parallel_for (intT i=0; i<n; i++) {
    if (triangs[i].initialized) {
      simplex t = simplex(&triangs[i],0);
      for (int j=0; j < 3; j++) {
	simplex a = t.across();
	if (a.valid()) {
	  vertex* v = a.rotClockwise().firstVertex();

          // Check that the neighbor is outside the triangle
	  if (!t.outside(v)) {
	    double vz = triAreaNormalized(t.t->vtx[(t.o+2)%3]->pt, 
					  v->pt, t.t->vtx[t.o]->pt);
	    //cout << "i=" << i << " vz=" << vz << endl;
	    // allow for small error
	    if (vz < -1e-10)  utils::writeMin(&insideOutError,i);
	  }

          // Check that the neighbor is not in circumcircle of the triangle
	  if (t.inCirc(v)) {
	    double vz = inCircleNormalized(t.t->vtx[0]->pt, t.t->vtx[1]->pt, 
					  t.t->vtx[2]->pt, v->pt);
	    //cout << "i=" << i << " vz=" << vz << endl;
	    // allow for small error
	    if (vz > 1e-10) utils::writeMin(&inCircleError,i);
	  }
	} else bcount[i]++;
	t = t.rotClockwise();
      }
    }
  }
  intT cnt = sequence::plusReduce(bcount,n);
  //if (boundarySize != cnt) {
  //cout << "delaunayCheck: wrong boundary size, should be " << boundarySize 
  //<< " is " << cnt << endl;
  //return 1;
  //}
  free(bcount);

  if (insideOutError < n) {
    cout << "delaunayCheck: neighbor inside triangle at triangle " 
	 << inCircleError << endl;
    return 1;
  }
  if (inCircleError < n) {
    cout << "In Circle Violation at triangle " << inCircleError << endl;
    return 1;
  }

  return 0;
}
