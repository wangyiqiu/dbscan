///////////////// CONFIGURATION ///////////////////
// Enable one out of four to select an algorithm
#define OUR_EXACT
// #define OUR_EXACT_QT
// #define OUR_APPROX
// #define OUR_APPROX_QT

// Enable or disable bucketing optimization
#define USE_BUCKETING
//////////////////////////////////////////////////

// Calling example for exact DBSCAN
// ./DBSCAN -eps 10 -minpts 10 -o clusters.txt dataset.txt

// Calling example for approximate DBSCAN
// ./DBSCAN -eps 10 -minpts 10 -rho 0.1 -o clusters.txt dataset.txt

/* Input dataset example, 3 dimensional dataset of 5 points
dim 3
10 20 30
5 2 1
3 3 5
2 20 10
22 10 30
*/

/* Output file meaning
<point index (starting 0)> <cluster id>
same cluster id means belong to the same cluster
*/

/* For best performance (benchmarking)
 - disable LONG_VARIABLE
 - enable USE_JEMALLOC (see more details on github page)
 in DBSCAN.h
 */

#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <numeric>
#include "DBSCANGeometry.h"
#include "geometry.h"
#include "DBSCAN.h"
#include "GetNeighborKDTree.h"
#include "Grids.h"
#include "unionFind.h"
#include "QuadTree.h"

//////////////// PREPROCESSORS ///////////////////
#ifdef USE_JEMALLOC
#include<jemalloc/jemalloc.h>
#include "je_allocator.h"
#endif

#if defined(OUR_EXACT)
#endif

#if defined(OUR_EXACT_QT)
#define ALLPT_QUADTREE
#endif

#if defined(OUR_APPROX)
#define APPROX
#define COREPT_QUADTREE
#endif

#if defined(OUR_APPROX_QT)
#define APPROX
#define ALLPT_QUADTREE
#define COREPT_QUADTREE
#endif

#if defined(ALLPT_QUADTREE)
#define QUADTREE_MARKCORE
#endif

#if defined(COREPT_QUADTREE)
#define QUADTREE_BCP
#endif

#if defined(ALLPT_QUADTREE) || defined(COREPT_QUADTREE) || defined(APPROX)

#if defined(LONG_VARIABLE)
//#error "long variables not yet supported for quadtree, turn off in DBSCAN.h"
#endif

#if defined(FLOAT_DATATYPE)
#error "quadtree does not support float datatype"
#endif

#endif

#if defined(ALLPT_QUADTREEL)
#error "all point quadtree under maintenance"
#endif

// for bcp grid ordering heuristic (bucketing)
struct gridDegComparator {
  Grids *m_grids;
  gridDegComparator(Grids *t_grids) : m_grids(t_grids) {}
  inline intT operator() (const GdLongDataType& t_g1, const GdLongDataType& t_g2) {
    return m_grids->GetNumPts(t_g1) > m_grids->GetNumPts(t_g2); // no obvious improvement
  }
};

struct pointDistComparator {
  pointDistComparator() {}
  inline intT operator() (const pair<PtIntermediateDataType, PtDataType*>& i, const pair<PtIntermediateDataType, PtDataType*>& j) {
    return i.first < j.first;
  }
};

inline bool PointGridDistanceLeq(intT dim, Grids *t_grids, GdLongDataType g, PtLongDataType p, PtIntermediateDataType t_threshold) {

  OtLongDataType gridDimNum[dim];
  t_grids->GetGridVec(g, gridDimNum);
  PtIntermediateDataType result = 0;
  PtIntermediateDataType dimSmall, dimLarge, pCoord, dimDist;
  for (intT d = 0; d < dim; ++ d) {
    dimSmall = gridDimNum[d] * t_grids->GetSize();
    dimLarge = (gridDimNum[d] + 1) * t_grids->GetSize();
    pCoord = t_grids->GetPtData(p)[d] - t_grids->GetLowCoord(d);
    dimDist = 0;

    if (pCoord >= dimSmall && pCoord <= dimLarge) {
      dimDist = 0;
    } else if (pCoord < dimSmall) {
      dimDist = dimSmall - pCoord;
    } else {
      dimDist = pCoord - dimLarge;
    }
    
    result += dimDist * dimDist;
  }
  
  if (result > t_threshold) {
    return false;
  }
  return true;
}

inline PtIntermediateDataType PointGridDistance(intT dim, Grids *t_grids, GdLongDataType g, PtLongDataType p) {

  OtLongDataType gridDimNum[dim];
  t_grids->GetGridVec(g, gridDimNum);
  PtIntermediateDataType result = 0;
  PtIntermediateDataType dimSmall, dimLarge, pCoord, dimDist;
  for (intT d = 0; d < dim; ++ d) {
    dimSmall = gridDimNum[d] * t_grids->GetSize();
    dimLarge = (gridDimNum[d] + 1) * t_grids->GetSize();
    pCoord = t_grids->GetPtData(p)[d] - t_grids->GetLowCoord(d);
    dimDist = 0;

    if (pCoord >= dimSmall && pCoord <= dimLarge) {
      dimDist = 0;
    } else if (pCoord < dimSmall) {
      dimDist = dimSmall - pCoord;
    } else {
      dimDist = pCoord - dimLarge;
    }
    
    result += dimDist * dimDist;
  }
  
  return result;
}

pair<PtLongDataType, intT *> MarkCorePointGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *t_isCoreGrid) {
  intT *corePointMarkers = newA(intT,t_points->n);
  parallel_for (PtLongDataType p = 0; p < t_points->n; ++ p) {
    corePointMarkers[p] = 0;
  }

  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    t_isCoreGrid[g] = false;
    GdPtLongDataType gNumPts = t_grids->GetNumPts(g);

    if (gNumPts >= t_params->m_minPts) {

      t_isCoreGrid[g] = true;
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t (p, t_grids->GetStartOffset(g), iterEnd, 2000, PtLongDataType, {
        corePointMarkers[t_grids->GetPt(p)] = 1;
      });

    } else {
      
      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      
      // compute potential more near points other than own grid by summing up that of
      // every neighboring grid, can just skip the loop when insufficient
      PtLongDataType potentialMoreNearPoints = 0;
      for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
        GdLongDataType n = neighbors->at(ni);
        if (n < 0 || n == g) continue;
        potentialMoreNearPoints += t_grids->GetNumPts(n);
      }

      if (potentialMoreNearPoints >= t_params->m_minPts - gNumPts) {
        PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
        parallel_for_1 (PtLongDataType p = t_grids->GetStartOffset(g); p < iterEnd; ++p) {
          PtLongDataType totalNearPoints = gNumPts;
          for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
            GdLongDataType n = neighbors->at(ni);
            if (n < 0 || n == g) continue;

            GdPtLongDataType nNumPts = t_grids->GetNumPts(n);

            PtLongDataType iterEndInner = t_grids->GetStartOffset(n) + nNumPts;
            for (PtLongDataType np = t_grids->GetStartOffset(n); np < iterEndInner; ++np) {
              if (PointDistNDSqrLeq(t_grids->GetPtData(p), t_grids->GetPtData(np), dim, t_params->m_epsilon * t_params->m_epsilon)) { // todo opt comp cast
                totalNearPoints ++;
                if (totalNearPoints >= t_params->m_minPts) {
                  goto decide_core;
                }
              }
            }
          }
          decide_core:
          if (totalNearPoints >= t_params->m_minPts) {
            corePointMarkers[t_grids->GetPt(p)] = 1;
            t_isCoreGrid[g] = true;
          } else {
            corePointMarkers[t_grids->GetPt(p)] = 0;
          }
        } // for pts in g
      } // if has potentially enough points

    }
  }

  PtLongDataType numCore;
  if (sizeof(PtLongDataType) > 4) {
    numCore = sequence::reduceIL(corePointMarkers, (PtLongDataType)t_points->n, utils::addF<intT>());
  } else {
    numCore = sequence::reduce(corePointMarkers, (PtLongDataType)t_points->n, utils::addF<intT>());
  }
  return make_pair(numCore, corePointMarkers);
}

#ifdef QUADTREE_H
pair<PtLongDataType, intT *> MarkCorePointGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, QuadTree *t_allPtQuadTree, bool *t_isCoreGrid) {

  intT *corePointMarkers = newA(intT,t_points->n);
  parallel_for (PtLongDataType p = 0; p < t_points->n; ++ p) {
    corePointMarkers[p] = 0;
  }

  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    t_isCoreGrid[g] = false;
    GdPtLongDataType gNumPts = t_grids->GetNumPts(g);

    // If core grid, just mark all.
    if (gNumPts >= t_params->m_minPts) {

      t_isCoreGrid[g] = true;
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t (p, t_grids->GetStartOffset(g), iterEnd, 100, PtLongDataType, {
        corePointMarkers[t_grids->GetPt(p)] = 1;
      });

    } else {
      
      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);

      PtLongDataType potentialMoreNearPoints = 0;
      for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
        GdLongDataType n = neighbors->at(ni);
        if (n < 0 || n == g) continue;
        potentialMoreNearPoints += t_grids->GetNumPts(n);
      }
      
      if (potentialMoreNearPoints >= t_params->m_minPts - gNumPts) {
        // Iterate point p in grid.
        PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
        parallel_for_1 (PtLongDataType p = t_grids->GetStartOffset(g); p < iterEnd; ++p) {
          PtLongDataType totalNearPoints = gNumPts;
          
  	      for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
            GdLongDataType n = neighbors->at(ni);
            if (n < 0 || n == g) continue;
            GdPtLongDataType nNumPts = t_grids->GetNumPts(n);
        	  // instead of itering points in neighbor, do a range query
        	  totalNearPoints += t_allPtQuadTree->RangeQuery(pointsData, pointsData + dim * t_grids->GetPt(p), dim, t_grids, n, t_params->m_epsilon * t_params->m_epsilon);
        	  if (totalNearPoints >= t_params->m_minPts) {
        	    goto decide_core;
        	  }	  
          }
          decide_core:
          if (totalNearPoints >= t_params->m_minPts) {
            corePointMarkers[t_grids->GetPt(p)] = 1;
            // numCore ++;
            t_isCoreGrid[g] = true;
          } else {
            corePointMarkers[t_grids->GetPt(p)] = 0;
          }
        } // for pts in g
      } // if has potentially enough points
    }
  }
  
  PtLongDataType numCore = sequence::reduce(corePointMarkers, t_points->n, utils::addF<intT>());
  return make_pair(numCore, corePointMarkers);
}

// quadtree version
bool edgeAddableBCPND(QuadTree *t_corePtQuadTree, PtDataType *pointsData, intT t_dim, DbscanParams *t_params, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  GdPtLongDataType g1NumPts = t_corePtQuadTree->GetGridNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_corePtQuadTree->GetGridNumPts(t_g2);

  if (g1NumPts == 0 || g2NumPts == 0) {
    return false;
  }

  // iter g1 core pts
  for (PtLongDataType p = t_corePtQuadTree->GetGridPtStartOffset(t_g1); p < t_corePtQuadTree->GetGridPtEndOffset(t_g1); ++ p) {
    PtLongDataType myPt = t_corePtQuadTree->GetPt(p); // point id in global pt array
    // exact binary range query on g2
    if(t_corePtQuadTree->BinRangeQuery(pointsData, pointsData + myPt * t_dim, t_dim, t_grids, t_g2, t_params->m_epsilon * t_params->m_epsilon)) {
      return true;
    }
  }
  
  return false;
}

#endif // ifdef quadtree h

// radius optimized version
bool edgeAddableBCPND(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  // t_g1 and t_g2 are compact-grid ids.
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);

  // PtDataType **g1Core = newA(PtDataType *,g1NumPts); // fuse alloc, and use jemalloc optionally
  // PtDataType **g2Core = newA(PtDataType*, g2NumPts);
#if defined(USE_JEMALLOC)
  PtDataType **g1Core = (PtDataType **)JE_MALLOC_FUNC(sizeof(PtDataType*) * (g1NumPts+g2NumPts));
  PtDataType **g2Core = g1Core + g1NumPts;
#else
  PtDataType **g1Core = (PtDataType **)malloc(sizeof(PtDataType*) * (g1NumPts+g2NumPts));
  PtDataType **g2Core = g1Core + g1NumPts;
#endif
  //PtIntermediateDataType threshold = static_cast<PtIntermediateDataType>(t_params->m_epsilon) * static_cast<PtIntermediateDataType>(t_params->m_epsilon);

  // Iterate through both lists, and eliminate points that are further than r from the other grid.
  // (This way there is no way the point being close to a point).
  
  GdPtLongDataType num1 = 0;
  PtLongDataType iterEnd = t_grids->GetStartOffset(t_g1) + g1NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g1); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      if (PointGridDistanceLeq(dim, t_grids, t_g2, p, t_radiusSqr1)) {
        g1Core[num1 ++] = t_grids->GetPtData(p);
      }
    }
  }
  
  GdPtLongDataType num2 = 0;
  iterEnd = t_grids->GetStartOffset(t_g2) + g2NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g2); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      //PtIntermediateDataType distToGridSqr = NDPointGridDistanceSqr(pointsData, dim, t_grids, t_g1, p);
      //if (distToGridSqr <= threshold) {
      if (PointGridDistanceLeq(dim, t_grids, t_g1, p, t_radiusSqr1)) {
        g2Core[num2 ++] = t_grids->GetPtData(p);
      }
    }
  }

  for (GdPtLongDataType i = 0; i < num1; ++ i) {
    for (GdPtLongDataType j = 0; j < num2; ++ j) {
      if (PointDistNDSqrLeq(g1Core[i], g2Core[j], dim, t_radiusSqr1)) {
	// free(g1Core);
	// free(g2Core);
#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
	return true;
      }
    }
  }
  
  // free(g1Core);
  // free(g2Core);
#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
  return false;

}

// radius optimized version with distance ordering
bool edgeAddableBCPNDWithSort(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  // t_g1 and t_g2 are compact-grid ids.
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);

#if defined(USE_JEMALLOC)
  pair<PtIntermediateDataType, PtDataType*> *g1Core = (pair<PtIntermediateDataType, PtDataType*>*)JE_MALLOC_FUNC(sizeof(pair<PtIntermediateDataType, PtDataType*>) * (g1NumPts+g2NumPts));
  pair<PtIntermediateDataType, PtDataType*>*g2Core = g1Core + g1NumPts;
#else
  pair<PtIntermediateDataType, PtDataType*>*g1Core = (pair<PtIntermediateDataType, PtDataType*>*)malloc(sizeof(pair<PtIntermediateDataType, PtDataType*>) * (g1NumPts+g2NumPts));
  pair<PtIntermediateDataType, PtDataType*>*g2Core = g1Core + g1NumPts;
#endif
  // Iterate through both lists, and eliminate points that are further than r from the other grid.
  // (This way there is no way the point being close to a point).
  
  GdPtLongDataType num1 = 0;
  PtLongDataType iterEnd = t_grids->GetStartOffset(t_g1) + g1NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g1); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      PtIntermediateDataType distToGridSqr = PointGridDistance(dim, t_grids, t_g2, p);
      if (distToGridSqr <= t_radiusSqr1) {
        g1Core[num1 ++] = pair<PtIntermediateDataType, PtDataType*>(distToGridSqr, t_grids->GetPtData(p));
      }
    }
  }
  
  GdPtLongDataType num2 = 0;
  iterEnd = t_grids->GetStartOffset(t_g2) + g2NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g2); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      PtIntermediateDataType distToGridSqr = PointGridDistance(dim, t_grids, t_g1, p);
      if (distToGridSqr <= t_radiusSqr1) {
        g2Core[num2 ++] = pair<PtIntermediateDataType, PtDataType*>(distToGridSqr, t_grids->GetPtData(p));
      }
    }
  }


  sampleSort(g1Core, num1, pointDistComparator());
  sampleSort(g2Core, num2, pointDistComparator());
  GdPtLongDataType pos1 = 0;
  GdPtLongDataType pos2 = 0;
  while (pos1 < num1 && pos2 < num2) {
    if (g1Core[pos1].first <= g2Core[pos2].first) {
      for (GdPtLongDataType i = pos2; i < num2; i++) {
        if (PointDistNDSqrLeq(g1Core[pos1].second, g2Core[pos2].second, dim, t_radiusSqr1)) {
#if defined(USE_JEMALLOC)
          JE_MALLOC_FREE(g1Core);
#else
          free(g1Core);
#endif
          return true;
        }
      }
      pos1++;
    } else {
      for (GdPtLongDataType i = pos1; i < num1; i++) {
        if (PointDistNDSqrLeq(g1Core[pos1].second, g2Core[pos2].second, dim, t_radiusSqr1)) {
#if defined(USE_JEMALLOC)
          JE_MALLOC_FREE(g1Core);
#else
          free(g1Core);
#endif
          return true;
        }
      }
      pos2++;
    }
  }
  
#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
  return false;

}

// radius optimized version with distance ordering only on grid g1 (imbalanced points)
bool edgeAddableBCPNDWithG1Sort(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  // t_g1 and t_g2 are compact-grid ids.
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);

#if defined(USE_JEMALLOC)
  pair<PtIntermediateDataType, PtDataType*> *g1Core = (pair<PtIntermediateDataType, PtDataType*>*)JE_MALLOC_FUNC(sizeof(pair<PtIntermediateDataType, PtDataType*>) * (g1NumPts+g2NumPts));
  pair<PtIntermediateDataType, PtDataType*>*g2Core = g1Core + g1NumPts;
#else
  pair<PtIntermediateDataType, PtDataType*>*g1Core = (pair<PtIntermediateDataType, PtDataType*>*)malloc(sizeof(pair<PtIntermediateDataType, PtDataType*>) * (g1NumPts+g2NumPts));
  pair<PtIntermediateDataType, PtDataType*>*g2Core = g1Core + g1NumPts;
#endif
  // Iterate through both lists, and eliminate points that are further than r from the other grid.
  // (This way there is no way the point being close to a point).
  
  GdPtLongDataType num1 = 0;
  PtLongDataType iterEnd = t_grids->GetStartOffset(t_g1) + g1NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g1); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      PtIntermediateDataType distToGridSqr = PointGridDistance(dim, t_grids, t_g2, p);
      if (distToGridSqr <= t_radiusSqr1) {
        g1Core[num1 ++] = pair<PtIntermediateDataType, PtDataType*>(distToGridSqr, t_grids->GetPtData(p));
      }
    }
  }
  
  GdPtLongDataType num2 = 0;
  iterEnd = t_grids->GetStartOffset(t_g2) + g2NumPts;
  for (PtLongDataType p = t_grids->GetStartOffset(t_g2); p < iterEnd; ++p) {
    if (PointGridDistanceLeq(dim, t_grids, t_g1, p, t_radiusSqr1)) {
      g2Core[num2 ++] = pair<PtIntermediateDataType, PtDataType*>(-1, t_grids->GetPtData(p));
    }
  }

  sampleSort(g1Core, num1, pointDistComparator());

  for (GdPtLongDataType i = 0; i < num1; ++ i) {
    for (GdPtLongDataType j = 0; j < num2; ++ j) {
      if (PointDistNDSqrLeq(g1Core[i].second, g2Core[j].second, dim, t_radiusSqr1)) {
#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
  return true;
      }
    }
  }
  
#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
  return false;

}

// simplest possible version, suitable for smaller number of points, say both < 100
bool edgeAddableBCPNDSimple(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);
  PtLongDataType end1 = t_grids->GetStartOffset(t_g1) + g1NumPts;
  PtLongDataType end2 = t_grids->GetStartOffset(t_g2) + g2NumPts;
  for(PtLongDataType p1 = t_grids->GetStartOffset(t_g1); p1 < end1; ++p1) {
    for(PtLongDataType p2 = t_grids->GetStartOffset(t_g2); p2 < end2; ++p2) {
      if (t_corePointMarkers[t_grids->GetPt(p1)] && t_corePointMarkers[t_grids->GetPt(p2)] && PointDistNDSqrLeq(t_grids->GetPtData(p1), t_grids->GetPtData(p2), dim, t_radiusSqr1)) {
	return true;
      }
    }
  }
  return false;
}

inline GdPtLongDataType GetCorePtsParallel(Grids *t_grids, intT *t_corePointMarkers, PtDataType **g1Core, GdLongDataType t_g1, GdLongDataType t_g2, PtIntermediateDataType t_radiusSqr1, intT dim) {
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);

#if defined(USE_JEMALLOC)
  PtLongDataType *coreArray = (PtLongDataType *)JE_MALLOC_FUNC(sizeof(PtLongDataType) * g1NumPts * 2);
#else
  PtLongDataType *coreArray = newA(PtLongDataType, g1NumPts * 2);
#endif
  PtLongDataType *g1CoreFlag = coreArray;
  PtLongDataType *g1CoreRaw = coreArray + g1NumPts;

  parallel_for(GdPtLongDataType p = 0; p < g1NumPts; ++p) {
    g1CoreFlag[p] = 0;
  }
  parallel_for(PtLongDataType p = t_grids->GetStartOffset(t_g1);
      p < t_grids->GetStartOffset(t_g1) + g1NumPts; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
      if (PointGridDistanceLeq(dim, t_grids, t_g2, p, t_radiusSqr1)) {
        g1CoreRaw[p - t_grids->GetStartOffset(t_g1)] = p;
        g1CoreFlag[p - t_grids->GetStartOffset(t_g1)] = 1;
      }
    }
  }
  
  GdPtLongDataType num1 = -1;
  if (sizeof(GdPtLongDataType) > 4) {
    num1 = (GdPtLongDataType)sequence::prefixSumL<PtLongDataType>(g1CoreFlag, 0, g1NumPts);
  } else {
    num1 = (GdPtLongDataType)sequence::prefixSum<PtLongDataType>(g1CoreFlag, 0, g1NumPts);
  }
  parallel_for(GdPtLongDataType i = 0; i < g1NumPts - 1; ++ i) {
    if (g1CoreFlag[i] != g1CoreFlag[i + 1]) {
      g1Core[g1CoreFlag[i]] = t_grids->GetPtData(g1CoreRaw[i]);
    }
  }
  GdPtLongDataType i = g1NumPts - 1;
  if (g1CoreFlag[i] != num1) {
    g1Core[g1CoreFlag[i]] = t_grids->GetPtData(g1CoreRaw[i]);
  }

#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(coreArray);
#else
  free(coreArray);
#endif
  return num1;
}

// radius optimized version - parallel
bool edgeAddableBCPNDParallel(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {

  const OtLongDataType parConst1 = 5000;
  const OtLongDataType parConst2 = 50000;
  const OtLongDataType grain = 1000;

  // t_g1 and t_g2 are compact-grid ids.
  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);

#if defined(USE_JEMALLOC)
  PtDataType **g1Core = (PtDataType **)JE_MALLOC_FUNC(sizeof(PtDataType*) * (g1NumPts+g2NumPts));
#else
  PtDataType **g1Core = (PtDataType **)malloc(sizeof(PtDataType*) * (g1NumPts+g2NumPts));
#endif
  PtDataType **g2Core = g1Core + g1NumPts;

  GdPtLongDataType num1 = 0;

  if (g1NumPts < parConst1) {
    PtLongDataType iterEnd = t_grids->GetStartOffset(t_g1) + g1NumPts;
    for (PtLongDataType p = t_grids->GetStartOffset(t_g1); p < iterEnd; ++p) {
      if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
        if (PointGridDistanceLeq(dim, t_grids, t_g2, p, t_radiusSqr1)) {
          g1Core[num1 ++] = t_grids->GetPtData(p);
        }
      }
    }
  } else {
    num1 = GetCorePtsParallel(t_grids, t_corePointMarkers, g1Core, t_g1, t_g2, t_radiusSqr1, dim);
  }

  if (num1 == 0) {
#if defined(USE_JEMALLOC)
    JE_MALLOC_FREE(g1Core);
#else
    free(g1Core);
#endif
    return false;
  }

  GdPtLongDataType num2 = 0;

  if (g2NumPts < parConst1) {
    PtLongDataType iterEnd = t_grids->GetStartOffset(t_g2) + g2NumPts;
    for (PtLongDataType p = t_grids->GetStartOffset(t_g2); p < iterEnd; ++p) {
      if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {
        if (PointGridDistanceLeq(dim, t_grids, t_g1, p, t_radiusSqr1)) {
          g2Core[num2 ++] = t_grids->GetPtData(p);
        }
      }
    }
  } else {
    num2 = GetCorePtsParallel(t_grids, t_corePointMarkers, g2Core, t_g2, t_g1, t_radiusSqr1, dim);
  }
  
  if (num2 == 0) {
#if defined(USE_JEMALLOC)
    JE_MALLOC_FREE(g1Core);
#else
    free(g1Core);
#endif
    return false;
  }

  /*
#if defined(USE_JEMALLOC)
    JE_MALLOC_FREE(g1Core);
#else
    free(g1Core);
#endif
    return false;
  */

  if (num1 < parConst2 || num2 < parConst2) {
    for (GdPtLongDataType i = 0; i < num1; ++ i) {
      for (GdPtLongDataType j = 0; j < num2; ++ j) {
        if (PointDistNDSqrLeq(g1Core[i], g2Core[j], dim, t_radiusSqr1)) {
  #if defined(USE_JEMALLOC)
          JE_MALLOC_FREE(g1Core);
  #else
          free(g1Core);
  #endif
          return true;
        }
      }
    }
  } else {
    GdPtLongDataType numIters = num1 / grain;
    GdPtLongDataType numRemain = num1 % grain;
    GdPtLongDataType toAlloc = -1;
    if (numIters == 0) {
      toAlloc = numRemain;
    } else {
      toAlloc = grain;
    }
    GdPtLongDataType numIters2 = num2 / grain;
    GdPtLongDataType numRemain2 = num2 % grain;
    GdPtLongDataType toAlloc2 = -1;
    if (numIters2 == 0) {
      toAlloc2 = numRemain2;
    } else {
      toAlloc2 = grain;
    }
    /*
    cout << "numIters1 " << numIters << endl;
    cout << "numRemain " << numRemain << endl;
    cout << "numIters2 " << numIters2 << endl;
    cout << "numRemain2 " << numRemain2 << endl;
    */
#if defined(USE_JEMALLOC)
    PtIntermediateDataType *distancesSqr = (PtIntermediateDataType *)JE_MALLOC_FUNC(sizeof(PtIntermediateDataType) * toAlloc * toAlloc2);
#else
    PtIntermediateDataType *distancesSqr = (PtIntermediateDataType *)malloc(sizeof(PtIntermediateDataType) * toAlloc * toAlloc2);
#endif
    //cout << "alloc " << toAlloc * toAlloc2 << endl;

    GdPtLongDataType i = 0;
    for (; i < numIters; ++ i) {
      GdPtLongDataType j = 0;
      for (; j < numIters2; ++ j) {
	parallel_for_1(GdPtLongDataType ii = 0; ii < grain; ++ ii) {
	  granular_for_t(jj, 0, grain, 3000, GdPtLongDataType, {
	      distancesSqr[ii * grain + jj] = PointDistNDSqr(g1Core[i * grain + ii], g2Core[j * grain + jj], dim);
	    });
	}
	PtIntermediateDataType minDistSqr = sequence::reduce(distancesSqr, (OtLongDataType) grain * grain, utils::minF<PtIntermediateDataType>());
	if (minDistSqr <= t_radiusSqr1) {
#if defined(USE_JEMALLOC)
	  JE_MALLOC_FREE(distancesSqr);
	  JE_MALLOC_FREE(g1Core);
#else
	  free(distancesSqr);
	  free(g1Core);
#endif
	  return true;
	}
      }

      if (numRemain2 > 0) {
	parallel_for_1(GdPtLongDataType ii = 0; ii < grain; ++ ii) {
	  granular_for_t(jj, 0, numRemain2, 3000, GdPtLongDataType, {
	      distancesSqr[ii * numRemain2 + jj] = PointDistNDSqr(g1Core[i * grain + ii], g2Core[j * grain + jj], dim);
	    });
	}
	PtIntermediateDataType minDistSqr = sequence::reduce(distancesSqr, (OtLongDataType) grain * numRemain2, utils::minF<PtIntermediateDataType>());
	if (minDistSqr <= t_radiusSqr1) {
#if defined(USE_JEMALLOC)
	  JE_MALLOC_FREE(distancesSqr);
	  JE_MALLOC_FREE(g1Core);
#else
	  free(distancesSqr);
	  free(g1Core);
#endif
	  return true;
	}
      }
    } // iters

    if (numRemain > 0) {
      GdPtLongDataType j = 0;
      for (; j < numIters2; ++ j) {
	parallel_for_1(GdPtLongDataType ii = 0; ii < numRemain; ++ ii) {
	  granular_for_t(jj, 0, grain, 3000, GdPtLongDataType, {
	      distancesSqr[ii * grain + jj] = PointDistNDSqr(g1Core[i * grain + ii], g2Core[j * grain + jj], dim);
	    });
	}
	PtIntermediateDataType minDistSqr = sequence::reduce(distancesSqr, (OtLongDataType) numRemain * grain, utils::minF<PtIntermediateDataType>());
	if (minDistSqr <= t_radiusSqr1) {
#if defined(USE_JEMALLOC)
	  JE_MALLOC_FREE(distancesSqr);
	  JE_MALLOC_FREE(g1Core);
#else
	  free(distancesSqr);
	  free(g1Core);
#endif
	  return true;
	}
      }

      if (numRemain2 > 0) {
	parallel_for_1(GdPtLongDataType ii = 0; ii < numRemain; ++ ii) {
	  granular_for_t(jj, 0, numRemain2, 3000, GdPtLongDataType, {
	      distancesSqr[ii * numRemain2 + jj] = PointDistNDSqr(g1Core[i * grain + ii], g2Core[j * grain + jj], dim);
	    });
	}
	PtIntermediateDataType minDistSqr = sequence::reduce(distancesSqr, (OtLongDataType) numRemain * numRemain2, utils::minF<PtIntermediateDataType>());
	if (minDistSqr <= t_radiusSqr1) {
#if defined(USE_JEMALLOC)
	  JE_MALLOC_FREE(distancesSqr);
	  JE_MALLOC_FREE(g1Core);
#else
	  free(distancesSqr);
	  free(g1Core);
#endif
	  return true;
	}
      } // rem2
    } // rem1


#if defined(USE_JEMALLOC)
    JE_MALLOC_FREE(distancesSqr);
#else
    free(distancesSqr);
#endif

  }

#if defined(USE_JEMALLOC)
  JE_MALLOC_FREE(g1Core);
#else
  free(g1Core);
#endif
  return false;

}

// tries to automatically choose the right version of bcp to call
// this is tailored for dealing with extremely larger grids

// it could be to call this in serial, for example in geom cases
/* other functions to use, but have not found use case
  return edgeAddableBCPNDWithSort(pointsData, dim, t_radiusSqr1, t_grids, t_corePointMarkers, g1, g2);
  return edgeAddableBCPNDWithG1Sort(pointsData, dim, t_radiusSqr1, t_grids, t_corePointMarkers, g1, g2);
*/
inline bool edgeAddableBCPNDSmart(PtDataType *pointsData, intT dim, PtIntermediateDataType t_radiusSqr1, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2) {
  const OtLongDataType thresh1 = 100;

  GdPtLongDataType g1NumPts = t_grids->GetNumPts(t_g1);
  GdPtLongDataType g2NumPts = t_grids->GetNumPts(t_g2);
  if (g1NumPts < thresh1 && g2NumPts < thresh1) {
    // this one helps on most datasets
    return edgeAddableBCPNDSimple(pointsData, dim, t_radiusSqr1, t_grids, t_corePointMarkers, t_g1, t_g2);
  }
  GdLongDataType g1 = g1NumPts > g2NumPts ? t_g1 : t_g2;
  GdLongDataType g2 = g1NumPts > g2NumPts ? t_g2 : t_g1;
  
  /* helps in geolife dataset, very skewed case, need to make edgeaddable calls sequential
  GdPtLongDataType g1NumPtsNew = t_grids->GetNumPts(g1);
  GdPtLongDataType g2NumPtsNew = t_grids->GetNumPts(g2);
  if (g1NumPtsNew > thresh2 || g2NumPtsNew > thresh2 ) {
  */
  return edgeAddableBCPNDParallel(pointsData, dim, t_radiusSqr1, t_grids, t_corePointMarkers, g2, g1);
    //  }

  //return edgeAddableBCPND(pointsData, dim, t_radiusSqr1, t_grids, t_corePointMarkers, g2, g1);
}

/**
 * Returns cluster id for all points, given boxes and core points.
 * @param t_points Sequence of input points.
 * @param t_params DBSCAN parameters.
 * @param t_grids  Grids data structure.
 * @param t_isCoreGrid A boolean vector indicating whether a grid has >=1 core points (core grid).
 * @param t_corePointMarkers pair<numCorePoints, corePointFlags>
 * @return A vector of flag of length t_points->n, containing clusterIds for core points, and -1 for non-assigned (border) points.
 */
GdLongDataType* ClusterCoreBCPGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid) {

  bool par_edgeadd = (t_points->n / t_grids->GetNumGrids() > 10000) || t_grids->GetNumGrids() < 1000;

  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  PtIntermediateDataType threshold = t_params->m_epsilon * t_params->m_epsilon;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());

  // Iterate through each grid (that has points / valid).
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {

    // If grid not core, no need to look.
    if (isCoreGrid[g]) {

      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
        GdLongDataType n = neighbors->at(ni);

        if (isCoreGrid[n] && n >= 0) {
          if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
	    if(edgeAddableBCPNDSmart(pointsData, dim, threshold, t_grids, corePointMarkersGrid.second, g, n)) {
	      bcpUF.link(n, g);
	    }
          }
        }
        
      }
    }
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }

  // Now that every core-grid should have a cluster id (root). Should assign corresponding core point to those ids.
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
        if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
          clusterIds[t_grids->GetPt(p)] = gridCluster;
        }
      });
    }
  }
  //PrintParser("timer3", t1.stop());

  return clusterIds;
}

GdLongDataType* ClusterCoreBCPGridParallelBucketingND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid) {
  cout << "clustercore bucketing" << endl;
   
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  PtIntermediateDataType threshold = t_params->m_epsilon * t_params->m_epsilon;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());

  GdLongDataType *tmp_sorted = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids() );
  GdLongDataType *ordering = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids() );
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    tmp_sorted[g]= g;
  }
  sampleSort(tmp_sorted, t_grids->GetNumGrids(), gridDegComparator(t_grids));
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    ordering[tmp_sorted[g]]= g;
  }
  free(tmp_sorted);

  GdLongDataType remGrids = t_grids->GetNumGrids();
  GdLongDataType blkSize = getWorkers() * 8;
  GdLongDataType gdOffset = 0;
  long avgPtsPerGrid = t_points->n/remGrids;
  bool serialBuckets = (remGrids <= blkSize) && (avgPtsPerGrid >= 2*blkSize);
  if(!serialBuckets) {
    while(remGrids > 0) {
      GdLongDataType loopEnd = min(remGrids, blkSize) + gdOffset;
      parallel_for_1 (GdLongDataType o = gdOffset; o < loopEnd; ++ o) {
	GdLongDataType g = ordering[o];
	if (isCoreGrid[g]) {

	  NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
	  parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
	    GdLongDataType n = neighbors->at(ni);
	    if (isCoreGrid[n] && n >= 0) {

	      if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
		if(edgeAddableBCPNDSmart(pointsData, dim, threshold, t_grids, corePointMarkersGrid.second, g, n)) {
		  bcpUF.link(n, g);
		}
	      }
	    
	    }
	  }
	}
      }
      gdOffset = loopEnd;
      remGrids = t_grids->GetNumGrids() - loopEnd;
      blkSize *= 2;
    }
  }
  else {
    for (GdLongDataType o = 0; o < remGrids; ++ o) {
      GdLongDataType g = ordering[o];
      if (isCoreGrid[g]) {

	NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
	for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
	  GdLongDataType n = neighbors->at(ni);
	  if (isCoreGrid[n] && n >= 0) {
	      
	    if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
	      if(edgeAddableBCPNDSmart(pointsData, dim, threshold, t_grids, corePointMarkersGrid.second, g, n)) {
		bcpUF.link(n, g);
	      }
	    }
	  }
	}
      }
    }
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }

  // Now that every core-grid should have a cluster id (root). Should assign corresponding core point to those ids.
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
        if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
          clusterIds[t_grids->GetPt(p)] = gridCluster;
        }
      });
    }
  }

  return clusterIds;
}

#ifdef QUADTREE_H

GdLongDataType* ClusterCoreBCPGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid, QuadTree *t_corePtQuadTree) {

  const OtLongDataType PAR_CONSTANT_BCP = 4294967295;
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());

  // Iterate through each grid (that has points / valid).
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    
    // If grid not core, no need to look.
    if (isCoreGrid[g]) {

      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
        GdLongDataType n = neighbors->at(ni);

        if (isCoreGrid[n] && n >= 0) {
          if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
            GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
            GdPtLongDataType nNumPts = t_grids->GetNumPts(n);
            if(edgeAddableBCPND(t_corePtQuadTree, pointsData, dim, t_params, t_grids, corePointMarkersGrid.second, g, n)) {
              bcpUF.link(n, g);
            }
          }
        }
      }
    }
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }

  // Now that every core-grid should have a cluster id (root). Should assign corresponding core point to those ids.
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
        if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
          clusterIds[t_grids->GetPt(p)] = gridCluster;
        }
      });
    }
  }

  return clusterIds;
}

// bucketing method for quadtree
GdLongDataType* ClusterCoreBCPGridParallelBucketingND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid, QuadTree *t_corePtQuadTree) {
  cout << "quadtree clustercore bucketing" << endl;
  
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());
  
  GdLongDataType *tmp_sorted = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids());
  GdLongDataType *ordering = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids());
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    tmp_sorted[g]= g;
  }
  sampleSort(tmp_sorted, t_grids->GetNumGrids(), gridDegComparator(t_grids));
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    ordering[tmp_sorted[g]]= g;
  }
  free(tmp_sorted);

  GdLongDataType remGrids = t_grids->GetNumGrids();
  GdLongDataType blkSize = getWorkers() * 8;
  GdLongDataType gdOffset = 0;
  long avgPtsPerGrid = t_points->n/remGrids;
  bool serialBuckets = (remGrids <= blkSize) && (avgPtsPerGrid >= 2*blkSize);
  if(!serialBuckets) { // parallel bucketing
    while(remGrids > 0) {
      GdLongDataType loopEnd = min(remGrids, blkSize) + gdOffset;
      parallel_for_1 (GdLongDataType o = gdOffset; o < loopEnd; ++ o) {
	GdLongDataType g = ordering[o];
	if (isCoreGrid[g]) {

	  NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
	  parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
	    GdLongDataType n = neighbors->at(ni);
	    if (isCoreGrid[n] && n >= 0) {
	      if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
		GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
		GdPtLongDataType nNumPts = t_grids->GetNumPts(n);
		if(edgeAddableBCPND(t_corePtQuadTree, pointsData, dim, t_params, t_grids, corePointMarkersGrid.second, g, n)) {
		  bcpUF.link(n, g);
		}
	      }
	    }
	  }
	}
      }
      gdOffset = loopEnd;
      remGrids = t_grids->GetNumGrids() - loopEnd;
      blkSize *= 2;
    }
  }
  else { // serial bucketing
    for (GdLongDataType o = 0; o < remGrids; ++ o) {
      GdLongDataType g = ordering[o];
      if (isCoreGrid[g]) {
	NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
	for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
	    GdLongDataType n = neighbors->at(ni);
	    if (isCoreGrid[n] && n >= 0) {
	      if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
		GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
		GdPtLongDataType nNumPts = t_grids->GetNumPts(n);
		if(edgeAddableBCPND(t_corePtQuadTree, pointsData, dim, t_params, t_grids, corePointMarkersGrid.second, g, n)) {
		  bcpUF.link(n, g);
		}
	      }
	    }
	  }
      }
    }
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }

  // mark core point labels
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
	  if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
	    clusterIds[t_grids->GetPt(p)] = gridCluster;
	  }
	});
    }
  }

  return clusterIds;
}

bool EdgeAddableApprox(PtDataType *pointsData, intT t_dim, PtIntermediateDataType t_radiusSqr1, PtIntermediateDataType t_radiusSqr2, Grids *t_grids, intT *t_corePointMarkers, GdLongDataType t_g1, GdLongDataType t_g2, QuadTree *t_Q) {
  // select the grid with fewer points to iterate on
  GdLongDataType queryingGrid;
  GdLongDataType queriedGrid;
  if (t_grids->GetNumPts(t_g1) > t_grids->GetNumPts(t_g2)) {
    queriedGrid = t_g1;
    queryingGrid = t_g2;
  } else {
    queriedGrid = t_g2;
    queryingGrid = t_g1;
  }

  PtLongDataType iterEnd = t_grids->GetStartOffset(queryingGrid) + t_grids->GetNumPts(queryingGrid);
  for (PtLongDataType p = t_grids->GetStartOffset(queryingGrid); p < iterEnd; ++p) {
    if (t_corePointMarkers[t_grids->GetPt(p)] == 1) {

      if (PointGridDistanceLeq(t_dim, t_grids, queriedGrid, p, t_radiusSqr1)) {
      	if (t_Q->BinApproxRangeQuery(pointsData, t_grids->GetPtData(p), t_dim, t_grids, queriedGrid, t_radiusSqr1, t_radiusSqr2)) {
      	  return true;
      	}
      }
    }
  }
  return false;
}

// approximate cluster core
GdLongDataType* ClusterCoreBCPGridApproxParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid, QuadTree *quadTree) {

  // each core point will query on the other box
  // once found a pair of core points that is close enough, return true
  PtIntermediateDataType radiusSqr1 = t_params->m_epsilon * t_params->m_epsilon;
  PtIntermediateDataType radiusSqr2 = t_params->m_epsilon * (1 + t_params->m_rho) * t_params->m_epsilon * (1 + t_params->m_rho);
  
  if (t_params->m_rho == -1) {
    cout << "error: rho not set for approx-dbscan, exiting" << endl;
    exit(1);
  }
  
  intT dim = t_points->A[0].m_dim;
  PtDataType *t_pointsData = t_points->A[0].m_data;

  const OtLongDataType PAR_CONSTANT_BCP = 4294967295;
  PtDataType *pointsData = t_points->A[0].m_data;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());

  parallel_for_1 (GdLongDataType o = 0; o < t_grids->GetNumGrids(); ++ o) {

    GdLongDataType g = o;

    if (isCoreGrid[g]) {
      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
        GdLongDataType n = neighbors->at(ni);

        if (isCoreGrid[n] && n >= 0) {
          if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
	    bool E = EdgeAddableApprox(pointsData, dim, radiusSqr1, radiusSqr2, t_grids, corePointMarkersGrid.second, g, n, quadTree);
            if(E) {
              bcpUF.link(g,n);
            }
#if defined(APPROX_VERI)
            bool E1 = edgeAddableBCPNDSmart(pointsData, dim, radiusSqr1, t_grids, corePointMarkersGrid.second, g, n);
            bool E2 = edgeAddableBCPNDSmart(pointsData, dim, radiusSqr2, t_grids, corePointMarkersGrid.second, g, n);
            if (!( E1 <= E && E <= E2 )) {
              cout << "error: E1 " << E1 << " E " << E << " E2 " << E2 << endl;   
              exit(1);
            }
#endif
          }
        }
      }
    }
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }
  
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
        if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
          clusterIds[t_grids->GetPt(p)] = gridCluster;
        }
      });
    }
  }
  
  return clusterIds;
}


GdLongDataType* ClusterCoreBCPGridApproxParallelBucketingND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, bool *isCoreGrid, pair<PtLongDataType, intT *> corePointMarkersGrid, QuadTree *quadTree) {
  cout << "clustercore approx bucketing" << endl;

  PtIntermediateDataType radiusSqr1 = t_params->m_epsilon * t_params->m_epsilon;
  PtIntermediateDataType radiusSqr2 = t_params->m_epsilon * (1 + t_params->m_rho) * t_params->m_epsilon * (1 + t_params->m_rho);
  
  if (t_params->m_rho == -1) {
    cout << "error: rho not set for approx-dbscan, exiting" << endl;
    exit(1);
  }
   
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  PtIntermediateDataType threshold = t_params->m_epsilon * t_params->m_epsilon;

  parallelUnionFind bcpUF = parallelUnionFind(t_grids->GetNumGrids());

  GdLongDataType *tmp_sorted = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids() );
  GdLongDataType *ordering = (GdLongDataType *) malloc(sizeof(GdLongDataType) * t_grids->GetNumGrids() );
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    tmp_sorted[g]= g;
  }
  sampleSort(tmp_sorted, t_grids->GetNumGrids(), gridDegComparator(t_grids));
  parallel_for (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    ordering[tmp_sorted[g]]= g;
  }
  free(tmp_sorted);

  GdLongDataType remGrids = t_grids->GetNumGrids();
  GdLongDataType blkSize = getWorkers() * 8;
  GdLongDataType gdOffset = 0;
  while(remGrids > 0) {
    GdLongDataType loopEnd = min(remGrids, blkSize) + gdOffset;
    parallel_for_1 (GdLongDataType o = gdOffset; o < loopEnd; ++ o) {
      GdLongDataType g = ordering[o];
      if (isCoreGrid[g]) {

	NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
        parallel_for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
          GdLongDataType n = neighbors->at(ni);
          if (isCoreGrid[n] && n >= 0) {
	    if (g > n && bcpUF.find(g) != bcpUF.find(n)) {
	      bool E = EdgeAddableApprox(pointsData, dim, radiusSqr1, radiusSqr2, t_grids, corePointMarkersGrid.second, g, n, quadTree);
	      if(E) {
		bcpUF.link(g,n);
	      }
#if defined(APPROX_VERI)
	      bool E1 = edgeAddableBCPNDSmart(pointsData, dim, radiusSqr1, t_grids, corePointMarkersGrid.second, g, n);
	      bool E2 = edgeAddableBCPNDSmart(pointsData, dim, radiusSqr2, t_grids, corePointMarkersGrid.second, g, n);
	      if (!( E1 <= E && E <= E2 )) {
		cout << "error: E1 " << E1 << " E " << E << " E2 " << E2 << endl;   
		exit(1);
	      }
#endif
	    }
          }
        }
      }
    }
    gdOffset = loopEnd;
    remGrids = t_grids->GetNumGrids() - loopEnd;
    blkSize *= 2;
  }

  GdLongDataType *clusterIds = newA(GdLongDataType,t_points->n);
  parallel_for (PtLongDataType i = 0; i < t_points->n; ++ i) {
    clusterIds[i] = -1;
  }

  // Now that every core-grid should have a cluster id (root). Should assign corresponding core point to those ids.
  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {
    if (isCoreGrid[g]) {
      GdLongDataType gridCluster = bcpUF.find(g);
      GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
      PtLongDataType iterEnd = t_grids->GetStartOffset(g) + gNumPts;
      granular_for_t(p, t_grids->GetStartOffset(g), iterEnd, 10000, PtLongDataType, {
        if (corePointMarkersGrid.second[t_grids->GetPt(p)] != 0) { // Is a core point.
          clusterIds[t_grids->GetPt(p)] = gridCluster;
        }
      });
    }
  }

  return clusterIds;
}


#endif // ifdef quadtree h

/**
 * Assigns border points after assigning core points to clusters.
 * @param t_points Input sequence, must be after calling assignBox (reordered).
 * @param t_params DBSCAN parameters.
 * @param t_grids  Grid data structure constructed earlier.
 * @param t_clusterIds A vector of cluster ids of all points. -1 means unmarked.
 * @param t_isCoreGrid A boolean vector indicating whether a grid has >=1 core points (core grid).
 * @return Number of border points newly assigned to clusters.
 */
void AssignBorderPointGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, GdLongDataType *t_clusterIds, intT *t_corePointMarkers, bool *isCoreGrid) {
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  PtIntermediateDataType epsSqr = static_cast<PtIntermediateDataType>(t_params->m_epsilon) * static_cast<PtIntermediateDataType>(t_params->m_epsilon);

  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {

    GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
    if (gNumPts < t_params->m_minPts) {

      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      PtLongDataType iterEndOuter = t_grids->GetStartOffset(g) + gNumPts;
      parallel_for (PtLongDataType p = t_grids->GetStartOffset(g); p < iterEndOuter; ++p) {

        if (t_corePointMarkers[t_grids->GetPt(p)] == 0) { // Is a noncore point.

          PtLongDataType nearCore = -1;
          // Let the initial near distance be not as good as qualifying a reachable point.
          PtIntermediateDataType nearDist = epsSqr + 1; // Hack.

          // First look at grid g itself.
          if (isCoreGrid[g]) { // If has no core points, no need to look
	          PtLongDataType iterEndInner = t_grids->GetStartOffset(g) + gNumPts;
            for (PtLongDataType np = t_grids->GetStartOffset(g); np < iterEndInner; ++np) {
              if (t_corePointMarkers[t_grids->GetPt(np)] != 0) { // Core point.
                PtIntermediateDataType myDistSqr = PointDistNDSqr(t_grids->GetPtData(p), t_grids->GetPtData(np), dim);
                if (myDistSqr < nearDist && myDistSqr <= epsSqr ) {
                  nearCore = t_grids->GetPt(np);
                  nearDist = myDistSqr;
                }
              }
            }
          }

          // Then iterate neighboring boxes, find the nearets core to p.
	        for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
            GdLongDataType n = neighbors->at(ni);
	          if (n < 0 || n == g) continue;

            if (isCoreGrid[n]) {

              GdPtLongDataType nNumPts = t_grids->GetNumPts(n);

              // Iterate through core points of neighboring grids to determine closest core point.
	            PtLongDataType iterEndInner = t_grids->GetStartOffset(n) + nNumPts;
              for (PtLongDataType np = t_grids->GetStartOffset(n); np < iterEndInner; ++np) {
                if (t_corePointMarkers[t_grids->GetPt(np)] != 0) { // Core point.
                  PtIntermediateDataType myDistSqr = PointDistNDSqr(t_grids->GetPtData(p), t_grids->GetPtData(np), dim);
                  if (myDistSqr < nearDist && myDistSqr <= epsSqr ) {
                    nearCore = t_grids->GetPt(np);
                    nearDist = myDistSqr;
                  }
                }
              }
            }
          }

          if (nearCore != -1) {
            t_clusterIds[t_grids->GetPt(p)] = t_clusterIds[nearCore];
          }
        }
      }
    }
  }
}

//void AssignBorderPointMultiGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, GdLongDataType *t_clusterIds, intT *t_corePointMarkers, bool *isCoreGrid, Table<hashBorderCluster, OtLongDataType> *t_borderHash) {
void AssignBorderPointMultiGridParallelND(_seq<pointNd> *t_points, DbscanParams *t_params, Grids *t_grids, GdLongDataType *t_clusterIds, intT *t_corePointMarkers, bool *isCoreGrid, ClusterContainer **t_borderHash) {
  intT dim = t_points->A[0].m_dim;
  PtDataType *pointsData = t_points->A[0].m_data;
  PtIntermediateDataType epsSqr = static_cast<PtIntermediateDataType>(t_params->m_epsilon) * static_cast<PtIntermediateDataType>(t_params->m_epsilon);

  parallel_for_1 (GdLongDataType g = 0; g < t_grids->GetNumGrids(); ++ g) {

    GdPtLongDataType gNumPts = t_grids->GetNumPts(g);
    if (gNumPts < t_params->m_minPts) {

      NbrContainer *neighbors = GetNeighborGridsKDTree(g, t_grids);
      PtLongDataType iterEndOuter = t_grids->GetStartOffset(g) + gNumPts;
      parallel_for (PtLongDataType p = t_grids->GetStartOffset(g); p < iterEndOuter; ++p) {

        if (t_corePointMarkers[t_grids->GetPt(p)] == 0) { // Is a noncore point.
	  bool isBorder = false;
	  ClusterContainer *myClusterContainer; // do not initialize yet

          // First look at grid g itself.
          if (isCoreGrid[g]) { // If has no core points, no need to look
	    PtLongDataType iterEndInner = t_grids->GetStartOffset(g) + gNumPts;
            for (PtLongDataType np = t_grids->GetStartOffset(g); np < iterEndInner; ++np) {
              if (t_corePointMarkers[t_grids->GetPt(np)] != 0) { // Core point.
                PtIntermediateDataType myDistSqr = PointDistNDSqr(t_grids->GetPtData(p), t_grids->GetPtData(np), dim);
                if (myDistSqr <= epsSqr ) { // in circle
		  if (!isBorder) {
		    isBorder = true;
		    myClusterContainer = new ClusterContainer();
		  }
		  myClusterContainer->insert(t_clusterIds[t_grids->GetPt(np)]);
                }
              }
            }
          }

          // Then iterate neighboring boxes, find the nearets core to p.
	  for (GdLongDataType ni = 0; ni < neighbors->size(); ++ ni) {
            GdLongDataType n = neighbors->at(ni);
	    if (n < 0 || n == g) continue;

            if (isCoreGrid[n]) {

              GdPtLongDataType nNumPts = t_grids->GetNumPts(n);

              // Iterate through core points of neighboring grids to determine closest core point.
	      PtLongDataType iterEndInner = t_grids->GetStartOffset(n) + nNumPts;
              for (PtLongDataType np = t_grids->GetStartOffset(n); np < iterEndInner; ++np) {
                if (t_corePointMarkers[t_grids->GetPt(np)] != 0) { // Core point.
                  PtIntermediateDataType myDistSqr = PointDistNDSqr(t_grids->GetPtData(p), t_grids->GetPtData(np), dim);
                  if (myDistSqr <= epsSqr ) {
		    if (!isBorder) {
		      isBorder = true;
		      myClusterContainer = new ClusterContainer();
		    }
		    myClusterContainer->insert(t_clusterIds[t_grids->GetPt(np)]);
                  }
                }
              }
            }
          }

          if (isBorder) {
	    t_borderHash[t_grids->GetPt(p)] = myClusterContainer;
            t_clusterIds[t_grids->GetPt(p)] = -2; // -1 means noise, -2 means border
          }
        }
      }
    }
  }
}

/**
 * Computes DBSCAN.
 * @param t_P Sequence of input points, typed point2d (_point2d<PtDataType>). The point order WILL be altered during the computation, and the returned value corresponds to the ALTERED points.
 * @param t_param DBSCAN parameter, defined in DBSCAN.h.
 * @return An array of positive integers indicating cluster ID for each point. -1 means no cluster assigned.
 */
#if defined CLOSEST_BORDER_CLUSTER
pair<intT*, GdLongDataType*> DBSCAN(_seq<pointNd> t_P, DbscanParams t_params)
#else
pair<GdLongDataType*, ClusterContainer** > DBSCAN(_seq<pointNd> t_P, DbscanParams t_params)
#endif
{

  PrintParser("serialalgo", "parallel");
  PrintParser("ndalgo", "nd");
#if defined(APPROX)
  PrintParser("algo", "approx");
#else
  PrintParser("algo", "exact");
#endif
  PrintParser("numprocs", getWorkers());

  int dim = t_P.A[0].m_dim;

#ifdef COMPONENT_TIMER
  startTime();
#endif
  Grids *grids = AssignGridParallelND(&t_P, &t_params, dim);
#ifdef COMPONENT_TIMER
  double assignGridTime = _tm.next();
  PrintParser("assigngridtime", assignGridTime);
#endif

  PrintParser("numgrids", grids->GetNumGrids());

#if defined(ALLPT_QUADTREE)
  QuadTree *allPtQuadTree = new QuadTree(&t_P, &t_params, grids);
#ifdef COMPONENT_TIMER
  double buildAllPtQuadTreeTime = _tm.next();
  cout << " build all pt quadtree: " << buildAllPtQuadTreeTime << endl;
#endif
#endif

  GetNeighborGridKDTreeInit(grids);
#ifdef COMPONENT_TIMER
  double buildGetNbrKdtTime = _tm.next();
  PrintParser("buildgetnbrkdt", buildGetNbrKdtTime);
#endif

  bool *isCoreGrid = newA(bool,grids->GetNumGrids());

#if defined(QUADTREE_MARKCORE)
  pair<PtLongDataType, intT *> corePointMarkersGrid = MarkCorePointGridParallelND(&t_P, &t_params, grids, allPtQuadTree, isCoreGrid);
#else
  pair<PtLongDataType, intT *> corePointMarkersGrid = MarkCorePointGridParallelND(&t_P, &t_params, grids, isCoreGrid);
#endif
  
#ifdef COMPONENT_TIMER
  double markcoreTime = _tm.next();
#if defined(ALLPT_QUADTREE)
  markcoreTime += buildAllPtQuadTreeTime;
#endif
  PrintParser("markcoretime", markcoreTime);
#endif
  
  PrintParser("numcore", corePointMarkersGrid.first);
  PrintParser("numpts", t_P.n);

#if defined(COREPT_QUADTREE) || defined(APPROX)
  QuadTree *corePtQuadTree = new QuadTree(&t_P, &t_params, grids, corePointMarkersGrid, isCoreGrid);

#ifdef COMPONENT_TIMER
  double coreptquadtreetime = _tm.next();
  PrintParser("buildcoreptquadtree", coreptquadtreetime);
#endif
#endif

#ifdef APPROX
  
#if defined (USE_BUCKETING)
  GdLongDataType *clusterIds = ClusterCoreBCPGridApproxParallelBucketingND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid, corePtQuadTree);
#else
  GdLongDataType *clusterIds = ClusterCoreBCPGridApproxParallelND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid, corePtQuadTree);
#endif
  
#else
  
#if defined(QUADTREE_BCP)
  
#if defined (USE_BUCKETING)
  GdLongDataType *clusterIds = ClusterCoreBCPGridParallelBucketingND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid, corePtQuadTree);
#else
  GdLongDataType *clusterIds = ClusterCoreBCPGridParallelND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid, corePtQuadTree);
#endif
  
#else

#if defined (USE_BUCKETING)
  GdLongDataType *clusterIds = ClusterCoreBCPGridParallelBucketingND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid);
#else
  GdLongDataType *clusterIds = ClusterCoreBCPGridParallelND(&t_P, &t_params, grids, isCoreGrid, corePointMarkersGrid);
#endif // bucketing
  
#endif
#endif

#ifdef COMPONENT_TIMER
  // Timer code in function
  double clusterCoreTime = _tm.next();
  PrintParser("clustercoretime", clusterCoreTime);
#endif

#if defined CLOSEST_BORDER_CLUSTER
  AssignBorderPointGridParallelND(&t_P, &t_params, grids, clusterIds, corePointMarkersGrid.second, isCoreGrid);
#else
  ClusterContainer **myBorderClusterTable = newA(ClusterContainer *, t_P.n);
  AssignBorderPointMultiGridParallelND(&t_P, &t_params, grids, clusterIds, corePointMarkersGrid.second, isCoreGrid, myBorderClusterTable);
#endif
  
#ifdef COMPONENT_TIMER
  double borderTime = _tm.next();
  PrintParser("assignbordertime", borderTime);
#endif

#if defined ALLPT_QUADTREE
  delete allPtQuadTree;
#endif
#if defined COREPT_QUADTREE
  delete corePtQuadTree;
#endif
  
  delete grids;
  delete g_kdtree;
  free(isCoreGrid);

#ifdef COMPONENT_TIMER
  double cleanupTime = _tm.next();
  PrintParser("cleanuptime", cleanupTime);
#endif

#if defined CLOSEST_BORDER_CLUSTER
  return make_pair(corePointMarkersGrid.second, clusterIds);
#else
  free(corePointMarkersGrid.second);
  return make_pair(clusterIds, myBorderClusterTable);
#endif
}
