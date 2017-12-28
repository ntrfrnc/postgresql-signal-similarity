#ifndef SIGNAL_SIMILARITY_H
#define SIGNAL_SIMILARITY_H

#include "postgres.h"
#include <math.h>
#include <string.h>
#include "catalog/pg_type.h"
#include "utils/builtins.h"
#include "utils/array.h"
#include "utils/varlena.h"
#include "utils/formatting.h"
#include "catalog/pg_collation.h"
#include "fmgr.h"

PG_MODULE_MAGIC;

#define ARRDOUBLEPTR(x) ((double*)ARR_DATA_PTR(x))
#define ARRNELEMS(x) ArrayGetNItems(ARR_NDIM(x), ARR_DIMS(x))
#define AARR_FREE_IF_COPY(array,n) \
    do { \
        if (!VARATT_IS_EXPANDED_HEADER(array)) \
            PG_FREE_IF_COPY(array, n); \
    } while (0)

#define min(x, y) ((x) < (y) ? (x) : (y))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define dist(x, y) ((x - y) * (x - y))

#define INF 1e20  // Pseudo Infitinte number for this code

PG_FUNCTION_INFO_V1(signal_similarity);

/// Data structure for similarity search results
typedef struct ssr {
  double distance;
  long long location;
} ssr;

/// Data structure for sorting the query
typedef struct Idx {
  double value;
  int index;
} Idx;

/// Data structure (circular array) for finding minimum and maximum for LB_Keogh
/// envolop
typedef struct deque {
  int* dq;
  int size, capacity;
  int f, r;
} deque;

/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void* a, const void* b);

/// Initial the queue at the begining step of envelop calculation
void init(deque* d, int capacity);

/// Destroy the queue
void destroy(deque* d);

/// Insert to the queue at the back
void push_back(struct deque* d, int v);

/// Delete the current (front) element from queue
void pop_front(struct deque* d);

/// Delete the last element from queue
void pop_back(struct deque* d);

/// Get the value at the current position of the circular queue
int front(struct deque* d);

/// Get the value at the last position of the circular queueint back(struct
/// deque *d)
int back(struct deque* d);

/// Check whether or not the queue is empty
int empty(struct deque* d);

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern
/// Recognition 42(9), 2009.
void lower_upper_lemire(double* t, int len, int r, double* l, double* u);

/// Calculate quick lower bound
/// Usually, LB_Kim take time O(m) for finding top,bottom,fist and last.
/// However, because of z-normalization the top and bottom cannot give
/// siginifant benefits.
/// And using the first and last points can be computed in constant time.
/// The prunning power of LB_Kim is non-trivial, especially when the query is
/// not long, say in length 128.
double lb_kim_hierarchy(double* t,
                        double* q,
                        int j,
                        int len,
                        double mean,
                        double std,
                        double bsf);

/// LB_Keogh 1: Create Envelop for the query
/// Note that because the query is known, envelop can be created once at the
/// begenining.
///
/// Variable Explanation,
/// order : sorted indices for the query.
/// uo, lo: upper and lower envelops for the query, which already sorted.
/// t     : a circular array keeping the current data.
/// j     : index of the starting location in t
/// cb    : (output) current bound at each position. It will be used later for
/// early abandoning in DTW.
double lb_keogh_cumulative(int* order,
                           double* t,
                           double* uo,
                           double* lo,
                           double* cb,
                           int j,
                           int len,
                           double mean,
                           double std,
                           double best_so_far);

/// LB_Keogh 2: Create Envelop for the data
/// Note that the envelops have been created (in main function) when each data
/// point has been read.
///
/// Variable Explanation,
/// tz: Z-normalized data
/// qo: sorted query
/// cb: (output) current bound at each position. Used later for early abandoning
/// in DTW.
/// l,u: lower and upper envelop of the current data
double lb_keogh_data_cumulative(int* order,
                                double* tz,
                                double* qo,
                                double* cb,
                                double* l,
                                double* u,
                                int len,
                                double mean,
                                double std,
                                double best_so_far);

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double* cb, int m, int r, double bsf);

/// Log errors
void error(int id);

/// Search for query in data, returns nearest neighbour distance and location
/// data,query: data and query, respectively
/// dl, ql: data and query lengths
/// ww : warping window (0-1)
/// return  : double* {distance, location}
ssr search_query(double* data, long long di, long long dl, double* query, long long qi, long long ql, float4 ww);

char** text_array_to_str_array(ArrayType *array);

Datum signal_similarity(PG_FUNCTION_ARGS);

#endif // SIGNAL_SIMILARITY_H
