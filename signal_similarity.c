#include "signal_similarity.h"

/// Sorting function for the query, sort by abs(z_norm(q[i])) from high to low
int comp(const void* a, const void* b) {
  Idx* x = (Idx*)a;
  Idx* y = (Idx*)b;
  return abs(y->value) - abs(x->value);  // high to low
}

/// Initial the queue at the begining step of envelop calculation
void init(deque* d, int capacity) {
  d->capacity = capacity;
  d->size = 0;
  d->dq = (int*)malloc(sizeof(int) * d->capacity);
  d->f = 0;
  d->r = d->capacity - 1;
}

/// Destroy the queue
void destroy(deque* d) {
  free(d->dq);
}

/// Insert to the queue at the back
void push_back(struct deque* d, int v) {
  d->dq[d->r] = v;
  d->r--;
  if (d->r < 0)
    d->r = d->capacity - 1;
  d->size++;
}

/// Delete the current (front) element from queue
void pop_front(struct deque* d) {
  d->f--;
  if (d->f < 0)
    d->f = d->capacity - 1;
  d->size--;
}

/// Delete the last element from queue
void pop_back(struct deque* d) {
  d->r = (d->r + 1) % d->capacity;
  d->size--;
}

/// Get the value at the current position of the circular queue
int front(struct deque* d) {
  int aux = d->f - 1;

  if (aux < 0)
    aux = d->capacity - 1;
  return d->dq[aux];
}

/// Get the value at the last position of the circular queueint back(struct
/// deque *d)
int back(struct deque* d) {
  int aux = (d->r + 1) % d->capacity;
  return d->dq[aux];
}

/// Check whether or not the queue is empty
int empty(struct deque* d) {
  return d->size == 0;
}

/// Finding the envelop of min and max value for LB_Keogh
/// Implementation idea is intoruduced by Danial Lemire in his paper
/// "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern
/// Recognition 42(9), 2009.
void lower_upper_lemire(double* t, int len, int r, double* l, double* u) {
  struct deque du, dl;

  init(&du, 2 * r + 2);
  init(&dl, 2 * r + 2);

  push_back(&du, 0);
  push_back(&dl, 0);

  for (int i = 1; i < len; i++) {
    if (i > r) {
      u[i - r - 1] = t[front(&du)];
      l[i - r - 1] = t[front(&dl)];
    }
    if (t[i] > t[i - 1]) {
      pop_back(&du);
      while (!empty(&du) && t[i] > t[back(&du)]) {
        pop_back(&du);
      }
    } else {
      pop_back(&dl);
      while (!empty(&dl) && t[i] < t[back(&dl)]) {
        pop_back(&dl);
      }
    }

    push_back(&du, i);
    push_back(&dl, i);

    if (i == 2 * r + 1 + front(&du)) {
      pop_front(&du);
    } else if (i == 2 * r + 1 + front(&dl)) {
      pop_front(&dl);
    }
  }

  for (int i = len; i < len + r + 1; i++) {
    u[i - r - 1] = t[front(&du)];
    l[i - r - 1] = t[front(&dl)];

    if (i - front(&du) >= 2 * r + 1) {
      pop_front(&du);
    }
    if (i - front(&dl) >= 2 * r + 1) {
      pop_front(&dl);
    }
  }

  destroy(&du);
  destroy(&dl);
}

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
                        double bsf) {
  /// 1 point at front and back
  double d, lb, x1, y1, x2, y2;
  double x0 = (t[j] - mean) / std;
  double y0 = (t[(len - 1 + j)] - mean) / std;
  lb = dist(x0, q[0]) + dist(y0, q[len - 1]);
  if (lb >= bsf) {
    return lb;
  }

  /// 2 points at front
  x1 = (t[(j + 1)] - mean) / std;
  d = min(dist(x1, q[0]), dist(x0, q[1]));
  d = min(d, dist(x1, q[1]));
  lb += d;
  if (lb >= bsf) {
    return lb;
  }

  /// 2 points at back
  y1 = (t[(len - 2 + j)] - mean) / std;
  d = min(dist(y1, q[len - 1]), dist(y0, q[len - 2]));
  d = min(d, dist(y1, q[len - 2]));
  lb += d;
  if (lb >= bsf) {
    return lb;
  }

  /// 3 points at front
  x2 = (t[(j + 2)] - mean) / std;
  d = min(dist(x0, q[2]), dist(x1, q[2]));
  d = min(d, dist(x2, q[2]));
  d = min(d, dist(x2, q[1]));
  d = min(d, dist(x2, q[0]));
  lb += d;
  if (lb >= bsf) {
    return lb;
  }

  /// 3 points at back
  y2 = (t[(len - 3 + j)] - mean) / std;
  d = min(dist(y0, q[len - 3]), dist(y1, q[len - 3]));
  d = min(d, dist(y2, q[len - 3]));
  d = min(d, dist(y2, q[len - 2]));
  d = min(d, dist(y2, q[len - 1]));
  lb += d;

  return lb;
}

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
                           double best_so_far) {
  double lb = 0;
  double x, d;

  for (int i = 0; i < len && lb < best_so_far; i++) {
    x = (t[(order[i] + j)] - mean) / std;
    d = 0;

    if (x > uo[i]) {
      d = dist(x, uo[i]);
    } else if (x < lo[i]) {
      d = dist(x, lo[i]);
    }

    lb += d;
    cb[order[i]] = d;
  }

  return lb;
}

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
                                double best_so_far) {
  double lb = 0;
  double uu, ll, d;

  for (int i = 0; i < len && lb < best_so_far; i++) {
    uu = (u[order[i]] - mean) / std;
    ll = (l[order[i]] - mean) / std;
    d = 0;

    if (qo[i] > uu) {
      d = dist(qo[i], uu);
    } else {
      if (qo[i] < ll) {
        d = dist(qo[i], ll);
      }
    }

    lb += d;
    cb[order[i]] = d;
  }

  return lb;
}

/// Calculate Dynamic Time Wrapping distance
/// A,B: data and query, respectively
/// cb : cummulative bound used for early abandoning
/// r  : size of Sakoe-Chiba warpping band
double dtw(double* A, double* B, double* cb, int m, int r, double bsf) {
  double* cost;
  double* cost_prev;
  double* cost_tmp;
  int i, j, k;
  double x, y, z, min_cost, final_dtw;

  /// Instead of using matrix of size O(m^2) or O(mr), we will reuse two array
  /// of size O(r).
  cost = (double*)malloc(sizeof(double) * (2 * r + 1));
  for (k = 0; k < 2 * r + 1; k++) {
    cost[k] = INF;
  }

  cost_prev = (double*)malloc(sizeof(double) * (2 * r + 1));
  for (k = 0; k < 2 * r + 1; k++) {
    cost_prev[k] = INF;
  }

  for (i = 0; i < m; i++) {
    k = max(0, r - i);
    min_cost = INF;

    for (j = max(0, i - r); j <= min(m - 1, i + r); j++, k++) {
      /// Initialize all row and column
      if ((i == 0) && (j == 0)) {
        cost[k] = dist(A[0], B[0]);
        min_cost = cost[k];
        continue;
      }

      if ((j - 1 < 0) || (k - 1 < 0)) {
        y = INF;
      } else {
        y = cost[k - 1];
      }

      if ((i - 1 < 0) || (k + 1 > 2 * r)) {
        x = INF;
      } else {
        x = cost_prev[k + 1];
      }

      if ((i - 1 < 0) || (j - 1 < 0)) {
        z = INF;
      } else {
        z = cost_prev[k];
      }

      /// Classic DTW calculation
      cost[k] = min(min(x, y), z) + dist(A[i], B[j]);

      /// Find minimum cost in row for early abandoning (possibly to use column
      /// instead of row).
      if (cost[k] < min_cost) {
        min_cost = cost[k];
      }
    }

    /// We can abandon early if the current cummulative distace with lower bound
    /// together are larger than bsf
    if (i + r < m - 1 && min_cost + cb[i + r + 1] >= bsf) {
      free(cost);
      free(cost_prev);
      return min_cost + cb[i + r + 1];
    }

    /// Move current array to previous array.
    cost_tmp = cost;
    cost = cost_prev;
    cost_prev = cost_tmp;
  }

  k--;

  /// the DTW distance is in the last cell in the matrix of size O(m^2) or at
  /// the middle of our array.
  final_dtw = cost_prev[k];
  free(cost);
  free(cost_prev);

  return final_dtw;
}

/// Log errors
void error(int id) {
  switch (id){
    case 1:
      elog(INFO, "ERROR : Memory can't be allocated!");
      break;
  }
}

/// Search for query in data, returns nearest neighbour distance and location
/// data,query: data and query, respectively
/// di, qi: data and query start indexes
/// dl, ql: data and query lengths
/// ww : warping window (0-1)
/// return  : ssr {double distance, long long location}
ssr search_query(double* data, long long di, long long dl, double* query, long long qi, long long ql, float4 ww) {
  double *t, *q;  /// data array and query array copies for manipulation (z
                  /// normalization etc.)
  double bsf;     /// best-so-far
  int* order;     /// new order of the query
  double *u, *l, *qo, *uo, *lo, *tz, *cb, *cb1, *cb2, *u_d, *l_d;

  double d;
  long long i, j;
  double ex, ex2, mean, std;
  int r = -1;
  long long loc = 0;
  int kim = 0, keogh = 0, keogh2 = 0;
  double dist = 0, lb_kim = 0, lb_k = 0, lb_k2 = 0;
  double *buffer, *u_buff, *l_buff;
  Idx* Q_tmp;

  /// For every EPOCH points, all cummulative values, such as ex (sum), ex2 (sum
  /// square), will be restarted for reducing the floating point error.
  int EPOCH = 100000;

  bool done = false;
  int it = 0, ep = 0, k = 0;
  long long I;  /// the starting index of the data in current chunk of size EPOCH

  ssr result;

  /// increment data length by start offset
  dl += di;

  /// read warping windows
  r = floor(ww * ql);

  /// malloc everything here
  q = (double*)malloc(sizeof(double) * ql);
  if (q == NULL)
    error(1);
  qo = (double*)malloc(sizeof(double) * ql);
  if (qo == NULL)
    error(1);
  uo = (double*)malloc(sizeof(double) * ql);
  if (uo == NULL)
    error(1);
  lo = (double*)malloc(sizeof(double) * ql);
  if (lo == NULL)
    error(1);

  order = (int*)malloc(sizeof(int) * ql);
  if (order == NULL)
    error(1);

  Q_tmp = (Idx*)malloc(sizeof(Idx) * ql);
  if (Q_tmp == NULL)
    error(1);

  u = (double*)malloc(sizeof(double) * ql);
  if (u == NULL)
    error(1);

  l = (double*)malloc(sizeof(double) * ql);
  if (l == NULL)
    error(1);

  cb = (double*)malloc(sizeof(double) * ql);
  if (cb == NULL)
    error(1);

  cb1 = (double*)malloc(sizeof(double) * ql);
  if (cb1 == NULL)
    error(1);

  cb2 = (double*)malloc(sizeof(double) * ql);
  if (cb2 == NULL)
    error(1);

  u_d = (double*)malloc(sizeof(double) * ql);
  if (u == NULL)
    error(1);

  l_d = (double*)malloc(sizeof(double) * ql);
  if (l == NULL)
    error(1);

  t = (double*)malloc(sizeof(double) * ql * 2);
  if (t == NULL)
    error(1);

  tz = (double*)malloc(sizeof(double) * ql);
  if (tz == NULL)
    error(1);

  buffer = (double*)malloc(sizeof(double) * EPOCH);
  if (buffer == NULL)
    error(1);

  u_buff = (double*)malloc(sizeof(double) * EPOCH);
  if (u_buff == NULL)
    error(1);

  l_buff = (double*)malloc(sizeof(double) * EPOCH);
  if (l_buff == NULL)
    error(1);

  /// Copy query array for manipulation, calculate sum and variance
  bsf = INF;
  j = 0;
  ex = ex2 = 0;
  for (i = qi; i < ql + qi; i++, j++) {
    ex += query[i];
    ex2 += query[i] * query[i];
    q[j] = query[i];
  }

  /// Do z-normalize the query, keep in same array, q
  mean = ex / ql;
  std = ex2 / ql;
  std = sqrt(std - mean * mean);
  for (i = 0; i < ql; i++) {
    q[i] = (q[i] - mean) / std;
  }

  /// Create envelop of the query: lower envelop, l, and upper envelop, u
  lower_upper_lemire(q, ql, r, l, u);

  /// Sort the query one time by abs(z-norm(q[i]))
  for (i = 0; i < ql; i++) {
    Q_tmp[i].value = q[i];
    Q_tmp[i].index = i;
  }
  qsort(Q_tmp, ql, sizeof(Idx), comp);

  /// also create another arrays for keeping sorted envelop
  for (i = 0; i < ql; i++) {
    int o = Q_tmp[i].index;
    order[i] = o;
    qo[i] = q[o];
    uo[i] = u[o];
    lo[i] = l[o];
  }
  free(Q_tmp);

  /// Initial the cummulative lower bound
  for (i = 0; i < ql; i++) {
    cb[i] = 0;
    cb1[i] = 0;
    cb2[i] = 0;
  }

  i = 0;  /// current index of the data in current chunk of size EPOCH
  j = 0;  /// the starting index of the data in the circular array, t
  ex = ex2 = 0;

  while (!done) {
    /// Read first m-1 points
    ep = 0;
    if (it == 0) {
      for (k = 0; k < ql - 1; k++) {
          if (di < dl) {
            buffer[k] = data[di];
            ++di;
          }
      }
    } else {
      for (k = 0; k < ql - 1; k++) {
        buffer[k] = buffer[EPOCH - ql + 1 + k];
      }
    }

    /// Read buffer of size EPOCH or when all data has been read.
    ep = ql - 1;
    while (ep < EPOCH && di < dl) {
      buffer[ep] = data[di];
      ++di;
      ep++;
    }

    /// Data are read in chunk of size EPOCH.
    /// When there is nothing to read, the loop is end.
    if (ep <= ql - 1) {
      done = true;
    } else {
      lower_upper_lemire(buffer, ep, r, l_buff, u_buff);

      /// Do main task here..
      ex = 0;
      ex2 = 0;
      for (i = 0; i < ep; i++) {
        /// A bunch of data has been read and pick one of them at a time to use
        d = buffer[i];

        /// Calcualte sum and sum square
        ex += d;
        ex2 += d * d;

        /// t is a circular array for keeping current data
        t[i % ql] = d;

        /// Double the size for avoiding using modulo "%" operator
        t[(i % ql) + ql] = d;

        /// Start the task when there are more than m-1 points in the current
        /// chunk
        if (i >= ql - 1) {
          mean = ex / ql;
          std = ex2 / ql;
          std = sqrt(std - mean * mean);

          /// compute the start location of the data in the current circular
          /// array, t
          j = (i + 1) % ql;
          /// the start location of the data in the current chunk
          I = i - (ql - 1);

          /// Use a constant lower bound to prune the obvious subsequence
          lb_kim = lb_kim_hierarchy(t, q, j, ql, mean, std, bsf);

          if (lb_kim < bsf) {
            /// Use a linear time lower bound to prune; z_normalization of t
            /// will be computed on the fly.
            /// uo, lo are envelop of the query.
            lb_k = lb_keogh_cumulative(order, t, uo, lo, cb1, j, ql, mean, std,
                                       bsf);
            if (lb_k < bsf) {
              /// Take another linear time to compute z_normalization of t.
              /// Note that for better optimization, this can merge to the
              /// previous function.
              for (k = 0; k < ql; k++) {
                tz[k] = (t[(k + j)] - mean) / std;
              }

              /// Use another lb_keogh to prune
              /// qo is the sorted query. tz is unsorted z_normalized data.
              /// l_buff, u_buff are big envelop for all data in this chunk
              lb_k2 = lb_keogh_data_cumulative(order, tz, qo, cb2, l_buff + I,
                                               u_buff + I, ql, mean, std, bsf);
              if (lb_k2 < bsf) {
                /// Choose better lower bound between lb_keogh and lb_keogh2 to
                /// be used in early abandoning DTW
                /// Note that cb and cb2 will be cumulative summed here.
                if (lb_k > lb_k2) {
                  cb[ql - 1] = cb1[ql - 1];
                  for (k = ql - 2; k >= 0; k--)
                    cb[k] = cb[k + 1] + cb1[k];
                } else {
                  cb[ql - 1] = cb2[ql - 1];
                  for (k = ql - 2; k >= 0; k--)
                    cb[k] = cb[k + 1] + cb2[k];
                }

                /// Compute DTW and early abandoning if possible
                dist = dtw(tz, q, cb, ql, r, bsf);

                if (dist < bsf) {  /// Update bsf
                  /// loc is the real starting location of the nearest neighbor
                  /// in the file
                  bsf = dist;
                  loc = (it) * (EPOCH - ql + 1) + i - ql + 1;
                }
              } else
                keogh2++;
            } else
              keogh++;
          } else
            kim++;

          /// Reduce obsolute points from sum and sum square
          ex -= t[j];
          ex2 -= t[j] * t[j];
        }
      }

      /// If the size of last chunk is less then EPOCH, then no more data and
      /// terminate.
      if (ep < EPOCH)
        done = true;
      else
        it++;
    }
  }

  free(q);
  free(u);
  free(l);
  free(uo);
  free(lo);
  free(qo);
  free(cb);
  free(cb1);
  free(cb2);
  free(tz);
  free(t);
  free(l_d);
  free(u_d);
  free(l_buff);
  free(u_buff);

//    elog(INFO, "Data scaned: %d", (it) * (EPOCH - ql + 1) + ep);
//    elog(INFO, "Pruned by LB_Kim: %d", ((double) kim / i)*100));
//    elog(INFO, "Pruned by LB_Keogh: %d", ((double) keogh / i)*100));
//    elog(INFO, "Pruned by LB_Keogh2: %d", ((double) keogh2 / i)*100));
//    elog(INFO, "DTW Calculation: %d", 100 - (((double)kim + keogh + keogh2) / i * 100));

  result.distance = sqrt(bsf); // distance
  result.location = loc; // location

  return result;
}

/// Convert 1D array of postgresql text type to c compatibile array of strings
/// array: 1D array of text type
char** text_array_to_str_array(ArrayType *array) {
  Datum* dimdatums;
  int ndim, i;
  char** result_array;
  char* s;

  Assert(ARR_ELEMTYPE(array) == TEXTOID);

  deconstruct_array(array, TEXTOID, -1, false, 'i', &dimdatums, NULL, &ndim);

  result_array = (char**)palloc(ndim * sizeof(char*));

  for (i = 0; i < ndim; i++) {
    s = TextDatumGetCString(dimdatums[i]);
    result_array[i] = s;
  }

  return result_array;
}

/// PostgreSQL function for evaluating similarity between multivariate signals
Datum signal_similarity(PG_FUNCTION_ARGS) {
  ArrayType* atqn = PG_GETARG_ARRAYTYPE_P(0);  /// query series names
  ArrayType* atqv = PG_GETARG_ARRAYTYPE_P(1);  /// query series values
  ArrayType* atdn = PG_GETARG_ARRAYTYPE_P(2);  /// data series names
  ArrayType* atdv = PG_GETARG_ARRAYTYPE_P(3);  /// data series values
  float4 ww = PG_GETARG_FLOAT4(4);  /// window warping (0-1)

  char **qn = text_array_to_str_array(atqn);  /// query series names
  double *qv = ARRDOUBLEPTR(atqv);  /// query series values
  char **dn = text_array_to_str_array(atdn);  /// data series names
  double *dv = ARRDOUBLEPTR(atdv);  /// data series values

  double similarity = 0;
  int coverage = 0; /// how many series in query have equivalent (by name) in data
  int i, j;

  ssr *results;
  Datum *datresult;
  ArrayType *result;

  /// get query and data series number
  int qs = ARR_DIMS(atqn)[0];
  int ds = ARR_DIMS(atdn)[0];

  /// get query and data series length
  int dl = ARR_DIMS(atdv)[1];
  int ql = ARR_DIMS(atqv)[1];

  results = (ssr*)malloc(sizeof(ssr)*qs);
  for (i = 0; i < qs; i++) {
    for (j = 0; j < ds; j++) {
      if (strcmp(qn[i], dn[j]) == 0) {
        results[i] = search_query(dv, j * dl, dl, qv, i * ql, ql, ww);
        similarity += results[i].distance;
        coverage++;
        break; /// ommit query series that already have match in data
      } else {
        results[i].distance = -1;
        results[i].location = -1;
      }
    }
  }

  similarity = (coverage == 0) ? 0 : 1 / (1 + (similarity / coverage));

  datresult = (Datum*)palloc(sizeof(Datum) * (qs * 2 + 1));
  datresult[0] = Float8GetDatum(similarity);
  j = 1;
  for (i = 0; i < qs; i++) {
    datresult[j++] = Float8GetDatum(results[i].distance);
  }
  for (i = 0; i < qs; i++) {
    datresult[j++] = Float8GetDatum((double)results[i].location);
  }

  result = construct_array(datresult, qs*2+1, FLOAT8OID, sizeof(float8), true, 'i');

  free(results);
  pfree(datresult);

  /* Avoid leaking memory when handed toasted input. */
  AARR_FREE_IF_COPY(atqn, 0);
  AARR_FREE_IF_COPY(atqv, 1);
  AARR_FREE_IF_COPY(atdn, 2);
  AARR_FREE_IF_COPY(atdv, 3);

  /**
   * Return array with structure as
   * {similarity, query distance for series 1, 2..., location (index) of matched query on series 1, 2...}
   */
  PG_RETURN_ARRAYTYPE_P(result);
}
