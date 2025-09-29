#pragma once

#include <algorithm>
#include <vector>
#include <stack>
#include <map>
#include <sstream>
#include <iterator>
#include <cstdio>
#include <mpi.h>

// NOTE: Do not use RUNTIME and UNSYNCNESS at the same time!
enum DbgMeasureMode { OFF=0, RUNTIME=1, UNSYNCNESS=2, ROOT_STAT=4, TOTAL_STAT=8, };

static std::string join(std::vector<std::string> const &strings, std::string delim) {
    std::stringstream ss;
    std::copy(strings.begin(), strings.end(),
        std::ostream_iterator<std::string>(ss, delim.c_str()));
    return ss.str();
}

struct CommInfo {
  int me, root, n_ranks;
  MPI_Comm comm;
  bool i_am_root() { return me == root; }
  static CommInfo getDefault(MPI_Comm comm = MPI_COMM_WORLD) {
    int me, n_ranks;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &n_ranks);
    return {me, 0, n_ranks, comm};
  }
};


struct TimeMeasurer {
  using strvec = std::vector<std::string>;

  FILE* fp = nullptr;
  size_t flagMask;
  CommInfo ci;

  strvec header;
  bool root_stat, total_stat;

  enum StatKind { MIN, AVG, MAX, SUM };
  std::stack<double> start_times;
  double diff_time, time_stat[3];
  std::vector<double[4][3]> total_time_stat;
  std::vector<size_t> n_timesteps;
  size_t col_cnt;
  
  TimeMeasurer(std::string fname, size_t flagMask, strvec header, MPI_Comm comm = MPI_COMM_WORLD)
      : flagMask{flagMask}, header{header}, ci{CommInfo::getDefault(comm)} {
    total_stat = flagMask & DbgMeasureMode::TOTAL_STAT;
    root_stat = (flagMask & DbgMeasureMode::ROOT_STAT) or total_stat;
    if (header.empty() and total_stat) { fprintf(stderr, "ERROR: incompatible: total_stat without header given\n"); exit(1); }
    if (total_stat) {
      total_time_stat = std::vector<double[4][3]>(header.size());
      for (size_t comp = 0; comp < total_time_stat.size(); ++comp)
        for (size_t row = StatKind::MIN; row <= StatKind::SUM; ++row)
          for (size_t col = StatKind::MIN; col <= StatKind::MAX; ++col)
            total_time_stat[comp][row][col] = (row == StatKind::MIN) ? 1.0 : 0.0;
      n_timesteps = std::vector<size_t>(header.size(), 0);
    }
    if ((not root_stat) or ci.i_am_root()) {
      fp = fopen(fname.c_str(), "w");
      fprintf(fp, "%s\n", join(header, " ").c_str());
    }
    col_cnt = 0;
  }
  TimeMeasurer(std::string fname, size_t flagMask, strvec header, CommInfo ci)
      : TimeMeasurer(fname, flagMask, header) {
    ci = ci;
  }

  virtual ~TimeMeasurer() {
    if (fp) {
      if (total_stat and ci.i_am_root()) {
        newline();

        for (size_t comp = 0; comp < total_time_stat.size(); ++comp) {
          fprintf(fp, "[n:%lu]  ", n_timesteps[comp]);
          for (size_t col = StatKind::MIN; col <= StatKind::MAX; ++col)
            total_time_stat[comp][StatKind::AVG][col] = total_time_stat[comp][StatKind::SUM][col] / n_timesteps[comp];
        } newline();

        for (size_t comp = 0; comp < total_time_stat.size(); ++comp) {
          printStat(total_time_stat[comp][StatKind::SUM], '[', ']');
        } newline();

        for (size_t row = StatKind::MIN; row <= StatKind::MAX; ++row) {
          for (size_t comp = 0; comp < total_time_stat.size(); ++comp) {
            printStat(total_time_stat[comp][row], '[', ']');
          } newline();
        }
      }
      fclose(fp);
      fp = nullptr;
    }
  }

  virtual void before() =0;
  virtual void after(std::string prefix) =0;

  double get() {
    return diff_time;
  }

  void write(std::string prefix) {
    if (not fp) return;
    if (header.empty()) {   // no header -> print name and value
      if (prefix != "") print(prefix + ":");
    } else {                // header -> pad if neccessary, then only print value
      std::string na = "   NA   ";
      if (root_stat) na = "(" + na + " " + na + " " + na + ")  ";
      while (prefix != "" and prefix != header.at(col_cnt++)) print(na);
    }
    if (root_stat) printStat(time_stat);
    else printDbl(diff_time);
  }

  void collect(std::string prefix) {
    MPI_Reduce(&diff_time, &time_stat[StatKind::MIN], 1, MPI_DOUBLE, MPI_MIN, ci.root, ci.comm);
    MPI_Reduce(&diff_time, &time_stat[StatKind::AVG], 1, MPI_DOUBLE, MPI_SUM, ci.root, ci.comm);
    MPI_Reduce(&diff_time, &time_stat[StatKind::MAX], 1, MPI_DOUBLE, MPI_MAX, ci.root, ci.comm);
    time_stat[StatKind::AVG] /= ci.n_ranks;
    if (total_stat and ci.i_am_root()) {
      // DONE?: make it work
      int comp = index_in_header(prefix);
      if (comp == -1) { fprintf(stderr, "ERROR: invalid prefix: '%s'\n", prefix.c_str()); exit(1); }
      for (size_t col = StatKind::MIN; col <= StatKind::MAX; ++col) {
        total_time_stat[comp][StatKind::MIN][col] = std::min(total_time_stat[comp][StatKind::MIN][col], time_stat[col]);
        total_time_stat[comp][StatKind::MAX][col] = std::max(total_time_stat[comp][StatKind::MAX][col], time_stat[col]);
        total_time_stat[comp][StatKind::SUM][col] += time_stat[col];
      }
      ++n_timesteps[comp];
    }
  }

  void print(std::string str)     { if (fp) fprintf(fp, "%s", str.c_str()); }
  void printLine(std::string str) { if (fp) fprintf(fp, "%s\n", str.c_str()); col_cnt = 0; }
  void printDbl(double d)         { if (fp) fprintf(fp, "%f ", d); }
  void printStat(double* arr, char ob='(', char cb=')') 
                                  { if (fp) fprintf(fp, "%c%f %f %f%c  ", ob, arr[0], arr[1], arr[2], cb); }
  void newline()                  { if (fp) fprintf(fp, "\n"); col_cnt = 0; }
  void space(size_t n=2)          { if (fp) for (size_t i = 0; i < n; ++i) fprintf(fp, " "); }
  void flush()                    { if (fp) fflush(fp); }

  int index_in_header(std::string s) {
    for (int i = 0; i < header.size(); ++i) if (s == header[i]) return i;
    return -1;
  }
};


struct RunTimeMeasurer : TimeMeasurer {
  RunTimeMeasurer(std::string fname, size_t flagMask, strvec header, CommInfo ci)
      : TimeMeasurer(fname, flagMask, header, ci) { }
  
  void before() override { start_times.push(MPI_Wtime()); }
  void after(std::string prefix) override {
    diff_time = MPI_Wtime() - start_times.top();
    start_times.pop();
    if (root_stat) collect(prefix);
  }
};

struct SyncTimeMeasurer : TimeMeasurer {
  SyncTimeMeasurer(std::string fname, size_t flagMask, strvec header, CommInfo ci)
      : TimeMeasurer(fname, flagMask, header, ci) { }
  
  void before() override { MPI_Barrier(ci.comm); }
  void after(std::string prefix) override {
    double start_time = MPI_Wtime();
    MPI_Barrier(ci.comm);
    diff_time = MPI_Wtime() - start_time;
    if (root_stat) collect(prefix);
  }
};

static TimeMeasurer* make_TimeMeasurer(std::string fname, size_t flagMask, std::vector<std::string> header, CommInfo ci) {
  if (flagMask & DbgMeasureMode::RUNTIME)         return new RunTimeMeasurer(fname, flagMask, header, ci);
  else if (flagMask & DbgMeasureMode::UNSYNCNESS) return new SyncTimeMeasurer(fname, flagMask, header, ci);
  else                                            return nullptr;
}

//#define DEBUG_MEASURE_ENABLED true
#if defined(DEBUG_MEASURE_ENABLED) && DEBUG_MEASURE_ENABLED == true

  #pragma message "DEBUG_MEASURE enabled!"

  #define DEBUG_MEASURE_ATTRS(tm_ptr) \
    TimeMeasurer* tm_ptr

  #define DEBUG_MEASURE_SETUP(tm_ptr, flagMask, ci, fname_base, fname_func, ... /* header list */) do { \
      char fname[128]; sprintf(fname, "%s-%s-%03i.log", fname_base, fname_func, ci.me); \
      tm_ptr = make_TimeMeasurer(fname, flagMask, __VA_ARGS__, ci); \
    } while (false)

  #define DEBUG_MEASURE_DESTROY(tm_ptr) do { \
      if (tm_ptr) delete tm_ptr; \
    } while (false)

  #define DEBUG_MEASURE_EOL(tm_ptr, flags) do { \
      if (tm_ptr and flags != DbgMeasureMode::OFF) tm_ptr->newline(); \
    } while (false)

  #define DEBUG_MEASURE(tm_ptr, flags, prefix, code_block) do { \
      if (tm_ptr and flags != DbgMeasureMode::OFF ) { \
        tm_ptr->before(); \
        { code_block } \
        tm_ptr->after(prefix); \
        tm_ptr->write(prefix); \
      } else { \
        code_block \
      } \
    } while (false)

// TODO: accumulative measure option

#else   // !defined(DEBUG_MEASURE_ENABLED) || DEBUG_MEASURE_ENABLED == false

  #pragma message "DEBUG_MEASURE disabled!"

  #define DEBUG_MEASURE_ATTRS(tm_ptr) ;
  #define DEBUG_MEASURE_SETUP(tm_ptr, flagMask, ci, fname_base, fname_func, ... /* header list */) {} 
  #define DEBUG_MEASURE_DESTROY(tm_ptr) {}
  #define DEBUG_MEASURE_EOL(tm_ptr, flags) {}
  #define DEBUG_MEASURE(tm_ptr, flags, prefix, code_block) code_block

#endif  // defined(DEBUG_MEASURE_ENABLED) && DEBUG_MEASURE_ENABLED == true

