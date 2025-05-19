#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cstdint>

constexpr int N_DIMS = 3;

std::string to_wdtstr(int i, int wdt) {
  std::stringstream sstr;
  sstr << std::setw(wdt) << std::setfill('0') << i;
  return sstr.str();
}

template<typename... Args>
void error_msg(const char* fmt, Args... args) {
  printf(fmt, args...);
  exit(1);
}

std::vector<std::string> toStrVec(int argc, char* argv[]) {
  std::vector<std::string> res(argc);
  std::generate(res.begin(), res.end(), [n=0]() mutable { return std::string(argv[n]); });
  return res;
}

enum class GridStyle {
  Tensor,
  Staggered
};

template<class T>
class Triple {
  T data[N_DIMS];
  T& x = data[0];
  T& y = data[1];
  T& z = data[2];
  Triple() : x(0), y(0), z(0) { }
  Triple(T x, T y, T z) : x{x}, y{y}, z{z} { }
  T operator[](int i) const { return data[i]; }
  T& operator[](int i) { return data[i]; }
  T sum() { return x + y + z; }
  T prod() { return x * y * z; }
};


struct Cell {
  double borders[2*N_DIMS];
  double& x_low  = borders[0];
  double& x_high = borders[1];
  double& y_low  = borders[2];
  double& y_high = borders[3];
  double& z_low  = borders[4];
  double& z_high = borders[5];
  std::string name;
  Cell(std::string name) : name{name} { }
  double operator[](int i) const { return data[i]; }
  double& operator[](int i) { return data[i]; }
};

struct Cuts {
  // interface: x, y, z
  // internal:  z, y, x
  bool zero_one_scale;
  Triple<int> dim;
  GridStyle grid;
  std::vector<double> z_cuts;
  std::vector<double> y_cuts;
  std::vector<double> x_cuts;
  Cuts(int x, int y, int z, GridStyle grid) : dim(x, y, z), grid{grid} {
    zero_one_scale = true;
    if (grid == GridStyle::Tensor) {
        z_cuts = std::vector<double>(z+1, 0);
        y_cuts = std::vector<double>(y+1, 0);
        x_cuts = std::vector<double>(x+1, 0);
    } else if (grid == GridStyle::Staggered) {
        z_cuts = std::vector<double>(z+1, 0);
        y_cuts = std::vector<double>(z * (y+1), 0);
        x_cuts = std::vector<double>(z * y * (x+1), 0);
    }
  }

  Cuts(int x, int y, int z, GridStyle grid, std::vector<double> zyx_cut_params, bool zero_one_scale = true) : Cuts(x, y, z, grid), zero_one_scale{zero_one_scale} {
    assert(z_cuts.size() + y_cuts.size() + x_cuts.size() == zyx_cut_params.size());
    std::copy(z_cuts.begin(), z_cuts.end(), zyx_cut_params.begin());
    std::copy(y_cuts.begin(), y_cuts.end(), zyx_cut_params.begin()+z_cuts.size());
    std::copy(x_cuts.begin(), x_cuts.end(), zyx_cut_params.begin()+z_cuts.size()+y_cuts.size());
  }

  template<class T>
  Cuts operator*(Triple<T> lens) {
    if (not zero_one_scale) error_msg("The Cuts are not on a 0-1 scale and can't be rescaled anymore!");
    Cuts res(*this);
    for (double& val : res.z_cuts) val *= lens.z;
    for (double& val : res.y_cuts) val *= lens.y;
    for (double& val : res.x_cuts) val *= lens.x;
    res.zero_one_scale = false;
    return res;
  }

  friend std::vector<double> complete_zero_one_scale(std::vector<double> zyx_cut_params, GridStyle grid) {
    // TODO: make a minimalistic cut-list complete
    //      otherwise only accept already complete cut lists
    return zyx_cut_params;
  }
};

struct System {
  Triple<int> procs, lens;
  Cuts cuts;
  std::vector<Cell> cells;      // len: procs.prod()

  // TODO: add cuts (optionally?)
  // TODO: remove args-parsing
  System(char* args[], double region_gap = 0) : procs(args), lens(args+N_DIMS) {
    // calculate ideal, arithmetic cuts (= borders of regions)
    for (int d = 0; d < N_DIMS; ++d) {
      DimCuts dc(procs[d]);
      std::generate(dc.data.begin(), dc.data.end(), [&,i=0]() mutable { return lens[d] * (i++ * 1.0 / procs[d]); } );
      dimcuts.push_back(dc);
    }
    // calculate borders for each region
    for (int reg = 0; reg < procs.prod(); ++reg) {
      Coord coord = idx_to_coord(reg);
      Cell cell("blk" + to_wdtstr(reg, 2));
      for (int d = 0; d < N_DIMS; ++d) {
        cell[2*d]     = dimcuts[d][coord[d]]     + region_gap/2;  // cell.x_low  = dimcuts[0][coord.x];
        cell[2*d + 1] = dimcuts[d][coord[d] + 1] - region_gap/2;  // cell.x_high = dimcuts[0][coord.x + 1];
      }
      cells.push_back(cell);
    }
  }

  System(Triple<int> procs, Triple<int> lens, Cuts cuts, double region_gap = 0) : procs{procs}, lens{lens}, cuts{cuts} {
    // scale cuts
    if (cuts.zero_one_scale) cuts = cuts * lens;
    // calculate cell-borders
    for (int reg = 0; reg < procs.prod(); ++reg) {
      Triple<int> coord = idx_to_coord(reg);
      Cell cell("reg" + to_wdtstr(reg, 2));
      cell.z_low  = cuts.z_cuts[coord.z]     + region_gap/2;
      cell.z_high = cuts.z_cuts[coord.z + 1] - region_gap/2;
      if (cuts.grid == GridStyle::Tensor) {
        cell.y_low  = cuts.y_cuts[coord.y]     + region_gap/2;
        cell.y_high = cuts.y_cuts[coord.y + 1] - region_gap/2;
        cell.x_low  = cuts.x_cuts[coord.x]     + region_gap/2;
        cell.x_high = cuts.x_cuts[coord.x + 1] - region_gap/2;
      } else if (cuts.grid == GridStyle::Staggered) {
        cell.y_low  = cuts.y_cuts[coord.z * (procs.y + 1) + coord.y]     + region_gap/2;
        cell.y_high = cuts.y_cuts[coord.z * (procs.y + 1) + coord.y + 1] - region_gap/2;
        cell.x_low  = cuts.x_cuts[coord.z * (procs.y + 1) * (procs.x + 1) + coord.y * (procs.x + 1) + coord.x]     + region_gap/2;
        cell.x_high = cuts.x_cuts[coord.z * (procs.y + 1) * (procs.x + 1) + coord.y * (procs.x + 1) + coord.x + 1] - region_gap/2;
      }
      cells.push_back(cell);
    }
  }

  std::string lmp_style_regions() {
    // "region  <region-name>  <kind=>block  <x-low> <x-high> <y-low> <y-high> <z-low> <z-high>"
    std::stringstream lmp_cmd;
    for (auto& cell : cells) {
      lmp_cmd << "region  " << cell.name << "  block  ";
      std::copy(cell.data, cell.data+6, std::ostream_iterator<double>(lmp_cmd, "  "));
      lmp_cmd << std::endl;
    }
    return lmp_cmd.str();
  }

  std::string lmp_style_create(int n_particles_per_region=1000, int seed=873984, bool inc_seed=true) {
    // "create_atoms  <atom_kind=1>  random  <n-atoms> <seed> <region-name>"
    std::stringstream lmp_cmd;
    for (auto& cell : cells) {
      lmp_cmd << "create_atoms  1  random  " << n_particles_per_region << " " << seed << " " << cell.name << std::endl;
      if (inc_seed) seed += 1;
    }
    return lmp_cmd.str();
  }

  int coord_to_idx(Coord c) {
    return c.z * (procs.y * procs.x) + c.y * procs.x + c.x;
  }

  Triple<int> idx_to_coord(int idx) {
    Triple<int> coord;
    coord.x = idx % procs.x;
    coord.y = (idx % (procs.x * procs.y)) / procs.x;
    coord.z = idx / (procs.x * procs.y);
    return coord;
  }
};




/**
 *  interface: ./atom_regions <grid>  <procs-x> <procs-y> <procs-z>  <len-x> <len-y> <len-z>  <region-gap>  <cut-parameters>
 *    grid:                       tensor or staggered (string)
 *    procs-x, procs-y, procs-z:  number of processors in direction x, y and z (int)
 *    len-x, len-y, len-z:        length of the system in dimension x, y and z (int|float)
 *    region-gap:                 gap between regions
 *    cut-parameters:             comma separated list of cuts in order of z, y, x (float in ]0,1[)
 *                                  for tensor   grid exactly  (procs-z - 1) + (procs-y - 1) + (procs-x - 1)  values
 *                                  for staggerd grid exactly  (procs-z - 1) + procs-z * (procs-y - 1) + procs-z * procs-y * (procs-x - 1)  values
 */
int main(int argc, char* argv[]) {
  // TODO: add option to choose between tensor and staggered grid
  // check number of arguments
  int argn = 1 + 2 * N_DIMS + 2;
  if (argc-1 != argn)
    error_msg("ERROR: This program requires %i arguments!\n" \
      "./atom_regions <grid>  <procs-x> <procs-y> <procs-z>  <len-x> <len-y> <len-z>  <region-gap>  <cut-params>\n" \
      "  grid:                      tensor or staggered (string)\n" \
      "  procs-x, procs-y, procs-z: number of processors in direction x, y and z (int)\n" \
      "  len-x, len-y, len-z:       length of the system in dimension x, y and z (int|float)\n" \
      "  region-gap:                gap between regions\n" \
      "  cut-parameters:            comma separated list of cuts in order of z, y, x (float in ]0,1[)\n" \
      "                               tensor   grid: exactly (procs-z - 1) + (procs-y - 1) + (procs-x - 1) values\n" \
      "                               staggerd grid: exactly (procs-z - 1) + procs-z * (procs-y - 1) + procs-z * procs-y * (procs-x - 1) values\n",
      argn);


  // parse arguments
  std::vector<std::string> args = toStrVec(argc-1, argv+1);
  int cnt = 0;
  GridStyle grid;
  if (args[cnt] == "tensor")         grid = GridStyle::Tensor;
  else if (args[cnt] == "staggered") grid = GridStyle::Staggered;
  else  error_msg("ERROR: Unknown grid style %s\n", argv);
  ++cnt;

  double region_gap = 0.0; // 0.001
  if (argc-1 == 2 * N_DIMS + 1) {
    region_gap = std::atof(argv[2*N_DIMS+1]);
  }

  // parse arguments
  System sys(argv+1, region_gap);

  std::cout << sys.lmp_style_regions() << std::endl;

  std::cout << sys.lmp_style_create() << std::endl;

  return 0;
}
