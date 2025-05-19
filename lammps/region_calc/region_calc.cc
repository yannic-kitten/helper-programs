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

struct Coord {
  int data[N_DIMS];
  int& x = data[0];
  int& y = data[1];
  int& z = data[2];
  Coord() { }
  Coord(char* args[]) { for (int d = 0; d < N_DIMS; ++d) data[d] = std::atoi(args[d]); }
  int operator[](int i) const { return data[i]; }
  int& operator[](int i) { return data[i]; }
};

struct Dim : Coord {
  int prod() { return x * y * z; }
  Dim(char* args[]) : Coord(args) { }
};

struct DimCuts {
  std::vector<double> data;
  DimCuts(int n) : data(n+1, 0) { }
  double operator[](int i) const { return data[i]; }
  double& operator[](int i) { return data[i]; }
};

struct Cell {
  double data[6];
  double& x_low  = data[0];
  double& x_high = data[1];
  double& y_low  = data[2];
  double& y_high = data[3];
  double& z_low  = data[4];
  double& z_high = data[5];
  std::string name;
  Cell(std::string name) : name{name} { }
  double operator[](int i) const { return data[i]; }
  double& operator[](int i) { return data[i]; }
};

struct System {
  Dim procs;
  Dim lens;
  std::vector<DimCuts> dimcuts; // len: N_DIMS
  std::vector<Cell> cells;      // len: procs.prod()

  System(char* args[], double border_space = 0) : procs(args), lens(args+N_DIMS) {
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
        cell[2*d]     = dimcuts[d][coord[d]]     + border_space;  // cell.x_low  = dimcuts[0][coord.x];
        cell[2*d + 1] = dimcuts[d][coord[d] + 1] - border_space;  // cell.x_high = dimcuts[0][coord.x + 1];
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

  Coord idx_to_coord(int idx) {
    Coord coord;
    coord.x = idx % procs.x;
    coord.y = (idx % (procs.x * procs.y)) / procs.x;
    coord.z = idx / (procs.x * procs.y);
    return coord;
  }
};



int main(int argc, char* argv[]) {
  // check number of arguments
  if (argc-1 < 2 * N_DIMS or argc-1 > 2 * N_DIMS + 1) {
    std::cout << "ERROR: This program requires " << (2 * N_DIMS) << " or " << (2 * N_DIMS + 1) << " arguments!" << std::endl;
    exit(1);
  }

  double border_space = 0.0; // 0.001
  if (argc-1 == 2 * N_DIMS + 1) {
    border_space = std::atof(argv[2*N_DIMS+1]);
  }

  // parse arguments
  System sys(argv+1, border_space);

  std::cout << sys.lmp_style_regions() << std::endl;

  std::cout << sys.lmp_style_create() << std::endl;

  return 0;
}
