#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <random>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>

using namespace std;

#define MIN_DISTANCE 0.025

/*
float rnd(float min=0.0, float max=1.0) {
  return (rand() * 1.0 / RAND_MAX) * (max - min) + min;
}
*/


int main(int argc, char* argv[]) {
  std::random_device rnd_dev;
  int seed = rnd_dev();

  int nargs = argc-1;
  if (nargs < 3) {
    std::cout << "ERROR: at least three parameters (proc grid dims) required (plus optional random seed)!" << std::endl; std::exit(1);
  } else if (nargs > 4) {
    std::cout << "ERROR: too many arguments: three or four parameters required!" << std::endl; std::exit(1);
  } else if (nargs == 4) {
    seed = std::atoi(argv[nargs]);
  }

  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> runif(MIN_DISTANCE, 1-MIN_DISTANCE);

  std::vector<std::string> comp_names {"x", "y", "z"};
  std::stringstream output;

  for (int comp = 0; comp < 3; ++comp) {
    output << comp_names[comp] << " ";
    int n = std::atoi(argv[comp+1]);
    std::vector<float> cuts(n-1);
    bool valid;
    do {
      std::generate(cuts.begin(), cuts.end(), [&](){ return runif(rng); });
      std::sort(cuts.begin(), cuts.end());
      // check that cuts are not too close ; otherwise regenerate last cut
      valid = true;
      for (int i = 0; i < cuts.size()-1; ++i)
        if (std::abs(cuts[i] - cuts[i+1]) < MIN_DISTANCE)
          valid = false;
    } while (!valid);
    std::copy(cuts.begin(), cuts.end(), std::ostream_iterator<float>(output, " "));
  }

  std::cout << output.str() << std::endl;

  return 0;
}
