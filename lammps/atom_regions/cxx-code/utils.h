#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cstdarg>

std::string to_wdtstr(int i, int wdt) {
  std::stringstream sstr;
  sstr << std::setw(wdt) << std::setfill('0') << i;
  return sstr.str();
}

void error_msg(std::string fmt, ...) {
  std::fprintf(stderr, "ERROR: ");
  std::va_list args;
  va_start(args, fmt);
  std::vfprintf(stderr, fmt.data(), args);
  va_end(args);
  std::exit(1);
}

std::vector<std::string> toStrVec(int argc, char* argv[]) {
  std::vector<std::string> res(argc);
  std::transform(argv, argv+argc, res.begin(), [](const char* s) { return std::string(s); });
  return res;
}
