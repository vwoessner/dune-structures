#ifndef DUNE_BLOCKLAB_UTILITIES_STRINGSPLIT_HH
#define DUNE_BLOCKLAB_UTILITIES_STRINGSPLIT_HH

#include<algorithm>
#include<string>
#include<vector>


namespace Dune::BlockLab {

  // String trimming - C++ standard library, you failed me!
  // Taken from here: https://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
  static inline void string_trim(std::string &s) {
	  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		  return !std::isspace(ch);
		  }));
	  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		  return !std::isspace(ch);
	      }).base(), s.end());
  }

  // String splitting - C++ standard library, you failed me!
  // Taken from here: http://www.martinbroadhurst.com/how-to-split-a-string-in-c.html
  std::vector<std::string> string_split(const std::string& str, char delim=',', bool trim=true)
  {
    if (str == "")
      return {};

    std::vector<std::string> cont;
    std::size_t current, previous = 0;
    current = str.find(delim);
    while (current != std::string::npos)
    {
      cont.push_back(str.substr(previous, current - previous));
      previous = current + 1;
      current = str.find(delim, previous);
    }
    cont.push_back(str.substr(previous, current - previous));

    // Trim strings if requested
    if(trim)
      for(auto& s : cont)
	string_trim(s);

    return cont;
  }

} // namespace Dune::BlockLab

#endif
