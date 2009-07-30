#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/config.hpp>

using namespace boost;

#define MIN_OVERLAP 3

int main()
{
  std::map<std::string, int> edge_counts;
  std::map<std::string, int> doc_node_set;
  std::map<std::string, std::vector<std::string> > doc_nodes;

  std::ifstream datafile("../../input/data.txt");
  //std::ifstream datafile("data.txt");

  if (!datafile) {
    std::cerr << "No data.txt file" << std::endl;
    return EXIT_FAILURE;
  }

  int j=0;

  for (std::string line; std::getline(datafile, line);) {
    if (line.rfind(":") == std::string::npos)
      continue;

    char_delimiters_separator < char >sep(false, "", ":");
    tokenizer <> line_toks(line, sep);
    tokenizer <>::iterator i = line_toks.begin();
    std::string user = *i++;
    std::string doc = *i;
    std::vector<std::string>::iterator it;
    std::vector<std::string> users = doc_nodes[doc];

    for (it = users.begin(); it < users.end(); it++) {
      std::stringstream edge;
      int cmp = user.compare(*it);
      if (!cmp)
        continue;
      if (cmp < 0) {
        edge << *it << ":" << user;
      }
      else {
        edge << user << ":" << *it;
      }

      edge_counts[edge.str().c_str()] ++;
    }

    users.push_back(user);
    doc_nodes[doc] = users;
    if ((++j) % 10000 == 0) {
      std::cerr << "Finished with line " << j << std::endl;
    }
  }

  for (std::map<std::string, int>::iterator it = edge_counts.begin();
       it != edge_counts.end(); it ++) {
    if (it->second < 3)
      continue;
    std::cout << it->first << ":" << it->second << std::endl;
  }

}
