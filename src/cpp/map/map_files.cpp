#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/config.hpp>
#include <cstring>
#include <sstream>

#include "generate_dims.h"

using namespace std;

void map_file(istream *datafile, DimInfo *dims);

int main(int argc, char **argv)
{
  ifstream datafile("../../../input/data.txt");

  /*
    if (argc < 2) {
    fprintf(stderr, "Usage: %s (-p|-r)\n", argv[0]);
    return 2;
    }*/

  DimInfo * dims;


  if (argc > 1 && !strncmp(argv[1], "-c", 25)) {
    dims = generate_dims(&datafile);
    ofstream outuserfile("../../../dat/usermap.dat");
    ofstream outrepofile("../../../dat/repomap.dat");
    ofstream outcountfile("../../../dat/usercount.dat");

    dimmap::iterator it;
    for (it = dims->user_to.begin(); it != dims->user_to.end(); it++) {
      outuserfile << (*it).first << ":" << (*it).second << endl;
    }

    for (it = dims->repos_to.begin(); it != dims->repos_to.end(); it++) {
      outrepofile << (*it).first << ":" << (*it).second << endl;
    }

    intmap::iterator it2;
    for (it2 = dims->user_counts.begin(); it2 != dims->user_counts.end(); it2++) {
      outcountfile << (*it2).first << ":" << (*it2).second << endl;
    }

    outcountfile.close();
    outuserfile.close();
    outrepofile.close();
  }
  else {
    dims = read_dims();
  }

  map_file(&cin, dims);

  destroy_dims(dims);
  return 0;
}


void map_file(istream *datafile, DimInfo *dims)
{
  for (string line; getline(*datafile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string user = *i++;
    string repoinfo = *i;

    boost::char_delimiters_separator < char >sep2(false, "", ",");
    boost::tokenizer <> line_toks2(repoinfo, sep2);
    boost::tokenizer <>::iterator i2;

    int size = 0;
    ostringstream repoids;

    for (i2 = line_toks2.begin(); i2 != line_toks2.end(); i2++) {
      int repoid = dims->repos_to[*i2];
      if (!dims->repos_to.count(*i2)) {
        continue;
      }
      if (!size) {
        repoids << repoid;
      }
      else {
        repoids << "," << repoid;
      }
      size = 1;
    }

    int exists = dims->user_to.count(user);
    int userid = dims->user_to[user];
    if (!exists) {
      fprintf(stderr, "Missing mapping for user %s\n", user.c_str());
      continue;
    }
    if (!size) {
      fprintf(stderr, "Unable to map repo data: %s\n", repoinfo.c_str());
      continue;
    }

    cout << userid << ":" << repoids.str() << endl;
  }
}
