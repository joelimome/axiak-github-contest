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
void map_file_reverse(istream *datafile, DimInfo *dims);
void help(const char * name);

int main(int argc, char **argv)
{
  ifstream datafile("../../../input/data.txt");

  /*
    if (argc < 2) {
    fprintf(stderr, "Usage: %s (-p|-r)\n", argv[0]);
    return 2;
    }*/

  DimInfo * dims;

  string arg = "";
  if (argc > 1) {
    arg = string(argv[1]);
  }

  if (!arg.compare("-h")) {
    help(argv[0]);
  }

  if (!arg.compare("-c")) {
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
  }
  else if (!arg.compare("-r")) {
    dims = read_dims();
    map_file_reverse(&cin, dims);
  }
  else {
    dims = read_dims();
    map_file(&cin, dims);
  }

  destroy_dims(dims);
  return 0;
}


void map_file(istream *datafile, DimInfo *dims)
{
  for (string line; getline(*datafile, line);) {
    if (line.rfind(":") == string::npos) {
      if (dims->user_to.count(line)) {
        cout << dims->user_to[line] << endl;
      }
      continue;
    }

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

void map_file_reverse(istream *datafile, DimInfo *dims)
{

  for (string line; getline(*datafile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    int userid = atoi((*i++).c_str());
    string repoinfo = *i;

    boost::char_delimiters_separator < char >sep2(false, "", ",");
    boost::tokenizer <> line_toks2(repoinfo, sep2);
    boost::tokenizer <>::iterator i2;

    int size = 0;
    ostringstream repoids;

    for (i2 = line_toks2.begin(); i2 != line_toks2.end(); i2++) {
      int repoid = atoi((*i2).c_str());
      if (!dims->repos_from.count(repoid)) {
        continue;
      }
      string repo = dims->user_from[repoid];
      if (!size) {
        repoids << repo;
      }
      else {
        repoids << "," << repo;
      }
      size = 1;
    }

    int exists = dims->user_from.count(userid);
    if (!exists) {
      fprintf(stderr, "Missing mapping for user %d\n", userid);
      continue;
    }
    if (!size) {
      fprintf(stderr, "Unable to map repo data: %s\n", repoinfo.c_str());
      continue;
    }

    cout << dims->user_from[userid] << ":" << repoids.str() << endl;
  }
}

void help(const char * name)
{
  fprintf(stderr, "Usage: %s [-c|-r]\n", name);
  fprintf(stderr, "  Default: Maps input file from username to userid.\n");
  fprintf(stderr, "  -c: Create all the map files.\n");
  fprintf(stderr, "  -r: Map input file with reverse maps.\n");
}
