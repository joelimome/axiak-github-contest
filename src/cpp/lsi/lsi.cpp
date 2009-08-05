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

int print_matrix(DimInfo *dims, ifstream *datafile, ofstream *outfile);

int main(int argc, char **argv)
{
  ifstream datafile("../../../input/data.txt");

  if (argc < 2) {
    fprintf(stderr, "Usage: %s (-p|-r)\n", argv[0]);
    return 2;
  }

  DimInfo * dims = generate_dims(&datafile);

  cout << dims->user_to.size() << endl;
  cout << dims->repos_to.size() << endl;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s (-p|-r)\n", argv[0]);
    return 2;
  }
  if (strncmp(argv[1], "-p", 20) == 0) {
    ifstream datafile2("../../../input/data.txt");
    ofstream outfile("matrix.dat", ofstream::binary);
    print_matrix(dims, &datafile2, &outfile);
  }

  destroy_dims(dims);

  return 0;
}

int print_matrix(DimInfo *dims, ifstream *datafile, ofstream *outfile) {
  dimmap user_to = dims->user_to;
  dimmap repos_to = dims->repos_to;
  intmap user_counts = dims->user_counts;

  dimmap matrixvals;

  int rows = user_to.size();
  int cols = repos_to.size();

  datafile->seekg(0);

  char *label = new char[50];

  for (string line; getline(*datafile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string user = *i++;
    string repo = *i;

    int userid = user_to[user];
    int repoid = repos_to[repo];

    sprintf(label, "%u,%u", userid, repoid);

    if (user_counts.count(userid) && user_counts[userid])
      matrixvals[label] = user_counts[userid];
  }

  
  for (int i=0; i < rows; i++) {
    int size = matrixvals.size();
    outfile->write((char *)&rows, 4);
    outfile->write((char *)&cols, 4);
    outfile->write((char *)&size, 4);

    stringstream oss;

    int num_nonzero = 0;
    int num;
    float x;
    for (int j=0; j < cols; j++) {
      sprintf(label, "%u,%u", i, j);
      if (matrixvals.count(label)) {
        num = matrixvals[label];
        if (num) {
          x = 1.0 / float(num);
          oss.write((char *)&j, 4);
          oss.write((char *)&x, 4);
          num_nonzero ++;
        }
      }
    }
    outfile->write((char *)&num_nonzero, 4);
    (*outfile) << oss.rdbuf();
  }

  printf("SIZE: %u\n", matrixvals.size());

  delete label;
  return 0;
}
