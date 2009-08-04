#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <vector>
#include <boost/tokenizer.hpp>
#include <boost/config.hpp>

#include <armadillo>

#include "generate_dims.h"

//fmat * construct_doc_matrix (DimInfo * dims);

using namespace std;
using namespace arma;

arma::fmat * construct_doc_matrix (DimInfo * dims, ifstream * datafile);

int main()
{
  ifstream datafile("../../../input/data.txt");

  DimInfo * dims = generate_dims(&datafile);

  cout << dims->user_to.size() << endl;
  cout << dims->repos_to.size() << endl;

  arma::fmat *A = construct_doc_matrix(dims, &datafile);

  destroy_dims(dims);

  return 0;
}

fmat * construct_doc_matrix (DimInfo * dims, ifstream * datafile) {

  dimmap user_to = dims->user_to;
  dimmap repos_to = dims->repos_to;
  intmap user_counts = dims->user_counts;

  arma::fmat A(user_to.size(), repos_to.size());

  A.zeros();

  datafile->seekg(0);

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

    A(userid, repoid) = 1.0 / user_counts[userid];
  }

  return &A;
}

