#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <boost/tokenizer.hpp>
#include <boost/config.hpp>

#include "generate_dims.h"

using namespace std;

DimInfo * generate_dims(ifstream *datafile)
/* Compute the number of distinct users
   and return maps.
*/
{
  DimInfo *result = new DimInfo;
  
  dimmap user_to; // = new user_to;
  dimmap repos_to; // = new dimmap;
  dimmaprev user_from; //= new dimmaprev;
  dimmaprev repos_from; // = new dimmaprev;

  intmap user_counts;

  int userid = 0, repoid = 0;

  datafile->seekg(0);

  for (string line; getline(*datafile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string user = *i++;
    string repo = *i;

    if (!user_to.count(user)) {
      user_to[user] = userid;
      user_counts[userid] ++;
      user_from[userid++] = user;
    }
    else {
      user_counts[user_to[user]] ++;
    }

    if (!repos_to.count(repo)) {
      repos_to[repo] = repoid;
      repos_from[repoid++] = repo;
    }
  }

  result->repos_to = repos_to;
  result->user_to = user_to;
  result->repos_from = repos_from;
  result->user_from = user_from;
  result->user_counts = user_counts;

  datafile->seekg(0);

  return result;

}

#ifndef DAT_USERMAP
#define DAT_USERMAP "../../../dat/usermap.dat"
#endif
#ifndef DAT_REPOMAP
#define DAT_REPOMAP "../../../dat/repomap.dat"
#endif
#ifndef DAT_USERCOUNT
#define DAT_USERCOUNT "../../../dat/usercount.dat"
#endif

DimInfo * read_dims()
{
  DimInfo *result = new DimInfo;
  
  dimmap user_to;
  dimmap repos_to;
  dimmaprev user_from;
  dimmaprev repos_from;

  intmap user_counts;

  int userid = 0, repoid = 0;

  ifstream userfile(DAT_USERMAP);
  for (string line; getline(userfile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string user = *i++;
    string x = *i;
    userid = atoi(x.c_str());
    user_to[user] = userid;
    user_from[userid] = user;
  }
  userfile.close();

  ifstream repofile(DAT_REPOMAP);
  for (string line; getline(repofile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string repo = *i++;
    string x = *i;
    repoid = atoi(x.c_str());
    repos_to[repo] = repoid;
    repos_from[repoid] = repo;
  }
  repofile.close();

  ifstream usercountfile(DAT_USERCOUNT);
  for (string line; getline(usercountfile, line);) {
    if (line.rfind(":") == string::npos)
      continue;

    boost::char_delimiters_separator < char >sep(false, "", ":");
    boost::tokenizer <> line_toks(line, sep);
    boost::tokenizer <>::iterator i = line_toks.begin();

    string userstring = *i++;
    string countstring = *i;
    user_counts[atoi(userstring.c_str())] = atoi(countstring.c_str());
  }
  usercountfile.close();

  result->repos_to = repos_to;
  result->user_to = user_to;
  result->repos_from = repos_from;
  result->user_from = user_from;
  result->user_counts = user_counts;

  return result;
}


void destroy_dims(DimInfo *diminfo)
{
  /*
  delete &diminfo->user_to;
  delete &diminfo->repos_from;
  delete &diminfo->user_from;
  delete &diminfo->repos_to;
  */
  delete diminfo;
}
