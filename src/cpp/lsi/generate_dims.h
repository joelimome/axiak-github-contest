#ifndef GENERATE_DIMS_H
#include <map>

using namespace std;

typedef map<string, int> dimmap;
typedef map<int, string> dimmaprev;
typedef map<int, int> intmap;

struct DimInfo {
        dimmap user_to;
        dimmap repos_to;
        dimmaprev user_from;
        dimmaprev repos_from;
        intmap user_counts;
};

DimInfo * generate_dims(ifstream *datafile);

#endif

