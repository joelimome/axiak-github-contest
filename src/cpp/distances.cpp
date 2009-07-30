#include <boost/tuple/tuple.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <map>

using namespace boost;

template <typename DistanceMap>
class d_recorder : public default_bfs_visitor
{
public:
  d_recorder(DistanceMap dist) : d(dist) { }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) const
  {
    typename graph_traits<Graph>::vertex_descriptor
      u = source(e, g);
    if (d[u] > 2) {
      return;
    }
    typename graph_traits<Graph>::vertex_descriptor
      v = target(e, g);
    d[v] = d[u] + 1;
  }
private:
    DistanceMap d;
};

// Convenience function
template < typename DistanceMap >
d_recorder<DistanceMap>
record_d(DistanceMap d)
{
  return d_recorder < DistanceMap > (d);
}


int
main()
{
  std::ifstream datafile("./edge.dat");
  if (!datafile) {
    std::cerr << "No ./edge.dat file" << std::endl;
    return EXIT_FAILURE;
  }

  typedef adjacency_list < vecS, vecS, undirectedS, property < vertex_name_t,
    std::string >, property < edge_name_t, std::string > > Graph;
  Graph g;

  typedef property_map < Graph, vertex_name_t >::type user_name_map_t;
  user_name_map_t user_name = get(vertex_name, g);
  typedef property_map < Graph, edge_name_t >::type movie_name_map_t;
  movie_name_map_t connecting_movie = get(edge_name, g);

  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::map < std::string, Vertex > NameVertexMap;
  NameVertexMap users;

  for (std::string line; std::getline(datafile, line);) {
    char_delimiters_separator < char >sep(false, "", ":");
    tokenizer <> line_toks(line, sep);
    tokenizer <>::iterator i = line_toks.begin();
    std::string users_name = *i++;
    NameVertexMap::iterator pos;
    bool inserted;
    Vertex u, v;
    tie(pos, inserted) = users.insert(std::make_pair(users_name, Vertex()));
    if (inserted) {
      u = add_vertex(g);
      user_name[u] = users_name;
      pos->second = u;
    } else
      u = pos->second;

    tie(pos, inserted) = users.insert(std::make_pair(*i, Vertex()));
    if (inserted) {
      v = add_vertex(g);
      user_name[v] = *i;
      pos->second = v;
    } else
      v = pos->second;

    graph_traits < Graph >::edge_descriptor e;
    tie(e, inserted) = add_edge(u, v, g);
  }


  for (NameVertexMap::const_iterator it = users.begin(); it != users.end(); ++it) {
    /* cout << "Who(key = first): " << it->first;
    cout << " Score(value = second): " << it->second << '\n';
    */
    std::vector < int >d(num_vertices(g));
    Vertex src = it->second;
    d[src] = 0;

    breadth_first_search(g, src,
                         visitor(record_d(&d[0])));

    std::cout << it->first;
    graph_traits < Graph >::vertex_iterator i, end;
    for (tie(i, end) = vertices(g); i != end; ++i) {
      if (d[*i] && d[*i] < 3) {
        std::cout << "," << user_name[*i] << ":"
                  << d[*i];
      }
    }
    std::cout << std::endl;
  }
  return 0;
}

