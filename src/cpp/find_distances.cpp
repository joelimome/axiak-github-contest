//=======================================================================
// Copyright 2001 Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <map>

using namespace boost;

template <typename DistanceMap>
class hitler_number_recorder : public default_bfs_visitor
{
public:
  hitler_number_recorder(DistanceMap dist) : d(dist) { }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) const
  {
    typename graph_traits<Graph>::vertex_descriptor
      u = source(e, g), v = target(e, g);
      d[v] = d[u] + 1;
  }
private:
    DistanceMap d;
};

// Convenience function
template < typename DistanceMap >
hitler_number_recorder<DistanceMap>
record_hitler_number(DistanceMap d)
{
  return hitler_number_recorder < DistanceMap > (d);
}


int
main(int argc, char **argv)
{
  typedef adjacency_list < vecS, vecS, undirectedS, property < vertex_name_t,
    std::string >, property < edge_name_t, std::string > > Graph;
  Graph g;

  typedef property_map < Graph, vertex_name_t >::type article_map_t;
  article_map_t article = get(vertex_name, g);

  typedef graph_traits < Graph >::vertex_descriptor Vertex;
  typedef std::map < std::string, Vertex > NameVertexMap;
  NameVertexMap articles;

  for (std::string line; std::getline(std::cin, line);) {
    char_delimiters_separator < char >sep(false, "", ";");
    tokenizer <> line_toks(line, sep);
    tokenizer <>::iterator i = line_toks.begin();
    std::string article_name = *i++;
    NameVertexMap::iterator pos;
    bool inserted;
    Vertex u, v;
    tie(pos, inserted) = articles.insert(std::make_pair(article_name, Vertex()));
    if (inserted) {
      u = add_vertex(g);
      article[u] = article_name;
      pos->second = u;
    } else
      u = pos->second;

    tie(pos, inserted) = articles.insert(std::make_pair(*i, Vertex()));
    if (inserted) {
      v = add_vertex(g);
      article[v] = *i;
      pos->second = v;
    } else
      v = pos->second;

    graph_traits < Graph >::edge_descriptor e;
    tie(e, inserted) = add_edge(u, v, g);
  }

  std::vector < int >hitler_number(num_vertices(g));

  Vertex src = articles["4273"];
  hitler_number[src] = 0;

  breadth_first_search(g, src,
                       visitor(record_hitler_number(&hitler_number[0])));

  graph_traits < Graph >::vertex_iterator i, end;
  for (tie(i, end) = vertices(g); i != end; ++i) {
    std::cout << article[*i] << " has a Hitler number of "
      << hitler_number[*i] << std::endl;
  }

  return 0;
}

