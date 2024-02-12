// Find cliques
// Meet-in-the-Middle solution O(n/2)
#include <assert.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>

using std::cin;
using std::cout;
using std::vector;

constexpr int kMaxVerticesCount = 62;

namespace BitManipulation {
constexpr bool IsSubmask(int64_t mask, int64_t submask) {
  return (submask & mask) == submask;
}
constexpr int64_t GetMaskByBitIndex(int vertex) { return 1llu << vertex; }
constexpr int64_t BitExtract(int64_t mask, int i) {
  return mask ^ GetMaskByBitIndex(i);
}
constexpr int64_t UnitCount(int64_t mask) {
  mask = (mask & 0x55555555u) + ((mask >> 1) & 0x55555555u);
  mask = (mask & 0x33333333u) + ((mask >> 2) & 0x33333333u);
  mask = (mask & 0x0f0f0f0fu) + ((mask >> 4) & 0x0f0f0f0fu);
  return mask;
}
}  // namespace BitManipulation

class MaxClique {
 public:
  // O(2^(n/2))
  using Graph64Vertices = int64_t;
  using Graph32Vertices = int32_t;

  void FillLeftAndRightGraphs(const vector<vector<bool>>& general_graph,
                              vector<vector<bool>>& left_graph,
                              vector<vector<bool>>& right_graph) {
    for (size_t i = 0; i < kN; ++i) {
      for (size_t j = 0; j < kN; ++j) {
        if (i < kN / 2 && j < kN / 2) {
          left_graph[i][j] = general_graph[i][j];
        }
        if (i >= kN / 2 && j >= kN / 2) {
          right_graph[GIndToRInd(i)][GIndToRInd(j)] = general_graph[i][j];
        }
      }
    }
  }

  MaxClique(const vector<vector<bool>>& general_graph)
      : kN(general_graph.size()) {
    assert(kN <= kMaxVerticesCount);
    vector<vector<bool>> left_graph(kN / 2, vector<bool>(kN / 2, false));
    vector<vector<bool>> right_graph(kN - kN / 2,
                                     vector<bool>(kN - kN / 2, false));
    FillLeftAndRightGraphs(general_graph, left_graph, right_graph);
    vector<bool> is_clique_in_leftgraph = FindCliques(left_graph);
    vector<Graph32Vertices> max_clique_by_right_subgraph =
        FindMaxCliquesByMatrixGraph(right_graph);
    vector<Graph32Vertices> neighbors_left_graph_in_right_graph =
        FindNeigboursLeftGraphInRightGraph(general_graph, left_graph,
                                           right_graph);
    max_clique_ = 0;
    for (Graph32Vertices l_subgraph = 0;
         l_subgraph < GetNextGraphAfterCliqueOnNVertices(left_graph.size());
         ++l_subgraph) {
      if (is_clique_in_leftgraph[l_subgraph]) {
        Graph64Vertices clique =
            l_subgraph |
            (static_cast<Graph64Vertices>(
                 max_clique_by_right_subgraph
                     [neighbors_left_graph_in_right_graph[l_subgraph]])
             << left_graph.size());
        if (VerticesCount(max_clique_) < VerticesCount(clique)) {
          max_clique_ = clique;
        }
      }
    }
  }

  // Output max clique
  friend std::ostream& operator<<(std::ostream& stream,
                                  const MaxClique& cliques) {
    for (size_t older = 0; older < cliques.kN; ++older) {
      if (IsSubgraph(cliques.max_clique_, GetGraphWithOneVertex(older))) {
        stream << older << ' ';
      }
    }

    return stream;
  }

 private:
  static int VerticesCount(Graph64Vertices mask) {
    return BitManipulation::UnitCount(mask);
  }

  static Graph64Vertices RemoveVertex(Graph64Vertices mask, int vertex) {
    return BitManipulation::BitExtract(mask, vertex);
  }

  static bool IsSubgraph(Graph64Vertices mask, Graph64Vertices submask) {
    return BitManipulation::IsSubmask(mask, submask);
  }

  static Graph64Vertices GetGraphWithOneVertex(int vertex) {
    return BitManipulation::GetMaskByBitIndex(vertex);
  }
  static Graph64Vertices GetNextGraphAfterCliqueOnNVertices(
      int vertices_count) {
    return BitManipulation::GetMaskByBitIndex(vertices_count);
  }

  static Graph64Vertices UnionGraph(Graph64Vertices mask1,
                                    Graph64Vertices mask2) {
    return mask1 | mask2;
  }

  static Graph64Vertices IntersectGraph(Graph64Vertices mask1,
                                        Graph64Vertices mask2) {
    return mask1 & mask2;
  }

  int RIndToGInd(int r_ind) { return kN / 2 + r_ind; }

  int GIndToRInd(int g_ind) { return g_ind - kN / 2; }

  vector<Graph32Vertices> FindNeigboursOfEachVertex(
      const vector<vector<bool>>& matrix_graph) {
    assert(matrix_graph.size() <= kMaxVerticesCount / 2);
    vector<Graph32Vertices> neighbor(matrix_graph.size(), 0);
    for (size_t i = 0; i < matrix_graph.size(); ++i) {
      for (size_t j = 0; j < matrix_graph.size(); ++j) {
        if (matrix_graph[i][j]) {
          neighbor[i] = UnionGraph(neighbor[i], GetGraphWithOneVertex(j));
        }
      }
    }
    return neighbor;
  }
  // returns right_graph masks O(2^(n/2))
  vector<Graph32Vertices> FindNeigboursLeftGraphInRightGraph(
      const vector<vector<bool>>& general_graph,
      const vector<vector<bool>>& left_graph,
      const vector<vector<bool>>& right_graph) {
    vector<Graph32Vertices> dp(
        GetNextGraphAfterCliqueOnNVertices(left_graph.size()), 0);

    dp[0] = GetNextGraphAfterCliqueOnNVertices(right_graph.size()) - 1;
    for (size_t left_vertex = 0; left_vertex < left_graph.size();
         ++left_vertex) {
      Graph32Vertices r_mask = 0;
      for (size_t right_vertex = 0; right_vertex < right_graph.size();
           ++right_vertex) {
        if (general_graph[left_vertex][RIndToGInd(right_vertex)]) {
          r_mask = UnionGraph(r_mask, GetGraphWithOneVertex(right_vertex));
        }
      }
      dp[GetGraphWithOneVertex(left_vertex)] = r_mask;
    }

    Graph32Vertices older = -1;
    for (size_t mask = 1; mask < dp.size(); ++mask) {
      if ((mask & (mask - 1)) == 0) {
        ++older;
      }
      dp[mask] =
          dp[RemoveVertex(mask, older)] & dp[GetGraphWithOneVertex(older)];
    }

    return dp;
  }
  // O(2^n)
  vector<bool> FindCliques(vector<vector<bool>>& matrix_graph) {
    assert(matrix_graph.size() <= kMaxVerticesCount / 2);
    vector<bool> dp(GetNextGraphAfterCliqueOnNVertices(matrix_graph.size()),
                    false);
    dp[0] = true;

    int older = -1;
    vector<Graph32Vertices> neighbor = FindNeigboursOfEachVertex(matrix_graph);
    for (size_t mask = 1; mask < dp.size(); ++mask) {
      if ((mask & (mask - 1)) == 0) {
        ++older;
      }
      dp[mask] = dp[RemoveVertex(mask, older)] &&
                 IsSubgraph(neighbor[older], RemoveVertex(mask, older));
    }

    return dp;
  }
  // O(2^n)
  vector<Graph32Vertices> FindMaxCliquesByMatrixGraph(
      vector<vector<bool>>& matrix_graph) {
    assert(matrix_graph.size() <= kMaxVerticesCount / 2);
    vector<Graph32Vertices> neighbors = FindNeigboursOfEachVertex(matrix_graph);
    vector<Graph32Vertices> dp(
        GetNextGraphAfterCliqueOnNVertices(matrix_graph.size()));
    dp[0] = 0;
    for (size_t older = 0; older < matrix_graph.size(); ++older) {
      dp[GetGraphWithOneVertex(older)] = GetGraphWithOneVertex(older);
    }

    int older = -1;
    for (size_t mask = 1; mask < dp.size(); ++mask) {
      if ((mask & (mask - 1)) == 0) {
        ++older;
      }
      Graph32Vertices clique_without_older = dp[RemoveVertex(mask, older)];
      Graph32Vertices clique_with_older =
          GetGraphWithOneVertex(older) | dp[neighbors[older] & mask];
      dp[mask] = clique_without_older;
      if (VerticesCount(clique_without_older) <
              VerticesCount(clique_with_older) ||
          (VerticesCount(clique_without_older) ==
               VerticesCount(clique_with_older) &&
           clique_without_older < clique_with_older)) {
        dp[mask] = clique_with_older;
      }
    }
    return dp;
  }

 private:
  const size_t kN;
  Graph64Vertices max_clique_;
};

int main() {
  int64_t n;
  cin >> n;
  vector<vector<bool>> graph(n, vector<bool>(n, false));
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      bool is_edge;
      cin >> is_edge;
      graph[i][j] = is_edge;
    }
  }
  cout << MaxClique(graph);
}