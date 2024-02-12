#include "../libs.hpp"

unordered_map<int, unordered_set<int>> g;
vector<char> type;
vector<bool> visited;

void DFS(int v) {
  if (visited[v]) {
    return;
  }
  visited[v] = true;

  for (int u : g[v]) {
    DFS(u);
    if (type[u] == 'P') {
      type[v] = 'N';
      return;
    }
    if (type[u] == 'D') {
      type[v] = 'D';
    }
  }
}

void RetrospectiveAnalysisDFS(int root) {
  type.assign(g.size(), 'P');
  visited.assign(g.size(), false);
  DFS(root);
}