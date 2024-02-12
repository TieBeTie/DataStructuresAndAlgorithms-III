#include "../libs.hpp"

unordered_map<int, unordered_set<int>> g;
vector<int> type;
vector<bool> visited;

int D = -1;

bool IsD(int v) { return type[v] == D; }

bool IsP(int v) { return type[v] % 2 == 0; }

bool IsN(int v) { return type[v] % 2 == 1; }

void DFS(int v) {
  if (visited[v]) {
    return;
  }
  visited[v] = true;

  int n_to_min_p = numeric_limits<int>::max();
  int p_to_max_n = 0;
  bool is_to_tie_exist = false;

  for (int u : g[v]) {
    DFS(u);
    if (IsD(u)) {
      is_to_tie_exist = true;
    } else if (IsP(u)) {
      n_to_min_p = min(n_to_min_p, type[u] + 1);
    } else {
      p_to_max_n = max(p_to_max_n, type[u] + 1);
    }
  }

  if (n_to_min_p != numeric_limits<int>::max()) {
    type[v] = n_to_min_p;
  } else if (is_to_tie_exist) {
    type[v] = D;
  } else {
    type[v] = p_to_max_n;
  }
}

void RetrospectiveAnalysisDFS(int root) {
  type.assign(g.size(), D);
  visited.assign(g.size(), false);
  DFS(root);
}