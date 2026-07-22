import math
import networkx as nx
import random
import numpy as np
import ndlib.models.ModelConfig as mc
import ndlib.models.epidemics as ep
from collections import deque


class Edge:
    def __init__(self, to, rev, cap):
        self.to = to
        self.rev = rev
        self.cap = cap


class Dinic:
    def __init__(self, n):
        self.size = n
        self.graph = [[] for _ in range(n)]

    def add_edge(self, fr, to, cap):
        forward = Edge(to, len(self.graph[to]), cap)
        backward = Edge(fr, len(self.graph[fr]), 0)
        self.graph[fr].append(forward)
        self.graph[to].append(backward)

    def bfs_level(self, s, t, level):
        q = deque()
        level[:] = [-1] * self.size
        level[s] = 0
        q.append(s)
        while q:
            v = q.popleft()
            for edge in self.graph[v]:
                if edge.cap > 0 and level[edge.to] == -1:
                    level[edge.to] = level[v] + 1
                    q.append(edge.to)
                    if edge.to == t:
                        return

    def dfs_flow(self, v, t, upTo, iter_, level):
        if v == t:
            return upTo
        for i in range(iter_[v], len(self.graph[v])):
            edge = self.graph[v][i]
            if edge.cap > 0 and level[v] < level[edge.to]:
                d = self.dfs_flow(edge.to, t, min(upTo, edge.cap), iter_, level)
                if d > 0:
                    edge.cap -= d
                    self.graph[edge.to][edge.rev].cap += d
                    return d
            iter_[v] += 1
        return 0

    def max_flow(self, s, t):
        flow = 0
        level = [-1] * self.size
        while True:
            self.bfs_level(s, t, level)
            if level[t] == -1:
                return flow
            iter_ = [0] * self.size
            while True:
                f = self.dfs_flow(s, t, float('inf'), iter_, level)
                if f == 0:
                    break
                flow += f
        return flow


def PCS_Loc(G, obs_seq):
    E = list(obs_seq.keys())
    sorted_obs = sorted(obs_seq.items(), key=lambda x: x[1])
    kb, _ = sorted_obs[0]
    kr, _ = sorted_obs[1]

    L = nx.shortest_path_length(G, kb, kr)
    R = L / 2

    C = set()
    all_nodes = list(G.nodes())
    dist_kb = nx.single_source_shortest_path_length(G, kb)
    dist_kr = nx.single_source_shortest_path_length(G, kr)
    for v in all_nodes:
        d1 = dist_kb[v]
        d2 = dist_kr[v]
        if abs(d1 - d2) <= L:
            C.add(v)

    Vc = set()
    for c in C:
        dist_c = nx.single_source_shortest_path_length(G, c)
        for node, d in dist_c.items():
            if d <= R:
                Vc.add(node)
    Vc = list(Vc)
    if len(Vc) == 0:
        return -1

    geo_score = {}
    for v in Vc:
        d1 = dist_kb[v]
        d2 = dist_kr[v]
        geo_score[v] = (d1 + d2) / 2
    sorted_Vc = sorted(geo_score.items(), key=lambda x: x[1])
    M_min = 5
    beta = 0.05
    M = max(M_min, math.ceil(beta * len(Vc)))
    VM = [node for node, score in sorted_Vc[:M]]

    node2id = {node: i for i, node in enumerate(G.nodes())}
    total_node_num = len(G.nodes())
    super_sink_id = total_node_num
    best_rho = -1
    pred_source = -1

    for cand in VM:
        dinic = Dinic(total_node_num + 1)
        for u, v in G.edges():
            dinic.add_edge(node2id[u], node2id[v], 1)
            dinic.add_edge(node2id[v], node2id[u], 1)
        for e in E:
            dinic.add_edge(node2id[e], super_sink_id, 1)
        rho = dinic.max_flow(node2id[cand], super_sink_id)
        if rho > best_rho:
            best_rho = rho
            pred_source = cand
    return pred_source


if __name__ == "__main__":
    file_path = r'soc-karate.csv'
    with open(file_path, 'r') as f:
        next(f)
        G = nx.parse_edgelist(f, delimiter=',', nodetype=int)

    total_error = 0
    correct_count = 0
    runs = 100
    obs_ratio = 0.2

    for run in range(1, runs + 1):
        source_node = random.choice(list(G.nodes()))
        model = ep.IndependentCascadesModel(G)
        config = mc.Configuration()
        config.add_model_initial_configuration("Infected", [source_node])
        p = 0.1
        for e in G.edges():
            config.add_edge_configuration("threshold", e, 1 - p)
        model.set_initial_status(config)

        infect_time = {}
        time_step = 0
        infect_time[source_node] = time_step
        all_infected = {source_node}
        max_iter = nx.diameter(G) + 5
        for _ in range(max_iter):
            iter_res = model.iteration()
            time_step += 1
            for node, status in iter_res['status'].items():
                if status == 1 and node not in all_infected:
                    infect_time[node] = time_step
                    all_infected.add(node)

        infect_list = list(all_infected)
        obs_num = max(2, int(len(infect_list) * obs_ratio))
        observers = random.sample(infect_list, obs_num)
        obs_seq = {e: infect_time[e] for e in observers}

        predicted_source = PCS_Loc(G, obs_seq)
        if predicted_source == -1:
            print(f"Run {run}: Failed to predict source, skipped.")
            continue

        dist_err = nx.shortest_path_length(G, source_node, predicted_source)
        total_error += dist_err
        if dist_err == 0:
            correct_count += 1

        print(f"Run {run}: Actual={source_node}, Predicted={predicted_source}, Distance Error={dist_err}")

    avg_error = total_error / runs
    accuracy = correct_count / runs

    print(f"\n===== PCS-Loc Results (100 runs) =====")
    print(f"Average Distance Error: {avg_error:.3f}")
    print(f"Accuracy (error=0): {accuracy * 100:.2f}%")