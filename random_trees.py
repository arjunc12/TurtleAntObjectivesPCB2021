import networkx as nx
from itertools import combinations
from map_network import *
from random import choice, shuffle, uniform, seed, randint
from collections import defaultdict
import argparse
from repeatability_utils import *
from itertools import product, combinations
from sys import argv, maxint

def exists_path(G, u, v):
    return nx.has_path(G, u, v) or nx.has_path(G, v, u)

def are_connected(G, u, v):
    for w in G.nodes():
        if exists_path(G, u, w) and exists_path(G, v, w):
            return True
    return False

def random_steiner(G, terminals=None, direction='out'):
    if terminals == None:
        terminals = G.graph['terminals']

    terminals = terminals[:]
    shuffle(terminals)

    DG = directed_subgraph(G, direction=direction, terminals=terminals)
    for u, v in G.edges():
        if G[u][v]['direction'] == direction:
            assert DG.has_edge(u, v)

    for u, v in DG.edges():
        if u not in terminals and v not in terminals:
            assert G[u][v]['direction'] == direction

    T = nx.DiGraph()
    root = terminals.pop(0)
    T.add_node(root)


    def path_prefix(T, path):
        prefix = []
        for u in path:
            prefix.append(u)
            if T.has_node(u):
                break

        return prefix

    def path_suffix(T, path):
        return path_prefix(T, reversed(path))

    while len(terminals) > 0:
        candidate_paths = None
        next_terminal = None
        for t in terminals:
            candidate_paths = set()
            for u in T.nodes():
                paths1 = nx.all_simple_paths(DG, t, u)
                paths2 = nx.all_simple_paths(DG, u, t)
                for path in paths1:
                    assert t == path[0]
                    candidate_paths.add(tuple(path_prefix(T, path)))

                for path in paths2:
                    assert t in path[-1]
                    candidate_paths.add(tuple(path_suffix(T, path)))

            if len(candidate_paths) > 0:
                next_terminal = t
                break

        assert candidate_paths != None
        assert next_terminal != None
        terminals.remove(next_terminal)

        path = choice(list(candidate_paths))
        for i in xrange(len(path) - 1):
            curr, succ = path[i], path[i + 1]
            T.add_edge(curr, succ)
            T.add_edge(succ, curr)

            for key, val in G[curr][succ].iteritems():
                T[curr][succ][key] = val
            for key, val in G[succ][curr].iteritems():
                T[succ][curr][key] = val

    T.graph['terminals used'] = terminals

    assert nx.is_strongly_connected(T)
    for terminal in terminals:
        assert T.has_node(terminal)
    return T


def random_steiner_lerw(G, terminals=None, euclidean=False):
    if terminals == None:
        terminals = G.graph['terminals'][:]
    shuffle(terminals)

    H = nx.DiGraph()

    root = terminals[0]

    in_tree = defaultdict(bool)
    in_tree[root] = True
    successor = {}
    successor[root] = None

    edges_removed = set()

    for u in terminals[1:]:
        curr = u
        curr_loop = []
        while not in_tree[curr]:
            candidates = list(G.neighbors(curr))
            succ = choice(candidates)
            successor[curr] = succ
            edges_removed.add((curr, succ))
            edges_removed.add((succ, curr))

            curr = successor[curr]

        curr = u
        while not in_tree[curr]:
            in_tree[curr] = True
            succ = successor[curr]

            for n in [curr, succ]:
                if not H.has_node(n):
                    H.add_node(n)
                    for key, val in G.node[n].iteritems():
                        H.node[n][key] = val

            H.add_edge(curr, succ)
            H.add_edge(succ, curr)

            edges_removed.remove((curr, succ))
            edges_removed.remove((succ, curr))

            for key, val in G[curr][succ].iteritems():
                H[curr][succ][key] = val
            for key, val in G[succ][curr].iteritems():
                H[succ][curr][key] = val

            curr = succ

    H.graph['terminals used'] = terminals[:]
    H.graph['edges removed'] = edges_removed

    return H

def directed_random_steiner_subgraph_lerw(G, terminals=None, euclidean=False, direction=None):
    '''
    modifed loop-erased random walk method to construct a steiner tree. picks directed paths randomly
    '''
    if direction != None:
        G = directed_subgraph(G, direction=direction)

    nodes_used = set()
    nodes_used.add(terminals[0])

    for terminal in terminals[1:]:
        curr = terminal
        curr_path = set()
        prev = None
        while curr not in nodes_used:
            curr_path.add(curr)
            predecessors = list(G.predecessors(curr))
            if prev in predecessors:
                predecessors.remove(prev)

            successors = list(G.successors(curr))
            if prev in successors:
                successors.remove(prev)

            print curr
            print list(G.neighbors(curr))
            succ = choice(list(G.neighbors(curr)))

            prev = curr
            curr = succ

        nodes_used.update(curr_path)

    T = G.subgraph(nodes_used)
    assert nx.is_strongly_connected(T)
    for terminal in terminals:
        assert T.has_node(terminal)
    return T.to_directed()

def random_steiner_subgraph_lerw(G, terminals=None, euclidean=False):
    if terminals == None:
        terminals = G.graph['terminals']

    nodes_used = set()
    edges_used = set()
    nodes_used.add(terminals[0])

    cycle_nodes = set()

    for terminal in terminals[1:]:
        if terminal in nodes_used:
            continue

        curr = terminal
        curr_path = set()
        walk = []
        while curr not in nodes_used:
            curr_path.add(curr)
            walk.append(curr)

            candidates = []
            for u in G.neighbors(curr):
                if u not in terminals or u in nodes_used:
                    candidates.append(u)

            successor = choice(candidates)
            curr = successor
        nodes_used.update(curr_path)

        walk.append(curr)
        for i in xrange(len(walk) - 1):
            edges_used.add((walk[i], walk[i + 1]))
            edges_used.add((walk[i + 1], walk[i]))

        H = nx.Graph()
        H.add_path(walk)
        assert walk[0] == terminal
        assert H.has_node(terminal)
        H.remove_node(terminal)
        for component in nx.connected_components(H):
            C = nx.Graph()
            C = H.subgraph(component)
            C = C.copy()

            remove_dead_ends(C)
            apoints = nx.articulation_points(C)
            cnodes = set(C.nodes())
            for ap in apoints:
                assert C.degree(ap) > 1
                if C.degree(ap) == 2:
                    cnodes.remove(ap)

            cycle_nodes.update(cnodes)

    #T = G.subgraph(nodes_used)
    T = G.edge_subgraph(edges_used)
    T.graph['cycle nodes'] = cycle_nodes

    assert nx.is_connected(T)
    for terminal in terminals:
        assert T.has_node(terminal)
    return T.to_directed()

def random_steiner_independent(G):
    trees = []
    for day in G.graph['days']:
        G2 = G.copy().to_undirected()
        for u, v, rupture_day in G.graph['edge ruptures']:
            if rupture_day == day:
                G2.remove_edge(u, v)

        for u, rupture_day in G.graph['node ruptures']:
            if rupture_day == day:
                G2.remove_node(u)

        D = day_subnetwork(G, day)
        terminals_used = D.graph['terminals used'][:]
        T = random_steiner_subgraph_lerw(G2, terminals=terminals_used)
        trees.append(T)

    return trees

def random_steiner_progression(G):
    trees = []
    prev_terminals = []
    prev_tree = nx.Graph()
    prev_tree.graph['terminals'] = []
    endpoints = set()

    for day in G.graph['days']:
        G2 = G.copy().to_undirected()

        for u, v, rupture_day in G.graph['edge ruptures']:
            if rupture_day == day:
                G2.remove_edge(u, v)

        for u, rupture_day in G.graph['node ruptures']:
            if rupture_day == day:
                G2.remove_node(u)


        D = day_subnetwork(G, day)
        T = prev_tree.copy()
        terminals_used = D.graph['terminals used']

        for terminal in prev_terminals:
            if terminal not in terminals_used:
                remove_terminal(G2, T, terminal, endpoints)

        shuffle(terminals_used)
        for terminal in terminals_used:
            if terminal not in prev_terminals:
                add_terminal(G2, T, terminal, endpoints)

        components = list(nx.connected_components(T))
        assert len(components) == 1
        assert T.number_of_nodes() >= 3

        trees.append(T.to_directed())
        prev_terminals = terminals_used
        prev_tree = T

    return trees

def random_walk_to_tree(G, T, u):
    if T.number_of_nodes() == 0:
       return [u]

    assert G.has_node(u)
    assert not T.has_node(u)

    curr = u
    walk = []
    while not T.has_node(curr):
        walk.append(curr)
        candidates =[]
        for u in G.neighbors(curr):
            if u not in G.graph['terminals'] or T.has_node(u):
                candidates.append(u)
        curr = choice(candidates)

    walk.append(curr)
    return walk

def add_terminal(G, T, terminal, endpoints):
    assert terminal not in T.graph['terminals']
    T.graph['terminals'].append(terminal)

    walk = random_walk_to_tree(G, T, terminal)
    T.add_path(walk)
    T.node[terminal]['walk'] = walk

    endpoint = walk[-1]
    assert T.has_node(endpoint)
    if len(walk) > 1 and endpoint:
        endpoints.add(endpoint)
        T.node[terminal]['endpoint'] = endpoint

def remove_terminal(G, T, terminal, endpoints):
    assert terminal in T.graph['terminals']
    T.graph['terminals'].remove(terminal)

    walk = T.node[terminal]['walk']
    walk_to_remove = None

    if len(walk) > 1:
        walk_to_remove = walk
        assert walk_to_remove[0] not in endpoints
    else:
        assert walk == [terminal]
        assert terminal in endpoints
        for term in T.graph['terminals']:
            if T.node[term]['endpoint'] == terminal:
                walk_to_remove = T.node[term]['walk']
                break
        endpoints.remove(terminal)
        assert walk != walk_to_remove

    remove_in_reverse = walk != walk_to_remove

    appearances = defaultdict(list)
    endpoint_appearances = []
    for i, u in enumerate(walk_to_remove):
        appearances[u].append(i)
        if u in endpoints:
            endpoint_appearances.append(i)

    assert len(endpoint_appearances) > 0 or len(T.graph['terminals']) == 1

    nodes_to_remove = []
    if len(endpoint_appearances) == 0:
        assert len(T.graph['terminals']) == 1
        term = T.graph['terminals'][0]
        assert walk_to_remove == T.node[term]['walk']
        assert walk_to_remove[0] == term
        nodes_to_remove = walk_to_remove[1:]
    else:
        for u, appearances in appearances.iteritems():
            if remove_in_reverse:
                if min(appearances) > endpoint_appearances[-1]:
                    nodes_to_remove.append(u)
            else:
                if max(appearances) < endpoint_appearances[0]:
                    nodes_to_remove.append(u)

    T.remove_nodes_from(nodes_to_remove)

def remove_ruptured_edges(G, T, ruptures):
    assert nx.is_connected(T)

    for u, v in ruptures:
        assert not G.has_edge(u, v)
        if T.has_edge(u, v):
            T.remove_edge(u, v)

    for u, v in ruptures:
        if nx.has_path(T, u, v):
            continue
        components = nx.connected_components(T)
        end_component = None
        for component in components:
            if v in component:
                assert u not in component
                end_component = component
                break

        T2 = T.subgraph(end_component)
        walk = random_walk_to_tree(G, T2, u)
        assert walk[0] == u
        assert T2.has_node(walk[-1])

        T.add_path(walk)

    assert nx.is_connected(T)

def remove_ruptured_nodes(G, T, ruptures):
    assert nx.is_connected(T)

    rupture_neighbors = set()
    for u in ruptures:
        assert not G.has_node(u)
        if T.has_node(u):
            rupture_neighbors.update(list(T.neighbors(u)))
            T.remove_node(u)

    for rupture_neighbor in rupture_neighbors:
        components = list(nx.connected_components(T))
        if len(components) == 1:
            break

        other_components = []
        for component in components:
            if rupture_neighbor not in component:
                other_components += component

        T2 = T.subgraph(other_components)
        walk = random_walk_to_tree(G, T2, rupture_neighbor)
        assert walk[0] == rupture_neighbor
        assert T2.has_node(walk[-1])

        T.add_path(walk)

    assert nx.is_connected(T)


def main():
    for colony in COLONIES:
        G, LG = read_network_paper(colony)
        trees = random_steiner_progression(G)
        trees2 = random_steiner_independent(G)

if __name__ == '__main__':
    main()
