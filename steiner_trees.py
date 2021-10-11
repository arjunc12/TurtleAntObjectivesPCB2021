import networkx as nx
from itertools import combinations
from map_network import *
from random import choice, shuffle, uniform, seed, randint
from collections import defaultdict
import argparse
from repeatability_utils import *
from itertools import product, combinations
from sys import argv, maxint

#seed(2)
#seed(3)
#seed(4)
#seed(33)
#seed(int(argv[1]))

REPEATABILITIES = [1, 2, 3, 4]

def path_weight(G, path, weight=None):
    wt = 0
    for i in xrange(len(path) - 1):
        if weight == None:
            wt += 1
        else:
            u, v = path[i], path[i + 1]
            wt += G[u][v][weight]
    return wt

def transitive_closure(G, weight='weight'):
    '''
    Assumes a connected graph
    '''
    assert nx.is_strongly_connected(G)
    TC = nx.DiGraph()
    paths = nx.shortest_path(G, weight=weight)
    for u, v in product(G.nodes(), G.nodes()):
        if u == v:
            TC.add_edge(u, v)
            TC[u][v][weight] = 0
        else:
            TC.add_edge(u, v)
            weight1 = path_weight(G, paths[u][v])
            TC[u][v][weight] = weight1

            TC.add_edge(v, u)
            weight2 = path_weight(G, paths[v][u])
            TC[v][u][weight] = weight2

    return TC

def steiner_partial_approx(TC, k, X, weight='weight'):
    d = float("inf")
    B = None
    for u, v in product(TC.nodes(), TC.nodes()):
        s = {}
        for a, b in X:
            s[(a, b)] = TC[a][u][weight] + TC[v][b][weight]
        s = sorted(s.iteritems(), key = lambda x: x[1])

        for p in xrange(1, k + 1):
            C = float(TC[u][v][weight])
            s2 = s[:p]
            for (a, b), wt in s2:
                C += wt

            if C / p <= d:
                d = C / p
                B = nx.DiGraph()
                B.add_edge(u, v)
                for (a, b), wt in s2:
                    B.add_edge(a, u)
                    B.add_edge(v, b)

    return B

def steiner_approx(G, weight='weight', terminals=None):
    ST = nx.DiGraph()
    if terminals == None:
        terminals = G.graph['terminals']

    if len(terminals) == 1:
        terminal = terminals[0]
        assert G.has_node(terminal)
        ST = nx.DiGraph()
        ST.add_node(terminal)
        return ST

    G2 = G.copy()
    TC = transitive_closure(G, weight=weight)

    X = list(product(terminals, terminals))
    X = filter(lambda (x, y): x != y, X)
    k = len(X)
    while k > 0:
        B = steiner_partial_approx(TC, k, X, weight=weight)
        covered = 0
        for u, v in B.edges():
            sp = nx.shortest_path(G, u, v, weight=weight)
            for i in xrange(len(sp) - 1):
                x, y = sp[i], sp[i + 1]
                ST.add_edge(x, y)
                for key, val in G[x][y].iteritems():
                    ST[x][y][key] = val

        for (a, b) in X:
            if ST.has_node(a) and ST.has_node(b) and nx.has_path(ST, a, b):
                X.remove((a, b))
                covered += 1

        assert covered <= k
        k -= covered

    ST.graph['terminals used'] = terminals[:]
    return ST

def preprocess_line_graph(G, LG, terminals):
    LG2 = LG.copy()
    for u, v in LG.nodes():
        for terminal in terminals:
            if G.has_edge(terminal, u) or G.has_edge(v, terminal):
                LG2.add_edge(terminal, (u, v))
                LG2[terminal][(u, v)]['transition index'] = 0

                LG2.add_edge((u, v), terminal)
                LG2[(u, v)][terminal]['transition index'] = 0

    return LG2

def postprocess_line_graph(G, ST):
    edgelist = []
    for x in ST.nodes():
        if len(x) == 2:
            u, v = x
            edgelist += [(u, v), (v, u)]
        else:
            assert x in G.graph['terminals']
            for u, v in ST.neighbors(x):
                if G.has_edge(x, u):
                    edgelist += [(x, u), (u, x)]
                elif G.has_edge(v, x):
                    edgelist += [(v, x), (x, v)]
    T = G.edge_subgraph(edgelist)
    assert nx.is_strongly_connected(T)
    assert T.number_of_edges() > 2
    return T

def average_steiner_approx(G, weight='weight', terminals=None):
    assert nx.is_directed(G)
    if terminals == None:
        terminals = G.graph['terminals']

    node_components = {}
    connected_components = {}
    for i, u in enumerate(G.nodes()):
        node_components[u] = i
        connected_components[i] = [u]

    G2 = G.to_undirected()
    sorted_edges = sorted(list(G2.edges()), key = lambda (x, y) : G[x][y][weight])
    i = 0
    T = nx.Graph()
    for u, v in sorted_edges:
        T.add_edge(u, v)
        wt = G[u][v][weight]
        T[u][v][weight] = wt

        comp1 = node_components[u]
        comp2 = node_components[v]
        if comp1 == comp2:
            continue

        min_comp = min(comp1, comp2)
        relabel_comp = max(comp1, comp2)
        for c in connected_components[relabel_comp]:
            connected_components[min_comp].append(c)
            node_components[c] = min_comp
        del connected_components[relabel_comp]

        min_terminal_comp = float("inf")
        max_terminal_comp = float("-inf")
        for terminal in terminals:
            min_terminal_comp = min(min_terminal_comp, node_components[terminal])
            max_terminal_comp = max(max_terminal_comp, node_components[terminal])

        if min_terminal_comp == max_terminal_comp:
            break

    for u, component in node_components.iteritems():
        if T.has_node(u) and component != node_components[terminals[0]]:
            T.remove_node(u)

    assert nx.is_connected(T)
    for terminal in terminals:
        assert T.has_node(terminal)

    total_weight = 0
    num_edges = 0
    for u, v in T.edges():
        total_weight += T[u][v][weight]
        num_edges += 1

    total_weight = float(total_weight)
    current_average = total_weight / num_edges

    sorted_edges = sorted(list(T.edges()), key = lambda (x, y) : T[x][y][weight])
    for u, v in reversed(sorted_edges):
        assert T.has_edge(u, v)
        wt = T[u][v][weight]
        if wt < current_average:
            break

        T.remove_edge(u, v)
        if nx.has_path(T, u, v):
            total_weight -= wt
            num_edges -= 1
            current_average = total_weight / num_edges
        else:
            T.add_edge(u, v)
            T[u][v][weight] = wt

        i -= 1

    assert nx.is_connected(T)

    edgelist = []
    for u, v in T.edges():
        edgelist += [(u, v), (v, u)]
    return G.edge_subgraph(edgelist)

def draw_random_trees(colony):
    G, LG = read_network_paper(colony)

    graphscale = 1
    for u, v in G.edges():
        length = G[u][v]['length']
        G[u][v]['drawing length'] = pylab.log(1 + length)
    pos = nx.kamada_kawai_layout(G, scale=graphscale)#, weight='drawing length')

    days = G.graph['days']
    trees = random_steiner_progression(G)
    assert len(days) == len(trees)

    fig = pylab.figure()

    def redraw(frame):
        print frame
        pylab.clf()
        nodelist = []
        node_color = []
        T = trees[frame]
        day = days[frame]
        D = day_subnetwork(G, day)

        for u in G.nodes():
            nodelist.append(u)

            ants = int(D.has_node(u))
            rand = 2 * int(T.has_node(u))

            status = ants + rand
            if status == 3:
                node_color.append('r')
            elif status == 2:
                node_color.append('r')
            elif status == 1:
                node_color.append('k')
            else:
                node_color.append('k')

        nx.draw(G, pos=pos, nodelist=nodelist, node_color=node_color, arrows=False,\
                with_labels=True, font_color='w', font_size=3)
        pylab.draw()

    ani = animation.FuncAnimation(fig, redraw, frames=len(days), interval=5000)
    #mywriter = animation.AVConvWriter()
    ani.save('network_changes/random-trees.mp4', dpi=300)
    pylab.close()

def opt_total_nodes(G, LG, terminals):
    return steiner_approx(G, weight='length', terminals=terminals)

def opt_total_length(G, LG, terminals):
    return steiner_approx(G, weight='length', terminals=terminals)

def opt_average_length(G, LG, terminals):
    return average_steiner_approx(G, weight='length', terminals=terminals)

def opt_total_transition_index(G, LG, terminals):
    LG2 = preprocess_line_graph(G, LG, terminals)
    ST = steiner_approx(LG2, weight='transition index', terminals=terminals)
    assert nx.is_strongly_connected(ST)
    assert ST.number_of_nodes() > 1
    return postprocess_line_graph(G, ST)

def opt_average_transition_index(G, LG, terminals):
    LG2 = preprocess_line_graph(G, LG, terminals)
    ST = average_steiner_approx(LG2, weight='transition index', terminals=terminals)
    assert nx.is_strongly_connected(ST)
    assert ST.number_of_nodes() > 1
    return postprocess_line_graph(G, ST)

def main():
    for colony in ['T500']:
        print colony
        G, LG = read_network_paper(colony)
        for day in G.graph['days']:
            print day
            D = day_subnetwork(G, day)
            terminals = D.graph['terminals used']
            LG2 = preprocess_line_graph(G, LG, terminals)
            ST = average_steiner_approx(LG2, terminals=terminals, weight='transition index')
            assert nx.is_strongly_connected(ST)
            assert ST.number_of_nodes() > 1
            T = postprocess_line_graph(G, ST)
            print T.number_of_nodes()

            ST2 = average_steiner_approx(G, terminals=terminals, weight='length')
            assert nx.is_strongly_connected(ST2)
            print ST2.number_of_nodes()

if __name__ == '__main__':
    main()
