import matplotlib as mpl
mpl.use('agg')
from matplotlib import animation
import pandas as pd
import networkx as nx
from sys import argv, exit
from random import sample
import pylab
from networkx.algorithms import approximation
from collections import OrderedDict, defaultdict
from itertools import combinations, product
import os
import argparse
#from steiner_approx import steiner_approx

DEFAULT_REPEATABILITY = 1
DEFAULT_LENGTH = 5 # cm

DATA_DIR = 'mapping_network/CAustin-Chamela-09-14-2018/'

#DATA_DIR = 'mapping_network/csv'
HOP_DRAWINGS_DIR = 'hop_neighborhoods'

COLONIES = ['turtle_hill', 'tejon7', 'T189', 'T500']
#COLONY_LABELS = {'turtle_hill' : 'EC460', 'tejon7' : 'T430', 'T189' : 'T189'}
COLONY_LABELS = {'turtle_hill' : 'Turtle Hill 460', 'tejon7' : 'Tejon 446',\
                 'T189' : 'Tejon 189', 'T500' : 'Tejon 500'}

RAW_DATA_DIR = 'paper_data/raw_data'
PLOT_DATA_DIR = 'paper_data/plot_data'

COLONY_RAW_DATA_FILES = {'turtle_hill' : 'Turtle_Hill_460_network_data2018.csv',\
                         'tejon7' : 'Tejon_446_network_data2018.csv',\
                         'T189' : 'Tejon_189_network_data2018.csv',\
                         'T500' : 'Tejon_500_network_data2019.csv'}


def similar_nodes(G):
    '''
    Go through a graph and print out nodes that have different labels but the
    same set of neighbors. These nodes are possibly the same node named
    differently
    '''
    for u, v in combinations(list(G.nodes()), 2):
        n1 = sorted(G.neighbors(u))
        n2 = sorted(G.neighbors(v))
        if n1 == n2:
            print u, v, n1

def remove_dead_ends(G):
    '''
    Go through a graph and remove all dead ends nodes and edges
    '''
    done = False
    while not done:
        bad_nodes = [] # nodes with degree <= 1
        for u in G.nodes():
            if G.degree(u) <= 1:
                bad_nodes.append(u)
        '''
        if there are bad nodes, get rid of them and repeat. Otherwise the
        algorithm is complete.
        '''
        if len(bad_nodes) > 0:
            for bad_node in bad_nodes:
                assert G.degree(bad_node) <= 1
                G.remove_node(bad_node)
        else:
            done = True

    for u in G.nodes():
        assert G.degree(u) > 1

def dfs(G, start, end):
    '''
    DFS subroutine for computing cycles in a graph

    Non-deterministically traverse all paths from start to end. Return a
    generator for all such paths
    '''
    fringe = [(start, [])]
    while fringe:
        state, path = fringe.pop()
        if len(path) > 2 and state == end:
            '''
            If a sufficiently long path is found, add it to the generator
            '''
            yield path
            continue
        for next_state in G.neighbors(state):
            '''
            try all possible ways to extend the path
            '''
            if next_state in path:
                continue
            fringe.append((next_state, path+[next_state]))

def get_cycles(G, terminals):
    '''
    Find all undirected cycles in a network. Exclude all nodes designated as
    terminals.
    '''
    G2 = G.to_undirected()

    # remove dead ends and terminals to reduce running time
    for terminal in terminals:
        G2.remove_node(terminal)
    remove_dead_ends(G2)

    # iterate through nodes, find all cycles starting at each node
    cycles = [[node]+path  for node in G2.nodes() for path in dfs(G2, node, node)]
    unique_cycles = set()

    # rotate cycles to start at the minimal node, remove duplicates
    for cycle in cycles:
        assert len(cycle) > 3
        assert cycle[0] == cycle[-1]
        cycle = cycle[:-1]
        start = min(cycle)
        start_pos = cycle.index(start)
        cycle = tuple(cycle[start_pos:] + cycle[:start_pos] + [cycle[start_pos]])
        #cycle = tuple(sorted(cycle))
        if tuple(reversed(cycle)) not in unique_cycles:
            unique_cycles.add(cycle)

    return list(unique_cycles)

def sort_key(u):
    '''
    Helper function for ordering nodes based on their name

    Node names have two parts: a leading integer, and a trailing
    sequence of letter/integer pairs i.e. 24A1B2C1, 15A2, etc.
    We sort nodes by the number of trailing letter number pairs, and then use
    the leading integer as a tiebreaker.
    '''
    pair = node_to_pair(u)
    if pair == None:
        return None
    else:
        x, y = pair
        return (len(y), x)

def get_contraction_map(G, contractions):
    '''
    Determine which nodes are equivalent. Within each group of equivalent
    nodes, determine which one should remain, and which nodes should be
    contracted into the remaining node.
    '''

    # Create a connecting pairs of equivalent nodes and find connected components
    contraction_graph = nx.Graph()
    for u1, u2 in contractions:
        contraction_graph.add_edge(u1, u2)

    contraction_map = {}
    for component in nx.connected_components(contraction_graph):
        component = list(component)
        assert len(component) >= 2

        '''
        sort nodes based on our lexicographic naming scheme, then pick
        the first node that the graph contains as the node taht will remain.
        '''
        component = sorted(component, key = lambda u : sort_key(u))
        final_label = None
        for u in component:
            if G.has_node(u):
                final_label = u
                break
        assert final_label != None

        contraction_map[final_label] = component

    return contraction_map

def copy_edge_data(G, u1, v1, u2, v2):
    '''
    Copies edge attributes from edge (u1, v1) to edge (u2, v2)

    We use this when contracting nodes to preserve data from edges that are
    removed as a result of the contraction
    '''
    for key, val in G[u1][v1].iteritems():
        G[u2][v2][key] = val

def contract_nodes(G, u1, u2):
    '''

    '''
    assert u1 != u2
    assert G.has_node(u1)
    assert G.has_node(u2)

    v1 = G.number_of_nodes()

    def all_days_used(u1, v1, u2, v2):
        days1 = set(G[u1][v1]['days used'])
        days2 = set(G[u2][v2]['days used'])
        return list(days1 | days2)

    for p in G.predecessors(u2):
        if not G.has_edge(p, u1):
            G.add_edge(p, u1)
            copy_edge_data(G, p, u2, p, u1)
        else:
            G[p][u1]['days used'] = all_days_used(p, u1, p, u2)

    for s in G.successors(u2):
        if not G.has_edge(u1, s):
            G.add_edge(u1, s)
            copy_edge_data(G, u2, s, u1, s)
        else:
            G[u1][s]['days used'] = all_days_used(u1, s, u2, s)

    G.remove_node(u2)

    v2 = G.number_of_nodes()
    assert v2 < v1

    if G.has_edge(u1, u1):
        G.remove_edge(u1, u1)

    return G

def contract_graph(G, contraction_map):
    v1 = G.number_of_nodes()
    for final_label, equivalent_nodes in contraction_map.iteritems():
        assert G.has_node(final_label)
        for e in equivalent_nodes:
            if G.has_node(e) and e != final_label:
                v1 = G.number_of_nodes()
                G = contract_nodes(G, final_label, e)
                v2 = G.number_of_nodes()
                assert v2 < v1
    v2 = G.number_of_nodes()
    assert v2 < v1
    return G

def infer_contractions(G, out_file=None):
    done = False
    all_contractions = []
    contraction_round = 0
    while not done:
        contraction_round += 1
        contractions = []
        for u1, u2 in combinations(list(G.nodes()), 2):
            n1 = directed_successors(G, u1, direction='out', terminals=None)
            n2 = directed_successors(G, u2, direction='out', terminals=None)
            if sorted(n1) == sorted(n2) and len(n1) > 1:
                contractions.append((u1, u2))

        all_contractions += contractions
        if len(contractions) == 0:
            done = True
        else:
            if out_file != None:
                out_file.write('-----\n')
                out_file.write('round %d\n' % contraction_round)
                out_file.write('-----\n')
                for u1, u2 in sorted(contractions, key = lambda (u1, u2): tuple(sorted((sort_key(u1), sort_key(u2))))):
                    out_file.write("equivalent nodes: %s, %s\n" % (u1, u2))
                    out_file.write("identical neighbors: %s\n" % str(list(G.neighbors(u1))))
                    out_file.write('\n')
                out_file.write('\n')

            G, contraction_map = contract_nodes(G, contractions)

    return all_contractions

def compare_contractions(predicted_contractions, contractions):
    false_positives = set()
    true_positives = set()
    false_negatives = set()
    for u in set(predicted_contractions) | set(contractions):
        x = 0
        if u in predicted_contractions:
            x += 1
        if u in contractions:
            x += 2

        assert x > 0

        if x == 1:
            false_positives.add(u)
        elif x == 2:
            false_negatives.add(u)
        else:
            true_positives.add(u)

def suspicious_directions(G):
    suspicious = []
    for u, v in G.edges():
        p1, p2 = node_to_pair(u), node_to_pair(v)
        if p1 == None or p2 == None:
            continue

        p1 = (p1[0], len(p1[1]))
        p2 = (p2[0], len(p2[1]))

        if p1 < p2 and G[u][v]['direction'] == 'in':
            suspicious.append((u, v))

    return suspicious

def read_network(network, contract=True):
    df = pd.read_csv('%s/%s/network.csv' % (DATA_DIR, network), skipinitialspace=True)
    G = nx.DiGraph()
    G.graph['name'] = network

    G.graph['days'] = []
    for column in list(df):
        column = column.lower()
        if 'used' in column:
            column = column.strip('used ')
            G.graph['days'].append(column)

    prev_edge = None
    for row in df.iterrows():
        row = row[1]
        nodes = row['Nodes']
        if pd.isnull(nodes):
            continue
        nodes = nodes.split('to')
        nodes = map(lambda x : x.strip(), nodes)
        assert len(nodes) >= 2

        used_map = {'yes' : True, 'y' : True, 'Yes' : True, 'Y' : True,\
                    'no' : False, 'n' : False, 'No' : False, 'N' : False}
        days_used = []
        for key, item in row.iteritems():
            key = key.lower()
            if 'used' in key:
                day = key.strip('used ')
                used = None
                if pd.isnull(item):
                    used = False
                else:
                    item = item.strip()
                    used = used_map[item]

                if used:
                    days_used.append(day)

        r1 = row['Repeatability index from nest out']
        r2 = row['Repeatability index toward nest']

        if pd.isnull(r1):
            r1 = DEFAULT_REPEATABILITY
        if pd.isnull(r2):
            r2 = DEFAULT_REPEATABILITY

        same_plant = row['same plant']
        if 'dead end' in nodes[1]:
            same_plant = True
        elif pd.isnull(same_plant):
            same_plant = False
        else:
            same_plant = used_map[same_plant]

        length = row['Distance in cm']
        if pd.isnull(length):
            length = DEFAULT_LENGTH

        for i in xrange(1, len(nodes)):
            n1, n2 = nodes[i - 1], nodes[i]

            n1 = n1.strip()
            n2 = n2.strip()

            pair1 = node_to_pair(n1)
            pair2 = node_to_pair(n2)

            if pair2 < pair1 and 'dead end' not in n2:
                pass #n1, n2 = n2, n1

            G.add_edge(n1, n2)
            G[n1][n2]['repeatability'] = r1
            G[n1][n2]['length'] = length
            G[n1][n2]['days used'] = days_used
            G[n1][n2]['nodes'] = 1
            G[n1][n2]['direction'] = 'out'
            G[n1][n2]['same plant'] = same_plant

            G.add_edge(n2, n1)
            G[n2][n1]['repeatability'] = r2
            G[n2][n1]['length'] = length
            G[n2][n1]['days used'] = days_used
            G[n2][n1]['nodes'] = 1
            G[n2][n1]['direction'] = 'in'
            G[n2][n1]['same plant'] = same_plant

            if prev_edge != None:
                m1, m2 = prev_edge
                if m2 == n1:
                    G.node[n1]['repeatability predecessor'] = {}
                    G.node[n1]['repeatability predecessor']['out'] = m1
                    G.node[n1]['repeatability predecessor']['in'] = n2

            prev_edge = (n1, n2)

    for u in G:
        if 'repeatability predecessor' not in G.node[u]:
            G.node[u]['repeatability predecessor'] = {}
            G.node[u]['repeatability predecessor']['out'] = None
            G.node[u]['repeatability predecessor']['in'] = None

    terminals = []
    food_nodes = set()
    nests = set()

    with open('%s/%s/terminals.csv' % (DATA_DIR, network)) as f:
        for line in f:
            line = line.strip('\n')
            node, terminal = line.split('to')
            node = node.strip()
            if not G.has_node(node):
                for final_label, equivalent_nodes in contraction_map.iteritems():
                    if node in equivalent_nodes:
                        node = final_label
                        break
            assert G.has_node(node)

            terminal = terminal.strip()
            terminals.append(terminal)

            days_used = []
            for day in G.graph['days']:
                for neighbor in G.neighbors(node):
                    if day in G[node][neighbor]['days used']:
                        days_used.append(day)

            G.add_edge(node, terminal)
            G.add_edge(terminal, node)
            G[node][terminal]['repeatability'] = DEFAULT_REPEATABILITY
            G[node][terminal]['length'] = DEFAULT_LENGTH
            G[node][terminal]['days used'] = days_used
            G[node][terminal]['nodes'] = 1
            G[node][terminal]['direction'] = 'in'
            G[node][terminal]['same plant'] = False

            G[terminal][node]['repeatability'] = DEFAULT_REPEATABILITY
            G[terminal][node]['length'] = DEFAULT_LENGTH
            G[terminal][node]['days used'] = days_used
            G[terminal][node]['nodes'] = 1
            G[terminal][node]['direction'] = 'out'
            G[terminal][node]['same plant'] = False

            if 'food' in terminal:
                food_nodes.add(terminal)
            elif 'nest' in terminal:
                nests.add(terminal)

    G.graph['nests'] = list(nests)
    G.graph['food nodes'] = list(food_nodes)
    G.graph['terminals'] = list(food_nodes) + list(nests)
    contractions = []
    with open('%s/%s/contract.csv' % (DATA_DIR, network)) as f:
        for line in f:
            line = line.strip('\n')
            line = line.split('=')
            line = map(lambda x : x.strip(), line)

            assert len(line) >= 2

            for i in xrange(1, len(line)):
                contractions.append((line[i - 1], line[i]))

    contraction_map = get_contraction_map(G, contractions)
    G.graph['contraction map'] = contraction_map
    for u in contraction_map:
        component = contraction_map[u]
        for c in component:
            if G.has_node(c):
                G.node[c]['equivalent nodes'] = component

    for u in G:
        if 'equivalent nodes' not in G.node[u]:
            G.node[u]['equivalent nodes'] = [u]

    if contract:
        G = contract_graph(G, contraction_map)

    if contract:
        G.graph['node ruptures'] = []
        ruptured_nodes_fname = '%s/%s/node_ruptures.csv' % (DATA_DIR, network)
        if os.path.exists(ruptured_nodes_fname):
            with open(ruptured_nodes_fname) as f:
                for line in f:
                    line = line.strip('\n')
                    node, day = line.split(', ')
                    assert G.has_node(node)
                    assert day in G.graph['days']

                    G.graph['node ruptures'].append((node, day))

        G.graph['edge ruptures'] = []
        ruptured_edges_fname = '%s/%s/edge_ruptures.csv' % (DATA_DIR, network)
        if os.path.exists(ruptured_edges_fname):
            with open(ruptured_edges_fname) as f:
                for line in f:
                    line = line.strip('\n')
                    u, v, day = line.split(', ')
                    assert G.has_edge(u, v)
                    assert G.has_edge(v, u)
                    assert day in G.graph['days']

                    G.graph['edge ruptures'].append((u, v, day))

        G.graph['repeatability changes'] = []
        repeatability_changes_fname = '%s/%s/repeatability_changes.csv' % (DATA_DIR, network)
        with open(repeatability_changes_fname) as f:
            for line in f:
                line = line.strip('\n')
                u, v, date, old, new = line.split(', ')
                G.graph['repeatability changes'].append((u, v, date, old, new))

        G.graph['barely'] = []
        barely_fname = '%s/%s/barely.csv' % (DATA_DIR, network)
        with open(barely_fname) as f:
            for line in f:
                line = line.strip('\n')
                day, u, v = line.split(', ')
                G.graph['barely'].append((day, u, v))

    resolutions_fname = '%s/%s/resolutions.csv' % (DATA_DIR, network)
    resolutions = []
    with open(resolutions_fname) as f:
        for line in f:
            x, y, z, r = line.split(', ')
            resolutions.append((x, y, z, int(r)))
    G.graph['resolutions'] = resolutions

    return G

def read_network_paper(colony):

    G = nx.DiGraph()
    G.graph['name'] = colony
    G.graph['nests'] = []
    G.graph['food nodes'] = []
    G.graph['terminals'] = []
    G.graph['days'] = set()
    G.graph['node ruptures'] = []
    G.graph['edge ruptures'] = []
    G.graph['repeatability changes'] = []

    LG = nx.DiGraph()
    LG.graph['name'] = colony

    fname = '%s/%s' % (RAW_DATA_DIR, COLONY_RAW_DATA_FILES[colony])
    with open(fname) as f:
        for line in f:
            if line[0] == '#':
                continue

            line = line.strip('\n')
            line = line.split(',')
            if len(line) <= 1:
                continue

            line = map(lambda x: x.strip(), line)

            if line[0] == 'edge':
                assert len(line) >= 4
                u, v = line[1], line[2]

                length = line[3]
                if length == '':
                    length = DEFAULT_LENGTH
                else:
                    length = float(length)
                days_used = []
                for day in line[4:]:
                    if day != '':
                        days_used.append(day)

                assert not G.has_edge(u, v)

                G.add_edge(u, v)
                G[u][v]['length'] = length
                G[u][v]['days used'] = days_used
                G[u][v]['direction'] = 'out'
                G[u][v]['nodes'] = 1

                G.add_edge(v, u)
                G[v][u]['length'] = length
                G[v][u]['days used'] = days_used
                G[v][u]['direction'] = 'in'
                G[v][u]['nodes'] = 1

                G.graph['days'].update(days_used)

            elif line[0] == 'transition':
                assert len(line) >= 5
                x, y, z = line[1], line[2], line[3]
                assert G.has_node(x)
                assert G.has_node(y)
                assert G.has_node(z)
                assert G.has_edge(x, y)
                assert G.has_edge(y, z)
                assert G[x][y]['direction'] == G[y][z]['direction'] == 'out'

                transition_index = line[4]
                ambiguous = False
                if line[4] == '':
                    transition_index = DEFAULT_REPEATABILITY
                    ambiguous = True
                else:
                    transition_index = int(transition_index)

                u = (x, y)
                v = (y, z)
                assert not LG.has_edge(u, v)

                LG.add_edge(u, v)
                LG[u][v]['transition index'] = transition_index
                LG[u][v]['direction'] = 'out'
                LG[u][v]['ambiguous'] = ambiguous

                LG.add_edge(v, u)
                LG[v][u]['transition index'] = transition_index
                LG[v][u]['direction'] = 'in'
                LG[u][v]['ambiguous'] = ambiguous

            elif line[0] == 'terminal':
                u, terminal = line[1], line[2]

                G.graph['terminals'].append(terminal)
                if 'nest' in terminal:
                    G.graph['nests'].append(terminal)
                elif 'food' in terminal:
                    G.graph['food nodes'].append(terminal)
                else:
                    raise ValueError("invalid type of terminal")

                days_used = set()
                for v in G.neighbors(u):
                    days_used.update(G[u][v]['days used'])
                days_used = sorted(list(days_used))

                G.add_edge(terminal, u)
                G[terminal][u]['days used'] = days_used
                G[terminal][u]['length'] = DEFAULT_LENGTH
                G[terminal][u]['repeatability'] = DEFAULT_REPEATABILITY
                G[terminal][u]['direction'] = 'out'


                G.add_edge(u, terminal)
                G[u][terminal]['days used'] = days_used
                G[u][terminal]['length'] = DEFAULT_LENGTH
                G[u][terminal]['repeatability'] = DEFAULT_REPEATABILITY
                G[u][terminal]['direction'] = 'in'

            elif line[0] == 'equivalence':
                assert len(line) >= 3
                u = line[1]
                assert G.has_node(u)
                equivalent_nodes = []
                for e in line[2:]:
                    if e != '':
                        assert not G.has_node(e)
                        equivalent_nodes.append(e)
                G.node[u]['equivalent nodes'] = equivalent_nodes

            elif line[0] == 'node rupture':
                assert len(line) >= 3
                u = line[1]
                assert G.has_node(u)

                for day in line[2:]:
                    G.graph['node ruptures'].append((u, day))

            elif line[0] == 'edge rupture':
                assert len(line) >= 4
                u, v = line[1], line[2]
                assert G.has_edge(u, v)
                for day in line[3:]:
                    if day != '':
                        G.graph['edge ruptures'].append((u, v, day))

        def sort_days_key(day):
            if 'morning' in day:
                return day.replace("morning", "am")
            elif 'afternoon' in day:
                return day.replace("afternoon", "pm")
            else:
                return day

        G.graph['days'] = sorted(list(G.graph['days']), key=sort_days_key)
        G.graph['terminals'] = list(set(G.graph['terminals']))

        for u in G.nodes():
            if 'equivalent nodes' not in G.node[u]:
                G.node[u]['equivalent nodes'] = []

        for graph in [G, LG]:
            for u, v in graph.edges():
                graph[u][v]['nodes'] = 1

        assert nx.is_strongly_connected(G)
        #assert nx.is_strongly_connected(LG)

        return G, LG


def day_subnetwork(G, day, barely=False):
    assert day in G.graph['days']
    G2 = G.copy()

    for u, v, rupture_day in G.graph['edge ruptures']:
        if rupture_day == day:
            G2.remove_edge(u, v)
            G2.remove_edge(v, u)

    for u, rupture_day in G.graph['node ruptures']:
        if rupture_day == day:
            G2.remove_node(u)

    edgelist = []
    for u, v in G2.edges():
        if day in G2[u][v]['days used']:
            edgelist.append((u, v))

    if barely:
        for day, u, v in G.graph['barely']:
            if (u, v) in edgelist:
                edgelist.remove((u, v))
            if (v, u) in edgelist:
                edgelist.remove((v, u))

    D = G2.edge_subgraph(edgelist)
    D.graph['terminals used'] = []
    for terminal in G2.graph['terminals']:
        if D.has_node(terminal):
            D.graph['terminals used'].append(terminal)

    D.graph['pseudo terminals'] = []
    if 'nest1' not in D.graph['terminals used']:
        closest = min(D.nodes(), key = lambda u : nx.shortest_path_length(G, u, 'nest1'))
        D.graph['pseudo terminals'].append(closest)

    for u, v, change_day, old, new in G.graph['repeatability changes']:
        if D.has_edge(u, v) and day >= change_day:
            D[u][v]['repeatability'] = int(new)

    D.graph['cycles'] = get_cycles(D, D.graph['terminals used'])

    return D

def hop_neighborhood(G, r, hops=2):
    neighborhood_nodes = set()
    queue = set()
    queue.add(r)
    for i in xrange(hops + 1):
        new_queue = set()
        for v in queue:
            neighborhood_nodes.add(v)
            for n in G.neighbors(v):
                if n not in neighborhood_nodes:
                    new_queue.add(n)
        queue = new_queue

    return G.subgraph(neighborhood_nodes)

def draw_hop_neighborhood(colony, r, G=None, L=None, hops=2):
    if G == None:
        G = read_network(colony)

    if L == None:
        L = line_network(G)

    neighborhood = hop_neighborhood(G, r, hops=hops)
    neighborhood = directed_subgraph(neighborhood, direction='out')
    neighborhood = neighborhood.copy()
    original_nodes = set(neighborhood.nodes())
    equivalent_nodes = defaultdict(list)

    for u in original_nodes:
        neighbors = list(neighborhood.neighbors(u))
        if 'equivalent nodes' in neighborhood.node[u]:
            for e in neighborhood.node[u]['equivalent nodes']:
                if e != u:
                    equivalent_nodes[u].append(e)

    for u in neighborhood:
        assert neighborhood.degree(u) > 0

    equivalence_graph = nx.Graph()
    for u, v in neighborhood.edges():
        for e in equivalent_nodes[u]:
            for e2 in [v] + equivalent_nodes[v]:
                neighborhood.add_edge(e, v)
                neighborhood[e][v]['repeatability'] = neighborhood[u][v]['repeatability']
            equivalence_graph.add_edge(u, e)

    pos = nx.kamada_kawai_layout(neighborhood, scale=1)

    edge_labels = {}
    for u, v in neighborhood.edges():
        edge_labels[(u, v)] = neighborhood[u][v]['repeatability']

    nodelist = []
    node_color = []
    for u in neighborhood:
        assert neighborhood.degree(u) > 0
        nodelist.append(u)
        if u in original_nodes:
            node_color.append('r')
        else:
            node_color.append('y')

    pylab.figure()
    nx.draw(neighborhood, pos=pos, with_labels=True, arrows=False, font_size=7, with_arrows=True, nodelist=nodelist, node_color=node_color)
    nx.draw_networkx_edge_labels(neighborhood, pos=pos, edge_labels=edge_labels, font_size=7)
    nx.draw_networkx_edges(equivalence_graph, pos, style='dashed')
    pylab.draw()

    outdir = '%s/%s' % (HOP_DRAWINGS_DIR, colony)
    os.system('mkdir -p %s' % outdir)
    pylab.savefig('%s/node-%s.pdf' % (outdir, r))

    pylab.close()

def reformat_day(day):
    afternoon = 'afternoon' in day
    day = day[:10]
    if afternoon:
        #day += ' pm'
        day += '*'
    return day

def node_to_pair(node):
    digits = []
    done = False
    digit_chars = map(str, range(10))
    i = 0
    while i < len(node) and not done:
        character = node[i]
        if character in digit_chars:
            digits.append(character)
            i += 1
        else:
            done = True

    digits = ''.join(digits)
    if digits == '':
        return None
    prefix = int(digits)
    suffix = node[i:]
    return (prefix, suffix)

def directed_subgraph(G, direction='out', terminals=None):
    if terminals == None:
        terminals = []

    subgraph_edges = []
    for u, v in G.edges():
        if (G[u][v]['direction'] == direction) or (u in terminals) or (v in terminals):
            subgraph_edges.append((u, v))

    return G.edge_subgraph(subgraph_edges)

def directed_predecessors(G, v, direction='out', terminals=None):
    directed_predecessors = []
    for p in G.predecessors(v):
        if G[p][v]['direction'] == direction and (terminals == None or p not in terminals):
            directed_predecessors.append(p)

    return directed_predecessors

def directed_successors(G, v, direction='out', terminals=None):
    directed_successors = []
    for s in G.successors(v):
        if G[v][s]['direction'] == direction and (terminals == None or s not in terminals):
            directed_successors.append(s)

    return directed_successors

def max_repeatability(G, v, direction='in'):
    preds = directed_predecessors(G, v1, direction=direction, terminals=None)
    max_repeatability = None
    for p in preds:
        repeatability = G[p][v1]['repeatability']
        if max_repeatability == None:
            max_repeatability = repeatability
        else:
            max_repeatability = max(max_repeatability, repeatability)
    return max_repeatability

def repeatability_confidence(G, x, y, z, direction='out'):
    assert G.has_edge(x, y)
    assert G.has_edge(y, z)

    if len(directed_predecessors(G, y, direction=direction)) == 1:
        return 3
    elif G.node[y]['repeatability predecessor'][direction] == x:
        return 3

    same1 = G[x][y]['same plant']
    same2 = G[y][z]['same plant']
    r1 = G[x][y]['repeatability']
    r2 = G[y][z]['repeatability']

    if same1 and same2 and r1 == r2 == 1:
        return 2
    elif not same1 and same2:
        return 1

    return 0

def guess_repeatability(G, G_unmerged, x, y, z, direction='out'):
    r1 = G[x][y]['repeatability']
    r2 = G[y][z]['repeatability']
    max_repeatability = max(r1, r2)

    confidence = repeatability_confidence(G, x, y, z, direction=direction)
    if confidence >= 2:
        return r2, confidence
    elif confidence == 1:
        return max_repeatability, confidence

    equiv1 = G.node[x]['equivalent nodes']
    equiv2 = G.node[y]['equivalent nodes']
    equiv3 = G.node[z]['equivalent nodes']

    max_confidence = 0
    max_confidence_repeatability = None
    for a, b, c in product(equiv1, equiv2, equiv3):
        if G_unmerged.has_edge(a, b) and G_unmerged.has_edge(b, c):
            confidence = repeatability_confidence(G_unmerged, a, b, c, direction=direction)
            if confidence > max_confidence:
                max_confidence = confidence
                max_confidence_repeatability = G_unmerged[b][c]['repeatability']

                if confidence == 3:
                    break

    assert 0 <= confidence <= 3

    if max_confidence_repeatability == None:
        max_confidence_repeatability = r2

    return max_confidence_repeatability, max_confidence

def line_network(G, G_unmerged, direction='out'):
    L = nx.line_graph(G)

    ambiguous_edges = defaultdict(list)
    same_plant_edges = defaultdict(list)
    different_plant_edges = defaultdict(list)
    resolution_edges = defaultdict(list)
    unambiguous_edges = defaultdict(list)

    unambiguous_repeatabilities = 0
    same_plant_guesses = 0
    different_plant_guesses = 0
    resolved_repeatabilities = 0
    ambiguous_guesses = 0

    unambiguous_repeatabilities_used = 0
    same_plant_guesses_used = 0
    different_plant_guesses_used = 0
    resolved_repeatabilities_used = 0
    ambiguous_guesses_used = 0

    for (u, v) in list(L.nodes())[:]:
        if u in G.graph['terminals'] or v in G.graph['terminals']:
            L.remove_node((u, v))

    for (x, y, z, r) in G.graph['resolutions']:
        assert G.has_edge(x, y)
        assert G.has_edge(y, z)
        for (u1, v1), (u2, v2) in [((x, y), (y, z))]:
            direction1 = G[u1][v1]['direction']
            direction2 = G[u2][v2]['direction']

            if not (direction1 == direction2 == 'out'):
                L.remove_edge((u1, v1), (u2, v2))
                continue

            L[(u1, v1)][(u2, v2)]['repeatability'] = r
            L[(u1, v1)][(u2, v2)]['confidence'] = 2

            used1 = set(G[u1][v1]['days used'])
            used2 = set(G[u2][v2]['days used'])
            days_used = list(used1 & used2)
            L[(u1, v1)][(u2, v2)]['days used'] = days_used
            length1 = G[u1][v1]['length']
            length2 = G[u2][v2]['length']
            L[(u1, v1)][(u2, v2)]['length'] = length1 + length2

            resolution_edges[y].append(((u1, v1), (u2, v2)))
            resolved_repeatabilities += 1
            if len(days_used) > 0:
                resolved_repeatabilities_used += 1

    for (u1, v1), (u2, v2) in sorted(list(L.edges())):
        assert v1 == u2
        assert G.has_edge(u1, v1)
        assert G.has_edge(u2, v2)
        if 'repeatability' in L[(u1, v1)][(u2, v2)]:
            continue

        direction1 = G[u1][v1]['direction']
        direction2 = G[u2][v2]['direction']

        if not (direction1 == direction2 == direction):
            L.remove_edge((u1, v1), (u2, v2))
            continue

        repeatability, confidence = guess_repeatability(G, G_unmerged, u1, v1, v2, direction=direction)
        assert 0 <= confidence <= 3
        L[(u1, v1)][(u2, v2)]['repeatability'] = repeatability
        L[(u1, v1)][(u2, v2)]['confidence'] = confidence
        length1 = G[u1][v1]['length']
        length2 = G[u2][v2]['length']
        L[(u1, v1)][(u2, v2)]['length'] = length1 + length2

        used1 = set(G[u1][v1]['days used'])
        used2 = set(G[u2][v2]['days used'])
        days_used = list(used1 & used2)
        L[(u1, v1)][(u2, v2)]['days used'] = days_used
        used = len(days_used) > 0

        if confidence == 0:
            ambiguous_edges[v1].append(((u1, v1), (u2, v2)))
            ambiguous_guesses += 1
            if len(days_used) > 0:
                ambiguous_guesses_used += 1
        elif confidence == 1:
            different_plant_edges[v1].append(((u1, v1), (u2, v2)))
            different_plant_guesses += 1
            if len(days_used) > 0:
                different_plant_guesses_used += 1
        elif confidence == 2:
            same_plant_edges[v1].append(((u1, v1), (u2, v2)))
            same_plant_guesses += 1
            if len(days_used) > 0:
                same_plant_guesses_used += 1
        else:
            assert confidence == 3
            unambiguous_edges[v1].append(((u1, v1), (u2, v2)))
            unambiguous_repeatabilities += 1
            if len(days_used) > 0:
                unambiguous_repeatabilities_used += 1

    L.graph['terminals'] = []
    for u, v in L.nodes():
        if u in G.graph['terminals'] or v in G.graph['terminals']:
            L.graph['terminals'].append((u, v))

    L.graph['ambiguous edges'] = ambiguous_edges
    L.graph['same plant edges'] = same_plant_edges
    L.graph['different plant edges'] = different_plant_edges
    L.graph['unambiguous edges'] = unambiguous_edges
    L.graph['resolution edges'] = resolution_edges

    L.graph['unambiguous repeatabilities'] = unambiguous_repeatabilities
    L.graph['ambiguous guesses'] = ambiguous_guesses
    L.graph['same plant guesses'] = same_plant_guesses
    L.graph['different plant guesses'] = different_plant_guesses
    L.graph['resolved repeatabilities'] = resolved_repeatabilities

    L.graph['unambiguous repeatabilities used'] = unambiguous_repeatabilities_used
    L.graph['ambiguous guesses used'] = ambiguous_guesses_used
    L.graph['same plant guesses used'] = same_plant_guesses_used
    L.graph['different plant guesses used'] = different_plant_guesses_used
    L.graph['resolved repeatabilities used'] = resolved_repeatabilities_used

    return L

def write_ambiguous_nodes_file(G, D, f, terminals=None, directions=['in', 'out'], verbose=False):
    ambiguous_nodes = 0
    for u in sorted(D.nodes(), key=lambda u : node_to_pair(u)):
        if terminals != None and u in terminals:
            continue

        for direction in directions:
            in_neighbors = directed_predecessors(G, u, direction=direction, terminals=terminals)

            if len(in_neighbors) > 1:
                out_neighbors = directed_successors(G, u, direction=direction, terminals=terminals)
                if len(out_neighbors) > 0:
                    ambiguous_nodes += 1
                    f.write('%s %s\n' % (u, direction))
                    for p, s in product(in_neighbors, out_neighbors):
                        r1 = G[p][u]['repeatability']
                        r2 = G[u][s]['repeatability']
                        f.write('%s -> %s -> %s: %d -> %d\n' % (p, u, s, r1, r2))
                    f.write('-----\n')

    return ambiguous_nodes

def write_ambiguous_nodes():
    f1 = open('ambiguous_nodes_used.txt', 'w')
    f2 = open('ambiguous_nodes.txt', 'w')
    f3 = open('ambiguous_nodes_unique.txt', 'w')
    directions = ['out']
    for network in COLONIES:
        print network
        f1.write('*****%s*****\n' % network)
        f2.write('*****%s*****\n' % network)
        f3.write('*****%s*****\n' % network)
        G = read_network(network)

        unique_colony_nodes = set()
        for day in G.graph['days']:
            print day
            f1.write('###%s###\n' % day)
            D = day_subnetwork(G, day)
            unique_colony_nodes.update(list(D.nodes()))
            #print network, day, D.graph['terminals used']
            ambiguous_nodes = write_ambiguous_nodes_file(G, D, f1, terminals=D.graph['terminals used'], directions=directions)
            f1.write('\n')
            print ambiguous_nodes, "/", D.number_of_nodes() * len(directions)

        colony_ambiguous_nodes = write_ambiguous_nodes_file(G, G, f2, terminals=G.graph['terminals'], directions=directions, verbose=True)
        print '>', colony_ambiguous_nodes, '/', G.number_of_nodes() * len(directions)

        D_unique = G.subgraph(unique_colony_nodes)
        unique_ambiguous_nodes = write_ambiguous_nodes_file(G, D_unique, f3, terminals=D.graph['terminals'], directions=directions)

        f1.write('\n\n\n')
        f2.write('\n\n\n')
        f3.write('\n\n\n')

    f1.close()
    f2.close()

def print_ambiguous_edges():
    for colony in COLONIES:
        print "*****", colony, "*****"
        G = read_network(colony, contract=True)
        G_unmerged = read_network(colony, contract=False)
        direction = 'out'
        L = line_network(G, G_unmerged, direction=direction)
        print_str = []
        ambiguous_out_nodes = 0
        for u in sorted(G.nodes(), key=node_to_pair):
            node_str = []
            node_str.append("ambiguous node %s" % u)
            if 'repeatability predecessor' in G.node[u]:
                node_str.append('path predecessor %s' % G.node[u]['repeatability predecessor'][direction])
            node_str.append('predecessors %s' % str(directed_predecessors(G, u, direction=direction)))

            ambiguous_edges_present = 0

            dict_keys = ['unambiguous edges', 'resolution edges', 'same plant edges', 'different plant edges', 'ambiguous edges']
            headers = ['unambiguous repeatabilities', 'resolved edges', 'same plant repeatabilities', 'different plant repeatabilities', 'ambiguous_repeatabilities']
            ambiguities = [True, True, True, True, False]
            for dict_key, header, ambiguity in zip(dict_keys, headers, ambiguities):
                node_str.append("-------------------")
                node_str.append(header)
                dictionary = L.graph[dict_key]
                for (u1, v1), (u2, v2) in sorted(dictionary[u], key = lambda (x, y) : len(L[x][y]['days used'])):
                    assert v1 == u2 == u
                    used_prefix = '>>>' * min(len(L[(u1, v1)][(u2, v2)]['days used']), 1)
                    if G[u1][v1]['direction'] == G[u2][v2]['direction'] == direction:
                        edge_str = '%s%s -> %s -> %s' % (used_prefix, u1, v1, v2)
                        if ambiguity:
                            repeatability = L[(u1, v1)][(u2, v2)]['repeatability']
                            edge_str += ': %d' % repeatability
                        else:
                            ambiguous_edges_present += 1
                        node_str.append(edge_str)

            if ambiguous_edges_present > 0:
                print_str.append('\n'.join(node_str))
                print_str.append("\n")
                ambiguous_out_nodes += 1
                draw_hop_neighborhood(colony, u, G=G, L=L, hops=3)

        print "unique ambiguous nodes", ambiguous_out_nodes, '/', G.number_of_nodes()

        unambiguous_repeatabilities = L.graph['unambiguous repeatabilities']
        resolved_repeatabilities = L.graph['resolved repeatabilities']
        same_plant_guesses = L.graph['same plant guesses']
        different_plant_guesses = L.graph['different plant guesses']
        educated_guesses = same_plant_guesses + different_plant_guesses

        ambiguous_guesses = L.graph['ambiguous guesses']

        unambiguous_repeatabilities_used = L.graph['unambiguous repeatabilities used']
        resolved_repeatabilities_used = L.graph['resolved repeatabilities used']
        same_plant_guesses_used = L.graph['same plant guesses used']
        different_plant_guesses_used = L.graph['different plant guesses used']
        ambiguous_guesses_used = L.graph['ambiguous guesses used']

        print "unambiguous repeatabilities: %d / %d" % (unambiguous_repeatabilities, L.number_of_edges())
        print "     %d used by the ants" % unambiguous_repeatabilities_used
        print "resolved repeatabilities: %d / %d" % (L.graph['resolved repeatabilities'], L.number_of_edges())
        print "     %d used by the ants" % resolved_repeatabilities_used
        print "educated repeatabilities guesses: %d / %d" %  (educated_guesses, L.number_of_edges())
        print "     same plant guesses: %d / %d" % (same_plant_guesses, L.number_of_edges())
        print "          %d used by the ants" % same_plant_guesses_used
        print "     different plant guesses: %d / %d" % (different_plant_guesses, L.number_of_edges())
        print "          %d used by the ants" % different_plant_guesses_used
        print "ambiguous repeatabilities: %d / %d" % (ambiguous_guesses, L.number_of_edges())
        print "     %d used by the ants" % ambiguous_guesses_used
        print "************************\n"
        print '\n'.join(print_str)

def node_edge_stats():
    vegetation_nodes = []
    vegetation_edges = []
    day_nodes = defaultdict(list)
    day_edges = defaultdict(list)
    for colony in COLONIES:
        G = read_network(colony)
        vegetation_nodes.append(G.number_of_nodes())
        vegetation_edges.append(G.number_of_edges())

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            day_nodes[colony].append(D.number_of_nodes())
            day_edges[colony].append(D.number_of_edges())

    for column, vals in zip(['vegetation nodes', 'vegetation edges'], [vegetation_nodes, vegetation_edges]):
        print column
        print 'min:', min(vals)
        print 'max:', max(vals)
        print 'mean (std):', pylab.mean(vals), pylab.std(vals, ddof=1)
        print '\n'

    nodes_used = []
    edges_used = []
    for column, dictionary, final_list in zip(['nodes used', 'edges used'], [day_nodes, day_edges], [nodes_used, edges_used]):
        print column
        for colony, counts in dictionary.iteritems():
            mu = pylab.mean(counts)
            sigma = pylab.std(counts, ddof=1)
            final_list.append(mu)
            print colony, ':', mu, '+/-', sigma

        print 'min:', min(final_list)
        print 'max:', max(final_list)
        print 'mean (std):', pylab.mean(final_list), pylab.std(final_list, ddof=1)
        print '\n'

def get_depth(node_name):
    depth = 0
    for char in node_name:
        if char in 'ABCDEF':
            depth += 1
    return depth

def aggregate_network_data():
    for colony in COLONIES:
        G = read_network(colony)
        G_unmerged = read_network(G, contract=False)
        L = line_network(G, G_unmerged, direction='out')

        colony_name = COLONY_LABELS[colony].replace(' ', '_')
        fname = '%s/%s_network_data2018.csv' % (RAW_DATA_DIR, colony_name)
        with open(fname,'w') as f:
            for u, v in G.edges():
                if G[u][v]['direction'] == 'out' and u not in G.graph['terminals'] and v not in G.graph['terminals']:
                    line_str = ['edge', u, v, str(G[u][v]['length'])]
                    for day in G[u][v]['days used']:
                        line_str.append(day)
                    line_str = ', '.join(line_str)
                    f.write('%s\n' % line_str)

            f.write('\n\n')
            for u, v in L.edges():
                a, b = u
                c, d = v
                assert b == c
                f.write('transition, %s, %s, %s, %d\n' % (a, c, d, L[u][v]['repeatability']))

            f.write('\n\n')
            for terminal in G.graph['terminals']:
                for n in G.neighbors(terminal):
                    f.write('terminal, %s, %s\n' % (n, terminal))

            f.write('\n\n')
            for u in G.nodes():
                line_str = ['equivalence', u]
                for e in G.node[u]['equivalent nodes']:
                    if e != u:
                        line_str.append(e)
                if len(line_str) > 2:
                    line_str = ', '.join(line_str)
                    f.write('%s\n' % line_str)

            f.write('\n\n')
            node_ruptures = defaultdict(list)
            for node, day in G.graph['node ruptures']:
                node_ruptures[node].append(day)

            for node, days in node_ruptures.iteritems():
                line_str = ['node rupture', node] + days
                line_str = ', '.join(line_str)
                f.write('%s\n' % line_str)

            f.write('\n\n')
            edge_ruptures = defaultdict(list)
            for u, v, day in G.graph['edge ruptures']:
                edge_ruptures[(u, v)].append(day)

            for edge, days in edge_ruptures.iteritems():
                u, v = edge
                line_str = ['edge rupture', u, v] + days
                line_str = ', '.join(line_str)
                f.write('%s\n' % line_str)

def print_ambiguous_transitions(G, LG):
    for x, y in sorted(LG.edges()):
        a, b = x
        c, d = y
        if b == c and LG[x][y]['ambiguous'] == True:
            print 'ambiguous transition: %s --> %s --> %s' % (a, b, d)

def print_missing_transitions(G, LG):
    terminals = G.graph['terminals']
    for u in G.nodes():
        for p in directed_predecessors(G, u, direction='out', terminals=terminals):
            for s in directed_successors(G, u, direction='out', terminals=terminals):
                x = (p, u)
                y = (u, s)
                if not LG.has_edge((p, u), (u, s)):
                    print "missing transition: %s -> %s -> %s" % (p, u, s)

def print_missing_data(colonies):
    if colonies == None or len(colonies) == 0:
        colonies = COLONIES

    for colony in colonies:
        print colony
        G, LG = read_network_paper(colony)
        assert nx.is_strongly_connected(G)
        assert nx.is_strongly_connected(LG)

        print_ambiguous_transitions(G, LG)
        print_missing_transitions(G, LG)
        print '\n'

def read_all_networks():
    for colony in COLONIES:
        G, LG = read_network_paper(colony)
        print G.graph['edge ruptures']

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--t1', action='store_true')
    parser.add_argument('-n', '--network_data', action='store_true')
    parser.add_argument('--missing', action='store_true')
    parser.add_argument('--colonies', nargs='*')
    parser.add_argument('-r', '--read', action='store_true')

    args = parser.parse_args()

    if args.t1:
        node_edge_stats()

    if args.network_data:
        aggregate_network_data()

    if args.missing:
        print_missing_data(args.colonies)

    if args.read:
        read_all_networks()

if __name__ == '__main__':
    main()
