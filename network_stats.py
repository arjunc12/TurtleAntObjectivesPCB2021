from map_network import *
import networkx as nx
import matplotlib as mpl
mpl.use('agg')
import pylab
import seaborn as sns
sns.set_style('white')
from collections import defaultdict
import os
from scipy.stats import sem, pearsonr, spearmanr, kendalltau, kruskal, mannwhitneyu
import pandas as pd
from collections import defaultdict
from itertools import combinations
from utils import node_order
import argparse

FIGS_DIR = 'network_stats'
MAX_REPEATABILITY = 4

CSV_DIR = 'paper_data/plot_data/'

def bootstrap_clustering(G):
    '''

    '''
    deviations = []
    for sample_size in xrange(10, G.number_of_nodes(), 10):
        vals = []
        for i in xrange(10000):
            nbunch = sample(G.nodes(), sample_size)
            vals.append(nx.average_clustering(G.subgraph(nbunch)))
        print sample_size
        print pylab.std(vals, ddof=1)

def plot_cycle_counts():
    day_lists = []
    count_lists = []
    retained_lists = []
    all_days = set()
    for colony in COLONIES:
        print colony
        G, LG = read_network_paper(colony)
        days = []
        counts = []
        retained = []
        prev_cycles = None
        for day in G.graph['days']:
            print day
            D = day_subnetwork(G, day)
            days.append(reformat_day(day))
            cycles = set(D.graph['cycles'])
            counts.append(len(cycles))

            if prev_cycles != None:
                retained_cycles = cycles.intersection(prev_cycles)
                retained.append(len(retained_cycles))

            prev_cycles = set(list(cycles)[:])

        day_lists.append(days)
        all_days.update(days)
        count_lists.append(counts)
        retained_lists.append(retained)

    all_days = sorted(all_days)
    #sns.set()

    for yvals, title_str, fname in zip([count_lists, retained_lists], ['unique cycles formed', 'cycles retained'], ['cycle_counts', 'cycles_retained']):
        pylab.figure()
        for colony, days, y in zip(COLONIES, day_lists, yvals):
            x = []
            for day in days:
                x.append(all_days.index(day))

            len_diff = len(x) - len(y)
            assert len_diff >= 0
            x = x[len_diff:]

            pylab.scatter(x, y, label='_nolegend_')
            pylab.plot(x, y, label=colony)

        pylab.title(title_str, size=20)

        pylab.xticks(range(len(all_days)), all_days)
        pylab.tick_params(axis='both', labelsize=20)
        pylab.tick_params(axis='x', rotation=-90)
        #pylab.xlabel('observation day', size=20)
        pylab.ylabel('cycles', size=20)

        '''
        leg = pylab.legend()
        ax = pylab.gca()
        pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
        leg.get_frame().set_linewidth(5)
        leg.get_frame().set_edgecolor('k')
        '''

        pylab.tight_layout()

        pylab.savefig('%s/%s.pdf' % (FIGS_DIR, fname), format='pdf')
        pylab.close()

def plot_multipath_repeatabilities():
    connectivity_df = pd.read_csv('local_connectivities_repeatability.csv', skipinitialspace=True)
    #connectivity_df['connectivity'] **= -1

    df_dict = {'colony' : [], 'day' : [], 'source' : [], 'target' : [], 'repeatability' : [], 'cycle nodes' : []}
    colony_order = []
    for colony in COLONIES:
        colony_label = COLONY_LABELS[colony]
        colony_order.append(colony_label)

        G, LG = read_network_paper(colony)

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            D2 = D.copy()
            for u, v in D.edges():
                if D[u][v]['direction'] == 'in':
                    D2.remove_edge(u, v)

            for cycle in D.graph['cycles']:
                for source, target in combinations(cycle[:-1], 2):
                    source, target = sorted([source, target])
                    paths = list(nx.all_simple_paths(D2, source, target)) + list(nx.all_simple_paths(D2, target, source))
                    viable_paths = []
                    for path in paths:
                        if len(path) >= 3:
                            viable_paths.append(path)
                    if len(viable_paths) >= 2:
                        repeatabilities = []
                        cycle_nodes = set()
                        for path in viable_paths:
                            repeatabilities.append(path_repeatability(L, path))
                            cycle_nodes.update(path)
                        repeatability = pylab.mean(repeatabilities)
                        #size = sum(sizes)
                        path_nodes = '-'.join(sorted(list(cycle_nodes)))
                        df_dict['colony'].append(colony_label)
                        df_dict['day'].append(reformat_day(day))
                        df_dict['source'].append(source)
                        df_dict['target'].append(target)
                        df_dict['repeatability'].append(repeatability)
                        df_dict['cycle nodes'].append(path_nodes)

    df = pd.DataFrame.from_dict(df_dict)
    df.drop_duplicates(subset=['colony', 'source', 'target', 'cycle nodes'])

    for colony, colony_group in df.groupby('colony'):
        repeatability = colony_group['repeatability']
        print colony, "loop repeatabilities", pylab.mean(repeatability), "+/-", pylab.std(repeatability, ddof=1)
        colony_connectivity = connectivity_df[connectivity_df['colony'] == colony]
        connectivity = colony_connectivity['connectivity']
        print mannwhitneyu(connectivity, repeatability)

    #sns.set()
    pylab.figure()
    #sns.boxenplot(x='day', y='size', hue='colony', data=df)
    sns.barplot(x='colony', y='repeatability', data=df)
    pylab.xlabel('')
    pylab.ylabel('cycle repeatability (index)', size=20)

    '''
    leg = pylab.legend()
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    leg.get_frame().set_linewidth(5)
    leg.get_frame().set_edgecolor('k')
    '''

    pylab.tick_params(axis='both', labelsize=20)
    #pylab.tick_params(axis='x', rotation=-90)
    pylab.title('Distribution of sizes of ant cycles', size=20)
    pylab.tight_layout()
    pylab.savefig('network_stats/multipath_sizes.pdf', format='pdf')
    pylab.close()

def plot_multipath_sizes():
    df_dict = {'colony' : [], 'day' : [], 'source' : [], 'target' : [], 'size' : [], 'cycle nodes' : []}
    colony_order = []
    for colony in COLONIES:
        colony_label = COLONY_LABELS[colony]
        colony_order.append(colony_label)

        G, LG = read_network_paper(colony)

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            D2 = D.copy()
            for u, v in D.edges():
                if D[u][v]['direction'] == 'in':
                    D2.remove_edge(u, v)

            for cycle in D.graph['cycles']:
                for source, target in combinations(cycle[:-1], 2):
                    source, target = sorted([source, target])
                    paths = list(nx.all_simple_paths(D2, source, target)) + list(nx.all_simple_paths(D2, target, source))
                    if len(paths) >= 2:
                        sizes = []
                        cycle_nodes = set()
                        for path in paths:
                            sizes.append(len(path))
                            cycle_nodes.update(path)
                        size = pylab.mean(sizes)
                        #size = sum(sizes)
                        path_nodes = '-'.join(sorted(list(cycle_nodes)))
                        df_dict['colony'].append(colony_label)
                        df_dict['day'].append(reformat_day(day))
                        df_dict['source'].append(source)
                        df_dict['target'].append(target)
                        df_dict['size'].append(size)
                        df_dict['cycle nodes'].append(path_nodes)

    df = pd.DataFrame.from_dict(df_dict)
    df.drop_duplicates(subset=['colony', 'source', 'target', 'cycle nodes'])

    for colony, colony_group in df.groupby('colony'):
        size = colony_group['size']
        print colony, "loop sizes", pylab.mean(size), "+/-", pylab.std(size, ddof=1)

    #sns.set()
    pylab.figure()
    #sns.boxenplot(x='day', y='size', hue='colony', data=df)
    sns.barplot(x='colony', y='size', data=df)
    pylab.xlabel('')
    pylab.ylabel('cycle size (nodes)', size=20)

    '''
    leg = pylab.legend()
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    leg.get_frame().set_linewidth(5)
    leg.get_frame().set_edgecolor('k')
    '''

    pylab.tick_params(axis='both', labelsize=20)
    #pylab.tick_params(axis='x', rotation=-90)
    pylab.title('Distribution of sizes of ant cycles', size=20)
    pylab.tight_layout()
    pylab.savefig('network_stats/multipath_sizes.pdf', format='pdf')
    pylab.close()

def plot_cycle_sizes():
    df_dict = {'colony' : [], 'day' : [], 'cycle size' : []}
    all_sizes = []
    for colony in COLONIES:
        print "*****", colony, "*****"
        colony_sizes = []
        G, LG = read_network_paper(colony)
        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            day_sizes = []
            cycles = D.graph['cycles']
            if len(cycles) > 0:
                for cycle in D.graph['cycles']:
                    cycle_size = len(cycle)
                    day_sizes.append(cycle_size)
                    df_dict['colony'].append(colony)
                    df_dict['day'].append(reformat_day(day))
                    df_dict['cycle size'].append(cycle_size)

                #print day, pylab.mean(day_sizes), "+/-", pylab.std(day_sizes, ddof=1)
                colony_sizes += day_sizes
        if len(colony_sizes) > 0:
            print "-" * len(colony)
            print colony, pylab.mean(colony_sizes), "+/-", pylab.std(colony_sizes, ddof=1)
            all_sizes += colony_sizes
    print "*****", "all colonies", "*****"
    print pylab.mean(all_sizes), "+/-", pylab.std(all_sizes, ddof=1)

    df = pd.DataFrame.from_dict(df_dict)
    #sns.set()
    pylab.figure()
    ax = sns.boxenplot(x='day', y='cycle size', hue='colony', data=df)
    pylab.xlabel('size of cycles')
    ax.tick_params(axis='both', labelsize=20)
    ax.tick_params(axis='x', rotation=-90)

    '''
    leg = pylab.legend()
    ax = pylab.gca()
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    leg.get_frame().set_linewidth(5)
    leg.get_frame().set_edgecolor('k')
    '''

    pylab.tight_layout()
    pylab.savefig('network_stats/cycle_sizes.pdf', format='pdf')
    pylab.close()

def path_transition_index(LG, path):
    if len(path) < 3:
        print path
    assert len(path) >= 3
    total_ti = 0
    transitions = 0

    for i in xrange(len(path) - 2):
        x, y, z = path[i], path[i + 1], path[i + 2]
        if not LG.has_edge((x, y), (y, z)):
            print path
            print x, y, z
        assert LG.has_edge((x, y), (y, z))
        total_ti += LG[(x, y)][(y, z)]['transition index']
        transitions += 1

    return total_ti / float(transitions)

def local_ti_connectivity(G, LG, D, u, n):
    assert G.has_edge(u, n)
    assert D.has_node(u)
    assert not D.has_node(n)

    G2 = G.copy()
    G2.remove_edge(u, n)
    min_ti_connectivity = None
    min_ti_node = None
    for v in G2.nodes():
        if nx.has_path(G2, n, v):
            for path in nx.all_simple_paths(G2, n, v):
                path = [u] + path
                ti_connectivity = path_transition_index(LG, path)
                if min_ti_connectivity == None or ti_connectivity < min_ti_connectivity:
                    min_ti_connectivity = ti_connectivity
                    min_ti_node = v

    return min_ti_connectivity, min_ti_node

def local_connectivity(G, D, u, n):
    assert G.has_edge(u, n)
    assert D.has_node(u)
    assert not D.has_node(n)
    G.remove_edge(u, n)
    min_connectivity = None
    min_connectivity_node = None
    for v in D.nodes():
        if v in G.graph['terminals']:
            continue
        if nx.has_path(G, n, v):
            sp = nx.shortest_path(G, n, v)
            #sp = [u] + sp
            connectivity = len(sp)
            if min_connectivity == None or connectivity < min_connectivity:
                min_connectivity = connectivity
                min_connectivity_node = v

    G.add_edge(u, n)
    return min_connectivity, min_connectivity_node

def local_connectivities(G, L, D, weight='nodes'):
    connectivity_func = None
    if weight == 'transition index':
        connectivity_func = lambda u, n : local_ti_connectivity(G, L, D, u, n)
    else:
        connectivity_func = lambda u, n : local_connectivity(G, D, u, n)
    connectivities = {}
    for u in D.nodes():
        if u in G.graph['terminals']:
            continue
        for n in G.neighbors(u):
            if not D.has_node(n):
                connectivity, connectivity_node = connectivity_func(u, n)
                if connectivity != None:
                    assert connectivity_node != None
                    pair = tuple(sorted([u, connectivity_node]))
                    connectivities[pair] = connectivity

    return connectivities

def write_connectivities_file(connectivity_fname, weight='nodes'):
    with open(connectivity_fname, 'w') as connectivity_file:
        connectivity_file.write('colony, day, node, connectivity\n')
        for colony in COLONIES:
            print colony
            G, LG = read_network_paper(colony)
            G2 = G.copy()

            colony_str = COLONY_LABELS[colony]

            for u, v in G.edges():
                if G[u][v]['direction'] == 'in':
                    G2.remove_edge(u, v)

            for u in G.nodes():
                if u in G.graph['terminals']:
                    G2.remove_node(u)

            for day in G.graph['days']:
                print day
                D = day_subnetwork(G, day)
                connectivities = local_connectivities(G2, LG, D, weight=weight)
                for (u, v), connectivity in sorted(connectivities.iteritems(), key = lambda ((x, y), c) : c):
                    assert u != None
                    assert v != None

                    connectivity_file.write('%s, %s, %s, %f\n' % (colony_str,\
                                                                  reformat_day(day),\
                                                                  u + '-' + v,\
                                                                  connectivity))

def plot_connectivity_distribution(transition=False):
    print "connectivity distribution"
    df_dict = {'colony' : [], 'day' : [], 'node' : [], 'connectivity' : []}

    connectivity_prefix = 'local_connectivities'
    connectivity_weight = 'length'
    if transition:
        connectivity_prefix += '_transition'
        connectivity_weight = 'transition index'

    connectivity_fname = '%s/%s.csv' % (CSV_DIR, connectivity_prefix)

    if not os.path.exists(connectivity_fname):
        write_connectivities_file(connectivity_fname, weight=connectivity_weight)

    df = pd.read_csv(connectivity_fname, skipinitialspace=True)
    #df['connectivity'] **= -1
    #df['colony'] = df['colony'].map(COLONY_LABELS)

    all_connectivities = []
    for colony, colony_group in df.groupby('colony'):
        colony_connectivity = colony_group['connectivity']
        #colony_connectivity **= -1
        print colony, "connectivity", pylab.mean(colony_connectivity), "+/-", pylab.std(colony_connectivity, ddof=1)
        all_connectivities.append(list(colony_connectivity))

    print "connectivity kruskal-wallis", kruskal(*all_connectivities)

    ax = sns.barplot(x='colony', y='connectivity', data=df, order=sorted(COLONY_LABELS.values()))

    pylab.tick_params(axis='both', labelsize=20)
    pylab.xlabel('')
    pylab.ylabel('local connectivity (1 / nodes)', size=20)
    pylab.savefig('%s/%s_distribution.pdf' % (FIGS_DIR, connectivity_prefix), bbox_inches='tight')
    pylab.close()

def degree_distribution(G):
    proportions = defaultdict(float)
    n = G.number_of_nodes()
    for u in G.nodes():
        degree = G.degree(u) / 2
        proportions[degree] += 1.0 / n

    return proportions

def plot_degree_distribution():
    df_dict = {'colony' : [], 'degree' : []}
    hue_order = []
    for colony in COLONIES:
        G, LG = read_network_paper(colony)
        G = G.to_undirected()

        #colony = colony.replace('_', ' ')
        colony = COLONY_LABELS[colony]
        hue_order.append(colony)

        degrees = []
        for u in G.nodes():
            pair = node_to_pair(u)
            if pair == None:
                continue

            node, level = pair
            if len(level) > 8:
                continue

            degree = len(list(G.neighbors(u)))
            if degree == 12:
                print u

            if degree == 1:
                continue

            df_dict['colony'].append(colony)

            df_dict['degree'].append(degree)
            degrees.append(degree)

        print colony, "degrees", pylab.mean(degrees), "+/-", pylab.std(degrees, ddof=1)


    df = pd.DataFrame.from_dict(df_dict)
    all_degrees = []
    for colony, colony_group in df.groupby('colony'):
        all_degrees.append(list(colony_group['degree']))
    print "degree kruskal-wallis", kruskal(*all_degrees)

    pylab.figure()
    df1 = df.groupby(['colony', 'degree']).size().reset_index(name='count')
    df2 = df1.groupby(['colony'], as_index=False).agg('sum').rename(columns={'count' : 'total'})
    df2.drop(columns='degree', inplace=True)
    df = pd.merge(df1, df2)
    df['proportion'] = df['count'] / df['total']
    ax = sns.barplot(x='degree', y='proportion', hue='colony', data=df, hue_order=hue_order)
    ax.get_legend().remove()
    pylab.legend(loc='upper right')
    pylab.tick_params(axis='both', labelsize=20)

    pylab.tight_layout()

    pylab.xlabel('Node degree', size=20)
    pylab.ylabel('Proportion of nodes\nwith this degree', size=20)
    pylab.tight_layout()
    pylab.savefig('%s/degree_distribution_colony.pdf' % FIGS_DIR)
    pylab.close()

def network_clustering(colony):
    G, LG = read_network_paper(colony)
    G = G.to_undirected()
    print colony, nx.average_clustering(G)

def nodes_length_correlation():
    pylab.figure()
    all_node_counts = []
    all_total_lengths = []
    for colony in os.listdir('mapping_network/csv'):
        node_counts = []
        total_lengths = []
        G, LG = read_network_paper(colony)
        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            node_count = D.number_of_nodes()
            total_length = 0
            for u, v in D.edges():
                total_length += D[u][v]['length']
            node_counts.append(node_count)
            total_lengths.append(total_length)
        print colony, pearsonr(node_counts, total_lengths), spearmanr(node_counts, total_lengths)
        pylab.scatter(node_counts, total_lengths, label=colony)
        all_node_counts += node_counts
        all_total_lengths += total_lengths
    print "all colonies", pearsonr(all_node_counts, all_total_lengths), spearmanr(all_node_counts, all_total_lengths)

    pylab.xlabel('total nodes')
    pylab.ylabel('total length')
    pylab.savefig('%s/nodes_length_correlation.pdf' % FIGS_DIR, format='pdf')

def objective_correlations():
    def mean_repeatability(G, L, direction='out'):
        H = directed_subgraph(G)
        H = nx.line_graph(H)
        total = 0
        transitions = 0
        for (u1, v1), (u2, v2) in H.edges():
            assert v1 == u2
            assert G[u1][v1]['direction'] == G[u2][v2]['direction'] == 'out'
            total += L[(u1, v1)][(u2, v2)]['transition index']
            transitions += 1

        return total / float(transitions)

    def total_length(G):
        total = 0.0
        for u, v in G.edges():
            total += G[u][v]['length']
        return total

    def mean_length(G):
        return total_length(G) / G.number_of_edges()

    def total_nodes(G):
        return float(G.number_of_nodes()) #/ (len(G.graph['terminals used']) ** 1)

    def length_repeatability(G):
        return mean_repeatability(G) * mean_length(G)

    def nodes_repeatability(G):
        return mean_repeatability(G) * G.number_of_nodes()

    def nodes_length(G):
        return mean_length(G) * G.number_of_nodes()

    objectives = {'length' : mean_length, 'total length' : total_length,\
                  'nodes' : total_nodes, 'repeatability' : None}


    objective_vals = {}
    for objective in objectives.keys():
        objective_vals[objective] = []

    for colony in COLONIES:
        print "*****" + colony + "*****"
        G, LG = read_network_paper(colony)

        def colony_mean_repeatability(D):
            D2 = D.copy()
            for t in G.graph['terminals']:
                if D2.has_node(t):
                    D2.remove_node(t)
            return mean_repeatability(D2, LG)
        objectives['repeatability'] = colony_mean_repeatability

        colony_objective_vals = {}
        for objective in objectives.keys():
            colony_objective_vals[objective] = []

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            for objective, objective_func in objectives.iteritems():
                val = objective_func(D)
                colony_objective_vals[objective].append(val)
                objective_vals[objective].append(val)

        for objective1, objective2 in combinations(colony_objective_vals.keys(), 2):
            print objective1, "vs", objective2
            vals1, vals2 = colony_objective_vals[objective1], colony_objective_vals[objective2]
            print spearmanr(vals1, vals2)

    print "*****all vals*****"
    for objective1, objective2 in combinations(objective_vals.keys(), 2):
        print objective1, "vs", objective2
        vals1, vals2 = objective_vals[objective1], objective_vals[objective2]
        print spearmanr(vals1, vals2)

def make_proportions_df(df, col):
    df1 = df.groupby(['colony', col]).size().reset_index(name='count')
    df2 = df1.groupby(['colony'], as_index=False).agg('sum').rename(columns={'count' : 'total'})
    df2.drop(columns=col, inplace=True)
    df = pd.merge(df1, df2)
    df['proportion'] = df['count'] / df['total']
    return df

def edge_length_distribution():
    #df_dict = {'colony' : [], 'length' : []}
    fname = '%s/%s' % (CSV_DIR, 'edge_lengths.csv')
    if not os.path.exists(fname):
        with open(fname, 'w') as f:
            f.write('colony, u, v, length\n')
            for colony in COLONIES:
                G, LG = read_network_paper(colony)
                for u, v in G.edges():
                    length = G[u][v]['length']
                    #df_dict['colony'].append(COLONY_LABELS[colony])
                    #df_dict['length'].append(length)
                    f.write('%s, %s, %s, %d\n' % (colony, u, v, length))

    #df = pd.DataFrame.from_dict(df_dict)
    df = pd.read_csv(fname, skipinitialspace=True)
    all_lengths = []
    for colony, colony_group in df.groupby('colony'):
        lengths = colony_group['length']
        print "edge lengths", colony, pylab.mean(lengths), "+/-", pylab.std(lengths, ddof=1)
        all_lengths.append(list(lengths))
    print "edge lengths kruskal-wallis", kruskal(*all_lengths)

    pylab.figure()
    sns.barplot(x='colony', y='length', data=df, order=sorted(COLONY_LABELS.values()))
    pylab.xlabel('')
    pylab.ylabel('edge length (cm)', size=20)
    pylab.tick_params(axis='both', labelsize=20)
    pylab.tight_layout()
    pylab.savefig('%s/edge_length_distribution.pdf' % FIGS_DIR, format='pdf', bbox_inches='tight')
    pylab.close()

def used_edge_length_distribution():
    print "lengths distribution"
    all_lengths = []
    all_weights = []
    labels = []
    for colony in COLONIES:
        colony_lengths = []
        colony_weights = []
        G, LG = read_network_paper(colony)

        lengths = []
        for u, v in G.edges():
            lengths.append(G[u][v]['length'])
        print colony, "lengths", pylab.mean(lengths), "+/-", pylab.std(lengths, ddof=1)

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            lengths = []
            for u, v in D.edges():
                lengths.append(D[u][v]['length'])

            colony_lengths += lengths

        all_lengths.append(colony_lengths)
        colony_weights = pylab.ones_like(colony_lengths) / float(len(colony_lengths))
        all_weights.append(colony_weights)

        labels.append(colony.replace('_', ' '))

    pylab.figure()
    pylab.hist(all_lengths, weights=all_weights, label=labels)

    pylab.xlabel('edge length (cm)', size=20)
    pylab.ylabel('proportion of edge lengths', size=20)
    pylab.tick_params(axis='both', labelsize=20)
    pylab.tight_layout()
    pylab.savefig('%s/edge_length_distribution.pdf' % FIGS_DIR, format='pdf', bbox_inches='tight')
    pylab.close()

def repeatability_distribution():
    print "Transition Index distribution"
    df_dict = {'colony' : [], 'transition index' : []}
    for colony in COLONIES:
        G, LG = read_network_paper(colony)

        repeatabilities = []
        for u, v in LG.edges():
            repeatability = LG[u][v]['transition index']
            repeatabilities.append(repeatability)
            df_dict['colony'].append(COLONY_LABELS[colony])
            df_dict['transition index'].append(repeatability)

    df = pd.DataFrame.from_dict(df_dict)
    all_repeatabilities = []
    for colony, colony_group in df.groupby('colony'):
        repeatabilities = colony_group['transition index']
        print colony, "transition index", pylab.mean(repeatabilities), "+/-", pylab.std(repeatabilities, ddof=1)
        all_repeatabilities.append(list(repeatabilities))
    print "transition index kruskal-wallis", kruskal(*all_repeatabilities)

    pylab.figure()
    sns.barplot(x='colony', y='transition index', data=df, order=sorted(COLONY_LABELS.values()))
    pylab.tick_params(axis='both', labelsize=20)
    pylab.xlabel('')
    pylab.ylabel('transition index', size=20)
    pylab.tight_layout()
    pylab.savefig('%s/repeatability_distribution_colony.pdf' % FIGS_DIR, format='pdf')
    pylab.close()

def used_repeatability_distribution():
    print "repeatability distribution"
    df_dict = {'colony' : [], 'day' : [], 'repeatability' : []}
    hue_order = []
    for colony in COLONIES:
        G, LG = read_network_paper(colony)

        colony = colony.replace('_', ' ')
        hue_order.append(colony)

        repeatabilities = []
        for u, v in L.edges():
            repeatability = L[u][v]['repeatability']
            repeatabilities.append(repeatability)
        print colony, "repeatability", pylab.mean(repeatabilities), "+/-", pylab.std(repeatabilities, ddof=1)
        repeatabilities = filter(lambda x : x <= 2, repeatabilities)


        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            day = day[:10]
            for (u1, v1), (u2, v2) in nx.line_graph(D).edges():
                assert v1 == u2
                if L.has_edge((u1, v1), (u2, v2)):
                    repeatability = L[(u1, v1)][(u2, v2)]['repeatability']
                    assert type(repeatability) is not str
                    #df_dict['colony'].append(colony.replace('_', ' '))
                    df_dict['colony'].append(colony)
                    df_dict['day'].append(day)
                    df_dict['repeatability'].append(repeatability)

    df = pd.DataFrame.from_dict(df_dict)

    days_order = sorted(list(df['day'].unique()))
    pylab.figure()
    ax = sns.barplot(x='day', y='repeatability', hue='colony', hue_order=COLONIES, data=df, order=days_order)
    ax.tick_params(axis='x', rotation=75, labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    pylab.ylabel('repeatability index', size=20)

    pylab.ylim(0, 1)
    pylab.tight_layout()
    pylab.savefig('%s/repeatability_distribution_day.pdf' % FIGS_DIR, format='pdf')
    pylab.close()

    pylab.figure()

    df1 = df.groupby(['colony', 'repeatability']).size().reset_index(name='count')
    df2 = df1.groupby(['colony'], as_index=False).agg('sum').rename(columns={'count' : 'total'})
    df2.drop(columns='repeatability', inplace=True)
    df = pd.merge(df1, df2)
    df['proportion'] = df['count'] / df['total']
    print df
    ax = sns.barplot(x='repeatability', y='proportion', hue='colony', data=df, hue_order=hue_order)
    ax.get_legend().remove()
    pylab.tick_params(axis='both', labelsize=20)
    pylab.xlabel('transition index', size=20)
    pylab.ylabel('proportion of transition indices', size=20)
    pylab.ylim(0, 1)

    pylab.tight_layout()
    pylab.savefig('%s/repeatability_distribution_colony.pdf' % FIGS_DIR, format='pdf')
    pylab.close()

def used_nodes_distribution():
    df_dict = {'colony' : [], 'day' : [], 'nodes' : []}
    for colony in COLONIES:
        G, LG = read_network_paper(colony)

        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            n = D.number_of_nodes()
            df_dict['colony'].append(COLONY_LABELS[colony])
            df_dict['day'].append(day)
            df_dict['nodes'].append(n)

    df = pd.DataFrame.from_dict(df_dict)
    all_nodes = []
    for colony, colony_group in df.groupby('colony'):
        n = colony_group['nodes']
        all_nodes.append(list(n))
        print "total nodes", colony, pylab.mean(n), "+/-", pylab.std(n, ddof=1)

    print "total nodes kruskal-wallis", kruskal(*all_nodes)

    pylab.figure()
    ax = sns.barplot(x='colony', y='nodes', data=df, order=sorted(COLONY_LABELS.values()))
    pylab.tick_params(axis='both', labelsize=20)
    pylab.xlabel('')
    pylab.ylabel('nodes in foraging network', size=20)
    pylab.tight_layout()
    pylab.savefig('%s/used_nodes_distribution_colony.pdf' % FIGS_DIR, format='pdf')
    pylab.close()

def depth_distribution():
    def get_depth(node_name):
        depth = 0
        for char in node_name:
            if char in 'ABCDEF':
                depth += 1
        return depth

    fname = '%s/depth_distribution.csv' % CSV_DIR
    if not os.path.exists(fname):
        with open(fname, 'w') as f:
            f.write('colony, day, node, depth\n')
            depth_counts = defaultdict(int)
            unique_nodes = set()
            for colony in COLONIES:
                G, LG = read_network_paper(colony)
                for day in G.graph['days']:
                    D = day_subnetwork(G, day)
                    for u in D.nodes():
                        if u in G.graph['terminals']:
                            continue
                        pair = (colony, u)
                        if pair in unique_nodes:
                            continue
                        else:
                            depth = get_depth(u)
                            f.write('%s, %s, %s, %d\n' % (colony, day, u, depth))
                            depth_counts[depth] += 1
                            unique_nodes.add(pair)

    '''
    total = float(sum(depth_counts.values()))
    df_dict = {'depth' : [], 'proportion' : []}
    for depth, count in depth_counts.iteritems():
        df_dict['depth'].append(depth)
        df_dict['proportion'].append(count / total)

    df = pd.DataFrame.from_dict(df_dict)
    '''
    df = pd.read_csv(fname, skipinitialspace=True)
    df = make_proportions_df(df, 'depth')
    df['colony'] = df['colony'].map(COLONY_LABELS)

    pylab.figure()
    ax = sns.barplot(x='depth', y='proportion', hue='colony', data=df, hue_order=sorted(COLONY_LABELS.values()))
    ax.get_legend().remove()
    pylab.xlabel('node depth', size=20)
    pylab.ylabel('proportion of nodes\nused by ants', size=20)
    pylab.tick_params(axis='both', labelsize=20)
    pylab.ylim(0, 1)
    pylab.savefig('%s/depth_distribution.pdf' % FIGS_DIR, format='pdf', bbox_inches='tight')
    pylab.close()

def paper_plots():
    plot_connectivity_distribution()
    edge_length_distribution()
    repeatability_distribution()
    plot_multipath_sizes()
    plot_degree_distribution()
    used_nodes_distribution()
    depth_distribution()

def table1_stats():
    all_files = []

    vegetation_lengths_file = open('%s/full_vegetation_edge_lengths2018-19.csv' % PLOT_DATA_DIR, 'w')
    vegetation_lengths_file.write('colony, u, v, length\n')
    all_files.append(vegetation_lengths_file)

    day_lengths_file = open('%s/day_by_day_edge_lengths2018-19.csv' % PLOT_DATA_DIR, 'w')
    day_lengths_file.write('colony, day, u, v, length\n')
    all_files.append(day_lengths_file)

    vegetation_ti_file = open('%s/full_vegetation_transition_indices2018-19.csv' % PLOT_DATA_DIR, 'w')
    vegetation_ti_file.write('colony, x, y, z, transition index\n')
    all_files.append(vegetation_ti_file)

    day_ti_file = open('%s/day_by_day_transition_indices2018-19.csv' % PLOT_DATA_DIR, 'w')
    day_ti_file.write('colony, day, x, y, z, transition index\n')
    all_files.append(day_ti_file)

    vegetation_nodes_file = open('%s/full_vegetation_nodes2018-19.csv' % PLOT_DATA_DIR, 'w')
    vegetation_nodes_file.write('colony, nodes\n')
    all_files.append(vegetation_nodes_file)

    day_nodes_file = open('%s/day_by_day_nodes2018-19.csv' % PLOT_DATA_DIR, 'w')
    day_nodes_file.write('colony, day, nodes\n')
    all_files.append(day_nodes_file)

    vegetation_depths_file = open('%s/full_vegetation_used_node_depths2018-19.csv' % PLOT_DATA_DIR, 'w')
    vegetation_depths_file.write('colony, u, depth\n')
    all_files.append(vegetation_depths_file)

    connectivities_fname = '%s/day_by_day_local_connectivities2018-19.csv' % PLOT_DATA_DIR
    if not os.path.exists(connectivities_fname):
        write_connectivities_file(connectivities_fname)
    connectivities_df = pd.read_csv(connectivities_fname, skipinitialspace=True)

    for colony in COLONIES:
        print '\n--------------------'
        colony_str = COLONY_LABELS[colony]
        print colony_str
        G, LG = read_network_paper(colony)

        colony_lengths = []
        colony_transition_indices = []

        for u, v in G.edges():
            length = G[u][v]['length']
            colony_lengths.append(length)
            vegetation_lengths_file.write('%s, %s, %s, %d\n' % (colony_str, u, v, length))

        for u, v in LG.edges():
            a, b = u
            c, d = v
            assert G.has_edge(a, b)
            assert G.has_edge(c, d)
            if b != c:
                continue

            ti = LG[u][v]['transition index']
            colony_transition_indices.append(ti)
            vegetation_ti_file.write('%s, %s, %s, %s, %d\n' % (colony_str, a, b, d, ti))

        print "available \ntotal length:", sum(colony_lengths)
        print "available average transition index:", pylab.mean(colony_transition_indices)
        print "available total nodes:", G.number_of_nodes()
        vegetation_nodes_file.write('%s, %d\n' % (colony_str, G.number_of_nodes()))

        day_mean_lengths = []
        day_mean_transition_indices = []
        day_nodes = []
        day_mean_connectivities = []

        unique_nodes_used = set()
        for day in G.graph['days']:
            D = day_subnetwork(G, day)
            unique_nodes_used.update(list(D.nodes()))
            day_lengths = []
            for u, v in D.edges():
                length = G[u][v]['length']
                day_lengths.append(length)
                day_lengths_file.write('%s, %s, %s, %s, %d\n' % (colony_str, day, u, v, length))
            day_mean_lengths.append(pylab.mean(day_lengths))

            day_transition_indices = []
            LD = nx.line_graph(D)
            for u, v in LD.edges():
                if LG.has_edge(u, v):
                    a, b = u
                    c, d = v
                    assert b == c

                    ti = LG[u][v]['transition index']
                    day_transition_indices.append(LG[u][v]['transition index'])
                    day_ti_file.write('%s, %s, %s, %s, %s, %d\n' % (colony_str, day, a, b, c, ti))
            day_mean_transition_indices.append(pylab.mean(day_transition_indices))

            day_nodes.append(D.number_of_nodes())
            day_nodes_file.write('%s, %s, %d\n' % (colony_str, day, D.number_of_nodes()))

        colony_connectivities_df = connectivities_df[connectivities_df['colony'] == colony_str]
        for day, day_group in colony_connectivities_df.groupby('day'):
            day_mean_connectivities.append(pylab.mean(day_group['connectivity']))

        for n in unique_nodes_used:
            depth = get_depth(n)
            vegetation_depths_file.write('%s, %s, %d\n' % (colony_str, n, depth))

        print "\nused average lengths:", pylab.mean(day_mean_lengths), '+/-', pylab.std(day_mean_lengths, ddof=1)
        print "used total nodes:", pylab.mean(day_nodes), "+/-", pylab.std(day_nodes, ddof=1)
        print "used average transition indices:", pylab.mean(day_mean_transition_indices), "+/-", pylab.std(day_mean_transition_indices, ddof=1)
        print "used average connectivities:", pylab.mean(day_mean_connectivities), "+/-", pylab.std(day_mean_connectivities, ddof=1)

        print "--------------------"

    for f in all_files:
        f.close()

def connecting_nodes(D, terminal1, terminal2):
    connecting_nodes = set()
    for path in nx.all_simple_paths(D, terminal1, terminal2):
        connecting_nodes.update(path[1:-1])
    return connecting_nodes

def nodes_to_terminal(D, terminal):
    anchor = None
    if D.has_node('nest1'):
        anchor = 'nest1'
    else:
        assert 'pseudo terminals' in D.graph
        anchor = D.graph['pseudo terminals'][0]

    return connecting_nodes(D, anchor, terminal)

def conserved_terminals(prev_day, D):
    prev_terminals = prev_day.graph['terminals used']
    terminals = D.graph['terminals used']

    conserved = set(prev_terminals) & set(terminals)
    conserved.discard('nest1')
    return conserved

def added_terminals(prev_day, D):
    prev_terminals = prev_day.graph['terminals used']
    terminals = D.graph['terminals used']

    added = set(terminals) - set(prev_terminals)
    added.discard('nest1')
    return added

def node_changes(prev_day, D):
    prev_nodes = set(prev_day.nodes())
    curr_nodes = set(D.nodes())

    added_nodes = curr_nodes - prev_nodes
    discarded_nodes = prev_nodes - curr_nodes

    conserved = conserved_terminals(prev_day, D)
    conservation_nodes = set()
    for terminal in conserved:
        conservation_nodes.update(nodes_to_terminal(D, terminal))
        conservation_nodes.update(nodes_to_terminal(prev_day, terminal))

    modification_nodes = added_nodes & conservation_nodes

    discarded_nodes = discarded_nodes & conservation_nodes

    new_path_nodes = added_nodes - conservation_nodes

    return map(len, [modification_nodes, discarded_nodes, new_path_nodes])

def supp_table1_stats():
    for colony in COLONIES:
        prev_network = None

        print colony
        print '-' * len(colony)
        G, LG = read_network_paper(colony)
        for i, day in enumerate(G.graph['days']):
            print '*' * len(day)
            print day
            print '*' * len(day)
            D = day_subnetwork(G, day)

            if prev_network != None:
                modification_nodes, discarded_nodes, new_path_nodes = node_changes(prev_network, D)

                print 'day %d | %d | %d | %d' % (i + 1, modification_nodes, discarded_nodes, new_path_nodes)

            prev_network = D

def test_new_function():
    depth_distribution()

def plot_connectivities():
    plot_connectivity_distribution()
    plot_connectivity_distribution(transition=True)
    #plot_multipath_sizes()
    #plot_multipath_repeatabilities()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--paper', action='store_true')
    parser.add_argument('-t', '--test', action='store_true')
    parser.add_argument('-c', '--connectivity', action='store_true')
    parser.add_argument('--table1', action='store_true')
    parser.add_argument('--supptable1', action='store_true')
    parser.add_argument('--correlations', action='store_true')

    args = parser.parse_args()

    if args.paper:
        paper_plots()

    if args.test:
        test_new_function()

    if args.connectivity:
        plot_connectivities()

    if args.table1:
        table1_stats()

    if args.supptable1:
        supp_table1_stats()

    if args.correlations:
        objective_correlations()

if __name__ == '__main__':
    main()
