import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import pylab
from map_network import *
import argparse
import os
from collections import defaultdict
import seaborn as sns
sns.set_style('white')
from random import uniform
from steiner_trees import *
from random_trees import *
from scipy.stats import f_oneway, spearmanr, sem, pearsonr, ttest_ind, percentileofscore, mannwhitneyu, kruskal, wilcoxon, friedmanchisquare
from statsmodels.stats.anova import anova_lm
from itertools import combinations, product
from sys import exit
from scikit_posthocs import posthoc_conover, posthoc_wilcoxon
#from statsmodels.formula.api import ols

MAX_TRANSITION_INDEX = 4
DEFAULT_NULL_MODEL = 'progression'
DEFAULT_ALGORITHM = 'ants'

UNITS = {'average transition index' : 'index', 'length' : 'cm',\
         'total repeatability' : 'index', 'total nodes' : 'nodes',\
         'total length' : 'cm'}

FILENAMES = {'average transition index' : 'transition', 'length' : 'length',\
             'total transition index' : 'total_transition',\
             'total length' : 'total_length', 'nodes' : 'nodes',}

LABEL_STR = {'average transition index' : 'Average transition index',\
             'total transition index' : 'Total transition index',\
             'average length' : 'Average edge length',\
             'total length' : 'Total edge length',\
             'total nodes' : 'Total nodes'}

OBJECTIVE_LABELS = {'total nodes' : 'Total nodes',\
                    'average transition index' : 'Average transition index',\
                    'average length' : 'Average edge length',\
                    'total length' : 'Total edge length'}

OBJECTIVE_ZORDERS = {'total nodes' : 3, 'average transition index' : 2,\
                     'average length' : 1, 'total transition index' : 4,\
                     'total length' : 5}

OBJECTIVE_COLORS = {'average transition index' : '#800080',\
                    'average length' : '#996515', 'total nodes' : 'c',\
                    'total length' : 'm', 'total transition index' : 'g'}

NULL_MODEL_NAMES = {'progression' : 'random progression subgraph',\
                    'independent' : 'random independent subgraph'}

NULL_MODEL_FUNCS = {'progression' : random_steiner_progression,\
                    'independent' : random_steiner_independent}

ALGORITHM_LABELS = {'ants' : 'Observed', 'heuristic' : 'Optimized'}

ANOVA_OBJECTIVES = ['total nodes', 'average transition index', 'average length']

HIST_PLOT_OBJECTIVES = ['total nodes', 'average transition index', 'average length', 'total length']

FIGS_DIR = 'objective_percentiles'

PAPER_DIR = 'paper_data/plot_data'
PAPER_FNAME = 'objective_percentiles_paper2018-19'
PAPER_VALUES_FNAME = 'objective_values_paper2018-19'

INCLUDE_LEGEND = False
FONT_FAMILY = 'Arial'
LEGEND_BORDER_AXES_PAD = 0.5
LEGEND_BORDER_PAD = 0.5
LEGEND_LABEL_SPACING = 0.5

BROKEN_COLONY_SHARE_X = False

def set_objectives(G):
    for u, v in G.edges():
        length = G[u][v]['length']
        repeatability = G[u][v]['repeatability']
        G[u][v]['nodes'] = 1
        G[u][v]['length*repeatability'] = length * repeatability
        G[u][v]['nodes*repeatability'] = repeatability
        G[u][v]['nodes*length'] = length
        G[u][v]['random'] = uniform(0, 1)

def remove_dead_ends(G, terminals=None):
    if terminals == None:
        terminals = []
    done = False
    while not done:
        bad_nodes = []
        for u in G.nodes():
            if G.degree(u) <= 1 and u not in terminals:
                bad_nodes.append(u)
        if len(bad_nodes) > 0:
            for bad_node in bad_nodes:
                assert G.degree(bad_node) <= 1
                assert bad_node not in terminals
                G.remove_node(bad_node)
        else:
            done = True

    for u in G.nodes():
        assert G.degree(u) > 1 or u in terminals

def total_length(G, LG, T, **kwargs):
    total = 0.0
    for u, v in T.edges():
        total += G[u][v]['length']
    return total

def mean_length(G, LG, T, **kwargs):
    return total_length(G, LG, T) / T.number_of_edges()

def total_nodes(G, LG, T, **kwargs):
    return float(T.number_of_nodes())

def get_transitions(G, LG, T, max_ti=MAX_TRANSITION_INDEX, **kwargs):
    transitions = []
    for u in T.nodes():
        for p in T.predecessors(u):
            for s in T.successors(u):
                x = (p, u)
                y = (u, s)
                if LG.has_edge(x, y):
                    transition_index = LG[x][y]['transition index']
                    if transition_index <= max_ti:
                        transitions.append(transition_index)

    return transitions

def mean_transition_index(G, LG, T, **kwargs):
    max_ti = None
    if 'max_ti' in kwargs:
        max_ti = kwargs['max_ti']
    else:
        max_ti = MAX_TRANSITION_INDEX

    transitions = get_transitions(G, LG, T, **kwargs)
    assert len(transitions) > 0
    return sum(transitions) / float(len(transitions))

def total_transition_index(G, LG, T, **kwargs):
    max_ti = None
    if 'max_ti' in kwargs:
        max_ti = kwargs['max_ti']
    else:
        max_ti = MAX_TRANSITION_INDEX

    transitions = get_transitions(G, LG, T, **kwargs)
    assert len(transitions) > 0
    return sum(transitions)

def write_objectives(random_trials=100, max_ti=MAX_TRANSITION_INDEX,\
                     null_model=DEFAULT_NULL_MODEL):
    values_fname = 'objective_values'
    values_fname += '.csv'

    first_time = not os.path.exists(values_fname)

    objectives = {'average length' : mean_length, 'total length' : total_length,\
                  'total nodes' : total_nodes,\
                  'average transition index' : mean_transition_index,\
                  'total transition index' : total_transition_index}

    heuristic_funcs = {'average length' : opt_average_length,\
                     'total length' : opt_total_length,\
                     'total nodes' : opt_total_nodes,\
                     'average transition index' : opt_average_transition_index,\
                     'total transition index' : opt_total_transition_index}

    null_model_name = NULL_MODEL_NAMES[null_model]
    null_model_func = NULL_MODEL_FUNCS[null_model]

    with open(values_fname, 'a') as f:

        if first_time:
            f.write('colony, day, max transition index, objective, model, value\n')

        for colony in COLONIES:
            print '*****', colony, '*****'
            G, LG = read_network_paper(colony)

            days = G.graph['days']

            if first_time:
                for day in days:
                    D = day_subnetwork(G, day)
                    for objective, objective_func in objectives.iteritems():
                        ant_value = objective_func(G, LG, D, max_ti=max_ti)
                        f.write('%s, %s, %d, %s, %s, %f\n' % (colony, day, max_ti, objective,\
                                                              'ants', ant_value))

                        heuristic_func = heuristic_funcs[objective]
                        heuristic_network = heuristic_func(G, LG, D.graph['terminals used'][:])
                        opt_value = objective_func(G, LG, heuristic_network, max_ti=max_ti)
                        f.write('%s, %s, %d, %s, %s, %f\n' % (colony, day, max_ti, objective,\
                                                              'heuristic', opt_value))

            for i in xrange(random_trials):
                null_trees = null_model_func(G)
                assert len(null_trees) == len(days)

                for day, null_tree in zip(days, null_trees):
                    for objective, objective_func in objectives.iteritems():
                        null_value = objective_func(G, LG, null_tree, max_ti=max_ti)
                        f.write('%s, %s, %d, %s, %s, %f\n' % (colony, day, max_ti,\
                                                              objective,\
                                                              null_model_name, null_value))

    write_percentiles()

def write_percentiles():
    values_fname = 'objective_values'
    values_fname += '.csv'
    df = pd.read_csv(values_fname, skipinitialspace=True)
    #df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df_dict = {'colony' : [], 'day' : [], 'max transition index' : [], 'objective' : [],\
            'null model' : [], 'algorithm' : [], 'value' : [], 'percentile' : []}

    algorithms = ['ants', 'heuristic']
    for name, group in df.groupby(['colony', 'day', 'max transition index', 'objective']):
        colony, day, max_ti, objective = name

        for model, model_group in group.groupby('model'):
            if model in algorithms:
                continue

            for algorithm in algorithms:
                alg_group = group[group['model'] == algorithm]
                alg_val = alg_group['value'].mean()

                null_vals = model_group['value']
                percentile = percentileofscore(null_vals, alg_val)

                df_dict['colony'].append(colony)
                df_dict['day'].append(day)
                df_dict['max transition index'].append(max_ti)
                df_dict['objective'].append(objective)
                df_dict['null model'].append(model)
                df_dict['algorithm'].append(algorithm)
                df_dict['value'].append(alg_val)
                df_dict['percentile'].append(percentile)

    anova_df = pd.DataFrame.from_dict(df_dict)
    percentiles_fname = values_fname.replace('values', 'percentiles')
    anova_df.to_csv(percentiles_fname, index=False)

def plot_changes(differences=False, paper=False, max_ti=MAX_TRANSITION_INDEX,\
                 null_model=DEFAULT_NULL_MODEL, algorithm=DEFAULT_ALGORITHM):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df = df[df['max transition index'] == max_ti]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['algorithm'] == algorithm]

    def date_key(date):
        month, day, year = date.split('/')
        etc = year[4:]
        year = year[:4]
        return (year, month, day, date[1], etc)

    df['day'] = df['day'].map(reformat_day)
    broken_axis_colonies = ['turtle_hill', 'tejon7', 'T189']
    broken_axis_df = df[df['colony'].isin(broken_axis_colonies)]
    broken_axis_days = sorted(list(broken_axis_df['day'].unique()), key = date_key)

    y_axis_colonies = ['T189', 'turtle_hill']
    x_axis_colonies = ['turtle_hill', 'T500']
    legend_colonies = ['tejon7']

    labelsize = 50
    ticklabelsize = 50
    labelpad = 75
    ticklabelpadding = 5
    xlabel_rotation = 0
    xlim_padding = 0.35


    outdir = '%s/%s' % (FIGS_DIR, null_model)
    os.system('mkdir -p %s' % outdir)
    for colony, colony_group in df.groupby('colony'):
        print colony
        broken_colony = colony in broken_axis_colonies
        y_axis_colony = colony in y_axis_colonies
        x_axis_colony = colony in x_axis_colonies
        legend_colony = colony in legend_colonies

        G, LG = read_network_paper(colony)
        days = G.graph['days']

        day_ruptures = defaultdict(int)
        for node, day in G.graph['node ruptures']:
            day_ruptures[reformat_day(day)] += 1

        for u, v, day in G.graph['edge ruptures']:
            day_ruptures[reformat_day(day)] += 1

        days_list = None
        if broken_colony and BROKEN_COLONY_SHARE_X:
            days_list = broken_axis_days
        else:
            days_list = sorted(list(colony_group['day'].unique()), key=date_key)

        day_index_cutoff = 0
        if broken_colony:
            day_index_cutoff += 1
        num_left_days = len(days_list) - day_index_cutoff
        num_right_days = day_index_cutoff

        pylab.rcParams['font.family'] = FONT_FAMILY
        fig, all_axes = None, None

        default_width = 6.4
        default_height = 4.8
        fig_width = 20
        fig_height= 10
        area_scale = (fig_width * fig_height) / (default_width * default_height)

        if colony in broken_axis_colonies:
            fig, (ax, ax2) = pylab.subplots(1, 2, sharey=True,\
                                            gridspec_kw = {'width_ratios': [num_left_days, num_right_days]},\
                                            figsize=(fig_width, fig_height))
            all_axes = [ax, ax2]
        else:
            fig, ax = pylab.subplots(1, 1, figsize=(fig_width, fig_height))
            all_axes = [ax]

        x = None
        percentile = None
        labels = []

        for objective, objective_group in colony_group.groupby('objective'):
            if objective not in ANOVA_OBJECTIVES:
                continue

            x = []
            labels = []
            percentile = []
            for day, day_group in objective_group.groupby('day'):
                idx = days_list.index(day) + 1
                x.append(idx)
                labels.append(str(idx))
                percentile.append(day_group['percentile'].mean())

            x = pylab.array(x)
            percentile = pylab.array(percentile)
            order = pylab.argsort(x)
            x = x[order]
            percentile = percentile[order]

            if differences:
                x = x[1:]
                labels = labels[1:]
                percentile = pylab.diff(percentile)

            objective_label = OBJECTIVE_LABELS[objective]
            color = OBJECTIVE_COLORS[objective]
            zorder = OBJECTIVE_ZORDERS[objective]
            width = 10
            if colony in ['T189', 'T500'] and objective == 'average transition index':
                width += 10

            for axis in all_axes:
                axis.plot(x, percentile, label='_nolegend_', zorder=zorder, linewidth=width)
                axis.scatter(x, percentile, s=750, label=objective_label, zorder=zorder)
                axis.axhline(y=50, c='k', linestyle='dashed', linewidth=4, zorder=0.5)

        for i, axis in enumerate(all_axes):
            N = len(days_list)

            tickmarks = range(1, N + 1)
            axis.set_xticks(tickmarks)
            axis.set_xticklabels(map(str, tickmarks))

            axis.tick_params(axis='both', labelsize=ticklabelsize)
            axis.tick_params(axis='both', which='major', pad=ticklabelpadding)

            xtick_color = 'k' if x_axis_colony else 'k'
            axis.tick_params(axis='x', rotation=-xlabel_rotation, colors=xtick_color)

            ytick_color = 'k' if y_axis_colony else 'w'
            axis.tick_params(axis='y', colors=ytick_color)

            for spine in ['right', 'top']:
                axis.spines[spine].set_visible(False)
            if i == 1:
                axis.spines['left'].set_visible(False)

            axis.spines['bottom'].set_linewidth(4)
            if i == 0:
                axis.spines['left'].set_linewidth(4)

            if i == 0:
                axis.yaxis.tick_left()

            if i == 1:
                axis.yaxis.tick_right()

            left_lim = None
            right_lim = None
            if i == 0:
                left_lim = 1 - xlim_padding
                right_lim = num_left_days +  xlim_padding
            else:
                left_lim = num_left_days + 1 - xlim_padding
                right_lim = num_left_days + num_right_days + xlim_padding

            axis.set_xlim(left=left_lim, right=right_lim)

        ymin = None
        if differences:
            ymin = -115
        else:
            ymin = -15

        ymax = None
        if INCLUDE_LEGEND:
            ymax = 200
        else:
            ymax = 105

        if ymax == None:
            pylab.ylim(bottom=ymin)
        else:
            pylab.ylim(bottom=ymin, top=ymax)

        ytick_color = 'k' if y_axis_colony else 'w'
        pylab.yticks([0, 50, 100], color=ytick_color)

        if len(all_axes) > 1:
            ax, ax2 = all_axes[0], all_axes[1]
            d = 0.010  # how big to make the diagonal lines in axes coordinates
            # arguments to pass to plot, just so we don't keep repeating them
            kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
            scale_factor = float(num_right_days + 1) / float(len(broken_axis_days))
            ax.plot((1 - d * scale_factor, 1 + d * scale_factor), (-d, d), zorder=float('inf'), **kwargs)        # top-left diagonal
            #ax.plot((1 - d * scale_factor, 1 + d * scale_factor), (1 - d, 1 + d), zorder=float('inf'), **kwargs)  # top-right diagonal

            kwargs.update(transform=ax2.transAxes)  # switch to the right axes
            ax2.plot((-d, d), (-d, d), zorder=float('inf'), **kwargs)  # bottom-left diagonal

        ax = all_axes[0]
        legend_scale = 2.5

        if INCLUDE_LEGEND:
            leg = ax.legend(frameon=True, loc=2, borderaxespad=LEGEND_BORDER_AXES_PAD * legend_scale,\
                                                 borderpad=LEGEND_BORDER_PAD * legend_scale,\
                                                 labelspacing=LEGEND_LABEL_SPACING * legend_scale)


            pylab.setp(ax.get_legend().get_texts(), fontsize=20 * legend_scale) # for legend text
            leg.get_frame().set_linewidth(4)
            leg.get_frame().set_edgecolor('k')
        elif legend_colony:
            leg = ax.legend(frameon=True)
            pylab.setp(ax.get_legend().get_texts(), fontsize=labelsize) # for legend text
            leg.get_frame().set_linewidth(4)
            leg.get_frame().set_edgecolor('k')

        final_axis = fig.add_subplot(111, frameon=False)
        final_axis.set_xticks([])
        final_axis.set_yticks([])

        xlabel_color = 'k' if x_axis_colony else 'w'
        final_axis.set_xlabel('Observation day', size=labelsize, labelpad=labelpad, color=xlabel_color)
        ylab = None
        if differences:
            ylab = 'difference in percentile of ants\ncompared to previous day'
        else:
            ylab = 'Comparison to\nrandom networks (%)'

        ylabel_color = 'k' if y_axis_colony else 'w'
        final_axis.set_ylabel(ylab, size=labelsize, labelpad=labelpad + 10, color=ylabel_color)

        # code to add annotations to the y-axis
        if y_axis_colony:
            yaxis_annotation_color = 'k' if y_axis_colony else 'w'
            final_axis.annotate('Equal', xy=(-0.3, 0.515), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

            final_axis.annotate('', xy=(-0.233, 0.17), xytext=(-0.233, 0.47),\
                                arrowprops={'color' : yaxis_annotation_color},\
                                annotation_clip=False,\
                                fontsize=labelsize)
            final_axis.annotate('Better', xy=(-0.3, 0.1), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

            final_axis.annotate('', xy=(-0.233, 0.9), xytext=(-0.233, 0.6),\
                                arrowprops={'color' : yaxis_annotation_color},\
                                annotation_clip=False,\
                                fontsize=labelsize)
            final_axis.annotate('Worse', xy=(-0.3, 0.93), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)


        fname = 'objective_'
        if differences:
            fname += 'differences_'
        else:
            fname += 'changes_'
        fname += colony + '.pdf'

        pylab.tight_layout()
        pylab.savefig('%s/%s' % (outdir, fname), format='pdf', bbox_inches='tight')
        pylab.close()


def plot_percentiles_distribution(paper=False, max_ti=MAX_TRANSITION_INDEX, \
                                  null_model=DEFAULT_NULL_MODEL,\
                                  algorithm=DEFAULT_ALGORITHM):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    percentile_df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    percentile_df = percentile_df[percentile_df['objective'].isin(HIST_PLOT_OBJECTIVES)]
    percentile_df = percentile_df[percentile_df['max transition index'] == max_ti]
    percentile_df = percentile_df[percentile_df['null model'] == NULL_MODEL_NAMES[null_model]]
    #percentile_df = percentile_df[percentile_df['algorithm'] == algorithm]

    percentile_df['objective'] = percentile_df['objective'].map(LABEL_STR)

    def reformat_objective(objective):
        objective = '\n'.join(objective.split())
        return objective

    def reformat_algorithm(algorithm):
        return ALGORITHM_LABELS[algorithm]

    percentile_df['colony'] = percentile_df['colony'].map(COLONY_LABELS)
    percentile_df['objective'] = percentile_df['objective'].apply(reformat_objective)
    percentile_df['algorithm'] = percentile_df['algorithm'].apply(reformat_algorithm)
    pylab.figure(figsize=(8, 4.8))
    pylab.rcParams['font.family'] = FONT_FAMILY

    ax = sns.barplot(x='objective', y='percentile', hue='algorithm', data=percentile_df,\
                     order=['Average\nedge\nlength', 'Average\ntransition\nindex', 'Total\nnodes', 'Total\nedge\nlength'])
    pylab.ylabel('Comparison to\nrandom networks (%)', size=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.xaxis.label.set_visible(False)
    for spine in ['left', 'right', 'top', 'bottom']:
        pass #ax.spines[spine].set_linewidth(4)

    pylab.yticks([0, 50, 100], fontsize=20)
    ylim_padding = 0.5
    pylab.ylim(-ylim_padding, 100 + ylim_padding + 10)

    pylab.axhline(y=50, lw=4, c='k', ls='dashed', zorder=0.5)

    leg = pylab.legend(frameon=True, loc="upper right", borderaxespad=LEGEND_BORDER_AXES_PAD,\
                                            borderpad=LEGEND_BORDER_PAD,\
                                            labelspacing=LEGEND_LABEL_SPACING)
    pylab.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
    #leg.get_frame().set_linewidth(4)
    leg.get_frame().set_edgecolor('k')

    final_axis = pylab.gca()
    labelsize = 20
    yaxis_annotation_color = 'k'
    final_axis.annotate('Equal', xy=(-1.8, 48.5), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

    final_axis.annotate('', xy=(-1.54, 7), xytext=(-1.54, 45),\
                        arrowprops={'color' : yaxis_annotation_color, 'width' : 1.25, 'headlength' : 4, 'headwidth' : 4},\
                        annotation_clip=False,\
                        fontsize=labelsize)
    final_axis.annotate('Better', xy=(-1.8, -2), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

    final_axis.annotate('', xy=(-1.54, 93), xytext=(-1.54, 55),\
                        arrowprops={'color' : yaxis_annotation_color, 'width' : 1.25, 'headlength' : 4, 'headwidth' : 4},\
                        annotation_clip=False,\
                        fontsize=labelsize)
    final_axis.annotate('Worse', xy=(-1.8, 97), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)


    y_diff = 3

    asterisk_y_low = 67
    asterisk_y_high = asterisk_y_low + y_diff

    ns_y_low = 53
    ns_y_high = ns_y_low + y_diff

    x_diff = 1
    x_pad = 0.05
    x_min = -0.2
    ns_x_min = x_min + x_diff

    asterisk_text_x_offset = 0.07
    asterisk_text_y_offset = 3

    ns_text_x_offset = 0.05
    ns_text_y_offset = -5


    final_axis.plot([x_min] * 2 + [x_min + x_diff] * 3 + [x_min + 2 * x_diff] * 2,\
                    [asterisk_y_low] + [asterisk_y_high] * 2 +  [asterisk_y_low] + [asterisk_y_high] * 2 + [asterisk_y_low],\
                    c='k')
    final_axis.annotate('*', (x_min + 2 * x_diff / 2.0 - asterisk_text_x_offset,
                               asterisk_y_low - asterisk_text_y_offset),\
                               fontsize=40)

    final_axis.plot([ns_x_min] * 2 + [ns_x_min + x_diff] * 2,
                    [ns_y_low] + [ns_y_high] * 2 +  [ns_y_low],\
                    c='k')
    final_axis.annotate('ns', ( ns_x_min + x_diff / 2.0 - ns_text_x_offset,
                                ns_y_low - ns_text_y_offset),\
                                fontsize=20)

    # final_axis.plot([-0.20, -0.20, 1.8,  1.8], [80, 83, 83, 80], c='k')
#     final_axis.annotate('*', (0.74, 77), fontsize=40)

    pylab.tight_layout()

    outdir = '%s/%s' % (FIGS_DIR, null_model)
    os.system('mkdir -p %s' % outdir)

    pylab.savefig('%s/objective_percentiles_distribution.pdf' % (outdir), format='pdf', bbox_inches='tight')
    pylab.close()

def plot_percentiles_swarm(paper=False, max_ti=MAX_TRANSITION_INDEX, \
                                  null_model=DEFAULT_NULL_MODEL,\
                                  algorithm=DEFAULT_ALGORITHM):
    values_fname = None
    if paper:
        values_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_VALUES_FNAME)
    else:
        values_fname = 'objective_percentiles.csv'

    values_df = pd.read_csv(values_fname, skipinitialspace=True)
    values_df = values_df[values_df['objective'].isin(HIST_PLOT_OBJECTIVES)]
    values_df = values_df[values_df['max transition index'] == max_ti]
    values_df = values_df[values_df['model'].isin(['ants', NULL_MODEL_NAMES[null_model]])]

    values_df['objective'] = values_df['objective'].map(LABEL_STR)

    def reformat_objective(objective):
        objective = '\n'.join(objective.split())
        return objective

    def reformat_algorithm(algorithm):
        return ALGORITHM_LABELS[algorithm]

    values_df['colony'] = values_df['colony'].map(COLONY_LABELS)
    values_df['objective'] = values_df['objective'].apply(reformat_objective)

    percentiles_df = pd.DataFrame()
    for objective, objective_group in values_df.groupby('objective'):
        for colony, colony_group in objective_group.groupby('colony'):
            for day, day_group in colony_group.groupby('day'):
                random_df = day_group[day_group['model'] == NULL_MODEL_NAMES[null_model]].copy()
                all_values = day_group['value']
                random_df = random_df.sample(n=25)
                def get_percentile(x):
                    return percentileofscore(all_values, x)

                ant_df = day_group[day_group['model'] == 'ants'].copy()
                ant_df['percentile'] = ant_df['value'].apply(get_percentile)


                random_df['percentile'] = random_df['value'].apply(get_percentile)

                percentiles_df = pd.concat([percentiles_df, ant_df, random_df])


    random_percentiles = percentiles_df[percentiles_df['model'] == NULL_MODEL_NAMES[null_model]]
    ant_percentiles = percentiles_df[percentiles_df['model'] == 'ants']

    pylab.figure(figsize=(11.5, 4.8))
    pylab.rcParams['font.family'] = FONT_FAMILY

    ax = sns.swarmplot(x='objective', y='percentile', data=random_percentiles,\
                       size=5, color='k', facecolors=None, alpha=0.1, label='Random',\
                       order=['Average\nedge\nlength', 'Average\ntransition\nindex', 'Total\nnodes', 'Total\nedge\nlength'])
    sns.swarmplot(x='objective', y='percentile', hue='colony', data=ant_percentiles, ax=ax,\
                  size=10, marker="X", dodge=True,\
                  order=['Average\nedge\nlength', 'Average\ntransition\nindex', 'Total\nnodes', 'Total\nedge\nlength'])

    pylab.ylabel('Comparison to\nrandom networks (%)', size=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.xaxis.label.set_visible(False)
    for spine in ['left', 'right', 'top', 'bottom']:
        pass #ax.spines[spine].set_linewidth(4)

    pylab.yticks([0, 50, 100], fontsize=20)
    ylim_padding = 0.5
    pylab.ylim(-ylim_padding, 100 + ylim_padding + 30)

    pylab.axhline(y=50, lw=4, c='k', ls='dashed', zorder=0.5)


    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[3:], labels[3:], frameon=True, loc="upper right",
                                                bbox_to_anchor=(1.35, 1),
                                                borderaxespad=LEGEND_BORDER_AXES_PAD,\
                                                borderpad=LEGEND_BORDER_PAD,\
                                                labelspacing=LEGEND_LABEL_SPACING,
                                                ncol=1)

    for handle in leg.legendHandles:
        handle._sizes = [100]

    pylab.setp(leg.get_texts(), fontsize=20) # for legend text
    #leg.get_frame().set_linewidth(4)
    leg.get_frame().set_edgecolor('k')

    final_axis = pylab.gca()
    labelsize = 20
    yaxis_annotation_color = 'k'
    final_axis.annotate('Equal', xy=(-1.75, 48.5), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

    final_axis.annotate('', xy=(-1.54, 7), xytext=(-1.54, 45),\
                        arrowprops={'color' : yaxis_annotation_color, 'width' : 1.25, 'headlength' : 4, 'headwidth' : 4},\
                        annotation_clip=False,\
                        fontsize=labelsize)
    final_axis.annotate('Better', xy=(-1.75, -2), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

    final_axis.annotate('', xy=(-1.54, 93), xytext=(-1.54, 55),\
                        arrowprops={'color' : yaxis_annotation_color, 'width' : 1.25, 'headlength' : 4, 'headwidth' : 4},\
                        annotation_clip=False,\
                        fontsize=labelsize)
    final_axis.annotate('Worse', xy=(-1.75, 97), annotation_clip=False, fontsize=labelsize, color=yaxis_annotation_color)

    y_diff = 3

    asterisk_y_low = 116
    asterisk_y_high = asterisk_y_low + y_diff

    ns_y_low = 103
    ns_y_high = ns_y_low + y_diff

    x_diff = 1
    x_pad = 0.05
    x_min = 0
    ns_x_min = x_min + x_diff

    asterisk_text_x_offset = 0.07
    asterisk_text_y_offset = 3

    ns_text_x_offset = 0.05
    ns_text_y_offset = -5

#     final_axis.plot([0, 0, 0.95,  0.95], [103, 105, 105, 103], c='k')
#     final_axis.annotate('*', (0.45, 98), fontsize=40)

    final_axis.plot([ns_x_min, ns_x_min, ns_x_min + x_diff,  ns_x_min + x_diff],\
                    [ns_y_low, ns_y_high, ns_y_high, ns_y_low],\
                    c='k')
    final_axis.annotate('ns', (ns_x_min + x_diff / 2.0 - ns_text_x_offset,\
                               ns_y_low - ns_text_y_offset), fontsize=20)

    final_axis.plot([x_min] * 2 + [x_min + x_diff] * 3 + [x_min + 2 * x_diff] * 2,\
                    [asterisk_y_low] + [asterisk_y_high] * 2 + [asterisk_y_low] + [asterisk_y_high] * 2 + [asterisk_y_low],\
                    c='k')
    final_axis.annotate('*', (x_min + 2 * x_diff / 2.0 - asterisk_text_x_offset,\
                              asterisk_y_low - asterisk_text_y_offset), fontsize=40)

    pylab.tight_layout()

    outdir = '%s/%s' % (FIGS_DIR, null_model)
    os.system('mkdir -p %s' % outdir)

    pylab.savefig('%s/objective_percentiles_swarmplot.pdf' % (outdir), format='pdf', bbox_inches='tight')
    pylab.close()

def test_percentiles(barely=False, max3=False, paper=False):
    percentiles_fname = 'objective_percentiles'
    if barely:
        percentiles_fname += '_barely'
    if max3:
        percentiles_fname += '_max3'
    if paper:
        percentiles_fname += '_paper'
    percentiles_fname += '.csv'
    anova_df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    anova_df = anova_df[anova_df['objective'].isin(ANOVA_OBJECTIVES)]
    for objective1, objective2 in combinations(ANOVA_OBJECTIVES, 2):
        df1 = anova_df[anova_df['objective'] == objective1]
        df2 = anova_df[anova_df['objective'] == objective2]

        percentiles1 = df1['percentile']
        percentiles2 = df2['percentile']

        print objective1, "vs.", objective2
        print mannwhitneyu(percentiles1, percentiles2)

    #sns.set()
    pylab.figure()
    ax = sns.violinplot(x='objective', y='percentile', hue='colony', data=anova_df)
    pylab.yticks(fontsize=20)
    pylab.ylabel('percentile', size=20)
    pylab.savefig('%s/objective_percentiles_distribution.pdf' % FIGS_DIR, format='pdf')
    pylab.close()
    print anova_df[['colony', 'objective', 'percentile', 'value']].groupby(['objective', 'colony']).agg([pylab.mean, sem])
    formula = 'percentile ~ C(objective) + C(colony) + C(objective):C(colony)'
    model = ols(formula, anova_df).fit()
    aov_table = anova_lm(model, typ=2)
    print aov_table

def one_way_anova(barely=False, max3=False, paper=False):
    percentiles_fname = 'objective_percentiles'
    if barely:
        percentiles_fname += '_barely'
    if max3:
        percentiles_fname += '_max3'
    if paper:
        percentiles_fname += '_paper'
    percentiles_fname += '.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df['day'] = df['day'].map(reformat_day)

    def anova(dframe, factor):
        groups = []
        for factor_val, factor_group in dframe.groupby(factor):
            groups.append(factor_group['percentile'])

        return kruskal(*groups)

    for factor in ['objective', 'colony']:
        print factor
        print "standard anova"
        print "--------------"
        test = anova(df, factor)
        print "H = %.15f, p = %.15f" % (test.statistic, test.pvalue)

        if not max3:
            print "\nleave one out"
            print "--------------"
            h_min = float('inf')
            h_max = float('-inf')

            p_min = float('inf')
            p_max = float('-inf')

            for name, group in df.groupby(['colony', 'day']):
                colony, day = name
                df2 = df[((df['colony'] != colony) | (df['day'] != day))]
                test = anova(df2, factor)

                h, p = test.statistic, test.pvalue

                h_min = min(h_min, h)
                h_max = max(h_max, h)

                p_min = min(p_min, p)
                p_max = max(p_max, p)

            print "%.15f < H < %.15f" % (h_min, h_max)
            print "%.15f < p < %.15f" % (p_min, p_max)

            print "\nall but august"
            print "----------------"
            df3 = df[df['day'] != '08/03']
            test = anova(df3, factor)
            print "H = %.15f, p = %.15f" % (test.statistic, test.pvalue)


def compare_heuristic(paper=False, max_ti=MAX_TRANSITION_INDEX,\
                      null_model=DEFAULT_NULL_MODEL,):
    fname = None
    if paper:
        fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        fname = 'objective_percentiles.csv'

    df = pd.read_csv(fname, skipinitialspace=True)
    df = df[df['objective'].isin(HIST_PLOT_OBJECTIVES)]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['max transition index'] == max_ti]

    for name, group in df.groupby('objective'):
        print "____________"
        print name
        print "____________"

        ant_percentiles = []
        heuristic_percentiles = []
        for name2, group2 in group.groupby(['colony', 'day']):
        #for name2, group2 in group.groupby(['colony']):
            percentiles = group2['percentile']
            ant_percentile = percentiles[group2['algorithm'] == 'ants'].mean()
            heuristic_percentile = percentiles[group2['algorithm'] == 'heuristic'].mean()

            ant_percentiles.append(ant_percentile)
            heuristic_percentiles.append(heuristic_percentile)

        ant_percentiles = pylab.array(ant_percentiles)
        heuristic_percentiles = pylab.array(heuristic_percentiles)

        print ant_percentiles
        print heuristic_percentiles

        print "ant percentile:", pylab.mean(ant_percentiles), "+/-", pylab.std(ant_percentiles, ddof=1)
        print "heuristic percentile:", pylab.mean(heuristic_percentiles), "+/-", pylab.std(heuristic_percentiles, ddof=1)

        random_differences = 50 - ant_percentiles
        heuristic_differences = ant_percentiles - heuristic_percentiles

        print random_differences
        print heuristic_differences

        print "ants vs random:", pylab.mean(random_differences), "+/-", pylab.std(random_differences, ddof=1)

        print "ants vs heuristic", pylab.mean(heuristic_differences), "+/-", pylab.std(heuristic_differences, ddof=1)

        print wilcoxon(heuristic_differences, random_differences)

def objectives_friedman(paper=False, max_ti=MAX_TRANSITION_INDEX, null_model=DEFAULT_NULL_MODEL):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df = df[df['max transition index'] == max_ti]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['algorithm'] == 'ants']

    samples = []
    for objective, objective_group in df.groupby(['objective']):
        objective_samples = []
        for colony, colony_group in objective_group.groupby('colony'):
            objective_samples.append(colony_group['percentile'].mean())
        samples.append(objective_samples)

    print friedmanchisquare(*samples)

def colonies_friedman(paper=False, max_ti=MAX_TRANSITION_INDEX, null_model=DEFAULT_NULL_MODEL):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df = df[df['max transition index'] == max_ti]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['algorithm'] == 'ants']

    samples = []
    for colony, colony_group in df.groupby(['colony']):
        colony_samples = []
        for objective, objective_group in colony_group.groupby('objective'):
            colony_samples.append(objective_group['percentile'].mean())
        samples.append(colony_samples)

    print friedmanchisquare(*samples)

def percentiles_conover(paper=False, max_ti=MAX_TRANSITION_INDEX, null_model=DEFAULT_NULL_MODEL):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df = df[df['max transition index'] == max_ti]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['algorithm'] == 'ants']

    conover_df = df[['objective', 'colony', 'percentile']]
    conover_df = conover_df.groupby(['objective', 'colony'], as_index=False).agg(pylab.mean)
    print posthoc_conover(conover_df, val_col='percentile', group_col='objective')

    for colony, colony_group in df.groupby('colony'):
        print '-' * len(colony)
        print colony
        print '-' * len(colony)

        print posthoc_conover(colony_group, val_col='percentile', group_col='objective')

def percentiles_wilcoxon(paper=False, max_ti=MAX_TRANSITION_INDEX, null_model=DEFAULT_NULL_MODEL):
    percentiles_fname = None
    if paper:
        percentiles_fname = '%s/%s.csv' % (PAPER_DIR, PAPER_FNAME)
    else:
        percentiles_fname = 'objective_percentiles.csv'

    df = pd.read_csv(percentiles_fname, skipinitialspace=True)
    df = df[df['objective'].isin(ANOVA_OBJECTIVES)]
    df = df[df['max transition index'] == max_ti]
    df = df[df['null model'] == NULL_MODEL_NAMES[null_model]]
    df = df[df['algorithm'] == 'ants']

    for colony, colony_group in df.groupby('colony'):
        print '-' * len(colony)
        print colony
        print '-' * len(colony)

        print posthoc_wilcoxon(colony_group, val_col='percentile', group_col='objective')


def sandbox():
    pass

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--write', action='store_true')
    parser.add_argument('-p', '--plot', action='store_true')
    parser.add_argument('-n', '--null', default=DEFAULT_NULL_MODEL)
    parser.add_argument('-a', '--algorithm', default=DEFAULT_ALGORITHM)
    parser.add_argument('-r', '--random_trials', type=int, default=0)
    parser.add_argument('-c', '--changes', action='store_true')
    parser.add_argument('-t', '--test', action='store_true')
    parser.add_argument('--max_ti', type=int, default=MAX_TRANSITION_INDEX)
    parser.add_argument('--a1', action='store_true')
    parser.add_argument('--paper', action='store_true')
    parser.add_argument('--loops', action='store_true')
    parser.add_argument('--correlations', action='store_true')
    parser.add_argument('--compare_heuristic', action='store_true')
    parser.add_argument('--friedman', action='store_true')
    parser.add_argument('--conover', action='store_true')
    parser.add_argument('--wilcoxon', action='store_true')
    parser.add_argument('--swarm', action='store_true')
    parser.add_argument('--sandbox', action='store_true')

    args = parser.parse_args()
    null_model = args.null
    algorithm = args.algorithm
    random_trials = args.random_trials
    paper = args.paper
    max_ti = args.max_ti
    loops = args.loops
    anova1 = args.a1
    correlations = args.correlations
    swarm = args.swarm

    if args.sandbox:
        sandbox()
        exit()

    if args.write:
        write_objectives(random_trials, max_ti=max_ti, null_model=null_model)
    if args.plot:
        plot_changes(differences=False, paper=paper, max_ti=max_ti,\
                     null_model=null_model, algorithm=algorithm)
        plot_percentiles_distribution(paper=paper, max_ti=max_ti,\
                                      null_model=null_model,\
                                      algorithm=algorithm)
    if args.test:
        test_percentiles(barely=barely, max3=max3, paper=paper)
    if anova1:
        one_way_anova(barely=barely, max3=max3, paper=paper)
    if correlations:
        objective_correlations()
    if args.compare_heuristic:
        compare_heuristic(paper=paper, max_ti=max_ti, null_model=null_model)
    if args.friedman:
        objectives_friedman(paper=paper, max_ti=max_ti, null_model=null_model)
        colonies_friedman(paper=paper, max_ti=max_ti, null_model=null_model)
    if args.conover:
        percentiles_conover(paper=paper, max_ti=max_ti, null_model=null_model)
    if args.wilcoxon:
        percentiles_wilcoxon(paper=paper, max_ti=max_ti, null_model=null_model)
    if args.swarm:
        plot_percentiles_swarm(paper=paper, max_ti=max_ti, null_model=null_model)


if __name__ == '__main__':
    main()
