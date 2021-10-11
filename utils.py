import networkx as nx

def node_order(G):
    terminals = sorted(G.graph['terminals'])
    terminal_distances = []
    for u in G.nodes():
        distances = []
        for terminal in terminals:
            distance = nx.shortest_path_length(G, u, terminal)
            distances.append(distance)
        distances.append(u)
        terminal_distances.append(tuple(distances))

    terminal_distances = sorted(terminal_distances)
    order = {}
    for i, distances in enumerate(terminal_distances):
        u = distances[-1]
        order[u] = i

    return order
