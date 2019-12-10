import networkx as nx
import numpy as np

from .genomics import MiBIGBGC

def _edges_dict_from_networkx(graph):
    return {'start': [x[0] for x in graph.edges()], 'end': [x[1] for x in graph.edges()]}

def _combined_edges_dict(g1, g2):
    e1 = _edges_dict_from_networkx(g1)
    e2 = _edges_dict_from_networkx(g2)
    e1['start'].extend(e2['start'])
    e1['end'].extend(e2['end'])
    return e1

def _layout_molnet(G, width, iterations, offset=(0, 0)):
    ncomp = nx.number_connected_components(G)
    nrows = np.ceil(np.sqrt(ncomp))
    cc = [G.subgraph(c) for c in nx.connected_components(G)]
    allpos = {}
    scale = (width / nrows) / 4
    for i, ccc in enumerate(cc):
        # compute center
        col_no = i % nrows
        row_no = i // nrows
        col_pos = width * (col_no / nrows) - (width / 2)
        row_pos = width * (row_no / nrows) - (width / 2)
        center = (offset[0] + row_pos, offset[1] + col_pos)
        pos = nx.spring_layout(ccc, iterations=iterations, threshold=1e-6, center=center, scale=scale)
        allpos.update(pos)
    return allpos

def create_metabolomics_graph(spectra, width=10, iterations=100, singletons=False, split_singletons=False):
    if not singletons and split_singletons:
        raise Exception('singletons must be True if split_singletons is True!')

    # if split_singletons is False, combine everything into a single graph
    if not split_singletons:
        print('Constructing single combined graph, singletons={}'.format(singletons))
        g = nx.Graph()
        for spec in spectra:
            if spec.family.family_id == -1 and not singletons:
                continue # ignore singletons

            g.add_node(spec.id)
            for e_iid, e_sid, e_cosine in spec.edges:
                if spectra[e_iid].family_id == -1 and not singletons:
                    continue 
                g.add_node(e_iid)
                g.add_edge(spec.id, e_iid)
        print('Graph: nodes={}, edges={}'.format(g.number_of_nodes(), g.number_of_edges()))
        return _layout_molnet(g, width, iterations), _edges_dict_from_networkx(g)
    else:
        print('Constructing split-graph, singletons={}'.format(singletons))
        gn = nx.Graph()
        gs = nx.Graph()
        
        for spec in spectra:
            g = gs if spec.family.family_id == -1 else gn
            g.add_node(spec.id)
            for e_iid, e_sid, e_cosine in spec.edges:
                g.add_node(e_iid)
                g.add_edge(spec.id, e_iid)
        
        print('Graph normal: nodes={}, edges={}'.format(gn.number_of_nodes(), gn.number_of_edges()))
        print('Graph singletons: nodes={}, edges={}'.format(gs.number_of_nodes(), gs.number_of_edges()))
        # dict of node ID: (x,y) entries
        gn_layout = _layout_molnet(gn, width, iterations)
        gs_layout = _layout_molnet(gs, width, iterations, offset=(0, -width * 1.1))
        gn_layout.update(gs_layout)

        return gn_layout, _combined_edges_dict(gn, gs)

def create_genomics_graph(bgcs, width=10, iterations=100, mibig=False, split_mibig=False):
    if not mibig and split_mibig:
        raise Exception('mibig must be True if split_mibig is True')
        
    if not split_mibig:
        print('Constructing single combined graph, mibig={}'.format(mibig))
        g = nx.Graph()
        for bgc in bgcs:
            if isinstance(bgc, MiBIGBGC) and not mibig:
                continue
            g.add_node(bgc.id)
            for obgc_id in bgc.edges:
                if isinstance(bgcs[obgc_id], MiBIGBGC) and not mibig:
                    continue
                g.add_node(obgc_id)
                g.add_edge(bgc.id, obgc_id)
        print('Graph: nodes={}, edges={}'.format(g.number_of_nodes(), g.number_of_edges()))
        return _layout_molnet(g, width, iterations), _edges_dict_from_networkx(g)
    else:
        print('Constructing split-graph, mibig={}'.format(mibig))
        gn = nx.Graph()
        gs = nx.Graph()
        
        for bgc in bgcs:
            g = gs if isinstance(bgc, MiBIGBGC) else gn
            g.add_node(bgc.id)
            for obgc_id in bgc.edges:
                g.add_node(obgc_id)
                g.add_edge(bgc.id, obgc_id)
            
        print('Graph normal: nodes={}, edges={}'.format(gn.number_of_nodes(), gn.number_of_edges()))
        print('Graph singletons: nodes={}, edges={}'.format(gs.number_of_nodes(), gs.number_of_edges()))
        # dict of node ID: (x,y) entries
        gn_layout = _layout_molnet(gn, width, iterations)
        gs_layout = _layout_molnet(gs, width, iterations, offset=(0, -width * 1.1))
        gn_layout.update(gs_layout)
        return gn_layout, _combined_edges_dict(gn, gs)
