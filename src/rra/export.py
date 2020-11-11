from typing import Dict, List, Set, Tuple

from pandas import DataFrame


class DirectedHyperGraph:
    def __init__(self,
            df:DataFrame,
            reac_col:str='reaction',
            reac_label_col:str='reaction',
            comp_col:str='compound',
            comp_label_col:str='compound',
            ispr_col:str='isproduct'
    ) -> Tuple[List[int], List[int], List[int], List[int], List[int]]:
    
        self.nodes = {}
        for c in df[comp_col].drop_duplicates().values:
            substance = df[df[comp_col]==c]
            self.nodes[c] = set([str(l) for l in substance[comp_label_col].values])

        self.pseudo_nodes = {}
        self.nondirected = []
        self.directed = []
        for r in df.index.drop_duplicates().values:
            reaction = df.loc[r,:]
            for s in reaction[~reaction[ispr_col]][comp_col].values:
                self.nondirected.append((s, r))
            for t in reaction[reaction[ispr_col]][comp_col].values:
                self.directed.append((r, t))
            self.pseudo_nodes[r] = set([str(l) for l in reaction.index.values])


def dot_reactions(df:DataFrame,
        std_df:DataFrame=None,
        reac_col:str='reaction',
        reac_label_col:str='reaction',
        comp_col:str='compound',
        comp_label_col:str='compound',
        ispr_col:str='isproduct',
        dot_file:str='.reactions.dot',
        png_file:str=None
) -> None:
    from pydotplus.graphviz import graph_from_dot_file

    dhg = DirectedHyperGraph(df, reac_col, reac_label_col, comp_col, comp_label_col, ispr_col)
    
    with open(dot_file, "w") as dot:
        dot.write("digraph {\n")
        for node, label in dhg.nodes.items():
            dot.write(f' {node} [label="{", ".join(label)}", color="blue"]\n')
        for node, label in dhg.pseudo_nodes.items():
            dot.write(f' "r:{node}" [label="{", ".join(label)}"]\n')
        dot.write(" subgraph Arrowd {\n")
        for (ds,dt) in dhg.directed:
            dot.write(f'  "r:{ds}" -> {dt}\n')
        dot.write(" }\n")
        dot.write(" subgraph Blunt {\n  edge [dir=none]\n")
        for (ns,nt) in dhg.nondirected:
            dot.write(f'  {ns} -- "r:{nt}"\n')
        dot.write(" }\n")
        
        if std_df is not None:
            dot.write(' subgraph Equival {\n  edge [dir=none, color="red"]\n')
            def add_standard_edge(row):
                cmpd = row['compound']
                if row.eqcl in dhg.nodes and cmpd in dhg.nodes:
                    dot.write(f'  {row.eqcl} -- {cmpd}\n')
            std_df.apply(add_standard_edge, axis=1)
            dot.write(" }\n")
        
        dot.write("}\n")
    
    graph = graph_from_dot_file(dot_file)
    png = graph.create_png()
    if png_file:
        with open(png_file, "wb") as write:
            write.write(png)
    return png
