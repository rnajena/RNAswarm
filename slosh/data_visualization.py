import nxviz as nv
import matplotlib.pyplot as plt

# We assume you have a graph G that is a NetworkX graph object.
# In this example, all nodes possess the "group" and "value" node attributes
# where "group" is categorical and "value" is continuous,
# and all edges have the "edge_value" node attribute as well.

ax = nv.circos(
    G,
    group_by="group",
    sort_by="value",
    node_color_by="group",
    edge_alpha_by="edge_value"
)

nv.annotate.circos_group(G, group_by="group")