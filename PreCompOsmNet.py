import osmnx
import pandas as pd

def simplify_network(G, tolerance=10):
    Gp=osmnx.projection.project_graph(G)
    # simplify
    G_simp=osmnx.simplification.consolidate_intersections(
            Gp, tolerance=tolerance, rebuild_graph=True, dead_ends=False, 
            reconnect_edges=True)
    print('\t Simplified from {} to {} edges and {} to {} nodes'.format(
        len(G.edges), len(G_simp.edges), len(G.nodes), len(G_simp.nodes)))
    # project back to WGS
    G_simp_wgs=osmnx.projection.project_graph(G_simp, to_crs='epsg:4326')
    return G_simp_wgs

def pre_compute_paths(G, weight_metric='travel_time', save_route_costs=False):
    print('\t Pre-computing paths')
    fw_pred, fw_dist=osmnx.graph.nx.floyd_warshall_predecessor_and_distance(
        G, weight=weight_metric)
    if save_route_costs:
    	return fw_pred, fw_dist
    else:
    	return fw_pred
 

class PreCompOSMNet():
    def __init__(self, G, pred):
        """
        Takes as input an osmnx network object, with a 'speed kph' attribute already added
        On intialisation, generates predecessors
        """
        for i_e, edge in enumerate(list(G.edges)):
            G.edges[edge]['edge_id']=i_e
        self.G=G
        self.predecessors=pred
        
    # def pre_compute_paths(self, save_travel_time=False):
    #     print('\t Pre-computing paths')
    #     fw_pred, fw_dist=osmnx.graph.nx.floyd_warshall_predecessor_and_distance(
    #         self.G, weight='travel_time')
    #     self.predecessors=fw_pred
    #     if save_dist:
    #     	self.travel_times=fw_dist
        
    def shortest_paths(self, nodes_a, nodes_b, imp_name=None):
        paths=[]
        for i in range(len(nodes_a)):
            paths.append(osmnx.graph.nx.algorithms.shortest_paths.dense.reconstruct_path(
                nodes_a[i], nodes_b[i], self.predecessors))
        return paths
        
    def get_node_ids(self, x_col, y_col):
        return osmnx.distance.get_nearest_nodes(self.G, x_col, y_col, method='kdtree')

    def get_nodes_df(self):
        return pd.DataFrame(data=[{'x': self.G.nodes[n]["x"], 
                    'y': self.G.nodes[n]["y"], 
                    'id': n} for n in self.G.nodes]).set_index('id')

    def get_path_link_attributes(self, path, attribute_names=['travel_time']):
        # https://github.com/gboeing/osmnx/blob/master/osmnx/plot.py for actual edge geometries
        output={attr: [] for attr in attribute_names}
        if len(path)>1:
            coordinates=[]
            edges=[]
            for u, v in zip(path[:-1], path[1:]):
#                 edges.append((u,v))
                # if there are parallel edges, select the shortest in length
                data = min(self.G.get_edge_data(u, v).values(), key=lambda d: d["length"])
                edges.append(data['edge_id'])
                coordinates.append([
                    self.G.nodes[u]["x"], 
                    self.G.nodes[u]["y"]])
                for attr in attribute_names:
                    output[attr].append(data[attr])
            coordinates.append([
                self.G.nodes[path[-1]]["x"], 
                self.G.nodes[path[-1]]["y"]])
            output['coordinates']=coordinates
            output['edges']=edges
            return output
        else:
            output['coordinates']=[]
            output['edges']=[]
            return output
