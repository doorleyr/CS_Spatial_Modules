import urllib
import geopandas as gpd
import osmnx
import PreCompOsmNet
import numpy as np
import requests
import json
import random
import pandas as pd

import OpenCity
from brix import Indicator, Handler

def mode_choice_model(all_trips_df):
#     all_trips_df['mode']=random.choices(['walk', 'drive'], k=len(all_trips_df))
    all_trips_df['route_distance']=all_trips_df.apply(lambda row: m.route_lengths[row['from_possible_nodes_drive'][0]]
                                                                        [row['to_possible_nodes_drive'][0]], axis=1)
    all_trips_df['mode']=all_trips_df.apply(lambda row: 'drive' if row['route_distance']>1500 else 'walk', axis=1)
    return all_trips_df

class Mobility_indicator(Indicator):
    def setup(self, state, table_name, model_radius, external_hw_tags=["motorway","motorway_link","trunk","trunk_link"]):
        self.N_max=2000
        self.external_hw_tags=external_hw_tags
        self.state=state
        self.table_name=table_name
        if not hasattr(self.state, 'od'):
            self.state.get_lodes_data(include=['od'])
        self.get_geogrid()
        self.get_overlap_geoids()
        self.state.geom['sim_area']=self.state.geom.index.isin(self.overlapping_geoids)
        self.state.subset_geom_by_distance([self.geogrid['x_centroid'].mean(), self.geogrid['y_centroid'].mean()],
                                           model_radius, 'model_area')
        simpop_df=state.lodes_to_pop_table(
            model_subset_name='model_area', sim_subset_name='sim_area')
        self.base_simpop_df=simpop_df.copy()
        self.build_mobsys() 
        print('Init simulation')
        grid_zones=self.create_grid_zones()
        model_zones=self.state.return_geometry('model_area')
        model_zones['grid_area']=False
        combined_zones=model_zones.append(grid_zones)
        sim_geoids=list(self.state.return_geometry('sim_area').index)+list(grid_zones.index)
        self.sim=OpenCity.Simulation(simpop_df, self.mob_sys, combined_zones, sim_geoids=sim_geoids)
        self.sim.set_choice_models(mode_chooser=mode_choice_model)
        self.create_zone_dist_mat()

    def create_grid_zones(self):
        grid_zones=self.geogrid.copy()
        centroids=grid_zones['geometry'].centroid
        grid_zones['x_centroid']=[c.x for c in centroids]
        grid_zones['y_centroid']=[c.y for c in centroids]
        cols_to_init=['total_pop_rac']+ [col for col in self.state.geom.columns if (
            ('emp' in col) or ('res' in col))]
        for col in cols_to_init:
            grid_zones[col]=0
        for area_col in ['model_area', 'sim_area', 'grid_area']:
            grid_zones[area_col]=1
        return grid_zones
        
    def get_geogrid(self):
        # repetition with proximity indicator
        get_url='https://cityio.media.mit.edu/api/table/'+self.table_name
        with urllib.request.urlopen(get_url+'/GEOGRID') as url:
            self.geogrid=gpd.read_file(url.read().decode())
        self.geogrid.set_crs("EPSG:4326")
        centroids=self.geogrid['geometry'].centroid
        self.geogrid['x_centroid']=[c.x for c in centroids]
        self.geogrid['y_centroid']=[c.y for c in centroids]
                       
    def get_overlap_geoids(self):
        # repetition with proximity indicator
        """
        find the geoids of the baseline zones which overlap with hthe geogrid
        
        """
        self.state.geom['copy_GEOID']=self.state.geom.index
        grid_intersect_zones=gpd.overlay(self.geogrid, self.state.geom, 'intersection')
        self.overlapping_geoids=grid_intersect_zones['copy_GEOID'].unique()
        
    def build_mobsys(self):
        print('Building Mobility System')
        print('\t getting graph')
        G_drive_sim = self.get_graph_buffered_to_hw_type(self.state.geom.loc[self.state.geom['sim_area']], 
                                                       self.external_hw_tags, 'drive')
        print('getting external roads')
        G_drive_model=osmnx.graph_from_polygon(self.state.geom.loc[self.state.geom['model_area']].unary_union,
                                 network_type='drive', custom_filter='["highway"~"{}"]'.format('|'.join(self.external_hw_tags)))
        G_drive_combined=osmnx.graph.nx.compose(G_drive_model,G_drive_sim)
#         for edge in list(G_walk.edges):
#             G_walk.edges[edge]['speed_kph']=4.8
        G_drive_combined=osmnx.add_edge_speeds(G_drive_combined)
        G_drive_combined=PreCompOsmNet.simplify_network(G_drive_combined)
        G_drive_combined=osmnx.add_edge_travel_times(G_drive_combined)
        
        for edge in list(G_drive_combined.edges):
            G_drive_combined.edges[edge]['travel_time_walk']=G_drive_combined.edges[edge]['length']/(4800/3600)

        G_drive_combined=osmnx.utils_graph.get_undirected(G_drive_combined)
#         fw_pred_drive=PreCompOsmNet.pre_compute_paths(G_drive_combined)

        # Note: this approach will assume the same route for each mode but different travel times
        # For different routes, would need to compute the fw result for each mode
        fw_all, self.route_lengths=PreCompOsmNet.pre_compute_paths(G_drive_combined, 
                                                         weight_metric='length', 
                                                         save_route_costs=True)
        pre_comp_drive=PreCompOsmNet.PreCompOSMNet(G_drive_combined, fw_all)
        networks={'drive': pre_comp_drive, 'walk': pre_comp_drive}
        
        drive_dict={
            'target_network_id': 'drive',
            'travel_time_metric': 'travel_time'}
        walk_dict={
            'target_network_id': 'walk',
            'travel_time_metric': 'travel_time_walk'}
        
        modes={'drive': OpenCity.Mode(drive_dict), 'walk': OpenCity.Mode(walk_dict)}

        self.mob_sys=OpenCity.MobilitySystem(modes=modes,
                              networks=networks)
        
    def create_zone_dist_mat(self):
        print('Computing zone-zone dist mat')
        zone_index=self.sim.zones.index
        self.dist_mat={}
        chosen_nodes_ls=[pn[0] for pn in self.sim.zones['possible_nodes_drive']]
        for i in range(len(self.sim.zones)):
            self.dist_mat[zone_index[i]]={}
            for j in range(len(self.sim.zones)):
                from_node=chosen_nodes_ls[i]
                to_node=chosen_nodes_ls[j]
                if from_node==to_node:
                    self.dist_mat[zone_index[i]][zone_index[j]]=50
                else:
                    self.dist_mat[zone_index[i]][zone_index[j]]=self.route_lengths[from_node][to_node] 
        
    def get_graph_buffered_to_hw_type(self,geom, external_hw_tags, network_type):
        for buffer in [i*200 for i in range(10)]:
            print('Buffer : {}'.format(buffer))
            geom_projected=osmnx.projection.project_gdf(geom)
            geom_projected_buffered=geom_projected.unary_union.buffer(buffer)

            geom_projected_buffered_gdf=gpd.GeoDataFrame(geometry=[geom_projected_buffered], crs=geom_projected.crs)
            geom_wgs_buffered_gdf=geom_projected_buffered_gdf.to_crs(geom.crs) 

            G_temp=osmnx.graph.graph_from_polygon(geom_wgs_buffered_gdf.iloc[0]['geometry'], network_type=network_type)
            all_hw_tags=[e[2]['highway'] for e in G_temp.edges(data=True)]

            if any([tag in all_hw_tags for tag in external_hw_tags]):
                return G_temp
            print('Doesnt contain external link types')
        print('Warning: internal network will not connect to external network')
        return G_temp 
        
        
    def routes_to_deckgl_trip(self, route_table):
        
        mode_id_map={'drive': "0", "cycle": "1", "walk": "2", "pt": "3"}
        profile_id_map={'u1250':"0", '1251-3333':"1",  '3333+':"2"}
        
        mode_split={}
        for mode in mode_id_map:
            mode_split[mode]=100*sum(route_table['mode']==mode)/len(route_table)
        attr_map={"mode": {
                    "0": {
                    "color": "#e41a1d",
                    "name": "Drive : {}%".format(mode_split['drive'])
                    },
                    "1": {
                    "color": "#377eb8",
                    "name": "Cycle : {}%".format(mode_split['cycle'])
                    },
                    "2": {
                    "color": "#4daf4a",
                    "name": "Walk : {}%".format(mode_split['walk'])
                    },
                    "3": {
                    "color": "#ffff33",
                    "name": "Public Transport : {}%".format(mode_split['pt'])
                    }},
                 "profile": {
                    "0": {
                    "color": "#7fc97f",
                    "name": "Low Income"
                    },
                    "1": {
                    "color": "#beaed4",
                    "name": "Med Income"
                    },
                    "2": {
                    "color": "#fdc086",
                    "name": "High Income"
                    }}}
                    #         earnings_to_ind={'u1250': 0, '1250to3333': 1, '3333plus': 2}
        trips=[]
        for ind, row in route_table.iterrows():
            coords=row['attributes']['coordinates']
            if len(coords)>1:
                travel_time_metric=self.mob_sys.modes[row['mode']].travel_time_metric
                cum_time=np.cumsum(row['attributes'][travel_time_metric])
                start_time=int(row['start_time'])
                timestamps=[int(start_time)] + [int(start_time)+ int(ct) for ct in cum_time]
                this_trip={'path': coords, 'timestamps': timestamps}
#                 for attr in ['earnings']:
#                     this_trip[attr]=row[attr]
#                 this_trip['earnings']=earnings_to_ind[row['earnings']]
                this_trip['mode']=mode_id_map[row['mode']],
                this_trip['profile']=profile_id_map[row['earnings']]
                trips.append(this_trip)
        return {'attr': attr_map, 'trips': trips}
    
    def geogrid_updates(self, geogrid_data):
        new_simpop=[]
        side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
        type_def=geogrid_data.get_type_info()
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            type_info=type_def[name]
            if type_info['NAICS'] is not None:
                height=cell['height']
                cell_area=side_length*side_length
                if isinstance(height, list):
                    height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:
                    sqm_pperson=50
                total_capacity=height*cell_area/sqm_pperson
                workers={code: int(type_info['NAICS'][code]*total_capacity) for code in  type_info['NAICS']}   
                for code in workers:
                    home_locations=self.sample_home_locations(i_c, '3333+', n=workers[code])
                    for i_w in range(workers[code]):
                        new_simpop.append({'work_geoid': i_c,'home_geoid': home_locations[i_w],
                                           'naics': code, 'earnings': '3333+',
                                          'age': '30-54'})
        return new_simpop
            
    def sample_home_locations(self, work_geoid, earnings, n):
        attraction=self.sim.zones['res_income_{}_rac'.format(earnings)]
        impedance=[self.dist_mat[hid][work_geoid] for hid in self.sim.zones.index]
        weights=np.divide(attraction,np.array(impedance))
        return np.random.choice(
            self.sim.zones.index, replace=True, p=weights/sum(weights), size=n)
        
    def simulate(self, simpop_df):
        print('Schedules and locations')
        simpop_df=self.sim.create_simple_HWH_schedules(simpop_df)
        print('Trip table')
        all_trips_df=self.sim.create_trip_table(simpop_df)
        print('Route table')
        route_table=self.sim.get_routes_table(all_trips_df)
        print('DeckGL')
        deckgl_trips=self.routes_to_deckgl_trip(route_table)
        return deckgl_trips
    
    def return_indicator(self, geogrid_data):
        print('Starting MM Update')
        new_simpop=self.geogrid_updates(geogrid_data)
        new_simpop_df=pd.DataFrame(new_simpop)
        combined_simpop=self.base_simpop_df.append(new_simpop_df)
        sample_simpop_df=combined_simpop.sample(min(self.N_max, len(combined_simpop)))
        deckgl_trips=self.simulate(sample_simpop_df)
        self.post_trips(deckgl_trips)
        print('Finished MM Update')
        return None
    
    def post_trips(self, deckgl_trips):
        post_url='https://cityio.media.mit.edu/api/table/update/'+self.table_name
        r = requests.post(post_url+'/ABM2', data = json.dumps(deckgl_trips))
        print('Post ABM: {}'.format(r))