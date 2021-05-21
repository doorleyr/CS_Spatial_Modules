from brix import Indicator, Handler
import OpenCity
import PreCompOsmNet
import statsmodels
from statsmodels.distributions.empirical_distribution import ECDF
import geopandas as gpd
import urllib
import json
import requests
import math
import osmnx
from shapely.geometry import Point, shape
import datetime
from scipy import spatial
import numpy as np
import sys
import pandas as pd
import copy

"""
"""

def aggregate_types_over_grid(geogrid_data, side_length, type_def):
    cell_area=side_length*side_length
    aggregated={}  
    for cell in geogrid_data:
        name=cell['name']
        type_info=type_def[name]
        height=cell['height']
        if isinstance(height, list):
            height=height[-1]
        if 'sqm_pperson' in type_info:
            sqm_pperson=type_info['sqm_pperson']
        else:
            sqm_pperson=50
        total_capacity=height*cell_area/sqm_pperson
        if name in aggregated:
            aggregated[name]+=total_capacity
        else:
            aggregated[name]=total_capacity
    return aggregated

def aggregate_attributes_over_grid(agg_types, attribute, side_length, type_def, digits=None):
    # TODO: eliminate repetition with previous function
    cell_area=side_length*side_length
    aggregated={}
    for cell_type in agg_types:
        type_info=type_def[cell_type]
        if type_info[attribute] is not None:
            total_capacity=agg_types[cell_type]
            for code in type_info[attribute]:
                if digits==None:
                    code_digits=code
                else:
                    code_digits=code[0:digits]
                attr_capacity=total_capacity*type_info[attribute][code]
                if code_digits in aggregated:
                    aggregated[code_digits]+= attr_capacity
                else:
                    aggregated[code_digits]= attr_capacity
    return aggregated

def get_central_nodes(geodf, G):
    """ takes a geodf and returns a list of nodes closest to their centroids
    returns both the nearest nodes and the distance
    """
    geodf_projected=osmnx.projection.project_gdf(geodf)
    projected_centroids=geodf_projected['geometry'].centroid.values
    projected_centroids_lst=[[c.x, c.y] for c in projected_centroids]
    G_proj=osmnx.projection.project_graph(G, geodf_projected.crs)
    G_proj_coordinates={node: [G_proj.nodes[node]['x'], G_proj.nodes[node]['y']] for node in G_proj.nodes}
    node_order=[node for node in G_proj_coordinates]
    nodes_kdtree=spatial.KDTree([G_proj_coordinates[node] for node in node_order])
    dist, ind_nearest=nodes_kdtree.query(projected_centroids_lst)
    nearest_nodes=[node_order[ind] for ind in ind_nearest]
    return nearest_nodes, dist


class Proximity_Indicator(Indicator):
    def setup(self, zones, geogrid, buffer=1200):
        print('Setting up Proximity Indicator')
        self.zone_to_node_tolerance=500
        self.grid_to_node_tolerance=100
        self.naics_codes=[col for col in zones.columns if 'naics' in col]
        self.buffer=buffer
        self.zones=zones
        self.indicator_type = 'hybrid'
        self.geogrid=geogrid
        self.overlapping_geoids=list(zones.loc[zones['sim_area']].index)
        self.all_site_ids=self.overlapping_geoids+list(range(len(geogrid)))  
        self.get_graph_reference_area()
        
        print('\t Getting central nodes')
        zones_nearest_nodes, zones_nearest_dist= get_central_nodes(self.zones, self.ref_G)
        self.zones['central_node']=zones_nearest_nodes
        self.zones['nearest_dist']=zones_nearest_dist
        
        grid_nearest_nodes, grid_nearest_dist= get_central_nodes(self.geogrid, self.ref_G)
        self.geogrid['central_node']=grid_nearest_nodes
        self.geogrid['nearest_dist']=grid_nearest_dist

        self.get_reachable_geoms_from_all()
        self.calculate_baseline_scores()
            
            
    def make_ego_graph_around_geometry(self, zone, tolerance):
        """
        For geometries within the buffered geogrid only.
        Returns the graph within a walkable distance of the centre of the zone
        """
        if zone['nearest_dist']<tolerance:
            sub_graph = osmnx.graph.nx.ego_graph(self.ref_G, zone['central_node'], radius=1200, distance='length')
        else:
            sub_graph = osmnx.graph.nx.Graph()
        return sub_graph
    
    def get_graph_reference_area(self):
        reference_zones=self.zones.loc[self.zones['reference_area']]
        print('\t Downloading graph for reference area')
        reference_zone_graph=self.get_network_around_geom_buffered(reference_zones)
        self.ref_G=reference_zone_graph
        
    def calculate_baseline_scores(self):
        # TODO: should use the get_reachable_geoms function?
        print('\t Calculating baseline scores')
#         self.base_scores={'walkable_{}'.format(x): [] for x in [
#             'employment', 'housing', 'healthcare', 'hospitality', 'shopping']}
        base_scores={}
        self.base_attributes={}
        self.score_ecdfs={}
        stats_to_aggregate=[col for col in self.zones.columns if (('res_' in col) or ('emp_' in col))]
        # get the baseline reachable attributes and scores for every zone
        for ind, row in self.zones.loc[self.zones['reference_area']].iterrows():
            reachable_zones=self.zone_to_reachable[ind]['zones']
            self.base_attributes[ind]=self.zones.loc[reachable_zones][stats_to_aggregate].sum().to_dict()
            self.base_attributes[ind]['source_res']=row['res_total']
            self.base_attributes[ind]['source_emp']=row['emp_total']
            # get scores for individual zones- weighting cancels out
            base_scores[ind]=self.attributes_to_scores([self.base_attributes[ind]])
            
        # Create the ECDFs for each score using only the zones (not the grid cells
        self.base_zones_scores=pd.DataFrame.from_dict(base_scores, orient='index')
        for score in self.base_zones_scores.columns:
            base_scores_no_nan=[x for x in self.base_zones_scores[score] if x==x]
            self.score_ecdfs[score]=ECDF(base_scores_no_nan)
            
        # get weighted scores across the simulation area zones 
        # (ignore the grid which is empty in reference and therefore would be weighted zero)
        ref_scores=self.attributes_to_scores([self.base_attributes[ind] for ind in self.overlapping_geoids])
        self.ref_ind=self.normalise_ind(ref_scores)
            
        # get the base reachable attributes for every grid cell location
        for i_c in range(len(self.geogrid)):
            reachable_zones=self.grid_to_reachable[i_c]['zones']
            self.base_attributes[i_c]=self.zones.loc[reachable_zones][stats_to_aggregate].sum().to_dict()
            
            
    def get_reachable_geoms(self, zone, tolerance):
        """
        find all grid cells and all zones reachable from a geometry
        """
        sub_graph=self.make_ego_graph_around_geometry(zone, tolerance)
        sub_graph_nodes=sub_graph.nodes(data=False)
        reachable_zones= list(self.zones.loc[
                ((self.zones['central_node'].isin(list(sub_graph_nodes)))&
                 (self.zones['nearest_dist']<self.zone_to_node_tolerance))
                ].index.values)
        reachable_grid_cells=list(self.geogrid.loc[
                ((self.geogrid['interactive']=='Web')&
                    (self.geogrid['central_node'].isin(list(sub_graph_nodes)))&
                    (self.geogrid['nearest_dist']<self.grid_to_node_tolerance))
                ].index.values)
        return {'zones': reachable_zones, 'cells': reachable_grid_cells}
    
    def get_reachable_geoms_from_all(self):
        """
        For every grid cell and every zone which intersects the grid:
        find the reachable zones and reachable grid cells
        """
        print('\t Finding all reachable geometries from each geometry')
        self.grid_to_reachable, self.zone_to_reachable={}, {}
        for ind, row in self.geogrid.iterrows():
            self.grid_to_reachable[ind]=self.get_reachable_geoms(row, self.grid_to_node_tolerance)
        for ind, row in self.zones.loc[self.zones['reference_area']].iterrows():
            self.zone_to_reachable[ind]=self.get_reachable_geoms(row, self.zone_to_node_tolerance)
        # create a reverse lookup to map from eacg grid cell to the cells from which it is reachable
        self.grid_to_reverse_reachable={}
        for i, row_i in self.geogrid.iterrows():
            self.grid_to_reverse_reachable[i]={'zones': [], 'cells': []}
            for j, row_j in self.geogrid.iterrows():
                if i in self.grid_to_reachable[j]['cells']:
                    self.grid_to_reverse_reachable[i]['cells'].append(j)
            for ind_z, row_z in self.zones.loc[self.zones['reference_area']].iterrows():
                if i in self.zone_to_reachable[ind_z]['cells']:
                    self.grid_to_reverse_reachable[i]['zones'].append(ind_z)                    
            
    def attributes_to_scores(self, attributes):        
        total_emp=sum([s['source_emp'] for s in attributes])
        total_res=sum([s['source_res'] for s in attributes])

        scores={}

        scores['walkable_housing']=sum([s['source_emp']*s['res_total'] for s in attributes])/total_emp
        scores['walkable_employment']=sum([s['source_res']*s['emp_total'] for s in attributes])/total_res
        scores['walkable_healthcare']=sum([s['source_res']*s['emp_naics_62'] for s in attributes])/total_res
        scores['walkable_hospitality']=sum([s['source_res']*s['emp_naics_72'] for s in attributes])/total_res
        scores['walkable_shopping']=sum([s['source_res']*s['emp_naics_44-45'] for s in attributes])/total_res
        return scores
    
    def return_indicator(self, geogrid_data):
        start_ind_calc=datetime.datetime.now()
        # make copy of base_scores
        new_attributes=self.get_new_reachable_attributes(geogrid_data)
        new_attributes_site= [new_attributes[ind] for ind in self.all_site_ids]
        new_scores=self.attributes_to_scores(new_attributes_site)
        new_ind=self.normalise_ind(new_scores)
        outputs=[]
        for ind_name in new_scores:
            outputs.append({'name': ind_name ,
                            'value': new_ind[ind_name],
                            'raw_value': new_scores[ind_name],
                            'ref_value': self.ref_ind[ind_name],
                            'viz_type': self.viz_type})
        end_ind_calc=datetime.datetime.now()
        
        new_attributes_site= [new_attributes[i_c] for i_c in range(len(geogrid_data))]
        heatmap=self.compute_heatmaps(new_attributes_site)
        
        end_hm_calc=datetime.datetime.now()
        
        print('Prox Ind: {}'.format(end_ind_calc-start_ind_calc))
        print('Prox HM: {}'.format(end_hm_calc-end_ind_calc))
        
        return {'heatmap':heatmap,'numeric':outputs}
    
    
    def get_new_reachable_attributes(self, geogrid_data):
        new_attributes=copy.deepcopy(self.base_attributes)
        side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
        cell_area=side_length*side_length
        type_def=geogrid_data.get_type_info()
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            height=cell['height']
            if isinstance(height, list):
                height=height[-1]
            type_info = type_def[name] 
            if 'sqm_pperson' in type_info:
                sqm_pperson=type_info['sqm_pperson']
            else:
                sqm_pperson=50
            total_capacity=height*cell_area/sqm_pperson
            agg_naics=aggregate_attributes_over_grid(
                {name: total_capacity}, 'NAICS', side_length, type_def, digits=2)
            agg_lbcs=aggregate_attributes_over_grid(
                {name: total_capacity}, 'LBCS', side_length, type_def, digits=1)

            added_attributes={}
            cell_employment=sum(agg_naics.values())
            new_attributes[i_c]['source_emp']=cell_employment
            added_attributes['emp_total']=cell_employment

            if '1' in agg_lbcs:
                cell_population=agg_lbcs['1']
            else:
                cell_population=0

            added_attributes['res_total']=cell_population
            new_attributes[i_c]['source_res']=cell_population

            for combined_code in self.naics_codes:
                naics_codes=combined_code.split('naics_')[1].split('-')
                for code in naics_codes:
                    if code in agg_naics:
                        if code in added_attributes:
                            added_attributes[combined_code]+=agg_naics[code]
                        else:
                            added_attributes[combined_code]=agg_naics[code]
            # the newly attributes are added to every zone and cell which can reach this cell   
            reverse_reachable=self.grid_to_reverse_reachable[i_c]
            for ind_z in reverse_reachable['zones']:
                for attr in added_attributes:
                    new_attributes[ind_z][attr]+=added_attributes[attr]
            for j_c in reverse_reachable['cells']:
                for attr in added_attributes:
                    new_attributes[j_c][attr]+=added_attributes[attr] 
        return new_attributes

    def normalise_ind(self, raw_ind):
        norm_ind={}
        for ind_name in raw_ind:
            norm_ind[ind_name]=self.score_ecdfs[ind_name](raw_ind[ind_name])       
        return norm_ind


    def compute_heatmaps(self, grid_reachable_area_stats):
        max_scores={score: self.base_zones_scores[score].max() for score in self.base_zones_scores}
        features=[]
        heatmap={'type': 'FeatureCollection',
                 'properties': ['housing', 'employment', 'healthcare', 'hospitality', 'shopping']}
        x_centroid_list, y_centroid_list=self.geogrid['x_centroid'], self.geogrid['y_centroid']
        for i_c, cell_stats in enumerate(grid_reachable_area_stats):
            features.append({
              "type": "Feature",
              "properties": [min((cell_stats['res_total']/max_scores['walkable_housing']), 1), 
                             min((cell_stats['emp_total']/max_scores['walkable_employment']), 1),
                             min((cell_stats['emp_naics_62']/max_scores['walkable_healthcare']), 1), 
                             min((cell_stats['emp_naics_72']/max_scores['walkable_hospitality']), 1), 
                             min((cell_stats['emp_naics_44-45']/max_scores['walkable_shopping']), 1)],
              "geometry": {
                "type": "Point",
                "coordinates": [
                  x_centroid_list[i_c],
                  y_centroid_list[i_c]
                ]
              }
            })
        heatmap['features']=features
        return heatmap
      
    def get_network_around_geom_buffered(self, geom):
        """
        Creates a buffer around the geometry and downloads the graph for this area
        """
        geom_projected=osmnx.projection.project_gdf(geom)
        geom_projected_buffered=geom_projected.unary_union.buffer(self.buffer)

        geom_projected_buffered_gdf=gpd.GeoDataFrame(geometry=[geom_projected_buffered], crs=geom_projected.crs)
        geom_wgs_buffered_gdf=geom_projected_buffered_gdf.to_crs(geom.crs) 
        
        return osmnx.graph.graph_from_polygon(geom_wgs_buffered_gdf.iloc[0]['geometry'], network_type='walk')

class Density_Indicator(Indicator):
    def setup(self, zones):
        self.zones=zones
        self.overlapping_geoids=list(zones.loc[zones['sim_area']].index)
        self.grid_cell_area=None
        self.compute_base_score_distribution()
                
    def compute_base_score_distribution(self):
        """
        computes ECDFs of the indicators across the baseline zones
        the ECDFs are later used to nromalise indicators by finding the percentile rank of the updated site wrt
        the baseline zones
        """
        self.score_ecdfs={}
        self.base_scores={}
        self.base_scores['res_density']=self.zones.apply(lambda row: self.res_density(row), axis=1)
        self.base_scores['emp_density']=self.zones.apply(lambda row: self.emp_density(row), axis=1)
        self.base_scores['live_work_score']=self.zones.apply(
            lambda row: self.get_live_work_score(row), axis=1)
        
        # Diversity
        industry_columns=[col for col in self.zones.columns if 'emp_naics' in col]
        res_income_columns=[col for col in self.zones.columns if 'res_income' in col]
        
        self.base_scores['industry_diversity']=self.zones.apply(
            lambda row: self.get_diversity(row, species_columns=industry_columns), axis=1)
        self.base_scores['income_diversity']=self.zones.apply(
            lambda row: self.get_diversity(row, species_columns=res_income_columns), axis=1)
        
        for score in self.base_scores:
#             plt.figure()
#             plt.hist(self.base_scores[score])
#             plt.title(score)
            base_scores_no_nan=[x for x in self.base_scores[score] if x==x]
            self.score_ecdfs[score]=ECDF(base_scores_no_nan)
        
    def combine_site_attributes(self, geogrid_data=None):
        """
        takes the attributes of the geogrid_data (new programming) and 
        the zones which overlap with the geogrid and (pre-existing programming)
        aggregates them together to get the updated site stats
        """
        stats_to_aggregate=['area']+[col for col in self.zones.columns if (('res_' in col) or ('emp_' in col))]
        temp_site_stats=dict(self.zones.loc[self.overlapping_geoids, 
                                                 stats_to_aggregate].sum())
        if geogrid_data is not None:
            side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
            type_def=geogrid_data.get_type_info()
            agg_types=aggregate_types_over_grid(geogrid_data, side_length=side_length, type_def=type_def)
            agg_naics=aggregate_attributes_over_grid(agg_types, 'NAICS', side_length=side_length, type_def=type_def, digits=2)
            agg_lbcs=aggregate_attributes_over_grid(agg_types, 'LBCS', side_length=side_length, type_def=type_def, digits=1)
            
            # update total residential and total employment
            add_emp=sum(agg_naics.values())
            if '1' in agg_lbcs:
                add_res=agg_lbcs['1']
            else:
                add_res=0    
            temp_site_stats['res_total']+=add_res
            temp_site_stats['emp_total']+=add_emp
            
            # update employment for each NAICS code
            for col in temp_site_stats:
                if 'naics' in col:
                    # if this is a double naics code column (eg. 44-45), add the new employment for both 44 and 45
                    col_naics_codes=col.split('naics_')[1].split('-')
                    for code in col_naics_codes:
                        if code in agg_naics:
                            temp_site_stats[col]+=agg_naics[code]  
                            
            # update residential types
            if 'Residential Low Income' in agg_types:
                temp_site_stats['res_income_low']+=agg_types['Residential Low Income']
            if 'Residential Med Income' in agg_types:
                temp_site_stats['res_income_mid']+=agg_types['Residential Med Income']
            if 'Residential High Income' in agg_types:
                temp_site_stats['res_income_high']+=agg_types['Residential High Income']          
        return temp_site_stats

    
    def calculate_indicators(self, site_stats):
        raw_ind={}
        raw_ind['res_density']=self.res_density(site_stats)
        raw_ind['emp_density']=self.emp_density(site_stats)
        raw_ind['live_work_score']=self.get_live_work_score(site_stats)
        
        industry_columns=[col for col in self.zones.columns if 'emp_naics' in col]
        res_income_columns=[col for col in self.zones.columns if 'res_income' in col]
        
        raw_ind['industry_diversity']=self.get_diversity(site_stats, species_columns=industry_columns)
        raw_ind['income_diversity']=self.get_diversity(site_stats, species_columns=res_income_columns)
               
        norm_ind={}
        for ind_name in raw_ind:
            norm_ind[ind_name]=self.score_ecdfs[ind_name](raw_ind[ind_name])       
        return {'raw': raw_ind, 'norm': norm_ind}
                  
    def return_indicator(self, geogrid_data):
        start_ind_calc=datetime.datetime.now()
        new_site_stats=self.combine_site_attributes(geogrid_data=geogrid_data)
        new_ind=self.calculate_indicators(new_site_stats)
        
        base_site_stats=self.combine_site_attributes(geogrid_data=None)
        base_ind=self.calculate_indicators(base_site_stats)
        
        outputs=[]
        for ind_name in new_ind['raw']:
            outputs.append({'name': ind_name.replace('_', ' ').title(),
                           'raw_value': new_ind['raw'][ind_name],
                           'value': new_ind['norm'][ind_name],
                           'ref_value': base_ind['norm'][ind_name]})
        end_ind_calc=datetime.datetime.now()
        print('Dens Ind: {}'.format(end_ind_calc-start_ind_calc))
#         print(outputs)
        return outputs
    
    @staticmethod
    def res_density(obj):
        """
        input can be a row if the baseline geodataframe
        or a dict representing a dynamic site
        """
        if obj['area']>0:
            return obj['res_total']/obj['area']
        return 0
    
    @staticmethod
    def emp_density(obj):
        """
        input can be a row if the baseline geodataframe
        or a dict representing a dynamic site
        """
        if obj['area']>0:
            return obj['emp_total']/obj['area'] 
        return 0
    
    @staticmethod
    def get_live_work_score(obj):
        if obj['emp_total']*obj['res_total']==0:
            return 0
        if obj['emp_total'] > obj['res_total']:
            return obj['emp_total']/obj['res_total']
        else:
            return obj['res_total']/obj['emp_total']
     
    @staticmethod
    def get_diversity(obj, species_columns):
        species_counts=[obj[col] for col in species_columns]
        diversity=0
        pop_size=sum(species_counts)
        if ((len(species_counts)>1) and (pop_size>0)):        
            for count in species_counts:
                pj=count/pop_size
                if not pj==0:
                    diversity+= -pj*math.log(pj)
            equitability=diversity/math.log(len(species_counts))
            return equitability
        else:
            return float('nan')

def mode_choice_model(all_trips_df):
#     all_trips_df['mode']=random.choices(['walk', 'drive'], k=len(all_trips_df))
    # all_trips_df['route_distance']=all_trips_df.apply(lambda row: m.route_lengths[row['from_possible_nodes_drive'][0]]
    #                                                                     [row['to_possible_nodes_drive'][0]], axis=1)
    options=['drive', 'cycle', 'walk', 'pt']
    all_trips_df['mode']='drive'
    ind_u1500=all_trips_df['route_distance']<1500
    ind_1500_5k=((all_trips_df['route_distance']>1500)&(all_trips_df['route_distance']<5000))
    ind_5k_10k=((all_trips_df['route_distance']>5000)&(all_trips_df['route_distance']<10000))
    ind_10kplus=all_trips_df['route_distance']>10000

    all_trips_df.loc[ind_u1500==True, 'mode']=np.random.choice(
        options, size=sum(ind_u1500), 
        replace=True, p=[0.1, 0.2, 0.6, 0.1])

    all_trips_df.loc[ind_1500_5k==True, 'mode']=np.random.choice(
        options, size=sum(ind_1500_5k), 
        replace=True, p=[0.25, 0.15, 0.3, 0.3])

    all_trips_df.loc[ind_5k_10k==True, 'mode']=np.random.choice(
        options, size=sum(ind_5k_10k), 
        replace=True, p=[0.6, 0.04, 0.01, 0.35])

    all_trips_df.loc[ind_10kplus==True, 'mode']=np.random.choice(
        options, size=sum(ind_10kplus), 
        replace=True, p=[0.9, 0.01, 0.01, 0.08])
    return all_trips_df

class Mobility_indicator(Indicator):
    def setup(self, zones, geogrid, table_name, simpop_df, external_hw_tags=["motorway","motorway_link","trunk","trunk_link"]):
        self.N_max=250
        self.geogrid=geogrid
        self.external_hw_tags=external_hw_tags
        self.zones=zones
        self.table_name=table_name
        self.base_simpop_df=simpop_df.copy()
        self.build_mobsys() 
        print('Init simulation')
        grid_zones=self.create_grid_zones()
        model_zones=zones.loc[zones['model_area']]
        model_zones['grid_area']=False
        combined_zones=model_zones.append(grid_zones)
        sim_geoids=list(zones.loc[zones['sim_area']].index)+list(grid_zones.index)
        self.sim=OpenCity.Simulation(simpop_df, self.mob_sys, combined_zones, sim_geoids=sim_geoids)
        self.sim.set_choice_models(mode_chooser=mode_choice_model)
        self.create_zone_dist_mat()

    def create_grid_zones(self):
        grid_zones=self.geogrid.copy()
        centroids=grid_zones['geometry'].centroid
        grid_zones['x_centroid']=[c.x for c in centroids]
        grid_zones['y_centroid']=[c.y for c in centroids]
        cols_to_init=[col for col in self.zones.columns if (
            ('emp_' in col) or ('res_' in col))]
        for col in cols_to_init:
            grid_zones[col]=0
        for area_col in ['model_area', 'sim_area', 'grid_area']:
            grid_zones[area_col]=1
        return grid_zones
        
    def build_mobsys(self):
        print('Building Mobility System')
        print('\t getting graph')
        G_drive_sim = self.get_graph_buffered_to_hw_type(self.zones.loc[self.zones['sim_area']], 
                                                       self.external_hw_tags, 'drive')
        print('getting external roads')
        G_drive_model=osmnx.graph_from_polygon(self.zones.loc[self.zones['model_area']].unary_union,
                                 network_type='drive', custom_filter='["highway"~"{}"]'.format('|'.join(self.external_hw_tags)))
        G_drive_combined=osmnx.graph.nx.compose(G_drive_model,G_drive_sim)
#         for edge in list(G_walk.edges):
#             G_walk.edges[edge]['speed_kph']=4.8
        G_drive_combined=osmnx.add_edge_speeds(G_drive_combined)
        G_drive_combined=PreCompOsmNet.simplify_network(G_drive_combined)
        G_drive_combined=osmnx.add_edge_travel_times(G_drive_combined)
        
        for edge in list(G_drive_combined.edges):
            G_drive_combined.edges[edge]['travel_time_walk']=G_drive_combined.edges[edge]['length']/(4800/3600)
            G_drive_combined.edges[edge]['travel_time_cycle']=G_drive_combined.edges[edge]['length']/(14000/3600)
            G_drive_combined.edges[edge]['travel_time_pt']=G_drive_combined.edges[edge]['length']/(25000/3600)

        G_drive_combined=osmnx.utils_graph.get_undirected(G_drive_combined)
#         fw_pred_drive=PreCompOsmNet.pre_compute_paths(G_drive_combined)

        # Note: this approach will assume the same route for each mode but different travel times
        # For different routes, would need to compute the fw result for each mode
        fw_all, self.route_lengths=PreCompOsmNet.pre_compute_paths(G_drive_combined, 
                                                         weight_metric='length', 
                                                         save_route_costs=True)
        pre_comp_drive=PreCompOsmNet.PreCompOSMNet(G_drive_combined, fw_all)
        networks={'drive': pre_comp_drive, 'walk': pre_comp_drive, 'cycle': pre_comp_drive, 'pt': pre_comp_drive}
        
        drive_dict={
            'target_network_id': 'drive',
            'travel_time_metric': 'travel_time'}
        walk_dict={
            'target_network_id': 'walk',
            'travel_time_metric': 'travel_time_walk'}
        cycle_dict={
            'target_network_id': 'cycle',
            'travel_time_metric': 'travel_time_cycle'}
        pt_dict={
            'target_network_id': 'pt',
            'travel_time_metric': 'travel_time_pt'}
        
        modes={'drive': OpenCity.Mode(drive_dict), 'walk': OpenCity.Mode(walk_dict),
        'cycle': OpenCity.Mode(cycle_dict), 'pt': OpenCity.Mode(pt_dict)}

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
        for buffer in [i*200 for i in range(1, 10)]:
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
        profile_id_map={'low':"0", 'mid':"1",  'high':"2"}
        
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
            coords=[[int(c[0]*1e5)/1e5, int(c[1]*1e5)/1e5] for c in row['attributes']['coordinates']]
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
        cols_to_zero= [col for col in self.sim.zones.columns if (
            ('emp_' in col) or ('res_' in col))]
        self.sim.zones.loc[self.sim.zones.grid_area==True, cols_to_zero]=0
        side_length=geogrid_data.get_geogrid_props()['header']['cellSize']
        type_def=geogrid_data.get_type_info()
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            type_info=type_def[name]
            if not name =='None':
                height=cell['height']
                cell_area=side_length*side_length
                if isinstance(height, list):
                    height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:
                    sqm_pperson=50
                total_capacity=height*cell_area/sqm_pperson
                # update where the residential capacity exists
                if 'Residential' in name:
                    self.sim.zones.loc[i_c, 'res_total']=total_capacity
                    if name=='Residential Low Income':
                        self.sim.zones.loc[i_c, 'res_income_low']=total_capacity
                    elif name=='Residential Med Income':
                        self.sim.zones.loc[i_c, 'res_income_mid']=total_capacity
                    elif name=='Residential High Income':
                        self.sim.zones.loc[i_c, 'res_income_high']=total_capacity
        for i_c, cell in enumerate(geogrid_data):
            name=cell['name']
            type_info=type_def[name]
            if not name =='None':
                height=cell['height']
                cell_area=side_length*side_length
                if isinstance(height, list):
                    height=height[-1]
                if 'sqm_pperson' in type_info:
                    sqm_pperson=type_info['sqm_pperson']
                else:
                    sqm_pperson=50
                total_capacity=height*cell_area/sqm_pperson                    
                # update the agents
                if type_info['NAICS'] is not None:
                    workers={code: int(type_info['NAICS'][code]*total_capacity) for code in  type_info['NAICS']}   
                    for code in workers:
                        home_locations=self.sample_home_locations(i_c, 'high', n=workers[code])
                        for i_w in range(workers[code]):
                            new_simpop.append({'work_geoid': i_c,'home_geoid': home_locations[i_w],
                                               'naics': code, 'earnings': 'high',
                                              'age': '30-54'})

        return new_simpop
            
    def sample_home_locations(self, work_geoid, earnings, n):
        attraction=self.sim.zones['res_income_{}'.format(earnings)]
        impedance=[self.dist_mat[hid][work_geoid] for hid in self.sim.zones.index]
        # weights=np.divide(attraction,np.array(impedance))
        weights=np.array(attraction)
        return np.random.choice(
            self.sim.zones.index, replace=True, p=weights/sum(weights), size=n)
        
    def simulate(self, simpop_df):
        print('Schedules and locations')
        simpop_df=self.sim.create_simple_HWH_schedules(simpop_df)
        print('Trip table')
        all_trips_df=self.sim.create_trip_table(simpop_df)
        all_trips_df['route_distance']=all_trips_df.apply(lambda row: self.route_lengths[row['from_possible_nodes_drive'][0]]
                                                                    [row['to_possible_nodes_drive'][0]], axis=1)
        all_trips_df=self.sim.mode_chooser(all_trips_df)
        print('Route table')
        route_table=self.sim.get_routes_table(all_trips_df)
        print('DeckGL')
        deckgl_trips=self.routes_to_deckgl_trip(route_table)
        by_mode=route_table.groupby('mode').size()
        ind=(by_mode['cycle']+by_mode['walk'])/by_mode.sum()
        return deckgl_trips, ind
    
    def return_indicator(self, geogrid_data):
        print('Starting MM Update')
        new_simpop=self.geogrid_updates(geogrid_data)
        new_simpop_df=pd.DataFrame(new_simpop)
        combined_simpop=self.base_simpop_df.append(new_simpop_df)
        sample_simpop_df=combined_simpop.sample(min(self.N_max, len(combined_simpop)))
        deckgl_trips, ind=self.simulate(sample_simpop_df)
        self.post_trips(deckgl_trips)
        print('Finished MM Update')
        return {'name': 'Active Mobility',
                           'raw_value': ind,
                           'value': ind,
                           'ref_value': 0.1,
                           'viz_type':'radar'}
    
    def post_trips(self, deckgl_trips):
        post_url='https://cityio.media.mit.edu/api/table/'+self.table_name
        r = requests.post(post_url+'/ABM2', data = json.dumps(deckgl_trips),
        	headers={'Content-Type': 'application/json'})
        print('Post ABM: {}'.format(r))


if __name__ == "__main__":
    table_name=sys.argv[1]

    # Load the saved data
    geogrid=gpd.read_file('tables/{}/geogrid.geojson'.format(table_name))
    zones=gpd.read_file('tables/{}/zones.geojson'.format(table_name))
    zones['GEOID']=zones['GEOID'].astype(int)
    zones=zones.set_index('GEOID')
    simpop_df=pd.read_csv('tables/{}/simpop_df.csv'.format(table_name))


    H=Handler(table_name=table_name)
    H.reset_geogrid_data()

    d=Density_Indicator(zones=zones)
    p=Proximity_Indicator(zones=zones, geogrid=geogrid)
    m=Mobility_indicator(zones, geogrid, table_name, simpop_df)
    H.add_indicators([d, p, m])

    H.listen()