#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:32:27 2020

@author: doorleyr
"""

import pandas as pd
import geopandas as gpd
import numpy as np
#from numpy.random import choice
import us
import math
import urllib.request as ur
from gzip import GzipFile
import random
#
#def add_centroid_xy_to_gdf(row):
#    centroid=row['geometry'].centroid
#    row['x_centroid']=centroid.x
#    row['y_centroid']=centroid.y
#    return row

def get_haversine_distance(point_1, point_2):
    """
    Calculate the distance between any 2 points on earth given as [lon, lat]
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [point_1[0], point_1[1], 
                                                point_2[0], point_2[1]])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a)) 
    r = 6371000 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

        
class US_State():
    def __init__(self, state_fips, year=2017, geom_type='block_group'):
        """
        geom_type should be one of ['block_group', 'block']
        """
        self.state_fips=str(state_fips).zfill(2)
        self.state=us.states.lookup(self.state_fips)
        self.geom_type=geom_type
        self.year=year
    def get_geometry(self, center_x_y=None, 
                     radius=None, target_crs="EPSG:4326"):
        """
        if center_lon_lat and radius are given, the returned geometry will be
        a subset of the zones in the state
        radius should be specified in meters
        """
        
        print('Getting geometry ({}) for state: {}'.format(self.geom_type, self.state.name))
        if self.geom_type=='block_group':
            self.geom=gpd.read_file('https://www2.census.gov/geo/tiger/TIGER{}/BG/tl_{}_{}_bg.zip'.format(
                    self.year, self.year, self.state_fips))
            self.geom=self.geom.set_index('GEOID')
        elif self.geom_type=='block': 
            self.geom=gpd.read_file('https://www2.census.gov/geo/tiger/TIGER{}/TABBLOCK/tl_{}_{}_tabblock10.zip'.format(
                    self.year, self.year, self.state_fips))
            self.geom=self.geom.set_index('GEOID10')
        else:
            print('Unrecoginised geometry type: {}'.format(self.geom_type))
        self.geoid_include=list(self.geom.index)
        self.geom=self.geom.to_crs(target_crs)
        self.geom['centroid']=self.geom['geometry'].centroid
        self.geom['x_centroid']=self.geom.apply(lambda row: row['centroid'].x, axis=1)
        self.geom['y_centroid']=self.geom.apply(lambda row: row['centroid'].y, axis=1)
        if ((center_x_y is not None) and (radius is not None)):
            print('\t Subsetting zones by distance')
            if target_crs=="EPSG:4326":
                dist_from_centre=self.geom.apply(lambda row: get_haversine_distance(
                        [row['x_centroid'], row['y_centroid']], center_x_y), axis=1)
            elif self.geom.crs.axis_info[0] == 'metre':
                print('Warning: distance calculation with projected coordinate system untested')
                dist_from_centre=np.sqrt(np.power(self.geom['x_centroid']-center_x_y[0], 2)+ 
                                         np.power(self.geom['y_centroid']-center_x_y[1], 2))
            self.geoid_include=list(self.geom.index[dist_from_centre<=radius].values)
            if len(self.geoid_include)==0:
                print('\t Warning: resulting geometry is empty')
                
    def get_bounds_included(self):
        bounds_all = self.geom.loc[self.geoid_include].total_bounds
        return bounds_all
    
    def get_geometries_included(self):
        subset_geom=self.geom.loc[self.geoid_include]
        return subset_geom
                
    def subset_geometries(self, included):
        """
        trip_end should be one of ['home', 'work']
        included should be a list of geometries
        
        """
        self.geoid_include=included
                
    def get_lodes_data(self, include=['wac', 'rac', 'od']):
        if 'wac' in include:
            self.wac=self.get_wac_data()
        if 'rac' in include:
            self.rac=self.get_rac_data()  
        if 'od' in include:
            self.od=self.get_od_data()

            
    def get_od_data(self):
        print('Getting OD data from https://lehd.ces.census.gov/data/lodes/LODES7/')
        req = ur.Request('https://lehd.ces.census.gov/data/lodes/LODES7/{}/od/{}_od_main_JT00_{}.csv.gz'.format(
                self.state.abbr.lower(), self.state.abbr.lower(), self.year)) 
        z_f = ur.urlopen(req)
        f = GzipFile(fileobj=z_f, mode="r")
        od = pd.read_csv(f)
        print('\t Formatting OD data')
        od=self.format_lodes_data(od)
        return od

    def get_rac_data(self):
        print('Getting RAC data from https://lehd.ces.census.gov/data/lodes/LODES7/')
        req = ur.Request('https://lehd.ces.census.gov/data/lodes/LODES7/{}/rac/{}_rac_S000_JT00_{}.csv.gz'.format(
                self.state.abbr.lower(), self.state.abbr.lower(), self.year)) 
        z_f = ur.urlopen(req)
        f = GzipFile(fileobj=z_f, mode="r")
        rac = pd.read_csv(f)
        print('\t Formatting RAC data')
        rac=self.format_lodes_data(rac)
        return rac
    
    def get_wac_data(self):
        print('Getting WAC data from https://lehd.ces.census.gov/data/lodes/LODES7/')
        req = ur.Request('https://lehd.ces.census.gov/data/lodes/LODES7/{}/wac/{}_wac_S000_JT00_{}.csv.gz'.format(
                self.state.abbr.lower(), self.state.abbr.lower(), self.year)) 
        z_f = ur.urlopen(req)
        f = GzipFile(fileobj=z_f, mode="r")
        wac = pd.read_csv(f)
        print('\t Formatting WAC data')
        wac=self.format_lodes_data(wac)
        return wac
        
    def format_block_id(self, block_id):
        return str(block_id).zfill(15)
        
    def format_lodes_data(self, block_df):
        block_cols=[c for c in ['h_geocode', 'w_geocode'] if c in block_df.columns]
        for col in block_cols:
            block_df[col]=block_df.apply(lambda row: self.format_block_id(row[col]), axis=1)
        if self.geom_type=='block':
            # use the existing geoid, just rename it
            cols_to_rename={col: col.replace('geocode', 'geoid') for col in block_cols}
            block_df=block_df.rename(columns=cols_to_rename)
#            for col in block_cols:
#                new_block_col=col.split('_')[0]+'_geoid'
#                block_df[new_block_col]=block_df[col]
            block_df=block_df.set_index(list(cols_to_rename.values()), drop=False)
            return block_df
        elif self.geom_type=='block_group':
            # aggregate blocks to block groups
            cols_not_to_sum=block_cols +['createdate']
            cols_to_sum=[col for col in block_df.columns if not col in cols_not_to_sum]
            block_group_cols=[]
            for col in block_cols:
                new_bg_col=col.split('_')[0]+'_geoid'
                block_df[new_bg_col]=block_df.apply(lambda row: row[col][0:12], axis=1)
                block_group_cols.append(new_bg_col)
            bg_df=block_df.groupby(block_group_cols, as_index=False)[cols_to_sum].agg('sum')
            bg_df=bg_df.set_index(block_group_cols, drop=False)
            return bg_df
        else:
            print('Geometry {} not recognised'.format(self.geom_type))
            
    def lodes_to_pop_table(self):
        od_subset=self.od
        od_subset=od_subset.loc[((od_subset['h_geoid'].isin(self.geoid_include))&
                                 (od_subset['w_geoid'].isin(self.geoid_include)))]

        print('Using {} of {} rows in OD data'.format(len(od_subset), len(self.od)))
        # convert counts by attribute to probabilities for each attribute
        attr_cols=['S000', 'SA01', 'SA02', 'SA03', 'SE01', 'SE02','SE03', 'SI01', 'SI02', 'SI03']
        probs_all=od_subset[attr_cols].div(od_subset.S000, axis=0)
        probs_all=probs_all.rename(columns={col: 'prob_'+col for col in attr_cols})
        od_subset_probs=pd.concat([od_subset, probs_all], axis=1)
#        for each row, sample the appropriate number of people
        sampled_age, sampled_earning, sampled_industry, sampled_home, sampled_work= [], [],[],[],[]
        count=0
        for ind, row in od_subset_probs.iterrows():
            if count%10000==0:
                print('{} of {}'.format(count, len(od_subset_probs)))
            count+=1
            n=row['S000']
            age_samples=np.random.choice(a=['u30','30to54','55plus'], size=n, 
                                         p=row[['prob_SA01','prob_SA02','prob_SA03']])
            earn_samples=np.random.choice(a=['u1250','1250to3333','3333plus'], size=n, 
                                          p=row[['prob_SE01','prob_SE02','prob_SE03']])
            indust_samples=np.random.choice(a=['goods_prod','trade_trans_util','other'], 
                                            size=n, p=row[['prob_SI01','prob_SI02','prob_SI03']])
        #     sampled_persons=np.column_stack([age_samples, earn_samples, indust_samples, 
        #                                      [row['h_geoid']]*n, [row['w_geoid']]*n ])
            sampled_age.extend(age_samples)
            sampled_earning.extend(earn_samples)
            sampled_industry.extend(indust_samples)
            sampled_home.extend([row['h_geoid']]*n)
            sampled_work.extend([row['w_geoid']]*n)
        self.simpop_df=pd.DataFrame.from_dict({'age':sampled_age, 'earnings': sampled_earning, 'industry': sampled_industry, 'home_geoid': sampled_home,
                                                   'work_geoid': sampled_work}, orient='columns')
           
class MobilitySystem():
    def __init__(self, modes, networks):
        """
        modes should be a dictionary with the identifier as key 
        and the mode object as value
        networks should be a dictionary where keys are network ids and values 
        are pandana network objects
        """
        assert all(modes[m].target_network_id in networks.keys() for m in modes), "A network ID doesn't exist"
        self.modes = modes
        self.networks = networks
        
    def update_modes(self,modes):
        """
        add new networks and/or update existing modes
        modes should be a dictionary with the identifier as key 
        and the mode object as value
        """
        self.modes.update(modes)
        
    def remove_mode(self,mode_id):
        del self.modes[mode_id]
        
    def get_mode_by_id(self, mode_id):
        return self.modes[mode_id]
    
    def update_networks(self, networks):
        """
        add new networks and/or update existing networks
        networks should be a dictionary where keys are network ids and values 
        are pandana network objects
        """
        self.modes.update(networks)
        
    def remove_network(self, network_id):
        del self.networks[network_id]
        
    def get_network_by_id(self, network_id):
        return self.networks[network_id]
    
    def get_one_route(self, mode_id, node_a, node_b, include_coordinates=False):
        target_network_id=self.modes[mode_id].target_network_id
        weight_metric=self.modes[mode_id].weight_metric
        node_path= self.networks[target_network_id].shortest_path(
                node_a, node_b, imp_name=weight_metric)
        output={'node_path': node_path}
        if include_coordinates:
            output['coordinates']=self.networks[target_network_id].nodes_df.loc[
                    node_path, ['x', 'y']].values
        return output       
    
    def get_many_routes(self, mode_id, nodes_a, nodes_b, include_coordinates=False):
        target_network_id=self.modes[mode_id].target_network_id
        weight_metric=self.modes[mode_id].weight_metric
        node_paths= self.networks[target_network_id].shortest_paths(
                nodes_a, nodes_b, imp_name=weight_metric)
        output={'node_paths': node_paths}
        if include_coordinates:
            output['coordinates']=[self.networks[target_network_id].nodes_df.loc[
                    np, ['x', 'y']].values for np in node_paths]
        return output 
    # TODO: add function to download OSM networks

class Mode():
    def __init__(self, mode_dict):
        """
        mode_dict must include at minimum:
            network_id: the network this mode uses
            weight_metric: the relevant weight metric of the target network (eg. distance)
        
        """
        self.mode_dict=mode_dict
        self.target_network_id=mode_dict['target_network_id']
        self.weight_metric=mode_dict['weight_metric']

class Zones():
    def __init__(self, zones_gdf, geoid_col=None):
        """
        If the index if not the geoid, then the geoid_col should be specified
        """       
        # TODO: support geojson?
        zones_gdf_centroid=zones_gdf["geometry"].centroid
        zones_gdf['x_centroid']=[c.x for c in zones_gdf_centroid]
        zones_gdf['y_centroid']=[c.y for c in zones_gdf_centroid]
        if geoid_col is not None:
            zones_gdf=zones_gdf.set_index(geoid_col)
        self.gdf=zones_gdf
    
        
#class SimPop:
#    def __init__(self, simpop_df):
#        """
#        pop_df should be a pandas df including columns home_geoid and work_geoid
#        """
#        self.simpop_df=simpop_df

        
    
class Simulation():
    def __init__(self, sim_pop, mobsys, zones, sample_N=None, sim_geoids=None):
        self.check_zones(sim_pop, zones)           
        sim_pop=sim_pop.reset_index(drop=True)
        self.mobsys=mobsys
        self.zones=zones
        if sim_geoids is not None:
            sim_pop=sim_pop.loc[((sim_pop['home_geoid'].isin(sim_geoids))|
                          (sim_pop['work_geoid'].isin(sim_geoids)))]
        if sample_N is not None:
            if len(sim_pop)>sample_N:
                sim_pop=sim_pop.sample(n=sample_N)
        self.sim_pop=sim_pop
        self.get_close_node_table()
        
    def check_zones(self, sim_pop, zones):
        all_zones_sim_pop=set(list(sim_pop['home_geoid'])+list(sim_pop['work_geoid']))
        assert all(zone in zones.index for zone in all_zones_sim_pop)
        
    def get_close_node_table(self):
        """
        Creates a lookup table for each network, mapping zone geoids to closest node ids
        returns an object where each key is a network id and value is a pandas dataframe,
        indexed by zone geoids
        """
        print('Finding closest nodes to every zone centroid')
        close_nodes_obj={}
        for network_id in self.mobsys.networks:
            close_nodes_this_net=pd.DataFrame(index=self.zones.index,
                                              data=self.mobsys.networks[network_id].get_node_ids(
                                                      x_col=self.zones['x_centroid'],
                                                      y_col=self.zones['y_centroid']))
            close_nodes_obj[network_id]=close_nodes_this_net            
        self.close_nodes=close_nodes_obj
        
    def default_scheduler(self, sim_pop_row):
        sim_pop_row['activities']=['H', 'W', 'H']
        sim_pop_row['start_times']=[0, 3600*7, 3600*18]
        return sim_pop_row
    
    def default_location_chooser(self, sim_pop_row, zones):
        locations=[]
        for activity in sim_pop_row['activities']:
            if activity=='H':
                locations.append(sim_pop_row['home_geoid'])
            elif activity=='W':
                locations.append(sim_pop_row['work_geoid'])
            else:
                locations.append(random.choice(zones.index))
        sim_pop_row['locations']=locations
        return sim_pop_row
        
    def set_activity_scheduler(self, scheduler=None, location_chooser=None):
        """
        If no activity scheduler is specified, a simple home-work-home 
        scheduler is used
        """
        if scheduler==None:
            self.scheduler=self.default_scheduler
        else:
            self.scheduler=scheduler
        if location_chooser==None:
            self.location_chooser=self.default_location_chooser
        else:
            self.location_chooser=location_chooser
            
    def create_activity_schedules(self):
        if self.scheduler is None:
            self.set_activity_scheduler()
        self.sim_pop=self.sim_pop.apply(lambda row: self.scheduler(row), axis=1)
        self.sim_pop=self.sim_pop.apply(lambda row: self.location_chooser(row, self.zones), axis=1)
        
    def create_trip_table(self):
        """
        For every person who passes through the simulation area:
            add all their trips to a trip table
        For every trip:
            find route
            add a deckgl path
        
        """
        
    

        
#        
#class Network():
#    def __int__(nodes_df, links_df,
#                node_x_col, node_y_col,
#                edge_from_col, edge_to_col, 
#                attributes):
#        self.pdna_net=.network.Network(node_x, node_y, edge_from, edge_to, edge_weights, twoway=True)
#        """
#        mode_dict should include at minimum the id, the 
#        """
#        
        
    
        
        