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
from shapely.geometry import Point, LineString

def get_haversine_distance(point_1, point_2):
    # TODO: vectorise for fast calculation of multiple queries
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

def coord_list_to_line_string(coordinates):
    if len(coordinates)>1:
        return LineString(coordinates)
    elif len(coordinates)>0:
        return Point(coordinates)
    else:
        return None
       
class US_State():
    def __init__(self, state_fips, year=2017, geom_type='block_group'):
        """
        geom_type should be one of ['block_group', 'block']
        """
        self.state_fips=str(state_fips).zfill(2)
        self.state=us.states.lookup(self.state_fips)
        self.geom_type=geom_type
        self.year=year
    def get_geometry(self, centre_x_y=None, 
                     radius=None, target_crs="EPSG:4326"):
        """
        if centre_lon_lat and radius are given, the returned geometry will be
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
        if ((centre_x_y is not None) and (radius is not None)):
            print('\t Subsetting zones by distance')
            if target_crs=="EPSG:4326":
                dist_from_centre=self.geom.apply(lambda row: get_haversine_distance(
                        [row['x_centroid'], row['y_centroid']], centre_x_y), axis=1)
            elif self.geom.crs.axis_info[0] == 'metre':
                print('Warning: distance calculation with projected coordinate system untested')
                dist_from_centre=np.sqrt(np.power(self.geom['x_centroid']-centre_x_y[0], 2)+ 
                                         np.power(self.geom['y_centroid']-centre_x_y[1], 2))
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
        return pd.DataFrame.from_dict({'age':sampled_age, 'earnings': sampled_earning, 'industry': sampled_industry, 'home_geoid': sampled_home,
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
        self.build_link_attribute_lookup()
        self.edge_ordering={}
        for network in networks:
            self.edge_ordering[network]=list(networks[network].edges_df.index)
        
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
        self.networks.update(networks)
        self.build_link_attribute_lookup()
        
    def remove_network(self, network_id):
        del self.networks[network_id]
        del self.attr_lookups[network_id]
        
    def get_network_by_id(self, network_id):
        return self.networks[network_id]

    def build_link_attribute_lookup(self):
        print("building link attribute lookup")
        self.attr_lookups={}
        for network_id in self.networks:
            edges_df=self.networks[network_id].edges_df.copy()
            if self.networks[network_id]._twoway:
                temp=edges_df.set_index(['to', 'from'], drop=False)
                temp['from']=list(edges_df['to'].values)
                temp['to']=list(edges_df['from'].values)
                edges_df=edges_df.append(temp)
                # eliminate duplicates, keeping the lowest of each weight metric
                # in general other attributes are the same for each duplicate
                # note that for other attributes which differ across links with the same a and b nodes (like link type)
                # this approach may give incorrect results.
            weight_metrics= self.networks[network_id].impedance_names
            agg_dict={weight: 'min' for weight in weight_metrics}
            for col in edges_df.columns:
                if col not in weight_metrics:
                    agg_dict.update({col: lambda x: x.iloc[0]})
            edges_df=edges_df.groupby(edges_df.index).agg(agg_dict)
            # add node coordinates
            edges_df=edges_df.join(self.networks[network_id].nodes_df, how='left', on='from').rename(
                columns={'x': 'a_node_x', 'y': 'a_node_y'})
            edges_df=edges_df.join(self.networks[network_id].nodes_df, how='left', on='to').rename(
                columns={'x': 'b_node_x', 'y': 'b_node_y'})
            self.attr_lookups[network_id]=edges_df.to_dict(orient='index')

    def get_one_route(self, mode_id, node_a, node_b, include_attributes=False):
        target_network_id=self.modes[mode_id].target_network_id
        weight_metric=self.modes[mode_id].weight_metric
        node_path= self.networks[target_network_id].shortest_path(
                node_a, node_b, imp_name=weight_metric)
        if include_attributes==True:
            attributes=self.get_path_link_attributes(node_path, target_network_id)
            return {'node_path': node_path, 'attributes': attributes}
        return node_path       
    
    def get_many_routes(self, mode_id, nodes_a, nodes_b, include_attributes=False):
        target_network_id=self.modes[mode_id].target_network_id
        weight_metric=self.modes[mode_id].weight_metric
        node_paths= self.networks[target_network_id].shortest_paths(
                nodes_a, nodes_b, imp_name=weight_metric)
        if include_attributes==True:
            attributes=[self.get_path_link_attributes(np, target_network_id) for np in node_paths]
            return {'node_path': node_paths, 'attributes': attributes}
        return node_paths 
    
    def get_path_link_attributes(self, path, network_id, attribute_names=None):
        if attribute_names==None:
            attribute_names=self.networks[network_id].impedance_names
        output={attr: [] for attr in attribute_names}
        if len(path)>1:
            coordinates=[]
            edges=[]
            for node_ind in range(len(path)-1):
                from_node=path[node_ind]
                to_node=path[node_ind+1]
                edges.append((from_node, to_node))
                this_link_attrs=self.attr_lookups[network_id][(from_node, to_node)]
                for attr in attribute_names:
                    output[attr].append(this_link_attrs[attr])
                coordinates+=[[this_link_attrs['a_node_x'], this_link_attrs['a_node_y']]]
            coordinates+=[[this_link_attrs['b_node_x'], this_link_attrs['b_node_y']]]
            output['coordinates']=coordinates
            output['edges']=edges
            return output
        else:
            output['coordinates']=[]
            output['edges']=[]
            return output
            

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

class Simulation():
    def __init__(self, sim_pop, mob_sys, zones,  sim_geoids=None, person_attributes=None,
                centre_x_y=None, vis_radius=None):
        self.scheduler=self.default_scheduler
        self.location_chooser=self.default_location_chooser
        self.mode_chooser=self.default_mode_chooser
        self.check_zones(sim_pop, zones)           
        self.mob_sys=mob_sys
        self.zones=zones
        self.sim_geoids=sim_geoids
        self.centre_x_y=centre_x_y
        self.vis_radius=vis_radius
        self.get_close_node_table()
        if person_attributes==None:
            person_attributes=[col for col in sim_pop.columns if col not in [
                        'home_geoid', 'work_geoid']]
        self.person_attributes=person_attributes
        if ((centre_x_y is not None) and (vis_radius is not None)):
            self.get_vis_nodes()
        
    def get_vis_nodes(self):
        print('Finding visible nodes')
        self.vis_nodes={}
        for net in self.mob_sys.networks:
            self.vis_nodes[net]=[ind for ind, row in self.mob_sys.networks[net].nodes_df.iterrows(
            ) if get_haversine_distance([row['x'], row['y']], self.centre_x_y)<self.vis_radius]
  
    def get_simpop_subset(self, sim_pop, sample_N=None):
        sim_pop=sim_pop.reset_index(drop=True)
        if self.sim_geoids is not None:
            sim_pop=sim_pop.loc[((sim_pop['home_geoid'].isin(self.sim_geoids))|
                          (sim_pop['work_geoid'].isin(self.sim_geoids)))]
        if sample_N is not None:
            if len(sim_pop)>sample_N:
                sim_pop=sim_pop.sample(n=sample_N)
        return sim_pop
        
    def get_od_flows(self, sim_pop):
        return sim_pop.groupby(['home_geoid', 'work_geoid']).size()
    
    def get_zone_coords(self):
        return self.zones[['x_centroid', 'y_centroid']]
    
    def check_zones(self, sim_pop, zones):
        all_zones_sim_pop=set(list(sim_pop['home_geoid'])+list(sim_pop['work_geoid']))
        assert all(zone in zones.index for zone in all_zones_sim_pop)
        
    def get_close_node_table(self):
        """
        Creates a lookup table for each network, mapping zone geoids to closest node ids
        returns a dataframe with a column for each network, indexed by zone geoid
        """
        print('Finding closest nodes to every zone centroid')
        
        close_nodes={}
        for net in self.mob_sys.networks:
            close_nodes['node_'+net]=self.mob_sys.networks[net].get_node_ids(
                x_col=self.zones['x_centroid'],y_col=self.zones['y_centroid'])            
        self.close_nodes_df=pd.DataFrame(index=self.zones.index, data=close_nodes)
        
        
    def default_scheduler(self, sim_pop_row):
        """
        gives everyone a schedule of [home, work, other, work, home]
        """
        sim_pop_row['activities']=['H', 'W', 'O','W','H']
        durations_hrs= [6+4*random.random(), # starts work between 6am and 1pm
                   3+2*random.random(), # works for 3-5 hours
                   0.5+1*random.random(), # break between 0.5 and 1.5 hours
                   3+2*random.random()] # works for 3-5 hours
        sim_pop_row['start_times']=[0] + list(3600*np.cumsum(durations_hrs))
        return sim_pop_row
    
    def default_mode_chooser(self, all_trips_df):
        """
        Input is entire trip table- modifies and returns the table
        """
        print("Choosing modes")
        all_trips_df['mode']='drive'
#         all_trips_df['net']='drive'
        
        return all_trips_df
    
    def default_location_chooser(self, sim_pop_row, zones):
        locations=[]
        for activity in sim_pop_row['activities']:
            if activity=='H':
                locations.append(sim_pop_row['home_geoid'])
            elif activity=='W':
                locations.append(sim_pop_row['work_geoid'])
            else:
                if self.sim_geoids is None:
                    locations.append(random.choice(zones.index))
                else:
                    locations.append(random.choice(self.sim_geoids))
        sim_pop_row['locations']=locations
        return sim_pop_row
        
    def set_choice_models(self, scheduler=None, location_chooser=None, mode_chooser=None):
        """
        Replace default choice models with user-specified models
        """
        if scheduler is not None:
            self.scheduler=scheduler
        if location_chooser is not None:
            self.location_chooser=location_chooser
        if mode_chooser is not None:
            self.mode_chooser=mode_chooser
            
    def create_activity_schedules(self, sim_pop):
        print("Scheduling activities")
        sim_pop=sim_pop.apply(lambda row: self.scheduler(row), axis=1)
        print("Chhosing locations for each activity")
        sim_pop=sim_pop.apply(lambda row: self.location_chooser(row, self.zones), axis=1)
        return sim_pop
        
    def create_trip_table(self, sim_pop):
        """
        For every person who passes through the simulation area:
            add all their trips to a trip table
        
        """
        all_trips=[]
        for ind, row in sim_pop.iterrows():
            activities=row['activities']
            start_times=row['start_times']
            locations=row['locations']
            for i in range(len(activities)-1):
                trip={'from_activity': activities[i],
                      'to_activity': activities[i+1],
                      'from_zone': locations[i],
                      'to_zone': locations[i+1],
                      'person_id': row.name,
                      'start_time': start_times[i+1]}
                for attr in self.person_attributes:
                    trip[attr]=row[attr]
                all_trips.append(trip)
        all_trips_df=pd.DataFrame(all_trips)
        all_trips_df=self.mode_chooser(all_trips_df)
        all_trips_df=all_trips_df.join(self.close_nodes_df, how='left', on='from_zone').rename(columns={
                'node_'+net: 'from_node_'+net for net in self.mob_sys.networks})
        all_trips_df=all_trips_df.join(self.close_nodes_df, how='left', on='to_zone').rename(columns={
                'node_'+net: 'to_node_'+net for net in self.mob_sys.networks})
        return all_trips_df
    
    # TODO: if I'm ultimately always using the traversal tables (rather than route_table)
    # then I should rethink this process- there's no need to get attributes/coordinates for every time a link
    # is traversed if I'm then going to groupby link and just take the first instance for each.
    # just get the link ids with the person attributes. After the grouby operations, then get the
    # cordinates etc.
    # Edge weight is needed however for building the node traversal table
    # For the node traversals
    
    def get_routes_table(self, all_trips_df):
        route_table=pd.DataFrame()
        for mode_id in self.mob_sys.modes:
            target_network_id=self.mob_sys.modes[mode_id].target_network_id
            route_table_this_mode=all_trips_df.loc[((all_trips_df['mode']==mode_id)&
                                          (~all_trips_df['from_node_'+target_network_id].isnull())&
                                          (~all_trips_df['to_node_'+target_network_id].isnull()))]
            routes=self.mob_sys.get_many_routes(
                mode_id, 
                [int(n) for n in route_table_this_mode['from_node_'+target_network_id].values],
                [int(n) for n in route_table_this_mode['to_node_'+ target_network_id].values], 
                include_attributes=True)
            route_table_this_mode['node_path']=routes['node_path']
            route_table_this_mode['attributes']=routes['attributes']
            route_table=route_table.append(route_table_this_mode)
        return route_table
    
    def get_edge_traversals(self, route_table):
        """
        Create a table with a row for every time someone travserses a link
        Include edge identifier, from coordinate, to coordinate and personal attributes
        """
        route_table['n_edges']=route_table.apply(lambda row: len(row['attributes']['edges']), axis=1)
        edge_traversal={'edge_id': [], 'from_coord': [], 'to_coord': []}
        for ind, row in route_table.iterrows():
            edge_traversal['edge_id'].extend(row['attributes']['edges'].copy())
            route_coords=row['attributes']['coordinates'].copy()
            edge_traversal['from_coord'].extend(route_coords[:-1])
            edge_traversal['to_coord'].extend(route_coords[1:])
        # person-level attributes- replicated for each edge in the route
        for p_attr in (self.person_attributes+['start_time', 'mode']):
            attr_as_list=list(route_table[p_attr])
            n_edges_by_route=list(route_table['n_edges'])
            this_attribute_by_edge=[]
            for t in range(len(attr_as_list)):
                this_attribute_by_edge.extend([attr_as_list[t]]*n_edges_by_route[t])
            edge_traversal[p_attr]=this_attribute_by_edge
        return pd.DataFrame(edge_traversal) 

    def get_node_traversals(self, route_table):
        route_table['n_nodes']=route_table.apply(lambda row: len(row['node_path']), axis=1)
        node_traversal={'coord': [], 'ts': [], 'node_id':[]}
        for ind, row in route_table.iterrows():
            mode_id=row['mode']
            target_network_id=self.mob_sys.modes[mode_id].target_network_id
            weight_metric=self.mob_sys.modes[mode_id].weight_metric
            if row['n_nodes']>1:
#                 print('{}_{}_{}'.format(i_r, len(route['coordinates']), len(node_paths[i_r])))
                node_traversal['coord'].extend(row['attributes']['coordinates'].copy())
                node_traversal['node_id'].extend(row['node_path'].copy())
                trip_start=int(row['start_time'])
                timestamps=[trip_start]+list((trip_start+np.cumsum(row['attributes'][weight_metric].copy())).astype(int))
                node_traversal['ts'].extend(timestamps)
        # person-level attributes- replicated for each edge in the route
        for p_attr in (self.person_attributes+['mode', 'person_id']):
            attr_as_list=list(route_table[p_attr])
            n_nodes_by_route=list(route_table['n_nodes'])
            this_attribute_by_node=[]
            for t in range(len(attr_as_list)):
                if n_nodes_by_route[t]>1:
                    this_attribute_by_node.extend([attr_as_list[t]]*n_nodes_by_route[t])
            node_traversal[p_attr]=this_attribute_by_node
#         for col in node_traversal:
#             print('{}_{}'.format(col, len(node_traversal[col])))
        return pd.DataFrame(node_traversal)  
    
    def get_edge_traffic(self, edge_traversal_table, mode):
        edge_traversal_table['count']=1
        agg_dict={'from_coord': lambda x: x.iloc[0],
                  'to_coord': lambda x: x.iloc[0],
                  'count': sum}
        agg_edge_traffic=edge_traversal_table.groupby('edge_id').agg(agg_dict)

        if ((self.centre_x_y is not None) and (self.vis_radius is not None)):
            agg_edge_traffic['in_scope']=agg_edge_traffic.apply(lambda row:get_haversine_distance(
                row['from_coord'], self.centre_x_y)<self.vis_radius, axis=1)
            return agg_edge_traffic[agg_edge_traffic['in_scope']]
        else:
            return agg_edge_traffic
        
    def node_traversal_to_trips_json(self, node_traversal_table):
        if ((self.centre_x_y is not None) and (self.vis_radius is not None)):
#             node_traversal_table['in_scope']=node_traversal_table.apply(lambda row:get_haversine_distance(
#                 row['coord'], centre_x_y)<vis_radius, axis=1)           
#             node_traversal_table=node_traversal_table[node_traversal_table['in_scope']]
            for mode in node_traversal_table['mode'].unique():
                net=self.mob_sys.modes[mode].target_network_id
                node_traversal_table.loc[node_traversal_table['mode']==mode, 'in_scope'
                                        ]=node_traversal_table.loc[node_traversal_table['mode']==mode, 'node_id'].isin(self.vis_nodes[net])
        node_traversal_table['coord_ts']=node_traversal_table.apply(lambda row: row['coord']+[row['ts']], axis=1)
        agg_dict={'coord_ts': lambda x: list(x),
#                   'coord': lambda x: list(x),
#                   'ts': lambda x: list(x),
                  'age': lambda x: x.iloc[0],
                  'earnings': lambda x: x.iloc[0],
                  'industry': lambda x: x.iloc[0]}
        node_traversal_per_person_mode=node_traversal_table.groupby(['person_id', 'mode']).agg(agg_dict)
        trips=node_traversal_per_person_mode.to_dict(orient='records')
        return trips       
  
    def route_table_to_geo(self, route_table):    
        route_table['line_string']=route_table.apply(lambda row: 
            coord_list_to_line_string(row['attributes']['coordinates']), axis=1)
        route_gdf=gpd.GeoDataFrame(route_table, geometry='line_string')
        return route_gdf
    
    def route_gdf_to_trips_geojson(self, route_gdf, start_day_time_stamp):
        geo_dict=route_gdf.__geo_interface__
        for i_f, feat in enumerate(geo_dict['features']):
            if feat['geometry'] is not None:
                coordinates=feat['geometry']['coordinates']
                time_metric=self.mob_sys.modes[feat['properties']['mode']].weight_metric
                timestamps=[feat['properties']['start_time']]+list(feat['properties']['start_time']+60*np.cumsum(feat['properties']['attributes'][time_metric]))
                timestamps=[start_day_time_stamp+int(t) for t in timestamps]
                new_coordinates=[[c[0], c[1], 0, timestamps[i]] for i,c in enumerate(coordinates)]
                geo_dict['features'][i_f]['geometry']['coordinates']=new_coordinates
            del geo_dict['features'][i_f]['properties']['attributes']
            del geo_dict['features'][i_f]['properties']['node_path']    
        return geo_dict
                

    
        
        