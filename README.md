# OpenCity
Easily create mobility simulations for any location in the USA using open data. 

## Usage
### Build geometry for state
Specify
- the state where the simulation will take place 
- the geometry type ('block_group' (default) or 'block')
- the year to use for the geometry and commuting data
- Optionally: a center coordinate and radius  to restrict the analysis to

```
fips=25
year=2017
geom_type='block'
center_x_y=[-71.059763, 42.359378]
radius=5000
state=US_State(state_fips=fips, year=year, geom_type='block')
state.get_geometry(center_x_y=center_x_y, radius=radius)
```
The state geometry can be accessed as a [geopandas](https://geopandas.org/) GeoDataFrame
```
zones=state.get_geometries_included()
```

### Get commuting data for state from LEHD
```
state.get_lodes_data()
```
The commuting data can be used to build a simulated population
```
simpop_df=state.lodes_to_pop_table()
```

### Build mobility system for same area
Get the road network(s). These are used to create [pandana](https://github.com/UDST/pandana) network objects.
```
networks={}

bbox=state.get_bounds_included()
drive_nodes_df,drive_edges_df=osmnet.load.network_from_bbox(lat_min=bbox[1], lng_min=bbox[0], lat_max=bbox[3], 
                          lng_max=bbox[2], bbox=None, network_type='drive', 
                          two_way=True, timeout=180, 
                          custom_osm_filter=None)
drive_net=pandana.Network(drive_nodes_df["x"], drive_nodes_df["y"], drive_edges_df["from"], drive_edges_df["to"],
                 drive_edges_df[["distance"]])

networks['drive']=drive_net
 ```
Define modes of transportation
```
drive_dict={
    'target_network_id': 'drive',
    'weight_metric': 'distance'}
```

Create mobility system using the pandana network(s) and mode definition(s)

```
mob_sys=MobilitySystem(modes=modes,
                      networks=networks)
```
### Create the Simulation Model
Optionally define "sim_geoids", a list of geoids to which the simulation will be restricted. Only individuals who live and/or work in one of these zones will be simulated. In order to reduct computational burden, the population can be further restricted to a maximum by using the "sample_N" parameter
```
sim=Simulation(simpop_df, mob_sys, zones, sim_geoids=sim_geoids)
simpop_df=sim.get_simpop_subset(simpop_df, sample_N=10000)
```

### Simulate trips and trajectories
simpop_df=sim.create_activity_schedules(simpop_df)
all_trips_df=sim.create_trip_table(attributes)
route_table=sim.get_routes_table(all_trips_df)

### Visualising Outputs
Get resulting trips as a GeoDataFrame of 'LineString's
```
route_gdf=sim.route_table_to_geo(route_table)
```
Get resulting trips as a geojsoon compatible with the [kepler.gl](https://kepler.gl/) [Trips Layer](https://deck.gl/docs/api-reference/geo-layers/trips-layer). 'start_day_time_stamp' should be the timestamp in epoch seconds format at midnight at the beginning of the day being simulated.

```
geo_dict=sim.route_gdf_to_trips_geojson(route_gdf, start_day_time_stamp)
```

