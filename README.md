# CS_Spatial_Modules
Create spatial indicators and mobility simulations for CityScope.

## OpenCity
A tool for creating simple mobility simulations and vizualisations with minimal data inputs.

### Data
The data required are:
- a shapefile of geometries in the study region
- an O-D matrix, optionally stratified by demographic attributes and/or industry. eg:

| Home GEOID | Work GEOID | N_total | N_income_low | ... | N_NAICS_44-45 | ...|
| --------- | --------- | --------- | --------- | --------- | --------- | --------- |
| 1  | 1 | 30  | 15  | ...  | 5  | ...  |
| 1  | 2 | 40  | 17  | ...  | 9  | ...  |
| ... | ... | ... | ... | ... | ... |...  |

For USA-based projects the data above can be easily obtained from public data sources using the tools in this repository.

### Usage
#### Build geometry for state
Specify
- the state where the simulation will take place 
- the geometry type ('block_group' (default) or 'block')
- the year to use for the geometry and commuting data
- Optionally: a center coordinate and radius  to restrict the analysis to

```
fips=25
year=2017
geom_type='block'
centre= {'lat': 42.352927, 'lon': -71.059435}
centre_x_y=[centre['lon'], centre['lat']]

model_area_radius=5000
sim_area_radius=2000

state=OpenCity.US_State(state_fips=fips, year=year, 
#                         geom_type=geom_type
                       )
state.get_geometry()

state.subset_geom_by_distance(centre_x_y, model_area_radius, 'model_area')
state.subset_geom_by_distance(centre_x_y, sim_area_radius, 'sim_area')
```
The state geometry can be accessed as a [geopandas](https://geopandas.org/) GeoDataFrame
```
all_zones=state.return_geometry()
sim_zones=state.return_geometry('sim_area')
model_zones=state.return_geometry('model_area')
```

#### Get commuting data for state from LEHD
```
state.get_lodes_data()
```
The commuting data can be used to build a simulated population
```
simpop_df=state.lodes_to_pop_table(model_subset_name='model_area',
                                  sim_subset_name='sim_area')
```

#### Build mobility system for same area
Get the road network(s). These are used to create [pandana](https://github.com/UDST/pandana) network objects.
```
import pandana
import osmnet

networks={}

bbox=state.get_bounds(subset_name='model_area')
drive_nodes_df,drive_edges_df=osmnet.load.network_from_bbox(lat_min=bbox[1], lng_min=bbox[0], lat_max=bbox[3], 
                          lng_max=bbox[2], bbox=None, network_type='drive', 
                          two_way=True, timeout=180, 
                          custom_osm_filter=None)
drive_edges_df['travel_time']=drive_edges_df['distance']*50000/3600

drive_net=pandana.Network(drive_nodes_df["x"], drive_nodes_df["y"], drive_edges_df["from"], drive_edges_df["to"],
                 drive_edges_df[["distance", "travel_time"]])

networks['drive']=OpenCity.PdnaNetwork(drive_net)
 ```
Define modes of transportation
```
drive_dict={
    'target_network_id': 'drive',
    'travel_time_metric': 'travel_time'}


modes={'drive': OpenCity.Mode(drive_dict)}
```

Create mobility system using the pandana network(s) and mode definition(s)

```
mob_sys=MobilitySystem(modes=modes,
                      networks=networks)
```
#### Create the Simulation Model
In order to reduct computational burden, the population can be further restricted to a maximum by using the "sample_N" parameter
```
sim=OpenCity.Simulation(simpop_df, mob_sys, model_zones)
simpop_df=sim.get_simpop_subset(simpop_df, sample_N=1000)

```

#### Simulate trips and trajectories
```
all_trips_df=sim.create_trip_table(simpop_df)
all_trips_df=sim.mode_chooser(all_trips_df)
route_table=sim.get_routes_table(all_trips_df)

```

#### Visualising Outputs
Get resulting trips as a GeoDataFrame of 'LineString's
```
route_gdf=sim.route_table_to_geo(route_table)
```
Get resulting trips as a geojsoon compatible with the [kepler.gl](https://kepler.gl/) [Trips Layer](https://deck.gl/docs/api-reference/geo-layers/trips-layer). 'start_day_time_stamp' should be the timestamp in epoch seconds format at midnight at the beginning of the day being simulated.

```
start_day_time_stamp=
geo_dict=sim.route_gdf_to_trips_geojson(route_gdf, start_day_time_stamp)
```

## CityScope Indicators
A set of modules which generate interactive indicators for CityScope project. Communication with [cityIO](https://github.com/CityScope/CS_CityIO) is handled via the [cs-brix](https://github.com/CityScope/CS_Brix) library.

### Data
The data required are:
- a shapefile of zone geometries in the study region. This should include columns indicating (i) counts of people living in each area by demographic categories and (ii) counts of people employed in each area by demographic categories.


| GEOID  | res_income_low | res_income_med | ... | emp_total | emp_income_low | ... | geometry |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |------------- |
| 1  | 1000 | 500  | 300 | ...  | 600  | ...  | POLYGON(...) |
| 2  | 1200 | 600  | 200  | ...  | 550 | ...  |POLYGON(...) |
|  ...   |  ...  |  ...   |  ...   | ...  |  ...  | ...  |POLYGON(...) |

- the GEOGRID of the CityScope table (can be created using [CityScopeJS](https://cityscope.media.mit.edu/CS_cityscopeJS/))
- the simulated population (can be created from an O-D matrix using OpenCity as outlined above)

### Usage
import CS_Indicators as CS
from brix import Indicator, Handler
import geopandas as gpd
import pandas as pd

geogrid=gpd.read_file(SAVED_GEOGRID_LOCATION)
zones=gpd.read_file(SAVED_ZONES_LOCATION)
simpop_df=pd.read_csv(SAVED_SIMPOP_LOCATION)

table_name='cs_course_volpe'

H=Handler(table_name=table_name)
H.reset_geogrid_data() # Optional

d=CS.Density_Indicator(zones=zones)
p=CS.Proximity_Indicator(zones=zones, geogrid=geogrid)
m=CS.Mobility_indicator(zones, geogrid, table_name, simpop_df)
H.add_indicators([
    d, 
    p, 
    m
])

H.listen()



