import xarray
import numpy as np

# Set up an ELM domain with polygons, based on Fengming's dataset from Wang et al 2019 JAMES paper
# 'This is for 7 land-units summarized from NGEE Site Area C with additional Arctic PFTS: 1-trough, 2-LCPcenter, 3-LCPtransition, 4-HCPcenter, 5-HCPtransition, 6-Rim, 7-average

landcover_types=[
  'trough',
  'LCPcenter',
  'LCPtransition',
  'HCPcenter',
  'HCPtransition',
  'Rim',
  'average',
]

# For now we are setting up each land cover type as a separate grid cell. 
# In the future we could do some of this using topo units within grid cells

# Hydrology uses the coastal wetland configuration setup for specifying a time series of hydrological boundary condition
# Hourly time series, just do one year and the model will repeat it. Currently using water level of zero. 
# Surface data set below will define ground surface height above drainage to give different hydrological conditions
ntimes=365*24
# Two grid cells (trough and high centered polygon)
num_grids=len(landcover_types)

tide_data_multicell=xarray.Dataset(
    data_vars={'tide_height':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids))+0.1,attrs={'units':'m'}),
               'tide_salinity':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids)),attrs={'units':'ppt'}),
               # Setting nonzero nitrate so leaching doesn't become a problem
               'tide_nitrate':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids))+0.3e-3,attrs={'units':'mol/L'}),
               },
    coords   ={'time':xarray.Variable(('time',),data=np.arange(ntimes),attrs={'units':'hours'}),
                'gridcell':np.arange(num_grids)},
    attrs    ={'Description':'Hydrological boundary conditions for trough/polygon grid cells'}
)



# Make new multi-grid cell configuration. Treating each land cover type as a separate grid cell
# This is easier for setting up tidal forcing, but won't work in a larger scale simulation where grid cells need to be spatially defined
# Long term solution is to use topo units, but will need to figure a way to do hydro forcing in that framework
# Here we start from the single grid cell configuration for the site and then multiply it into multiple grid cells
site='beo'
# domain_onecol=xarray.open_dataset(f'/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_{site}/domain.nc')
# surfdata_onecol=xarray.open_dataset(f'/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_{site}/surfdata.nc')
# landuse_onecol=xarray.open_dataset(f'/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_{site}/surfdata.pftdyn.nc')

# BUT there is another set of domain files here which includes a gridded setup for BEO. Need to ask Fengming about best starting point.
# 7 point polygon transect: /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_7x1pt_beo-NGEE_areaC_simyr1850_c20230320-arctic.nc
# arctic pfts one has a lot of PFTs, need to ask Fengming which parameter file is appropriate for that
# surfdata_multicell_arcticpfts=xarray.open_dataset('/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_7x1pt_beo-NGEE_areaC_simyr1850_c20230320-arctic.nc')
surfdata_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_7x1pt_beo-NGEE_areaC_simyr1850_c20230320.nc')
domain_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/share/domains/domain.clm/domain.lnd.7x1pt_beo-NGEE_areaC_navy.nc')

domain_multicell.to_netcdf('BEO_domain_wetlandtype_gridcells.nc')

# Add new surface data fields specific to gridded hydrological forcing
# Let's define land surface heights relative to the trough
# ht_above_stream in meters units
surfdata_multicell['ht_above_stream'] = surfdata_multicell['TOPO']-surfdata_multicell['TOPO'][0]
# Here we specify the height of the polygon relative to the "zero" point in the hydrological forcing (in meters)
# Hubbard et al. (2013) has a detailed analysis of surface heights at BEO using Lidar https://link.springer.com/article/10.1007/s10040-012-0939-y
# Distance proportional to height? Or constant?
surfdata_multicell['dist_from_stream'] = surfdata_multicell['ht_above_stream']*0.0 + 1.0

# We can change soil texture including organic content which is important for hydrological and thermal properties
# that contribute to active layer thickness
# BEO surface data already includes different organic matter percentages
# surfdata_multicell['ORGANIC'][...] = surfdata_multicell['ORGANIC'][...]*3.0
# surfdata_multicell['PCT_SAND'][...] = 75.0
# surfdata_multicell['PCT_CLAY'][...] = 15.0
# fdrain controls the relationship between water table depth and topographic drainage
# Higher number -> drainage declines faster as WT goes below surface
surfdata_multicell['fdrain']=surfdata_multicell['dist_from_stream']*0.0 + 100.0

surfdata_multicell['SLOPE'][:]=0.05

# Change PFTs to match different ecosystems. This is the list of default ELM PFTs
pftnames = [s.strip() for s in [
  "not_vegetated                           ",
  "needleleaf_evergreen_temperate_tree     ",
  "needleleaf_evergreen_boreal_tree        ",
  "needleleaf_deciduous_boreal_tree        ",
  "broadleaf_evergreen_tropical_tree       ",
  "broadleaf_evergreen_temperate_tree      ",
  "broadleaf_deciduous_tropical_tree       ",
  "broadleaf_deciduous_temperate_tree      ",
  "broadleaf_deciduous_boreal_tree         ",
  "broadleaf_evergreen_shrub               ",
  "broadleaf_deciduous_temperate_shrub     ",
  "broadleaf_deciduous_boreal_shrub        ",
  "c3_arctic_grass                         ",
  "c3_non-arctic_grass                     ",
  "c4_grass                                ",
  "c3_crop                                 ",
  "c3_irrigated                            ",
  "corn                                    ",
  "irrigated_corn                          ",
  "spring_temperate_cereal                 ",
  "irrigated_spring_temperate_cereal       ",
  "winter_temperate_cereal                 ",
  "irrigated_winter_temperate_cereal       ",
  "soybean                                 ",
  "irrigated_soybean                       " ]
]

pftnames_arctic = [s.strip() for s in [
  "not_vegetated                           ",
  "arctic_lichen                           ",
  "arctic_bryophyte                        ",
  "arctic_evergreen_shrub_dwarf            ",
  "arctic_evergreen_shrub_tall             ",
  "arctic_deciduous_shrub_dwarf            ",
  "arctic_deciduous_shrub_low              ",
  "arctic_deciduous_shrub_tall             ",
  "arctic_deciduous_shrub_alder            ",
  "arctic_forb                             ",
  "arctic_dry_graminoid                    ",
  "arctic_wet_graminoid                    " 
]]

# The current setup for BEO is 33% nonvegetated, 2% evergreen boreal tree, 41% arctic shrub, and 24% arctic grass
# We can distribute these differently across land cover types, for example assume that troughs are mostly unvegetated and shrubs are mostly in polygon centers
# surfdata_multicell['PCT_NAT_PFT'][:]=0.0
# surfdata_multicell['PCT_NAT_PFT'][pftnames.index('c3_non-arctic_grass'),:,landcover_types.index('Fresh marsh')]=100.0

# Create a new version of the surface data with the number of PFTs matching the Arctic ones
surfdata_multicell_arcticpfts = surfdata_multicell.copy(deep=True).isel(natpft=slice(0,len(pftnames_arctic)),lsmpft=slice(0,len(pftnames_arctic)))
surfdata_multicell_arcticpfts['PCT_NAT_PFT'][1:,:]=0.0
surfdata_multicell_arcticpfts['PCT_NAT_PFT'][pftnames_arctic.index('arctic_deciduous_shrub_low'),:]=surfdata_multicell['PCT_NAT_PFT'][pftnames.index('broadleaf_deciduous_boreal_shrub'),:]
surfdata_multicell_arcticpfts['PCT_NAT_PFT'][pftnames_arctic.index('arctic_evergreen_shrub_dwarf'),:]=surfdata_multicell['PCT_NAT_PFT'][pftnames.index('broadleaf_evergreen_shrub'),:]

flooded=surfdata_multicell_arcticpfts['ht_above_stream']<tide_data_multicell['tide_height'].mean()
surfdata_multicell_arcticpfts['PCT_NAT_PFT'][pftnames_arctic.index('arctic_dry_graminoid'),~flooded]=surfdata_multicell['PCT_NAT_PFT'][pftnames.index('c3_arctic_grass'),~flooded]
surfdata_multicell_arcticpfts['PCT_NAT_PFT'][pftnames_arctic.index('arctic_wet_graminoid'),flooded]=surfdata_multicell['PCT_NAT_PFT'][pftnames.index('c3_arctic_grass'),flooded]

# PFT percents are required to sum to 100 in each grid cell or the model will crash
if (surfdata_multicell['PCT_NAT_PFT'].sum(dim='natpft')!=100).any():
    raise ValueError('PFTs do not all add up to 100')
if (surfdata_multicell_arcticpfts['PCT_NAT_PFT'].sum(dim='natpft')!=100).any():
    raise ValueError('Arctic PFTs do not all add up to 100')

tide_data_multicell.to_netcdf('BEO_hydro_BC_multicell.nc')
surfdata_multicell.to_netcdf('BEO_surfdata_multicell.nc')
surfdata_multicell_arcticpfts.to_netcdf('BEO_surfdata_multicell_arcticpfts.nc')
domain_multicell.to_netcdf('BEO_domain_multicell.nc')

import matplotlib.pyplot as plt
f,a=plt.subplots(num='Water heights',clear=True,nrows=1)
a.fill_between(np.arange(len(landcover_types)),np.zeros(len(landcover_types))-.1,surfdata_multicell['ht_above_stream'],ls='-',color='brown',label='Soil surface',step='mid')
a.axhspan(-0.1,tide_data_multicell['tide_height'].mean(),color='b',alpha=0.5,label='Water level')
plt.xticks(ticks=np.arange(len(landcover_types)),labels=landcover_types)
a.set_ylabel('Height (m)')
a.set(xlim=(0,len(landcover_types)-1.5),ylim=(-0.1,0.23),title='Polygon landform levels')
# a.legend()

plt.show()