import xarray
import numpy as np

#BAM editing to reflect changes after discussion w/Neslihan and evaluation of exsiting BEO data
# Set up an ELM domain with polygons, based on Fengming's dataset from Wang et al 2019 JAMES paper
# 'This is for 8 land-units summarized from NGEE Site Area C with additional Arctic PFTS: 1-LCPtrough, 2-LCPcenter, 3-LCPrim, 4-FCPtrough, 5-FCPcenter, 6-FCPrim, 7-HCPtrough, 8-HCPcenter, 9-average
#BAM: jk for first testing just keeping the 7-cell structure and changing individual variables, will figure out how to add the other 2 gridcells later
#once I make sure that this works (so dropping 3-LCPrim, 9-average for now)
#BAM 9/17/24 edit: replacing FCPrim w/LCPrim, no difference in hydrological regime between FCPcenter and FCPrim so makes more sense to include LCPrim if 
#sticking w/the 7-cell structure for now

landcover_types=[
  'LCPtrough',
  'LCPcenter',
  'FCPtrough',
    'FCPcenter',
    'LCPrim',
  'HCPtrough',
  'HCPcenter',
]

# For now we are setting up each land cover type as a separate grid cell. 
# In the future we could do some of this using topo units within grid cells

# Hydrology uses the coastal wetland configuration setup for specifying a time series of hydrological boundary condition
# Hourly time series, just do one year and the model will repeat it. Currently using water level of zero. 
# Surface data set below will define ground surface height above drainage to give different hydrological conditions
ntimes=365*24
# 7 grid cells 
num_grids=len(landcover_types)
#BAM: sets up array w/dimensions corresponding to number of grid cells to be used. Think can drop for this go around too, since using the 
#files already created by Ben as a starting point, just edit the tide height variable the same way editing organic matter and veg fractional cover
#only editing tide height, leaving salinity and nitrate as they are
#tide_data_multicell=xarray.Dataset(
#    data_vars={'tide_height':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids))+0.1,attrs={'units':'m'}),
#               'tide_salinity':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids)),attrs={'units':'ppt'}),
               # Setting nonzero nitrate so leaching doesn't become a problem
#               'tide_nitrate':xarray.Variable(('time','gridcell'),data=np.zeros((ntimes,num_grids))+0.3e-3,attrs={'units':'mol/L'}),
#               },
#    coords   ={'time':xarray.Variable(('time',),data=np.arange(ntimes),attrs={'units':'hours'}),
#                'gridcell':np.arange(num_grids)},
#    attrs    ={'Description':'Hydrological boundary conditions for polygon grid cells'}
#)





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
##surfdata_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/lnd/clm2/surfdata_map/surfdata_7x1pt_beo-NGEE_areaC_simyr1850_c20230320.nc')
##domain_multicell=xarray.open_dataset('/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/share/domains/domain.clm/domain.lnd.7x1pt_beo-NGEE_areaC_navy.nc')

#BAM:just pull in the domain and surface files that Ben already set up for BEO w/this code, don't need to start from Fengming's original stuff again
#also changing to do locally for testing, downloaded the files from CADES
surfdata_multicell=xarray.open_dataset('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_surfdata_multicell.nc')
domain_multicell=xarray.open_dataset('~/GitHub/NGEE_ELM_BGC_Bailey/BEO_domain_multicell.nc')

#BAM: not changing anything in the domain file yet so can just use the one thats already been created, don't have to run this again
##domain_multicell.to_netcdf('BEO_domain_wetlandtype_gridcells.nc')

# Add new surface data fields specific to gridded hydrological forcing
# Let's define land surface heights relative to the trough
# ht_above_stream in meters units
# Here we specify the height of the polygon relative to the "zero" point in the hydrological forcing (in meters)
# Hubbard et al. (2013) has a detailed analysis of surface heights at BEO using Lidar https://link.springer.com/article/10.1007/s10040-012-0939-y
#BAM: this is what determines the soil height to water height ratio, aka the inundation level for the different polygon types
#surfdata_multicell['ht_above_stream'] = surfdata_multicell['TOPO']-surfdata_multicell['TOPO'][0]
#BAM: DO I ALSO NEED TO CHANGE THINGS IN THE HYDRO FILE...? OR CAN IT JUST BE SET UP HERE? CHECK!!!!!
new_inundation_values = [[0], #gridcell 1 (LCPtrough)
                        [0.02991], #gridcell 2 (LCPcenter)
                        [0.008905], #gridcell 3 (FCPtrough)
                        [0.1228581], #gridcell 4 (FCPcenter)
                        [0.11896484], #gridcell 5 (LCPrim)
                        [0.01781], #gridcell 6 (HCPtrough)
                        [0.21580621]] #gridcell 7 (HCPcenter)

# Convert list of new inundation values to a numpy array
new_inundation_values = np.array(new_inundation_values)

# Update the ht_above_stream variable with new values
surfdata_multicell['ht_above_stream'][:] = new_inundation_values


# Distance proportional to height? Or constant?
surfdata_multicell['dist_from_stream'] = surfdata_multicell['ht_above_stream']*0.0 + 1.0
#BAM: dist_from_stream is set as 1 for everything

# We can change soil texture including organic content which is important for hydrological and thermal properties
# that contribute to active layer thickness
# BEO surface data already includes different organic matter percentages
# surfdata_multicell['ORGANIC'][...] = surfdata_multicell['ORGANIC'][...]*3.0
# surfdata_multicell['PCT_SAND'][...] = 75.0
# surfdata_multicell['PCT_CLAY'][...] = 15.0

new_organic_values = [
    [102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24],  # Values for gridcell 1 (LCPtrough)
    [102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24],  # Values for gridcell 2 (LCPcenter)
    [101.452, 101.452, 101.452, 29.237, 28.886, 50.895, 13.61285714, 13.624, 5.27475, 5.148],  # Values for gridcell 3 (FCPtrough)
    [101.452, 101.452, 101.452, 29.237, 28.886, 50.895, 13.61285714, 13.624, 5.27475, 5.148],  # Values for gridcell 4 (FCPcenter)
    [102.973, 56.927, 68.51975, 43.667, 42.3306, 45.864, 22.02525, 6.24, 6.24, 6.24],  # Values for gridcell 5 (LCPrim)
    [49.621, 49.621, 49.621, 80.314, 62.16166667, 73.53666667, 52.71825, 10.647, 6.097, 6.097],  # Values for gridcell 6 (HCPtrough)
    [120.081, 120.081, 11.895, 59.0525, 62.16166667, 73.53666667, 52.71825, 10.647, 6.56175, 4.706],  # Values for gridcell 7 (HCPcenter)
]

# Convert list of new organic matter density values to a numpy array
new_organic_values = np.array(new_organic_values)

# Update the ORGANIC variable with new values
surfdata_multicell['ORGANIC'][:, :] = new_organic_values

# fdrain controls the relationship between water table depth and topographic drainage
# Higher number -> drainage declines faster as WT goes below surface
#BAM: set as = 100 for all polygon types
surfdata_multicell['fdrain']=surfdata_multicell['dist_from_stream']*0.0 + 100.0

surfdata_multicell['SLOPE'][:]=0.05

#BAM: change soil order from 2 (andisols) to 3 (gelisols) 
#see https://passel2.unl.edu/view/lesson/2eafec8dd762/3#:~:text=This%20lesson%20will%20examine%20each,Histosols%2C%20Aridisols%2C%20and%20Vertisols. 
surfdata_multicell['ORDER'][:]=3

#BAM: dropping all the PFT changing code for now since I'm pulling in the surface files Ben already did this for
#just changing the veg fractional coverage values the same way I'm changing the organic matter density values (above)
#...also ignoring the arctic PFT versions of things for now

# The current setup for BEO is 33% nonvegetated, 2% evergreen boreal tree, 41% arctic shrub, and 24% arctic grass
# We can distribute these differently across land cover types, for example assume that troughs are mostly unvegetated and shrubs are mostly in polygon centers
# surfdata_multicell['PCT_NAT_PFT'][:]=0.0
# surfdata_multicell['PCT_NAT_PFT'][pftnames.index('c3_non-arctic_grass'),:,landcover_types.index('Fresh marsh')]=100.0

#BAM: these are veg fractional cover for default PFTs (so 17 PFT options)
#we still are working w/the version w/default arctic PFTs only
new_veg_cover_values = [
    [5.241799121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 94.75820088, 0, 0, 0, 0],  # Values for gridcell 1 (LCPtrough)
    [43.04960196, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 56.95039804, 0, 0, 0, 0],  # Values for gridcell 2 (LCPcenter)
    [2.781795927, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 97.21820407, 0, 0, 0, 0],  # Values for gridcell 3 (FCPtrough)
    [27.09212377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 72.90787623, 0, 0, 0, 0],  # Values for gridcell 4 (FCPcenter)
    [11.87843147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.040826416, 87.08074211, 0, 0, 0, 0],  # Values for gridcell 5 (LCPrim)
    [15.46722699, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 84.53277301, 0, 0, 0, 0],  # Values for gridcell 6 (HCPtrough)
    [45.68772594, 0, 0, 0, 0, 0, 0, 0, 0, 11.4764446, 0, 0, 42.83582946, 0, 0, 0, 0],  # Values for gridcell 7 (HCPcenter)
]

# Convert list of new veg fractional cover values to a numpy array
new_veg_cover_values = np.array(new_veg_cover_values)

# Update PCT_NAT_PFT with new values
surfdata_multicell['PCT_NAT_PFT'][:, :] = new_veg_cover_values


# PFT percents are required to sum to 100 in each grid cell or the model will crash
if (surfdata_multicell['PCT_NAT_PFT'].sum(dim='natpft')!=100).any():
    raise ValueError('PFTs do not all add up to 100')

#BAM: leaving PCT_SAND and PCT_CLAY alone right now, waiting on data from Neslihan

#BAM: save everything as new versions so I don't get confuuuuused
#tide_data_multicell.to_netcdf('BEO_hydro_BC_multicell.nc') #didn't actually change anything here, just use the original file Ben created
surfdata_multicell.to_netcdf('BEO_surfdata_multicell_BAM.nc')
#surfdata_multicell_arcticpfts.to_netcdf('BEO_surfdata_multicell_arcticpfts.nc')
domain_multicell.to_netcdf('BEO_domain_multicell_BAM.nc')

import matplotlib.pyplot as plt
f,a=plt.subplots(num='Water heights',clear=True,nrows=1)
a.fill_between(np.arange(len(landcover_types)),np.zeros(len(landcover_types))-.1,surfdata_multicell['ht_above_stream'],ls='-',color='brown',label='Soil surface',step='mid')
a.axhspan(-0.1,tide_data_multicell['tide_height'].mean(),color='b',alpha=0.5,label='Water level')
plt.xticks(ticks=np.arange(len(landcover_types)),labels=landcover_types)
a.set_ylabel('Height (m)')
a.set(xlim=(0,len(landcover_types)-1.5),ylim=(-0.1,0.23),title='Polygon landform levels')
# a.legend()

plt.show()