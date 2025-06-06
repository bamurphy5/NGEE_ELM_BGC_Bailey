#4/28/25
#version of the job submission script for use w/CADES baseline setup (bails not 9v6)
#Ben recently updated the master tidal branch to incorporate the surface gas flux changes and to use freshwater instead of saltwater freezing point to resolve issues w/Teri lowering the tidal forcing freezing point that potentially were resulting in the lack of freezing at depth we were seeing
#this uses the short spin up for testing (200 years instead of 600 years)
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_1 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 151 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc \
--parm_file $HOME/NGEE_ELM/parms_BEO.txt

#5/22/25
#editing to use parms_BEO_may.txt which has higher rsub_global_max 
#also changed BEO_hydro_BC_multicell.nc so there's not nitrate in the water being forced...set all = 0
#version of the job submission script for use w/CADES baseline setup (bails not 9v6)
#Ben recently updated the master tidal branch to incorporate the surface gas flux changes and to use freshwater instead of saltwater freezing point to resolve issues w/Teri lowering the tidal forcing freezing point that potentially were resulting in the lack of freezing at depth we were seeing
#this uses the short spin up for testing (200 years instead of 600 years)
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_2 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 151 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc \
--parm_file $HOME/NGEE_ELM/parms_BEO_may.txt


#5/23/25
#editing to use parms_BEO_may.txt which has higher rsub_global_max (this time inc slightly from last time)
#version of the job submission script for use w/CADES baseline setup (bails not 9v6)
#Ben recently updated the master tidal branch to incorporate the surface gas flux changes and to use freshwater instead of saltwater freezing point to resolve issues w/Teri lowering the tidal forcing freezing point that potentially were resulting in the lack of freezing at depth we were seeing
#this uses the short spin up for testing (200 years instead of 600 years)
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_5 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 151 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc \
--parm_file $HOME/NGEE_ELM/parms_BEO_may.txt


site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro_higherfdrain.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_6 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 151 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc \
--parm_file $HOME/NGEE_ELM/parms_BEO_may.txt

#6/2/25 trying w/fdrain amped back up to 100 for all polygons and w/o passing the parm file, also shortening transient to make things quicker
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro_higherfdrain.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_8 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 150 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--alquimia_ad $HOME/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming_adspinup.in \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc
#--parm_file $HOME/NGEE_ELM/parms_BEO_may.txt


