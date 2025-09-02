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
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_hydro_higherfdrain2.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate,H2OSFC_TIDE,ALT,\
FCH4,RAIN,TSA,FSAT,ZWT_PERCH,TBOT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_arctic_BAM_9 \
--nyears_ad_spinup 200 --nyears_final_spinup 200 --tstep 1 --nyears_transient 100 \
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

#6/27/25: trying running ELM w/o PFLOTRAN turned on to just look at the default CH4 model output
#7/14/25 added outputs to have the vertically resolved CH4 and O2 for the default model as well, using surface file w/differences by poylgon type set up
#see notes in "testing_drainage_rates_2.R"
#edited to drop SIF as an output variable to save during transient simulations, not a valid output variable
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_polygon_diff.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,TSOI,H2OSFC_TIDE,ALT,\
FCH4,FCH4TOCO2,CH4PROD,RAIN,TSA,FSAT,ZWT_PERCH,TBOT,\
FINUNDATED,CH4_SURF_DIFF_SAT,CH4_SURF_DIFF_UNSAT,CH4_EBUL_TOTAL_SAT,CH4_EBUL_TOTAL_UNSAT,CH4_SURF_EBUL_SAT,\
CH4_SURF_EBUL_UNSAT,CH4_SURF_AERE_SAT,CH4_SURF_AERE_UNSAT,CONC_CH4_SAT,CONC_CH4_UNSAT,CH4_OXID_DEPTH_SAT,CH4_OXID_DEPTH_UNSAT,CONC_O2_SAT,CONC_O2_UNSAT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_defaultCH4_arctic_BAM_2 \
--nyears_ad_spinup 300 --nyears_final_spinup 400 --tstep 1 --nyears_transient 165 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc
#--parm_file $HOME/NGEE_ELM/parms_BEO_may.txt

#creating a version of the above to restart from transient since spin up worked it was just the transient period that failed b/c SIF wasn't a valid output variable to save
#tags --noad --nofnsp say to do no ad spinup and no final spin up
#--finidat tag is the absolute path of the restart file (*.elm.r.*.nc) to use to initialize the simulation
#don't need to specify --run_startyear b/c I'm starting the whole transient period over, would need to specify if only starting part of the transient period (like from 1980 on or something), it'll default to start in 1850 if not specified
#some examples used --finitfile tag to point to restart file and some used --finidat...not sure which is right, maybe --finitfile is when completely restarting transient period and --finidat is for restarting partway through?
#I think it also needs to be able to find the .cpl and .rho files to restart, not sure what tag to point to these, Fengmings example for the workshop just copies these files over from the directory where we uploaded them
#but maybe thats just needed when restarting partway through transient (we were trying to restart from year 2000), in his example he deleted the --finitfile command and copied those files over instead then used a CONTINUE_RUN=TRUE command to restart (https://github.com/dmricciuto/OLMT/commit/616629b8e7bafcfd13141309ce79a76c80fc34d5)
#hmm got error site_fullrun.py: error: no such option: --finitfile, maybe thats a docker specific command? Searching site_fullrun.py the description for --finitfile is 'initial ELM data file to start/restart' but its grouped under options for surface data, whereas --finidat is grouped under CASE options and described as 'Full path of ELM restart file to use (for transient only)'
#NEED TO MANUALLY EDIT --finidat FOR OTHER RUN CASES!!!!!!!!!!!!!!!!!!!!!!!!!!!!
site=beo
metdir=/gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_polygon_diff.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,TSOI,H2OSFC_TIDE,ALT,\
FCH4,FCH4TOCO2,CH4PROD,RAIN,TSA,FSAT,ZWT_PERCH,TBOT,\
FINUNDATED,CH4_SURF_DIFF_SAT,CH4_SURF_DIFF_UNSAT,CH4_EBUL_TOTAL_SAT,CH4_EBUL_TOTAL_UNSAT,CH4_SURF_EBUL_SAT,\
CH4_SURF_EBUL_UNSAT,CH4_SURF_AERE_SAT,CH4_SURF_AERE_UNSAT,CONC_CH4_SAT,CONC_CH4_UNSAT,CONC_O2_SAT,CONC_O2_UNSAT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_defaultCH4_arctic_BAM_2 \
--noad --nofnsp --tstep 1 --nyears_transient 165 \
--cpl_bypass --machine cades-baseline --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--trans_varlist $varlist \
--finidat /gpfs/wolf2/cades/cli185/scratch/bails/Alaska_defaultCH4_arctic_BAM_2_AK-BEO_ICB1850CNRDCTCBC/run/Alaska_defaultCH4_arctic_BAM_2_AK-BEO_ICB1850CNRDCTCBC.elm.r.0401-01-01-00000.nc \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc
#--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
#--finitfile $ELM_USER_DATA/OLMT_${site_code}_ICB1850CNPRDCTCBC.elm.r.0601-01-01-00000.nc \
#--finidat /inputdata/lnd/clm2/initdata/20230315_ARW_ICB20TRCNPRDCTCBC.elm.r.1980-01-01-00000.nc \

#9/2/25
#starting another round of the default sims but adding in a few output variables and switching the met forcing, so can't restart from transient need to go from the beginning, reducing final spin up to 300 years to make it a little quicker
#also edited soil water content during spin up per Katrina’s suggestion 
#replaced the --gswp3 tag with --daymet, didn't see a tag for ERA5 but maybe this will work since I'm using daymet corrected ERA5 and it says the --daymet tag is for daymet corrected met forcing? https://github.com/bamurphy5/OLMT/blob/master/site_fullrun.py 
site=BEO
metdir=/gpfs/wolf2/cades/cli185/world-shared/e3sm/inputdata/atm/datm7/Daymet_ERA5_ngee4/cpl_bypass_$site
domain=$HOME/NGEE_ELM/BEO_domain_multicell.nc
surf=$HOME/NGEE_ELM/BEO_surfdata_multicell_arcticpfts_polygon_diff.nc
paramfile=$HOME/NGEE_ELM/clm_params_arctic_updated.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
HR,ER,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,TSOI,H2OSFC_TIDE,ALT,\
FCH4,FCH4TOCO2,CH4PROD,RAIN,TSA,FSAT,ZWT_PERCH,TBOT,FSDS,EFLX_LH_TOT,FSH,\
FINUNDATED,CH4_SURF_DIFF_SAT,CH4_SURF_DIFF_UNSAT,CH4_EBUL_TOTAL_SAT,CH4_EBUL_TOTAL_UNSAT,CH4_SURF_EBUL_SAT,\
CH4_SURF_EBUL_UNSAT,CH4_SURF_AERE_SAT,CH4_SURF_AERE_UNSAT,CONC_CH4_SAT,CONC_CH4_UNSAT,CH4_OXID_DEPTH_SAT,CH4_OXID_DEPTH_UNSAT,CONC_O2_SAT,CONC_O2_UNSAT"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_defaultCH4_arctic_BAM_3 \
--nyears_ad_spinup 300 --nyears_final_spinup 300 --tstep 1 --nyears_transient 173 \
--cpl_bypass --machine cades-baseline --no_dynroot --nofire --nopftdyn --nopointdata \
--model_root $HOME/ELM-alquimia/E3SM --ccsm_input /gpfs/wolf2/cades/cli185/proj-shared/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 24 --maxpatch_pft 12 \
--mod_parm_file $paramfile \
--caseroot ~/cases --runroot /gpfs/wolf2/cades/cli185/scratch/bails/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--marsh --tide_forcing_file $HOME/NGEE_ELM/BEO_hydro_BC_multicell.nc




















