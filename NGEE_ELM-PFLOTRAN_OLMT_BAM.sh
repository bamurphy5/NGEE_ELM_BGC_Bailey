# Multi column short test alquimia simulation, BAM version w/updated organic matter density, inundation levels, and veg fractional coverage
site=beo
metdir=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=/home/9v6/NGEE_ELM/BEO_domain_multicell.nc
surf=/home/9v6/NGEE_ELM/BEO_surfdata_multicell.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_test  \
--noad --tstep 1 --nyears_final_spinup 15 --nyears_transient 15 \
--cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root /home/9v6/ELM-alquimia/E3SM --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 3 \
--caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/9v6/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia /home/9v6/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in 


# Multi column short test alquimia simulation, version WITHOUT updated organic matter density, inundation levels, and veg fractional coverage (so using Ben's original surface data files)
site=beo
metdir=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_domain_multicell.nc
surf=/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_surfdata_multicell.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_test  \
--noad --tstep 1 --nyears_final_spinup 15 --nyears_transient 15 \
--cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root /home/9v6/ELM-alquimia/E3SM --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 3 \
--caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/9v6/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia /home/9v6/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in 

# trying w/specifying param file, still Ben's original surf files
site=beo
metdir=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_domain_multicell.nc
surf=/nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/BEO_surfdata_multicell.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_test  \
--noad --tstep 1 --nyears_final_spinup 15 --nyears_transient 15 \
--cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
--model_root /home/9v6/ELM-alquimia/E3SM --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
--domainfile $domain \
--surffile $surf --np 7 --walltime 3 \
--caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/9v6/  --mpilib openmpi --pio_version 2 \
--hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
--trans_varlist $varlist \
--alquimia /home/9v6/ELM-alquimia/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in \
--parm_file /nfs/data/ccsi/proj-shared/b0u/NGEE_ELM/parms_BEO.txt
