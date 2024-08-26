# Multi column short test alquimia simulation
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
                       --model_root /home/b0u/models/E3SM-NGEE --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
                       --domainfile $domain \
                       --surffile $surf --np 7 --walltime 3 \
                       --caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/b0u/  --mpilib openmpi3 --pio_version 2 \
                       --hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
                       --trans_varlist $varlist \
                       --alquimia /home/b0u/models/PFLOTRAN/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in

# Single column simulation
site=beo
domain=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site/domain.nc
surf=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site/surfdata.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,QFLX_ADV,\
QFLX_LAT_AQU,QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,QDRAI_VR,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate"
python site_fullrun.py --site AK-BEO --sitegroup NGEEArctic --caseidprefix Alaska_alquimia_test_onecol  \
                       --noad --tstep 1 --nyears_final_spinup 15 --nyears_transient 15 \
                       --cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
                       --model_root /home/b0u/models/E3SM-NGEE --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
                       --domainfile $domain \
                       --surffile $surf --np 1 --walltime 3 \
                       --caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/b0u/  --mpilib openmpi3 --pio_version 2 \
                       --hist_nhtfrq_trans -1 --hist_mfilt_trans 8760 --hist_mfilt_spinup 0 --hist_nhtfrq_spinup 12 --cn_only \
                       --trans_varlist $varlist \
                       --alquimia /home/b0u/models/PFLOTRAN/REDOX-PFLOTRAN/ELM_decks/CTC_alquimia_forELM_O2consuming.in


# Test IM4 simulation
# I don't think this code includes the capability for reading in Fengming's site meteorology but I can check with him
# git checkout ngee-arctic-IM4/merged
site=kougarok
metdir=/nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata/atm/datm7/GSWP3_daymet/cpl_bypass_$site
domain=/home/b0u/models/E3SM-NGEE/components/elm/cime_config/testdefs/testmods_dirs/elm/usrpft_arctic_I1850CNPRDCTCBC/domain.lnd.6x1pt_kougarok-NGEE_TransA_navy.nc
surf=/home/b0u/models/E3SM-NGEE/components/elm/cime_config/testdefs/testmods_dirs/elm/usrpft_arctic_I1850CNPRDCTCBC/surfdata_6x1pt_kougarok-NGEE_TransA_simyr1850_c20201008-sub12.nc
parmfile=/home/b0u/models/E3SM-NGEE/components/elm/cime_config/testdefs/testmods_dirs/elm/usrpft_arctic_I1850CNPRDCTCBC/clm_params_c180524-sub12_updated20240201.nc
varlist="TOTVEGC,TOTSOMC,TOTLITC,SOIL1C_vr,SOIL2C_vr,SOIL3C_vr,SOIL4C_vr,LITR1C_vr,LITR2C_vr,LITR3C_vr,LEAFC,\
soil_O2,HR,GPP,NEE,NPP,SMINN,SMINN_TO_PLANT,DIC_vr,SIC_vr,H2OSOI,H2OSFC,SOILLIQ,SOILICE,ZWT,\
QFLX_EVAP_TOT,QVEGT,watsat,chem_dt,soil_salinity,soil_pH,DOC_vr,DIC_vr,DOC_RUNOFF,DIC_RUNOFF,SMIN_NO3_RUNOFF,\
soil_sulfate,soil_sulfide,CH4_vr,CH4FLUX_ALQUIMIA,QDRAI,TSOI,soil_Fe2,soil_FeOxide,soil_FeS,soil_acetate"
varlist_pft="LEAFC,FROOTC,WOODC,LIVECROOTC,TLAI,CPOOL,STORVEGC,STORVEGN,LEAFN,FROOTN,LIVECROOTN,LIVESTEMC,LIVESTEMN,DEADCROOTC,DEADCROOTN,DEADSTEMC,DEADSTEMN,TOTVEGC,TOTVEGN,GPP,\
NPP,HTOP,AGNPP,BGNPP,LEAF_MR,FROOT_MR,LIVESTEM_MR,LIVECROOT_MR,MR,VCMAX25TOP,GR,AVAILC,PLANT_CALLOC,EXCESS_CFLUX,XSMRPOOL_RECOVER,XSMRPOOL,CPOOL,LEAFC_XFER_TO_LEAFC,FROOTC_XFER_TO_FROOTC,DOWNREG,\
INIT_GPP,SMINN_TO_NPOOL,PLANT_PDEMAND,PLANT_NDEMAND,SMINP_TO_PPOOL"
python site_fullrun.py --site AK-K64 --sitegroup NGEEArctic --caseidprefix Kougarok_IM4_Arcticpfts \
                       --nyears_ad_spinup 200 --tstep 1 --nyears_final_spinup 600 --walltime 18 \
                       --cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
                       --model_root /home/b0u/models/E3SM-NGEE --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
                       --domainfile $domain \
                       --surffile $surf --np=6 --maxpatch_pft=12 --var_soilthickness --mod_parm_file $parmfile \
                       --caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/b0u/  --mpilib openmpi3 --pio_version 2 \
                       --dailyvars --var_list_pft $varlist_pft


python site_fullrun.py --site AK-K64 --sitegroup NGEEArctic --caseidprefix Kougarok_IM4_E3SMpfts \
                       --nyears_ad_spinup 200 --tstep 1 --nyears_final_spinup 600 --walltime 18 \
                       --cpl_bypass --machine cades --no_dynroot --gswp3 --nofire --nopftdyn --nopointdata \
                       --model_root /home/b0u/models/E3SM-NGEE --ccsm_input /nfs/data/ccsi/proj-shared/E3SM/pt-e3sm-inputdata \
                       --exeroot /lustre/or-scratch/cades-ccsi/b0u/Kougarok_IM4_test_AK-K64_ICB1850CNRDCTCBC_ad_spinup/bld/ \
                       --domainfile $domain \
                       --surffile /home/b0u/Kougarok_param_edits/param_files/surfdata_Kougarok_defaultPFTs_all-shrubs-decid-boreal.nc \
                       --np=6 --var_soilthickness \
                       --caseroot ~/cases --runroot /lustre/or-scratch/cades-ccsi/b0u/  --mpilib openmpi3 --pio_version 2 \
                       --dailyvars --var_list_pft $varlist_pft

