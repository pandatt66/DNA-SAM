
SHELL := /bin/bash

USER = qb

#directories
ifeq ($(USER),xtong3)
PYTHON_DIR = /usr/local/anaconda/bin
else ifeq ($(USER),qb)
PYTHON_DIR =/usr/local/packages/python/2.7.7-anaconda/bin
endif

NUCLEOTIDES = dAMP dCMP dGMP dTMP

.PHONY: all \
dNMP_traj \
dNMP_wall_mindist \
equil_time\
ads_des_stats

all: dNMP_traj \
dNMP_wall_mindist \
equil_time\
ads_des_stats


#####################################################
#dNMP trajectory
#####################################################
dNMP_traj:$(foreach n,$(NUCLEOTIDES),$(eval $(call NOREST_TPR,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call NONWATER_TPR,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call NONWATER_INDEX,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call dNMP_COORD,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call dNMP_COORD_FIXED,$(n)))) \

define NOREST_TPR
x_norest_$(1).tpr:x_temp.mdp
	gmx grompp -f x_temp.mdp -c md_x1_$(1).gro -r ref_throw_$(1).pdb -p ch3_md_$(1).top -n index_ch3_$(1).ndx -o x_norest_$(1).tpr -maxwarn 1 

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call NOREST_TPR,$(n))))

define NONWATER_TPR

x_nonwater_$(1).tpr:x_norest_$(1).tpr
	echo 7|gmx convert-tpr -s x_norest_$(1).tpr -o x_nonwater_$(1).tpr

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call NONWATER_TPR,$(n))))


define NONWATER_INDEX

x_nonwater_$(1).ndx:x_nonwater_$(1).tpr
	source make_nonwater_ndx.bash $(1)
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call NONWATER_INDEX,$(n))))


define dNMP_COORD
$(1)_1_coord.xvg:x_nonwater_$(1).ndx md_x1_$(1).xtc
	echo 5| gmx traj -f md_x1_$(1).xtc -s x_nonwater_$(1).tpr -n x_nonwater_$(1).ndx -com -ox $(1)_1_coord.xvg
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call dNMP_COORD,$(n))))


define dNMP_COORD_FIXED
$(1)_1_coord_fixed.dat:$(1)_1_coord.xvg md_x1_$(1).gro read_xvg.py fix_dNMP_traj.py
	$(PYTHON_DIR)/ipython fix_dNMP_traj.py md_x1_$(1).gro $(1)_1_coord.xvg $(1)_1_coord_fixed

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call dNMP_COORD_FIXED,$(n))))


#################################################
#Plot rmsd of monolayers, detect equilibrium time 
#################################################
equil_time:$(foreach n,$(NUCLEOTIDES),$(eval $(call PBC_MOL,$(n))))\
$(foreach n,$(NUCLEOTIDES),$(eval $(call PBC_NOJUMP,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call FIRST_FRAME,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call RMSD,$(n)))) \
$(foreach n,$(NUCLEOTIDES),$(eval $(call EQUIL_TIME,$(n))))\

define PBC_MOL

md_x1_$(1)_mol.xtc: md_x1_$(1).xtc x_nonwater_$(1).tpr
	echo 0| gmx trjconv -f md_x1_$(1).xtc  -s x_nonwater_$(1).tpr -o md_x1_$(1)_mol.xtc -pbc mol

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PBC_MOL,$(n))))

define PBC_NOJUMP

md_x1_$(1)_nojump.xtc:md_x1_$(1)_mol.xtc
	gmx trjconv -f md_x1_$(1)_mol.xtc -o md_x1_$(1)_nojump.xtc -pbc nojump 

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PBC_NOJUMP,$(n))))

define FIRST_FRAME

md_x1_$(1)_1stfram.gro:md_x1_$(1)_nojump.xtc
	echo 0|gmx trjconv -f md_x1_$(1)_nojump.xtc -s x_nonwater_$(1).tpr -n x_nonwater_$(1).ndx -o md_x1_$(1)_1stfram.gro -dump 0

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call FIRST_FRAME,$(n))))

define RMSD

rmsd_$(1).xvg:md_x1_$(1)_1stfram.gro
	echo 2 2|gmx rms -f md_x1_$(1)_nojump.xtc -s md_x1_$(1)_1stfram.gro -n x_nonwater_$(1).ndx -o rmsd_$(1).xvg
	
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call RMSD,$(n))))


#########################################################################
#Have created python envioronment named my_root (source activate my_root)
#########################################################################
define EQUIL_TIME

equiltime_$(1).dat: rmsd_$(1).xvg equiltime.py
	python equiltime.py rmsd_$(1).xvg equiltime_$(1)

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call EQUIL_TIME,$(n))))

#################################################
#Adsorption, desorption statistics
#################################################
ads_des_stats:$(foreach n,$(NUCLEOTIDES),$(eval $(call TIME_MIN,$(n))))\
$(foreach n,$(NUCLEOTIDES),$(eval $(call NUM_CONT,$(n))))\
$(foreach n,$(NUCLEOTIDES),$(eval $(call ADS_DES_STATS,$(n))))\

#################################################
#dNMP-ch3 minimum distance
#################################################

#define TIME_MIN

#time_$(1):= $(shell cat equiltime_$(1).dat)

#endef

#$(foreach n,$(NUCLEOTIDES),$(eval $(call TIME_MIN,$(n))))

define NUM_CONT

$(1)_1_numcont.xvg:md_x1_$(1).xtc
#echo 2 5| gmx mindist -f  md_x1_$(1).xtc -s x_nonwater_$(1).tpr -n x_nonwater_$(1).ndx -od $(1)_1_mindist.xvg -on $(1)_1_numcont.xvg -b $(time_$(1)) -d 0.35#

	echo 2 5| gmx mindist -f  md_x1_$(1).xtc -s x_nonwater_$(1).tpr -n x_nonwater_$(1).ndx -od $(1)_1_mindist.xvg -on $(1)_1_numcont.xvg -b 50000 -d 0.35

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call NUM_CONT,$(n))))

define ADS_DES_STATS

$(1)_1_ads_times.dat: read_xvg.py ads_des_times.py ads_des_stats.py $(1)_1_coord_fixed.dat $(1)_1_numcont.xvg
	$(PYTHON_DIR)/ipython ads_des_stats.py $(1)_1_mindist.xvg $(1)_1_coord_fixed.dat $(1)_1 -ic $(1)_1_numcont.xvg -mt 0.3

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call ADS_DES_STATS,$(n))))

#############
#PMF
#####################
define NDX_BASE
index_ch3_$(1).ndx:make_ndx.bash
	source make_ndx.bash throw_$(1).gro $(1) ch3
	source monolayer_selections_to_index_file.bash $(1) CH3 SH 35 ch3 
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call NDX_BASE,$(n))))

define CTR_SH
$(1)_ctr_SH.trr:index_ch3_$(1).ndx
	echo 9 0|gmx trjconv -f md_x1_$(1).trr -s md_x1_$(1).tpr -o $(1)_ctr_SH.trr -pbc mol -n index_ch3_$(1).ndx -center

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call CTR_SH,$(n))))

define PMF_SH

pmf_$(1)_SH.xvg:$(1)_ctr_SH.trr
	echo 11|gmx traj -f $(1)_ctr_SH.trr -s md_x1_$(1).tpr -n index_ch3_$(1).ndx -com -ox -nox -noy pmf_$(1)_SH.xvg

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PMF_SH,$(n))))


define PMF_BASE
pmf_$(1)_base.xvg:$(1)_ctr_SH.trr pmf_$(1)_SH.xvg
	echo 10|gmx traj -f $(1)_ctr_SH.trr -s md_x1_$(1).tpr -n index_ch3_$(1).ndx -com -ox -nox -noy pmf_$(1)_base.xvg
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PMF_BASE,$(n))))

define FIND_TIME

pmf_$(1)_time.dat:pmf_$(1)_base.xvg
	ipython find_pmf_time.py pmf_$(1)_SH.xvg pmf_$(1)_base.xvg $(1) 15
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call FIND_TIME,$(n))))

##############################################################
#THe following:
#i)dump all the selectd time for configuration; 
#ii)translate system to put solvent in center of box 
#iii)generate tpr files
############################################################
define EXT_CONFIG
ref_pmf_$(1).pdb:trajconv_dump.bash pmf_$(1)_time.dat
	source trajconv_dump.bash $(1) pmf_$(1)_time.dat -3.523 ch3
endef
$(foreach n,$(NUCLEOTIDES),$(eval $(call EXT_CONFIG,$(n))))

##################################
#mdrun for all windowson cluster(....)
###################################
define PMF_PROF
$(1)_profile.xvg:$(1)_tpr-files.dat $(1)_pullf-files.dat
	gmx wham -it $(1)_tpr-files.dat -if $(1)_pullf-files.dat -o $(1)_profile.xvg -hist $(1)_histo.xvg -unit kT
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PMF_PROF,$(n))))


#ll dTMP*.xvg , get rid of details:cat 'dTMP_pullf-files.dat'| awk {'print $9'}
#find reference position in pmf: 
#gmx traj -f dCMP_x1_2.29056.trr -s dCMP_x1_2.29056.tpr -n index_ch3_dCMP.ndx -com -ox -nox -noy dCMP_SH_2.29056.xvg
#gmx traj -f dCMP_x1_2.29056.trr -s dCMP_x1_2.29056.tpr -n index_ch3_dCMP.ndx -com -ox -nox -noy dCMP_SH_2.29056.xvg
#python shift_pmf.py dCMP_profile.xvg dCMP dCMP_SH_2.29056.xvg dCMP_edgrp_2.29056.xvg

define PMF_SHFT
pmf_shift_$(1).dat: $(1)_profile.xvg
	python shift_pmf.py $(1)_profile.xvg $(1)
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PMF_SHFT,$(n))))


define CTR_SOL
$(1)_ctr_SOL.trr:
	echo 0|trjconv -f $(1)_ctr_SH.trr -s md_x1_$(1).tpr -pbc mol -n index_ch3_$(1).ndx -trans 0 0 -3.552 -o $(1)_ctr_SOL.trr

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call CTR_SOL,$(n))))

define SH1
pmf_$(1)_SH1.xvg:$(1)_ctr_SOL.trr
	echo 11|gmx traj -f $(1)_ctr_SOL.trr -s md_x1_$(1).tpr -n index_ch3_$(1).ndx -com -ox -nox -noy pmf_$(1)_SH1.xvg
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call SH1,$(n))))

define SH2
pmf_$(1)_SH2.xvg:$(1)_ctr_SOL.trr pmf_$(1)_SH1.xvg
	echo 12|gmx traj -f $(1)_ctr_SOL.trr -s md_x1_$(1).tpr -n index_ch3_$(1).ndx -com -ox -nox -noy pmf_$(1)_SH2.xvg
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call SH2,$(n))))

define PMF_CTR_BASE
pmf_$(1)_base_ctr.xvg:$(1)_ctr_SOL.trr pmf_$(1)_SH2.xvg
	echo 10|gmx traj -f $(1)_ctr_SOL.trr -s md_x1_$(1).tpr -n index_ch3_$(1).ndx -com -ox -nox -noy pmf_$(1)_base_ctr.xvg
endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PMF_CTR_BASE,$(n))))


define PRBLITY
$(1)_pmf_prob.dat:pmf_$(1)_base_ctr.xvg pmf_$(1)_SH1.xvg pmf_$(1)_SH2.xvg 
	ipython pmf_prob.py pmf_shift_$(1).dat pmf_$(1)_base_ctr.xvg pmf_$(1)_SH1.xvg pmf_$(1)_SH2.xvg $(1)

endef

$(foreach n,$(NUCLEOTIDES),$(eval $(call PRBLITY,$(n))))











