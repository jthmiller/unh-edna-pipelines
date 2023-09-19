#!/bin/bash

## Set dada2 denoise-paired parameters for reference database, references taxa, denoise, truncation and trim parameters
if [ -z "$primer" ] ; then

     echo " primer should be set to 16s, 18s, mifish, or klymus" && exit

elif [ $primer == "rbcl" ] ; then
    
    fw1="AGGTGAAGTAAAAGGTTCWTACTTAAA"
    fw2="AGGTGAAGTTAAAGGTTCWTAYTTAAA"
    fw3="AGGTGAAACTAAAGGTTCWTACTTAAA"
 
    rv1="CCTTCTAATTTACCWACWACTG"
    rv2="CCTTCTAATTTACCWACAACAG"

    cutadapt_config="--p-front-f $fw1 --p-front-f $fw2 --p-front-f $fw3 --p-front-r $rv1 --p-front-r $rv2"


    polyg_len=150

    ## denoise
    ## trunc
    trunclenr=200
    trunclenf=200
    ## trim
    trimleftf=0
    trimleftr=0

elif [ $primer == "mifish" ] ; then

    ## denoise
    maxaccepts=1
    query_cov=0.4 
    perc_identity=0.85 

    ## trunc
    trunclenr=150
    trunclenf=150

    ## trim
    trimleftf=0
    trimleftr=0

    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/MiFish/MiFish-reference-taxonomy-2021.9.fasta.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/MiFish/MiFish-reference-taxonomy-2021.9.tax.qza}

elif [ $primer == "mitohelper" ] ; then


    fw=GTCGGTAAAACTCGTGCCAGC	
    rv=CATAGTGGGGTATCTAATCCCAGTTTG
   cutadapt_config="--p-front-f $fw --p-front-r $rv"
  
    ## denoise
    polyg_len=100

    ## trunc
    trunclenf=100
    trunclenr=90

    ## trim
    trimleftf=0
    trimleftr=0

    ## taxonomy
    maxaccepts=20
    query_cov=0.8 
    perc_identity=0.95
    weak_id=0.80 
    tophit_perc_identity=0.90

    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/MiFish/12S-seqs-derep-uniq.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/MiFish/12S-tax-derep-uniq.qza}
    blastdb=${blastdb:-/home/unhAW/jtmiller/watts/ref-database/MiFish/mitohelper/QIIME-compatible/blast/12S-seqs-derep-db/12S-seqs-derep-db}

elif [ $primer == "18s" ] ; then
	
    #fw=GTACACACCGCCCGTC	
    #rv=TGATCCTTCTGCAGGTTCACCTAC
    
    rv=GTACACACCGCCCGTC
    fw=TGATCCTTCTGCAGGTTCACCTAC

    cutadapt_config="--p-front-f $fw --p-front-r $rv"

    echo $cutadapt_config

    ## denoise
    polyg_len=100

    ## taxonomy
    maxaccepts=10
    query_cov=0.8 
    perc_identity=0.90 
    weak_id=0.80
    
    ## trunc
    trunclenf=100 
    trunclenr=100

    ## trim
    trimleftf=0
    trimleftr=0

    ## reference
    refreads=${refreads:-/home/share/databases/SILVA_databases/silva-138-99-seqs.qza}
    reftax=${reftax:-/home/share/databases/SILVA_databases/silva-138-99-tax.qza}
    blastdb=${blastdb:-/home/unhAW/jtmiller/watts/ref-database/SILVA/silva-138-99/blast/18S-db}

elif [ $primer == "16s_V4-V5" ] ; then

    fw=GTGYCAGCMGCCGCGGTAA	
    rv=CCGYCAATTYMTTTRAGTTT
    
    reads=single

    ## denoise
    polyg_len=150
    
    ##taxonomy
    #maxaccepts=5
    #query_cov=0.4 
    #perc_identity=0.7 

    ## trunc
    trunclenf=225
    trunclenr=225

    ## trim
    trimleftf=0
    #trimleftr=20

    ## reference
    refreads=${refreads:-/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.qza}
    reftax=${reftax:-/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/majority_taxonomy_all_levels.qza}

elif [ $primer == "klymus" ] ; then
    
    ## denoise
    maxaccepts=1
    query_cov=0.4 
    perc_identity=0.85 

    ## Devin used: trunclenf=175 and trunclenr=175
    ## I used 150 forward and reverse
    trunclenf=150
    trunclenr=150
    
    ## trim
    trimleftf=19
    trimleftr=17

    ## reference
    refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/klymus/length_filtered_klymus_full_trimmed_final.fasta.qza}
    reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/klymus/length_filtered_klymus_full_trimmed_final.tax.qza}

elif [ $primer == "leray_co1" ] ; then

	## primer sequence
	fw=GGWACWGGWTGAACWGTWTAYCCYCC
	rv=TANACYTCNGGRTGNCCRAARAAYCA

	cutadapt_config="--p-front-f $fw --p-front-r $rv"
	##fastp param
	polyg_len=175

	## denoise
	## trunc
	trunclenr=185
	trunclenf=185
	## trim
	trimleftf=21
	trimleftr=27
	## taxonomy
	maxaccepts=10
	query_cov=0.80 
	perc_identity=0.90
	weak_id=0.80 
	tophit_perc_identity=0.90
	
	refreads=${refreads:-/home/unhAW/jtmiller/watts/ref-database/CO1/bold_derep/bold_derep1_seqs.qza}
	reftax=${reftax:-/home/unhAW/jtmiller/watts/ref-database/CO1/bold_derep/bold_derep1_taxa.qza}


else

    echo " primer should be 16s, 18s, mifish, or klymus - all lower case" && exit 

fi

## Adapted from Devins workflows
## /home/genome/devin/qiime2_tools/16S_workflow.py: trunc_len_f=230, trunc_len_r=210,
## /home/genome/devin/qiime2_tools/18S_workflow.py: trunc_len_f=175, trunc_len_r=175,
## /home/genome/devin/qiime2_tools/MiFish_workflow.py: trunc_len_f=150, trunc_len_r=150,
## /home/genome/devin/qiime2_tools/16S_workflow.py: trim_left_f=19, trim_left_r=20,
## /home/genome/devin/qiime2_tools/18S_workflow.py: trim_left_f=22, trim_left_r=31,
## /home/genome/devin/qiime2_tools/MiFish_workflow.py: trim_left_f=21, trim_left_r=27,


## /home/genome/devin/qiime2_tools/16S_workflow.py: reference_reads="/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.qza",
## /home/genome/devin/qiime2_tools/16S_workflow.py: reference_taxonomy="/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/majority_taxonomy_all_levels.qza",
## /home/genome/devin/qiime2_tools/18S_workflow.py: reference_reads="/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/rep_set/rep_set_all/99/silva132_99.qza",
## /home/genome/devin/qiime2_tools/18S_workflow.py: reference_taxonomy="/home/share/databases/SILVA_databases/SILVA_132_QIIME_release/taxonomy/taxonomy_all/99/majority_taxonomy_all_levels.qza",


