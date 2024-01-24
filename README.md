# unh-edna-pipelines for UNH MBON
Code for eDNA metabarcoding on RON and Premise @UNH

To do:  
- Upload mifish classifier and reference taoxnomy/sequences
- Share the feature tableand ASVs of water samples to fecal group


### Qiime workflow
[MBON Workflow for MiFish](sh/MBON_Workflow.sh)

### Refernce database info  
[Mitofish](https://mitofish.aori.u-tokyo.ac.jp/download/) for all Mifish primers  

### feature-classifier info  
We use qiime's 'classify-hybrid-vsearch-sklearn', which may be more robust to parameter changes  

### Comparison to fecal
- Level of bacterial contam? (what percent of )

### Cross reference of GBIF for each 

### Example outputs:    
BIOM Format  
![BIOM Format](images/qiime-biom.png)  


Custom ASV output format  
![custom format](images/custom.png)  