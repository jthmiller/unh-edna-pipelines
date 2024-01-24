# unh-edna-pipelines
Code for eDNA metabarcoding on RON and Premise @UNH

[MBON Workflow for MiFish](sh/MBON_Workflow.sh)

### Refernce database info  
[Mitofish](https://mitofish.aori.u-tokyo.ac.jp/download/) for all Mifish primers  

### feature-classifier info  
We use qiime's 'classify-hybrid-vsearch-sklearn', which may be more robust to parameter changes  

### Example outputs:    
BIOM Format  
![BIOM Format](images/qiime-biom.png)  


Custom ASV output format  
![custom format](images/custom.png)  