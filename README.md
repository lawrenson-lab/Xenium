Code used to explore data from the Multi-ethnic analysis of endometrioma project

Working pipeline for pembro subset
1. script/banksy_domains.R only for samples with targeted panel. Output: Banksy matrix
2. script/banksy_domains_pembro.R only for samples with targeted panel. 
   Output: Domains across these samples with coordinates for concatenated cores. 
   This may need improvements 
3. notebook/pembro_alt_celltypes_prote.Rmd only for samples with targeted panel.
   Output: Cell type labels. This absolutely needs improvements
   
4. script/banksy_domains_5k_pembro.R only for samples with 5k panel. 
   Output: Domains across these samples with coordinates for concatenated cores.
5. script/bansky.R only for samples with 5k panel. Output: Cell types these samples.
   This absolutely needs improvements
6. notebook/pembro_banksyD_prote.Rmd
   Output: Response association, markers and correlation between targeted and 5k domains
   This may need improvements 
7. notebook/pembro_cell_composition.Rmd
   Output: Domain composition, markers 
   
8. notebook/pembro_PD1summary.Rmd
   Output: PD1, PDL1, PanCK, CD45 summary
    This absolutely needs improvements
