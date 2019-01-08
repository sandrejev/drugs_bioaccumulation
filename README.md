---
output:
  html_document: default
  pdf_document: default
---
## Data files *

* **drug_map.csv**
* **bug_map.csv**

* **170511_processingUPLCpeaks_data.clean_aftergrowthmerge_andfixingget.datacleanfunction.RData** - 
depletion degradation of drugs, experiment #1, SANKEY LEFT

* **curves.rel_annotation_2016-11-28.tab** - 
growth data influenced by drugs, experiment #1, SANKEY RIGHT 
(From /g/patil/Martina/Documents/computational stuff/Drug Growth Curve analysis Sergej/161128_growthcurves anno/)

* **data.depletionmodeassay_long.csv**, experiment #2, SANKEY LEFT (color)

#
#
#

We use p-value 0.05 (instead 0.01) and bioaccumulation fold threshold 0.7 (if MedianDiff < -0.1) and 0.1 delta threshold

        species.short      drug.long  padjst.super  padjst.total  MedianDiff.super  MedianDiff.total        interaction          Activity
             C. comes    Montelukast  !0.645849342    0.03935458        -0.3310672        -0.1328072  Biotransformation   Bioaccumulation
           E. rectale    Montelukast  !0.280667695    0.00527550        -0.1189074        -0.1343513  Biotransformation   Bioaccumulation
         F. nucleatum     Ranitidine  !0.172219085    0.00527550         0.1607430        -0.2877953  Biotransformation       No activity
            L. lactis    Roflumilast   0.005275500    0.00527550       !-0.6122437       !-0.2909171    Bioaccumulation       No activity
                 <NA>     Duloxetine   0.010252040    0.00527550        -1.0000000        -1.0000000  Biotransformation       No activity
                 <NA>  Sulfasalazine   0.005275500    0.00527550        -1.0000000        -1.0000000  Biotransformation       No activity

                   old
new                 Bioaccumulation Biotransformation No activity Previously known
  Bioaccumulation                19                 0           1                0
  Biotransformation               2                17           3                9
  No activity                     0                 0          15                0
  
  
Ranitidine+L. gasseri not in original file

##############
# Some drugs/species are never tested for degradation
##############