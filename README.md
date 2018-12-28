# Xena_wrapper

This is a Python class for wrapping xenaPython

Copyright of Xena : The Regents of the University of California, Santa Cruz.
[UCSC Xena](https://xena.ucsc.edu/)

# Example :
```Python
import xenaPython as xena
from TCGAlib.TCGAlib import load_TCGA

### Host (GDC TCGA)
host = xena.PUBLIC_HUBS['gdcHub']
### Project name
prjID = 'TCGA-SKCM'
### Cohort name
cohort = ['GDC TCGA Melanoma (SKCM)']

tcga = load_TCGA(host, cohort, prjID)

# Input shoud be EntrezID
input_list=['9480', '367', '2137']

# Getting expression(FPKM)
tcga_expr = tcga.get_TCGA_expr(input_list)
# Getting CNV
tcga_cnv = tcga.get_TCGA_cnv(input_list)
# Getting Mutation
tcga_mut = tcga.get_TCGA_mut(input_list)
# Getting Survival
tcga_surv = tcga.get_TCGA_surv()
```
