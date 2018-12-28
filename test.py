import xenaPython as xena
from TCGAlib.TCGAlib import load_TCGA


host = xena.PUBLIC_HUBS['gdcHub']

"""
#####finding cohort here.
prj_arr = [x['name'].split('/')[0] for x in xena.all_datasets(host) if x['name'].split('/')[0].find('TCGA')!=-1]
prj_arr = list(set(prj_arr))
prj_arr_el = [x.split('-') for x in prj_arr]
cohort_arr = [y for i,item in enumerate(prj_arr_el) for y in xena.all_cohorts(host, exclude=['genomicVector']) if y.find(item[0])>-1 and y.find(item[1])>-1]
"""


host = xena.PUBLIC_HUBS['gdcHub']
prjID = 'TCGA-SKCM'
cohort = ['GDC TCGA Melanoma (SKCM)']

tcga = load_TCGA(host, cohort, prjID)

# return df, FPKM(NOT FPKM-UQ)
tcga_expr = tcga.get_TCGA_expr(input_list=['9480', '367', '2137'])

# return df, there are duplicated samples, but take the first.
tcga_cnv = tcga.get_TCGA_cnv(input_list=['9480', '367', '2137'])

# return array with df as elements, not merged df, because it has duplication.
tcga_mut = tcga.get_TCGA_mut(input_list=['9480', '367', '2137'])

# return df, Survival data
tcga_surv = tcga.get_TCGA_surv()
