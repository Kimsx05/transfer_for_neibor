"""Observed-expression review rules, using frozen v2 marker frameworks.
These are annotation evidence requirements, not cell-retention filters.
"""
MAJOR={
'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT5','KRT14','UPK1A','UPK2'],
'B':['CD79A','CD79B','MS4A1','CD19','CD22','CD37'],
'Plasma':['MZB1','JCHAIN','SDC1','XBP1','DERL3','TNFRSF17'],
'T/NK':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','NKG7','GNLY','KLRD1','KLRF1'],
'Myeloid':['LST1','LYZ','TYROBP','FCER1G','CSF1R','CD68','AIF1','FCGR3A'],
'Fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','RGS5','ACTA2','TAGLN','CSPG4','PDGFRB','MYH11'],
'Endothelial':['PECAM1','VWF','EMCN','KDR','CLDN5','RAMP2','ESAM','FLT1','PROX1'],
'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2']}
SUBTYPE={
'Naive B':['IGHD','TCL1A','FCER2'], 'Memory B':['CD27','TNFRSF13B','AIM2'], 'GC B':['BCL6','AICDA','RGS13'],
'Classical Mono':['FCN1','VCAN','S100A8','S100A9','S100A12','CCR2'],
'FOLR2 TAM':['FOLR2','LYVE1','SELENOP','F13A1'], 'C1QC TAM':['C1QA','C1QC','APOE'], 'SPP1 TAM':['SPP1','GPNMB','FN1'],
'cDC1':['CLEC9A','XCR1','CADM1','BATF3'], 'cDC2':['CD1C','FCER1A','CLEC10A'], 'LAMP3 DC':['LAMP3','CCR7','FSCN1'], 'pDC':['LILRA4','CLEC4C','GZMB','IL3RA'], 'Neutrophil':['CSF3R','FCGR3B','CXCR2'],
'Treg':['FOXP3','IL2RA','CTLA4','IKZF2'], 'CD8 Exhausted':['CD8A','CD8B','TOX','ENTPD1','LAG3','HAVCR2','ITGAE'],
'CD8 GZMK/Effector-memory-like':['GZMK','CD8A','CD8B','CCL5'], 'CD4 Helper/Conventional':['IL7R','KLRB1','CD40LG','RORA','IL17A','IL26','IL23R'],
'CD4 Naive/Memory-like':['CCR7','TCF7','LEF1','SELL','KLF2','IL7R'], 'Tfh/CXCL13-like':['CXCL13','IL21','TOX2','MAF','CD4','PDCD1','ICOS'], 'NK':['FGFBP2','FCGR3A','GNLY','GZMB','PRF1','KLRD1','NKG7'],
'iCAF':['PI16','DPT','PDGFRA','CXCL12','CFD','APOD','PLA2G2A'], 'Matrix_CAF':['COL11A1','COL10A1','MMP11','POSTN','CTHRC1','LRRC15','FAP'], 'rCAF':['CCL19','CCL21'], 'Pericyte':['RGS5','CSPG4','PDGFRB','MCAM','NOTCH3'], 'SMC':['MYH11','CNN1','MYOCD','LMOD1','DES','ACTA2','TAGLN'], 'vCAF':['MCAM','NOTCH3','COL18A1'],
'Venous EC':['ACKR1','SELE','VCAM1'], 'Capillary EC':['CA4','CD36','RGCC'], 'Arterial EC':['GJA5','EFNB2','HEY1'], 'Lymphatic EC':['PROX1','LYVE1','CCL21']}
# Additional identity anchors prevent generic activation/cytotoxic genes alone from supporting a fine identity.
ANCHORS={'GC B':['BCL6','AICDA','RGS13'],'Naive B':['IGHD','TCL1A','FCER2'],'Memory B':['CD27','TNFRSF13B','AIM2'],'Tfh/CXCL13-like':['CXCL13','IL21','TOX2'],'CD8 Exhausted':['TOX','ENTPD1','LAG3','HAVCR2'],'CD8 GZMK/Effector-memory-like':['GZMK'],'Treg':['FOXP3','IL2RA'],'NK':['GNLY','KLRD1','FGFBP2'],'SPP1 TAM':['SPP1'],'FOLR2 TAM':['FOLR2','LYVE1'],'C1QC TAM':['C1QA','C1QB','C1QC'],'vCAF':['COL18A1'],'Pericyte':['RGS5','CSPG4','PDGFRB'],'SMC':['MYH11','CNN1','MYOCD','LMOD1'],'Venous EC':['ACKR1'],'pDC':['LILRA4','CLEC4C'],'Neutrophil':['CSF3R','FCGR3B','CXCR2']}
