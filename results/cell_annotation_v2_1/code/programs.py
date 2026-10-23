# Marker panels used as interpretable expression evidence, not a universal AND veto.
MAJOR={
'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT5','KRT14','UPK1A','UPK2'],
'T/NK':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','NKG7','GNLY','KLRD1'],
'B':['CD79A','CD79B','MS4A1','CD19','CD22','CD37'],
'Plasma':['MZB1','JCHAIN','DERL3','SDC1','TNFRSF17','XBP1'],
'Myeloid':['LST1','LYZ','TYROBP','FCER1G','CSF1R','CD68','AIF1','C1QA','C1QB','C1QC'],
'Fibroblast':['DCN','LUM','COL1A1','COL1A2','COL3A1','PDGFRA','DPT','COL6A1'],
'Mural':['RGS5','CSPG4','PDGFRB','MCAM','NOTCH3','MYH11','CNN1','MYOCD','LMOD1','ACTA2','TAGLN'],
'Endothelial':['PECAM1','VWF','EMCN','KDR','CLDN5','RAMP2','ESAM','FLT1','PROX1'],
'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],
'Erythroid':['HBB','HBA1','HBA2','ALAS2','GYPA','KLF1'],
'Neural_glial':['SOX10','S100B','PLP1','MPZ','SLC1A3','GFAP']}
# Shared ECM, secretory and contractile markers carry less lineage-specific weight.
WEIGHTS={'COL1A1':.15,'COL1A2':.3,'COL3A1':.4,'COL6A1':.3,'XBP1':.2,'SDC1':.4,'CD37':.3,'LYZ':.6,'ACTA2':.35,'TAGLN':.35,'MCAM':.5,'NOTCH3':.5,'RAMP2':.5,'FLT1':.5}
STATE={'Cycling':['MKI67','TOP2A','UBE2C','TYMS','PCNA','STMN1'],'Stress':['FOS','JUN','JUNB','HSPA1A','DDIT3','HSPA1B'],'Interferon':['ISG15','IFIT1','IFIT3','MX1','ISG20'],
'Epithelial_basal':['KRT5','KRT14','KRT17','TP63'],'Epithelial_luminal':['UPK1A','UPK1B','UPK2','UPK3A','KRT20','GATA3','FOXA1'],'Epithelial_squamous':['IVL','DSG3','DSC3','KRT6A','KRT16'],
'Plasma_IgA':['IGHA1','IGHA2','JCHAIN'],'Plasma_IgG':['IGHG1','IGHG2','IGHG3','IGHG4'],'Mast_activation':['HDC','LTC4S','HPGDS','IL1RL1'],
'C1Q_program':['C1QA','C1QB','C1QC'],'SPP1_program':['SPP1','GPNMB','APOC1'],'FOLR2_program':['FOLR2','LYVE1','SELENOP','F13A1']}

# S100A8/A9 can be inflammatory or epithelial background; not sufficient alone for monocyte identity.
STATE['Myeloid_S100A8_A9_program']=['S100A8','S100A9','S100A12']
