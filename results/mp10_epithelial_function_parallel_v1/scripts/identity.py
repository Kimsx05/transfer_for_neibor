from public_common import *
def annotate(m):
 dm=pd.read_csv(C/'DepMap_24Q4_v1/Model.csv')
 # exact stripped-name matches; do not infer identity for ambiguous names.
 dm=dm[~dm.StrippedCellLineName.duplicated(keep=False)]
 m=m.merge(dm[['StrippedCellLineName','ModelID','PatientID','OncotreeCode']],how='left',left_on='name',right_on='StrippedCellLineName',validate='many_to_one')
 urinary=m.lineage.eq('urinary_tract')
 assert m.loc[urinary,'OncotreeCode'].notna().all()
 m['eligible']=urinary&m.OncotreeCode.eq('BLCA')
 m['unit']=m.PatientID.fillna(m.name)
 return m

