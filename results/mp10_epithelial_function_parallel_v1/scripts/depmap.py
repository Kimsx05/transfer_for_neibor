from public_common import *
import re
D=C/'DepMap_24Q4_v1'
def matrixpath(name):return D/name if (D/name).exists() else D/(name+'.xz')
ge=pd.read_csv(matrixpath('CRISPRGeneEffect.csv'),usecols=lambda c:c in ['Unnamed: 0','GPX4 (2879)'] or c=='ModelID')
assert len(ge.columns)==2
ge.columns=['ModelID','response'];assert ge.ModelID.str.startswith('ACH-').all() and ge.ModelID.is_unique
meta=pd.read_csv(D/'Model.csv')
s=[]
for df in pd.read_csv(matrixpath('OmicsExpressionProteinCodingGenesTPMLogp1.csv'),index_col=0,chunksize=48):
 assert df.index.str.startswith('ACH-').all()
 df.columns=df.columns.str.replace(r' \(\d+\)$','',regex=True)
 mat=df.T.groupby(level=0).mean();rank=mat.rank(pct=True);panel=rank.reindex(top);s.append(pd.DataFrame(dict(score=panel.mean(),background=mat.notna().sum(),top_covered=panel.notna().sum())))
sc=pd.concat(s);assert sc.index.is_unique
m=meta.merge(sc,left_on='ModelID',right_index=True,how='left').merge(ge,on='ModelID',how='left')
m['eligible']=m.OncotreeCode.eq('BLCA');m['endpoint']='GPX4 (2879)';m['unit']=m.PatientID.fillna(m.ModelID);m['lineage']=m.OncotreeLineage;m['name']=m.StrippedCellLineName
out=m[m.lineage.eq('Bladder/Urinary Tract')][['ModelID','unit','name','OncotreeSubtype','OncotreeCode','eligible','background','top_covered','score','response']];output(out,'tables/DepMap_urinary_crosswalk.tsv')
a=out[out.eligible];b=a.dropna(subset=['score','response']);output(pd.DataFrame([dict(eligible=len(a),no_expression=a.score.isna().sum(),no_effect=a.response.isna().sum(),both=len(b),units=b.unit.nunique(),score_min=b.score.min(),score_max=b.score.max(),effect_min=b.response.min(),effect_max=b.response.max())]),'tables/DepMap_intersections.tsv')
associations(m,'DepMap')

