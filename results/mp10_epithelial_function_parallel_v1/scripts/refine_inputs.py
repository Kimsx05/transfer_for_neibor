from pathlib import Path
o=Path(__file__).resolve().parent
p=o/'project.R';s=p.read_text()
s=s.replace("ww=w[source==m & target %in% rownames(mm)];s=apply(mm[ww$target,,drop=FALSE],2,function(v)coef(summary(lm(v~ww$mor)))[2,3])","ww=w[source==m & target %in% rownames(mm)];act=decoupleR::run_ulm(mm,ww,.source='source',.target='target',.mor='mor',minsize=10);s=act$score[match(p$profile_id,act$condition)]")
marker="scores=list(); genes=list()"
new="""targets=unique(c(unlist(sets),'SREBF2','SCAP','CHMP4B','CHMP2A'))
selected=intersect(targets,rownames(obj$counts))
node_counts=as.data.table(as.table(as.matrix(obj$detected_cells[selected,pr$profile_id,drop=FALSE])))
setnames(node_counts,c('gene','profile_id','n_detected')); node_counts=merge(node_counts,pr,by='profile_id')
node_counts[,det_fraction:=n_detected/n_cells];saveT(node_counts,'tables/candidate_profile_detection.tsv.gz')
"""
if new not in s:s=s.replace(marker,new+marker)
p.write_text(s)
p=o/'depmap.py';s=p.read_text()
if 'def matrixpath' not in s:
 s=s.replace("ge=pd.read_csv(D/'CRISPRGeneEffect.csv'","def matrixpath(name):return D/name if (D/name).exists() else D/(name+'.xz')\nge=pd.read_csv(matrixpath('CRISPRGeneEffect.csv')").replace("in pd.read_csv(D/'OmicsExpressionProteinCodingGenesTPMLogp1.csv',","in pd.read_csv(matrixpath('OmicsExpressionProteinCodingGenesTPMLogp1.csv'),")
p.write_text(s)

