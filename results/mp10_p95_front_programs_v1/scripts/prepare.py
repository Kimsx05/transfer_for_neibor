from pathlib import Path
import pandas as pd, shutil
R=Path.cwd();O=R/'results/mp10_p95_front_programs_v1';F=R/'mp10_paper_distance_full_v1'
g=pd.read_csv(F/'tables/gene_sets.tsv.gz',sep='\t'); mp=set('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split())
p30='SERPINE1 TGFBI MMP10 LAMC2 P4HA2 PDPN ITGA5 LAMA3 CDH13 TNC MMP2 EMP3 INHBA LAMB3 VIM SEMA3C PRKCDBP ANXA5 DHRS7 ITGB1 ACTN1 CXCR7 ITGB6 IGFBP7 THBS1 PTHLH TNFRSF6B PDLIM7 CAV1 DKK3'.split()
sets={'PURAM_PEMT_TOP30':p30,'FLORES_CRISPR_PRO_INVASION_10':'ARHGEF19 RWDD1 PRKAR2A EDIL3 LRRC41 ITGA6 PRKD2 NLK CEP19 DNMT3L'.split()}
p=O/'sources/pemt_signature_100_original.txt'
if p.exists():
 hundred=p.read_text().split();assert len(hundred)==100 and hundred[:30]==p30;sets['PURAM_PEMT_TOP100_SENSITIVITY']=hundred
alias={'PRKCDBP':'CAVIN3','CXCR7':'ACKR3','DFNA5':'GSDME','LEPREL1':'P3H2'}
rows=[]
for n,a in sets.items():
 for x in a:rows.append(dict(gs_name=n,original_symbol=x,gene_symbol=alias.get(x,x),source='Puram Cell 2017; Zenodo 3260758' if n.startswith('PURAM') else 'Flores DOI 10.1038/s41420-026-03202-y; user supplied complete positive ten',direction='positive',panel='sensitivity' if '100' in n else 'main',db_version='literature'))
for n,a in g.groupby('gs_name'):
 for x in a.gene_symbol.unique():rows.append(dict(gs_name=n,original_symbol=x,gene_symbol=x,source='MSigDB 2026.1.Hs frozen full_v1',direction='positive',panel='main',db_version='2026.1.Hs'))
 if set(a.gene_symbol)&mp:
  for x in sorted(set(a.gene_symbol)-mp):rows.append(dict(gs_name=n+'__NO_MP10',original_symbol=x,gene_symbol=x,source='Sensitivity: frozen official set minus MP10 Top20',direction='positive',panel='sensitivity',db_version='2026.1.Hs'))
t=pd.DataFrame(rows);t['MP10_overlap']=t.gene_symbol.isin(mp);t.to_csv(O/'tables/gene_sets.tsv',sep='\t',index=False)
t[~t.gs_name.isin(g.gs_name)].to_csv(O/'tables/new_gene_sets.tsv',sep='\t',index=False)
shutil.copyfile(F/'tables/input_manifest.tsv',O/'tables/input_manifest.tsv');shutil.copyfile(R/'neighborhood_niche_visualization_v1/data_tables/niche_palette.tsv',O/'tables/niche_palette.tsv');shutil.copyfile(F/'tables/source_R12_profiles.tsv',O/'tables/source_R12_profiles.tsv')
print(t.groupby('panel').gs_name.nunique(),len(sets))
