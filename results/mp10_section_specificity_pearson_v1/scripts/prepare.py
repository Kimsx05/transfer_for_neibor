import os,sys
from pathlib import Path
R=Path(__file__).resolve().parents[3];sys.path.insert(0,str(R/'spatial_ecology_top20_v2/runtime'))
import pandas as pd,numpy as np,shutil
O=R/'results/mp10_section_specificity_pearson_v1';OLD=R/'results/mp10_p95_front_programs_v1';F=R/'mp10_paper_distance_full_v1'
g=pd.read_csv(OLD/'tables/gene_sets.tsv',sep='\t');sets={
'BUDDEN_TG_CATABOLISM_24':'ABHD5 CAV1 FABP1 FABP12 FABP2 FABP3 FABP4 FABP5 FABP6 FABP7 FABP9 GPD2 LIPE MGLL PLIN1 PLIN3 PNPLA4 PNPLA5 PPP1CA PPP1CB PPP1CC PRKACA PRKACB PRKACG',
'BUDDEN_CHOLESTEROL_SYNTHESIS_27':'ACAT2 ARV1 CYP51A1 DHCR24 DHCR7 EBP FDFT1 FDPS GGPS1 HMGCR HMGCS1 HSD17B7 IDI1 IDI2 LBR LSS MSMO1 MVD MVK NSDHL PLPP6 PMVK SC5D SQLE SREBF1 SREBF2 TM7SF2',
'BUDDEN_S1P_RESPONSE_18':'CCN2 CD274 CXCL2 EGR1 ETS1 FN1 FOS IL6 JUN JUNB MMP9 NCAM1 PIEZO2 PTGS2 RXRA SERPINE1 TAGLN THBS1',
'BUDDEN_TNF_EMT_RESPONSE_15':'CXCL2 SQSTM1 PRSS2 IER2 IL1B IL1A FERMT2 ZFP36 JUN JUNB AREG DUSP1 FOSB SERPINE1 FOS',
'BUDDEN_CHOLESTEROL_SYNTHESIS_noMP10_18':'ARV1 CYP51A1 DHCR24 DHCR7 EBP GGPS1 HSD17B7 IDI2 LBR LSS MVK NSDHL PLPP6 PMVK SC5D SREBF1 SREBF2 TM7SF2'}
mp=set('HMGCS1 MSMO1 INSIG1 FDFT1 SQLE FDPS IDI1 ACSS2 ELOVL5 LDLR MVD SCD HMGCR ACAT2 SLC20A1 ACSL3 ACSL1 ADIPOR2 MXD1 NPC1'.split());a=[]
for name,gs in sets.items():
 for gene in gs.split():a.append(dict(gs_name=name,original_symbol=gene,gene_symbol=gene,source='Budden Nature Metabolism 2026; Zenodo 18790615 original scripts',direction='positive',panel='sensitivity' if 'noMP10' in name else 'main',db_version='literature',MP10_overlap=gene in mp))
g=pd.concat([g,pd.DataFrame(a)],ignore_index=True);assert g.gs_name.nunique()==85;assert len(set(sets['BUDDEN_CHOLESTEROL_SYNTHESIS_27'].split())&mp)==9
g.to_csv(O/'tables/gene_sets.tsv',sep='\t',index=False);pd.DataFrame(a).to_csv(O/'tables/new5_gene_sets.tsv',sep='\t',index=False)
m=pd.read_csv(F/'tables/input_manifest.tsv',sep='\t');m.to_csv(O/'tables/input_manifest.tsv',sep='\t',index=False)
for sec in m.section_id:
 d=pd.read_parquet(F/'spots'/f'{sec}.parquet');means=d.filter(regex='^mean__');valid=(means.ge(0)&np.isfinite(means)).all(1)&means.sum(1).gt(0)&d.center_qc&d.lesion_location.ne('quality_excluded');pd.DataFrame({'barcode':d.barcode,'valid_tissue':valid}).to_csv(O/'scores'/f'{sec}_eligibility.tsv',sep='\t',index=False)
factors=sorted(c[6:] for c in d if c.startswith('mean__') and c!='mean__Epithelial');assert len(factors)==29
pd.DataFrame({'candidate_type':factors,'original_mean_column':['mean__'+c for c in factors],'fraction_denominator':'all30 posterior mean'}).to_csv(O/'tables/candidate_mapping.tsv',sep='\t',index=False)
shutil.copyfile(OLD/'tables/niche_palette.tsv',O/'tables/niche_palette.tsv');shutil.copyfile(OLD/'tables/source_R12_profiles.tsv',O/'tables/source_R12_profiles.tsv')
print(g.groupby('panel').gs_name.nunique())
v=pd.read_csv(O/'tables/source_R12_profiles.tsv',sep='\t');pal=pd.read_csv(O/'tables/niche_palette.tsv',sep='\t').set_index('accepted_label').color;key=[]
for k in range(9):
 z=v[(v.cluster==k)&v.ring.eq('R1')&v.cell_type.ne('Epithelial')].nlargest(3,'mean_relative_composition');key.append(dict(formal_label=k,formal_column='accepted_label (source SAMPLE_R12_LEIDEN)',old_display='R12_niche'+str(k+1),color=pal[k],R1_largest_nonEpi_factors='; '.join(f'{r.cell_type}: {r.mean_relative_composition:.3f}' for r in z.itertuples()),interpretation='original R1 profile; descriptive, not a relabelling or center posterior composition'))
pd.DataFrame(key).to_csv(O/'tables/niche_label_key.tsv',sep='\t',index=False)
