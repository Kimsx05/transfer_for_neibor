from common import *
import shutil,zipfile,xml.etree.ElementTree as ET
src=Path('/home/data/t070721/.codex/attachments/b1ceba5e-6188-4736-ad25-0bc96cbf74f0');paths=[src/'niche_celltype_biology_review_v6.xlsx',src/'niche_celltype_review_suggestions_270.tsv']
for path in paths:shutil.copy2(path,OUT/'00_approval'/path.name)
s=pd.read_csv(paths[1],sep='\t');p=pairs();chosen=s[s.review_tier.isin(TIERS)];assert len(s)==270 and len(chosen)==139 and s.candidate_id.is_unique
j=p.merge(chosen[['candidate_id','niche','accepted_label','celltype','review_tier']],left_on='master_candidate_id',right_on='candidate_id',suffixes=('_p','_s'),validate='one_to_one');assert len(j)==139
for c in ['niche','accepted_label','celltype','review_tier']:assert j[c+'_p'].eq(j[c+'_s']).all(),c
# Fully parse every workbook sheet; preserve all sheet contents as TSV for auditable review.
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};workbook=[]; workbook_pairs={}
with zipfile.ZipFile(paths[0]) as z:
 shared=[]
 if 'xl/sharedStrings.xml' in z.namelist():shared=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)]
 wb=ET.fromstring(z.read('xl/workbook.xml'));rels=ET.fromstring(z.read('xl/_rels/workbook.xml.rels'));targets={r.attrib['Id']:r.attrib['Target'] for r in rels};allcandidates=set()
 for sh in wb.find('m:sheets',ns):
  rid=sh.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'];target=targets[rid];target=target.lstrip('/') if target.startswith('/') else 'xl/'+target
  root=ET.fromstring(z.read(target));rows=[]
  for row in root.findall('.//m:sheetData/m:row',ns):
   out={}
   for c in row.findall('m:c',ns):
    name=re.sub(r'\d','',c.attrib['r']) if False else ''.join(ch for ch in c.attrib['r'] if ch.isalpha());col=0
    for ch in name:col=col*26+ord(ch)-64
    val=c.find('m:v',ns);v=val.text if val is not None else ''
    if c.attrib.get('t')=='s':v=shared[int(v)]
    elif c.attrib.get('t')=='inlineStr':v=''.join(c.find('m:is',ns).itertext())
    out[col-1]=v
    if str(v).startswith('V6_'):allcandidates.add(v)
   rows.append(out)
   if str(out.get(0,'')).startswith('niche') and out.get(1) in TYPES and out.get(2) in ['核心正向','宽口径正向','补充探索','负向对照']:
    workbook_pairs[(out[0],out[1])]=out[2]
  data=pd.DataFrame(rows).reindex(columns=range(max((max(r,default=0) for r in rows),default=0)+1));safe='sheet_'+str(len(workbook)+1);tab('00_approval/workbook_parsed/'+safe+'.tsv',data);workbook.append(dict(sheet_name=sh.attrib['name'],n_rows=len(rows),n_columns=len(data.columns),parsed_file=safe+'.tsv'))
 translation=dict(POSITIVE_CORE='核心正向',POSITIVE_BROAD='宽口径正向',EXPLORATORY='补充探索',DEPLETION_CONTROL='负向对照')
 for row in p.itertuples():assert workbook_pairs.get((row.niche,row.celltype))==translation[row.review_tier],(row.niche,row.celltype,workbook_pairs.get((row.niche,row.celltype)))
tab('00_approval/workbook_sheet_inventory.tsv',workbook);tab('00_approval/supplemental_input_manifest.tsv',[input_record(p,'user supplemental review source') for p in paths]);js('00_approval/supplemental_crosscheck.json',dict(status='PASS',all_270_tsv_rows_read=True,approved_139_match_all_keys_and_tiers=True,workbook_all_sheets_read=True,workbook_all_139_niche_celltype_tiers_match=True,workbook_ID_column_present=False,original_Gate_remains_unchanged=True,authorization='current explicit user request, not blank decision fields in recommendation table'))
print(json.dumps(workbook,ensure_ascii=False));print('Supplemental tables agree with all 139 approved appendix pairs and tiers')
