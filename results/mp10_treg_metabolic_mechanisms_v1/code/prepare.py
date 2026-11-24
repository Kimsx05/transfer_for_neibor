"""Decode GEO records and preserve the exact accession, probe and abundance joins."""
from pathlib import Path
import csv,gzip,tarfile,io,re,collections,hashlib
out=Path(__file__).resolve().parents[1]; tab=out/'tables'
root=Path('cache/mp10_public_resources_v1')
def inp(acc,pattern):return next((root/acc/'snapshot_20261124').glob(pattern))
def write(name,header,rows):
    p=tab/name
    with (gzip.open(p,'wt') if name.endswith('.gz') else p.open('w')) as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
def bounds(start,end):return ['GSM'+str(i) for i in range(start,end+1)]
cases={'SCAP':('GSE149573',bounds(4505100,4505103),bounds(4505104,4505106)),
       'FASN':('GSE149573',bounds(4505107,4505109),bounds(4505110,4505112)),
       'TCR':('GSE165257',bounds(5029123,5029127),bounds(5029128,5029132)),
       'CD36':('GSE139325',bounds(4138181,4138183),bounds(4138184,4138186)),
       'SREBF2':('GSE243021',bounds(7777514,7777516),bounds(7777517,7777519))}
platforms={'GSE149573':'GPL23038','GSE165257':'GPL16570'}
csv.field_size_limit(20000000)
def parse_gpl(acc):
    res={}
    with gzip.open(inp(acc,'*self.soft.gz'),'rt',errors='replace') as f:
        for line in f:
            if line.strip()=='!platform_table_begin':break
        for row in csv.DictReader(f,delimiter='\t'):
            key=row['ID']
            if key.startswith('!'):break
            field='mrna_assignment' if acc=='GPL23038' else 'gene_assignment'
            text=row[field];syms=set()
            for hit in text.split(' /// '):
                parts=hit.split(' // ')
                if acc=='GPL23038' and len(parts)>2 and parts[1] in ('RefSeq','UCSC Genes'):
                    # Standard symbol directly before transcript/RNA punctuation.
                    m=re.findall(r'\(([A-Za-z0-9][A-Za-z0-9_.-]*)\)(?=\s*[,.;]|$)',parts[2])
                    syms.update(m)
                elif acc=='GPL16570' and len(parts)>1 and parts[1]!='---':syms.add(parts[1].strip())
            res[key]=(next(iter(syms)) if len(syms)==1 else '', ';'.join(sorted(syms)),field)
    return res
metaout=[]; joins=[]; audits=[]
for acc in dict.fromkeys(c[0] for c in cases.values()):
    with gzip.open(inp(acc,'*series_matrix.txt.gz'),'rt') as f:
        md=collections.defaultdict(list); reader=csv.reader(f,delimiter='\t')
        for r in reader:
            if not r:continue
            if r[0]=='!series_matrix_table_begin':break
            if r[0].startswith('!Sample_'):md[r[0][8:]].append(r[1:])
        head=next(reader);arr=[]
        for r in reader:
            if not r:continue
            if r[0]=='!series_matrix_table_end':break
            arr.append(r)
    gsms=md['geo_accession'][0];titles=dict(zip(gsms,md['title'][0]));gpl=platforms.get(acc)
    if gpl:
        assert head[1:]==gsms and (len(arr),len(gsms))==({'GPL23038':(28846,25),'GPL16570':(41345,10)}[gpl])
        annot=parse_gpl(gpl)
        for r in arr:
            sm,opts,field=annot.get(r[0],('','','absent_platform_ID'))
            joins.append((acc,r[0],sm,opts,field,'used' if sm else 'ambiguous' if opts else 'unresolved'))
    sourcecols={}
    for key,(series,a,b) in cases.items():
        if series!=acc:continue
        ids=a+b;rows=[];nraw=0
        if gpl:
            pos=[head.index(x) for x in ids]
            for r in arr:
                sym=annot.get(r[0],('',))[0]
                if sym:rows.append([sym,*[r[j] for j in pos]])
            sourcecols.update({x:x for x in ids});nraw=len(arr)
        elif key=='CD36':
            d={}; names=collections.defaultdict(set)
            with tarfile.open(inp(acc,'*.tar')) as tar:
                for gsm in ids:
                    members=[n for n in tar.getnames() if n.startswith(gsm+'_')];assert len(members)==1
                    with gzip.open(tar.extractfile(members[0]),'rt') as f:
                        rd=csv.reader(f,delimiter='\t');h=next(rd);assert h[:2]==['gene_id','gene_name']
                        sourcecols[gsm]=members[0]+'|'+h[2]
                        for gid,symbol,val in rd:
                            gid=gid.strip();names[gid].add(symbol.strip());dest=d.setdefault(gid,{})
                            assert gsm not in dest
                            dest[gsm]=val
            nraw=len(d)
            for gid,vals in d.items():
                sym=next(iter(names[gid])) if len(names[gid])==1 else ''
                usable=sym and sym!='-' and len(vals)==len(ids)
                joins.append((acc,gid,sym,';'.join(sorted(names[gid])),'gene_id/gene_name','used' if usable else 'missing_or_ambiguous'))
                if usable:rows.append([sym,*[vals[i] for i in ids]])
        else:
            with gzip.open(inp(acc,'*.csv.gz'),'rt') as f:
                rd=csv.DictReader(f)
                cols=[titles[g]+'.rpkm' for g in ids];assert set(cols)=={k for k in rd.fieldnames if k.startswith('VAT.')}
                sourcecols.update(dict(zip(ids,cols)))
                for r in rd:
                    nraw+=1;sym=r['SYMBOL'].strip();ok=sym not in ('','NA','---')
                    joins.append((acc,r['ENTREZID'],sym,sym,'ENTREZID/SYMBOL','used' if ok else 'missing_symbol'))
                    if ok:rows.append([sym,*[r[j] for j in cols]])
        write('expression_'+key+'.tsv.gz',['mouse_symbol',*ids],rows)
        audits.append((key,acc,nraw,len(rows),len(set(r[0] for r in rows)),len(a),len(b)))
    for j,g in enumerate(gsms):
        match=[(k,'control' if g in a else 'case') for k,(series,a,b) in cases.items() if g in a+b]
        k,arm=match[0] if match else ('excluded','')
        chars=' | '.join(v[j] for v in md['characteristics_ch1'] if v[j])
        if k in ('SCAP','FASN'):
            assert 'tumor infiltrating Treg' in chars and ('genotype: WT' if arm=='control' else ('Scap KO' if k=='SCAP' else 'Fasn KO')) in chars
        if k=='TCR':assert titles[g].startswith('0h' if arm=='control' else '8h')
        if k=='CD36':assert 'melanoma' in chars and ('genotype: WT' if arm=='control' else 'genotype: CD36 KO') in chars
        if k=='SREBF2':assert titles[g].startswith('VAT.WT' if arm=='control' else 'VAT.KO')
        metaout.append((acc,g,k,arm,titles[g],md['source_name_ch1'][0][j],chars,' | '.join(v[j] for v in md['description']),sourcecols.get(g,''), 'HOLD spleen' if acc=='GSE243021' and not match else 'included' if match else 'out of scope'))
write('sample_mapping.tsv',['accession','GSM','contrast','arm','title','source','characteristics','description','input_column','decision'],metaout)
write('input_gene_mapping.tsv.gz',['accession','source_id','mouse_symbol','candidate_symbols','annotation_field','status'],joins)
write('platform_coverage.tsv',['contrast','accession','rows_submitted','mapped_features','unique_symbols','n_control','n_case'],audits)
groups=collections.defaultdict(lambda:collections.defaultdict(set))
with inp('MGI','*.rpt').open() as f:
    for r in csv.DictReader(f,delimiter='\t'):
        groups[r['DB Class Key']][r['NCBI Taxon ID']].add((r['Symbol'],r['EntrezGene ID']))
pairs=[]; usage=collections.defaultdict(set)
for group,org in groups.items():
    a=org['10090'];b=org['9606']
    for x in a:usage['M',x[0]].add(group)
    for x in b:usage['H',x[0]].add(group)
for group,org in groups.items():
    a=org['10090'];b=org['9606']
    ok=len(a)==len(b)==1 and len(usage['M',next(iter(a))[0]])==len(usage['H',next(iter(b))[0]])==1
    for m in a or {('','')}:
        for h in b or {('','')}:
            pairs.append((group,*m,*h,len(a),len(b),'one_to_one' if ok else 'missing_species' if not (a and b) else 'non_unique'))
write('gene_mapping.tsv.gz',['MGI_class','mouse','mouse_entrez','human','human_entrez','n_mouse','n_human','status'],pairs)
print(audits)
