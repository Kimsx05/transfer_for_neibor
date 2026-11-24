from fetch import *
base='https://pmc-oa-opendata.s3.us-east-1.amazonaws.com/'
rows=[]
for pm in ('PMC8432967','PMC7997658','PMC5380221','PMC9395194','PMC10407908','PMC10515159','PMC6521885'):
    f=fetch(pm,'article.xml',[base+pm+'.1/'+pm+'.1.xml',base+pm+'.1/'+pm+'.1.xml?download=1','https://s3.amazonaws.com/pmc-oa-opendata/'+pm+'.1/'+pm+'.1.xml'],'xml')
    if not f:continue
    root=ET.parse(f).getroot()
    for s in root.iter():
        if s.tag not in ('supplementary-material','media'):continue
        caption=' '.join(''.join(s.itertext()).split())[:350]
        for el in s.iter():
            href=el.attrib.get('{http://www.w3.org/1999/xlink}href','')
            if href.endswith(('.xlsx','.xls','.pdf','.docx','.zip')):rows.append((pm,href,caption))
with (P/'tables/literature_supplement_index.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(('PMCID','file','caption'));w.writerows(sorted(set(rows)))
for pm,name,cap in set(rows):
    if pm in ('PMC10407908','PMC10515159','PMC9395194') and name.endswith('.xlsx'):
        fetch(pm,name,[base+pm+'.1/'+name],'xlsx')

