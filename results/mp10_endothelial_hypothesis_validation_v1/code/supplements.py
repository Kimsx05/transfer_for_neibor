from fetch import *
import subprocess
base='https://pmc-oa-opendata.s3.us-east-1.amazonaws.com/'
targets={'PMC8432967':['mmc1.pdf'],'PMC5380221':['NIHMS66220-supplement-1.pdf'],'PMC6521885':['NIHMS1511426-supplement-1.pdf','NIHMS1511426-supplement-2.pdf'],'PMC10515159':['pnas.2218150120.sapp.pdf'],'PMC10407908':['mmc1.xlsx']}
def parsed(f):
 if f.suffix=='.pdf':
    r=subprocess.run(['pdftotext','-layout',str(f),'-'],capture_output=True,check=True);txt=r.stdout.decode(errors='replace')
    return len(txt), ' | '.join(s.strip()[:100] for s in txt.splitlines() if 'Fig' in s or 'Table' in s)[:480]
 import xml.etree.ElementTree as ET
 with zipfile.ZipFile(f) as z:
    root=ET.fromstring(z.read('xl/workbook.xml'));sheets=[s.attrib['name'] for s in root.iter() if s.tag.endswith('}sheet')]
    return len(sheets),' | '.join(sheets)
for pm,names in targets.items():
 for name in names:
    url=base+pm+'.1/'+name;fetch(pm,name,[url+'?download=1',url],name.split('.')[-1])
rows=[]
for pm in list(targets)+['PMC9395194','eLife.63326']:
 for f in (C/pm).glob('*'):
    if f.suffix in ('.xlsx','.pdf'):
        try:n,s=parsed(f);rows.append((pm,f.name,n,s))
        except Exception as e:rows.append((pm,f.name,0,repr(e)))
with (P/'tables/supplement_readout.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(('accession','file','text_chars_or_worksheet_count','figures_or_sheets'));w.writerows(rows)

