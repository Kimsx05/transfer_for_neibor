"""Extract published figures without using pixels to derive spot labels."""
from pathlib import Path
import zipfile,xml.etree.ElementTree as ET
O=Path(__file__).resolve().parents[1]
p=O/'sources/SpatialBC_supplement.docx'
with zipfile.ZipFile(p) as z:
 for name in ['image1.png','image9.png']:(O/'sources'/name).write_bytes(z.read('word/media/'+name))
 tree=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 (O/'sources/SpatialBC_supplement_text.txt').write_text('\n'.join(''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in tree.findall('.//w:p',ns)))
