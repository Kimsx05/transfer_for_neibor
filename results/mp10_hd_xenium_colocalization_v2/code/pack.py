from common import *
from zipfile import ZipFile,ZIP_DEFLATED
import html,os
def pack():
 targets=[];root=OUT/'review';archive=root/'first_batch.zip';tmp=root/'first_batch.tmp.zip'
 with ZipFile(tmp,'w',compression=ZIP_DEFLATED,compresslevel=4) as f:
  f.write(root/'README.md','README.md')
  f.writestr('index.html','<meta charset="utf-8"><h2>HD 与 Xenium 盲审 · human pending</h2><p>先看局部计数与组织，再打开复核表填写。图像隐去所有代谢状态。</p><ul>'+''.join('<li><a href="'+s+'/index.html">'+s+'</a></li>' for s in HD+XE)+'</ul>')
  for s in HD+XE:
   q=paths(s);hu=pd.read_csv(q/'human_review.tsv',sep='\t');w=pd.read_csv(q/'window_index.tsv',sep='\t');f.write(q/'human_review.tsv',s+'/human_review.tsv')
   body='<meta charset="utf-8"><h2>'+s+' · pending</h2><a href="human_review.tsv">填写表</a><ol>'
   for i in hu.id:
    for ext in ['png','pdf']:
     file=q/'cards'/(i+'.'+ext);f.write(file,s+'/cards/'+file.name)
    body+='<li><a href="cards/'+i+'.png">'+html.escape(i)+'</a> · <a href="cards/'+i+'.pdf">PDF</a>'
    if s in HD:
     file=q/'histology'/(i+'.png');f.write(file,s+'/histology/'+file.name);body+=' · <a href="histology/'+i+'.png">原图</a>'
    body+='</li>'
   f.writestr(s+'/index.html',body+'</ol>');targets.append((s,len(hu),len(w)))
 with ZipFile(tmp) as f:assert f.testzip() is None
 os.replace(tmp,archive);put(pd.DataFrame(targets,columns=['sample','human_units','all_indexed_units']),root/'packet_manifest.tsv')
if __name__=='__main__':pack()
