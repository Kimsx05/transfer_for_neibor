from overview import combine_pdf
from report import *
files=[]
for sec in ['GSE319536_BC12','GSE319536_BC13','GSE319536_BC18']:
 subprocess.run(['/usr/bin/python3.12',str(O/'scripts/redraw.py'),'--section',sec,'--candidate','vCAF','--program','BUDDEN_TNF_EMT_RESPONSE_15','--range','common_support_0_1000'],check=True)
 files.append(O/'figures/redrawn'/f'{sec}__vCAF__BUDDEN_TNF_EMT_RESPONSE_15__AddModuleScore_SCT_data__common_support_0_1000.pdf')
combine_pdf(files,O/'figures/discussion_GSE319536_vCAF_TNF.pdf')
