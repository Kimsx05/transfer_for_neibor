from base6 import *
import subprocess,pyarrow.parquet as pq
s=sys.argv[1];o=R/s/'native_full';p=o/'in/square_002um/spatial/tissue_positions.parquet';cmd=js(o/'command.json');tr=js(N/s/'metadata/transform.json')
if any(float(cmd[cmd.index(k)+3])<0 for k in ['--cellpose-x-transform','--cellpose-y-transform']):
 t=pq.read_table(p)
 for k in ['pxl_row_in_fullres','pxl_col_in_fullres']:
  j=t.schema.get_field_index(k);t=t.set_column(j,k,pa.array(t.column(j).to_numpy()+10/tr['native_mpp']))
 pq.write_table(t,p)
 for k in ['--cellpose-x-transform','--cellpose-y-transform']:
  i=cmd.index(k)+3;cmd[i]=str(float(cmd[i])+10)
 (o/'native_offset_um.txt').write_text('10');(o/'command.json').write_text(json.dumps(cmd,indent=1))
with (o/'native.log').open('w') as f:sys.exit(subprocess.run(['/usr/bin/time','-v',*cmd],stdout=f,stderr=subprocess.STDOUT).returncode)
