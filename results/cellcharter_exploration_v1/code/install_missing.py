import sys,subprocess,importlib.metadata as md
from pathlib import Path
from packaging.requirements import Requirement
root=Path(__file__).resolve().parents[1];target=root/'.packages';sys.path.insert(0,str(target))
for step in range(10):
 todo=[];visited=set();queue=['cellcharter','squidpy']
 while queue:
  name=queue.pop()
  if name.lower() in visited:continue
  visited.add(name.lower())
  try:d=md.distribution(name)
  except md.PackageNotFoundError:todo.append(name);continue
  for item in d.requires or []:
   r=Requirement(item)
   if r.marker and not r.marker.evaluate({'extra':''}):continue
   try:v=md.version(r.name)
   except md.PackageNotFoundError:todo.append(str(r));continue
   if r.specifier and v not in r.specifier:todo.append(str(r))
   else:queue.append(r.name)
 if not todo:print('DEPS_COMPLETE');break
 todo=sorted(set(todo));print('INSTALL',todo,flush=True);subprocess.run([sys.executable,'-m','pip','install','--target',str(target),'--upgrade','--no-deps',*todo],check=True);md.MetadataPathFinder.invalidate_caches()
