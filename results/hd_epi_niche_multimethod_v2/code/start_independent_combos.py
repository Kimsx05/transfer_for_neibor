from night_supervise import run,ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=3) as pool:
 futures=[pool.submit(run,f'fit_combo_{ref}_{r}',offset) for ref,r,offset in [('blca',200,4),('hpca',50,8),('hpca',200,12)]]
 for f in futures:f.result()
