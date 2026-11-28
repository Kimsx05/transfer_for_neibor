"""Add headers to stage resource measurements and summarize a completed run."""
from pathlib import Path
import csv
from datetime import timedelta
def summarize(directory):
    columns = ['stage','exit_code','elapsed_seconds','sampled_tasks_peak','sampled_tree_rss_bytes_peak','sampled_cgroup_bytes_peak','service_kernel_tasks_peak','service_kernel_memory_bytes_peak']
    with (directory/'logs/resources.tsv').open() as f:
        rows = list(csv.DictReader(f, fieldnames=columns, delimiter='\t'))
    with (directory/'tables/resources_by_stage.tsv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=columns,delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    seconds=sum(float(r['elapsed_seconds']) for r in rows)
    mem=max(int(r['service_kernel_memory_bytes_peak']) for r in rows)
    tasks=max(int(r['service_kernel_tasks_peak']) for r in rows)
    rss=max(int(r['sampled_tree_rss_bytes_peak']) for r in rows)
    assert mem<=200000000000 and tasks<=10
    msg=f"按序执行的 {len(rows)} 个阶段用时合计 {timedelta(seconds=round(seconds))}。内核记录的任务线程峰值 {tasks}，cgroup 内存峰值 {mem:,} 字节，0.4 秒采样的进程树 RSS 合计峰值 {rss:,} 字节。RSS 合计可重复计算共享页。"
    out=directory/'REPORT.md'
    if out.exists():
        section='\n## 实测资源\n'
        report=out.read_text().split(section)[0]
        out.write_text(report+section+'\n'+msg+' [逐阶段资源记录](tables/resources_by_stage.tsv)\n')
if __name__=='__main__':
    summarize(Path(__file__).resolve().parents[1])

