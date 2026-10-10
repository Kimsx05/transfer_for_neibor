# Execution notes

- Original v5 scope/run are read-only. New work is diagnostic v6, following the explicitly invoked attachment.
- User later explicitly authorized resource-aware parallel model execution. This overrides the attachment's single-GPU sequential scheduling sentence. At most3 model processes were used, with existing task detection and resource checks. No statistical settings changed.
- The scheduler initially matched only absolute command paths. It was corrected to recognize both relative/absolute paths and existing tags, then only the scheduler was restarted; running fits were not stopped or duplicated. The maximum remains1 NB+5mapping fits.
- Old quantities require weighted moments before nonlinear RMSE/calibration. An initial analysis draft averaging sample slopes/RMSE did not reproduce v5; the code was corrected, and all44 historical reproduction checks passed. This did not alter counts, thresholds or models. Final tables use the corrected method.
- Pandas3 nullable/string assignment and read-only NumPy array behavior required local compatibility corrections in audit/evaluation scripts only. No environment reinstall.
- Official GEO deposited feature annotations were available; official10x probe CSV requests were unavailable (403/connection closure), and the download-page request returned429. These failures are recorded rather than interpreted as assay absence. No bulk spatial expression was downloaded.
- Newly generated control counts are not claimed to numerically reproduce old fullDEV predictions, even though the same cell bags are reused.
- Intermediate plots existed before all S30 branches were finished and explicitly marked unavailable panels. The final automated sequence regenerates all plots from completed outputs before delivery validation.
- Final outcome interpretation is in08_review/INTERPRETATION.json; report/package generation cannot begin without passing delivery checks. G1 and real spatial remain not run by design.
