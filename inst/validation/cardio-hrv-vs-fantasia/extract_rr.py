"""Extract a resting RR-interval series from a PhysioNet Fantasia record.

Reads the beat annotations of Fantasia young-subject record f1y03, forms the
RR-interval series (normal beats), and writes a 600-beat resting segment with
provenance to artifacts/ for the R analysis.

Portable: no absolute paths. Point FANTASIA_F1Y03 at a local record prefix (.../f1y03),
else the annotations are fetched from PhysioNet (https://physionet.org/content/fantasia/).
"""
import os, csv, json
import numpy as np
import wfdb

here = os.path.dirname(os.path.abspath(__file__))
art = os.path.join(here, "artifacts"); os.makedirs(art, exist_ok=True)
local = os.environ.get("FANTASIA_F1Y03")
ann = wfdb.rdann(local, "ecg") if local else wfdb.rdann("f1y03", "ecg", pn_dir="fantasia")
fs = 250.0
rpk = np.asarray(ann.sample)[np.asarray(ann.symbol) == "N"]
rr = np.diff(rpk) / fs * 1000.0                 # RR intervals (ms)
rr = rr[(rr > 400) & (rr < 1400)]               # physiological
seg = rr[500:1100]                              # a 600-beat resting segment

with open(os.path.join(art, "rr.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["rr_ms"]); w.writerows([[float(v)] for v in seg])
json.dump(dict(dataset="Fantasia Database", subject="f1y03 (young healthy)",
               fs=fs, n_beats=int(len(seg)),
               rr_mean_ms=round(float(seg.mean()), 1), rr_sd_ms=round(float(seg.std()), 1)),
          open(os.path.join(art, "real_meta.json"), "w"), indent=2)
print(f"f1y03: {len(seg)} beats, mean {seg.mean():.0f} ms, sd {seg.std():.1f} ms")
