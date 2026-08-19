"""Extract an eyes-closed resting occipital EEG channel from PhysioNet eegmmidb.

Downloads run R02 (eyes closed) of subject S001 of the EEG Motor Movement/Imagery
Database, takes an occipital channel (O1) and writes it, with the sampling rate and
provenance, to artifacts/ for the R analysis. A cross-check alpha summary
(dominant frequency, relative alpha power via Welch) is also recorded.

Portable: no absolute paths. Point EEGMMIDB_S001R02 at a local S001R02.edf, else it
is fetched from PhysioNet (https://physionet.org/content/eegmmidb/).
"""
import os, csv, json, tempfile, urllib.request
import numpy as np
from scipy import signal
import mne

here = os.path.dirname(os.path.abspath(__file__))
art = os.path.join(here, "artifacts"); os.makedirs(art, exist_ok=True)
edf = os.environ.get("EEGMMIDB_S001R02") or os.path.join(tempfile.gettempdir(), "S001R02.edf")
if not os.path.exists(edf):
    url = "https://physionet.org/files/eegmmidb/1.0.0/S001/S001R02.edf"
    urllib.request.urlretrieve(url, edf)

raw = mne.io.read_raw_edf(edf, preload=True, verbose="ERROR")
fs = float(raw.info["sfreq"])
ch = "O1.."                                         # occipital channel (10-10 label in eegmmidb)
x = raw.get_data(picks=ch)[0] * 1e6                 # microvolts

def welch(sig_, fs, seg=8):
    f, P = signal.welch(sig_, fs, nperseg=int(seg * fs))
    return f, P
f, P = welch(x, fs)
band = (f >= 1) & (f <= 40); alpha = (f >= 8) & (f < 13)
rel_alpha = float(P[alpha].sum() / P[band].sum())
dominant = float(f[band][np.argmax(P[band])])

with open(os.path.join(art, "real_eeg.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["eeg_uV"]); w.writerows([[float(v)] for v in x])
json.dump(dict(dataset="EEG Motor Movement/Imagery Database (eegmmidb)",
               subject="S001", run="R02 (eyes closed, baseline)", channel="O1",
               fs=fs, n_samples=int(len(x)),
               real_dominant_hz=round(dominant, 2), real_rel_alpha=round(rel_alpha, 3)),
          open(os.path.join(art, "real_meta.json"), "w"), indent=2)
print(f"O1 @ {fs:.0f} Hz, {len(x)} samples: dominant {dominant:.2f} Hz, rel-alpha {rel_alpha:.2f}")
