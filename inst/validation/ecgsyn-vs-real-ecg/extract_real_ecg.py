import wfdb, numpy as np, json, os
here = os.path.dirname(os.path.abspath(__file__))
# MIT-BIH record 100 (PhysioNet mitdb, https://physionet.org/content/mitdb/).
# Point MITDB_100 at the record prefix (.../100), or place 100.dat/.hea/.atr under data/.
rec = os.environ.get("MITDB_100", os.path.join(here, "data", "100"))
r = wfdb.rdrecord(rec)
ann = wfdb.rdann(rec, "atr")
fs = r.fs
sig = r.p_signal[:, 0]                      # lead MLII
# beat R-peak samples: keep normal beats 'N'
rpk = np.array([s for s, sym in zip(ann.sample, ann.symbol) if sym == "N"])
rpk = rpk[(rpk > int(0.3*fs)) & (rpk < len(sig)-int(0.5*fs))]
rr = np.diff(rpk) / fs                       # RR intervals (s)
rr = rr[(rr > 0.3) & (rr < 1.5)]             # physiological
hr = 60.0 / np.mean(rr)
sdnn = np.std(rr)
# average beat aligned on R, window -0.25..+0.45 s, resampled to 200 pts
pre, post, N = int(0.25*fs), int(0.45*fs), 200
beats = []
for s in rpk:
    if s-pre >= 0 and s+post < len(sig):
        seg = sig[s-pre:s+post]
        beats.append(np.interp(np.linspace(0, len(seg)-1, N), np.arange(len(seg)), seg))
avg = np.mean(np.array(beats), axis=0)
avg = (avg - avg.min()) / (avg.max() - avg.min())      # normalise 0..1
out = dict(record="MIT-BIH 100 (MLII)", fs=int(fs), n_beats=int(len(rpk)),
           hr=float(hr), sdnn=float(sdnn), rr_mean=float(np.mean(rr)),
           avg_beat=[float(x) for x in avg])
os.makedirs(os.path.join(here, "artifacts"), exist_ok=True)
json.dump(out, open(os.path.join(here, "artifacts", "ecg_real.json"), "w"))
print(f"record 100: n_beats={len(rpk)} HR={hr:.1f} SDNN={sdnn*1000:.0f}ms rr_mean={np.mean(rr):.3f}")
print(f"avg beat: R at idx {np.argmax(avg)}/200, range OK")
