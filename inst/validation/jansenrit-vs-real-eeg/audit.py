import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
import json, csv
ok = True
def check(l, c):
    global ok; print(f"  [{'PASS' if c else 'FAIL'}] {l}"); ok = ok and bool(c)

case = json.load(open("case.json"))
by   = {cl["id"]: cl for cl in case["claims"]}
s    = {r["metric"]: float(r["value"]) for r in csv.DictReader(open("artifacts/validation_summary.csv"))}
psd  = list(csv.DictReader(open("artifacts/psd.csv")))

# --- claims trace to artifacts ----------------------------------------------
check("real_alpha traces",    abs(s["real_dominant_hz"]             - by["real_alpha"]["value"])   < 0.3)
check("model_default traces", abs(s["model_dominant_default"]       - by["model_default"]["value"]) < 0.3)
check("model_person traces",  abs(s["model_dominant_personalized"]  - by["model_person"]["value"])  < 0.3)
check("real_rel_alpha traces",  abs(s["real_rel_alpha"]  - by["real_rel_alpha"]["value"])  < 5e-3)
check("model_rel_alpha traces", abs(s["model_rel_alpha"] - by["model_rel_alpha"]["value"]) < 5e-3)
check("model_rel_alpha_noaperiodic traces",
      abs(s["model_rel_alpha_noaperiodic"] - by["model_rel_alpha_noaperiodic"]["value"]) < 5e-3)
check("aperiodic_scale traces", abs(s["aperiodic_scale"] - by["aperiodic_scale"]["value"]) < 5e-3)

# --- validation-against-reality: the two requirements -----------------------
check("personalised Jansen-Rit reproduces the real alpha frequency (within 1.5 Hz)",
      abs(s["model_dominant_personalized"] - s["real_dominant_hz"]) < 1.5)
check("personalisation IMPROVES on the generic model",
      abs(s["model_dominant_personalized"] - s["real_dominant_hz"]) <= abs(s["model_dominant_default"] - s["real_dominant_hz"]))
check("alpha is the dominant band in BOTH the real EEG and the model",
      int(s["real_alpha_is_top"]) == 1 and int(s["model_alpha_is_top"]) == 1)
# the model PSD really does peak in the alpha band (8-13 Hz), pointwise on the shared grid
alpha_pts = [p for p in psd if 8 <= float(p["freq"]) < 13]
peak_freq = max(psd, key=lambda p: float(p["model_psd"]))["freq"]
check("model PSD peaks inside the alpha band", 8 <= float(peak_freq) < 13)

# --- spectral shape: the fitted 1/f background CLOSES the relative-alpha gap ------
check("the bare column really did over-concentrate alpha (relative alpha > 0.9)",
      s["model_rel_alpha_noaperiodic"] > 0.9)
check("the 1/f background brings the model's relative alpha within 0.1 of the real EEG",
      abs(s["model_rel_alpha"] - s["real_rel_alpha"]) < 0.1)
check("the 1/f background IMPROVES the spectral-shape match (vs the bare column)",
      abs(s["model_rel_alpha"] - s["real_rel_alpha"]) < abs(s["model_rel_alpha_noaperiodic"] - s["real_rel_alpha"]))

# --- honest scope -----------------------------------------------------------
note = case["validation"]["note"].lower()
check("note distinguishes validation-vs-reality from synthetic verification",
      "verif" in note and "reality" in note)
check("note is honest about single-channel / rhythm-not-waveform, and that the bare column over-concentrated alpha before the 1/f fix",
      "single channel" in note and "waveform" in note and "over-concentrat" in note)
check("note states the over-concentration is now RESOLVED by a 1/f background that is CALIBRATED (fit) to the record",
      "resolv" in note and "1/f" in note and "calibrat" in note)

print("\nRESULT:", "PASS -- mechanistic EEG model validated against a real alpha rhythm" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
