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
check("real_sd traces",     abs(s["real_sd_ms"]      - by["real_sd"]["value"])     < 0.2)
check("model_sd traces",    abs(s["model_sd_ms"]     - by["model_sd"]["value"])    < 0.5)
check("real_lf_hf traces",  abs(s["real_lf_hf"]      - by["real_lf_hf"]["value"])  < 0.05)
check("model_lf_hf traces", abs(s["model_lf_hf"]     - by["model_lf_hf"]["value"]) < 0.2)
check("lf_peak traces",     abs(s["real_lf_peak_hz"] - by["lf_peak"]["value"])     < 5e-3)
check("hf_peak traces",     abs(s["real_hf_peak_hz"] - by["hf_peak"]["value"])     < 5e-3)

# --- validation-against-reality: the requirements ---------------------------
check("the fitted model reproduces the real RR SD (within 15%)",
      abs(s["model_sd_ms"] - s["real_sd_ms"]) < 0.15 * s["real_sd_ms"])
check("the fitted model reproduces the real LF/HF autonomic balance (within 20%)",
      abs(s["model_lf_hf"] - s["real_lf_hf"]) < 0.20 * s["real_lf_hf"])
check("the real record has a separated Mayer (LF) and respiratory (HF) rhythm",
      s["real_hf_peak_hz"] > s["real_lf_peak_hz"] + 0.03 and s["real_lf_peak_hz"] < 0.15)
check("the fit is stable and physiological (interior damping, not railed)",
      0.30 < s["fit_rho"] < 0.995)
# the model PSD carries power in both the LF and HF bands (not a single-band collapse)
lf = sum(float(p["model_psd"]) for p in psd if 0.04 <= float(p["freq"]) < 0.15)
hf = sum(float(p["model_psd"]) for p in psd if 0.15 <= float(p["freq"]) < 0.40)
check("the fitted-model spectrum has power in both LF and HF bands", lf > 0 and hf > 0)

# --- honest scope -----------------------------------------------------------
note = case["validation"]["note"].lower()
check("note distinguishes validation-vs-reality from synthetic verification",
      "verif" in note and "reality" in note)
check("note is honest about second-order/spectral, linear idealisation, free amplitudes, single record",
      ("second-order" in note or "spectral" in note) and "linear" in note
      and "free amplitud" in note and "single record" in note)

print("\nRESULT:", "PASS -- closed-loop cardiovascular model validated against real HRV" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
