import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))
import json, csv
ok = True
def check(l, c):
    global ok; print(f"  [{'PASS' if c else 'FAIL'}] {l}"); ok = ok and bool(c)

case = json.load(open("case.json"))
by   = {cl["id"]: cl for cl in case["claims"]}
s    = {r["metric"]: float(r["value"]) for r in csv.DictReader(open("artifacts/validation_summary.csv"))}
bt   = list(csv.DictReader(open("artifacts/avg_beat.csv")))

# --- claims trace to artifacts ----------------------------------------------
check("morph_person traces", abs(s["morphology_cor_personalized"] - by["morph_person"]["value"]) < 5e-3)
check("morph_default traces", abs(s["morphology_cor_default"]      - by["morph_default"]["value"]) < 5e-3)
check("hr_match traces",      abs(s["syn_hr"]                      - by["hr_match"]["value"])      < 0.2)
check("sdnn_syn traces",      abs(s["syn_sdnn_ms"]                 - by["sdnn_syn"]["value"])      < 2)
check("n_beats traces",       int(s["real_n_beats"])              == by["n_beats"]["value"])

# --- validation-against-reality: the three requirements ---------------------
check("personalised ECGSYN reproduces the beat morphology (r > 0.85)", s["morphology_cor_personalized"] > 0.85)
check("personalisation IMPROVES on the generic model", s["morphology_cor_personalized"] > s["morphology_cor_default"] + 0.1)
check("synthetic heart rate matches the real (within 2 bpm)", abs(s["syn_hr"] - s["real_hr"]) < 2)
check("synthetic HRV is in the real ballpark (within 15 ms)", abs(s["syn_sdnn_ms"] - s["real_sdnn_ms"]) < 15)
# the personalised template really is closer to the real beat, pointwise
import statistics
def rms(pred): return (sum((float(b["real"]) - float(b[pred]))**2 for b in bt) / len(bt))**0.5
check("personalised template has lower error than default (pointwise)", rms("ecgsyn_personalized") < rms("ecgsyn_default"))

# --- honest scope -----------------------------------------------------------
check("note distinguishes validation-vs-reality from synthetic verification",
      "verif" in case["validation"]["note"].lower() and "reality" in case["validation"]["note"].lower())
check("note is honest about single-lead / by-construction HRV",
      "single lead" in case["validation"]["note"].lower() and "construction" in case["validation"]["note"].lower())

print("\nRESULT:", "PASS -- mechanistic ECG model validated against real data" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
