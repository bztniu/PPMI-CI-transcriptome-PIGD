# -*- coding: utf-8 -*-
# Verify every row of Table 1 (Table1.xlsx) from raw PPMI files.
# Grouping: figures_pc1/data/pc1_scores.csv (trusted PC1 median split).
import pandas as pd, numpy as np, glob, os
from scipy.stats import mannwhitneyu, fisher_exact

ROOT = r"E:/PPMI帕金森数据库专用"
FD   = ROOT + r"/结果/新结果/实验/movementdisorders/new/revision_v2/figures_pc1/data"

def rd(p, **kw):
    for enc in ("utf-8-sig", "gbk", "latin-1"):
        try:
            return pd.read_csv(p, encoding=enc, low_memory=False, **kw)
        except UnicodeDecodeError:
            continue
    raise IOError(p)

pc = rd(FD + "/pc1_scores.csv")
pc["group"] = pc["group"].astype(str)
assert pc["group"].value_counts().to_dict() == {"Low": 196, "High": 197}, pc["group"].value_counts()
print(f"N={len(pc)} (Low={sum(pc.group=='Low')} High={sum(pc.group=='High')})")

def bl(df):
    d = df[df["EVENT_ID"].isin(["BL", "SC"])]
    return d.drop_duplicates("PATNO")

def num(s):
    return pd.to_numeric(s, errors="coerce")

# Sex
demo = rd(ROOT + "/Demographics_03Feb2026.csv")
sx = demo.set_index("PATNO")["SEX"].reindex(pc["PATNO"])
print("SEX values:", sx.dropna().unique()[:6])
pc["MALE"] = sx.isin(["M", "Male", 1, "1"]).astype(float).where(sx.notna()).to_numpy()

# Age (BL age at visit)
age = rd(ROOT + "/Age_at_visit_17Mar2025.csv")
age = age[age["EVENT_ID"] == "BL"].drop_duplicates("PATNO")
pc["AGE"] = num(age.set_index("PATNO")["AGE_AT_VISIT"].reindex(pc["PATNO"])).to_numpy()

# Education
edu = rd(ROOT + "/ppmi数据表/Subject_Characteristics/Socio-Economics_05Apr2026.csv")
edu = edu[edu["EVENT_ID"].isin(["SC", "BL"])].drop_duplicates("PATNO")
pc["EDUC"] = num(edu.set_index("PATNO")["EDUCYRS"].reindex(pc["PATNO"])).to_numpy()

# Part III (BL/SC dedup)
p3 = bl(rd(ROOT + "/运动症状数据/MDS-UPDRS_Part_III_29Jan2026.csv")).set_index("PATNO")
pc["NP3TOT"] = num(p3["NP3TOT"].reindex(pc["PATNO"])).to_numpy()
pc["NHY"]    = num(p3["NHY"].reindex(pc["PATNO"])).to_numpy()
sub = p3.reindex(pc["PATNO"]).reset_index(drop=True)
pc["PIGD"] = num(sub["NP3GAIT"]).fillna(0) + num(sub["NP3PSTBL"]).fillna(0) + num(sub["NP3FRZGT"]).fillna(0)
pc["TotTrem"] = sum(num(sub[c]).fillna(0) for c in
    ["NP3PTRMR","NP3PTRML","NP3KTRMR","NP3KTRML","NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL"])
pc["Rest4"] = sum(num(sub[c]).fillna(0) for c in ["NP3RTARU","NP3RTALU","NP3RTARL","NP3RTALL"])
pc["Rest5"] = pc["Rest4"] + num(sub["NP3RTALJ"]).fillna(0)

# Disease duration (BL exam date - PDDXDT)/12, v2 construction
dx = bl(rd(ROOT + "/CI_影像学分析/data/longitudinal_imaging/raw_inputs/PD_Diagnosis_History_05Apr2026.csv")).set_index("PATNO")
pddx = dx["PDDXDT"].reindex(pc["PATNO"]).reset_index(drop=True)
def my(s):
    if not isinstance(s, str) or "/" not in s: return np.nan
    a = s.split("/")
    if len(a) != 2: return np.nan
    try: return int(a[0]) + 12*int(a[1])
    except: return np.nan
b = exam = sub["EXAMDT"].map(my); d0 = pddx.map(my)
pc["DUR"] = np.where((b.notna()) & (d0.notna()) & (b >= d0), (b - d0)/12.0, np.nan)

# UPDRS I / II
p1 = bl(rd(ROOT + "/运动症状数据/MDS-UPDRS_Part_I_31Jan2026.csv")).set_index("PATNO")
pc["NP1RTOT"] = num(p1["NP1RTOT"].reindex(pc["PATNO"])).to_numpy()
p2 = bl(rd(ROOT + "/运动症状数据/MDS_UPDRS_Part_II__Patient_Questionnaire_31Jan2026.csv")).set_index("PATNO")
pc["NP2PTOT"] = num(p2["NP2PTOT"].reindex(pc["PATNO"])).to_numpy()

# MoCA / RBDSQ / UPSIT
moca = bl(rd(ROOT + "/运动症状数据/Montreal_Cognitive_Assessment__MoCA__31Jan2026.csv")).set_index("PATNO")
pc["MoCA"] = num(moca["MCATOT"].reindex(pc["PATNO"])).to_numpy()
rb = bl(rd(ROOT + "/ppmi数据表/Non-motor_Assessments/REM_Sleep_Behavior_Disorder_Screening_Questionnaire_05Apr2026.csv")).set_index("PATNO").reindex(pc["PATNO"]).reset_index(drop=True)
items = ["PTCGBOTH","DRMVIVID","DRMAGRAC","DRMNOCTB","SLPLMBMV","SLPINJUR","DRMVERBL","DRMFIGHT","DRMUMV","DRMOBJFL"]
rbi = pd.concat([num(rb[c]) for c in items], axis=1)
pc["RBDSQ"] = rbi.sum(axis=1).where(rbi.notna().all(axis=1))
ups = bl(rd(ROOT + "/ppmi数据表/Non-motor_Assessments/University_of_Pennsylvania_Smell_Identification_Test_UPSIT_05Apr2026.csv")).set_index("PATNO").reindex(pc["PATNO"]).reset_index(drop=True)
ucols = [c for c in ups.columns if c.startswith("SCENT_") and c.endswith("_CORRECT")]
print("UPSIT correct columns:", len(ucols))
ui = pd.concat([num(ups[c]) for c in ucols], axis=1)
pc["UPSIT"] = ui.sum(axis=1).where(ui.notna().all(axis=1))

# SAA
saa = rd(ROOT + "/运动症状数据/SAA_Biospecimen_Analysis_Results_23Feb2026.csv")
saa = saa[(saa["CLINICAL_EVENT"]=="BL") & (saa["TYPE"]=="Cerebrospinal Fluid") & (saa["COHORT"]=="PD")].drop_duplicates("PATNO")
pc["SAA"] = saa.set_index("PATNO")["SAA_Status"].reindex(pc["PATNO"]).to_numpy()

# RIN
meta = rd(ROOT + "/metaDataIR3.csv")
mm = meta.set_index("Specimen Bar Code")["RIN Value"]
pc["RIN"] = num(mm.reindex(pc["SAMPLE_ID"]).values)

def cont(nm, k):
    lo = pc.loc[pc.group=="Low", k].dropna(); hi = pc.loc[pc.group=="High", k].dropna()
    p = mannwhitneyu(lo, hi, alternative="two-sided").pvalue
    print(f"{nm:<28}| {lo.mean():.2f} ({lo.std(ddof=1):.2f}) n={len(lo):<3}| {hi.mean():.2f} ({hi.std(ddof=1):.2f}) n={len(hi):<3}| P={p:.3f}")

print("\n--- continuous rows (Wilcoxon rank-sum) ---")
for nm, k in [("Age, years","AGE"),("Education, years","EDUC"),("Symptom duration, years","DUR"),
              ("Hoehn-Yahr stage","NHY"),("MDS-UPDRS II total","NP2PTOT"),("MDS-UPDRS III total","NP3TOT"),
              ("PIGD score","PIGD"),("Total tremor score","TotTrem"),("Rest tremor (4-limb)","Rest4"),
              ("Rest tremor (5-item+jaw)","Rest5"),("MDS-UPDRS I total","NP1RTOT"),("MoCA total","MoCA"),
              ("RBDSQ total","RBDSQ"),("UPSIT total","UPSIT"),("RIN","RIN")]:
    cont(nm, k)

print("\n--- Sex male ---")
lo = pc.loc[pc.group=="Low","MALE"].dropna(); hi = pc.loc[pc.group=="High","MALE"].dropna()
tab = [[int(lo.sum()), int((1-lo).sum())],[int(hi.sum()), int((1-hi).sum())]]
print(f"Low: {int(lo.sum())} ({100*lo.mean():.1f}%) of {len(lo)} | High: {int(hi.sum())} ({100*hi.mean():.1f}%) of {len(hi)} | Fisher P={fisher_exact(tab)[1]:.3f}")

print("\n--- CSF SAA ---")
lo = pc.loc[pc.group=="Low","SAA"].dropna(); hi = pc.loc[pc.group=="High","SAA"].dropna()
print("Low :", lo.value_counts().to_dict())
print("High:", hi.value_counts().to_dict())
l2 = lo[lo.isin(["Positive","Negative"])]; h2 = hi[hi.isin(["Positive","Negative"])]
t22 = [[int((l2=="Positive").sum()), int((l2=="Negative").sum())],
       [int((h2=="Positive").sum()), int((h2=="Negative").sum())]]
odds, p22 = fisher_exact(t22)
print(f"2x2 Pos-vs-Neg (inconclusive excluded): Low {t22[0][0]}/{sum(t22[0])} ({100*t22[0][0]/sum(t22[0]):.1f}%) | High {t22[1][0]}/{sum(t22[1])} ({100*t22[1][0]/sum(t22[1]):.1f}%) | Fisher P={p22:.3f}")
