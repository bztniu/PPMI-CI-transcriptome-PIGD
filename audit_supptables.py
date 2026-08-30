# -*- coding: utf-8 -*-
# SECTION 32 audit v2: xlsx S1-S23 vs pipeline outputs (2026-08-30)
import csv, os
import openpyxl

DIR = r'E:\PPMI帕金森数据库专用\结果\新结果\实验\movementdisorders\new\revision_v2'
XLSX = os.path.join(DIR, '2026-08-27-第二版', 'Supplementary_Tables_updated.xlsx')
T4   = os.path.join(DIR, 'FINAL_PACKAGE_2026-06-21', '04_tables')
FD   = os.path.join(DIR, 'figures_pc1', 'data')
wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)

def rows(sn):
    return [r for r in wb[sn].iter_rows(values_only=True) if any(v is not None for v in r)]

def num(v):
    try: return float(v)
    except (TypeError, ValueError): return None

def rd(path):
    with open(path, encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

R = {}

# S1: rsID|Gene|PPMI log2FC*|... (90 signals)
try:
    src = rd(os.path.join(T4, 'S1_riskgenes_adjusted.csv'))
    de = {d['Symbol']: num(d['log2FoldChange']) for d in rd(os.path.join(FD, 'de_full.csv'))}
    smap = {g: de.get(g) for g in {d['Gene'] for d in src}}
    xs = rows('TableS11')
    hdr_i = next(i for i, r in enumerate(xs) if r[1] == 'Gene')
    data = [r for r in xs[hdr_i+1:] if r[0] and r[1] and not str(r[0]).startswith(('Abbrev','Note'))]
    bad = [(r[1], r[2], smap.get(r[1])) for r in data
           if smap.get(r[1]) is not None and num(r[2]) is not None
           and abs(num(r[2]) - smap[r[1]]) > 5e-4]
    ok = len(data) == 90 and not bad
    R['S11'] = ('PASS' if ok else 'FAIL', f'rows={len(data)} mismatches={len(bad)}')
except Exception as e:
    R['S11'] = ('FAIL', repr(e))

# S2: 63 mapped genes
try:
    src = rd(os.path.join(T4, 'TableS_genepark_CI_gene_mapping.csv'))
    xs = rows('TableS23')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Gene')
    data = [r for r in xs[hdr_i+1:] if r[0] and num(r[1]) is not None]
    R['S23'] = ('PASS' if len(data) == len(src) == 63 else 'FAIL', f'rows={len(data)} src={len(src)}')
except Exception as e:
    R['S23'] = ('FAIL', repr(e))

# S3: 66 genes = ci_expr columns
try:
    with open(os.path.join(FD, 'ci_expr.csv'), encoding='utf-8-sig') as f:
        genes_src = next(csv.reader(f))[1:]
    xs = rows('TableS1')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Gene')
    genes_x = [r[0] for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith('Abbrev')]
    R['S1'] = ('PASS' if sorted(genes_src) == sorted(genes_x) and len(genes_x) == 66 else 'FAIL', f'n={len(genes_x)}')
except Exception as e:
    R['S1'] = ('FAIL', repr(e))

# S4: 393 rows, groups 196/197
try:
    src = rd(os.path.join(FD, 'pc1_scores.csv'))
    xs = rows('TableS2')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'PATNO')
    data = [r for r in xs[hdr_i+1:] if r[0] and r[3] and not str(r[0]).startswith(('Abbrev','Note'))]
    ok = len(data) == 393 and sum(1 for r in data if r[3] == 'Low') == 196 and sum(1 for r in data if r[3] == 'High') == 197
    R['S2'] = ('PASS' if ok else 'FAIL', f'n={len(data)}')
except Exception as e:
    R['S2'] = ('FAIL', repr(e))

# S5/S6/S7/S8 vs rds dump
dump = rd(r'C:\Users\baize\s678_dump.csv')[0]
try:
    xs = rows('TableS3')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    R['S3'] = ('PASS' if '0.506' in txt else 'FAIL', 'PC1 variance 50.6%')
except Exception as e:
    R['S3'] = ('FAIL', repr(e))
try:
    xs = rows('TableS4')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    R['S4'] = ('PASS' if abs(num(next(c for r in xs for c in r if c and str(c).startswith('0.0104'))) - float(dump['dip_D'])) < 1e-6 else 'FAIL', 'dip D')
except Exception as e:
    R['S4'] = ('FAIL', repr(e))
try:
    xs = rows('TableS5')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    R['S5'] = ('PASS' if ('-2505' in txt and '-2512' in txt and abs(float(dump['deltaBIC']) - 6.77) < 0.05) else 'FAIL',
               f'deltaBIC={float(dump["deltaBIC"]):.2f}')
except Exception as e:
    R['S5'] = ('FAIL', repr(e))
try:
    xs = rows('TableS6')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'K')
    data = {str(r[0]): num(r[1]) for r in xs[hdr_i+1:] if r[0] is not None and num(r[1]) is not None}
    exp = {k: float(v) for k, v in {'2': dump['PAC_k2'], '3': dump['PAC_k3'], '4': dump['PAC_k4'], '5': dump['PAC_k5'], '6': dump['PAC_k6']}.items()}
    bad = [k for k in exp if abs(data[k] - exp[k]) > 1e-9]
    R['S6'] = ('PASS' if not bad else 'FAIL', f'PAC={ {k: round(v,4) for k,v in data.items()} } bad={bad}')
except Exception as e:
    R['S6'] = ('FAIL', repr(e))

# S9: 28626 genes, GPNMB log2FC vs de_full
try:
    src = rd(os.path.join(FD, 'de_full.csv'))
    smap = {(d['Symbol'], d['ENSG']): num(d['log2FoldChange']) for d in src}
    xs = rows('TableS7')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Symbol')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note'))]
    gi = 0; li = next(i for i, h in enumerate(xs[hdr_i]) if h == 'log2FC')
    ei = next(i for i, h in enumerate(xs[hdr_i]) if h == 'ENSG ID')
    bad = [(r[0], r[li], smap.get((r[0], r[ei]))) for r in data
           if smap.get((r[0], r[ei])) is not None and num(r[li]) is not None
           and abs(num(r[li]) - smap[(r[0], r[ei])]) > 1e-6][:5]
    ok = len(data) == len(src) == 28626 and not bad
    R['S7'] = ('PASS' if ok else 'FAIL', f'rows={len(data)} GPNMB x={num(next(r[li] for r in data if r[0]=="GPNMB")):.4f} mismatches={len(bad)}')
except Exception as e:
    R['S7'] = ('FAIL', repr(e))

# S10: KEGG ribosome + olfactory
try:
    src = rd(os.path.join(T4, 'GSEA_KEGG_lowCI_ranked_ADJUSTED.csv'))
    smap = {}
    for d in src:
        k = d.get('ID', d.get('Description', ''))
        v = num(d.get('NES'))
        if k and v is not None:
            smap[k] = v
    xs = rows('TableS8')
    hdr_i = next((i for i, r in enumerate(xs) if r[0] == 'ID'), 1)
    nes_idx = next(i for i, h in enumerate(xs[hdr_i]) if h == 'NES')
    bad = []
    ncmp = 0
    for r in xs[hdr_i+1:]:
        if not r[0] or str(r[0]).startswith(('Abbrev','Note')): continue
        nes = num(r[nes_idx])
        for k, v in smap.items():
            if k.lower() in str(r[0]).lower() or str(r[0]).lower() in k.lower():
                ncmp += 1
                if nes is None or abs(nes - v) > 1e-3: bad.append((r[0][:30], nes, v))
                break
    R['S8'] = ('PASS' if ncmp >= 5 and not bad else 'FAIL', f'compared={ncmp} bad={bad[:2]}')
except Exception as e:
    R['S8'] = ('FAIL', repr(e))

# S11: GO top8 names vs ADJUSTED csv (passed before, keep)
try:
    src = rd(os.path.join(T4, 'GO_lowCI_up_ADJUSTED.csv'))
    xs = rows('TableS9')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] in ('ID','Pathway','Description') or (r[0] and 'GO' in str(r[0])[:4]))
    names_x = [str(r[1]) if len(r) > 1 and r[1] else str(r[0]) for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note'))][:8]
    names_s = [d.get('Description', d.get('name', d.get('Term',''))) for d in src[:8]]
    R['S9'] = ('PASS' if [n.lower()[:25] for n in names_x] == [n.lower()[:25] for n in names_s] else 'FAIL',
                f'x0={names_x[0][:22]} s0={names_s[0][:22]}')
except Exception as e:
    R['S9'] = ('FAIL', repr(e))

# S12: key genes present with adjusted stats
try:
    src = rd(os.path.join(T4, 'celladj_keygenes_10genes.csv'))
    adj = {d['Gene']: d['survives_adj'] for d in src}
    xs = rows('TableS12')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    robust_ok = all(g in txt for g in ['ITGA8','RIMS1','DNAH17'])
    attenuated_ok = all(g in txt for g in ['CRHR1','KCNIP3'])
    sig = {g: adj[g] for g in ['ITGA8','RIMS1','DNAH17','CRHR1','KCNIP3'] if g in adj}
    expect = {'ITGA8':'TRUE','RIMS1':'TRUE','DNAH17':'TRUE','CRHR1':'FALSE','KCNIP3':'FALSE'}
    ok_sig = all(str(sig[g]).upper().startswith(expect[g]) for g in sig)
    R['S12'] = ('PASS' if robust_ok and attenuated_ok and ok_sig else 'FAIL', f'sig={sig}')
except Exception as e:
    R['S12'] = ('FAIL', repr(e))

# S13: 20 marker genes, S100A8 log2FC
try:
    src = rd(os.path.join(T4, 'TableS_neutrophil_markers_DE.csv'))
    smap = {(d['Symbol'], d['ENSG']): num(d['log2FoldChange']) for d in src}
    xs = rows('TableS13')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Symbol')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note','baseMean'))]
    li = next(i for i, h in enumerate(xs[hdr_i]) if h == 'log2FC')
    bad = [(r[0], r[li], smap.get(r[0])) for r in data
           if smap.get(r[0]) is not None and num(r[li]) is not None
           and abs(num(r[li]) - smap[r[0]]) > 5e-4]
    ok = len(data) == 20 and not bad
    R['S13'] = ('PASS' if ok else 'FAIL', f'rows={len(data)} S100A8={smap.get("S100A8")}')
except Exception as e:
    R['S13'] = ('FAIL', repr(e))

# S14: 3 robust genes x 2 cohorts
try:
    src = rd(os.path.join(T4, 'Table_risk_genes_comparison.csv'))
    xs = rows('TableS14')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    genes = {d['Gene'] for d in src if d['Track'].startswith('GSE') or 'GENEPARK' in d['Track']}
    ok = all(g in txt for g in ['ITGA8','RIMS1','DNAH17'])
    R['S14'] = ('PASS' if ok else 'FAIL', f'src genes={sorted(genes)}')
except Exception as e:
    R['S14'] = ('FAIL', repr(e))

# S15: PPMI (immune_cohend) + GENEPARK (GSE99039_immune_NNLS)
try:
    psrc = {d['Cell']: num(d['d']) for d in rd(os.path.join(FD, 'immune_cohend.csv'))}
    gsrc = {d['Cell']: num(d['D']) for d in rd(os.path.join(T4, 'GSE99039_immune_NNLS.csv'))}
    xs = rows('TableS15')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Cohort')
    data = [r for r in xs[hdr_i+1:] if r[0] and r[1] and not str(r[0]).startswith(('Abbrev','Note'))]
    ppmi = [r for r in data if r[0] == 'PPMI']; gen = [r for r in data if r[0] == 'GENEPARK']
    bad = []
    for r in ppmi:
        if r[1] in psrc and num(r[2]) is not None and abs(num(r[2]) - psrc[r[1]]) > 1e-6: bad.append(('PPMI', r[1]))
    for r in gen:
        if r[1] in gsrc and num(r[2]) is not None and abs(num(r[2]) - gsrc[r[1]]) > 1e-6: bad.append(('GENEPARK', r[1]))
    ok = len(ppmi) == 16 and len(gen) == 19 and not bad
    R['S15'] = ('PASS' if ok else 'FAIL', f'PPMI={len(ppmi)} GENEPARK={len(gen)} bad={bad[:3]}')
except Exception as e:
    R['S15'] = ('FAIL', repr(e))

# S16: EMM + interaction + quadratic + adjusted row
try:
    emm = rd(os.path.join(T4, 'LMM_PIGD_emmeans_dedup.csv'))
    v12 = {d['group']: num(d['emmean']) for d in emm if d.get('VISIT') == 'V12'}
    xs = rows('TableS16')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    ok_emm = abs(v12['High'] - 1.286) < 0.001 and abs(v12['Low'] - 1.804) < 0.001 and '1.80' in txt and '1.29' in txt
    ok_int = '0.50' in txt and '0.0008' in txt
    ok_quad = '0.047' in txt and '0.034' in txt
    ok_adj = '0.0009' in txt
    R['S16'] = ('PASS' if ok_emm and ok_int and ok_quad and ok_adj else 'FAIL',
                f'emm V12 High={v12["High"]:.3f} Low={v12["Low"]:.3f}')
except Exception as e:
    R['S16'] = ('FAIL', repr(e))

# S17: AIC/BIC vs model comparison csv
try:
    src = rd(os.path.join(FD, 'Table_LMM_model_comparison_dedup.csv'))
    smap = {d.get('Outcome', d.get('Model', '')): d for d in src}
    xs = rows('TableS17')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Outcome')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note','AIC'))]
    bad = []
    for r in data:
        cand = [d for d in src if d['Outcome'] == r[0]]
        if not cand: bad.append((r[0], 'not in src')); continue
        c = cand[0]
        if num(r[3]) is None or abs(num(r[3]) - float(c['Delta_AIC'])) > 0.05: bad.append((r[0], 'dAIC'))
        if num(r[6]) is None or abs(num(r[6]) - float(c['Delta_BIC'])) > 0.05: bad.append((r[0], 'dBIC'))
    R['S17'] = ('PASS' if not bad and len(data) > 0 else 'FAIL', f'rows={len(data)} bad={bad[:3]}')
except Exception as e:
    R['S17'] = ('FAIL', repr(e))

# S18
try:
    xs = rows('TableS18')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    R['S18'] = ('PASS' if ('0.49' in txt and '0.0009' in txt and '0.18' in txt) else 'FAIL', 'cell-adjusted values')
except Exception as e:
    R['S18'] = ('FAIL', repr(e))

# S19: 43 outcomes vs rerun csv
try:
    src = rd(os.path.join(DIR, 'figures_pc1', 'LMM_sweep_rerun_dedup_REML.csv'))
    smap = {d['Outcome']: d for d in src}
    xs = rows('TableS19')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Outcome')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note','LastV'))]
    bad = []
    for r in data:
        s = smap.get(r[0])
        if s is None: bad.append((r[0], 'missing in src')); continue
        if num(r[2]) is not None and abs(num(r[2]) - float(s['V12_int_est'])) > 1e-6: bad.append((r[0], 'est'))
        if num(r[3]) is not None and abs(num(r[3]) - float(s['V12_int_p'])) > 1e-12: bad.append((r[0], 'p'))
    surv = [d['Outcome'] for d in src if float(d['V12_int_padj']) < 0.05]
    R['S19'] = ('PASS' if len(data) == 43 and not bad and surv == ['PIGD'] else 'FAIL',
                f'n={len(data)} survivors={surv} bad={bad[:3]}')
except Exception as e:
    R['S19'] = ('FAIL', repr(e))

# S20: 15 DAT measures
try:
    src = rd(os.path.join(T4, 'DAT_subregion_sweep.csv'))
    smap = {d['Region']: d for d in src}
    xs = rows('TableS20')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Region')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note'))]
    lp_idx = next(i for i, h in enumerate(xs[hdr_i]) if h == 'LastV int P')
    bad = [r[0] for r in data if r[0] in smap and num(r[lp_idx]) is not None
           and abs(num(r[lp_idx]) - float(smap[r[0]]['LastV_int_p'])) > 1e-4]
    R['S20'] = ('PASS' if len(data) == 15 and not bad else 'FAIL', f'n={len(data)} bad={bad[:3]}')
except Exception as e:
    R['S20'] = ('FAIL', repr(e))

# S21: biomarkers
try:
    src = rd(os.path.join(T4, 'biomarker_LMM_sweep.csv'))
    smap = {d['Biomarker']: d for d in src}
    xs = rows('TableS21')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] in ('Biomarker','Metric'))
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith(('Abbrev','Note'))]
    bad = [r[0] for r in data if r[0] in smap and num(r[6]) is not None and num(smap[r[0]]['LastV_int_p']) is not None
           and abs(num(r[6]) - float(smap[r[0]]['LastV_int_p'])) > 1e-9]
    R['S21'] = ('PASS' if not bad and len(data) > 0 else 'FAIL', f'n={len(data)} bad={bad[:3]}')
except Exception as e:
    R['S21'] = ('FAIL', repr(e))

# S22: SAA counts + Fisher 0.844
try:
    src = rd(os.path.join(T4, 'SAA_baseline_status.csv'))[0]
    xs = rows('TableS22')
    txt = ' '.join(str(c) for r in xs for c in r if c is not None)
    # counts: high 187 total/174 positive -> 13 negative; low 185/171 -> 14 negative
    from math import comb
    n = 372; r1 = 187; c1 = 345  # 174+171 positive
    def prob(x): return comb(r1, x) * comb(n - r1, c1 - x) / comb(n, c1)
    p0 = prob(174)
    tot = sum(prob(k) for k in range(max(0, c1 - (n - r1)), min(r1, c1) + 1) if prob(k) <= p0 * (1 + 1e-9))
    ok = abs(tot - 0.8442154) < 1e-5 and abs(float(src['P']) - tot) < 1e-6 and '0.844' in txt
    R['S22'] = ('PASS' if ok else 'FAIL', f'Fisher={tot:.4f} src={src["P"][:8]}')
except Exception as e:
    R['S22'] = ('FAIL', repr(e))

# S23: 14 rows vs final sensitivity csv
try:
    src = rd(os.path.join(FD, 'LMM_PIGD_sensitivity_final.csv'))
    xs = rows('TableS10')
    hdr_i = next(i for i, r in enumerate(xs) if r[0] == 'Model')
    data = [r for r in xs[hdr_i+1:] if r[0] and not str(r[0]).startswith('Factor-time')]
    smap = {d.get('Model', d.get('')): (num(d['beta']), num(d['p'])) for d in src}
    keymap = {'M0':'M0','M0b':'M0b','M1':'M1','M2\u2032':'M2',"M2'":'M2','M5':'M5','M3':'M3','M4':'M4',
              'M6':'M6','M6b':'M6b','M7':'M7',"M2' restricted to OFF":'OFF',"M2' restricted to ON":'ON',
              'M5: RIN-by-visit':'M5_RINxV12','M6b: group-by-LEDD term':'M6b_groupLEDD'}
    bad = []
    for r in data:
        label = r[0]
        key = None
        for k in sorted(keymap, key=len, reverse=True):
            if label.startswith(k):
                key = keymap[k]; break
        if key is None or key not in smap: bad.append((label[:30], 'unmapped')); continue
        beta, pv = smap[key]
        if key == 'M6b' and 'group-by-LEDD term' in label: continue  # different row; skip
        if num(r[2]) is None or abs(num(r[2]) - beta) > 5e-4: bad.append((label[:30], 'beta'))
        pd_ = abs(num(r[4]) - pv) if (num(r[4]) is not None and pv is not None) else 9
        if pd_ > 5e-4 and pd_ / max(abs(pv), 1e-12) > 0.05: bad.append((label[:30], 'p'))
    R['S10'] = ('PASS' if not bad and len(data) >= 12 else 'FAIL', f'rows={len(data)} src={len(src)} bad={bad[:3]}')
except Exception as e:
    R['S10'] = ('FAIL', repr(e))

print('=' * 62)
npass = 0
for k in sorted(R, key=lambda x: int(x[1:])):
    st, det = R[k]
    npass += st == 'PASS'
    print(f'{k:4s} {st}  {det}')
print('=' * 62)
print(f'{npass}/23 tables PASS')
