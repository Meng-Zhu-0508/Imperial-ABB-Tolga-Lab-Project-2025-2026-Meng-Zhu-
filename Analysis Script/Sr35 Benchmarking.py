#!/usr/bin/env python3
"""
af3_benchmark.py  ——  AF3 模型 vs 实验结构 的一站式打分脚本
 
计算: 全局/分链 Cα RMSD、TM-score、DockQ + fnat/fnonnat/iRMSD/LRMSD (每个界面 + 平均)
      pTM/ipTM (优先读 summary_confidences.json；没有则从 full_data 的 PAE 估算)
      pLDDT (per-residue = Cα 原子)
输出: Markdown 表 + CSV + 每个 job 的详细 txt
 
用法:
    python af3_benchmark.py jobs.json -o results
jobs.json 格式见脚本末尾 EXAMPLE_CONFIG。
 
依赖: numpy, scipy      (pip install numpy scipy)
"""
import sys, os, json, argparse, itertools, collections
import numpy as np
from scipy.spatial import cKDTree
 
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}
BACKBONE = ('N','CA','C','O')
 
# ---------------------------------------------------------------- mmCIF 解析
def load_cif(path):
    """返回 res[chain][resnum] = {atom: xyz}, comp[(chain,resnum)], plddt[chain][resnum] (Cα B-factor)"""
    cols, rows, inloop = [], [], False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('_atom_site.'):
                cols.append(line.strip().split('.')[1]); inloop = True; continue
            if inloop:
                if line.startswith(('ATOM', 'HETATM')):
                    rows.append(line.split())
                elif line.startswith('#') and rows:
                    break
    if not rows:
        raise SystemExit(f'{path}: 没找到 _atom_site 记录')
    I = {c: i for i, c in enumerate(cols)}
    res = collections.defaultdict(dict); comp = {}; bf = collections.defaultdict(dict)
    het = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        if r[I['type_symbol']] == 'H':
            continue
        cm = r[I['label_comp_id']]
        ch = r[I['auth_asym_id']]
        xyz = np.array([float(r[I['Cartn_x']]), float(r[I['Cartn_y']]), float(r[I['Cartn_z']])])
        if cm not in AA3:
            if cm != 'HOH':
                het[ch][cm].append(xyz)
            continue
        n = int(r[I['auth_seq_id']])
        res[ch].setdefault(n, {})[r[I['label_atom_id']]] = xyz
        comp[(ch, n)] = cm
        if r[I['label_atom_id']] == 'CA':
            bf[ch][n] = float(r[I['B_iso_or_equiv']])
    return dict(res), comp, dict(bf), {k: dict(v) for k, v in het.items()}
 
def seq_of(res, comp, ch):
    ks = sorted(res[ch]); return ks, ''.join(AA3[comp[(ch, k)]] for k in ks)
 
# ---------------------------------------------------------------- 编号对齐
def find_offset(mres, mcomp, mch, nres, ncomp, nch, span=1500):
    """暴力扫整数偏移，使 model 残基号 + off = native 残基号 的序列一致度最高"""
    nmap = {k: AA3[ncomp[(nch, k)]] for k in nres[nch]}
    best = (0, -1, 0)
    for off in range(-span, span + 1):
        hit = tot = 0
        for k in mres[mch]:
            t = k + off
            if t in nmap:
                tot += 1
                if nmap[t] == AA3[mcomp[(mch, k)]]:
                    hit += 1
        if tot >= 30 and hit > best[1]:
            best = (off, hit, tot)
    return best  # (offset, identical, overlap)
 
# ---------------------------------------------------------------- 几何
def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    V, S, Wt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(V @ Wt))
    return V @ np.diag([1, 1, d]) @ Wt, pc, qc
 
def transform(X, R, pc, qc): return (X - pc) @ R + qc
def rmsd(P, Q): return float(np.sqrt(((P - Q) ** 2).sum(1).mean()))
 
def tm_score(P, Q, Lnorm):
    """标准 TM-score：多起点片段 + 迭代扩展，取最大值"""
    d0 = max(1.24 * (Lnorm - 15) ** (1 / 3) - 1.8, 0.5)
    L = len(P); best = 0.0
    for frag in [L, L // 2, L // 4, L // 8, 4]:
        frag = max(int(frag), 4)
        for start in range(0, L - frag + 1, max(1, frag // 2)):
            idx = np.arange(start, start + frag)
            for _ in range(30):
                R, pc, qc = kabsch(P[idx], Q[idx])
                d = np.linalg.norm(transform(P, R, pc, qc) - Q, axis=1)
                score = float((1 / (1 + (d / d0) ** 2)).sum() / Lnorm)
                best = max(best, score)
                cut = d0
                new = np.where(d < cut)[0]
                while len(new) < 4:
                    cut += 0.5; new = np.where(d < cut)[0]
                if len(new) == len(idx) and np.all(new == idx):
                    break
                idx = new
    return best
 
# ---------------------------------------------------------------- DockQ
def residue_contacts(res, c1, c2, cutoff=5.0):
    A = [(n, x) for n, at in res[c1].items() for x in at.values()]
    B = [(n, x) for n, at in res[c2].items() for x in at.values()]
    if not A or not B: return set()
    ta = cKDTree(np.array([p[1] for p in A])); tb = cKDTree(np.array([p[1] for p in B]))
    out = set()
    for i, js in enumerate(ta.query_ball_tree(tb, cutoff)):
        for j in js:
            out.add((A[i][0], B[j][0]))
    return out
 
def interface_residues(res, c1, c2, cutoff=10.0):
    A = [(n, x) for n, at in res[c1].items() for x in at.values()]
    B = [(n, x) for n, at in res[c2].items() for x in at.values()]
    ta = cKDTree(np.array([p[1] for p in A])); tb = cKDTree(np.array([p[1] for p in B]))
    s1, s2 = set(), set()
    for i, js in enumerate(ta.query_ball_tree(tb, cutoff)):
        if js: s1.add(A[i][0])
        for j in js: s2.add(B[j][0])
    return s1, s2
 
def paired_backbone(mres, nres, mch, nch, resnums):
    P, Q = [], []
    for n in sorted(resnums):
        for a in BACKBONE:
            if a in mres[mch].get(n, {}) and a in nres[nch].get(n, {}):
                P.append(mres[mch][n][a]); Q.append(nres[nch][n][a])
    return np.array(P), np.array(Q)
 
def dockq_formula(fnat, irms, lrms):
    return (fnat + 1 / (1 + (irms / 1.5) ** 2) + 1 / (1 + (lrms / 8.5) ** 2)) / 3
 
def capri_class(dq):
    return 'High' if dq >= 0.80 else 'Medium' if dq >= 0.49 else 'Acceptable' if dq >= 0.23 else 'Incorrect'
 
def dockq_pair(model, native, m1, n1, m2, n2, fnat_cut=5.0, iface_cut=10.0):
    com1 = set(model[m1]) & set(native[n1]); com2 = set(model[m2]) & set(native[n2])
    NC = {(a, b) for a, b in residue_contacts(native, n1, n2, fnat_cut) if a in com1 and b in com2}
    MC = {(a, b) for a, b in residue_contacts(model, m1, m2, fnat_cut) if a in com1 and b in com2}
    if not NC: return None
    fnat = len(NC & MC) / len(NC)
    fnonnat = len(MC - NC) / len(MC) if MC else 0.0
    i1, i2 = interface_residues(native, n1, n2, iface_cut)
    P1, Q1 = paired_backbone(model, native, m1, n1, i1 & com1)
    P2, Q2 = paired_backbone(model, native, m2, n2, i2 & com2)
    P = np.vstack([P1, P2]); Q = np.vstack([Q1, Q2])
    R, pc, qc = kabsch(P, Q); irms = rmsd(transform(P, R, pc, qc), Q)
    if len(com1) >= len(com2): rm, rn, lm, ln, rc, lc = m1, n1, m2, n2, com1, com2
    else:                      rm, rn, lm, ln, rc, lc = m2, n2, m1, n1, com2, com1
    Pr, Qr = paired_backbone(model, native, rm, rn, rc)
    Pl, Ql = paired_backbone(model, native, lm, ln, lc)
    R, pc, qc = kabsch(Pr, Qr); lrms = rmsd(transform(Pl, R, pc, qc), Ql)
    dq = dockq_formula(fnat, irms, lrms)
    return dict(DockQ=dq, fnat=fnat, fnonnat=fnonnat, iRMSD=irms, LRMSD=lrms,
                nat_contacts=len(NC), model_contacts=len(MC), CAPRI=capri_class(dq))
 
# ---------------------------------------------------------------- 链映射
def pair_offsets(mres, mcomp, mchains, nres, ncomp, nchains):
    """对每个 (model chain, native chain) 求最佳整数偏移与序列一致度"""
    cand = {}
    for mc in mchains:
        cand[mc] = {}
        for nc in nchains:
            off, hit, tot = find_offset(mres, mcomp, mc, nres, ncomp, nc)
            cand[mc][nc] = (off, hit, tot, hit / tot if tot else 0.0)
    return cand
 
def best_chain_mapping(mres, mchains, nres, nchains, cand, min_ident=0.90, min_overlap=30):
    """只在序列同源的链之间穷举映射，取全局 Cα RMSD 最小者"""
    allowed = {mc: [nc for nc in nchains
                    if cand[mc][nc][2] >= min_overlap and cand[mc][nc][3] >= min_ident]
               for mc in mchains}
    for mc in mchains:
        if not allowed[mc]:                      # 放宽
            allowed[mc] = sorted(nchains, key=lambda nc: -cand[mc][nc][3])[:1]
    def rms_for(mapping):
        P, Q = [], []
        for mc in mchains:
            nc = mapping[mc]; off = cand[mc][nc][0]
            ks = [k + off for k in mres[mc] if 'CA' in mres[mc][k]]
            ks = [k for k in ks if k in nres[nc] and 'CA' in nres[nc][k]]
            if len(ks) < min_overlap: return None, 0
            P += [mres[mc][k - off]['CA'] for k in ks]
            Q += [nres[nc][k]['CA'] for k in ks]
        P, Q = np.array(P), np.array(Q)
        R, pc, qc = kabsch(P, Q)
        return rmsd(transform(P, R, pc, qc), Q), len(P)
    best = (None, 1e9, 0)
    combos = itertools.product(*[allowed[mc] for mc in mchains])
    for i, combo in enumerate(combos):
        if i > 5040: break                        # 安全上限
        if len(set(combo)) != len(combo): continue
        mp = dict(zip(mchains, combo))
        r, n = rms_for(mp)
        if r is not None and r < best[1]: best = (mp, r, n)
    return best
 
# ---------------------------------------------------------------- 置信度
def read_confidences(job):
    """优先 summary_confidences.json；否则从 full_data 的 PAE 估算 pTM/ipTM"""
    out = {'pTM': None, 'ipTM': None, 'source': '—'}
    if job.get('ptm') is not None or job.get('iptm') is not None:
        out['pTM'] = job.get('ptm'); out['ipTM'] = job.get('iptm')
        out['source'] = '配置文件手填'
        return out
    p = job.get('summary_confidences')
    if p and os.path.exists(p):
        d = json.load(open(p))
        for k in ('ptm', 'pTM'):
            if k in d: out['pTM'] = float(d[k])
        for k in ('iptm', 'ipTM'):
            if k in d: out['ipTM'] = float(d[k])
        out['source'] = 'summary_confidences.json'
        return out
    p = job.get('full_data')
    if p and os.path.exists(p):
        d = json.load(open(p))
        if 'pae' in d and 'token_chain_ids' in d:
            pae = np.asarray(d['pae'], float); ch = np.asarray(d['token_chain_ids'])
            N = len(ch); d0 = max(1.24 * (N - 15) ** (1 / 3) - 1.8, 0.5)
            f = 1 / (1 + (pae / d0) ** 2)
            out['pTM'] = float(f.mean(1).max())
            inter = ch[None, :] != ch[:, None]
            if inter.any():
                rowsum = (f * inter).sum(1); rowcnt = inter.sum(1)
                with np.errstate(invalid='ignore'):
                    out['ipTM'] = float(np.nanmax(np.where(rowcnt > 0, rowsum / np.maximum(rowcnt, 1), np.nan)))
            out['source'] = 'PAE 估算 (近似，非官方值)'
    return out
 
# ---------------------------------------------------------------- 主流程
def run_job(job):
    name = job['name']
    nres, ncomp, _, nhet = load_cif(job['native'])
    mres, mcomp, mbf, mhet = load_cif(job['model'])
    nchains = job.get('native_chains') or sorted(nres)
    mchains = job.get('model_chains') or sorted(mres)
    # 1) 逐 (model, native) 链对求编号偏移
    cand = pair_offsets(mres, mcomp, mchains, nres, ncomp, nchains)
    # 2) 链映射（只在同源链之间）
    mapping, gRMSD, nCA = best_chain_mapping(mres, mchains, nres, nchains, cand)
    offsets = {mc: cand[mc][mapping[mc]][0] for mc in mchains}
    detail = [f'  model {mc} -> native {mapping[mc]}: 偏移 {offsets[mc]:+d}  '
              f'(序列一致 {cand[mc][mapping[mc]][1]}/{cand[mc][mapping[mc]][2]})' for mc in mchains]
    M  = {mc: {k + offsets[mc]: v for k, v in mres[mc].items()} for mc in mchains}
    MB = {mc: {k + offsets[mc]: v for k, v in mbf[mc].items()} for mc in mchains if mc in mbf}
    # 3) 分链 RMSD + TM-score
    P, Q, perchain = [], [], []
    for mc in mchains:
        nc = mapping[mc]
        ks = sorted(set(M[mc]) & set(nres[nc]))
        ks = [k for k in ks if 'CA' in M[mc][k] and 'CA' in nres[nc][k]]
        p = np.array([M[mc][k]['CA'] for k in ks]); q = np.array([nres[nc][k]['CA'] for k in ks])
        R, pc, qc = kabsch(p, q)
        perchain.append((mc, nc, len(ks), rmsd(transform(p, R, pc, qc), q),
                         tm_score(p, q, len(nres[nc]))))
        P.append(p); Q.append(q)
    Pall, Qall = np.vstack(P), np.vstack(Q)
    R, pc, qc = kabsch(Pall, Qall)
    global_rmsd = rmsd(transform(Pall, R, pc, qc), Qall)
    # 4) DockQ（所有有接触的 native 链对）
    inv = {v: k for k, v in mapping.items()}
    ifaces = []
    for a, b in itertools.combinations([mapping[mc] for mc in mchains], 2):
        if len(residue_contacts(nres, a, b, 5.0)) < 10: continue
        r = dockq_pair(M, nres, inv[a], a, inv[b], b)
        if r: r['pair'] = f'{inv[a]}{inv[b]} -> {a}{b}'; ifaces.append(r)
    mean = {}
    if ifaces:
        for k in ('DockQ', 'fnat', 'fnonnat', 'iRMSD', 'LRMSD'):
            mean[k] = float(np.mean([x[k] for x in ifaces]))
        mean['CAPRI'] = capri_class(mean['DockQ'])
    # 5) pLDDT
    allp = np.array([v for mc in MB for v in MB[mc].values()])
    pl = dict(mean=float(allp.mean()), median=float(np.median(allp)), max=float(allp.max()),
              b90=float(100 * (allp >= 90).mean()), b70=float(100 * ((allp >= 70) & (allp < 90)).mean()),
              b50=float(100 * ((allp >= 50) & (allp < 70)).mean()), b0=float(100 * (allp < 50).mean()))
    conf = read_confidences(job)
    return dict(name=name, native_label=job.get('native_label'), mapping=mapping, offsets=offsets, offset_detail=detail,
                global_rmsd=global_rmsd, n_ca=len(Pall), perchain=perchain,
                interfaces=ifaces, mean=mean, plddt=pl, conf=conf,
                native=job['native'], model=job['model'],
                native_het=nhet, model_het=mhet)
 
def fmt(x, n=3):
    return '—' if x is None else (f'{x:.{n}f}' if isinstance(x, float) else str(x))
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('config', help='jobs.json')
    ap.add_argument('-o', '--out', default='af3_benchmark_out')
    args = ap.parse_args()
    jobs = json.load(open(args.config))
    os.makedirs(args.out, exist_ok=True)
    results = [run_job(j) for j in jobs]
 
    hdr = ['Cryo-EM', 'AF3 Model', 'pTM', 'ipTM', 'RMSD (Å)', '分链 RMSD (Å)', 'TM-score',
           'DockQ', 'CAPRI', 'fnat', 'LRMSD (Å)', 'iRMSD (Å)', 'mean pLDDT', 'nCA']
    rows = []
    for r in results:
        m = r['mean']
        per_r = ' / '.join(f'{p[3]:.2f}' for p in r['perchain'])
        per_t = ' / '.join(f'{p[4]:.3f}' for p in r['perchain'])
        if len(r['perchain']) > 3:                      # 同源多聚体，收敛成范围
            rr = [p[3] for p in r['perchain']]; tt = [p[4] for p in r['perchain']]
            per_r = f'{min(rr):.2f}–{max(rr):.2f}'; per_t = f'{min(tt):.3f}–{max(tt):.3f}'
        rows.append([
            r.get('native_label') or os.path.splitext(os.path.basename(r['native']))[0],
            r['name'], fmt(r['conf']['pTM'], 2), fmt(r['conf']['ipTM'], 2),
            f"{r['global_rmsd']:.2f}", per_r, per_t,
            fmt(m.get('DockQ')), m.get('CAPRI', '—'), fmt(m.get('fnat')),
            fmt(m.get('LRMSD'), 2), fmt(m.get('iRMSD'), 2),
            f"{r['plddt']['mean']:.1f}", str(r['n_ca'])])
    md = ['| ' + ' | '.join(hdr) + ' |', '|' + '|'.join(['---'] * len(hdr)) + '|']
    md += ['| ' + ' | '.join(row) + ' |' for row in rows]
    md = '\n'.join(md)
    foot = ('\n\n注:\n'
            '- RMSD = 全部指定链一起做刚体叠合后的 Cα RMSD；分链 RMSD = 每条链各自单独叠合。\n'
            '  异源复合物两条链误差差一个数量级时，全局 RMSD 没有意义，请看分链值。\n'
            '- TM-score 归一到 native 链长，>0.5 才算同一折叠。\n'
            '- DockQ/fnat/iRMSD/LRMSD 为所有有接触的链对的算术平均（与 DockQ v2 的\n'
            '  "Total DockQ over N native interfaces" 定义一致）；主链原子 = N/CA/C/O，\n'
            '  接触阈值 5.0 Å，界面阈值 10.0 Å。\n'
            '- pLDDT 为 per-residue（Cα 原子）。AF3 的 mmCIF 是逐原子 pLDDT，取 Cα 才是惯例值。\n'
            '- pTM/ipTM 优先读 summary_confidences.json；没有时从 PAE 估算，为近似值，论文请用官方值。\n')
    open(os.path.join(args.out, 'summary_table.md'), 'w').write(md + foot)
    with open(os.path.join(args.out, 'summary_table.csv'), 'w') as fh:
        fh.write(','.join(hdr) + '\n')
        for row in rows: fh.write(','.join(f'"{c}"' for c in row) + '\n')
 
    with open(os.path.join(args.out, 'details.txt'), 'w') as fh:
        for r in results:
            fh.write('=' * 72 + f"\n{r['name']}\n  native: {r['native']}\n  model : {r['model']}\n")
            fh.write('  编号偏移:\n' + '\n'.join(r['offset_detail']) + '\n')
            fh.write(f"  链映射 (model->native): {r['mapping']}\n")
            fh.write(f"  全局 Cα RMSD: {r['global_rmsd']:.3f} Å over {r['n_ca']} CA\n")
            for mc, nc, n, rr, tt in r['perchain']:
                fh.write(f'    {mc}->{nc}: n={n:5d}  单链叠合 RMSD {rr:6.2f} Å   TM-score {tt:.3f}\n')
            if r['interfaces']:
                fh.write('  DockQ 逐界面:\n')
                for i in r['interfaces']:
                    fh.write(f"    {i['pair']:<16s} DockQ {i['DockQ']:.3f} ({i['CAPRI']:<10s}) "
                             f"fnat {i['fnat']:.3f} fnonnat {i['fnonnat']:.3f} "
                             f"iRMSD {i['iRMSD']:6.2f} LRMSD {i['LRMSD']:6.2f} "
                             f"[nat {i['nat_contacts']} / model {i['model_contacts']} 接触对]\n")
                fh.write(f"    平均: DockQ {r['mean']['DockQ']:.3f} ({r['mean']['CAPRI']})\n")
            else:
                fh.write('  DockQ: 无界面（单链 job）\n')
            p = r['plddt']
            fh.write(f"  pLDDT(Cα): mean {p['mean']:.1f} median {p['median']:.1f} max {p['max']:.1f} | "
                     f">90 {p['b90']:.1f}%  70-90 {p['b70']:.1f}%  50-70 {p['b50']:.1f}%  <50 {p['b0']:.1f}%\n")
            fh.write(f"  pTM/ipTM 来源: {r['conf']['source']}\n")
            if r['native_het']: fh.write(f"  native 配体: {r['native_het'] and {k:list(v) for k,v in r['native_het'].items()}}\n")
            if r['model_het']:  fh.write(f"  model  配体: {r['model_het'] and {k:list(v) for k,v in r['model_het'].items()}}\n")
            fh.write('\n')
    print(md); print(f'\n-> {args.out}/summary_table.md, summary_table.csv, details.txt')
 
EXAMPLE_CONFIG = """
[
 {"name": "Sr35 Pentamer", "native_label": "7XE0",
  "native": "7XE0 CryoEM.cif", "native_chains": ["A","C","E","G","I"],
  "model": "fold_7xe0_..._model_0.cif", "model_chains": ["A","B","C","D","E"],
  "summary_confidences": "fold_7xe0_..._summary_confidences_0.json"},
 {"name": "Sr35 + AvrSr35", "native_label": "7XVG",
  "native": "7XVG CyroEM.cif", "model": "fold_7xvg_af3_model_0.cif",
  "summary_confidences": "fold_7xvg_af3_summary_confidences_0.json",
  "full_data": "fold_7xvg_af3_full_data_0.json"},
 {"name": "AvrSr35 Alone", "native_label": "7XDS",
  "native": "7XDS AvrSr35 Alone.cif", "model": "fold_7xvg_avrsr35_alone_model_0.cif"}
]
"""
if __name__ == '__main__':
    if len(sys.argv) == 1:
        print(__doc__); print('jobs.json 示例:'); print(EXAMPLE_CONFIG); sys.exit(0)
    main()
 