from pymol import cmd, stored
import math
import os

print("\n=============================================================")
print("  SS-Manager v16.0: The Ultimate Disulfide Suite")
print("  Integrates V11 (Compare), V12 (Transfer), V15 (CSV/Snap)")
print("=============================================================\n")

# --- 1. 通用辅助函数 ---
AA_MAP = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'CYX':'C',
    'GLN':'Q', 'GLU':'E', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L',
    'LYS':'K', 'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S', 'THR':'T',
    'TRP':'W', 'TYR':'Y', 'VAL':'V'
}

def get_dist(c1, c2):
    return math.sqrt(sum((c1[i]-c2[i])**2 for i in range(3)))

def get_bond_midpoint(c1, c2):
    return [(c1[i]+c2[i])/2 for i in range(3)]

def get_seq_context(obj, chain, resi, window=3):
    """提取CYS残基周围的序列指纹"""
    resi_int = int(resi); start = resi_int - window; end = resi_int + window
    seq_str = ""; myspace = {'residues': []}
    try:
        cmd.iterate_state(1, f"({obj}) and chain {chain} and resi {start}-{end} and name CA", 
                          "residues.append(resn)", space=myspace)
    except: return "X" * (window*2 + 1)
    for r in myspace['residues']: seq_str += AA_MAP.get(r, 'X')
    return seq_str

def calculate_seq_score(seq1, seq2):
    if not seq1 or not seq2: return 0.0
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    length = max(len(seq1), len(seq2))
    return matches / length if length > 0 else 0.0

# --- 2. 核心扫描逻辑 (Topology vs Geometry) ---
def scan_true_topology_bonds(obj_name):
    """(Ref) 从 PyMOL 内部拓扑提取真实的键 (V11核心)"""
    model = cmd.get_model(f"({obj_name}) and resn CYS+CYX and name SG")
    idx_to_atom = {i: at for i, at in enumerate(model.atom)}
    bonds_found = []; processed = set()
    if hasattr(model, 'bond'):
        for bd in model.bond:
            at1 = idx_to_atom[bd.index[0]]; at2 = idx_to_atom[bd.index[1]]
            if at1.name != 'SG' or at2.name != 'SG': continue
            k1 = (at1.chain, at1.resi); k2 = (at2.chain, at2.resi)
            pair_key = tuple(sorted([k1, k2]))
            if pair_key not in processed:
                processed.add(pair_key)
                s1 = get_seq_context(obj_name, at1.chain, at1.resi)
                s2 = get_seq_context(obj_name, at2.chain, at2.resi)
                bonds_found.append({
                    'obj': obj_name, 'atom1': {'chain': at1.chain, 'resi': at1.resi, 'coord': at1.coord}, 
                    'atom2': {'chain': at2.chain, 'resi': at2.resi, 'coord': at2.coord},
                    'midpoint': get_bond_midpoint(at1.coord, at2.coord),
                    'label': f"{at1.chain}:{at1.resi}-{at2.chain}:{at2.resi}",
                    'seqs': tuple(sorted([s1, s2]))
                })
    return bonds_found

def scan_existing_bonds_by_distance(obj_name, cutoff=3.0):
    """(Target) 基于几何距离扫描键 (V11核心)"""
    sg_atoms = []; cmd.iterate_state(1, f"({obj_name}) and resn CYS+CYX and name SG",
                         'sg_atoms.append({"chain":chain, "resi":resi, "coord":(x,y,z), "id":ID})',
                         space={'sg_atoms': sg_atoms})
    found_bonds = []; processed_pairs = set(); cutoff_float = float(cutoff)
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1 = sg_atoms[i]; a2 = sg_atoms[j]
            if a1['chain'] == a2['chain'] and abs(int(a1['resi']) - int(a2['resi'])) <= 1: continue
            if get_dist(a1['coord'], a2['coord']) < cutoff_float:
                pair_key = tuple(sorted([(a1['chain'], int(a1['resi'])), (a2['chain'], int(a2['resi']))]))
                if pair_key not in processed_pairs:
                    processed_pairs.add(pair_key)
                    s1 = get_seq_context(obj_name, a1['chain'], a1['resi'])
                    s2 = get_seq_context(obj_name, a2['chain'], a2['resi'])
                    found_bonds.append({
                        'obj': obj_name, 'atom1': a1, 'atom2': a2,
                        'midpoint': get_bond_midpoint(a1['coord'], a2['coord']),
                        'label': f"{a1['chain']}:{a1['resi']}-{a2['chain']}:{a2['resi']}",
                        'seqs': tuple(sorted([s1, s2]))
                    })
    return found_bonds

# ==============================================================================
# 命令 1: ss_compare (来自 V11 - 权威比对)
# ==============================================================================
def ss_compare(ref_obj, target_obj, match_tolerance=3.0):
    print(f"\n[SS_Compare] {ref_obj} (Ref, Topology) vs {target_obj} (Target, Dist)")
    ref_bonds = scan_true_topology_bonds(ref_obj)
    if len(ref_bonds) == 0:
        print("[WARNING] Ref has no topology. Falling back to distance scan.")
        ref_bonds = scan_existing_bonds_by_distance(ref_obj, cutoff=3.0)
    tgt_bonds = scan_existing_bonds_by_distance(target_obj, cutoff=3.0)
    print(f"  > Ref Bonds: {len(ref_bonds)} | Tgt Bonds: {len(tgt_bonds)}")
    
    candidates = []
    for r_idx, r_bond in enumerate(ref_bonds):
        for t_idx, t_bond in enumerate(tgt_bonds):
            dist = get_dist(r_bond['midpoint'], t_bond['midpoint'])
            if dist < float(match_tolerance):
                s1 = (calculate_seq_score(r_bond['seqs'][0], t_bond['seqs'][0]) + calculate_seq_score(r_bond['seqs'][1], t_bond['seqs'][1]))/2.0
                s2 = (calculate_seq_score(r_bond['seqs'][0], t_bond['seqs'][1]) + calculate_seq_score(r_bond['seqs'][1], t_bond['seqs'][0]))/2.0
                candidates.append({'r_idx': r_idx, 't_idx': t_idx, 'dist': dist, 'seq_score': max(s1, s2)})
    
    candidates.sort(key=lambda x: (-x['seq_score'], x['dist']))
    ref_matched = set(); tgt_matched = set(); final_matches = []
    for c in candidates:
        if c['r_idx'] not in ref_matched and c['t_idx'] not in tgt_matched:
            ref_matched.add(c['r_idx']); tgt_matched.add(c['t_idx']); final_matches.append(c)
            t_bond = tgt_bonds[c['t_idx']]; n = f"shared_{c['r_idx']}"
            cmd.distance(n, f"({target_obj}) and chain {t_bond['atom1']['chain']} and resi {t_bond['atom1']['resi']} and name SG",
                            f"({target_obj}) and chain {t_bond['atom2']['chain']} and resi {t_bond['atom2']['resi']} and name SG")
            cmd.color("green", n); cmd.hide("label", n)

    print("-" * 100); print(f"{'STATUS':<15} | {'Ref Bond (Crystal)':<22} | {'Target Bond (Model)':<22} | {'SeqSim':<6} | {'Dist'}")
    print("-" * 100)
    final_matches.sort(key=lambda x: x['r_idx'])
    for m in final_matches:
        r = ref_bonds[m['r_idx']]; t = tgt_bonds[m['t_idx']]
        print(f"{'SHARED':<15} | {r['label']:<22} | {t['label']:<22} | {m['seq_score']:.2f}   | {m['dist']:.1f} A")
    
    for i, r in enumerate(ref_bonds):
        if i not in ref_matched:
            n = f"missing_{i}"; p1=f"p1_{i}"; p2=f"p2_{i}"
            cmd.pseudoatom(p1, pos=r['atom1']['coord']); cmd.pseudoatom(p2, pos=r['atom2']['coord'])
            cmd.distance(n, p1, p2); cmd.color("red", n); cmd.set("dash_gap", 0.5, n)
            cmd.hide("nonbonded", p1); cmd.hide("nonbonded", p2); cmd.hide("label", n)
            print(f"{'MISSING':<15} | {r['label']:<22} | {'---':<22} | -      | -")

    for i, t in enumerate(tgt_bonds):
        if i not in tgt_matched:
            n = f"new_{i}"; cmd.distance(n, f"({target_obj}) and chain {t['atom1']['chain']} and resi {t['atom1']['resi']} and name SG",
                            f"({target_obj}) and chain {t['atom2']['chain']} and resi {t['atom2']['resi']} and name SG")
            cmd.color("blue", n); cmd.hide("label", n)
            print(f"{'NEW/ARTIFACT':<15} | {'---':<22} | {t['label']:<22} | -      | -")
    print("-" * 100)

# ==============================================================================
# 命令 2: autobond_ss (来自 V11 - 快速连线)
# ==============================================================================
def autobond_ss(obj_name, cutoff=3.0):
    print(f"\n[AutoBond] Scanning '{obj_name}' (Cutoff: {cutoff} A)...")
    sg_atoms = []; cmd.iterate_state(1, f"({obj_name}) and resn CYS+CYX and name SG",
                      'sg_atoms.append({"id":ID, "chain":chain, "resi":resi, "coord":(x,y,z)})',
                      space={'sg_atoms': sg_atoms})
    count = 0; cutoff_float = float(cutoff)
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1 = sg_atoms[i]; a2 = sg_atoms[j]
            if a1['chain'] == a2['chain'] and abs(int(a1['resi']) - int(a2['resi'])) <= 1: continue
            if get_dist(a1['coord'], a2['coord']) < cutoff_float:
                try:
                    cmd.bond(f"({obj_name}) and id {a1['id']}", f"({obj_name}) and id {a2['id']}")
                    cmd.color("yellow", f"({obj_name}) and id {a1['id']}+{a2['id']}")
                    count += 1
                except: pass
    cmd.rebuild(); print(f"[AutoBond] Created {count} bonds.\n")

# ==============================================================================
# 命令 3: ss_transfer (来自 V12 - 拓扑克隆)
# ==============================================================================
def ss_transfer(source_obj, target_obj):
    """Copies topology regardless of distance."""
    print(f"\n[Transfer] Cloning topology from '{source_obj}' -> '{target_obj}'...")
    model = cmd.get_model(f"({source_obj}) and resn CYS+CYX and name SG")
    idx_to_atom = {i: at for i, at in enumerate(model.atom)}
    bonds = []; processed = set()
    if hasattr(model, 'bond'):
        for bd in model.bond:
            at1 = idx_to_atom[bd.index[0]]; at2 = idx_to_atom[bd.index[1]]
            if at1.name != 'SG' or at2.name != 'SG': continue
            key = tuple(sorted([(at1.chain, at1.resi), (at2.chain, at2.resi)]))
            if key not in processed: processed.add(key); bonds.append(key)
    
    count = 0; cmd.set("stick_radius", 0.15)
    for b in bonds:
        c1, r1 = b[0]; c2, r2 = b[1]
        sel1 = f"({target_obj}) and chain {c1} and resi {r1} and name SG"
        sel2 = f"({target_obj}) and chain {c2} and resi {r2} and name SG"
        if cmd.count_atoms(sel1) and cmd.count_atoms(sel2):
            try:
                cmd.bond(sel1, sel2); d = cmd.get_distance(sel1, sel2)
                if d > 4.5: cmd.color("red", sel1); cmd.color("red", sel2)
                elif d > 3.0: cmd.color("magenta", sel1); cmd.color("magenta", sel2)
                else: cmd.color("yellow", sel1); cmd.color("yellow", sel2)
                count += 1
            except: pass
    print(f"[Transfer] Cloned {count} bonds to '{target_obj}'.")

# ==============================================================================
# 命令 4-8: V15 新功能 (Export, Import, Snap, Save, Check)
# ==============================================================================
def ss_export(obj_name, filename="bonds.csv"):
    print(f"\n[Export] Scanning '{obj_name}' -> '{filename}'...")
    sg_atoms = []; cmd.iterate_state(1, f"({obj_name}) and resn CYS+CYX and name SG",
                         'sg_atoms.append({"chain":chain, "resi":resi, "coord":(x,y,z)})', space={'sg_atoms': sg_atoms})
    bonds = []; processed = set()
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1 = sg_atoms[i]; a2 = sg_atoms[j]
            if a1['chain'] == a2['chain'] and abs(int(a1['resi']) - int(a2['resi'])) <= 1: continue
            if get_dist(a1['coord'], a2['coord']) < 3.2:
                key = tuple(sorted([(a1['chain'], a1['resi']), (a2['chain'], a2['resi'])]))
                if key not in processed: processed.add(key); bonds.append(key)
    try:
        with open(filename, 'w') as f:
            f.write("Chain1,Resi1,Chain2,Resi2,Note\n")
            for b in bonds: f.write(f"{b[0][0]},{b[0][1]},{b[1][0]},{b[1][1]},Detected\n")
        print(f"[Export] Saved {len(bonds)} bonds.")
    except Exception as e: print(e)

def ss_import(obj_name, filename="bonds.csv"):
    print(f"\n[Import] Applying '{filename}' to '{obj_name}'...")
    try:
        with open(filename, 'r') as f: lines = f.readlines()
        cmd.set("stick_radius", 0.15)
        for line in lines:
            if not line.strip() or line.startswith("Chain") or line.startswith("#"): continue
            parts = line.split(','); c1, r1, c2, r2 = parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[3].strip()
            sel1 = f"({obj_name}) and chain {c1} and resi {r1} and name SG"; sel2 = f"({obj_name}) and chain {c2} and resi {r2} and name SG"
            if cmd.count_atoms(sel1) and cmd.count_atoms(sel2):
                cmd.bond(sel1, sel2); d = cmd.get_distance(sel1, sel2)
                if d > 4.5: cmd.color("red", sel1); cmd.color("red", sel2)
                elif d > 3.0: cmd.color("magenta", sel1); cmd.color("magenta", sel2)
                else: cmd.color("yellow", sel1); cmd.color("yellow", sel2)
        print(f"[Import] Bonds applied visually.")
        cmd.show("sticks", f"({obj_name}) and resn CYS+CYX")
    except Exception as e: print(e)

def ss_snap(fixed_atom_sel, moving_atom_sel, moving_object_sel):
    print(f"\n[Snap] Moving '{moving_object_sel}' towards '{fixed_atom_sel}'...")
    try:
        fixed_model = cmd.get_model(fixed_atom_sel); moving_model = cmd.get_model(moving_atom_sel)
        if len(fixed_model.atom) == 0 or len(moving_model.atom) == 0: print("[Error] Atoms not found."); return
        p1 = fixed_model.atom[0].coord; p2 = moving_model.atom[0].coord
        vec = [p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]]
        cmd.translate(vec, moving_object_sel, camera=0)
        print(f"[Snap] Success! Distance near 0.0 A.")
    except Exception as e: print(f"[Error] Snap failed: {e}")

def ss_check_dist(obj_name, csv_filename="bonds.csv"):
    print(f"\n[Check] Measuring distances in '{obj_name}'...")
    try:
        with open(csv_filename, 'r') as f: lines = f.readlines()
        print("-" * 65); print(f"{'Bond (Chain:Resi)':<25} | {'Dist (A)':<10} | {'Status'}"); print("-" * 65)
        for line in lines:
            if not line.strip() or line.startswith("Chain") or line.startswith("#"): continue
            parts = line.split(','); c1, r1, c2, r2 = parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[3].strip()
            sel1 = f"({obj_name}) and chain {c1} and resi {r1} and name SG"; sel2 = f"({obj_name}) and chain {c2} and resi {r2} and name SG"
            if cmd.count_atoms(sel1) and cmd.count_atoms(sel2):
                d = cmd.get_distance(sel1, sel2); status = "OK (Green)" if d < 2.5 else "STRETCHED (Yellow)" if d < 4.5 else "BROKEN (Red)"
                print(f"{c1}:{r1} - {c2}:{r2}             | {d:<10.2f} | {status}")
            else: print(f"{c1}:{r1} - {c2}:{r2}             | ---        | [MISSING ATOM]")
        print("-" * 65)
    except Exception as e: print(e)

def ss_save_pdb(obj_name, csv_filename, output_pdb):
    print(f"\n[Save] Saving '{obj_name}' with SSBOND to '{output_pdb}'...")
    cmd.save(output_pdb, obj_name); ssbond_lines = []
    try:
        with open(csv_filename, 'r') as f: lines = f.readlines()
        serial = 1
        for line in lines:
            if not line.strip() or line.startswith("Chain") or line.startswith("#"): continue
            parts = line.split(','); c1, r1, c2, r2 = parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[3].strip()
            ssbond_lines.append(f"SSBOND {serial:>3} CYS {c1:>1} {r1:>4}    CYS {c2:>1} {r2:>4}                          1555\n"); serial += 1
        with open(output_pdb, 'r') as f: content = f.readlines()
        clean = [l for l in content if not l.startswith("SSBOND")]
        with open(output_pdb, 'w') as f: f.writelines(ssbond_lines + clean)
        print(f"[Save] Success!")
    except Exception as e: print(f"[Error] {e}")

# Register All
cmd.extend("ss_compare", ss_compare)
cmd.extend("autobond_ss", autobond_ss)
cmd.extend("ss_transfer", ss_transfer)
cmd.extend("ss_export", ss_export)
cmd.extend("ss_import", ss_import)
cmd.extend("ss_snap", ss_snap)
cmd.extend("ss_check_dist", ss_check_dist)
cmd.extend("ss_save_pdb", ss_save_pdb)
