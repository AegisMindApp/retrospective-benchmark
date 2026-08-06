"""Prove the resume path: cached compounds must NOT be re-docked."""
import json, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))   # portable: no hardcoded paths
sys.path.insert(0, HERE)
os.chdir(HERE)
import run_pilot_fast as rpf
import vina_cache as vc

# ProcessPoolExecutor forks: an in-memory list in the child is invisible to the
# parent. Every real dock appends to a file the parent reads back afterwards.
LEDGER = os.path.join(tempfile.mkdtemp(), "docked.txt")
def fake_dock_one(a):
    with open(LEDGER, "a") as fh:
        fh.write(a[1] + "\n")
    return -7.0
rpf._dock_one = fake_dock_one
rpf.N_WORKERS = 2

names, smiles = [f"C{i}" for i in range(10)], ["CCO"] * 10
cache = os.path.join(tempfile.mkdtemp(), "c.json")
json.dump({f"bench_C{i}": -10.0 - i for i in range(4)}, open(cache, "w"))  # C0..C3 done

scores = rpf.parallel_vina_scores(smiles, "rec.pdbqt", [0,0,0], [20,20,20],
                                  names, exhaustiveness=8, cache_path=cache)

print("\n--- assertions ---")
ok = True
DOCKED = sorted(open(LEDGER).read().split()) if os.path.exists(LEDGER) else []
redocked = [n for n in DOCKED if n in {"C0","C1","C2","C3"}]
print(f"1 cached re-docked   : {redocked or 'NONE'}  {'OK' if not redocked else 'FAIL'}")
ok &= not redocked
want = [f"C{i}" for i in range(4,10)]
print(f"2 actually docked    : {DOCKED}  {'OK' if DOCKED==want else 'FAIL want '+str(want)}")
ok &= DOCKED == want
exp = [-10.0,-11.0,-12.0,-13.0] + [-7.0]*6
print(f"3 positional order   : {scores}  {'OK' if scores==exp else 'FAIL'}")
ok &= scores == exp
final = json.load(open(cache))
missing = [f"bench_C{i}" for i in range(4,10) if f"bench_C{i}" not in final]
print(f"4 new persisted      : {len(final)} keys, missing={missing or 'NONE'}  {'OK' if not missing else 'FAIL'}")
ok &= not missing

# 5-7: the protocol guard must refuse to resume across a changed protocol.
pc = os.path.join(tempfile.mkdtemp(), "p.json")
proto = {"exhaustiveness": 32, "receptor": "r.pdbqt", "flexible_ligand": True}
print(f"5 fresh cache        : {vc.check_protocol(pc, proto)[0]}  "
      f"{'OK' if vc.check_protocol(pc, proto)[0] in ('new','ok') else 'FAIL'}")
ok &= vc.check_protocol(pc, proto)[0] == "ok"
same = vc.check_protocol(pc, proto)[0]
print(f"6 same protocol      : {same}  {'OK' if same=='ok' else 'FAIL'}")
ok &= same == "ok"
bad = vc.check_protocol(pc, dict(proto, exhaustiveness=8))[0]
print(f"7 exh 32->8 refused  : {bad}  {'OK' if bad=='mismatch' else 'FAIL'}")
ok &= bad == "mismatch"
rigid = vc.check_protocol(pc, dict(proto, flexible_ligand=False))[0]
print(f"8 rigid vs flex      : {rigid}  {'OK' if rigid=='mismatch' else 'FAIL'}")
ok &= rigid == "mismatch"
# 9: scores present but no sidecar -> 'legacy', never silently 'ok'
lc = os.path.join(tempfile.mkdtemp(), "l.json")
json.dump({"bench_X": -8.0}, open(lc, "w"))
leg = vc.check_protocol(lc, proto)[0]
print(f"9 legacy detected    : {leg}  {'OK' if leg=='legacy' else 'FAIL'}")
ok &= leg == "legacy"

print("\nRESULT:", "PASS - resume skips/orders/merges and refuses protocol mixing" if ok else "FAIL")
sys.exit(0 if ok else 1)
