#%% fix_ck_numbers.py
from pathlib import Path
import re

#%%

mech_in = Path("path_to_chem_file.inp.txt").read_text(errors="ignore")
therm_in_lines = Path("C:/Users/raubry/Desktop/Combustor-HX/combustion_chemistry/gasoline_surrogate/gasoline_surrogate_therm.dat.txt").read_text(errors="ignore").splitlines()

# (a) Fix kinetics numbers like 0.0000+03 -> 0.0000E+03
mech_fixed = re.sub(r'(\d+\.\d+)([+-]\d{2,3})(?=(\s|$))', r'\1E\2', mech_in)
Path("path_to_fixed_chem_file.inp.txt").write_text(mech_fixed)

# (b) Grab the SPECIES list from the mechanism (to avoid removing anything you actually use)
m = re.search(r'(?ms)^\s*SPECIES\s*$([\s\S]*?)^\s*END\s*$', mech_fixed)
species_declared = set()
if m:
    block = "\n".join([ln for ln in m.group(1).splitlines() if not ln.strip().startswith("!")])
    species_declared = set(re.findall(r'([A-Za-z0-9\-\+\(\),]+)', block))

# (c) De-duplicate thermo: define a "header" as any line whose last non-space char is '1'
header_idxs = [i for i,l in enumerate(therm_in_lines) if l.rstrip().endswith('1')]
seen = {}
del_ranges = []

for k, h in enumerate(header_idxs):
    name = therm_in_lines[h].strip().split()[0]
    if name in seen:
        end = header_idxs[k+1] if k+1 < len(header_idxs) else len(therm_in_lines)
        del_ranges.append((h, end))
    else:
        seen[name] = h

skip = [False]*len(therm_in_lines)
for a,b in del_ranges:
    for j in range(a,b): skip[j] = True

therm_dedup = [ln for j,ln in enumerate(therm_in_lines) if not skip[j]]

# (d) Optional: sanity check that all declared species exist in the deduped thermo
names_dedup = {ln.strip().split()[0] for ln in therm_dedup if ln.rstrip().endswith('1')}
missing = sorted(species_declared - names_dedup)
if missing:
    print(f"WARNING: {len(missing)} species declared in the mechanism have no thermo entry:", missing[:10], "...")

Path("gasoline_surrogate_therm_dedup3.dat.txt").write_text("\n".join(therm_dedup))
print("Wrote Chem323_fixed.inp.txt and gasoline_surrogate_therm_dedup3.dat.txt")

#%%

# patch_nasa4.py
# only for reduced mechanism where the transport file needs to be adapted with less species

p = Path("C:/Users/raubry/Desktop/Combustor-HX/combustion_chemistry/gasoline_surrogate/gasoline_surrogate_therm.dat.txt")
lines = p.read_text(errors="ignore").splitlines()

def is_numline(s):
    # heuristic: at least 4 Fortran floats on the line
    return len(re.findall(r'[ \-+]\d\.\d+(?:[EeDd][\+\-]\d{2,3})', ' '+s)) >= 4

fixed = []
i = 0
while i < len(lines):
    fixed.append(lines[i])
    # detect a header-like line for a NASA7 block (next 3 non-comment numeric lines exist)
    j = i + 1
    num = []
    raw = []
    while j < len(lines) and len(num) < 3 and (j - i) <= 10:
        s = lines[j]
        if s.strip() and not s.strip().startswith('!'):
            raw.append((j, s))
            if is_numline(s):
                num.append((j, s))
        j += 1
    if len(num) == 3:
        # we just copied lines[i]; now ensure the 3 numeric lines end with 2/3/4 markers
        for k,(idx,s) in enumerate(num, start=2):
            # pad to col 80 with the correct marker if missing/wrong
            s2 = s.rstrip('\n')
            if len(s2) >= 80:
                s2 = s2[:79] + str(k)
            else:
                s2 = s2 + ' ' * (80 - len(s2) - 1) + str(k)
            fixed.append(s2)
        # skip over the original numeric lines we just re-emitted
        i = num[-1][0] + 1
        continue
    i += 1

Path("C:/Users/raubry/Desktop/Combustor-HX/combustion_chemistry/gasoline_surrogate/gasoline_surrogate_therm_patched.dat.txt").write_text("\n".join(fixed) + "\n")
print("Wrote gasoline_surrogate_therm_patched.txt")

#%%