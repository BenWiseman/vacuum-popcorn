#!/usr/bin/env python3
"""Enforce STYLE.md on the letter: image counts within their caps, nothing from the banned lists.
Usage: python3 style_check.py vacuum_popcorn.md        (exit 1 on any breach)
       python3 style_check.py --selftest               (plants breaches and checks they are caught)"""
import re, sys
CAPS = {r"\brecipes?\b": 2, r"\bingredients?\b": 3, r"\bpips?\b": 3, r"\bpopcorn\b": 4}
COUNT_ONLY = {r"\bpops?\b": None}
BANNED = [r"\bcosts?\b", r"\bbuys?\b", r"\bpays?\b", r"\bowes?\b", r"\bbudget", r"\bprices?\b", r"\bcheap",
          r"\bexpensive\b", r"\bworth\b", r"\bspend", r"\bearns?\b", r"\bafford", r"\bbargain", r"trade-off",
          r"recipe for disaster", r"food for thought", r"half-baked", r"proof of the pudding", r"boils? down",
          r"cook(?:s|ed|ing)? up", r"\bsimmer", r"back burner", r"spill the beans", r"cherry-pick", r"sweet spot",
          r"icing on", r"bread and butter", r"piece of cake", r"melting pot", r"a taste of", r"whet",
          r"\bflavou?rs?\b", r"\bkernels?\b", r"popcorn mechanism", "—", r"\ba read\b", r"\ba tell\b",
          r"through-line"]
def strip(t):
    t = re.sub(r"<!--.*?-->", " ", t, flags=re.S)          # provenance comments
    t = re.sub(r"\$\$.*?\$\$", " ", t, flags=re.S); t = re.sub(r"\$[^$\n]*\$", " ", t)
    t = t.split("\n## References")[0]                       # titles of cited papers are not ours
    t = re.sub(r"https?://\S+", " ", t)                     # a web address is not an image (the repository is named vacuum-popcorn)
    return t
def check(text):
    t = strip(text); bad = []
    for pat, cap in CAPS.items():
        n = len(re.findall(pat, t, flags=re.I))
        if n > cap: bad.append(f"image over its cap: /{pat}/ used {n} times (cap {cap})")
    for pat in BANNED:
        for m in re.finditer(pat, t, flags=re.I):
            s = max(0, m.start() - 40); bad.append(f"banned /{pat}/: ...{t[s:m.end() + 40]!r}...")
    counts = {pat: len(re.findall(pat, t, flags=re.I)) for pat in list(CAPS) + list(COUNT_ONLY)}
    return bad, counts
if __name__ == "__main__":
    if "--selftest" in sys.argv:
        clean = "A recipe with ingredients. Pinch a pip. A pan of popcorn. One pop."
        assert check(clean)[0] == [], check(clean)[0]
        for planted in ["This costs nothing.", "It boils down to mass.", "pip pip pip pip.", "a flavour of",
                        "a gap — here", "popcorn popcorn popcorn popcorn"]:
            assert check(clean + " " + planted)[0], f"missed: {planted}"
        print("selftest: 6 planted breaches caught, clean text passes"); sys.exit(0)
    text = open(sys.argv[1], encoding="utf-8").read()
    bad, counts = check(text)
    print("image counts:", ", ".join(f"{re.sub(r'[^a-z]', '', k.replace('\\b', ''))}={v}" for k, v in counts.items()))
    for b in bad: print("  " + b)
    print("PASS" if not bad else f"FAIL: {len(bad)} breach(es)")
    sys.exit(1 if bad else 0)
