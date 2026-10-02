#!/usr/bin/env python3
"""Every in-text citation (Surname year) must have a reference entry, every entry must be cited, and
every entry's author list must match its INSPIRE record (refs_inspire.json).
Usage: python3 cite_check.py vacuum_popcorn.md     |    python3 cite_check.py --selftest"""
import re, sys
def surname_key(s):
    s = s.replace("Bödeker", "Bodeker").replace("Niţu", "Nitu")
    return s.split()[0].strip(",").lower()
def parse(text):
    body, refs = text.split("\n## References", 1)
    body = re.sub(r"<!--.*?-->", "", body, flags=re.S)
    body = re.sub(r"\s+", " ", body)                       # citations wrap across line breaks
    cites = set()
    # "Kawana, Lu & Xie 2022", "Carenza et al. (2024)", "(Pierre Auger Collaboration 2017)", "Lu, Kawana & Kusenko 2023"
    for m in re.finditer(r"([A-Z][A-Za-zöţü\-]+(?: [A-Z][A-Za-z\-]+)?)(?:,? [A-Z][A-Za-zöţü\-]+)*(?: & [A-Z][A-Za-zöţü\-]+(?: [A-Z][A-Za-z\-]+)?)?(?: et\s+al\.)?\s*\(?((?:19|20)\d\d[a-z]?)\)?", body):
        first, year = m.group(1), m.group(2)
        if first in ("May", "Figure", "Section", "Telescope", "Pierre", "Auger") and "Collaboration" not in m.group(0): pass
        cites.add((surname_key(first), year))
    for m in re.finditer(r"(Telescope Array|Pierre Auger) Collaboration\s*\(?((?:19|20)\d\d[a-z]?)", body):
        cites.add((m.group(1).split()[0].lower(), m.group(2)))
    entries = set()
    for line in refs.splitlines():
        m = re.match(r"\s*([A-Z][^\s,]+)[^\n]*?,\s*((?:19|20)\d\d[a-z]?),", line)
        if m: entries.add((surname_key(m.group(1)), m.group(2)))
    return cites, entries
def authors_check(text, inspire="refs_inspire.json"):
    """Each entry's author list against its INSPIRE record, and the in-text form against the count.

    The first version of this file checked surnames and years only, so two entries went out one
    author short (Lewicki et al. 2023 without Veermäe, Niţu et al. 2021 without Scaife) and the
    second was cited in the text as a three-author paper. One, two or three authors are written
    out; four or more are "et al.".
    """
    import json, os
    import unicodedata
    plain = lambda x: "".join(c for c in unicodedata.normalize("NFD", x) if not unicodedata.combining(c))
    recs = {r["arxiv"]: r for r in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), inspire)))}
    body, refs = text.split("\n## References", 1)
    body = re.sub(r"\s+", " ", re.sub(r"<!--.*?-->", "", body, flags=re.S))
    problems, checked, forms = [], 0, {}
    for line in refs.splitlines():
        m = re.match(r"\s*(.+?),\s*((?:19|20)\d\d[a-z]?),.*?arXiv:([^\s)]+)", line)
        if not m or "Collaboration" in m.group(1): continue
        names, year, ax = m.group(1), m.group(2), m.group(3)
        if ax not in recs: problems.append("no INSPIRE record for arXiv:%s" % ax); continue
        checked += 1
        want = [a.split(",")[0] for a in recs[ax]["authors"]]; n = recs[ax]["n_authors"]
        if "et al" not in names:
            listed = len(re.findall(r"(?:^|,\s*)[A-Z][^,]*?\s(?:[A-Z]\.\s?-?)+", names))
            if listed != n:
                problems.append("%s %s: entry lists %d authors, INSPIRE has %d (%s)"
                                % (want[0], year, listed, n, ", ".join(want)))
        form = (want[0] if n == 1 else "%s & %s" % (want[0], want[1]) if n == 2 else
                "%s, %s & %s" % (want[0], want[1], want[2]) if n == 3 else "%s et al." % want[0])
        forms.setdefault((plain(want[0]), year), set()).add(plain(form))
    seen = 0
    for (first, year), ok in sorted(forms.items()):      # two papers may share a first author and a year
        pat = r"(?<!, )(?<!& )" + re.escape(first) + r"(?:, [A-Z][^\s,()]+(?: [A-Z][a-z]+)?| & [A-Z][^\s,()]+(?:-[A-Z][a-z]+)?| et al\.)*\s*\(?" + year   # (?<!, )(?<!& ): not a later name in another paper's list
        for c in re.finditer(pat, plain(body)):
            seen += 1
            got = re.sub(r"\s*\(?" + year + "$", "", c.group(0))
            if got not in ok:
                problems.append("cited as '%s %s', should be '%s %s'" % (got, year, "' or '".join(sorted(ok)), year))
    return problems, checked, seen
def check(text):
    cites, entries = parse(text)
    real = {c for c in cites if c in entries or c[0] in {e[0] for e in entries}}
    missing = sorted(c for c in real if c not in entries)        # surname known, year differs
    unknown = sorted(c for c in cites if c[0] not in {e[0] for e in entries})
    uncited = sorted(e for e in entries if e not in cites)
    return missing, uncited, unknown
if __name__ == "__main__":
    if "--selftest" in sys.argv:
        t = open("vacuum_popcorn.md").read()
        m, u, _ = check(t)
        planted = t.replace("Arnold\n1993", "Arnold\n1994", 1).replace("Arnold 1993", "Arnold 1994", 1)
        m2, u2, _ = check(planted)
        planted2 = t.replace("(Carlson et al. 2005)", "(somewhere)", 1)
        m3, u3, _ = check(planted2)
        assert (m2 or u2) and u3, (m2, u2, u3)
        print("selftest: a wrong year and a dropped citation are both caught")
        p0, n0, c0 = authors_check(t)
        drop = t.replace("Kawana K., Lu P., Xie K.-P., 2022", "Kawana K., Xie K.-P., 2022", 1)
        form = t.replace("Carenza et al. (2024)", "Carenza, Eby & Iarygina (2024)", 1)
        assert drop != t and form != t
        p1 = authors_check(drop)[0]; p2 = authors_check(form)[0]
        assert n0 >= 20 and c0 >= 20 and len(p1) > len(p0) and len(p2) > len(p0), (n0, p0, p1, p2)
        print("selftest: a dropped author and a wrong in-text form are both caught (%d entries, %d in-text citations checked)" % (n0, c0))
        sys.exit(0)
    m, u, unk = check(open(sys.argv[1]).read())
    for c in m: print("cited with a year that has no entry:", c)
    for e in u: print("entry never cited:", e)
    ap, n, nc = authors_check(open(sys.argv[1]).read())
    for a in ap: print("authors:", a)
    print("%d entries and %d in-text citations checked against INSPIRE author lists" % (n, nc))
    print("PASS" if not (m or u or ap) else "FAIL")
    sys.exit(1 if (m or u or ap) else 0)
