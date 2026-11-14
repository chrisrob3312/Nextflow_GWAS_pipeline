#!/usr/bin/env python3
"""
Static consistency audit for the Nextflow code in this repository.

Checks, across main.nf, subworkflows/**/*.nf and modules/**/*.nf:
  * every `include { X } from '...'` points at a file that exists and defines X
  * every call `X( ... )` of a known process / workflow passes as many
    top-level arguments as X declares inputs / takes
  * every `X.out.<name>` names an emit that X actually declares

Usage:  python3 tests/nf_audit.py          (exit 1 when problems are found)
It is a structural check only; it does not run Nextflow.
"""
import glob
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

files = glob.glob("main.nf") + glob.glob("subworkflows/**/*.nf", recursive=True) + glob.glob("modules/**/*.nf", recursive=True)
src = {f: open(f).read() for f in files}


def strip_comments(s):
    # keep line numbers: block comments become the same number of newlines
    s = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group(0).count("\n"), s, flags=re.S)
    s = re.sub(r"^\s*//.*$", "", s, flags=re.M)
    return re.sub(r"(?<=[^:'\"])//[^\n]*", "", s)


defs = {}  # name -> dict(kind, file, n_in, takes, emits)
for f, s in src.items():
    s2 = strip_comments(s)
    for m in re.finditer(r"^(process|workflow)\s+(\w+)\s*\{", s2, flags=re.M):
        kind, name = m.group(1), m.group(2)
        i, depth = m.end(), 1
        while i < len(s2) and depth:
            if s2[i] == "{":
                depth += 1
            elif s2[i] == "}":
                depth -= 1
            i += 1
        body = s2[m.end():i - 1]
        if kind == "process":
            inp = re.search(r"^\s*input:\s*\n(.*?)^\s*(output|when|script|shell|exec):", body, flags=re.S | re.M)
            n_in = len([l for l in (inp.group(1).split("\n") if inp else []) if re.match(r"\s*(tuple|path|val|env|stdin|each)\b", l)])
            emits = re.findall(r"emit:\s*(\w+)", body)
            defs[name] = dict(kind=kind, file=f, n_in=n_in, emits=emits, takes=[])
        else:
            take = re.search(r"^\s*take:\s*\n(.*?)^\s*main:", body, flags=re.S | re.M)
            takes = [l.strip().split()[0] for l in (take.group(1).split("\n") if take else []) if l.strip()]
            emit = re.search(r"^\s*emit:\s*\n(.*)$", body, flags=re.S | re.M)
            emits = re.findall(r"^\s*(\w+)\s*=", emit.group(1), flags=re.M) if emit else []
            defs[name] = dict(kind=kind, file=f, n_in=len(takes), emits=emits, takes=takes)


def top_level_args(text):
    args, depth, cur, q = [], 0, "", None
    for ch in text:
        if q:
            if ch == q:
                q = None
            cur += ch
            continue
        if ch in "\"'":
            q = ch
            cur += ch
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            args.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        args.append(cur)
    return [a for a in args if a.strip()]


problems = []
for f, s in src.items():
    s2 = strip_comments(s)
    for m in re.finditer(r"include\s*\{\s*([\w\s,]+?)\s*\}\s*from\s*'([^']+)'", s2):
        names = [n.strip().split(" as ")[0] for n in m.group(1).split(",")]
        rel = m.group(2)
        p = os.path.normpath(os.path.join(os.path.dirname(f), rel))
        tgt = p + "/main.nf" if os.path.isdir(p) else (p if p.endswith(".nf") else p + ".nf")
        if not os.path.exists(tgt):
            problems.append(f"{f}: include target missing: {rel}")
            continue
        body = strip_comments(open(tgt).read())
        for n in names:
            if not re.search(r"^(process|workflow)\s+" + n + r"\b", body, flags=re.M):
                problems.append(f"{f}: '{n}' not defined in {tgt}")
    for name, d in defs.items():
        for m in re.finditer(r"(?<![\w.])" + name + r"\s*\(", s2):
            line_start = s2.rfind("\n", 0, m.start()) + 1
            line = s2[line_start:m.start()]
            if re.search(r"(process|workflow)\s*$", line) or line.strip().startswith(("include", "def ")):
                continue
            i, depth = m.end(), 1
            while i < len(s2) and depth:
                if s2[i] == "(":
                    depth += 1
                elif s2[i] == ")":
                    depth -= 1
                i += 1
            args = top_level_args(s2[m.end():i - 1])
            if len(args) != d["n_in"]:
                ln = s2[:m.start()].count("\n") + 1
                extra = f" takes={d['takes']}" if d["takes"] else ""
                problems.append(f"{f}:{ln}: {name}() called with {len(args)} args, {d['kind']} declares {d['n_in']}{extra}")
    for m in re.finditer(r"(\w+)\.out\.(\w+)", s2):
        name, em = m.groups()
        if name in defs and em not in defs[name]["emits"]:
            ln = s2[:m.start()].count("\n") + 1
            problems.append(f"{f}:{ln}: {name}.out.{em} not among emits {defs[name]['emits']}")

n_proc = sum(1 for d in defs.values() if d["kind"] == "process")
n_wf = sum(1 for d in defs.values() if d["kind"] == "workflow")
print(f"Scanned {len(files)} files, {n_proc} processes, {n_wf} workflows")
print(f"{len(problems)} problems")
for p in sorted(set(problems)):
    print(" -", p)
sys.exit(1 if problems else 0)
