"""Prepare the note for pandoc -> LaTeX (pdflatex, as arXiv compiles it): rewrite the Unicode symbols used in the prose
as raw LaTeX, leaving math segments ($...$, $$...$$), code spans and fenced code untouched. Symbols become
\\ensuremath{...} (raw LaTeX passes through pandoc; $...$ would not when a digit follows the closing dollar).
Writes <name>_tex.md next to the source; build with paper/build_tex.sh.
Usage: python paper/md2tex.py paper/exact_graphstate_sampling.md"""
import sys, re, pathlib
src = pathlib.Path(sys.argv[1]); text = src.read_text(encoding="utf-8").replace("\r\n", "\n")
SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-"}
SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6", "₇": "7", "₈": "8", "₉": "9"}
MATH = {"−": "-", "·": "\\cdot", "±": "\\pm", "⟩": "\\rangle", "⟨": "\\langle", "Σ": "\\Sigma", "½": "\\tfrac12",
        "χ": "\\chi", "√": "\\surd", "⌊": "\\lfloor", "⌋": "\\rfloor", "≤": "\\le", "≈": "\\approx", "σ": "\\sigma",
        "ψ": "\\psi", "Π": "\\Pi", "∈": "\\in", "φ": "\\varphi", "≥": "\\ge", "†": "\\dagger", "→": "\\to",
        "×": "\\times", "θ": "\\theta", "↦": "\\mapsto", "∎": "\\blacksquare", "⊗": "\\otimes", "′": "'",
        "≠": "\\ne", "∪": "\\cup", "δ": "\\delta", "π": "\\pi"}
TEXT = {"–": "--", "—": "---", "…": "\\ldots{}", "§": "\\S{}"}
def mathify(s):
    """Prose written with TeX-style sub/superscripts (2^{35}, r_k, Z^{v}, A_ij, log2 Σ_k 2^{r_k}) becomes real math.
    Braced groups may contain one level of nested braces (they already hold \\ensuremath{...} from the symbol pass)."""
    grp = r"\{(?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*\}"      # braces nested up to two levels (symbol pass adds \ensuremath{})
    s = re.sub(r"(\|[^|{}\n]{1,12}\|)\^(\w)\b", lambda m: "\\ensuremath{" + m.group(1) + "^{" + m.group(2) + "}}", s)   # |n|^w
    s = re.sub(r"(?<![\\\w{])(\d+(?:\.\d+)?|[A-Za-z])((?:\^|_)" + grp + r"(?:(?:\^|_)(?:" + grp + r"|\w))?)", lambda m: "\\ensuremath{" + m.group(1) + m.group(2) + "}", s)
    s = re.sub(r"(?<![\\\w{])([A-Za-z])_([A-Za-z0-9]{1,3})\b", lambda m: "\\ensuremath{" + m.group(1) + "_{" + m.group(2) + "}}", s)
    s = re.sub(r"(?<![\\\w{])(\d+|[A-Za-z])\^([A-Za-z0-9])\b", lambda m: "\\ensuremath{" + m.group(1) + "^{" + m.group(2) + "}}", s)   # Z^v, 2^r
    s = re.sub(r"(\((?:[^()\n]|\\ensuremath\{[^{}]*\})*\))\^(" + grp + r"|\w)", lambda m: "\\ensuremath{" + m.group(1) + "^" + m.group(2) + "}", s)   # (-1)^{...}
    s = re.sub(r"(?<![\\\w{])(\d+)\^(\d+(?:\.\d+)?)\b", lambda m: "\\ensuremath{" + m.group(1) + "^{" + m.group(2) + "}}", s)   # 2^16.7 in tables
    # a symbol that became \ensuremath{...} followed by its own sub/superscript, e.g. Σ_k or Π_j: fold the script inside
    s = re.sub(r"\\ensuremath\{([^{}]*)\}((?:\^|_)(?:" + grp + r"|\w))", lambda m: "\\ensuremath{" + m.group(1) + m.group(2) + "}", s)
    return s
def convert_text(s):
    s = re.sub("[" + "".join(SUP) + "]+", lambda m: "\\textsuperscript{" + "".join(SUP[c] for c in m.group(0)) + "}", s)
    s = re.sub("[" + "".join(SUB) + "]+", lambda m: "\\textsubscript{" + "".join(SUB[c] for c in m.group(0)) + "}", s)
    out = []
    for ch in s:
        if ch in MATH: out.append("\\ensuremath{" + MATH[ch] + "}")
        elif ch in TEXT: out.append(TEXT[ch])
        else: out.append(ch)
    return mathify("".join(out))
def convert(block):
    parts = re.split(r"(```.*?```)", block, flags=re.S); res = []
    for p in parts:
        if p.startswith("```"): res.append(p); continue
        segs = re.split(r"(\$\$.*?\$\$|\$[^$\n]+\$|`[^`\n]+`)", p, flags=re.S)
        res.append("".join(seg if (seg.startswith("$") or seg.startswith("`")) else convert_text(seg) for seg in segs))
    return "".join(res)
out = convert(text)
# pandoc needs a blank line before a list; the source sometimes starts one right after a paragraph line.
lines = out.split("\n"); fixed = []
for i, l in enumerate(lines):
    if i > 0 and re.match(r"^\s*([-*]|\d+\.)\s+", l) and lines[i - 1].strip() and not re.match(r"^\s*([-*]|\d+\.|\|)", lines[i - 1]) and not lines[i - 1].startswith("  "):
        fixed.append("")
    fixed.append(l)
out = "\n".join(fixed)
# Title block -> pandoc metadata (the first heading is the title; the next non-empty line is the author line).
lines = out.split("\n"); title = lines[0].lstrip("# ").strip(); rest = lines[1:]
while rest and not rest[0].strip(): rest.pop(0)
author_line = rest.pop(0).strip()          # "**Name** — affiliation — ORCID ... — email"
fields = [f.strip().strip("*") for f in author_line.split("---")]
author = fields[0]; affil = " \\ensuremath{\\cdot} ".join(fields[1:])
out = ("---\ntitle: \"" + title.replace('"', "'") + "\"\nauthor: \"" + author.replace('"', "'") + "\"\ndate: \"\"\n---\n"
       + (affil[0].upper() + affil[1:] + "\n\n" if affil else "") + "\n".join(rest))
# Images: drop the alt text so pandoc does not add its own numbered caption; the note's own caption paragraph follows.
out = re.sub(r"!\[[^\]]*\]\(", "![](", out)
dst = src.with_name(src.stem + "_tex.md"); dst.write_text(out, encoding="utf-8")
left = sorted({c for c in out if ord(c) > 127 and c not in "éüő"})
print("wrote", dst, "| non-ASCII left:", "".join(left) or "none")
