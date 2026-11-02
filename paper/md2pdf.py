"""Render the note (Markdown) to HTML and PDF with headless Chrome. No pandoc or LaTeX on this machine.
Handles: #/##/### headings, paragraphs, bullet lists, numbered lists, pipe tables, images, bold, italic,
inline code, fenced code, blockquotes, horizontal rules. Usage: python paper/md2pdf.py paper/exact_graphstate_sampling.md"""
import sys, re, html, subprocess, os, pathlib
src = pathlib.Path(sys.argv[1]); text = src.read_text(encoding="utf-8").replace("\r\n", "\n")

def inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<i>\1</i>", s)
    s = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)", r'<img src="\2" alt="\1">', s)
    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', s)
    s = re.sub(r"(?<![\"'>])(https?://[^\s<)]+)", r'<a href="\1">\1</a>', s)
    return s

out = []; lines = text.split("\n"); i = 0; para = []
def flush():
    global para
    if para: out.append("<p>" + inline(" ".join(para)) + "</p>"); para = []
while i < len(lines):
    l = lines[i]
    if l.startswith("```"):
        flush(); j = i + 1; buf = []
        while j < len(lines) and not lines[j].startswith("```"): buf.append(lines[j]); j += 1
        out.append("<pre>" + html.escape("\n".join(buf)) + "</pre>"); i = j + 1; continue
    m = re.match(r"^(#{1,3})\s+(.*)", l)
    if m: flush(); k = len(m.group(1)); out.append(f"<h{k}>{inline(m.group(2))}</h{k}>"); i += 1; continue
    if l.strip() == "---": flush(); out.append("<hr>"); i += 1; continue
    if l.startswith("|"):
        flush(); rows = []
        while i < len(lines) and lines[i].startswith("|"): rows.append(lines[i]); i += 1
        cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
        cells = [c for c in cells if not all(re.fullmatch(r":?-+:?", x) for x in c)]
        t = "<table><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells[0]) + "</tr>"
        for r in cells[1:]: t += "<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>"
        out.append(t + "</table>"); continue
    if re.match(r"^\s*[-*]\s+", l) or re.match(r"^\s*\d+\.\s+", l):
        # A list item is a sequence of blocks: paragraphs and pipe tables. Continuation lines are indented by at least
        # two spaces; a blank line inside an item starts a new paragraph (so a numbered list separated by blank lines
        # stays one list and does not restart at 1), and an indented table inside an item renders as a table.
        flush(); ordered = bool(re.match(r"^\s*\d+\.", l)); marker = r"^\s*\d+\.\s+" if ordered else r"^\s*[-*]\s+"
        items = []   # each item: list of blocks, a block being ["p", text] or ["t", [rows]]
        while i < len(lines):
            cur = lines[i]
            if re.match(marker, cur): items.append([["p", re.sub(marker, "", cur)]]); i += 1; continue
            if cur.strip() == "":
                j = i + 1
                while j < len(lines) and lines[j].strip() == "": j += 1
                if j < len(lines) and items and (re.match(marker, lines[j]) or lines[j].startswith("   ")):
                    items[-1].append(["p", ""]); i = j; continue
                break
            if items and cur.startswith("  "):
                s = cur.strip()
                if s.startswith("|"):
                    if items[-1][-1][0] == "t": items[-1][-1][1].append(s)
                    else: items[-1].append(["t", [s]])
                else:
                    if items[-1][-1][0] == "p": items[-1][-1][1] = (items[-1][-1][1] + " " + s).strip()
                    else: items[-1].append(["p", s])
                i += 1; continue
            break
        def table_html(rows):
            cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
            cells = [c for c in cells if not all(re.fullmatch(r":?-+:?", x) for x in c)]
            t = "<table><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells[0]) + "</tr>"
            for r in cells[1:]: t += "<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>"
            return t + "</table>"
        def item_html(blocks):
            parts = []
            for kind, body in blocks:
                if kind == "t": parts.append(table_html(body))
                elif body: parts.append(inline(body) if not parts else f"<p>{inline(body)}</p>")
            return "".join(parts)
        tag = "ol" if ordered else "ul"
        out.append(f"<{tag}>" + "".join(f"<li>{item_html(x)}</li>" for x in items) + f"</{tag}>"); continue
    if l.startswith(">"):
        flush(); buf = []
        while i < len(lines) and lines[i].startswith(">"): buf.append(lines[i][1:].strip()); i += 1
        out.append("<blockquote>" + inline(" ".join(buf)) + "</blockquote>"); continue
    if l.strip() == "": flush(); i += 1; continue
    para.append(l.strip()); i += 1
flush()
css = """body{font-family:'Georgia','Times New Roman',serif;font-size:10.5pt;line-height:1.38;max-width:17cm;margin:1.6cm auto;color:#111}
h1{font-size:17pt;line-height:1.25;margin-bottom:4pt}h2{font-size:13pt;margin-top:18pt;border-bottom:1px solid #999;padding-bottom:2pt}h3{font-size:11pt}
table{border-collapse:collapse;margin:8pt 0;font-size:9.5pt}th,td{border:1px solid #888;padding:2pt 6pt;text-align:left}th{background:#eee}
code{font-family:Consolas,monospace;font-size:9pt;background:#f3f3f3;padding:0 2px}pre{font-family:Consolas,monospace;font-size:8.5pt;background:#f6f6f6;padding:6pt;white-space:pre-wrap}
img{max-width:100%;display:block;margin:8pt auto}blockquote{border-left:3px solid #bbb;margin:6pt 0;padding-left:8pt;color:#333}
li{margin:2pt 0}a{color:#1a4a8a;text-decoration:none}@page{size:A4;margin:1.8cm}"""
htmlp = src.with_suffix(".html"); pdfp = src.with_suffix(".pdf")
htmlp.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>{html.escape(src.stem)}</title><style>{css}</style></head><body>" + "\n".join(out) + "</body></html>", encoding="utf-8")
chrome = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
before = pdfp.stat().st_mtime if pdfp.exists() else None
r = subprocess.run([chrome, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={pdfp.resolve()}", str(htmlp.resolve().as_uri())], capture_output=True, text=True, timeout=120)
written = "bytes written to file" in r.stderr and pdfp.exists() and pdfp.stat().st_mtime != before
print("html:", htmlp, "pdf:", pdfp, "bytes:", pdfp.stat().st_size if written else "NOT WRITTEN")
if not written:   # Chrome reports a locked or unwritable target only on stderr and exits 0; do not report a stale file as built
    print(r.stderr[-400:]); sys.exit(1)
