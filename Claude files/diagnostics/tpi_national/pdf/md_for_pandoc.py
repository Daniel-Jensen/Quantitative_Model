# PREPARE A GFM REPORT FOR PANDOC'S OWN MARKDOWN READER (SO PIPE-TABLE CELLS WRAP).
#   python3 md_for_pandoc.py in.md out.md
# 1. horizontal rules '---' become '* * *' (a dash line can open a pandoc multiline table)
# 2. a blank line goes before every list that follows a paragraph line (pandoc needs it, GFM does not)
# 3. each pipe table's separator gets dashes in proportion to its columns' content, so the relative widths
#    pandoc derives from the separator follow the content
import re
import sys

LIST = re.compile(r"^\s*([-*+]|\d+\.)\s+\S")
SEP = re.compile(r"^\|(\s*:?-{3,}:?\s*\|)+\s*$")


def cells(line):
    # THE CELLS OF ONE PIPE-TABLE ROW (ESCAPED PIPES KEPT INSIDE CELLS).
    body = line.strip().strip("|")
    return [c.strip() for c in re.split(r"(?<!\\)\|", body)]


def widths(rows):
    # RELATIVE COLUMN WEIGHTS: THE LONGEST CELL, SOFT-CAPPED SO ONE LONG CELL CANNOT TAKE THE PAGE.
    n = max(len(r) for r in rows)
    w = [0] * n
    for r in rows:
        for i, c in enumerate(r):
            text = re.sub(r"[*`]", "", c)
            w[i] = max(w[i], len(text))
    return [max(4, min(x, 70)) ** 0.85 for x in w]


def fix_tables(lines):
    # REWRITE EVERY PIPE TABLE'S SEPARATOR LINE.
    out, i = [], 0
    while i < len(lines):
        if i + 1 < len(lines) and lines[i].lstrip().startswith("|") and SEP.match(lines[i + 1].strip()):
            j = i + 2
            while j < len(lines) and lines[j].lstrip().startswith("|"):
                j += 1
            rows = [cells(lines[i])] + [cells(l) for l in lines[i + 2:j]]
            seps = cells(lines[i + 1])
            w = widths(rows)
            tot = sum(w)
            new = []
            for k, s in enumerate(seps):
                d = "-" * max(3, round(60 * w[k] / tot))
                new.append((":" if s.startswith(":") else "") + d + (":" if s.endswith(":") else ""))
            out.append(lines[i])
            out.append("|" + "|".join(new) + "|")
            out.extend(lines[i + 2:j])
            i = j
            continue
        out.append(lines[i])
        i += 1
    return out


def main():
    # READ, TRANSFORM, WRITE.
    lines = open(sys.argv[1]).read().split("\n")
    lines = ["* * *" if l.strip() == "---" else l for l in lines]
    fixed = []
    for k, l in enumerate(lines):
        prev = fixed[-1] if fixed else ""
        if LIST.match(l) and prev.strip() and not LIST.match(prev) and not prev.startswith((" ", "\t")) \
                and not prev.lstrip().startswith("|"):
            fixed.append("")
        fixed.append(l)
    open(sys.argv[2], "w").write("\n".join(fix_tables(fixed)))


if __name__ == "__main__":
    main()
