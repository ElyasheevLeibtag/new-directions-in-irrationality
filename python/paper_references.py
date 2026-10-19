"""Build a transparent label-to-number map for the supplied paper snapshot.

The section-local theorem counter is shared by the theorem, lemma, remark,
definition, corollary and proposition environments. This is the convention in
the project's supplied preamble snapshot. Main introduction theorems use mthm.
An actual LaTeX .aux file, if available, overrides these derived numbers.
This is a reference utility, not a general TeX interpreter.
"""
from pathlib import Path
import argparse
import json
import re

ROOT = Path(__file__).resolve().parent
ENVIRONMENTS = {
    'theorem': 'Theorem', 'thm': 'Theorem', 'mthm': 'Theorem',
    'lemma': 'Lemma', 'lem': 'Lemma', 'remark': 'Remark', 'rem': 'Remark',
    'prop': 'Proposition', 'proposition': 'Proposition', 'defn': 'Definition',
    'cor': 'Corollary', 'conj': 'Conjecture', 'klem': 'Key Lemma',
}

def derive_references(source):
    """Count active sections and theorem environments, ignoring TeX comments."""
    source = re.sub(r'(?<!\\)%[^\n]*', '', source)
    tokens = re.compile(r'\\appendix|\\(section|subsection|subsubsection)(\*)?\{|'
                        r'\\begin\{(' + '|'.join(ENVIRONMENTS) + r')\}')
    section = sub = subsub = theorem = main = 0
    appendix = False
    refs = {}
    for token in tokens.finditer(source):
        if token.group() == r'\appendix':
            appendix, section = True, 0
            continue
        level, star, env = token.groups()
        if level:
            if star:
                continue
            if level == 'section':
                section += 1
                sub = subsub = theorem = 0
            elif level == 'subsection':
                sub += 1
                subsub = 0
            else:
                subsub += 1
            prefix = chr(64 + section) if appendix else str(section)
            number = prefix + (f'.{sub}' if level != 'section' else '')
            if level == 'subsubsection':
                number += f'.{subsub}'
            # Headings can contain nested braces, e.g. \texorpdfstring.
            end, depth = token.end(), 1
            while depth and end < len(source):
                depth += (source[end] == '{') - (source[end] == '}')
                end += 1
            title = source[token.end():end-1]
            match = re.match(r'\s*(?:\\label\{[^}]+\}\s*)+', source[end:])
            labels = re.findall(r'\\label\{([^}]+)\}', match.group()) if match else []
            kind = 'Appendix' if appendix else 'Section'
        else:
            if env == 'mthm':
                main += 1
                number = str(main)
            else:
                theorem += 1
                prefix = chr(64 + section) if appendix else str(section)
                number = f'{prefix}.{theorem}'
            end = source.find(r'\end{' + env + '}', token.end())
            body = source[token.end():end]
            labels = re.findall(r'\\label\{([^}]+)\}', body)
            heading = re.match(r'\s*\[([^\]]*)\]', body)
            title = heading.group(1) if heading else ''
            kind = ENVIRONMENTS[env]
        for label in labels:
            refs[label] = {'kind': kind, 'number': number, 'title': title,
                           'display': f'{kind} {number}',
                           'number_source': 'derived from source and project numbering convention'}
    return refs

def build(aux=None):
    refs = derive_references((ROOT/'reference/full_paper_snapshot.tex').read_text())
    if aux:
        for label, number in re.findall(r'\\newlabel\{([^}]+)\}\{\{([^}]+)\}', Path(aux).read_text()):
            if label in refs:
                refs[label].update(number=number, display=f"{refs[label]['kind']} {number}",
                                   number_source='LaTeX aux')
    return refs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--aux', type=Path, help='Optional compiled paper .aux file')
    args = parser.parse_args()
    refs = build(args.aux)
    (ROOT/'reference/paper_labels.json').write_text(json.dumps(refs, indent=2, ensure_ascii=False)+'\n')
    rows = ['# Paper reference numbers', '',
            'Numbers below refer to the supplied manuscript snapshot. They are derived using the',
            'section-based shared theorem counter in the existing project preamble.',
            'No compiled `.aux` for this exact draft was supplied. Exact source labels remain',
            'the stable identifiers; rerun with `--aux path/to/main.aux` to use compiled numbers.', '',
            '| Reference | Source label | Title |', '|---|---|---|']
    for label, row in refs.items():
        if row['number'].startswith(('6.', 'B.')) or label in ('thm:intro-zeta34','thm:intro-zeta23','sec:proof_strategy'):
            title = ' '.join(row['title'].split()).replace('|', '\\|')
            rows.append(f"| {row['display']} | `{label}` | {title} |")
    (ROOT/'PAPER_REFERENCES.md').write_text('\n'.join(rows)+'\n')
    print(f'Wrote {len(refs)} label mappings.')

if __name__ == '__main__':
    main()
