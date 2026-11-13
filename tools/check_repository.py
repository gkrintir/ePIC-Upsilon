#!/usr/bin/env python3
"""Check source portability, notebook hygiene, links and recorded source hashes."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
EXCLUDE={'.git','.venv','workspaces','work','results','__pycache__','.pytest_cache'}


def check_analysis_snapshot(root, provenance):
    copied = {item['path']: item['sha256'] for item in provenance['copied_files']}
    if len(copied) != len(provenance['copied_files']):
        raise ValueError('copied_files contains duplicate paths')

    records = provenance.get('amended_files', [])
    if not isinstance(records, list):
        raise ValueError('amended_files must be a list')
    amendments = {}
    required_fields = {'path', 'baseline_sha256', 'new_sha256', 'reason'}
    for record in records:
        if not isinstance(record, dict) or set(record) != required_fields:
            raise ValueError('Each amendment requires exactly path, baseline_sha256, new_sha256 and reason')
        path = record['path']
        if not isinstance(path, str) or path not in copied:
            raise ValueError('Unknown amended path: ' + str(path))
        if path in amendments:
            raise ValueError(path + ': duplicate amendment')
        if not isinstance(record['reason'], str) or not record['reason'].strip():
            raise ValueError(path + ': amendment requires a nonempty reason')
        for field in ('baseline_sha256', 'new_sha256'):
            digest = record[field]
            if not isinstance(digest, str) or re.fullmatch(r'[0-9a-f]{64}', digest) is None:
                raise ValueError(path + ': invalid ' + field)
        if record['baseline_sha256'] != copied[path]:
            raise ValueError(path + ': amendment baseline differs from copied_files')
        amendments[path] = record

    counts = {
        'preserved_analysis_files': len(copied) - len(amendments),
        'amended_analysis_files': len(amendments),
    }
    errors = []
    for path, baseline in copied.items():
        expected = baseline
        reference = 'preserved analysis baseline'
        if path in amendments:
            expected = amendments[path]['new_sha256']
            reference = 'recorded analysis amendment'
        source = root / path
        if not source.is_file() or hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            errors.append(path + ': differs from ' + reference)
    return counts, errors


def main():
    errors=[];counts={'python':0,'notebooks':0,'preserved_analysis_files':0,'amended_analysis_files':0}
    sources=[p for p in ROOT.rglob('*') if p.is_file() and not any(x in EXCLUDE for x in p.relative_to(ROOT).parts)]
    for p in sources:
        rel=p.relative_to(ROOT)
        if p.suffix=='.py':
            try:ast.parse(p.read_text(),filename=str(rel),feature_version=(3,9));counts['python']+=1
            except (SyntaxError,UnicodeError) as exc:errors.append(str(rel)+': '+str(exc))
        if p.suffix=='.ipynb':
            counts['notebooks']+=1
            try:
                nb=json.loads(p.read_text())
                for cell in nb['cells']:
                    if cell.get('cell_type')=='code' and (cell.get('outputs') or cell.get('execution_count') is not None):
                        errors.append(str(rel)+': contains execution output');break
            except (ValueError,KeyError) as exc:errors.append(str(rel)+': '+str(exc))
        if p.suffix.lower() in {'.root','.zip','.tar','.pyc'}:errors.append(str(rel)+': generated/binary data in source tree')
        if p.suffix in {'.py','.md','.json','.toml','.ipynb','.sh','.C','.h'}:
            text=p.read_text(errors='replace')
            if re.search(r'\b(?:olp_|ghp_|github_pat_)[A-Za-z0-9_]{25,}',text):errors.append(str(rel)+': credential-shaped token; inspect privately')
    provenance=json.loads((ROOT/'provenance/analysis_snapshot.json').read_text())
    try:
        source_counts, source_errors = check_analysis_snapshot(ROOT, provenance)
        counts.update(source_counts)
        errors.extend(source_errors)
    except ValueError as exc:
        errors.append('provenance/analysis_snapshot.json: ' + str(exc))
    docs=[ROOT/'README.md']+list((ROOT/'docs').glob('*.md'))
    for p in docs:
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',p.read_text()):
            target=target.split('#',1)[0]
            if not target or '://' in target or target.startswith('mailto:'):continue
            if not (p.parent/target).exists():errors.append(str(p.relative_to(ROOT))+': missing link '+target)
    for key,value in counts.items():print(key+': '+str(value))
    if errors:
        for err in errors:print('FAIL '+err,file=sys.stderr)
        return 1
    print('PASS repository checks')
    return 0

if __name__=='__main__':raise SystemExit(main())
