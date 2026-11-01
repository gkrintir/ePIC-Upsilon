#!/usr/bin/env python3
"""Check source portability, notebook hygiene, links and preserved baseline hashes."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
EXCLUDE={'.git','.venv','workspaces','work','results','__pycache__','.pytest_cache'}

def main():
    errors=[];counts={'python':0,'notebooks':0,'preserved_analysis_files':0}
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
    for item in provenance['copied_files']:
        p=ROOT/item['path']
        if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:
            errors.append(item['path']+': differs from preserved analysis baseline')
        counts['preserved_analysis_files']+=1
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
