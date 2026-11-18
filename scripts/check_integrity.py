"""Verify distribution checksums immediately after extraction, before rerunning."""
from pathlib import Path
import csv, hashlib, sys
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((root/'SHA256SUMS.csv').open(encoding='utf-8')))
bad=[]
for row in rows:
    p=root/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:bad.append(row['path'])
if bad:
    print('Missing or changed files:',*bad,sep='\n');sys.exit(1)
print('Distribution integrity passed for',len(rows),'files.')
