"""Read-only verification of fixed-value Excel -> six raw analysis CSVs."""
from pathlib import Path
import csv, json, math
import openpyxl

root=Path(__file__).resolve().parents[1]
wb=openpyxl.load_workbook(root/'data/source_snapshots/six_outcomes_fixed_values_20261110.xlsx',read_only=True,data_only=True)
mapping=dict(zip(['CAIT','FAAM_ADL','Sports','Anterior','Posteromedial','Posterolateral'],['CAIT','FAAM-ADL','运动功能','前向平衡','后内侧平衡','后外侧平衡']))
summary=[]; cells=[];components=[]
for key,zh in mapping.items():
    sheet=zh+'_五列汇总'; rows=[]
    for rowno,row in enumerate(wb[sheet].iter_rows(min_row=2,values_only=True),2):
        if isinstance(row[0],str) and row[0].startswith('RPT') and all(isinstance(row[i],(int,float)) for i in [2,3,4]):
            rows.append(row[:5])
            cells.append(dict(outcome=key,sheet=sheet,row=rowno,id=row[0],t=row[1],mean=row[2],sd=row[3],n=row[4],sd_note=row[6] if len(row)>6 else None,scope_note=row[7] if len(row)>7 else None,source_note=row[8] if len(row)>8 else None))
    old=list(csv.DictReader((root/f'data/analysis_inputs/primary/{key}_arms.csv').open(encoding='utf-8')))
    assert len(old)==len(rows),(key,'row count')
    maxdiff=0
    for expected,got in zip(rows,old):
        assert expected[:2]==(got['id'],got['t']),(key,expected,got)
        for i,col in [(2,'mean'),(3,'sd'),(4,'n')]:
            diff=abs(expected[i]-float(got[col]));maxdiff=max(diff,maxdiff)
            assert math.isclose(expected[i],float(got[col]),rel_tol=1e-12,abs_tol=1e-12),(key,col,expected,got)
    summary.append(dict(outcome=key,source_rows=len(rows),numeric_cells=3*len(rows),max_absolute_difference=maxdiff,passed=True))
    sheet=zh+'_干预成分'; cc=list(wb[sheet].iter_rows(values_only=True));head=cc[0][1:]
    for row in cc[1:]:
        if not isinstance(row[0],str):continue
        if any(v is None for v in row[1:len(head)+1]):continue
        record={'outcome':key,'treatment':row[0]}
        for name,value in zip(head,row[1:]):
            assert int(name in row[0].split('+'))==int(value),(key,row[0],name,value)
            record[name]=value
        components.append(record)
out=root/'results/validation';out.mkdir(parents=True,exist_ok=True)
for filename,records in [('source_workbook_checks.csv',summary),('source_cell_map.csv',cells),('source_component_matrix.csv',components)]:
    with (out/filename).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
print('Workbook matches all',sum(s['source_rows'] for s in summary),'raw arm rows; component coding also matches.')
