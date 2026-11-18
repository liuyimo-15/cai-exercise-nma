from pathlib import Path
import csv, json, math, shutil

root=Path(__file__).resolve().parents[1]
got=json.loads((root/'work/cinema_20261112/assessment.json').read_text())
ref=json.loads((root/'expected/conditional_grade_assessment.json').read_text())
fields=['grade','original_grade','bias_down','imp_down','het_down','inc_down','merged_down','total_down','reporting','indirectness']
key=lambda r:(r['outcome'],r['comparison'])
a={key(r):r for r in ref['rows']}; b={key(r):r for r in got['rows']}
assert a.keys()==b.keys()
rows=[]
for k in a:
    assert all(a[k][f]==b[k][f] for f in fields), (k,'rating differs')
    for f in ['estimate','lower','upper','pi_lower','pi_upper','bias_score','delta']:
        assert math.isclose(a[k][f],b[k][f],rel_tol=1e-5,abs_tol=1e-7),(k,f)
    row=b[k].copy()
    for f in ['sensitivity']: row[f]=json.dumps(row[f],ensure_ascii=False)
    rows.append(row)
out=root/'results/tables/conditional_grade';out.mkdir(parents=True,exist_ok=True)
shutil.copyfile(root/'work/cinema_20261112/assessment.json',out/'assessment.json')
with (out/'comparison_ratings.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(root/'results/validation/grade_verification.json').write_text(json.dumps({'comparisons':len(rows),'all_ratings_match':True,'status':'Conditional user-specified RoB2/modified-CINeMA scenario, not evidence-validated final GRADE'},indent=2))
print('Conditional grade matches',len(rows),'comparison ratings.')
