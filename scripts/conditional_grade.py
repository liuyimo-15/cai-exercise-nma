from pathlib import Path
from collections import defaultdict,Counter
import csv,json,math,hashlib,re
from pypdf import PdfReader

BASE=Path('work/cinema_20261112');OUT=Path('outputs/2026-11-12_CINeMA网状GRADE与临床阈值')
names=dict(zip(['CAIT','FAAM_ADL','Sports','Anterior','Posteromedial','Posterolateral'],['CAIT','FAAM-ADL','运动功能','前向平衡','后内侧平衡','后外侧平衡']))
def read(p):return list(csv.DictReader(open(p,encoding='utf-8-sig')))
def num(x):
 try:return float(x) if str(x) not in ['NA','','None'] else None
 except:return None
def risklevel(score):return 0 if score<1.5 else 1 if score<=2.5 else 2
def precision(est,lo,hi,delta):
 if lo>=0 or hi<=0 or (lo>=-delta and hi<=delta):return 0
 if (est>=0 and lo< -delta) or (est<0 and hi>delta):return 2
 return 1
def certainty(down):return ['高','中','低','极低'][min(3,down)]
judges=['无明显顾虑','有些顾虑','严重顾虑']
threshold_scale=json.load(open(BASE/'threshold_scale.json'))
deltas={'CAIT':3.,'FAAM_ADL':8.,'Sports':threshold_scale['Sports_delta'],'Anterior':.5,'Posteromedial':.5,'Posterolateral':.5}
ratings={(r['outcome'],r['id']):r for r in read(BASE/'ratings.csv')}
for key in ['Anterior','Posteromedial','Posterolateral']:
 ratings[(key,'RPT092_306')]={'outcome':key,'id':'RPT092_306','author':'Reyes 2024 / Foreman 2024（合并输入）','D1':'Some concerns','D2':'Low risk','D3':'Some concerns','D4':'Some concerns','D5':'Some concerns','Overall':'Some concerns','Original_Overall':'Some concerns','supplement':True}
network={r['outcome']:{k:(num(v) if k!='outcome' else v) for k,v in r.items()} for r in read(BASE/'network_summary.csv')}
scores={'Low risk':1,'Some concerns':2,'High risk':3}
allrows=[];contribrows=[];scope=[];side_validation=[]
pdf=PdfReader('data/source_snapshots/R050_primary_analysis_reference.pdf')
for ix,(key,name) in enumerate(names.items()):
 od=BASE/'diagnostics'/key
 contrib=defaultdict(list)
 for c in read(od/'study_contributions.csv'):
  c['contribution']=float(c['contribution']);contrib[c['comparison']].append(c)
  rr=ratings[(key,c['study'])]
  contribrows.append([name,c['comparison'],c['study'],rr['author'],rr['Overall'],rr['Original_Overall'],c['contribution'],bool(rr.get('supplement'))])
 direct={r['comparison']:r for r in read(od/'direct.random.csv')}
 indirect={r['comparison']:r for r in read(od/'indirect.random.csv')}
 local={r['comparison']:r for r in read(od/'compare.random.csv')}
 pairs=read(od/'pairwise_input.csv')
 for row in read(od/'all_estimates.csv'):
  a,b=row['treatment'],row['reference'];cmp=':'.join(sorted([a,b]));delta=deltas[key]
  weights=contrib[cmp];assert abs(sum(x['contribution'] for x in weights)-1)<1e-7,(key,cmp)
  shares=Counter();origshares=Counter();supplement=0
  for c in weights:
   rr=ratings[(key,c['study'])];shares[rr['Overall']]+=c['contribution'];origshares[rr['Original_Overall']]+=c['contribution']
   if rr.get('supplement'):supplement+=c['contribution']
  score=sum(scores[g]*v for g,v in shares.items());bias=risklevel(score);origscore=sum(scores[g]*v for g,v in origshares.items());oldbias=risklevel(origscore)
  est,lo,hi,pl,pu=[float(row[x]) for x in ['estimate','lower','upper','pi_lower','pi_upper']]
  imp=precision(est,lo,hi,delta);pirisk=precision(est,pl,pu,delta);het=pirisk if pirisk>imp else 0
  lp=num(local.get(cmp,{}).get('p'));gp=network[key]['global_p'];inc=0
  if lp is None:inc=2 if gp is None or gp<.05 else 1 if gp<=.10 else 0
  elif lp<=.10:
   de=num(direct.get(cmp,{}).get('TE'));ie=num(indirect.get(cmp,{}).get('TE'))
   inc=2 if de is not None and ie is not None and de*ie<0 and abs(de)>delta and abs(ie)>delta else 1
  merged=max(imp,het,inc);total=bias+merged;grade=certainty(total)
  directpairs=[p for p in pairs if {p['treat1'],p['treat2']}=={a,b}]
  dk=len(directpairs);dn=sum(float(p['n1'])+float(p['n2']) for p in directpairs)
  sensitivity=[]
  variants=([2.,3.,4.] if key=='CAIT' else [6.,8.,10.] if key=='FAAM_ADL' else [.2,.5,delta,.8] if key=='Sports' else [.2,.5,.8])
  for dd in sorted(set(variants)):
   ii=precision(est,lo,hi,dd);pp=precision(est,pl,pu,dd);hh=pp if pp>ii else 0
   sensitivity.append({'delta':dd,'grade':certainty(bias+max(ii,hh,inc))})
  reason=f'贡献加权RoB分数{score:.3f}，低/中/高贡献={shares["Low risk"]:.1%}/{shares["Some concerns"]:.1%}/{shares["High risk"]:.1%}。δ={delta:.3f}；95%CI[{lo:.3f},{hi:.3f}]；95%PI[{pl:.3f},{pu:.3f}]。'
  reason+=f'局部SIDE P={lp:.3f}。' if lp is not None else f'本比较局部SIDE不可检验，按全局P={gp:.3f}暂评。'
  reason+=f'后三域取max({imp},{het},{inc})={merged}，不叠加。'
  flags=[]
  if key in ['Anterior','Posteromedial','Posterolateral']:flags+=['分方向锚定MCID未确立；SMD±0.5为操作性阈值']
  if supplement>1e-8:flags+=[f'合并ID RPT092_306补充暂评贡献{supplement:.1%}；未训练侧及共享队列待核']
  if key!='CAIT':flags+=['RPT023/RPT184潜在共享队列按R050原样保留']
  allrows.append(dict(outcome=key,outcome_zh=name,comparison=f'{a} vs {b}',treatment=a,reference=b,primary=b=='NAT',sm=row['sm'],estimate=est,lower=lo,upper=hi,pi_lower=pl,pi_upper=pu,delta=delta,bias_score=score,within=judges[bias],reporting='不降级（指定口径）',indirectness='不降级（指定口径）',imprecision=judges[imp],heterogeneity=judges[het],incoherence=judges[inc],bias_down=bias,imp_down=imp,het_down=het,inc_down=inc,merged_down=merged,total_down=total,grade=grade,reason=reason,sensitivity=sensitivity,original_grade=certainty(oldbias+merged),rob_low=shares['Low risk'],rob_some=shares['Some concerns'],rob_high=shares['High risk'],supplement_weight=supplement,direct_k=dk,direct_n=int(dn),global_p=gp,local_p=lp,flags='；'.join(flags),status='条件性初评：指定RoB2情景'+('；含补充暂评' if supplement>0 else '')))
 # Verify the PDF SIDE display against reconstructed diagnostics.
 txt=pdf.pages[ix*10+7].extract_text().replace('−','-');labels=re.findall(r'^([A-Z+]+:[A-Z+]+)$',txt,re.M);vals=re.findall(r'^(-?\d+\.\d+) \[(-?\d+\.\d+), (-?\d+\.\d+)\]$',txt,re.M)
 assert len(labels)==len(vals),(key,len(labels),len(vals))
 differences=[]
 for label,v in zip(labels,vals):
  r=local[label];diff=max(abs(float(r[k])-float(z)) for k,z in zip(['TE','lower','upper'],v));differences.append(diff)
 assert max(differences)<.00501,(key,max(differences))
 side_validation.append(dict(outcome=key,contrasts=len(vals),max_difference=max(differences)))
 arms=read(od/'analysed_arms.csv');active={x['id'] for x in arms};inputids={r['id'] for r in read(BASE/(key+'_arms.csv'))}
 scope.append(dict(outcome=name,k=int(network[key]['k']),nodes=int(network[key]['nodes']),active_ids=sorted(active),not_in_pdf_network=sorted(inputids-active),scope_method='按森林图、Q/df及SIDE的数值匹配重建R050分析子集；排除原因需原始R配置确认'))

summary=[]
for key,name in names.items():
 rows=[r for r in allrows if r['outcome']==key];primary=[r for r in rows if r['primary']];summary.append(dict(outcome=key,name=name,sm=rows[0]['sm'],network=network[key],comparisons=len(rows),primary_comparisons=len(primary),grades=dict(Counter(r['grade'] for r in primary)),all_grades=dict(Counter(r['grade'] for r in rows)),threshold=deltas[key],source_scope=next(x for x in scope if x['outcome']==name)))

references=[
 ['C01','Wright CJ, Linens SW, Cain MS.','Establishing the Minimal Clinical Important Difference and Minimal Detectable Change for the Cumberland Ankle Instability Tool.','Arch Phys Med Rehabil. 2017;98(9):1806-1811.','10.1016/j.apmr.2017.01.003','https://pubmed.ncbi.nlm.nih.gov/28137476/','CAI 50人；锚定MCID≥3分；MDC95=3.08。支持CAIT临床阈值3分。'],
 ['C02','Martin RL, Irrgang JJ, Burdett RG, et al.','Evidence of validity for the Foot and Ankle Ability Measure (FAAM).','Foot Ankle Int. 2005;26(11):968-983.','10.1177/107110070502601113','https://pubmed.ncbi.nlm.nih.gov/16309613/','MCID：ADL 8、Sport 9个百分点；源于混合足踝疾病人群约4周，向CAI外推需说明。MDC不是MCID。'],
 ['C03','Carcia CR, Martin RL, Drouin JM.','Validity of the Foot and Ankle Ability Measure in athletes with chronic ankle instability.','J Athl Train. 2008;43(2):179-183.','10.4085/1062-6050-43.2.179','https://pmc.ncbi.nlm.nih.gov/articles/PMC2267323/','支持CAI构念效度；FADI-Sport与FAAM-Sport的8项运动分量表相同。未为CAI单独建立MCID。'],
 ['C04','Munro AG, Herrington LC.','Between-session reliability of the star excursion balance test.','Phys Ther Sport. 2010;11(4):128-132.','10.1016/j.ptsp.2010.07.002','https://pubmed.ncbi.nlm.nih.gov/21055706/','22名健康运动者；归一化SEBT约6%–8%的改变可超出测量误差。属于SDD/MDC，不能当作CAI的分方向锚定MCID。'],
 ['M01','Nikolakopoulou A, Higgins JPT, Papakonstantinou T, et al.','CINeMA: An approach for assessing confidence in the results of a network meta-analysis.','PLoS Med. 2020;17(4):e1003082.','10.1371/journal.pmed.1003082','https://pmc.ncbi.nlm.nih.gov/articles/PMC7122720/','六域、贡献矩阵、临床等效区；对相关顾虑避免重复降级。研究内偏倚阈值为本次操作规则，非软件唯一官方算法。'],
 ['M02','Cochrane Handbook, Chapter 15.','Interpreting results and drawing conclusions.','Current online version; accessed 2026-11-12.','','https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-15','SMD 0.2/0.5/0.8为一般小/中/大效应解释规则，不能自动视为患者MCID；可通过参考SD重表达原量表单位。'],
 ['M03','netmeta development team.','netmeta reference manual: netcontrib, netsplit and netmeta.','netmeta 3.3-1 used in this analysis.','','https://cran.rstudio.com/web/packages/netmeta/refman/netmeta.html','使用shortestpath、study=TRUE、纳入多臂相关性的研究贡献权重；SIDE用Back-calculation，固定与主分析相同τ。'],
 ['L01','Reyes MC et al.','Cross-education effects of balance training in individuals with chronic ankle instability.','J Bodyw Mov Ther. 2024.','10.1016/j.jbmt.2024.07.042','https://doi.org/10.1016/j.jbmt.2024.07.042','原PDF Methods/Results：未训练/未受累侧；30人中11人无后测，采用多重插补；属于较大项目的一部分。合并ID补评及队列复核。'],
 ['L02','Foreman NN et al.','Cross-Education Effects of Hip Strengthening in Individuals with Chronic Ankle Instability.','Transl J Am Coll Sports Med. 2024.','10.1249/TJX.0000000000000252','https://doi.org/10.1249/TJX.0000000000000252','原PDF Methods/Results：未训练侧；封闭信封但隐藏细节不足，评估者未盲；30人中10人缺失、MI；较大项目的次要分析。'],
 ['L03','Suttmiller AMB et al.','Comparing the Effects of Progressive Balance and Hip Strengthening Rehabilitation in Individuals With Chronic Ankle Instability.','J Sport Rehabil. 2025.','10.1123/jsr.2024-0049','https://doi.org/10.1123/jsr.2024-0049','原PDF Methods/flow：BAL/HIP/CON各15人，共45；用于核对合并ID RPT092_306是否共享队列及患侧/非患侧差异。']
]
thresholds=[
 dict(outcome='CAIT',scale='0–30分，越高越好',clinical='3分',clinical_value=3,analysis_scale='MD',delta=3,equivalence='[-3, +3]',type='CAI锚定MCID',refs='C01',recommendation='建议主阈值3分；敏感性2/4分',limits='个体变化MCID借用作组间等效区边界为综述假设；不能等同个体应答率。'),
 dict(outcome='FAAM-ADL',scale='0–100%，越高越好',clinical='8个百分点',clinical_value=8,analysis_scale='MD',delta=8,equivalence='[-8, +8]',type='锚定MCID，混合足踝疾病人群',refs='C02；C03',recommendation='建议主阈值8个百分点；敏感性6/10',limits='来自约4周混合足踝疾病；不是CAI专属MCID；0–84原始分须先换成百分制。'),
 dict(outcome='运动功能',scale='FAAM-Sport/FADI-Sport 0–100%，越高越好',clinical='9个百分点',clinical_value=9,analysis_scale='SMD (Hedges g)',delta=deltas['Sports'],equivalence=f'[-{deltas["Sports"]:.3f}, +{deltas["Sports"]:.3f}]',type='原量表MCID＋本数据参考SD换算',refs='C02；C03；M02',recommendation=f'主δ=9/12.2726={deltas["Sports"]:.3f}；另以SMD 0.2/0.5/0.8敏感性',limits='0.733不是文献直接报告的MCID；参考SD是19份报告试验内合并变化SD的中位数。仅计算参考SD时将RPT184由0–32转为0–100，未改NMA。'),
 *[dict(outcome=n,scale='SEBT/YBT，原文含归一化距离或其他单位；越高越好',clinical='尚未找到可直接采用的CAI分方向锚定MCID',clinical_value=None,analysis_scale='SMD (Hedges g)',delta=.5,equivalence='[-0.5, +0.5]',type='操作性阈值，非验证MCID',refs='C04；M02',recommendation='暂用SMD 0.5；必须报告0.2/0.8阈值敏感性',limits='0.5仅为一般中等效应约定。健康人SEBT 6%–8%是测量误差界限，不是CAI临床重要差异；不可用4cm伤害风险截点代替MID。') for n in ['前向平衡','后内侧平衡','后外侧平衡']]
]
rules=[
 ['评价单位','每个结局×每个治疗节点对比；主表为48个对NAT比较，六个结局分表共219个全部网络比较。不能用每结局单一等级替代逐比较GRADE。'],
 ['等级性质','最新RoB2来自用户按比例设置的情景，不是新增原文证据重评。本次分级是该情景下的条件性初评，未声称已完成双人独立评价或CINeMA网站认证。'],
 ['起点与总公式','RCT证据从高起评；总降级=研究内偏倚降级+max(不精确、异质性、不一致性降级)。报告偏倚与间接性按要求固定0。总降级0/1/2/≥3对应高/中/低/极低。'],
 ['报告偏倚与间接性','固定不降级为用户指定假设，不等于已证明无发表偏倚、无间接性；漏斗图或局部检验不显著也不能证明无偏倚。'],
 ['研究内偏倚贡献','用netmeta 3.3-1 netcontrib shortestpath且考虑多臂相关性，按真实研究贡献加权，而非按RoB2图的报告数量百分比推断。每个比较权重和为1。'],
 ['偏倚操作规则','每报告×结局用图对应的Overall：Low=1、Some=2、High=3。贡献加权分数<1.5降0；1.5–2.5降1；>2.5降2。这是供复核的操作规则，不是CINeMA唯一官方阈值。'],
 ['不精确判据','本次采用CINeMA相对宽松的方向性判读：95%CI完全在0同侧，或完全在[-δ,+δ]内则0；跨0但不越过点估计相反方向δ为1；越过相反方向δ为2。不因单纯跨临床有益阈值、而全部仍在0同侧就自动降级。'],
 ['异质性判据','以相同δ比较95%CI和95%PI；若PI使方向性判断从无顾虑转为有些/严重顾虑，分别1/2；若CI已显示同样的不确定性，则不重复记异质性。报告τ²、Q、I²供解释，不用I²>50%机械降级。'],
 ['预测区间','PI=估计±1.96×sqrt(SE²+τ²)，与提供的Gaussian REML分析框架一致；未计入τ²估计和模型选择不确定性，不代表未来个体变化范围。'],
 ['不一致性判据','局部SIDE P>0.10暂不降级；≤0.10再看直接/间接效应临床方向与区间。局部不可检验时按同τ的全局设计-治疗Wald P：>0.10为0，0.05–0.10为1，<0.05为2；不是以“不显著”等同完全一致。'],
 ['主分析重建','MD用Gaussian REML；SMD用Hedges g、试验内共同合并SD、N作方差分母，多臂相关0.5。6个网络森林图、Q/df、77个SIDE显示均匹配PDF至两位小数；netmeta辅助拟合与重建拟合误差<1e-6。'],
 ['重建范围限制','未取得原始R050配置/模型对象；根据已给PDF数值匹配分析子集。重建可复现并不确认排除原因，须按原R配置复核。未改动源PDF或源数据。'],
 ['合并ID补评','RPT092_306仅在三个平衡网络；图中无该ID。暂补D1 Some、D2 Low、D3 Some、D4 Some、D5 Some，Overall Some。依据两文原文随机化、非盲、约1/3缺失及多重插补；该补评仍需独立复核。'],
 ['队列与测试侧','RPT092_306的未训练侧结果与RPT141等可能来自同一45人项目；RPT023/RPT184可能共享主动臂。R050中按原样保留，不能视为均为独立RCT。表中列出贡献权重；解决共享队列和患侧/非患侧后需重跑，当前等级不能作无条件发表定稿。'],
 ['原RoB2敏感性','以原始2026-11-10 Overall替代用户修订Overall、其余统计量和阈值不变，重算GRADE供对照；不覆盖用户选择的主情景。'],
 ['阈值敏感性','CAIT δ2/3/4，FAAM-ADL δ6/8/10；运动功能δ0.2/0.5/0.733/0.8；三方向平衡δ0.2/0.5/0.8。阈值为本次事后选定，需在报告中说明，不能冒充预注册阈值。'],
 ['成分网状边界','这些是传统治疗节点NMA比较的CINeMA分级；不能直接套用于CNMA成分系数、交互作用项或P-score排序。CNMA需另评模型加和性、选模不确定性与成分专属贡献。'],
 ['临床解释','CAIT 3、FAAM-ADL 8、Sport 9为原量表临床锚；把个体变化MCID用作组间临床等效区是额外假设。三方向0.5 SMD不是已验证CAI MCID。']
]
review=[
 ['优先1','RoB2来源','本次输入是按比例重新选位的RoB2情景。','回到原文和信号问题形成真实RoB2共识，再形成正式GRADE。','所有比较','目前均为条件性初评'],
 ['优先1','重复队列和测试侧','RPT092_306含未训练侧；RPT141为患侧。两文共同研究人员、样本及项目结构提示共享队列。','核对试验身份/原始注册及侧别；避免当独立试验重复纳入。','三个平衡网络','未按要求固定0的间接性栏抹去此问题'],
 ['优先1','共享干预臂','RPT023/RPT184主动平衡臂数值重合。','依据队列关联选择报告或正确处理共享臂相关性。','FAAM-ADL、运动功能、三个平衡网络','R050按原样保留'],
 ['优先2','缺少原R050配置','数值匹配显示部分报告未进入PDF分析，原因未写在PDF。','核对原始分析代码及筛选日志。','CAIT RPT284；三个平衡RPT008/128/284/289；后内侧另RPT030','子集由全套数值交叉验证重建'],
 ['优先2','平衡MCID','无可直接采用的CAI分方向锚定MCID。','用操作性SMD阈值并展示0.2/0.5/0.8敏感性，不声称0.5为验证MCID。','三个平衡网络','阈值支持限于一般效应大小规则'],
 ['优先2','模型层次','R050同时报告NMA与CNMA。','不要把节点NMA的等级复制给成分、交互项或治疗排名。','所有结局','当前分级覆盖219个节点比较']
]
rob_input=[]
for key,name in names.items():
 active=next(x for x in scope if x['outcome']==name)['active_ids']
 for id in active:
  r=ratings[(key,id)];rob_input.append([name,id,r['author'],r['D1'],r['D2'],r['D3'],r['D4'],r['D5'],r['Overall'],r['Original_Overall'],'原文补充暂评' if r.get('supplement') else '最新RoB2图对应用户比例情景'])
data=dict(rows=allrows,summary=summary,thresholds=thresholds,references=references,rules=rules,review=review,contributions=contribrows,rob_input=rob_input,scope=scope,side_validation=side_validation,threshold_scale=threshold_scale)
(BASE/'assessment.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
for s in summary:print(s['name'],s['network']['k'],'NAT',s['grades'],'ALL',s['all_grades'])
print('N comparisons',len(allrows),'primary',sum(x['primary'] for x in allrows),'contrib rows',len(contribrows),'SIDE matches',sum(x['contrasts'] for x in side_validation),'overall',dict(Counter(x['grade'] for x in allrows)))
