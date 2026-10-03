"""Import verified 2026 Georges Bank haddock summaries; retain incomplete status."""
import csv,io,re,sys
from pathlib import Path
root=Path(sys.argv[1])/'analysis/comp_assessments'
cache=root/'source_cache/nefsc_haddock_georges_bank_2026'
assessment='nefsc_haddock_georges_bank_2026'
stock_id='nefsc_haddock_georges_bank'
portal='https://apps-nefsc.fisheries.noaa.gov/saw/sasi_files.php?year=2026&species_id=5&stock_id=1&review_type_id=6&info_type_id=-1&map_type_id=&filename=Georges_Bank_haddock_Update_2026_06_17_100805.801712.pdf'
review='https://www.fisheries.noaa.gov/s3/2026-07/june-2026-management-track-peer-review-panel-report_508_20260710.pdf'
framework='https://www.fisheries.noaa.gov/s3/2023-06/2022-RT-Haddock-GB-EGB-summary-ONLY-6-23-23-nefsc.pdf'
text=(cache/'2026_assessment_text.txt').read_text(encoding='utf-8')
final=(cache/'2026_peer_review_text.txt').read_text(encoding='utf-8')
assert '29,037' in final and '26,615' in final and '2DAR1' in final
inputs,outputs,assumptions=[],[],[]
def values(label):
 line=next(line for line in text.splitlines() if line.startswith(label))
 nums=[float(x.replace(',','')) for x in re.findall(r'\d[\d,]*(?:\.\d+)?',line[len(label):])]
 assert len(nums)==8,(label,nums)
 return nums
for year,value in zip(range(2018,2026),values('Catch for Assessment')):
 inputs.append(dict(assessment_id=assessment,type='catch',measure='total_biomass',basis='biomass',fleet='Combined commercial removals',year=year,year_basis='calendar_year',value=value,unit='t',source_type='official_table',source_reference=portal+'; Table 1; Catch for Assessment',notes='Published modeled total, including US landings/discards and Canadian catch. Native age composition not yet recovered; source total retained without manufacturing catch-at-age.'))
for age in range(1,10):
 inputs.append(dict(assessment_id=assessment,type='M',measure='natural_mortality_at_age',basis='per_year',year=1931,year_basis='calendar_year',age=age,value=.2,unit='per year',source_type='official_table',source_reference=review+'; Georges Bank haddock TOR 3',notes='Fixed M=0.2 for all modeled years and ages; constant vector stored once. Age 9 is 9+.'))
for measure,label,unit,kind in [('SSB','Spawning Stock Biomass','t','biomass'),('Fbar','F¯','per year','mortality'),('recruitment','Recruits (age 1)','thousand fish','recruitment')]:
 for year,value in zip(range(2018,2026),values(label)):
  row=dict(assessment_id=assessment,type=kind,measure=measure,year=year,value=value,unit=unit,source_type='official_table',source_reference=portal+'; Table 1',notes='Accepted WHAM historical estimate; no retrospective adjustment. Rounded published values.')
  if measure=='Fbar': row['age_group']='5–7'
  if measure=='recruitment': row['age']=1
  if year==2025 and measure=='SSB': row.update(lwr=16839,upr=50070,notes=row['notes']+' Published 95% limits; SE not reconstructed.')
  if year==2025 and measure=='Fbar': row.update(lwr=.114,upr=.388,notes=row['notes']+' Published 95% limits; SE not reconstructed.')
  outputs.append(row)
for component,setting,value,notes,source in [
 ('assessment','accepted_model','WHAM 2026 Management Track','Final June 2026 peer review accepts this update as BSIA.',review),
 ('population','modeled_years','1931–2025','Report figure time series; exact native dimensions pending.',portal),
 ('population','ages','1–9; plus group 9+','Nine NAA panels in diagnostic bundle and base framework settings; current native inputs pending.',framework),
 ('population','recruitment_age','1','Published recruitment definition.',portal),
 ('catch','fleet_structure','One combined commercial-removals fleet','US landings/discards and Canadian catches; native composition still needed.',review),
 ('catch','composition_likelihood','Logistic normal','Unchanged 2022 Research Track model approach.',review),
 ('index','active_streams','NEFSC spring BTS; NEFSC fall BTS; DFO BTS','All three require numerical input extraction.',review),
 ('index','excluded_observations','NEFSC spring 2023 and NEFSC fall 2025','Accepted exclusions, not missing values to impute.',review),
 ('index','zero_observations','Treated as missing','Final TOR 3 model description.',review),
 ('index','sampling_time','unknown','Do not borrow timing from the separate Eastern Georges Bank assessment.',review),
 ('M','treatment','Fixed 0.2 for all years and ages','Constant input represented once.',review),
 ('N','process','2DAR1: correlated across age and year','Final TOR 3 confirms unchanged framework process.',review),
 ('assessment','retrospective_adjustment','None','No rho adjustment for stock status or projections.',portal),
 ('biology','projection_weights','Gaussian Markov random field','Projection method; not a replacement for historical supplied weights.',portal),
 ('q','parameter_structure','unknown','Current native configuration not yet available.',review),
 ('F','selectivity_structure','unknown','Current native configuration not yet available.',review)]:
 assumptions.append(dict(assessment_id=assessment,component=component,setting=setting,value=value,notes=notes,source_reference=source))
stock=dict(stock_id=stock_id,charbonneau_id='',authority='NOAA-NEFSC',authority_stock_id='Georges Bank haddock',scientific_name='Melanogrammus aeglefinus',common_name='Georges Bank haddock',area='Georges Bank',region='Northeast US',ocean='North Atlantic',notes='Not found in checked Charbonneau PLOS 2026 catalogue. NEFSC-GARMIII_5Y_Melanogrammus_Aeglefinus is Gulf of Maine haddock and is not this stock. Separate from DFO Eastern Georges Bank assessment.')
record=dict(assessment_id=assessment,stock_id=stock_id,assessment_year=2026,terminal_year=2025,estimate_terminal_year=2025,assessment_type='annual',model_family='WHAM',model_version='2026 management track',framework_year=2022,is_current='TRUE',is_applied='TRUE',assessment_url=portal,framework_url=framework,data_url='https://apps-nefsc.fisheries.noaa.gov/saw/sasi.php',assumptions_status='partial',inputs_status='partial',outputs_status='partial',notes='Accepted run verified against final June 2026 review, including survey exclusions and terminal SSB/F. Only published 2018–2025 totals/summaries and fixed M entered. Full catch-at-age, survey, biology and N/F-at-age inputs/outputs remain unresolved. Source review documents cached sources and searches.')
for filename,new,key in [('stocks.csv',[stock],'stock_id'),('assessments.csv',[record],'assessment_id'),('inputs.csv',inputs,'assessment_id'),('outputs.csv',outputs,'assessment_id'),('assumptions.csv',assumptions,'assessment_id')]:
 path=root/'database'/filename
 content=path.read_bytes().decode('utf-8-sig'); reader=csv.DictReader(io.StringIO(content)); fields=reader.fieldnames; old=list(reader)
 lines=content.splitlines(keepends=True); assert len(lines)==len(old)+1
 with path.open('w',newline='',encoding='utf-8') as f:
  f.writelines([lines[0]]+[line for line,row in zip(lines[1:],old) if row[key]!=new[0][key]])
  csv.DictWriter(f,fields).writerows({field:r.get(field,'') for field in fields} for r in new)
print(f'Imported {len(inputs)} inputs, {len(outputs)} outputs, {len(assumptions)} assumptions; all statuses partial.')
