// Profile-only AA normalization: synthetic scores, no network or source-data changes.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
function element(tag = 'div') {
  return { tag, attributes: {}, children: [], textContent: '', innerHTML: '', hidden: false, value: '',
    dataset: {}, classList: { toggle() {} }, handlers: {},
    addEventListener(type, handler) { this.handlers[type] = handler; },
    setAttribute(key, value) { this.attributes[key] = value; },
    append(child) { this.children.push(child); },
    cloneNode() { return Object.assign(element(this.tag), { attributes: { ...this.attributes }, innerHTML: this.innerHTML }); },
  };
}
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, element());
  return nodes.get(selector);
}
function serialize(el) {
  return `<${el.tag} ${Object.entries(el.attributes).map(([k, v]) => `${k}="${v}"`).join(' ')}>${el.textContent}${el.innerHTML}${el.children.map(serialize).join('')}</${el.tag}>`;
}
let sharedURL;
const context = vm.createContext({ console, URLSearchParams, Blob,
  document: { querySelector: node, querySelectorAll: () => [], addEventListener() {}, createElementNS: (_, tag) => element(tag) },
  location: { search: '', pathname: '/', hash: '' },
  history: { replaceState: (_, __, url) => { sharedURL = url; } },
  XMLSerializer: class { serializeToString(el) { return serialize(el); } },
  proteinTairIds: () => [],
});
const source = fs.readFileSync(path.join(root, 'js', 'app.js'), 'utf8');
vm.runInContext(source.replace(/\binit\(\);\s*$/, ''), context);
const run = text => vm.runInContext(text, context);
const plain = value => JSON.parse(JSON.stringify(value));
const close = (a, b) => assert.ok(Math.abs(a - b) < 1e-10, `${a} should equal ${b}`);
assert.equal(run('state.profileNormalize'), false, 'Normalization is off by default');
run(`
function fixture(id, pocket, score, extra = {}) {
  return { uniprot_id:id, protein:'AF-'+id+'-F1-model_v4', pocket, rank:Number(pocket),
    vina_affinity:score, sfct_score:score, vina_sfct_combined:score, vina_sfct_combined_50:score,
    vina_status:'success', sfct_status:'success', probability:.9, mean_pocket_plddt:95,
    residue_ids:'A_1 A_2', center_x:1, center_y:2, center_z:3, ...extra };
}
function resetFixture() {
  state.rawByAA.clear(); state.rankingCache.clear();
  Object.assign(state, { aa:'ALA', metric:'vina_affinity', p2rank:.7, plddt:90,
    profileValue:'raw_vina', profileNormalize:true, profileOrder:'score_asc', profilePocket:null, profilePocketAnchor:null,
    selectedProtein:'P1', search:'', top:100, aaRank:20, minDelta:null, maxCompetitors:19,
    requireLPreference:false, pocketMode:'best' });
}
resetFixture();
state.rawByAA.set('ALA', [fixture('P1','1',-10), fixture('P1','2',-9),
  fixture('P2','1',-30,{probability:.1}), fixture('P2','2',-8),
  fixture('P3','1',-100,{mean_pocket_plddt:40}), fixture('P3','2',-6),
  fixture('FAILED','1',-200,{vina_status:'failed'}), fixture('NONFINITE','1',NaN)]);
state.rawByAA.set('GLY', [fixture('P1','1',10),fixture('P2','2',12),fixture('P3','2',14)]);
state.rawByAA.set('LEU', [fixture('P1','1',4),fixture('P2','2',2),fixture('P3','2',0)]);
state.rawByAA.set('DALA', [fixture('P1','1',-1000),fixture('P2','1',1000)]);
const beforeComparison = JSON.stringify(getComparison('P1'));
const beforeRanking = JSON.stringify(getRanking());
const beforeTable = profileTableToDelimited(getProteinProfile('P1'), '\\t');
`);

// QC first, one best pocket per protein; each AA has its own mean and population SD.
assert.equal(run('Object.keys(PROFILE_VALUES).length'), 6);
assert.equal(run('PROFILE_VALUES.aa_zscore'), undefined, 'Z is a checkbox, not a value option');
assert.deepEqual(plain(run(`getProfileScoreSummary('ALA','vina_affinity')`)), { count:3, mean:-8, sd:Math.sqrt(8/3), top5Score:-10 });
close(run(`getProfilePlotData('P1').find(e=>e.code==='ALA').plotValue`), -Math.sqrt(1.5));
close(run(`getProfilePlotData('P1').find(e=>e.code==='GLY').plotValue`), -Math.sqrt(1.5));
close(run(`getProfilePlotData('P1').find(e=>e.code==='LEU').plotValue`), Math.sqrt(1.5));
assert.equal(run(`getProfilePlotData('P1').some(e=>e.code.startsWith('D'))`), false);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity') === getProfileScoreSummary('ALA','vina_affinity')`), true);

// Table filters do not truncate the reference population, even if no candidates remain.
run(`const baselinePlot = JSON.stringify(getProfilePlotData('P1'));
  Object.assign(state, {search:'NO_MATCH',top:1,aaRank:1,minDelta:100,maxCompetitors:0,requireLPreference:true,pocketMode:'all'});`);
assert.equal(run('filterRows().length'), 0);
assert.equal(run(`JSON.stringify(getProfilePlotData('P1')) === baselinePlot`), true);
assert.equal(run(`JSON.stringify(getComparison('P1')) === beforeComparison`), true);
assert.equal(run('JSON.stringify(getRanking()) === beforeRanking'), true);
assert.equal(run(`profileTableToDelimited(getProteinProfile('P1'), '\\t') === beforeTable`), true);

// An inspected pocket uses the unchanged best-pocket distribution, never all-pocket scores.
close(run(`getProfilePlotData('P1','2')[0].plotValue`), -1/Math.sqrt(8/3));
assert.equal(run(`getProfilePlotData('P1','2')[0].normalization.mean`), -8);
close(run(`getProfilePlotData('P1','2')[0].top5Z`), -Math.sqrt(1.5));
assert.equal(run(`getProfilePlotData('P1','2').length`), 1, 'No substitution for absent same-pocket AA scores');

// Changes in QC/source invalidate summaries; an existing population can reuse its summary.
run(`const oldSummary=getProfileScoreSummary('ALA','vina_affinity'); state.p2rank=0;`);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').mean`), (-10-30-6)/3);
run('state.p2rank=.7; state.plddt=0;');
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').mean`), (-10-8-100)/3);
run('state.plddt=90;');
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity') === oldSummary`), true);
run(`state.rawByAA.set('ALA',[fixture('P1','1',-20),fixture('P2','1',-10)]);`);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').mean`), -15);

// The checkbox follows all four plot scores, independently of Explorer's Ranking score.
run(`resetFixture(); state.rawByAA.set('ALA',[
  fixture('P1','1',-10,{sfct_score:2,vina_sfct_combined:4,vina_sfct_combined_50:-3}),
  fixture('P2','1',-8,{sfct_score:4,vina_sfct_combined:5,vina_sfct_combined_50:-2}),
  fixture('P3','1',-6,{sfct_score:8,vina_sfct_combined:10,vina_sfct_combined_50:9})]);`);
for (const [metric, values] of Object.entries({ vina_affinity:[-10,-8,-6], sfct_score:[2,4,8],
  vina_sfct_combined:[4,5,10], vina_sfct_combined_50:[-3,-2,9] })) {
  context.testMetric=metric;
  run(`state.profileValue=Object.keys(PROFILE_VALUES).find(key=>!PROFILE_VALUES[key].percentile && PROFILE_VALUES[key].metric===testMetric);
    state.metric=testMetric==='sfct_score'?'vina_affinity':'sfct_score';`);
  const mean=values.reduce((a,b)=>a+b)/3, sd=Math.sqrt(values.reduce((s,v)=>s+(v-mean)**2,0)/3);
  close(run(`getProfilePlotData('P1')[0].plotValue`), (values[0]-mean)/sd);
  close(run(`getProfilePlotData('P1')[0].top5Z`), (values[0]-mean)/sd);
  assert.equal(run('getProfileValueConfig().metric'), metric);
  assert.ok(run('getProfileValueConfig().label').includes(run('METRICS[testMetric].label')));
  run('state.profileNormalize=false;');
  close(run(`getProfilePlotData('P1')[0].plotValue`), values[0]);
  run('state.profileNormalize=true;');
}
run(`state.rawByAA.set('ALA',[fixture('VINA_ONLY','1',-10,{sfct_status:'failed'}),
  fixture('SFCT_ONLY','1',-20,{vina_status:'failed'}),fixture('BOTH','1',-5)]);`);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').mean`), -7.5);
assert.equal(run(`getProfileScoreSummary('ALA','sfct_score').mean`), -12.5);

// Empirical boundaries match the existing Top 5% tier, not a common Gaussian Z cutoff.
run(`resetFixture();
  state.rawByAA.set('ALA',Array.from({length:100},(_,i)=>fixture('P'+String(i).padStart(3,'0'),'1',i)));
  state.rawByAA.set('GLY',Array.from({length:100},(_,i)=>fixture('P'+String(i).padStart(3,'0'),'1',i*i)));
  const empiricalProfile=getProfilePlotData('P099');
`);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').top5Score`), 4);
assert.equal(run(`getProfileScoreSummary('GLY','vina_affinity').top5Score`), 16);
assert.equal(run(`getRanking().filter(row=>row.proteome_percentile<=5).length`), 5);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').top5Score === getDistributionData().cutoffs.find(c=>c.tier===5).value`), true);
for (const aa of ['ALA','GLY']) {
  context.testAA=aa;
  const summary=run(`getProfileScoreSummary(testAA,'vina_affinity')`);
  close(run('empiricalProfile.find(e=>e.code===testAA).top5Z'), (summary.top5Score-summary.mean)/summary.sd);
}
assert.notEqual(run('empiricalProfile[0].top5Z'),run('empiricalProfile[1].top5Z'));
const empiricalSVG=run('profileChartSVG(empiricalProfile)');
const segments = svg => [...svg.matchAll(/<line class="profile-top5-cutoff" data-aa="([A-Z]+)" x1="([^"]+)" x2="([^"]+)" y1="([^"]+)" y2="([^"]+)"/g)];
const cutoffLabel = svg => [...svg.matchAll(/<text class="profile-top5-label" data-aa="([A-Z]+)" text-anchor="([^"]+)" x="([^"]+)" y="([^"]+)">Top 5%<\/text>/g)];
assert.equal(segments(empiricalSVG).length,2);
assert.equal(cutoffLabel(empiricalSVG).length,1);
assert.equal(cutoffLabel(empiricalSVG)[0][1],segments(empiricalSVG).at(-1)[1]);
close(Number(cutoffLabel(empiricalSVG)[0][4]),Number(segments(empiricalSVG).at(-1)[4])-6);
assert.notEqual(segments(empiricalSVG)[0][4], segments(empiricalSVG)[1][4]);
for (const [, , x1, x2, y1, y2] of segments(empiricalSVG)) {
  assert.ok(Number(x1)>=48 && Number(x2)<=808 && Number(x2)>Number(x1));
  assert.ok(Number(y1)>18 && Number(y1)<243, 'Cutoff included in axis extent, even when all displayed scores are poor');
  assert.equal(y1,y2);
}
run(`state.rawByAA.set('ALA',Array.from({length:100},(_,i)=>fixture('P'+String(i).padStart(3,'0'),'1',i<10?0:i)));`);
assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').top5Score`),0);
assert.equal(run(`getRanking().filter(row=>row.proteome_percentile<=5).length`),5,'Ties do not expand rank-based tiers');
assert.equal(run(`getProfilePlotData('P005').find(e=>e.code==='ALA').proteome_percentile > 5`),true);
assert.match(run(`profileChartSVG(getProfilePlotData('P005'))`),/Boundary ties follow rank-based tier colors/);

// Missing, single-protein and constant populations cannot turn into zero Z or invalid SVG.
run(`resetFixture();
  state.rawByAA.set('ALA',[fixture('P1','1',1)]);
  state.rawByAA.set('GLY',[fixture('P1','1',.1),fixture('P2','1',.1),fixture('P3','1',.1)]);
  state.rawByAA.set('LEU',[fixture('P1','1',NaN),fixture('P2','1',-6)]);`);
assert.equal(run(`getProfileScoreSummary('EMPTY','vina_affinity').count`), 0);
assert.ok(Number.isNaN(run(`getProfileScoreSummary('EMPTY','vina_affinity').mean`)));
assert.equal(run(`getProfileScoreSummary('GLY','vina_affinity').sd`), 0);
assert.equal(run(`getProfilePlotData('P1').every(e=>Number.isNaN(e.plotValue))`), true);
assert.equal(run(`getProfilePlotData('P1').every(e=>Number.isNaN(e.top5Z))`), true);
assert.match(run(`profileChartSVG(getProfilePlotData('P1'))`), /No AA-normalized Z-scores available/);
run(`state.rawByAA.set('LEU',[fixture('P1','1',-10),fixture('P2','1',-6)]);`);
let svg=run(`profileChartSVG(getProfilePlotData('P1'))`);
assert.doesNotMatch(svg, /NaN|Infinity|meets cutoff/);
assert.match(svg, /<title>Z = 0<\/title>/);
assert.doesNotMatch(svg, /AA-specific mean/);
assert.match(svg, /AA-normalized Z-score/);
assert.match(svg, /2 reference proteins/);
assert.match(svg, /population SD 2.000/);
assert.equal((svg.match(/<circle /g)||[]).length, 1);
assert.equal(segments(svg).length, 1, 'No cutoff segment for missing/undefined Z');
assert.equal(cutoffLabel(svg)[0][2],'start','Single-AA label stays inside the plot');
assert.match(run(`getProfileNormalizationNote(getProfilePlotData('P1'))`), /Missing\/unavailable Z: ALA/);

// Score order is ascending Z, not raw score; canonical order and tier colors are preserved.
run(`state.rawByAA.set('ALA',[fixture('P1','1',-1),fixture('P2','1',-3)]);
  state.rawByAA.set('GLY',[fixture('P1','1',10),fixture('P2','1',20)]);
  state.rawByAA.set('LEU',[fixture('P1','1',-8),fixture('P2','1',-10),fixture('P3','1',-6)]);`);
const labels = svg => [...svg.matchAll(/<text class="profile-label"[^>]*>([^<]+)<\/text>/g)].map(match=>match[1]);
assert.deepEqual(labels(run(`profileChartSVG(getProfilePlotData('P1'))`)), ['GLY','LEU','ALA']);
assert.deepEqual(segments(run(`profileChartSVG(getProfilePlotData('P1'))`)).map(m=>m[1]),['GLY','LEU','ALA']);
assert.equal(cutoffLabel(run(`profileChartSVG(getProfilePlotData('P1'))`))[0][1],'ALA');
assert.match(run(`profileChartSVG(getProfilePlotData('P1'))`), /profile-dot tier-top1/);
run('state.profileOrder="aa";');
assert.deepEqual(labels(run(`profileChartSVG(getProfilePlotData('P1'))`)), ['ALA','GLY','LEU']);
assert.deepEqual(segments(run(`profileChartSVG(getProfilePlotData('P1'))`)).map(m=>m[1]),['ALA','GLY','LEU']);
assert.equal(cutoffLabel(run(`profileChartSVG(getProfilePlotData('P1'))`))[0][1],'LEU');

// Existing six values, raw table, selectivity summaries and the Vina axis title are unaffected.
for (const mode of ['raw_vina','raw_sfct','raw_combined','raw_combined_50','percentile','percentile_50']) {
  context.mode=mode;
  run('state.profileValue=mode; state.profileNormalize=false;');
  assert.equal(run(`getProfilePlotData('P1').every(e=>e.plotValue===(PROFILE_VALUES[mode].percentile?e.proteome_percentile:e[PROFILE_VALUES[mode].metric]))`), true);
  assert.equal(run(`getProfileNormalizationNote(getProfilePlotData('P1'))`), '');
  assert.doesNotMatch(run(`profileChartSVG(getProfilePlotData('P1'))`), /AA-specific mean|reference proteins|profile-top5-cutoff|profile-top5-label/);
}
run('state.profileValue="raw_vina";');
assert.match(run(`profileChartSVG(getProfilePlotData('P1'))`), /Vina score \(predicted affinity, kcal\/mol\)/);
run('renderPockets=()=>{}; renderProtein();');
const tableBefore=node('#profile-aa-body').innerHTML, selectivityBefore=node('#profile-z').textContent;
assert.equal(node('#profile-normalization-note').hidden, true);
assert.equal(node('#profile-top5-key').hidden, true);
run('state.profileNormalize=true; renderProtein();');
assert.equal(node('#profile-normalization-note').hidden, false);
assert.equal(node('#profile-top5-key').hidden, false);
assert.equal(node('#profile-threshold-key').innerHTML, '<i></i>Z = 0');
assert.match(node('#profile-normalization-note').textContent, /AutoDock Vina \(selected Value type\)/);
assert.equal(node('#profile-aa-body').innerHTML, tableBefore);
assert.equal(node('#profile-z').textContent, selectivityBefore);

// Both controls respond to changes; percentile selection clears/disables the checkbox.
for (const name of ['bindProteinStyleEvents','bindProteinColorEvents','bindPocketSticksEvents',
  'bindLigandViewerEvents','bindPocketCloudEvents','bindGeneDescriptionEvents','bindAnalysisEvents',
  'bindControlQCEvents','bindGOEvents']) context[name]=()=>{};
run('bindProfileProteinSearchEvents=()=>{}; bindEvents();');
const changeValue = value => node('#profile-value-select').handlers.change({target:{value}});
const changeNormalize = checked => node('#profile-normalize').handlers.change({target:{checked}});
changeNormalize(false);
assert.equal(run('state.profileNormalize'),false);
assert.equal(node('#profile-normalize').checked,false);
assert.equal(node('#profile-top5-key').hidden,true);
changeNormalize(true);
assert.equal(node('#profile-normalize').checked,true);
for (const value of ['percentile','percentile_50']) {
  changeValue(value);
  assert.equal(run('state.profileNormalize'),false);
  assert.equal(node('#profile-normalize').checked,false);
  assert.equal(node('#profile-normalize').disabled,true);
  assert.match(node('#profile-normalize-control').title,/unavailable/);
  assert.equal(node('#profile-top5-key').hidden,true);
  assert.doesNotMatch(node('#profile-chart').innerHTML,/profile-top5-cutoff/);
  changeNormalize(true);
  assert.equal(run('state.profileNormalize'),false,'Cannot normalize percentiles even via a synthetic event');
  changeValue('raw_vina');
  assert.equal(node('#profile-normalize').disabled,false);
  assert.equal(node('#profile-normalize').checked,false,'Returning to raw does not silently restore normalization');
  changeNormalize(true);
}
changeValue('raw_sfct');
assert.equal(node('#profile-normalize').checked,true,'Changing raw scores preserves the checkbox');
assert.match(node('#profile-normalization-note').textContent,/OnionNet-SFCT/);
changeValue('raw_vina');

// URL state and downloaded SVG retain the plot metric even when Explorer uses a different score.
run('state.metric="sfct_score"; updateURL();');
assert.match(sharedURL, /profile_value=raw_vina/);
assert.match(sharedURL, /profile_z=1/);
assert.match(sharedURL, /metric=sfct_score/);
context.location.search=sharedURL.slice(1);
run('state.profileValue="raw_combined_50"; state.profileNormalize=false; state.metric="vina_affinity"; readURLState();');
assert.equal(run('state.profileValue'), 'raw_vina');
assert.equal(run('state.profileNormalize'), true);
assert.equal(run('state.metric'), 'sfct_score');
assert.equal(run('getProfileValueConfig().metric'), 'vina_affinity');
const chart=element('svg');
chart.innerHTML=run(`profileChartSVG(getProfilePlotData('P1'))`).replace(/^<svg[^>]*>|<\/svg>$/g,'');
nodes.set('#profile-chart svg',chart);
const exported=run('buildProfilePlotExportSVG();');
assert.match(exported, /AA-normalized Z-score · AutoDock Vina/);
assert.match(exported, /<desc>.*P2Rank.*Not a binding probability/);
assert.match(exported, /<title>Z = 0<\/title>/);
assert.match(exported, /AA-specific Top 5% cutoff/);
assert.match(exported, /\.profile-top5-cutoff \{ stroke: #8e65a5;/);
assert.equal(segments(exported).length,3);
assert.equal(cutoffLabel(exported).length,1,'In-plot cutoff label is also exported');
assert.doesNotMatch(exported, /NaN|Infinity/);
run('let lastDownload; downloadBlob=(filename,blob)=>{lastDownload={filename,blob};};');
Promise.resolve(run('downloadProfilePlot("svg")')).then(() => {
  assert.match(run('lastDownload.filename'), /aa_zscore_Vina_best_pockets_20aa_profile\.svg$/);
  const html=fs.readFileSync(path.join(root,'index.html'),'utf8');
  assert.doesNotMatch(html, /<option value="aa_zscore">/);
  assert.equal((html.match(/type="checkbox" id="profile-normalize"/g)||[]).length,1);
  assert.ok(html.indexOf('id="profile-normalize-control"') < html.indexOf('class="plot-value-control"'), 'Checkbox is before Value type');
  assert.match(html, /id="profile-normalization-note"[^>]*hidden/);
  // Legacy links restore the originally intended normalized score, now as a base + checkbox.
  for (const [metric,value] of Object.entries({vina_affinity:'raw_vina',sfct_score:'raw_sfct',
    vina_sfct_combined:'raw_combined',vina_sfct_combined_50:'raw_combined_50'})) {
    context.location.search=`?metric=${metric}&profile_value=aa_zscore`;
    run('readURLState(); syncProfileValueControls(); updateURL();');
    assert.equal(run('state.profileValue'),value);
    assert.equal(node('#profile-normalize').checked,true);
    assert.equal(run('getProfileValueConfig().metric'),metric);
    assert.match(sharedURL,/profile_z=1/);
    assert.doesNotMatch(sharedURL,/profile_value=aa_zscore/);
  }
  context.location.search='?profile_value=percentile&profile_z=1';
  run('readURLState(); syncProfileValueControls(); updateURL();');
  assert.equal(node('#profile-normalize').disabled,true);
  assert.equal(node('#profile-normalize').checked,false);
  assert.doesNotMatch(sharedURL,/profile_z/);
  context.location.search='?profile_value=raw_vina';
  run('readURLState(); syncProfileValueControls();');
  assert.equal(node('#profile-normalize').disabled,false);
  assert.equal(node('#profile-normalize').checked,false);
  run('renderExplorer=()=>{}; state.profileNormalize=true;');
  node('#reset-filters').handlers.click();
  assert.equal(run('state.profileNormalize'),false);
  assert.equal(run('state.profileValue'),'raw_combined_50');
  assert.equal(node('#profile-normalize').checked,false);
  console.log('Profile AA-normalized Z-score checks passed.');
}).catch(error=>{console.error(error);process.exitCode=1;});
