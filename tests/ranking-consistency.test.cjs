// Shared QC ranking, inspected-pocket identity, histogram tiers and versioned-model checks.
// No network requests or source-data writes.
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const assert = require('node:assert/strict');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
const node = selector => {
  if (!nodes.has(selector)) nodes.set(selector, { value: '', textContent: '', innerHTML: '', hidden: false,
    dataset: {}, handlers: {}, attributes: {}, classList: { toggle() {} },
    setAttribute(key, value) { this.attributes[key] = value; },
    setCustomValidity() {}, addEventListener(type, fn) { (this.handlers[type] ??= []).push(fn); } });
  return nodes.get(selector);
};
const warnings = [];
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) }, URLSearchParams,
  document: { querySelector: node, querySelectorAll: () => [], addEventListener() {} },
  location: { pathname: '/', hash: '', search: '' }, history: { replaceState() {} },
  window: { scrollTo() {} }, setTimeout() {}, clearTimeout() {} });
for (const file of ['ligand-viewer.js', 'pocket-electrostatics.js', 'pocket-cloud.js', 'protein-viewer.js', 'app.js', 'analysis.js', 'go-analysis.js', 'control-qc.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
const plain = value => JSON.parse(JSON.stringify(value));
run(`
function fixture(id, pocket, score, probability = .9) {
  return { uniprot_id: id, protein: 'AF-' + id + '-F1-model_v4', pocket, rank: Number(pocket),
    vina_affinity: score, sfct_vina_score: score, sfct_score: score,
    vina_sfct_combined: score, vina_sfct_combined_50: score,
    vina_status: 'success', sfct_status: 'success', probability, mean_pocket_plddt: 95,
    residue_ids: 'A_1 A_2', center_x: 1, center_y: 2, center_z: 3 };
}
function resetFixture() {
  state.rawByAA.clear(); state.rankingCache.clear(); state.metadata.clear();
  Object.assign(state, { aa: 'ALA', metric: 'vina_affinity', p2rank: .7, plddt: 90, top: 100,
    maxCompetitors: 19, aaRank: 20, minDelta: null, requireLPreference: false, pocketMode: 'best',
    search: '', selectedProtein: null, selectedPocket: null, profilePocket: null, profilePocketAnchor: null });
}
`);

// Best-pocket selection and rank/percentile order agree in every consumer and score type.
for (const metric of ['vina_affinity', 'sfct_score', 'vina_sfct_combined', 'vina_sfct_combined_50']) {
  context.testMetric = metric;
  const check = plain(run(`(() => {
    resetFixture(); state.metric = testMetric;
    const rows = [fixture('Z','1',-8), fixture('A','1',-8), fixture('ALT','1',-20,.1),
      fixture('ALT','2',-9), fixture('FAILQC','1',-30,.1),
      { ...fixture('FAILED','1',-50), vina_status:'failed', sfct_status:'failed' }];
    state.rawByAA.set('ALA',rows); state.rawByAA.set('LEU',rows.map(row=>({...row})));
    const qc={metric:testMetric,p2:.7,plddt:90};
    const go={...goState, aa:'ALA', metric:testMetric, p2rank:.7, plddt:90, top:100};
    const ids=rows=>rows.map(row=>row.uniprot_id);
    return { explorer:ids(filterRows()), go:ids(goProteinSets(go).selected),
      compare:qualityRanking('ALA',qc).rows.map(row=>row.uniprot_id),
      overlap:[...topHitOverlap({...qc,aas:['ALA','LEU'],tier:100}).sets[0]],
      pockets:getRanking().map(row=>row.pocket), hist:getDistributionData().values,
      total:ids(goProteinSets({...go,background:'total'}).background),
      inspected:getRanking().map(row=>[row.uniprot_id,row.proteome_percentile,
        getProteinProfile(row.uniprot_id,testMetric,row.pocket)[0].proteome_percentile]),
      independentGO:ids(goProteinSets({...go,p2rank:0,plddt:0}).selected) };
  })()`));
  assert.deepEqual(check.explorer, ['ALT', 'A', 'Z']);
  for (const key of ['go', 'compare', 'overlap']) assert.deepEqual(check[key], check.explorer, `${metric}/${key}`);
  assert.deepEqual(check.pockets, ['2', '1', '1']);
  assert.deepEqual(check.hist, [-9, -8, -8]);
  assert.ok(check.total.includes('FAILQC'));
  assert.ok(check.independentGO.includes('FAILQC'), 'GO uses its own thresholds, not Explorer state');
  check.inspected.forEach(([, best, inspected]) => assert.equal(best, inspected));
}

// Compare inspected pockets only against OTHER proteins; preserve stable ties for every score.
run(`resetFixture(); state.rawByAA.set('ALA', [fixture('P1','1',-10),fixture('P1','2',-9),fixture('P2','1',-8)]);`);
assert.equal(run(`addPocketProteomePosition(state.rawByAA.get('ALA')[1], 'ALA','vina_affinity').proteome_percentile`), 0);
for (let i = 0; i < 30; i++) {
  context.testIndex = i;
  assert.equal(run(`(() => {
    const rows=Array.from({length:30},(_,n)=>fixture('P'+String(n).padStart(2,'0'),'1',-10+n%5));
    state.rawByAA.set('ALA',rows);
    const row={...rows[testIndex],pocket:'2',vina_affinity:-10+(testIndex%7)/2};
    const pop=rankedPopulation('ALA','vina_affinity');
    const expected=100*pop.rows.filter(other=>other.uniprot_id!==row.uniprot_id && compareProteinScores(other,row,'vina_affinity')<0).length/29;
    return pocketProteomePosition(pop,row,'vina_affinity').proteome_percentile-expected;
  })()`), 0);
}
run(`resetFixture(); state.rawByAA.set('ALA',[fixture('ONE','1',-8)]);`);
assert.equal(run(`getProteinProfile('ONE','vina_affinity','1')[0].proteome_percentile`), 0);
assert.ok(run(`addPocketProteomePosition(fixture('NONE','1',-8),'ALA','vina_affinity') === null`));

// Tied scores straddling tiers have the same membership in table and histogram, not all Top 1%.
run(`resetFixture(); state.rawByAA.set('ALA',Array.from({length:100},(_,i)=>fixture('P'+String(i).padStart(3,'0'),'1',-8)).reverse());`);
for (const tier of [1, 5, 10]) {
  context.testTier = tier;
  assert.equal(run('state.top=testTier; filterRows().length'), tier);
}
const histogram = run(`(() => {const d=getDistributionData();return histogramSVG(d.values,d.cutoffs,d.referenceValues,d.percentiles);})()`);
for (const [label, count] of [['Top 1%', 1], ['Top 1\u20135%', 4], ['Top 5\u201310%', 5], ['Outside top 10%', 90]]) {
  assert.ok(histogram.includes(`${label}: ${count} protein`));
}
run(`state.distributionSource='filtered'; state.top=100; state.search='P099'; filterRows();`);
assert.ok(run(`(() => {const d=getDistributionData();return histogramSVG(d.values,d.cutoffs,d.referenceValues,d.percentiles);})()`)
  .includes('Outside top 10%: 1 protein'));
run(`state.distributionSource='all';`);

// Cross-AA inspected pockets require model, pocket center and residue-set agreement.
run(`resetFixture();
  state.rawByAA.set('ALA',[fixture('P1','1',-10),fixture('P1','2',-7)]);
  state.rawByAA.set('LEU',[fixture('P1','1',-11),{...fixture('P1','2',-6),residue_ids:'A_2 A_1'}]);
  state.rawByAA.set('VAL',[{...fixture('P1','2',-8),protein:'AF-P1-F2-model_v4'}]);
  state.rawByAA.set('ARG',[{...fixture('P1','2',-8),center_x:300}]);
  state.rawByAA.set('SER',[{...fixture('P1','2',-8),residue_ids:'A_9 A_10'}]);
  state.rawByAA.set('THR',[fixture('P1','2',-8,.1)]);
`);
assert.deepEqual(Array.from(run(`getProteinProfile('P1','vina_affinity','2').map(row=>row.code)`)), ['ALA', 'LEU']);
run(`const originalSwitchView=switchView; switchView=()=>{}; selectProtein('P1','2');`);
assert.equal(run('state.profilePocket'), '2');
assert.equal(run('state.profilePocketAnchor'), run(`state.rawByAA.get('ALA')[1]`));
assert.deepEqual(Array.from(run(`getProteinProfile('P1','vina_affinity',state.profilePocket).map(row=>row.value)`)), [-7, -6]);
run('switchView=originalSwitchView;');

// L>D uses QC-first L and D populations, including in GO; ties/missing/mismatched sites cannot pass.
run(`state.rawByAA.set('DALA',[fixture('P1','1',-20,.1),fixture('P1','2',-5)]);`);
assert.equal(run(`getStereoControl('P1','ALA','vina_affinity').lPreferred`), true);
assert.equal(run(`getStereoControl('P1','ALA','vina_affinity','2').lPreferred`), true);
run(`state.rawByAA.set('DALA',[{...fixture('P1','2',-5),center_z:900}]);`);
assert.equal(run(`getStereoControl('P1','ALA','vina_affinity','2')`), null);
run(`state.rawByAA.set('DALA',[fixture('P1','2',-5,.1)]);`);
assert.equal(run(`getStereoControl('P1','ALA','vina_affinity')`), null);
assert.equal(run(`goProteinSets({...goState,aa:'ALA',metric:'vina_affinity',p2rank:.7,plddt:90,top:100,requireLPreference:true}).selected.length`), 0);

// Control QC and profile use the same percentile helper and exact matched-pocket anchor.
run(`resetFixture();
  state.rawByAA.set('ALA',[{...fixture('P1','1',-10),residue_ids:'A_9'},fixture('P1','2',-9),fixture('P2','1',-8)]);
  const controls={controls:[{aa:'ALA',protein:'P1',pdb:'1ABC',positions:[1,2],source:{}}],excluded:0};
  const qcResult=calculateControlQC(controls,{...controlQCState,aa:'ALA',pocketMode:'matched'});`);
assert.equal(run('qcResult.rows[0].scoreRows[0].pocket'), '2');
assert.equal(run('qcResult.rows[0].percentiles[0]'), 0);
assert.equal(run(`getProteinProfile('P1','vina_affinity','2')[0].proteome_percentile`), 0);

// Changing data or thresholds cannot reuse a different population's cached ranking.
run(`state.rawByAA.set('ALA',[fixture('NEW','1',-10)]);`);
assert.equal(run(`qualityRanking('ALA',{metric:'vina_affinity',p2:.7,plddt:90}).rows[0].uniprot_id`), 'NEW');
assert.equal(run(`qualityRanking('ALA',{metric:'vina_affinity',p2:1,plddt:90}).rows.length`), 0);

// Versioned model URLs, caching by model rather than accession, and explicit failure with no fallback.
assert.equal(run(`dockingModelReference(fixture('P12345','1',-8)).url`), 'https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.cif');
assert.throws(() => run(`dockingModelReference({...fixture('P1','1',-8),protein:'AF-P2-F1-model_v4'})`));
assert.throws(() => run(`dockingModelReference({...fixture('P1','1',-8),protein:'unversioned_model'})`));
const urls = [], presets = [];
context.window.molstar = { lib: { structure: { StructureElement: {
  Loci: { fromSchema: (structure, schema) => ({ structure, schema }), isEmpty: () => false },
  Bundle: { fromLoci: loci => loci },
} }, plugin: { StateTransforms: { Representation: { OverpaintStructureRepresentation3DFromBundle: {} } } } } };
const viewer = { plugin: { async clear() {},
  state: { data: { updateCellState() {} } },
  build: () => ({ delete: () => ({ async commit() {} }),
    to: () => ({ applyOrUpdateTagged() {}, update() {} }), async commit() {} }),
  managers: { camera: { focusLoci() {} } }, builders: {
  data: { async download({ url }) { urls.push(url); return {}; } },
  structure: { async parseTrajectory() { return {}; },
    async tryCreateComponent() { return { ref: 'pocket-component' }; },
    representation: { async addRepresentation() { return { ref: 'pocket-sticks' }; } },
    hierarchy: { async applyPreset(_, __, options) {
      presets.push(options);
      const structure = {};
      return { structureProperties: { ref: 'protein', obj: { data: structure } },
        representation: { components: { polymer: { ref: 'polymer' } },
          representations: { polymer: { ref: 'cartoon', obj: { data: { sourceData: { root: structure } } } } } } };
    } },
  },
} }, async structureInteractivity() {} };
context.testViewer = viewer;
run('getMolstarViewer=async()=>testViewer;');
(async () => {
  for (const version of [4, 5, 5]) {
    context.testVersion = version;
    await run(`(() => {const row={...fixture('P12345','1',-8),protein:'AF-P12345-F1-model_v'+testVersion};
      const request={protein:'P12345',row};molstarLatestRequest=request;return performMolstarUpdate(request);})()`);
  }
  assert.deepEqual(urls, ['https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.cif', 'https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v5.cif']);
  assert.ok(presets.every(options => options.representationPreset === 'polymer-cartoon'));
  assert.ok(presets.every(options => options.representationPresetParams.theme.globalName === 'plddt-confidence'));
  assert.ok(presets.every(options => options.structure.name === 'model'));
  viewer.plugin.builders.data.download = async ({ url }) => { urls.push(url); throw Error('Model unavailable'); };
  run(`state.selectedProtein='P12345'; state.currentView='protein'; updateMolstarPocket({...fixture('P12345','1',-8),protein:'AF-P12345-F1-model_v6'});`);
  await run('molstarUpdateQueue');
  assert.equal(node('#structure-viewer-message').textContent, 'Exact docking model unavailable');
  assert.ok(node('#structure-viewer-detail').textContent.includes('No alternate version'));
  assert.equal(urls.length, 3);
  assert.equal(warnings.length, 1, 'The intentional missing-model case reports its failure');
  console.log('Ranking consistency passed: all four metrics, QC-first populations, GO isolation/backgrounds, inspected ranks, ties, histogram tiers, pocket context/geometry, L/D, Control QC and pinned model versions.');
})().catch(error => { console.error(error); process.exitCode = 1; });
