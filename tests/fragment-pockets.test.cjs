// Repeated pocket names across model fragments: loading, selection and protein-level science.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const zlib = require('node:zlib');
const { createHash, webcrypto } = require('node:crypto');
const root = path.resolve(__dirname, '..');
const hash = text => createHash('sha256').update(text).digest('hex');
const files = new Map(), nodes = new Map(), downloads = [];
let pocketNodes = [];
const node = selector => {
  if (!nodes.has(selector)) nodes.set(selector, { value: '', innerHTML: '', textContent: '', dataset: {},
    classList: { toggle() {} }, setAttribute() {}, setCustomValidity() {}, handlers: {},
    addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
};
const context = vm.createContext({ console, URLSearchParams, TextEncoder, TextDecoder, crypto: webcrypto,
  Blob, Response, DecompressionStream, location: { pathname: '/', hash: '', search: '' },
  history: { replaceState() {} }, window: { scrollTo() {} }, setTimeout() {}, clearTimeout() {},
  document: { querySelector: node, querySelectorAll: selector => selector === '#pockets-body tr' ? pocketNodes : [], addEventListener() {} },
  fetch: async url => new Response(files.get(url.split('?')[0]) || '', { status: files.has(url.split('?')[0]) ? 200 : 404 }),
  captureDownload: (name, text) => downloads.push({ name, text }),
});
for (const file of ['organisms.js', 'ligand-viewer.js', 'compact-data.js', 'app.js', 'analysis.js', 'go-analysis.js', 'control-qc.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
const plain = value => JSON.parse(JSON.stringify(value));
run(`
function fixture(fragment, vina, sfct) {
  return { uniprot_id:'A2VEC9', protein:'AF-A2VEC9-F'+fragment+'-model_v6', pocket:'pocket3', rank:3,
    probability:.9, mean_pocket_plddt:95, residue_ids:'A_1 A_2', center_x:1, center_y:2, center_z:3,
    vina_affinity:vina, vina_ala_affinity:vina, vina_status:'success',
    sfct_vina_score:sfct, sfct_score:sfct, vina_sfct_combined:sfct, sfct_best_pose:0, sfct_n_poses:9,
    sfct_status:Number.isFinite(sfct)?'success':'missing' };
}
const sourceRows=[fixture(1,-5,-4),fixture(17,-12,NaN),fixture(20,-8,-9),fixture(5,-7,-6),fixture(6,-6,-5)];
Object.assign(state,{aa:'ALA',metric:'vina_affinity',p2rank:.7,plddt:90});
proteinTableIdentity=id=>id;
poseElectrostaticsCell=()=>'<td>—</td>';
poseElectrostaticsForRow=()=>({value:NaN,status:'unavailable',pose:null,error:'Missing optional dataset'});
downloadText=captureDownload;
switchView=()=>{};
let viewedPocket=null;
updateMolstarPocket=row=>{viewedPocket=row;};
`);

function putTable(file, fields, rows) {
  const text = fields.join('\t') + '\n' + rows.map(row => fields.map(field => Number.isNaN(row[field]) ? '' : row[field] ?? '').join('\t')).join('\n') + '\n';
  files.set('Ec_results/L/' + file, zlib.gzipSync(text));
  return { file, rows: rows.length, sha256: hash(text) };
}
const raw = run('sourceRows').map((row, i) => ({ ...row, pocket_id: String(i + 1) }));
const metadataFields = ['pocket_id', ...run('COMPACT_METADATA_FIELDS')];
const scoreFields = ['pocket_id', ...run('COMPACT_SCORE_FIELDS')];
const pockets = putTable('pockets.tsv.gz', metadataFields, raw);
const scores = putTable('scores_ala.tsv.gz', scoreFields, raw);
files.set('Ec_results/L/manifest.json', JSON.stringify({ format:'aa-proteome-atlas-compact', version:1,
  pockets, ligands:[{ ...scores, code:'ALA' }] }));

(async () => {
  // L only, exactly the user's five fragments, including independent Vina success for missing SFCT.
  const data = await run('loadCompactData()');
  context.compact = data;
  assert.equal(data.byAA.size, 1);
  await run(`loadCompactResultRows({code:'ALA'},compact).then(rows=>state.rawByAA.set('ALA',rows))`);
  assert.equal(run('state.rawByAA.get("ALA").length'), 5);
  assert.equal(run('getRanking().length'), 1);
  assert.equal(run('getRanking()[0].protein'), 'AF-A2VEC9-F17-model_v6');
  assert.equal(run('getRanking()[0].proteome_percentile'), 0);
  assert.equal(run(`getRanking('ALA','sfct_score')[0].protein`), 'AF-A2VEC9-F20-model_v6');
  assert.ok(run('Number.isNaN(state.rawByAA.get("ALA")[1].sfct_score)'));
  assert.equal((await run(`loadCompactResultRows({code:'DALA'},compact)`)).length, 0);
  const legacy = run(`mergeResultRows({code:'ALA',file:'ala'},sourceRows,[...sourceRows].reverse())`);
  assert.deepEqual(plain(legacy.map(row => [row.protein,row.vina_affinity,row.sfct_score])),
    plain(run('state.rawByAA.get("ALA").map(row=>[row.protein,row.vina_affinity,row.sfct_score])')));
  assert.throws(() => run(`mergeResultRows({code:'ALA',file:'ala'},[...sourceRows,sourceRows[0]],[])`), /Duplicate/);
  assert.throws(() => run(`mergeResultRows({code:'ALA',file:'ala'},sourceRows,[{...sourceRows[0],center_x:99}])`), /metadata mismatch/);
  const scoreEntry = data.byAA.get('ALA').ligands.get('ALA');
  Object.assign(scoreEntry, putTable('scores_ala.tsv.gz', scoreFields, [...raw, raw[0]]));
  await assert.rejects(run(`loadCompactResultRows({code:'ALA'},compact)`), /Duplicate compact/);
  Object.assign(scoreEntry, putTable('scores_ala.tsv.gz', scoreFields, raw));

  // An unqualified pocket name is ambiguous; exact identity works regardless of row order.
  assert.equal(run(`profilePocketAnchor('A2VEC9','pocket3')`), null);
  run(`selectProtein('A2VEC9','pocket3')`);
  assert.equal(run('state.selectedProtein'), null);
  run(`const f20=state.rawByAA.get('ALA')[2]; selectProtein('A2VEC9',pocketKey(f20));`);
  assert.equal(run('state.profilePocketAnchor.protein'), 'AF-A2VEC9-F20-model_v6');
  run(`state.rawByAA.set('ALA',[...state.rawByAA.get('ALA')].reverse());
    state.rawByAA.set('LEU',[...state.rawByAA.get('ALA')].map(row=>({...row})));
    state.rawByAA.set('VAL',[{...f20,protein:'AF-A2VEC9-F20-model_v4'}]);
    state.rawByAA.set('ARG',[{...f20,center_x:99}]);
    state.rawByAA.set('DALA',[{...sourceRows[0],vina_affinity:-30}]);`);
  assert.deepEqual(Array.from(run(`getProteinProfile('A2VEC9','vina_affinity',state.profilePocket).map(row=>row.code)`)), ['ALA','LEU']);
  assert.equal(run(`getStereoControl('A2VEC9','ALA','vina_affinity',state.profilePocket)`), null);
  assert.equal(run(`getStereoControl('A2VEC9').dRow.protein`), 'AF-A2VEC9-F1-model_v6');
  run(`state.rawByAA.get('DALA').push({...f20,vina_affinity:-6});`);
  assert.equal(run(`getStereoControl('A2VEC9','ALA','vina_affinity',state.profilePocket).delta`), 2);

  // Scientific units remain explicit: one protein, five model-pockets (four with SFCT).
  assert.deepEqual(plain(run(`(()=>{const s=macroStatistics({...analysisState.statistics,metric:'vina_affinity'});
    return [s.proteins,s.pockets,s.records.length,s.pocketsPerPair];})()`)), [1,5,1,[5]]);
  assert.equal(run(`macroStatistics({...analysisState.statistics,metric:'sfct_score',unit:'all'}).pockets`), 4);
  assert.equal(run(`getProfileScoreSummary('ALA','vina_affinity').count`), 1);
  assert.equal(run(`goProteinSets({...goState,aa:'ALA',metric:'vina_affinity',top:100}).background.length`), 1);
  run(`state.pocketMode='all'; renderResults(filterRows()); renderPockets('A2VEC9',f20);`);
  assert.equal(run('viewedPocket.protein'), 'AF-A2VEC9-F20-model_v6');
  assert.equal((node('#pockets-body').innerHTML.match(/pocket-best-badge/g) || []).length, 1);
  assert.match(node('#pockets-body').innerHTML, /pocket3 · F20 v6/);
  assert.match(node('#pocket-detail-list').innerHTML, /AF-A2VEC9-F20-model_v6/);
  for (const fragment of [1,17,20,5,6]) {
    assert.ok(node('#results-body').innerHTML.includes(`AF-A2VEC9-F${fragment}-model_v6`));
  }
  pocketNodes = raw.map(row => ({ dataset:{ pocket:run(`pocketKey(sourceRows.find(row=>row.protein===${JSON.stringify(row.protein)}))`) },
    classList:{ toggle(name, active) { this.active = active; } } }));
  run('renderPocketDetail(f20)');
  assert.deepEqual(Array.from(pocketNodes,row=>row.classList.active), [false,false,true,false,false]);

  // Exercise the actual Explorer and Inspect button handlers with exact keys.
  for (const name of ['bindProteinStyleEvents','bindProteinColorEvents','bindPocketSticksEvents',
    'bindLigandViewerEvents','bindPocketCloudEvents','bindGeneDescriptionEvents','bindAnalysisEvents',
    'bindControlQCEvents','bindGOEvents']) context[name]=()=>{};
  run('bindProfileProteinSearchEvents=()=>{}; bindEvents(); renderProtein=()=>{};');
  const key = run('pocketKey(sourceRows[4])');
  node('#results-body').handlers.click({target:{closest:()=>({dataset:{protein:'A2VEC9',pocket:key}})}});
  assert.equal(run('state.profilePocketAnchor.protein'), 'AF-A2VEC9-F6-model_v6');
  node('#pockets-body').handlers.click({target:{closest:()=>({dataset:{pocket:run('pocketKey(f20)')}})}});
  assert.equal(run('state.profilePocketAnchor.protein'), 'AF-A2VEC9-F20-model_v6');

  // Downloads carry exact models, and internal selection keys never leak into filenames.
  const parse = text => { const [header,...rows]=text.split('\n').map(line=>line.split('\t'));
    return rows.map(row=>Object.fromEntries(header.map((field,i)=>[field,row[i]]))); };
  const exported = parse(run(`rowsToDelimited(filterRows(),'\t')`));
  assert.equal(new Set(exported.map(row=>row.model)).size, 5);
  assert.ok(exported.every(row=>row.d_control_model==='AF-A2VEC9-F1-model_v6'));
  run(`downloadProfile(); downloadProfileTable('\t');`);
  assert.ok(downloads.every(item=>item.name.includes('AF-A2VEC9-F20-model_v6_pocket3')));
  assert.ok(parse(downloads[0].text).every(row=>row.model==='AF-A2VEC9-F20-model_v6'));
  assert.equal(parse(downloads[1].text).find(row=>row.AA.startsWith('ALA')).Model, 'AF-A2VEC9-F20-model_v6');
  assert.ok(run(`pairedAAResults({...analysisState.compare,metric:'vina_affinity'}).rows.every(row=>row.modelX && row.modelY)`));

  // Scores from non-F1 models remain usable; site overlap does not assume a numbering offset.
  assert.equal(run('controlQCPocket(f20)'), null);
  run(`const controls={controls:[{aa:'ALA',protein:'A2VEC9',pdb:'1ABC',positions:[1,2],source:{}}],excluded:0};
    const qc=calculateControlQC(controls,{...controlQCState,aa:'ALA'});`);
  assert.equal(run('qc.rows[0].eligible'), true);
  assert.equal(run('qc.rows[0].scoreOverlaps[0]'), null);
  run(`controlQCResult={result:qc,options:{...controlQCState,aa:'ALA'}}; downloadControlQC();`);
  const qcExport = parse(downloads.at(-1).text)[0];
  assert.equal(qcExport.vina_affinity_model, 'AF-A2VEC9-F17-model_v6');
  assert.equal(qcExport.sfct_score_model, 'AF-A2VEC9-F20-model_v6');
  assert.equal(qcExport.vina_affinity_d_model, 'AF-A2VEC9-F1-model_v6');
  run(`goResult={options:{...goState,aa:'ALA',metric:'vina_affinity',top:100},
    data:{proteins:new Map()},background:getRanking()}; downloadGOProteins(true);`);
  assert.equal(parse(downloads.at(-1).text)[0].model, 'AF-A2VEC9-F17-model_v6');
  assert.equal(run('dockingModelReference(f20).url'), 'https://alphafold.ebi.ac.uk/files/AF-A2VEC9-F20-model_v6.cif');
  console.log('Fragment pockets passed: L-only compact/legacy loading, exact model selection, duplicate/geometry rejection, rankings, statistics, GO, L/D, downloads and QC mapping limits.');
})().catch(error=>{console.error(error);process.exitCode=1;});
