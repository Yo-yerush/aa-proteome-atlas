// Explorer qphi annotation: read-only file validation, independent sorting and download checks.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const zlib = require('node:zlib');
const { createHash, webcrypto } = require('node:crypto');
const root = path.resolve(__dirname, '..');
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const nodes = new Map();
const node = selector => {
  if (!nodes.has(selector)) nodes.set(selector, { value: '', textContent: '', innerHTML: '', hidden: false, disabled: false,
    dataset: {}, handlers: {}, setAttribute() {}, classList: { toggle() {} },
    addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
};
const warnings = [], downloads = [], toasts = [], fetches = [];
let files = new Map(), gates = new Map();
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) },
  document: { querySelector: node, querySelectorAll: () => [], addEventListener() {} },
  crypto: webcrypto, TextEncoder, TextDecoder, Blob, Response, DecompressionStream, URLSearchParams,
  location: { pathname: '/', hash: '', search: '' }, history: { replaceState() {} },
  fetch: async url => {
    const file = url.split('?')[0]; fetches.push(file);
    if (gates.has(file)) await gates.get(file);
    return new Response(files.get(file) || '', { status: files.has(file) ? 200 : 404 });
  },
  downloadText: (name, text) => downloads.push({ name, text }),
  showToast: text => toasts.push(text),
});
for (const file of ['pocket-electrostatics.js', 'pose-electrostatics.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
context.captureDownload = (name, text) => downloads.push({ name, text });
context.captureToast = text => toasts.push(text);
run(`downloadText=captureDownload;showToast=captureToast;proteinTableIdentity=id=>id;
function scoreRow(id, protein, rank, score, aa='ALA') {
  return { pocket_id:String(id),uniprot_id:protein,protein:'AF-'+protein+'-F1-model_v6',pocket:'pocket'+rank,rank,
    probability:.9,mean_pocket_plddt:95,vina_affinity:score,sfct_vina_score:score,sfct_score:score,
    vina_sfct_combined:score,vina_sfct_combined_50:score,vina_status:'success',sfct_status:'success',sfct_best_pose:3,
    _compactDirectory:'At_results/L',_compactAA:aa,_compactPocketHash:'a'.repeat(64) };
}
function scoreFixture() {
  state.rawByAA.clear();state.rankingCache.clear();
  for(const aa of ['ALA','MET','LEU']) state.rawByAA.set(aa,[
    scoreRow(1,'P1',1,-9,aa),scoreRow(2,'P1',2,-8,aa),scoreRow(3,'P2',1,-7,aa),
    scoreRow(4,'P3',1,-6,aa),scoreRow(5,'P4',1,-5,aa),scoreRow(6,'P5',1,-4,aa),scoreRow(7,'P6',1,-3,aa)]);
  Object.assign(state,{aa:'ALA',metric:'vina_sfct_combined_50',pocketMode:'best',sortKey:'percentile',sortDirection:'asc',
    search:'',top:100,p2rank:.7,plddt:90,maxCompetitors:19,aaRank:20,minDelta:null,requireLPreference:false,pageSize:10,page:1});
  filterRows();
}`);
const columns = ['pocket_id', 'vina_pose', 'n_atoms', 'phi_mean', 'phi_min', 'phi_max', 'qphi_kT', 'status', 'error'];
const dataRows = [
  ['1','1','13','','','','-2.1234','success',''],
  ['2','1','13','','','','-9','success',''],
  ['3','1','13','','','','0','success',''],
  ['4','1','13','','','','-999','unreliable_clash','Overlap, exceeds cutoff'],
  ['5','','','','','','','receptor_error','Bad receptor <script>'],
  ['6','1','13','','','','','success',''],
  ['7','1','13','','','','3.3333','success',''],
];
const tableText = rows => columns.join('\t') + '\n' + rows.map(row => row.join('\t')).join('\n') + '\n';
const compact = Buffer.from(JSON.stringify({ format: 'aa-proteome-atlas-compact', version: 1, pockets: { sha256: 'a'.repeat(64) } }));
const manifest = { format: 'aa-pocket-electrostatics', version: 1,
  poses: { columns, pose: 'Vina MODEL 1', qphi_units: 'kT' },
  bundles: [{ manifest_sha256: hash(compact), codes: ['ALA', 'MET', 'LEU'], pose_files: {} }] };
function reset() {
  files = new Map(); gates = new Map(); fetches.length = 0; downloads.length = 0; toasts.length = 0;
  for (const aa of ['ALA', 'MET', 'LEU']) {
    const rows = dataRows.map(row => [...row]);
    if (aa !== 'ALA') rows[0][6] = aa === 'MET' ? '8.5' : '9.5';
    const bytes = zlib.gzipSync(tableText(rows)), file = `pose_electrostatics/electrostatics_${aa.toLowerCase()}.tsv.gz`;
    manifest.bundles[0].pose_files[aa] = { file, rows: rows.length, sha256: hash(bytes) };
    files.set('At_results/L/' + file, bytes);
  }
  files.set('At_results/electrostatics/manifest.json.gz', zlib.gzipSync(JSON.stringify(manifest)));
  files.set('At_results/L/manifest.json', compact);
  run('electrostaticsManifestPromise=null;poseElectrostaticsBundlePromise=null;poseElectrostaticsTables.clear();explorerPoseObservation=null;scoreFixture();');
}
const request = () => run('prepareExplorerPoseElectrostatics()');
const value = id => run(`poseElectrostaticsForRow(state.rawByAA.get(state.aa).find(row=>row.pocket_id==='${id}'))`);
const waitFor = async predicate => {
  for (let i=0;i<1000;i++) { if (predicate()) return; await new Promise(resolve => setTimeout(resolve, 1)); }
  throw new Error('Test event timeout');
};

(async () => {
  reset();
  assert.equal(fetches.length, 0);
  const before = run('JSON.stringify(state.rawByAA.get("ALA"))');
  const baseline = run('JSON.stringify(filterRows().map(row=>[row.uniprot_id,row.pocket,row.proteome_rank,row.proteome_percentile,row.comparison]))');
  run('renderExplorer()');
  const loading = request();
  assert.equal(node('#download-tsv').disabled, true);
  assert.equal(value(1).status, 'loading');
  const entry = await loading;
  assert.equal(entry.kind, 'ready');
  assert.equal(fetches.length, 3, 'Only shared manifest, compact manifest and selected AA pose file');
  assert.ok(fetches.every(file => !file.includes('/points_') && !file.includes('pocket_summary') && !file.includes('aa_positions')));
  assert.equal(node('#download-tsv').disabled, false);
  assert.equal(node('#explorer-qphi-status').hidden, true);
  assert.equal(value(1).value, -2.1234);
  assert.equal(value(3).value, 0, 'Valid zero stays numeric');
  assert.ok(Number.isNaN(value(4).value), 'A numeric value on a clashing row is never used');
  assert.equal(value(4).status, 'unreliable_clash');
  assert.equal(value(5).status, 'receptor_error');
  assert.equal(value(6).status, 'missing_value');
  assert.equal(run('JSON.stringify(state.rawByAA.get("ALA"))'), before, 'No mutation of docking records');
  assert.equal(run('JSON.stringify(filterRows().map(row=>[row.uniprot_id,row.pocket,row.proteome_rank,row.proteome_percentile,row.comparison]))'), baseline);
  assert.ok(node('#results-body').innerHTML.includes('-2.123</td>'));
  assert.ok(node('#results-body').innerHTML.includes('Vina pose 1'));
  assert.ok(node('#results-body').innerHTML.includes('qφ (kT)'));
  const page = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
  assert.match(page, /data-sort="qphi_kT">qφ \(kT\)/, 'Display notation does not change the sorting key');
  assert.ok(page.includes('<th>qφ (kT)</th>'), 'Methods use the same notation');
  assert.ok(page.includes('positive potential favors –COO⁻'));
  assert.ok(page.includes('negative potential favors –NH₃⁺'));
  assert.ok(node('#results-body').innerHTML.includes('saved SFCT/Combined pose: Vina MODEL 4 (SFCT index 3)'));
  assert.ok(node('#results-body').innerHTML.includes('&lt;script&gt;'));
  assert.ok(!node('#results-body').innerHTML.includes('<script>'));
  assert.equal((node('#results-body').innerHTML.match(/class="numeric qphi-cell"/g) || []).length, 6);

  // Sorting affects the displayed rows, not best-pocket selection, ranks or cutoffs.
  run("state.sortKey='qphi_kT';state.sortDirection='asc';filterRows();");
  assert.deepEqual(Array.from(run('state.filtered.map(row=>row.pocket_id)')), ['1','3','7','4','5','6']);
  run("state.sortDirection='desc';filterRows();");
  assert.deepEqual(Array.from(run('state.filtered.map(row=>row.pocket_id)')), ['7','3','1','4','5','6']);
  assert.equal(run('state.filtered.find(row=>row.uniprot_id==="P1").pocket'), 'pocket1', 'Never choose the lower-qphi pocket');
  run("state.pocketMode='all';state.sortDirection='asc';filterRows();");
  assert.equal(run('state.filtered[0].pocket_id'), '2', 'All-pockets rows use their own pocket ID');
  assert.equal(run('state.filtered[0].isBestPocket'), false);
  for (const metric of ['vina_affinity','sfct_score','vina_sfct_combined','vina_sfct_combined_50']) {
    context.metric = metric; run('state.metric=metric;filterRows();'); assert.equal(value(1).value, -2.1234);
  }
  await request(); assert.equal(fetches.length, 3, 'No refetch when sorting/filtering');

  // Exports retain full precision and explicit pose/status, across the whole filtered list.
  run('state.pageSize=1;state.page=2;renderResults(state.filtered);');
  for (const delimiter of ['\t', ',']) {
    context.delimiter = delimiter;
    await run('downloadFiltered(delimiter)');
    const download = downloads.at(-1);
    const lines = download.text.split('\n');
    assert.equal(lines.length, 8);
    assert.ok(lines[0].endsWith(['qphi_kT','qphi_status','qphi_vina_pose','qphi_error'].join(delimiter)));
    assert.ok(download.text.includes('-2.1234'));
    assert.ok(download.text.includes([0,'success',1,''].join(delimiter)));
    assert.ok(!download.text.includes('-999'));
    assert.ok(download.text.includes('unreliable_clash'));
    if (delimiter === ',') assert.ok(download.text.includes('"Overlap, exceeds cutoff"'));
  }
  assert.equal(run('state.page'), 2, 'Loading/downloading does not reset pagination');

  run("state.aa='MET';filterRows();"); await request();
  assert.equal(value(1).value, 8.5, 'AA key prevents matching another ligand');
  run("state.aa='LEU';filterRows();"); await request();
  assert.equal(run('poseElectrostaticsTables.size'), 2);
  assert.equal(fetches.filter(file => file.endsWith('manifest.json.gz')).length, 1);
  const tableLoads = fetches.filter(file => file.includes('/pose_electrostatics/')).length;
  run("state.aa='ALA';filterRows();"); await request();
  assert.equal(fetches.filter(file => file.includes('/pose_electrostatics/')).length, tableLoads + 1, 'Evicted AA can reload');

  // No D/legacy fallback, checksum bypass, arbitrary remote paths, or incompatible pose units.
  for (const directory of ['At_results/D', 'At_results']) {
    context.directory = directory;
    assert.equal(run('poseElectrostaticsSource({...state.rawByAA.get("ALA")[0],_compactDirectory:directory},"ALA")'), null);
  }
  assert.equal(run('poseElectrostaticsForRow({...state.rawByAA.get("ALA")[0],_compactAA:"MET"}).status'), 'unavailable');
  reset();
  run("state.rawByAA.get('ALA')[0]._compactPocketHash='b'.repeat(64)");
  assert.equal((await request()).kind, 'unavailable');
  assert.match(node('#explorer-qphi-status').textContent, /export does not match/);
  reset(); files.set('At_results/L/manifest.json', Buffer.from(compact.toString()+'\n'));
  assert.equal((await request()).kind, 'unavailable');
  assert.match(node('#explorer-qphi-status').textContent, /current compact L bundle/);
  reset();
  const aaFile = 'At_results/L/pose_electrostatics/electrostatics_ala.tsv.gz';
  const goodBytes = files.get(aaFile);
  files.set(aaFile, zlib.gunzipSync(goodBytes));
  assert.equal((await request()).kind, 'unavailable');
  assert.match(node('#explorer-qphi-status').textContent, /checksum mismatch/);
  assert.equal(node('#explorer-qphi-retry').hidden, false);
  assert.equal(node('#download-tsv').disabled, false, 'A broken optional file cannot disable docking downloads');
  files.set(aaFile, goodBytes);
  run('retryExplorerPoseElectrostatics()'); await request();
  assert.equal(value(1).value, -2.1234);
  reset(); files.delete(aaFile);
  assert.equal((await request()).kind, 'unavailable'); assert.match(node('#explorer-qphi-status').textContent, /404/);
  assert.equal(run('filterRows().length'), 6, 'Optional missing data never excludes docking rows');
  reset(); context.crypto = undefined;
  assert.equal((await request()).kind, 'unavailable'); assert.match(node('#explorer-qphi-status').textContent, /HTTPS or localhost/);
  context.crypto = webcrypto;
  for (const poses of [{...manifest.poses,pose:'SFCT pose'}, {...manifest.poses,qphi_units:'kcal/mol'}]) {
    reset(); files.set('At_results/electrostatics/manifest.json.gz', zlib.gzipSync(JSON.stringify({...manifest,poses})));
    assert.equal((await request()).kind, 'unavailable');
    assert.match(node('#explorer-qphi-status').textContent, /format, units or pose/);
  }
  for (const badRows of [
    [dataRows[0],dataRows[0]],
    [[...dataRows[0].slice(0,1),'2',...dataRows[0].slice(2)]],
    [[...dataRows[0].slice(0,6),'NaN',...dataRows[0].slice(7)]],
    [[...dataRows[0].slice(0,6),'Infinity',...dataRows[0].slice(7)]],
  ]) {
    context.bad = tableText(badRows); context.count = badRows.length;
    assert.throws(() => run('parsePoseElectrostatics(bad,count)'));
  }
  context.bad = tableText(dataRows); assert.throws(() => run('parsePoseElectrostatics(bad,2)'), /row count/);

  // Delayed ALA requests cannot overwrite MET rows or produce a mixed-AA download.
  reset();
  let release;
  gates.set(aaFile,new Promise(resolve=>{release=resolve;}));
  const pending = request();
  await waitFor(()=>fetches.includes(aaFile));
  const pendingDownload = run("downloadFiltered('\t')");
  run("state.aa='MET';filterRows();"); await request();
  const metHTML = node('#results-body').innerHTML;
  release(); await pending; await pendingDownload;
  assert.equal(node('#results-body').innerHTML, metHTML);
  assert.equal(value(1).value, 8.5);
  assert.equal(downloads.length, 0); assert.match(toasts.at(-1), /target AA changed/);

  if (process.argv.includes('--real-data')) {
    const gunzip = file => zlib.gunzipSync(fs.readFileSync(path.join(root,file))).toString();
    const realManifest = JSON.parse(gunzip('At_results/electrostatics/manifest.json.gz'));
    const realCompact = JSON.parse(fs.readFileSync(path.join(root,'At_results/L/manifest.json')));
    run('electrostaticsManifestPromise=null;poseElectrostaticsBundlePromise=null;poseElectrostaticsTables.clear();explorerPoseObservation=null;');
    context.fetch = async url => {
      const file=path.resolve(root,url.split('?')[0]); assert.ok(file.startsWith(root+path.sep));
      return new Response(fs.readFileSync(file));
    };
    let checked=0,successful=0,unavailable=0;
    for (const aa of realManifest.bundles[0].codes) {
      context.realSource = { aa,directory:'At_results/L',pocketHash:realCompact.pockets.sha256,
        key:`At_results/L|${aa}|${realCompact.pockets.sha256}` };
      const entry=await run('ensurePoseElectrostaticsTable(realSource,realSource.aa).promise');
      assert.equal(entry.kind,'ready',`${aa}: ${entry.error}`);
      assert.equal(entry.table.size,realManifest.bundles[0].pose_files[aa].rows);
      for (const record of entry.table.values()) {
        if (record.status==='success') { assert.ok(Number.isFinite(record.value));assert.equal(record.pose,1);successful++; }
        else { assert.ok(Number.isNaN(record.value));unavailable++; }
      }
      checked++;
    }
    console.log(`Real pose electrostatics: ${checked} AA files with matching checksums/row counts; ${successful} successful values, ${unavailable} unavailable values kept missing.`);
  }
  console.log('Pose electrostatics passed: exact bundle/AA/pocket joins, pose identity, lazy cache, validation/failure isolation, zero/missing handling, sortable best/all pockets, unchanged rankings and safe full-precision exports.');
})().catch(error=>{console.error(error);process.exitCode=1;});
