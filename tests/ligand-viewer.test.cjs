// Dependency-free coordinate, lifecycle and lazy-loading checks. No result-file writes.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const zlib = require('node:zlib');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, { textContent: '', hidden: false, checked: true,
    dataset: {}, handlers: {}, addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
}
const warnings = [];
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) },
  document: { querySelector: node, querySelectorAll: () => [] }, Blob, Response, TextDecoder, DecompressionStream });
for (const file of ['ligand-viewer.js', 'compact-data.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
const plain = value => JSON.parse(JSON.stringify(value));
const positions = 'C:10.125,-2.750,3.000;N:11.500,-2.700,3.000;H:12.000,-2.000,3.000';
const selectedPositions = 'C:-3.645,0.088,7.896;S:-4.940,-1.152,8.264';
context.positionText = positions;
context.selectedPositionText = selectedPositions;
assert.equal(run('parseLigandPositionTable("pocket_id\\tpose_id\\tpositions\\r\\n1\\t1\\t" + positionText + "\\r\\n1\\t4\\t\\r\\n").size'), 2);
assert.equal(run('parseLigandPositionTable("\\uFEFFpocket_id\\tpose_id\\tpositions\\n1\\t4\\t" + positionText).get("1:4")'), positions);
assert.equal(run('parseLigandAtoms("").length'), 0);
for (const bad of ['pocket_id\tother\n1\tX', 'pocket_id\tpositions\n1\tC:1,2,3',
  'pocket_id\tpose_id\tpositions\n1\t1\t\n1\t1\t', 'pocket_id\tpose_id\tpositions\n0\t1\t',
  'pocket_id\tpose_id\tpositions\n1\t0\t', 'pocket_id\tpose_id\tpositions\n1\t-1\t',
  'pocket_id\tpose_id\tpositions\n1\t1.5\t', 'pocket_id\tpose_id\tpositions\n1\t01\t',
  'pocket_id\tpose_id\tpositions\n1\t9007199254740992\t', 'pocket_id\tpose_id\tpositions\n1\t1',
  'pocket_id\tpose_id\tpositions\n1\t1\tx\textra']) {
  context.badText = bad;
  assert.throws(() => run('parseLigandPositionTable(badText)'));
}
for (const bad of ['C:NaN,1,2', 'C:Infinity,1,2', 'C:,1,2', 'X:1,2,3', 'C:1,2', 'C:1,2,3;', 'C:100000,2,3', 'H:1,2,3']) {
  context.badText = bad;
  assert.throws(() => run('parseLigandAtoms(badText)'));
}
const pdb = run('ligandPDB(parseLigandAtoms(positionText))');
const atomLines = pdb.split('\n').filter(line => line.startsWith('HETATM'));
assert.equal(atomLines.length, 3);
atomLines.forEach((line, i) => {
  const [element, xyz] = positions.split(';')[i].split(':');
  assert.equal(line.length, 80);
  assert.equal(line.slice(17, 20), 'LIG');
  assert.equal(line.slice(21, 22), 'Z');
  assert.equal(line.slice(76, 78).trim(), element);
  assert.deepEqual([30, 38, 46].map(start => Number(line.slice(start, start + 8))), xyz.split(',').map(Number));
});

run(`
const testHash = 'a'.repeat(64);
function requestFor(aa='ALA', id='1', metric='vina_affinity') {
  const request={aa,protein:'P12345',ligandMetric:metric,row:{_compactDirectory:'At_results/L',_compactAA:aa,_compactPocketHash:testHash,
    uniprot_id:'P12345',protein:'AF-P12345-F1-model_v6',pocket:'pocket'+id,pocket_id:id,
    vina_status:'success',sfct_status:'success',sfct_best_pose:3,vina_affinity:-4,sfct_score:3,
    vina_sfct_combined:1.6,vina_sfct_combined_50:-.5,center_x:11,center_y:-3,center_z:3}};
  state.currentView='protein';state.selectedProtein=request.protein;state.aa=aa;state.showLigand=true;
  molstarLatestRequest=request;
  return request;
}
`);
assert.equal(run('ligandPositionSource(requestFor())'), `At_results/L/aa_positions/positions_ala.tsv.gz?v=${'a'.repeat(64)}&poses=2`);
for (const mutation of ["r.row._compactDirectory='At_results/D'", "r.row._compactAA='LEU'", "r.aa='DALA'", "delete r.row._compactPocketHash", "r.row.pocket_id='../x'"]) {
  assert.equal(run(`(()=>{const r=requestFor();${mutation};return ligandPositionSource(r)})()`), null);
}
assert.equal(run('ligandPoseSelection(requestFor()).poseId'), 1);
for (const metric of ['sfct_score', 'vina_sfct_combined', 'vina_sfct_combined_50']) {
  context.testMetric = metric;
  assert.equal(run('ligandPoseSelection(requestFor("ALA","1",testMetric)).poseId'), 4);
}
for (const index of [NaN, null, undefined, -1, 1.5, '3', Infinity, Number.MAX_SAFE_INTEGER]) {
  context.testIndex = index;
  assert.equal(run('savedSfctPoseId({sfct_best_pose:testIndex})'), null);
}
assert.equal(run('savedSfctPoseId({sfct_best_pose:0})'), 1);
assert.equal(run('savedSfctPoseId({sfct_best_pose:9})'), 10);
assert.equal(run('state.ligandPoseMode'), 'auto');
assert.equal(run('ligandPoseSelection(requestFor("ALA","1","sfct_score"),"vina").poseId'), 1);
assert.equal(run('ligandPoseSelection(requestFor(),"sfct").poseId'), 4);
assert.equal(run('(()=>{const r=requestFor("ALA","1","sfct_score");r.row.sfct_status="failed";return ligandPoseSelection(r,"vina").poseId})()'), 1);
assert.equal(run('(()=>{const r=requestFor();r.row.vina_status="failed";return ligandPoseSelection(r,"sfct").poseId})()'), 4);

const cells = new Set(['receptor']);
const added = [], deleted = [], focused = [], fetches = [];
let serial = 0, failingRepresentation = false, rawGate = null;
const viewer = { plugin: {
  build: () => ({ delete: ref => ({ commit: async () => { assert.notEqual(ref, 'receptor'); deleted.push(ref); cells.delete(ref); } }) }),
  builders: {
    data: { rawData: async params => {
      const ref = `ligand-${++serial}`;
      cells.add(ref); added.push({ ref, ...params });
      if (rawGate) await rawGate;
      return { ref };
    } },
    structure: {
      parseTrajectory: async (_, format) => { assert.equal(format, 'pdb'); return {}; },
      createModel: async () => ({}),
      createStructure: async (_, options) => { assert.equal(options.name, 'model'); return {}; },
      representation: { addRepresentation: async (_, options) => {
        assert.equal(options.type, 'ball-and-stick'); assert.equal(options.typeParams.ignoreHydrogens, true);
        assert.equal(options.color, 'element-symbol');
        assert.equal(options.colorParams.carbonColor.name, 'element-symbol');
        assert.equal(options.colorParams.lightness, 0);
        if (failingRepresentation) throw Error('Representation failed');
      } },
    },
  },
  managers: { camera: { focusSpheres: (spheres) => focused.push(plain(spheres)) } },
} };
context.testViewer = viewer;
let fetchFailure = false, fetchGate = null;
context.fetch = async url => {
  fetches.push(url);
  if (fetchGate) await fetchGate;
  if (fetchFailure) return new Response('', { status: 404 });
  return new Response(zlib.gzipSync(`pocket_id\tpose_id\tpositions\n1\t1\t${positions}\n1\t4\t${selectedPositions}\n2\t1\t\n3\t1\tC:NaN,1,2\n4\t1\t${positions}\n`));
};
async function display(expression = 'requestFor()') {
  return run(`updateMolstarLigand(testViewer,${expression})`);
}

(async () => {
  assert.equal(fetches.length, 0, 'No coordinate downloads at startup');
  await display();
  assert.equal(node('#structure-ligand-status').dataset.state, 'ready');
  assert.ok(node('#structure-ligand-status').textContent.includes('2 heavy atoms'));
  assert.equal(node('#structure-ligand-legend').hidden, false);
  assert.equal(fetches.length, 1);
  assert.equal(cells.size, 2);
  assert.equal(added[0].data, pdb);
  assert.ok(added[0].label.includes('Vina MODEL 1 · AutoDock Vina'));
  assert.equal(focused.length, 1);
  const selectedPDB = run('ligandPDB(parseLigandAtoms(selectedPositionText))');
  for (const metric of ['sfct_score', 'vina_sfct_combined', 'vina_sfct_combined_50']) {
    await display(`requestFor("ALA","1","${metric}")`);
    assert.equal(node('#structure-ligand-status').dataset.state, 'ready');
    assert.equal(added.at(-1).data, selectedPDB, 'Saved selection loads MODEL 4, not MODEL 3 or MODEL 1');
    assert.match(added.at(-1).label, /Vina MODEL 4/);
    assert.equal(cells.size, 2, 'Only one ligand subtree remains after changing score');
    assert.equal(fetches.length, 1, 'Both poses share the AA coordinate cache');
  }
  await display('requestFor()');
  assert.equal(added.at(-1).data, pdb, 'Returning to Vina restores MODEL 1 unchanged');
  await display('(()=>{const r=requestFor("ALA","1","sfct_score");r.row.sfct_best_pose=0;return r})()');
  assert.equal(added.at(-1).data, pdb, 'Shared MODEL 1 requires only one exported row');
  await display('requestFor("ALA","4","sfct_score")');
  assert.equal(node('#structure-ligand-status').dataset.state, 'missing');
  assert.match(node('#structure-ligand-status').textContent, /MODEL 4 coordinates unavailable; no other pose substituted/);
  assert.equal(cells.size, 1, 'Absent selected MODEL 4 must never fall back to existing MODEL 1');
  await display('(()=>{const r=requestFor("ALA","1","sfct_score");r.row.vina_status="failed";return r})()');
  assert.equal(node('#structure-ligand-status').dataset.state, 'ready', 'Successful saved SFCT pose is independent of Vina score status');
  for (const mutation of ['r.row.sfct_status="failed"', 'r.row.sfct_score=NaN', 'delete r.row.sfct_best_pose', 'r.row.sfct_best_pose=-1', 'r.row.sfct_best_pose=1.5']) {
    const count = added.length;
    await display(`(()=>{const r=requestFor("ALA","1","sfct_score");${mutation};return r})()`);
    assert.equal(node('#structure-ligand-status').dataset.state, 'missing');
    assert.equal(added.length, count, 'Invalid saved selection is not rounded or substituted');
    assert.equal(cells.size, 1);
  }
  await display('requestFor("ALA","2")');
  assert.equal(node('#structure-ligand-status').dataset.state, 'missing');
  assert.equal(cells.size, 1, 'Empty pocket removes the old ligand, not the receptor');
  assert.equal(fetches.length, 1, 'Other pockets reuse the AA coordinate cache');
  await display('requestFor("ALA","3")');
  assert.equal(node('#structure-ligand-status').dataset.state, 'error');
  assert.equal(cells.size, 1);
  await display('(()=>{const r=requestFor();r.row.vina_status="failed";return r})()');
  assert.equal(fetches.length, 1);
  assert.equal(node('#structure-ligand-status').dataset.state, 'missing');
  await display('requestFor("DALA")');
  assert.equal(fetches.length, 1, 'Never fetch D coordinates from the L bundle');
  await display('(()=>{const r=requestFor();state.showLigand=false;return r})()');
  assert.equal(fetches.length, 1, 'Hidden ligand makes no coordinate request');
  assert.equal(node('#structure-ligand-status').dataset.state, 'hidden');

  await display('requestFor("LEU")');
  await display('requestFor("MET")');
  assert.equal(run('ligandPositionCache.size'), 2);
  const beforeEvicted = fetches.length;
  await display();
  assert.equal(fetches.length, beforeEvicted + 1, 'Evicted AA is fetched again');
  failingRepresentation = true;
  await display();
  assert.equal(cells.size, 1, 'Partial ligand subtree removed after representation failure');
  assert.equal(node('#structure-ligand-status').dataset.state, 'error');
  failingRepresentation = false;

  fetchFailure = true;
  await display('requestFor("VAL")');
  assert.equal(node('#structure-ligand-status').dataset.state, 'error');
  const beforeRetry = fetches.length;
  fetchFailure = false;
  await display('requestFor("VAL")');
  assert.equal(fetches.length, beforeRetry + 1, 'Failed downloads can retry');
  assert.equal(node('#structure-ligand-status').dataset.state, 'ready');

  // Out-of-order requests must not add a pose after another protein/AA is selected.
  let releaseFetch;
  fetchGate = new Promise(resolve => { releaseFetch = resolve; });
  const pendingFetch = display('requestFor("SER")');
  await new Promise(resolve => setImmediate(resolve));
  const beforeStale = added.length;
  run('requestFor("THR")');
  releaseFetch(); fetchGate = null;
  await pendingFetch;
  assert.equal(added.length, beforeStale);
  assert.equal(cells.size, 1);

  let releaseRaw;
  rawGate = new Promise(resolve => { releaseRaw = resolve; });
  const pendingRaw = display('requestFor("ALA")');
  await new Promise(resolve => setImmediate(resolve));
  run('state.showLigand=false;molstarLatestRequest=null');
  releaseRaw(); rawGate = null;
  await pendingRaw;
  assert.equal(cells.size, 1, 'Cancelled raw-data creation cleans its subtree');

  // A score change for the same protein/AA invalidates a pending old pose too.
  rawGate = new Promise(resolve => { releaseRaw = resolve; });
  const pendingVina = display('requestFor("ALA","1","vina_affinity")');
  await new Promise(resolve => setImmediate(resolve));
  run('requestFor("ALA","1","sfct_score")');
  releaseRaw(); rawGate = null;
  await pendingVina;
  assert.equal(cells.size, 1, 'Stale MODEL 1 does not survive a same-AA score switch');
  await run('updateMolstarLigand(testViewer,molstarLatestRequest)');
  assert.equal(added.at(-1).data, selectedPDB);
  await display('(()=>{const r=requestFor();state.showLigand=false;return r})()');

  // Pose and visibility controls touch only the ligand, preserving the structure request and camera.
  run('molstarViewer=testViewer; molstarLoadedModel="AF-P12345-F1-model_v6"; requestFor("ALA","1","sfct_score"); bindLigandViewerEvents();');
  cells.add('pocket-cloud'); cells.add('pocket-sticks');
  const controlRequest = run('molstarLatestRequest');
  const changePose = value => node('#ligand-pose-select').handlers.change({ target: { value } });
  const showLigand = checked => node('#show-pocket-ligand').handlers.change({ target: { checked } });
  const flushOverlay = () => run('molstarUpdateQueue');
  const fixedFocusCount = focused.length;
  const fixedFetchCount = fetches.length;
  const fixedScores = plain(run('({metric:state.metric,value:state.profileValue,row:molstarLatestRequest.row})'));
  assert.equal(node('#ligand-pose-select').value, 'auto');
  assert.equal(node('#ligand-pose-select').disabled, false);
  assert.equal(node('#ligand-sfct-pose-option').textContent, 'Saved SFCT/Combined - MODEL 4');
  changePose('vina'); await flushOverlay();
  assert.equal(added.at(-1).data, pdb);
  assert.match(node('#structure-ligand-status').textContent, /MODEL 1.*manual/);
  assert.match(node('#structure-ligand-note').textContent, /Different pose from OnionNet-SFCT/);
  assert.match(node('#structure-ligand-note').textContent, /scores and plots are unchanged/);
  assert.match(node('#structure-ligand-note').textContent, /qφ \(kT\) remains Vina MODEL 1/);
  changePose('sfct'); await flushOverlay();
  assert.equal(added.at(-1).data, selectedPDB);
  assert.match(node('#structure-ligand-status').textContent, /MODEL 4.*manual/);
  assert.doesNotMatch(node('#structure-ligand-note').textContent, /Different pose/);
  changePose('auto'); await flushOverlay();
  assert.equal(added.at(-1).data, selectedPDB);
  assert.doesNotMatch(node('#structure-ligand-status').textContent, /manual/);
  assert.equal(focused.length, fixedFocusCount, 'Dropdown never resets the camera');
  assert.equal(fetches.length, fixedFetchCount, 'Manual poses reuse the same AA file');
  assert.equal(run('molstarLatestRequest'), controlRequest, 'Pending cloud/structure requests are not invalidated');
  assert.deepEqual(plain(run('({metric:state.metric,value:state.profileValue,row:molstarLatestRequest.row})')), fixedScores);
  assert.ok(cells.has('pocket-cloud') && cells.has('pocket-sticks'));
  assert.equal(cells.size, 4, 'Protein, cloud, sticks and one ligand remain');
  const validMode = run('state.ligandPoseMode');
  changePose('invalid'); await flushOverlay();
  assert.equal(run('state.ligandPoseMode'), validMode);

  // Fast dropdown changes must cancel stale pose creation even though the structure request is unchanged.
  rawGate = new Promise(resolve => { releaseRaw = resolve; });
  changePose('vina');
  await new Promise(resolve => setImmediate(resolve));
  const beforeRapid = added.length;
  changePose('sfct');
  releaseRaw(); rawGate = null;
  await flushOverlay();
  assert.equal(added.length, beforeRapid + 1);
  assert.equal(added.at(-1).data, selectedPDB);
  assert.equal(cells.size, 4);
  assert.equal(focused.length, fixedFocusCount);
  assert.equal(run('molstarLatestRequest'), controlRequest);

  showLigand(false); await flushOverlay();
  assert.equal(cells.size, 3, 'Hiding preserves the other three layers');
  assert.equal(node('#structure-ligand-status').dataset.state, 'hidden');
  const hiddenFetches = fetches.length, hiddenAdds = added.length;
  changePose('vina'); await flushOverlay();
  assert.equal(fetches.length, hiddenFetches);
  assert.equal(added.length, hiddenAdds, 'Changing selection while hidden does not draw or fetch a ligand');
  showLigand(true); await flushOverlay();
  assert.equal(added.at(-1).data, pdb);
  assert.equal(focused.length, fixedFocusCount, 'Visibility toggles also preserve the camera');

  // Labels and availability are specific to the current pocket, not a remembered model number.
  run('molstarLatestRequest.row.sfct_best_pose=0; syncLigandPoseControls();');
  assert.equal(node('#ligand-sfct-pose-option').textContent, 'Saved SFCT/Combined - MODEL 1');
  changePose('sfct'); await flushOverlay();
  assert.equal(added.at(-1).data, pdb, 'Both choices can legitimately select the same exported row');
  run('molstarLatestRequest.row.sfct_best_pose=8; syncLigandPoseControls();');
  assert.equal(node('#ligand-sfct-pose-option').textContent, 'Saved SFCT/Combined - MODEL 9');
  changePose('sfct'); await flushOverlay();
  assert.equal(node('#structure-ligand-status').dataset.state, 'missing');
  assert.match(node('#structure-ligand-status').textContent, /MODEL 9 coordinates unavailable/);
  assert.equal(cells.size, 3, 'Missing manual pose is not replaced with MODEL 1');
  run('molstarLatestRequest.row.sfct_status="failed"; syncLigandPoseControls();');
  assert.equal(node('#ligand-sfct-pose-option').disabled, true);
  assert.equal(node('#ligand-sfct-pose-option').textContent, 'Saved SFCT/Combined - unavailable');
  run('syncLigandPoseControls(null);');
  assert.equal(node('#ligand-pose-select').disabled, true);
  changePose('auto'); await flushOverlay();
  cells.delete('pocket-cloud'); cells.delete('pocket-sticks');

  // The real request builder follows Value type, not Explorer's independent metric.
  context.setPocketCloudStatus = () => {};
  run('performMolstarUpdate=async()=>{};');
  for (const [value, metric, model] of [
    ['raw_vina', 'vina_affinity', 1], ['raw_sfct', 'sfct_score', 4],
    ['raw_combined', 'vina_sfct_combined', 4], ['raw_combined_50', 'vina_sfct_combined_50', 4],
    ['percentile', 'vina_sfct_combined', 4], ['percentile_50', 'vina_sfct_combined_50', 4],
  ]) {
    for (const normalized of [false, true]) {
      run(`state.profileValue='${value}';state.profileNormalize=${normalized};state.metric='vina_sfct_combined_50';
        updateMolstarPocket(requestFor().row);`);
      await run('molstarUpdateQueue');
      assert.equal(run('molstarLatestRequest.ligandMetric'), metric);
      assert.equal(run('ligandPoseSelection(molstarLatestRequest).poseId'), model);
    }
  }

  // Validate compact provenance using the actual loader, with a tiny in-memory bundle.
  run(`readCompactTable=async()=>[{pocket_id:'1',vina_status:'success',sfct_status:'failed',vina_affinity:-4}];`);
  const compact = await run(`loadCompactResultRows({code:'ALA'},{byAA:new Map([['ALA',{
    directory:'At_results/L',pocketHash:testHash,ligands:new Map([['ALA',{}]]),
    pockets:new Map([['1',{uniprot_id:'P12345',pocket:'pocket1'}]])}]])})`);
  assert.equal(compact[0]._compactAA, 'ALA');
  assert.equal(compact[0]._compactDirectory, 'At_results/L');
  assert.equal(compact[0]._compactPocketHash, 'a'.repeat(64));
  assert.equal(compact[0].vina_affinity, -4);

  if (process.argv.includes('--real-data')) {
    const readGzip = file => zlib.gunzipSync(fs.readFileSync(path.join(root, file))).toString('utf8');
    const pocketIDs = new Set(readGzip('At_results/L/pockets.tsv.gz').trimEnd().split(/\r?\n/).slice(1).map(line => line.split('\t')[0]));
    let coordinateCount = 0, sampleAtoms = 0, missingVina = 0, missingSfct = 0, fileCount = 0;
    for (const file of fs.readdirSync(path.join(root, 'At_results/L/aa_positions')).filter(file => file.endsWith('.tsv.gz'))) {
      fileCount++;
      context.realText = readGzip(`At_results/L/aa_positions/${file}`);
      const table = run('parseLigandPositionTable(realText)');
      const scoreLines = readGzip(`At_results/L/${file.replace('positions_', 'scores_')}`).trimEnd().split(/\r?\n/);
      const headers = scoreLines.shift().split('\t');
      const scores = new Map(scoreLines.map(line => {
        const columns = line.split('\t');
        const row = Object.fromEntries(headers.map((key, i) => [key, columns[i]]));
        return [row.pocket_id, row];
      }));
      for (const [key, value] of table) {
        const [id, pose] = key.split(':');
        assert.ok(pocketIDs.has(id), `${file}/${key}`);
        const saved = scores.get(id);
        assert.ok(saved, `${file}/${key}: known score row`);
        assert.ok(pose === '1' || (saved.sfct_best_pose.trim() !== '' && Number(pose) === Number(saved.sfct_best_pose) + 1), `${file}/${key}: only MODEL 1 or the saved selection`);
        if (value) coordinateCount++;
      }
      for (const [id, score] of scores) {
        if (score.vina_status === 'success' && !table.get(`${id}:1`)) missingVina++;
        if (score.sfct_status === 'success' && !table.get(`${id}:${Number(score.sfct_best_pose) + 1}`)) missingSfct++;
      }
      if (['positions_cys.tsv.gz', 'positions_asp.tsv.gz'].includes(file)) {
        assert.ok(table.get('10460:1') && table.get('10460:4'), `${file}: Q38933 exports both MODEL 1 and MODEL 4`);
        assert.notEqual(table.get('10460:1'), table.get('10460:4'));
        context.realRow = Object.fromEntries(Object.entries(scores.get('10460')).map(([key, value]) => [key, key.endsWith('_status') ? value : Number(value)]));
        for (const [metric, expected] of [['vina_affinity', 1], ['sfct_score', 4], ['vina_sfct_combined', 4], ['vina_sfct_combined_50', 4]]) {
          context.realMetric = metric;
          const selected = run('ligandPoseSelection({row:addDerivedScores(realRow),ligandMetric:realMetric})');
          assert.equal(selected.poseId, expected);
          context.realPositions = table.get(`10460:${selected.poseId}`);
          const rendered = run('ligandPDB(parseLigandAtoms(realPositions))').split('\n').filter(line => line.startsWith('HETATM'));
          const coordinates = context.realPositions.split(';').map(token => token.split(':')[1].split(',').map(Number));
          rendered.forEach((line, i) => assert.deepEqual([30,38,46].map(start => Number(line.slice(start,start+8))), coordinates[i]));
        }
      }
      const samples = [...table.values()].filter(Boolean);
      for (const value of [...samples.slice(0, 10), ...samples.slice(-10)]) {
        context.realPositions = value;
        const atomCount = run('parseLigandAtoms(realPositions).length');
        assert.ok(atomCount > 0); sampleAtoms += atomCount;
        const samplePDB = run('ligandPDB(parseLigandAtoms(realPositions))');
        assert.equal(samplePDB.split('\n').filter(line => line.startsWith('HETATM')).length, atomCount);
      }
    }
    console.log(`Real coordinates: ${fileCount} AA files; ${coordinateCount} nonempty poses with matching pocket/pose IDs; ${sampleAtoms} sampled atoms parsed and converted. Successful scores without coordinates: Vina ${missingVina}, saved SFCT ${missingSfct}. Q38933 CYS/ASP select MODEL 1 for Vina and MODEL 4 for saved scores without changing coordinates.`);
  }
  assert.ok(deleted.length > 0 && cells.has('receptor'));
  console.log('Ligand viewer passed: automatic/manual pose selection, camera/layer/score preservation, rapid switches, pocket/pose keys, profile controls, coordinates/PDB, provenance, lazy cache/retry, missing poses and cleanup.');
})().catch(error => { console.error(error); process.exitCode = 1; });
