// Read-only fixtures for potential provenance, lazy loading and isolated Mol* color updates.
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
  if (!nodes.has(selector)) nodes.set(selector, { value: '', hidden: false, textContent: '', dataset: {}, handlers: {},
    addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
};
const warnings = [];
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) },
  document: { querySelector: node, querySelectorAll: () => [] },
  crypto: webcrypto, TextEncoder, TextDecoder, Blob, Response, DecompressionStream });
for (const file of ['ligand-viewer.js', 'pocket-electrostatics.js', 'pocket-cloud.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
const model = 'AF-P1-F1-model_v6';
const points = [[0, -0, 1.125], [2, 3, 4], [5, 6, 7]];
const fingerprint = hash('0.000,0.000,1.125;2.000,3.000,4.000;5.000,6.000,7.000');
const compact = Buffer.from(JSON.stringify({ format: 'aa-proteome-atlas-compact', version: 1, pockets: { sha256: 'a'.repeat(64) } }));
const scale = { minimum: -5, midpoint: 0, maximum: 5, colors: ['#d73027', '#ffffff', '#4575b4'] };
const manifest = {
  format: 'aa-pocket-electrostatics', version: 1, bundles: [{ manifest_sha256: hash(compact) }], proteins: {},
  points: { columns: ['pocket_rank', 'point_index', 'phi_kT_per_e'], join: ['protein', 'pocket_rank'], point_index_base: 0,
    order: 'ATOM/HETATM encounter order within the original P2Rank pocket rank',
    coordinate_fingerprint: 'SHA256 of ASCII x,y,z;x,y,z;...; fixed 3 decimals, -0.000 normalized to 0.000; no newline',
    file_template: 'points_{protein}.tsv.gz', potential_units: 'kT/e', color_scale: scale },
  summary: { file: 'pocket_summary.tsv.gz', columns: ['protein', 'pocket_rank', 'n_points', 'points_sha256',
    'phi_mean', 'phi_min', 'phi_max', 'fraction_positive', 'fraction_negative', 'status', 'error'] },
};
let files, fetches, gate, httpFailure;
function resetFiles({ status = 'success', potential = '-10', fingerprintValue = fingerprint, missing = false } = {}) {
  const summary = zlib.gzipSync(manifest.summary.columns.join('\t') + '\n'
    + [model, '1', '3', fingerprintValue, '', '', '', '', '', status, status === 'success' ? '' : 'Receptor preparation failed'].join('\t') + '\n'
    + [model, '2', '3', fingerprintValue, '', '', '', '', '', 'success', ''].join('\t') + '\n');
  const potentials = zlib.gzipSync('pocket_rank\tpoint_index\tphi_kT_per_e\n'
    + `1\t0\t${potential}\n1\t1\t${missing ? '' : '0'}\n1\t2\t8\n2\t0\t1\n2\t1\t2\n2\t2\t3\n`);
  manifest.summary.sha256 = hash(summary);
  manifest.proteins[model] = { points_file: `points_${model}.tsv.gz`, points_file_sha256: hash(potentials) };
  files = new Map([
    ['At_results/electrostatics/manifest.json.gz', zlib.gzipSync(JSON.stringify(manifest))],
    ['At_results/electrostatics/pocket_summary.tsv.gz', summary],
    [`At_results/electrostatics/points_${model}.tsv.gz`, potentials],
    ['At_results/L/manifest.json', compact],
  ]);
  fetches = []; gate = null; httpFailure = false;
  run('electrostaticsManifestPromise = null; pocketElectrostaticsMetadata = null; pocketPotentialFiles.clear();');
}
context.points = points;
context.scale = scale;
context.request = { protein: 'P1', aa: 'ALA', row: { uniprot_id: 'P1', protein: model, pocket: 'pocket1', rank: 1,
  pocket_id: '1', _compactDirectory: 'At_results/L', _compactAA: 'ALA', _compactPocketHash: 'a'.repeat(64) } };
context.fetch = async url => {
  const file = url.split('?')[0];
  fetches.push(file);
  if (file.includes('/points_') && gate) await gate;
  return new Response(files.get(file) || '', { status: !httpFailure && files.has(file) ? 200 : 404 });
};
const load = () => run('loadValidatedPocketPotentials(request, points)');
const failLoad = async pattern => { await assert.rejects(load(), pattern); };
const waitFor = async check => {
  for (let i = 0; i < 1000; i++) {
    if (check()) return;
    await new Promise(resolve => setTimeout(resolve, 1));
  }
  throw new Error('Timed out waiting for test event');
};

(async () => {
  resetFiles();
  assert.equal(fetches.length, 0, 'No preload at startup');
  const result = await load();
  assert.deepEqual(Array.from(result.values), [-10, 0, 8]);
  assert.equal(result.missing, 0);
  assert.equal(fetches.length, 4, 'Only shared metadata, L manifest and selected model points');
  await load();
  context.request.row.rank = 2; context.request.row.pocket = 'pocket2';
  assert.deepEqual(Array.from((await load()).values), [1, 2, 3]);
  assert.equal(fetches.length, 4, 'Same protein and another pocket reuse the point file');
  context.request.row.rank = 1; context.request.row.pocket = 'pocket1';
  assert.equal(run('pocketCoordinateFingerprintText([[0,-0,1.125],[-0.0001,2,3]])'), '0.000,0.000,1.125;0.000,2.000,3.000');
  for (const [value, expected] of [[-10, 0xd73027], [-5, 0xd73027], [0, 0xffffff], [5, 0x4575b4], [10, 0x4575b4],
    [null, 0x9aa0a6], [NaN, 0x9aa0a6]]) {
    context.value = value;
    assert.equal(run('pocketPotentialColor(value, scale)'), expected);
  }
  context.value = -2.5;
  assert.equal(run('pocketPotentialColor(value, scale)'), 0xeb9893, 'Linear red-white interpolation');

  // Validation rejects wrong export/model/rank/count/coordinate order before using any colors.
  context.request.row._compactPocketHash = 'b'.repeat(64); await failLoad(/Pocket export checksum/);
  context.request.row._compactPocketHash = 'a'.repeat(64);
  context.request.row.uniprot_id = 'OTHER'; await failLoad(/supported, versioned/);
  context.request.row.uniprot_id = 'P1';
  context.request.row.pocket = 'pocket2'; await failLoad(/model\/rank/);
  context.request.row.pocket = 'pocket1';
  context.request.row._compactDirectory = 'At_results/D'; await failLoad(/model\/rank/);
  context.request.row._compactDirectory = 'At_results/L';
  context.points = points.slice(1); await failLoad(/point count/);
  context.points = [...points].reverse(); await failLoad(/coordinate\/order checksum/);
  context.points = points;
  resetFiles({ fingerprintValue: 'b'.repeat(64) }); await failLoad(/coordinate\/order checksum/);
  assert.ok(!fetches.some(file => file.includes('/points_')));
  resetFiles({ status: 'receptor_error' }); await failLoad(/receptor_error.*Receptor preparation/);
  assert.ok(!fetches.some(file => file.includes('/points_')), 'Failed calculations do not load/paint samples');
  resetFiles({ missing: true });
  assert.deepEqual(Array.from((await load()).values), [-10, null, 8]);
  assert.equal((await load()).missing, 1, 'Blank is missing, not neutral zero');

  for (const text of ['1\t0\tNaN', '1\t0\tInfinity', '1\t1\t3', '1\t0\t1\n1\t0\t2', '1\t0\t1\n1\t2\t2',
    '1\t0\t0x12', '0\t0\t1', '1\t0\t1\textra']) {
    context.badText = 'pocket_rank\tpoint_index\tphi_kT_per_e\n' + text + '\n';
    assert.throws(() => run('parsePocketPotentials(badText)'));
  }
  context.badManifest = { ...manifest, points: { ...manifest.points, point_index_base: 1 } };
  assert.throws(() => run('validatePocketPotentialManifest(badManifest)'), /Unsupported/);
  context.badManifest = { ...manifest, points: { ...manifest.points, color_scale: { ...scale, maximum: 10 } } };
  assert.throws(() => run('validatePocketPotentialManifest(badManifest)'), /color scale/);

  resetFiles();
  const pointPath = `At_results/electrostatics/points_${model}.tsv.gz`;
  const goodBytes = files.get(pointPath);
  files.set(pointPath, zlib.gunzipSync(goodBytes));
  await failLoad(/checksum mismatch/, 'HTTP transparent decompression must not bypass compressed-byte checksums');
  files.set(pointPath, goodBytes);
  await load();
  assert.equal(fetches.filter(file => file === pointPath).length, 2, 'Failed point-file checks can retry');
  resetFiles();
  files.set('At_results/L/manifest.json', Buffer.from(compact.toString() + '\n'));
  await failLoad(/different compact L bundle/);
  resetFiles(); files.set('At_results/electrostatics/pocket_summary.tsv.gz', zlib.gzipSync('bad'));
  await failLoad(/checksum mismatch/);
  resetFiles(); context.crypto = undefined;
  await failLoad(/HTTPS or localhost/); assert.equal(fetches.length, 0);
  context.crypto = webcrypto;
  resetFiles(); httpFailure = true; await failLoad(/HTTP 404/); httpFailure = false; await load();

  // A bounded LRU caches only models specifically requested; never fetch a proteome directory.
  resetFiles();
  for (const id of ['P2', 'P3', 'P4']) {
    const other = `AF-${id}-F1-model_v6`;
    const entry = { points_file: `points_${other}.tsv.gz`, points_file_sha256: hash(files.get(pointPath)) };
    files.set(`At_results/electrostatics/${entry.points_file}`, files.get(pointPath));
    context.otherModel = other; context.entry = entry;
    await run('loadProteinPocketPotentials(otherModel, entry)');
  }
  assert.equal(run('pocketPotentialFiles.size'), 2);
  assert.equal(fetches.length, 3);

  // Mock a reordered Mol* point carrier; only its existing representation may be updated.
  const unit = { elements: [0, 1, 2], model: { atomicHierarchy: { atoms: {
    label_atom_id: { value: element => ['P002', 'P000', 'P001'][element] },
  } } } };
  const structureRoot = { units: [unit], elementCount: 3 };
  structureRoot.root = structureRoot;
  const structure = { obj: { data: structureRoot } };
  const repr = { ref: 'only-cloud', params: { colorTheme: { name: 'uniform', params: { value: 0xffd700 } },
    type: { name: 'spacefill', params: { alpha: .3, sizeFactor: 1 } }, sizeTheme: { name: 'uniform', params: { value: .55 } } } };
  const providers = [], updates = [];
  const viewer = { plugin: {
    representation: { structure: { themes: { colorThemeRegistry: { add: provider => providers.push(provider) } } } },
    get builders() { throw new Error('Must not rebuild geometry'); },
    get managers() { throw new Error('Must not touch the camera'); },
    build() {
      let callback;
      const builder = { to(target) { assert.equal(target, repr, 'Protein/sticks/ligand are never targeted'); return builder; },
        update(fn) { callback = fn; return builder; },
        async commit(options) { assert.equal(options.doNotUpdateCurrent, true); callback(repr.params); updates.push(repr.params.colorTheme); } };
      return builder;
    },
  } };
  Object.assign(context, { viewer, structure, repr, unit });
  const setCloud = () => run(`
    state.currentView='protein'; state.selectedProtein='P1'; state.aa='ALA';
    state.showPocketCloud=true; state.pocketPointColor='gold'; molstarLatestRequest=request;
    setPocketCloudStatus('ready','3 exported points');
    registerPocketPointColorCloud(viewer,request,pocketCloudGeneration,points,structure,repr);
  `);
  resetFiles(); setCloud();
  assert.equal(fetches.length, 0, 'Gold does not load electrostatics');
  await run("state.pocketPointColor='potential'; refreshPocketPointColor()");
  assert.equal(repr.params.colorTheme.name, 'atlas-pocket-electrostatics');
  assert.equal(node('#pocket-potential-legend').hidden, false);
  assert.equal(node('#pocket-potential-missing-key').hidden, true);
  assert.equal(node('#structure-cloud-legend').hidden, true);
  assert.equal(providers.length, 1);
  const theme = providers[0].factory({ structure: structureRoot }, {});
  assert.equal(theme.granularity, 'group');
  assert.deepEqual([0, 1, 2].map(element => theme.color({ kind: 'element-location', unit, element })), [0x4575b4, 0xd73027, 0xffffff]);
  assert.deepEqual(repr.params.type, { name: 'spacefill', params: { alpha: .3, sizeFactor: 1 } });
  assert.deepEqual(repr.params.sizeTheme, { name: 'uniform', params: { value: .55 } });
  await run("state.pocketPointColor='gold'; refreshPocketPointColor()");
  assert.equal(repr.params.colorTheme.params.value, 0xffd700);
  assert.equal(node('#structure-cloud-legend').hidden, false);
  assert.equal(node('#pocket-potential-legend').hidden, true);
  await run("state.pocketPointColor='potential'; refreshPocketPointColor()");
  assert.equal(providers.length, 1, 'One provider per viewer');
  assert.equal(fetches.length, 4, 'Color toggles reuse validated files');

  resetFiles({ missing: true });
  await run('refreshPocketPointColor()');
  assert.equal(node('#pocket-potential-status').dataset.state, 'partial');
  assert.equal(node('#pocket-potential-missing-key').hidden, false);
  resetFiles({ status: 'coverage_error' });
  await run('refreshPocketPointColor()');
  assert.equal(node('#pocket-potential-status').dataset.state, 'unavailable');
  assert.match(node('#pocket-potential-status').textContent, /coverage_error/);
  assert.equal(repr.params.colorTheme.params.value, 0x9aa0a6);
  assert.equal(node('#pocket-potential-legend').hidden, true);
  assert.equal(node('#pocket-potential-retry').hidden, false);

  // Rapid Gold switching and pocket changes invalidate delayed potential colors.
  resetFiles();
  let release;
  gate = new Promise(resolve => { release = resolve; });
  const pending = run('refreshPocketPointColor()');
  await waitFor(() => fetches.includes(pointPath));
  await run('molstarUpdateQueue'); // Resolves even while the point download is blocked.
  await run("state.pocketPointColor='gold'; refreshPocketPointColor()");
  const beforeRelease = updates.length;
  release(); gate = null; await pending;
  assert.equal(updates.length, beforeRelease, 'Late potential cannot override Gold');
  assert.equal(repr.params.colorTheme.params.value, 0xffd700);
  resetFiles(); gate = new Promise(resolve => { release = resolve; });
  const stale = run("state.pocketPointColor='potential'; refreshPocketPointColor()");
  await waitFor(() => fetches.includes(pointPath));
  run('molstarLatestRequest={...request}; resetPocketPointColorCloud();');
  const beforePocketChange = updates.length;
  release(); gate = null; await stale;
  assert.equal(updates.length, beforePocketChange, 'Late potential cannot paint a different pocket');
  setCloud();
  run('state.showPocketCloud=false; bindPocketPointColorEvents();');
  const beforeHidden = fetches.length;
  node('#pocket-point-color').handlers.change({ target: { value: 'potential' } });
  assert.equal(fetches.length, beforeHidden, 'Hidden cloud does not initiate downloads');

  // Optional integration check with the actual compressed v1 files. No data writes.
  if (process.argv.includes('--real-data')) {
    run('electrostaticsManifestPromise=null; pocketElectrostaticsMetadata=null; pocketPotentialFiles.clear();');
    fetches = [];
    context.fetch = async url => {
      fetches.push(url.split('?')[0]);
      const local = path.resolve(root, url.split('?')[0]);
      assert.ok(local.startsWith(root + path.sep));
      return new Response(fs.readFileSync(local));
    };
    const gunzip = file => zlib.gunzipSync(fs.readFileSync(path.join(root, file))).toString();
    const meta = await run('loadPocketElectrostaticsMetadata()');
    const compactManifest = JSON.parse(fs.readFileSync(path.join(root, 'At_results/L/manifest.json')));
    const lines = gunzip('At_results/L/pockets.tsv.gz').trimEnd().split(/\r?\n/);
    const header = lines.shift().split('\t');
    const rows = lines.map(line => Object.fromEntries(line.split('\t').map((value, i) => [header[i], value])));
    context.coordinateText = gunzip('At_results/pocket_points.tsv.gz');
    const coordinateTable = run('parsePocketPointTable(coordinateText)');
    let checked = 0, failed = 0;
    const seen = new Set();
    for (const row of rows) {
      const summary = meta.summary.get(`${row.protein}|${row.rank}`);
      if (summary.status !== 'success' ? failed >= 2 : checked >= 6 || seen.has(row.protein)) continue;
      context.request = { protein: row.uniprot_id, aa: 'ALA', row: { ...row, rank: Number(row.rank),
        _compactDirectory: 'At_results/L', _compactAA: 'ALA', _compactPocketHash: compactManifest.pockets.sha256 } };
      context.pointText = coordinateTable.get(row.pocket_id);
      context.points = run('parsePocketPoints(pointText)');
      if (summary.status !== 'success') { await failLoad(/Calculation unavailable/); failed++; }
      else { assert.equal((await load()).values.length, context.points.length); checked++; seen.add(row.protein); }
      if (checked >= 6 && failed >= 2) break;
    }
    assert.equal(checked, 6); assert.equal(failed, 2);
    assert.equal(meta.summary.size, rows.length);
    assert.equal(fetches.filter(file => file.includes('/points_')).length, 6);
    console.log(`Real electrostatics: ${meta.summary.size} summary rows, ${checked} models with exact coordinate/order and compressed-file checksums; ${failed} failed pocket results rejected.`);
  }
  console.log('Pocket electrostatics passed: gzip/checksums, bundle/model/pocket/order validation, missing vs zero, fixed scale, lazy bounded caching, isolated recoloring and stale-request guards.');
})().catch(error => { console.error(error); process.exitCode = 1; });
