// Point-cloud validation, lazy download and serialized-layer checks; no source-data writes.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const zlib = require('node:zlib');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, { textContent: '', hidden: false, checked: true, dataset: {}, handlers: {},
    addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
}
const warnings = [];
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) },
  document: { querySelector: node, querySelectorAll: () => [] }, Blob, Response, DecompressionStream, TextDecoder });
for (const file of ['ligand-viewer.js', 'pocket-electrostatics.js', 'pocket-cloud.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
const coords = '1.125,-2.500,3.000;2.125,-2.500,3.000';
context.coords = coords;
const table = `pocket_id\tpoints\r\n1\t${coords}\r\n2\t\r\n3\tNaN,1,2\r\n`;
context.table = table;
assert.equal(run('parsePocketPointTable(table).size'), 3);
for (const value of ['wrong\tpoints\n1\t', 'pocket_id\tpoints\n1\t\n1\t', 'pocket_id\tpoints\n0\t', 'pocket_id\tpoints\n1\tx\textra']) {
  context.bad = value; assert.throws(() => run('parsePocketPointTable(bad)'));
}
for (const value of ['1,2', 'NaN,1,2', '1,,2', '1,Infinity,2', '999999,0,0', '1,2,3;', Array(10000).fill('1,2,3').join(';')]) {
  context.bad = value; assert.throws(() => run('parsePocketPoints(bad)'));
}
const pdb = run('pocketPointsPDB(parsePocketPoints(coords))');
pdb.split('\n').filter(line => line.startsWith('HETATM')).forEach((line, index) => {
  assert.equal(line.length, 80);
  assert.equal(line.slice(17, 20), 'PNT'); assert.equal(line.slice(21, 22), 'Q');
  assert.deepEqual([30, 38, 46].map(start => Number(line.slice(start, start + 8))), coords.split(';')[index].split(',').map(Number));
});
assert.ok(!pdb.includes('CONECT'), 'Points have no bonds');
run(`
function requestFor(id='1', aa='ALA') {
  const request={protein:'P1',aa,row:{pocket_id:id,pocket:'pocket'+id,protein:'AF-P1-F1-model_v6',
    _compactDirectory:'At_results/L',_compactAA:aa,_compactPocketHash:'a'.repeat(64),sas_points:3,
    center_x:1,center_y:-2,center_z:3}};
  Object.assign(state,{selectedProtein:'P1',aa,currentView:'protein',showPocketCloud:true});
  molstarLatestRequest=request;
  return request;
}
`);
for (const change of ["r.row._compactDirectory='At_results/D'", "r.aa='DALA'", "r.row._compactAA='LEU'", "delete r.row._compactPocketHash"]) {
  assert.equal(run(`(()=>{const r=requestFor();${change};return pocketPointSource(r)})()`), null);
}
const fetches = [], representations = [], additions = [], deleted = [], camera = [];
const cells = new Set(['protein', 'sticks', 'ligand']);
let fetchGate = null, rawGate = null, httpError = false, failRepr = false, serial = 0;
context.fetch = async url => {
  fetches.push(url);
  if (fetchGate) await fetchGate;
  return httpError ? new Response('', { status: 404 }) : new Response(zlib.gzipSync(table));
};
const viewer = { plugin: {
  build: () => ({ delete: ref => ({ commit: async () => {
    assert.ok(!['protein', 'sticks', 'ligand'].includes(ref)); cells.delete(ref); deleted.push(ref);
  } }) }),
  builders: {
    data: { rawData: async params => {
      const ref = 'cloud-' + ++serial; cells.add(ref); additions.push(params);
      if (rawGate) await rawGate;
      return { ref, count: params.data.split('\n').filter(line => line.startsWith('HETATM')).length };
    } },
    structure: {
      parseTrajectory: async (data, format) => { assert.equal(format, 'pdb'); return data; },
      createModel: async data => data,
      createStructure: async (data, options) => { assert.equal(options.name, 'model'); return { obj: { data: { elementCount: data.count } } }; },
      representation: { addRepresentation: async (_, options) => {
        if (failRepr) throw Error('Intentional render failure');
        representations.push(options);
        return { ref: 'cloud-representation' };
      } },
    },
  },
  managers: { camera: { focusSpheres: spheres => camera.push(spheres) } },
} };
context.viewer = viewer;
async function show(expression = 'requestFor()') {
  await run(`prepareMolstarPocketCloud(viewer,${expression})`);
  await run('molstarUpdateQueue');
}
(async () => {
  assert.equal(fetches.length, 0);
  await show();
  assert.equal(fetches.length, 1);
  assert.equal(node('#structure-cloud-status').dataset.state, 'ready');
  assert.ok(node('#structure-cloud-status').textContent.includes('2 exported points (metadata: 3)'));
  assert.equal(cells.size, 4);
  assert.equal(additions[0].data, pdb);
  const repr = representations[0];
  assert.equal(repr.type, 'spacefill'); assert.equal(repr.typeParams.alpha, 0.3);
  assert.equal(repr.colorParams.value, 0xffd700); assert.equal(repr.sizeParams.value, 0.55);
  assert.equal(repr.typeParams.ignoreHydrogens, false, 'Pseudoatom carrier must not be hidden as hydrogen');
  assert.deepEqual(Array.from(repr.typeParams.visuals), ['element-sphere']);
  await show('requestFor("1","MET")');
  assert.equal(fetches.length, 1, 'The point table is shared across AAs');
  assert.equal(cells.size, 4, 'Layers do not accumulate');
  await show('requestFor("2")');
  assert.equal(node('#structure-cloud-status').dataset.state, 'missing'); assert.equal(cells.size, 3);
  await show('requestFor("3")');
  assert.equal(node('#structure-cloud-status').dataset.state, 'error');
  await show('requestFor("1","DALA")');
  assert.equal(fetches.length, 1, 'D IDs cannot load the L cloud');
  await show('(()=>{const r=requestFor();state.showPocketCloud=false;return r})()');
  assert.equal(node('#structure-cloud-status').dataset.state, 'hidden'); assert.equal(fetches.length, 1);

  // A long shared download does not occupy the structure queue, and old requests cannot draw.
  run('pocketPointTableCache=null');
  let releaseFetch;
  fetchGate = new Promise(resolve => { releaseFetch = resolve; });
  const pending = run('prepareMolstarPocketCloud(viewer,requestFor())');
  let independentUpdate = false;
  await run('molstarUpdateQueue').then(() => { independentUpdate = true; });
  assert.ok(independentUpdate);
  const before = additions.length;
  run('requestFor("2")');
  releaseFetch(); fetchGate = null;
  await pending; await run('molstarUpdateQueue');
  assert.equal(additions.length, before);

  run('pocketPointTableCache=null'); httpError = true;
  await show(); assert.equal(node('#structure-cloud-status').dataset.state, 'error');
  const count = fetches.length; httpError = false;
  await show(); assert.equal(fetches.length, count + 1, 'Failed download can retry');
  failRepr = true;
  await show(); assert.equal(cells.size, 3, 'Failed rendering removes only its cloud subtree');
  failRepr = false;

  let releaseRaw;
  rawGate = new Promise(resolve => { releaseRaw = resolve; });
  await run('prepareMolstarPocketCloud(viewer,requestFor())');
  const drawing = run('molstarUpdateQueue');
  await new Promise(resolve => setImmediate(resolve));
  run('molstarLatestRequest=null');
  releaseRaw(); rawGate = null;
  await drawing; assert.equal(cells.size, 3, 'Cancelled render cleans up partially created data');
  await show();
  run('molstarViewer=viewer;bindPocketCloudEvents()');
  node('#show-pocket-cloud').handlers.change({ target: { checked: false } });
  await run('molstarUpdateQueue');
  assert.equal(cells.size, 3, 'Cloud checkbox preserves protein, pocket sticks and ligand');
  assert.equal(node('#structure-cloud-legend').hidden, true);
  assert.ok(deleted.length > 0 && camera.length > 0);

  if (process.argv.includes('--real-data')) {
    const gz = file => zlib.gunzipSync(fs.readFileSync(path.join(root, file))).toString('utf8');
    context.realTable = gz('At_results/pocket_points.tsv.gz');
    const records = run('parsePocketPointTable(realTable)');
    const pocketRows = gz('At_results/L/pockets.tsv.gz').trimEnd().split(/\r?\n/);
    const header = pocketRows.shift().split('\t');
    const ids = new Map(pocketRows.map(line => { const fields = line.split('\t'); return [fields[0], Number(fields[header.indexOf('sas_points')])]; }));
    let total = 0, countDifferences = 0;
    for (const [id, value] of records) {
      assert.ok(ids.has(id)); context.realPoints = value;
      const points = run('parsePocketPoints(realPoints)'); total += points.length;
      if (points.length !== ids.get(id)) countDifferences++;
      if (Number(id) % 101 === 1) {
        const exportPDB = run('pocketPointsPDB(parsePocketPoints(realPoints))');
        const lines = exportPDB.split('\n').filter(line => line.startsWith('HETATM'));
        assert.equal(lines.length, points.length);
        assert.equal(new Set(lines.map(line => line.slice(12, 16))).size, points.length);
        lines.forEach((line, i) => assert.deepEqual([30, 38, 46].map(start => Number(line.slice(start, start + 8))), Array.from(points[i])));
      }
    }
    console.log(`Real cloud data: ${records.size} matched pocket IDs, ${total} valid points, ${countDifferences} explicitly reported metadata-count differences.`);
  }
  console.log('Pocket cloud passed: exact positions, no sampling/bonds, transparency, shared lazy gzip cache, nonblocking loading, missing/errors/retry, cancellation and isolated visibility.');
})().catch(error => { console.error(error); process.exitCode = 1; });
