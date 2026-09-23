// Startup progress without network downloads, parsing real result files or browser dependencies.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(path.join(__dirname, '..', 'js', 'organisms.js'), 'utf8') + '\n' + fs.readFileSync(path.join(__dirname, '..', 'js', 'app.js'), 'utf8').replace(/\binit\(\);\s*$/, '');
function fixture(options = {}) {
  const nodes = new Map(), values = [], errors = [], loaded = [];
  let inFlight = 0, maximumInFlight = 0, hidden = false;
  const screen = { innerHTML: '', classList: { add(name) {
    assert.equal(name, 'hidden'); assert.equal(values.at(-1), 100); hidden = true;
  } } };
  nodes.set('#loading-screen', screen);
  nodes.set('#atlas-loading-progress', { set value(value) { values.push(value); } });
  function node(selector) {
    if (screen.innerHTML && selector.startsWith('#atlas-loading-')) return null;
    if (!nodes.has(selector)) nodes.set(selector, { textContent: '' });
    return nodes.get(selector);
  }
  const context = vm.createContext({ console: { ...console, error: error => errors.push(error) },
    document: { querySelector: node }, location: { protocol: options.file ? 'file:' : 'http:', hash: '' },
    requestAnimationFrame: fn => fn(),
    loadCompactData: async () => options.legacy ? null : { byAA: new Map((options.all
      ? vm.runInContext('DATA_LIGANDS.map(aa=>aa.code)', context) : ['ALA', 'MET', 'DALA']).map(code => [code, {}])) },
    loadCompactResultRows: readRows,
    async stubAnnotations() { if (options.failAnnotations) throw Error('Annotation failure'); },
    async stubRows(aa) { return readRows(aa); },
    stubRender() {
      assert.ok(values.at(-1) < 100, 'Do not claim completion until view preparation succeeds');
      if (options.failRender) throw Error('View failure');
    },
  });
  async function readRows(aa) {
    inFlight++; maximumInFlight = Math.max(maximumInFlight, inFlight);
    await Promise.resolve();
    inFlight--;
    if (aa.code === options.failAA) throw Error('Dataset <failure>');
    loaded.push(aa.code);
    return [{ uniprot_id: 'P-' + aa.code }];
  }
  vm.runInContext(source, context);
  vm.runInContext(`initializeOrganismUI=()=>{};readURLState=()=>{};populateAASelects=()=>{};syncControls=()=>{};bindEvents=()=>{};
    renderExplorer=stubRender;loadUniProtAnnotations=stubAnnotations;loadResultRows=stubRows;`, context);
  return { run: code => vm.runInContext(code, context), values, loaded, errors, nodes, screen,
    get hidden() { return hidden; }, get maximumInFlight() { return maximumInFlight; } };
}
(async () => {
  for (const options of [{}, { all: true }, { legacy: true }]) {
    const f = fixture(options);
    await f.run('init()');
    assert.equal(f.errors.length, 0);
    assert.equal(f.hidden, true);
    assert.equal(f.values[0], 0);
    assert.equal(f.values.at(-1), 100);
    assert.equal(f.values.filter(value => value === 100).length, 1);
    assert.ok(f.values.every((value, i, values) => value >= 0 && value <= 100 && (!i || value >= values[i - 1])));
    assert.equal(f.maximumInFlight, 2, 'Keep the existing two-dataset memory limit');
    assert.equal(f.run('state.rawByAA.size'), f.loaded.length);
    assert.equal(f.loaded.length, options.all ? 39 : options.legacy ? 21 : 3);
    assert.equal(f.nodes.get('#atlas-loading-percent').textContent, '100%');
  }
  const dataOnly = fixture();
  await dataOnly.run('loadData()');
  assert.deepEqual(dataOnly.values, [0, 17, 33, 50, 67, 83], 'Advance per AA, then annotations; leave a step for rendering');
  assert.equal(dataOnly.hidden, false);

  for (const options of [{ failAA: 'ALA' }, { failAnnotations: true }, { failRender: true }, { file: true }]) {
    const f = fixture(options);
    await f.run('init()');
    assert.equal(f.hidden, false);
    assert.equal(f.errors.length, 1);
    assert.ok(f.screen.innerHTML.includes('Atlas data could not load'));
    assert.ok(!f.values.includes(100));
    if (options.failAA) assert.ok(f.screen.innerHTML.includes('&lt;failure&gt;'), 'Keep escaped error details');
    assert.doesNotThrow(() => f.run('updateAtlasLoadingProgress(90)'), 'Late completions cannot replace or break the error screen');
  }
  const bounds = fixture();
  for (const value of ['-20', '130', 'NaN', 'Infinity', '12.4']) bounds.run(`updateAtlasLoadingProgress(${value})`);
  assert.deepEqual(bounds.values, [0, 100, 0, 0, 12]);
  console.log('Loading progress passed: compact/legacy populations, individual completions, bounded concurrency, annotations/view completion, monotonic percentages, error preservation and late updates.');
})().catch(error => { console.error(error); process.exitCode = 1; });
