// Protein-only style switching: preserve all other layers, camera and loaded coordinates.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const nodes = new Map();
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, { textContent: '', hidden: false, dataset: {},
    attributes: {}, setAttribute(name, value) { this.attributes[name] = value; } });
  return nodes.get(selector);
}
const buttons = ['cartoon', 'surface'].map(style => ({ dataset: { proteinStyle: style }, attributes: {}, handlers: {},
  setAttribute(name, value) { this.attributes[name] = value; }, addEventListener(type, handler) { this.handlers[type] = handler; } }));
const warnings = [], built = [], visibility = [], deleted = [];
const paints = new Map(), paintTransform = {}, proteinRoot = {};
let gate = null, fail = false, failPaint = false;
const viewer = { plugin: {
  state: { data: { updateCellState(ref, update) {
    assert.ok(/^(cartoon|surface)-/.test(ref), 'Only receptor representations can change visibility');
    visibility.push({ ref, ...update });
  } } },
  build: () => {
    const pending = [];
    return {
      delete: ref => ({ commit: async () => { deleted.push(ref); } }),
      to: representation => ({ applyOrUpdateTagged(tag, transform, params) {
        assert.ok(/^(cartoon|surface)-/.test(representation.ref));
        assert.equal(tag, 'atlas-pocket-color'); assert.equal(transform, paintTransform);
        pending.push([representation.ref, params.layers]);
      }, update() {} }),
      async commit() {
        if (failPaint) throw Error('Paint unavailable');
        for (const [ref, layers] of pending) paints.set(ref, layers);
      },
    };
  },
  builders: { structure: { representation: { async addRepresentation(component, options, params) {
    if (fail) throw Error('Surface unavailable');
    built.push({ component, options, params });
    if (gate) await gate;
    return { ref: 'surface-' + component.ref, obj: { data: { sourceData: { root: proteinRoot } } } };
  } } } },
  managers: { camera: { focusLoci() { throw Error('Camera must not move'); }, focusSpheres() { throw Error('Camera must not move'); } } },
} };
const context = vm.createContext({ console: { ...console, warn: (...args) => warnings.push(args) }, viewer,
  proteinRoot, window: { molstar: { lib: {
    plugin: { StateTransforms: { Representation: { OverpaintStructureRepresentation3DFromBundle: paintTransform } } },
    structure: { StructureElement: {
      Loci: { fromSchema: (root, schema) => { assert.equal(root, proteinRoot); return { root, schema }; }, isEmpty: () => false },
      Bundle: { fromLoci: loci => ({ schema: loci.schema }) },
    } },
  } } },
  document: { querySelector: node, querySelectorAll: selector => selector === '[data-protein-style]' ? buttons : [] },
  fetch() { throw Error('Switching style must not download data'); } });
const root = path.resolve(__dirname, '..');
for (const file of ['ligand-viewer.js', 'protein-viewer.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
run(`
function loadFixture(id='1') {
  resetProteinRepresentations();
  registerProteinRepresentations({representation:{components:{polymer:{ref:id}},representations:{polymer:{ref:'cartoon-'+id,obj:{data:{sourceData:{root:proteinRoot}}}}}}});
  state.currentView='protein';state.aa='ALA';state.selectedProtein='P1';molstarViewer=viewer;
  molstarLatestRequest={protein:'P1',aa:'ALA'};
  return molstarLatestRequest;
}
loadFixture();bindProteinStyleEvents();
`);
(async () => {
  assert.equal(run('state.proteinStyle'), 'cartoon');
  assert.equal(buttons[0].attributes['aria-pressed'], 'true');
  assert.equal(built.length, 0, 'Surface is not built at startup');
  await run('setProteinPocketHighlight(viewer,{items:[{auth_asym_id:"A",auth_seq_id:42}]})');
  assert.equal(paints.get('cartoon-1')[0].color, 0x39ff14);
  buttons[1].handlers.click();
  await run('molstarUpdateQueue');
  assert.equal(built.length, 1);
  assert.equal(built[0].options.type, 'molecular-surface');
  assert.equal(built[0].options.color, 'plddt-confidence');
  assert.equal(built[0].options.typeParams.alpha, 0.45);
  assert.equal(built[0].params.initialState.isHidden, true, 'Surface starts hidden until ready');
  assert.equal(run('molstarProteinStyle'), 'surface');
  assert.equal(buttons[1].attributes['aria-pressed'], 'true');
  assert.deepEqual(visibility.slice(-2), [{ ref: 'cartoon-1', isHidden: true }, { ref: 'surface-1', isHidden: false }]);
  assert.equal(paints.get('surface-1')[0].color, 0x39ff14, 'Lazy surface receives green pocket color before display');
  assert.equal(paints.get('surface-1')[0].bundle.schema.items[0].auth_seq_id, 42);
  buttons[0].handlers.click(); await run('molstarUpdateQueue');
  assert.equal(run('molstarProteinStyle'), 'cartoon');
  buttons[1].handlers.click(); await run('molstarUpdateQueue');
  assert.equal(built.length, 1, 'Previously computed surface is reused');

  run('loadFixture("2")');
  assert.equal(run('molstarPocketColorSchema'), null, 'New protein cannot inherit a different pocket highlight');
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(built.length, 2, 'New model gets its own surface');
  assert.equal(run('molstarProteinStyle'), 'surface', 'Chosen style survives model changes');

  run('loadFixture("3")'); fail = true;
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(run('state.proteinStyle'), 'cartoon', 'Surface failure retains the working cartoon');
  assert.ok(node('#structure-style-status').textContent.includes('Previous view kept'));
  fail = false;
  buttons[1].handlers.click(); await run('molstarUpdateQueue');
  assert.equal(run('molstarProteinStyle'), 'surface', 'Retry works');

  run('loadFixture("4")');
  let release;
  gate = new Promise(resolve => { release = resolve; });
  const pending = run('applyProteinStyle(viewer,molstarLatestRequest)');
  await new Promise(resolve => setImmediate(resolve));
  const count = visibility.length;
  run('state.proteinStyle="cartoon"');
  release(); gate = null; await pending;
  assert.equal(visibility.length, count, 'Late surface result cannot override a newer Cartoon choice');
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(run('molstarProteinStyle'), 'cartoon');
  assert.equal(visibility.at(-1).isHidden, true);

  run('loadFixture("5");state.proteinStyle="surface"');
  await run('setProteinPocketHighlight(viewer,{items:[{auth_asym_id:"B",auth_seq_id:7}]})');
  failPaint = true;
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(run('state.proteinStyle'), 'cartoon', 'Failed surface paint leaves cartoon displayed');
  assert.equal(paints.has('surface-5'), false);
  const builtBeforeRetry = built.length;
  failPaint = false;
  buttons[1].handlers.click(); await run('molstarUpdateQueue');
  assert.equal(built.length, builtBeforeRetry, 'Paint retry reuses cached surface geometry');
  assert.equal(run('molstarProteinStyle'), 'surface');
  assert.equal(paints.get('surface-5')[0].bundle.schema.items[0].auth_seq_id, 7, 'Retry paints cached surface');
  await run('setProteinPocketHighlight(viewer,null)');
  assert.equal(paints.get('cartoon-5').length, 0);
  assert.equal(paints.get('surface-5').length, 0);
  assert.equal(warnings.length, 2);
  assert.equal(deleted.length, 0, 'No protein/pocket/ligand/cloud subtrees are removed by normal toggles');
  console.log('Protein style passed: lazy surface green coloring, independent layers, pLDDT/opacity, cached reuse, persistent preference, camera preservation, paint failure/retry and late results.');
})().catch(error => { console.error(error); process.exitCode = 1; });
