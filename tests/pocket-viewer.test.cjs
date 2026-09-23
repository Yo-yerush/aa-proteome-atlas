// Pocket representation/lifecycle checks using the public Mol* API contract.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, { textContent: '', hidden: false, dataset: {}, handlers: {},
    attributes: {}, setAttribute(name, value) { this.attributes[name] = value; },
    addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
}
const receptor = { ref: 'receptor', obj: { data: { atoms: [
  { chain: 'A', residue: 1 }, { chain: 'A', residue: 2 }, { chain: 'B', residue: 1 },
] } } };
const schemas = [], components = [], representations = [], deleted = [], focused = [], actions = [];
const paints = new Map(), visibility = [], paintTransform = {};
const cells = new Set(['receptor', 'ligand']);
let serial = 0, gate = null, fail = false;
const StructureElement = {
  Loci: {
    fromSchema: (structure, schema) => {
      assert.equal(structure, receptor.obj.data, 'Select only on the recorded receptor');
      schemas.push(schema);
      return { atoms: structure.atoms.filter(atom => schema.items.some(item => item.auth_asym_id === atom.chain && item.auth_seq_id === atom.residue)) };
    },
    isEmpty: loci => !loci.atoms.length,
  },
  Bundle: { fromLoci: loci => ({ atoms: loci.atoms }) },
};
const viewer = {
  async structureInteractivity(options) { actions.push(options.action); },
  plugin: {
    build: () => ({
      delete: ref => ({ commit: async () => {
        assert.ok(!['receptor', 'ligand'].includes(ref)); deleted.push(ref); cells.delete(ref);
      } }),
      to: representation => ({ applyOrUpdateTagged(tag, transform, params) {
        assert.ok(['cartoon', 'surface'].includes(representation.ref), 'Paint only receptor representations');
        assert.equal(tag, 'atlas-pocket-color'); assert.equal(transform, paintTransform);
        paints.set(representation.ref, params.layers);
      }, update() {} }),
      async commit() {},
    }),
    state: { data: { updateCellState(ref, update) {
      assert.ok(ref.startsWith('sticks-'), 'Sticks visibility never touches other layers');
      visibility.push({ ref, ...update });
    } } },
    builders: { structure: {
      async tryCreateComponent(protein, params, key) {
        assert.equal(protein, receptor);
        assert.equal(params.type.name, 'bundle');
        assert.equal(key, 'atlas-selected-pocket');
        const ref = `pocket-${++serial}`; cells.add(ref);
        components.push({ ref, ...params });
        if (gate) await gate;
        return { ref };
      },
      representation: { async addRepresentation(component, options, params) {
        if (fail) throw Error('Representation failed');
        representations.push({ component, options, params });
        return { ref: 'sticks-' + component.ref };
      } },
    } },
    managers: { camera: { focusLoci: loci => focused.push(loci) } },
  },
};
const context = vm.createContext({ console,
  document: { querySelector: node, querySelectorAll: () => [] },
  window: { molstar: { lib: { structure: { StructureElement },
    plugin: { StateTransforms: { Representation: { OverpaintStructureRepresentation3DFromBundle: paintTransform } } },
  } } },
  viewer, receptor,
});
for (const file of ['organisms.js', 'ligand-viewer.js', 'protein-viewer.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
run(`molstarProteinStructure=receptor;
molstarViewer=viewer;state.currentView='protein';state.selectedProtein='P1';
molstarProteinRepresentations=Object.fromEntries(['cartoon','surface'].map(ref=>[ref,{ref,obj:{data:{sourceData:{root:receptor.obj.data}}}}]));
bindPocketSticksEvents();
function showPocket(ids) {
  const request={protein:'P1',aa:state.aa,row:{pocket:'pocket1',residue_ids:ids}};
  molstarLatestRequest=request;
  return focusMolstarPocket(viewer,request.row,request);
}`);
(async () => {
  assert.equal(run('state.showPocketSticks'), false);
  assert.equal(node('#show-pocket-sticks').checked, false);
  const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');
  assert.ok(!html.match(/<input\b[^>]*id="show-pocket-sticks"[^>]*>/)[0].includes('checked'));
  await run('showPocket("A_1 A_1 B_1")');
  assert.equal(schemas[0].items.length, 2, 'Duplicate residue IDs do not duplicate atoms');
  assert.deepEqual(components[0].type.params.atoms, [receptor.obj.data.atoms[0], receptor.obj.data.atoms[2]]);
  const options = representations[0].options;
  assert.equal(options.type, 'ball-and-stick');
  assert.deepEqual(Array.from(options.typeParams.visuals), ['intra-bond', 'inter-bond']);
  assert.equal(options.typeParams.ignoreHydrogens, true);
  assert.equal(options.color, 'uniform');
  assert.equal(options.colorParams.value, 0x39ff14);
  assert.equal(cells.size, 3);
  assert.equal(focused.length, 1);
  assert.ok(node('#structure-pocket-status').textContent.includes('residues mapped'));
  for (const layers of paints.values()) {
    assert.equal(layers.length, 1);
    assert.equal(layers[0].color, 0x39ff14);
    assert.deepEqual(layers[0].bundle.atoms, components[0].type.params.atoms);
  }
  assert.equal(representations[0].params.initialState.isHidden, true);
  assert.equal(visibility.at(-1).isHidden, true, 'Sticks default to hidden');
  const paintBeforeToggle = paints.get('cartoon');
  node('#show-pocket-sticks').handlers.change({ target: { checked: false } });
  await run('molstarUpdateQueue');
  assert.equal(visibility.at(-1).isHidden, true);
  assert.equal(paints.get('cartoon'), paintBeforeToggle, 'Hiding sticks leaves green protein untouched');
  await run('showPocket("A_2")');
  assert.equal(cells.size, 3, 'Pocket changes replace rather than accumulate sticks');
  assert.deepEqual(components[1].type.params.atoms, [receptor.obj.data.atoms[1]]);
  assert.equal(visibility.at(-1).isHidden, true, 'New pocket respects hidden-sticks preference');
  for (const layers of paints.values()) {
    assert.equal(layers.length, 1);
    assert.deepEqual(layers[0].bundle.atoms, [receptor.obj.data.atoms[1]], 'Old green region is replaced');
  }
  node('#show-pocket-sticks').handlers.change({ target: { checked: true } });
  await run('molstarUpdateQueue');
  assert.equal(visibility.at(-1).isHidden, false);
  await run('showPocket("")');
  assert.equal(cells.size, 2);
  assert.ok([...paints.values()].every(layers => !layers.length), 'Missing pocket restores pLDDT');
  assert.ok(node('#structure-pocket-status').textContent.includes('no residue list'));
  await run('showPocket("Z_999")');
  assert.equal(cells.size, 2, 'Unknown residues never select the ligand');
  assert.ok(node('#structure-pocket-status').textContent.includes('not found'));
  fail = true;
  await assert.rejects(run('showPocket("A_1")'), /Representation failed/);
  assert.equal(cells.size, 2, 'Failed pocket layers are cleaned up');
  assert.ok([...paints.values()].every(layers => layers.length === 1), 'Protein stays green even if sticks fail');
  fail = false;
  let release;
  gate = new Promise(resolve => { release = resolve; });
  const pending = run('showPocket("A_1")');
  await new Promise(resolve => setImmediate(resolve));
  run('molstarLatestRequest=null');
  const before = representations.length;
  release(); gate = null;
  await pending;
  assert.equal(cells.size, 2, 'Stale asynchronous pocket components are removed');
  assert.equal(representations.length, before);
  assert.ok(actions.every(action => action === 'select'), 'Selection markers are only cleared, not left on the cartoon');
  assert.ok(deleted.length > 0 && cells.has('receptor') && cells.has('ligand'));
  console.log('Pocket viewer passed: receptor-only green coloring, independent Show sticks, layer replacement, missing residues, failures, cancellation and ligand preservation.');
})().catch(error => { console.error(error); process.exitCode = 1; });
