// Receptor-only coloring and independent pocket override. No network or data-file changes.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const nodes = new Map();
function node(selector) {
  if (!nodes.has(selector)) nodes.set(selector, { value: '', textContent: '', innerHTML: '', hidden: false,
    attributes: {}, setAttribute(name, value) { this.attributes[name] = value; },
    handlers: {}, addEventListener(type, handler) { this.handlers[type] = handler; } });
  return nodes.get(selector);
}
const proteinRoot = { atoms: ['ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'ALA', 'UNK',
  'ASN', 'CYS', 'GLN', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'MSE', 'SEP'].map((name, i) => ({
  name, chain: 'A', residue: i + 1,
})) };
const themes = new Map(), paints = new Map(), visibility = [], builds = [], warnings = [];
const paintTransform = {};
let failNext = false, gate = null, commits = 0;
const clone = value => JSON.parse(JSON.stringify(value));
function representation(ref, type = 'cartoon', params = {}) {
  themes.set(ref, { type: { name: type, params }, colorTheme: { name: 'plddt-confidence', params: {} },
    sizeTheme: { name: 'uniform', params: { value: 1 } } });
  return { ref, obj: { data: { sourceData: { root: proteinRoot } } } };
}
const viewer = { plugin: {
  build() {
    const changes = [], layers = [];
    return {
      to(repr) {
        assert.ok(/^(cartoon|surface)-/.test(repr.ref), 'Color changes cannot touch sticks, ligand or cloud');
        return {
          update(callback) {
            const draft = clone(themes.get(repr.ref));
            callback(draft);
            changes.push([repr.ref, draft]);
          },
          applyOrUpdateTagged(tag, transform, params) {
            assert.equal(tag, 'atlas-pocket-color'); assert.equal(transform, paintTransform);
            layers.push([repr.ref, params.layers]);
          },
        };
      },
      async commit() {
        if (gate) await gate;
        if (failNext) { failNext = false; throw Error('Intentional paint failure'); }
        for (const [ref, params] of changes) themes.set(ref, params);
        for (const [ref, params] of layers) paints.set(ref, params);
        commits++;
      },
      delete() { throw Error('Color changes must not remove structure layers'); },
    };
  },
  state: { data: { updateCellState(ref, update) { visibility.push({ ref, ...update }); } } },
  builders: { structure: { representation: { async addRepresentation(component, options) {
    builds.push(options);
    return representation('surface-' + component.ref, options.type, options.typeParams);
  } } } },
  managers: { camera: { focusLoci() { throw Error('Colors must not move camera'); }, focusSpheres() { throw Error('Colors must not move camera'); } } },
} };
const StructureElement = {
  Loci: {
    fromSchema(root, schema) {
      assert.equal(root, proteinRoot, 'Selections refer only to the receptor root');
      return { atoms: root.atoms.filter(atom => schema.items.some(item => item.label_comp_id
        ? item.label_comp_id === atom.name
        : item.auth_asym_id === atom.chain && item.auth_seq_id === atom.residue)) };
    },
    isEmpty: loci => !loci.atoms.length,
  },
  Bundle: { fromLoci: loci => ({ atoms: loci.atoms }) },
};
const context = vm.createContext({ viewer, proteinRoot, representation, console: { ...console, warn: (...args) => warnings.push(args) },
  document: { querySelector: node, querySelectorAll: () => [] },
  window: { molstar: { lib: { structure: { StructureElement },
    plugin: { StateTransforms: { Representation: { OverpaintStructureRepresentation3DFromBundle: paintTransform } } },
  } } },
  fetch() { throw Error('Color controls must not download anything'); },
});
for (const file of ['ligand-viewer.js', 'protein-viewer.js', 'app.js']) {
  vm.runInContext(fs.readFileSync(path.join(__dirname, '..', 'js', file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = code => vm.runInContext(code, context);
run(`function loadFixture(id) {
  resetProteinRepresentations();
  registerProteinRepresentations({representation:{components:{polymer:{ref:id}},representations:{polymer:representation('cartoon-'+id)}}});
  state.currentView='protein';state.selectedProtein=id;molstarViewer=viewer;
  molstarLatestRequest={protein:id,aa:state.aa};
}
loadFixture('P1');bindProteinColorEvents();bindPocketSticksEvents();`);
const color = value => node('#protein-color').handlers.change({ target: { value } });
const highlight = checked => node('#highlight-pocket-residues').handlers.change({ target: { checked } });
const settle = () => run('molstarUpdateQueue');
function shownColor(ref, name) {
  const atom = proteinRoot.atoms.find(atom => atom.name === name);
  let color = themes.get(ref).colorTheme.params.value;
  for (const layer of paints.get(ref) || []) if (layer.bundle.atoms.includes(atom)) color = layer.color;
  return color;
}
(async () => {
  assert.equal(run('state.showPocketSticks'), false, 'Sticks default to off');
  assert.equal(node('#show-pocket-sticks').checked, false);
  // Exercise independent coloring with sticks explicitly enabled by the user.
  node('#show-pocket-sticks').handlers.change({ target: { checked: true } }); await settle();
  assert.equal(node('#protein-color').value, 'plddt');
  assert.equal(node('#highlight-pocket-residues').checked, true);
  assert.equal(node('#structure-charge-note').hidden, true);
  assert.equal(node('#structure-polarity-note').hidden, true);
  await run('setProteinPocketHighlight(viewer,{items:[{auth_asym_id:"A",auth_seq_id:1}]})');
  color('charge'); await settle();
  assert.equal(node('#structure-charge-note').hidden, false);
  assert.ok(node('#structure-protein-legend').innerHTML.includes('HIS (variable)'));
  assert.equal(shownColor('cartoon-P1', 'ASP'), 0x39ff14, 'Pocket override wins over charge');
  assert.equal(shownColor('cartoon-P1', 'GLU'), 0xd1495b);
  assert.equal(shownColor('cartoon-P1', 'LYS'), 0x3675c5);
  assert.equal(shownColor('cartoon-P1', 'ARG'), 0x3675c5);
  assert.equal(shownColor('cartoon-P1', 'HIS'), 0xb88732);
  for (const name of ['ALA', 'UNK']) assert.equal(shownColor('cartoon-P1', name), 0xb8bdc5);
  assert.equal(paints.get('cartoon-P1').length, 4);

  highlight(false); await settle();
  assert.equal(shownColor('cartoon-P1', 'ASP'), 0xd1495b, 'Disabling green reveals underlying charge');
  assert.equal(paints.get('cartoon-P1').length, 3);
  assert.equal(run('state.showPocketSticks'), true, 'Override does not disable sticks');
  assert.ok(node('#structure-protein-legend').innerHTML.includes('Pocket sticks'));
  assert.equal(visibility.length, 0, 'Color/override events never change any layer visibility');

  // Newly built and cached surfaces follow the same colors even without a green pocket.
  await run('state.proteinStyle="surface";applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(shownColor('surface-P1', 'ASP'), 0xd1495b);
  assert.equal(themes.get('surface-P1').type.params.alpha, 0.45);
  assert.equal(builds.length, 1);
  const visibilityCount = visibility.length;
  highlight(true); await settle();
  for (const ref of ['cartoon-P1', 'surface-P1']) assert.equal(shownColor(ref, 'ASP'), 0x39ff14);
  assert.equal(visibility.length, visibilityCount);
  highlight(false); await settle();
  // Polarity covers all 20 canonical residues; unrecognized/modified names stay gray.
  color('polarity'); await settle();
  const polar = ['ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','LYS','SER','THR','TYR'];
  const nonpolar = ['ALA','ILE','LEU','MET','PHE','PRO','TRP','VAL'];
  assert.deepEqual(Array.from(run('RESIDUE_POLARITY_GROUPS.flatMap(group=>group.names)')).sort(),
    Array.from(run('AMINO_ACIDS.map(aa=>aa.code)')).sort(), 'Exactly one class per canonical AA');
  for (const ref of ['cartoon-P1', 'surface-P1']) {
    for (const name of polar) assert.equal(shownColor(ref,name),0x4d86b8,`${ref}/${name}: polar`);
    for (const name of nonpolar) assert.equal(shownColor(ref,name),0xc49a6c,`${ref}/${name}: nonpolar`);
    for (const name of ['UNK','MSE','SEP']) assert.equal(shownColor(ref,name),0xb8bdc5,`${ref}/${name}: unclassified`);
    assert.equal(paints.get(ref).length,2,'Polarity replaces charge layers');
  }
  assert.equal(themes.get('surface-P1').type.params.alpha,.45);
  assert.equal(node('#structure-charge-note').hidden,true);
  assert.equal(node('#structure-polarity-note').hidden,false);
  assert.equal(node('#protein-color').attributes['aria-describedby'],'structure-polarity-note');
  for (const label of ['Polarity','Polar','Nonpolar','Other/unknown']) assert.ok(node('#structure-protein-legend').innerHTML.includes(label));
  assert.ok(node('#structure-protein-legend').innerHTML.includes(polar.join(', ')),'Legend exposes the exact convention');
  highlight(true); await settle();
  for (const ref of ['cartoon-P1','surface-P1']) {
    assert.equal(shownColor(ref,'ASP'),0x39ff14,'Green override still wins');
    assert.equal(shownColor(ref,'GLU'),0x4d86b8);
  }
  highlight(false); await settle();
  assert.equal(shownColor('surface-P1','ASP'),0x4d86b8);
  assert.equal(visibility.length,visibilityCount,'Polarity does not toggle sticks, ligand or cloud');
  assert.equal(run('state.showPocketSticks && state.showLigand && state.showPocketCloud'),true);
  color('none'); await settle();
  for (const ref of ['cartoon-P1', 'surface-P1']) {
    assert.equal(paints.get(ref).length, 0, 'None clears charge and disabled pocket layers');
    for (const atom of proteinRoot.atoms) assert.equal(shownColor(ref, atom.name), 0xb8bdc5);
    assert.equal(themes.get(ref).sizeTheme.params.value, 1);
  }
  assert.equal(node('#structure-charge-note').hidden, true);
  assert.equal(node('#structure-polarity-note').hidden, true);
  assert.ok(node('#structure-protein-legend').innerHTML.includes('Uniform'));
  highlight(true); await settle();
  assert.equal(shownColor('surface-P1', 'ASP'), 0x39ff14, 'None supports an independent green override');
  color('plddt'); await settle();
  for (const ref of ['cartoon-P1', 'surface-P1']) {
    assert.equal(themes.get(ref).colorTheme.name, 'plddt-confidence');
    assert.equal(paints.get(ref).length, 1);
  }

  // Missing pockets retain charge coloring; changing a pocket while green is off is remembered.
  color('charge'); highlight(false); await settle();
  await run('setProteinPocketHighlight(viewer,null)');
  assert.equal(paints.get('surface-P1').length, 3);
  await run('setProteinPocketHighlight(viewer,{items:[{auth_asym_id:"A",auth_seq_id:3}]})');
  highlight(true); await settle();
  assert.equal(shownColor('surface-P1', 'LYS'), 0x39ff14);
  assert.equal(shownColor('surface-P1', 'ASP'), 0xd1495b, 'Old pocket does not remain green');
  highlight(false); await settle();
  node('#show-pocket-sticks').handlers.change({ target: { checked: false } }); await settle();
  assert.ok(!node('#structure-protein-legend').innerHTML.includes('confidence-pocket'), 'No green legend if neither green layer is enabled');

  // Session preferences survive a protein change but old residue selections do not.
  run('loadFixture("P2")');
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(run('state.proteinColor'), 'charge');
  assert.equal(run('state.showPocketHighlight'), false);
  assert.equal(run('molstarPocketColorSchema'), null);
  assert.equal(shownColor('surface-P2', 'ASP'), 0xd1495b);
  assert.equal(shownColor('surface-P2', 'LYS'), 0x3675c5);

  // Failed changes restore prior state and colors and can be retried.
  failNext = true; color('polarity'); await settle();
  assert.equal(run('state.proteinColor'), 'charge');
  assert.equal(shownColor('surface-P2', 'ASP'), 0xd1495b);
  assert.equal(node('#structure-color-status').hidden, false);
  color('polarity'); await settle();
  assert.equal(run('state.proteinColor'), 'polarity');
  assert.equal(shownColor('surface-P2','ASP'),0x4d86b8);
  assert.equal(node('#structure-color-status').hidden, true);
  assert.equal(warnings.length, 1);

  // Rapid changes and requests superseded before execution cannot leave stale colors.
  let release;
  gate = new Promise(resolve => { release = resolve; });
  color('charge'); await new Promise(resolve => setImmediate(resolve));
  color('plddt'); highlight(true);
  release(); gate = null; await settle();
  assert.equal(themes.get('surface-P2').colorTheme.name, 'plddt-confidence');
  assert.equal(run('molstarAppliedProteinColors.highlight'), true);
  color('polarity'); highlight(false); await settle();
  run('loadFixture("P3")');
  await run('applyProteinStyle(viewer,molstarLatestRequest)');
  assert.equal(run('state.proteinColor'),'polarity','Polarity preference persists between proteins');
  assert.equal(shownColor('surface-P3','SER'),0x4d86b8,'New surfaces inherit polarity');
  assert.equal(shownColor('surface-P3','VAL'),0xc49a6c);
  assert.equal(run('molstarPocketColorSchema'),null);
  const before = commits;
  color('polarity'); run('molstarLatestRequest=null'); await settle();
  assert.equal(commits, before, 'Superseded request does not mutate an old protein');
  assert.equal(builds.length, 3, 'One surface per protein; color controls do not rebuild surfaces');
  const html=fs.readFileSync(path.join(__dirname,'..','index.html'),'utf8');
  assert.equal((html.match(/<option value="polarity">Residue polarity<\/option>/g)||[]).length,1);
  console.log('Protein colors passed: charge/polarity classes, unknown residues, four modes, optional green override, legends, cached/new surfaces, isolated layers, failures and rapid changes.');
})().catch(error => { console.error(error); process.exitCode = 1; });
