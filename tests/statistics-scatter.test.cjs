// Dependency-free checks for full-population Statistics scatter rendering.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const root = path.resolve(__dirname, '..');
const nodes = new Map();
const node = selector => {
  if (!nodes.has(selector)) nodes.set(selector, { innerHTML: '', textContent: '', hidden: false });
  return nodes.get(selector);
};
const drawn = [];
const drawing = {
  setTransform() {}, beginPath() {}, fill() {},
  arc(x, y, radius) { assert.ok(Number.isFinite(x) && Number.isFinite(y)); drawn.push({ x, y, radius }); },
};
const canvas = { handlers: {}, getContext: () => drawing,
  addEventListener(type, callback) { this.handlers[type] = callback; },
  getBoundingClientRect: () => ({ left: 0, top: 0, width: 600, height: 360 }) };
const tooltip = { hidden: true, textContent: '' };
node('#statistics-scatter').querySelector = selector => selector === 'canvas' ? canvas : tooltip;
const context = vm.createContext({ console, URLSearchParams, devicePixelRatio: 2,
  document: { querySelector: node, querySelectorAll: () => [] } });
for (const file of ['js/app.js', 'js/analysis.js']) {
  vm.runInContext(fs.readFileSync(path.join(root, file), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
}
const run = source => vm.runInContext(source, context);
run(`
  const scatterRow = (protein, pocket, score) => addDerivedScores({ uniprot_id: protein, pocket,
    vina_affinity: score, sfct_score: 2 * score + .5, sfct_vina_score: score,
    vina_sfct_combined: 1.8 * score + .4, vina_status: 'success', sfct_status: 'success',
    probability: .9, mean_pocket_plddt: 95 });
  const scatterRows = Array.from({ length: 2603 }, (_, i) =>
    [scatterRow('P' + i, '1', -10 + i / 1000), scatterRow('P' + i, '2', -9 + i / 1000)]).flat();
  // A late extreme must affect the full-population plot, even beyond the old cap.
  scatterRows[scatterRows.length - 1].sfct_score = 12345;
  state.rawByAA.set('ALA', [...scatterRows,
    { ...scatterRow('FAILED', '1', -20), vina_status: 'failed', sfct_status: 'failed' },
    { ...scatterRow('LOW_QC', '1', -20), probability: .1 },
    { ...scatterRow('NO_SFCT', '1', -20), sfct_score: NaN }]);
  state.rawByAA.set('LEU', scatterRows.map(row => ({ ...row, vina_affinity: row.vina_affinity + 1 })));
  state.rawByAA.set('DALA', scatterRows.slice(0, -6).map(row => ({ ...row, vina_affinity: row.vina_affinity + .5 })));
  analysisState.statistics.metric = 'vina_affinity';
`);

for (const unit of ['best', 'all']) {
  for (const aa of ['ALA', 'ALL']) {
    context.testUnit = unit;
    context.testAA = aa;
    run(`Object.assign(analysisState.statistics, { aa: testAA, unit: testUnit });
      analysisResults.statistics = macroStatistics(analysisState.statistics);`);
    const ids = Array.from(run('analysisResults.statistics.correlations.map(row => row.id)'));
    for (const id of ids) {
      context.testScatter = id;
      drawn.length = 0;
      run('analysisState.statistics.scatter = testScatter; renderStatisticsScatter();');
      const count = run('analysisResults.statistics.correlations.find(row => row.id === testScatter).n');
      assert.ok(count > 2000, `${aa}/${unit}/${id} exercises the old cap`);
      assert.equal(drawn.length, count, `${aa}/${unit}/${id}: every eligible pair is drawn`);
      assert.ok(node('#statistics-scatter').innerHTML.includes(`${count.toLocaleString()} of ${count.toLocaleString()} points shown`));
      assert.ok(!node('#statistics-scatter').innerHTML.includes('<circle'), 'No per-point SVG elements');
      assert.ok(node('#statistics-scatter-note').textContent.includes('All eligible paired observations'));
      if (id === 'vina_sfct') {
        assert.equal(count, 2603 * (unit === 'all' ? 2 : 1) * (aa === 'ALL' ? 2 : 1));
        assert.ok(!run(`statisticsCorrelationPairs(analysisResults.statistics,
          analysisResults.statistics.correlations.find(row => row.id === testScatter))
          .some(pair => ['FAILED', 'LOW_QC', 'NO_SFCT'].includes(pair.protein))`));
        if (unit === 'all') {
          assert.equal(run(`Math.max(...statisticsCorrelationPairs(analysisResults.statistics,
            analysisResults.statistics.correlations.find(row => row.id === testScatter)).map(pair => pair.y))`), 12345);
        }
      }
    }
  }
}
assert.equal(canvas.width, 1200);
assert.equal(canvas.height, 720);
const point = drawn[drawn.length - 1];
canvas.handlers.pointermove({ clientX: point.x, clientY: point.y });
assert.equal(tooltip.hidden, false);
assert.ok(tooltip.textContent.includes('P2602'));
canvas.handlers.pointerleave();
assert.equal(tooltip.hidden, true);

// Empty and singleton populations remain safe, with undefined correlations retained.
run(`analysisResults.statistics = { records: [], stereo: [], correlations: [
  { id: 'vina_sfct', xKey: 'vina_affinity', yKey: 'sfct_score', xLabel: 'Vina', yLabel: 'SFCT', n: 0, r: NaN, rho: NaN }
] }; analysisState.statistics.scatter = 'vina_sfct'; renderStatisticsScatter();`);
assert.ok(node('#statistics-scatter').innerHTML.includes('No paired observations'));
drawn.length = 0;
run(`analysisResults.statistics.records = [{ aa: 'ALA', row: scatterRow('SINGLE', '1', -8) }];
  analysisResults.statistics.correlations[0].n = 1; renderStatisticsScatter();`);
assert.equal(drawn.length, 1);
assert.ok(Number.isNaN(run('pairStatistics([{ x: 1, y: 1 }]).r')));
// Compare AAs keeps its existing display behavior; this change targets Statistics only.
assert.equal((run('analysisScatter(scatterRows.map(row => ({ x: row.vina_affinity, y: row.sfct_score })), "X", "Y")')
  .match(/<circle /g) || []).length, 2000);
console.log('Statistics scatter checks passed: all points beyond 2,000, all variable pairs, pooled/best/all-pocket modes, L/D, QC/missing values, hover, empty/singleton data.');
