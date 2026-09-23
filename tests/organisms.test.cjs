// Real local exports, separate document contexts and forced optional-file failures.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const { webcrypto } = require('node:crypto');
const root = path.resolve(__dirname, '..');
const scripts = [...fs.readFileSync(path.join(root, 'index.html'), 'utf8').matchAll(/<script src="(js\/[^?]+)\?/g)].map(m => m[1]);

function fixture(search = '', missing = () => false) {
  const nodes = new Map(), requests = [], downloads = [], navigations = [];
  const node = selector => {
    if (!nodes.has(selector)) nodes.set(selector, { value: '', textContent: '', innerHTML: '', hidden: false,
      dataset: {}, classList: { add() {}, remove() {}, toggle() {} }, setAttribute() {},
      addEventListener() {}, remove() {}, click() { downloads.push(this.download); } });
    return nodes.get(selector);
  };
  const location = { search, protocol: 'http:', pathname: '/atlas/', hash: '#explorer',
    href: `http://localhost/atlas/${search}#explorer`, assign: url => navigations.push(url), reload() {} };
  const context = vm.createContext({ console, URL, URLSearchParams, Blob, Response, TextEncoder, TextDecoder,
    DecompressionStream, AbortController, crypto: webcrypto, setTimeout, clearTimeout, location,
    document: { querySelector: node, querySelectorAll: () => [], createElement: () => node('anchor'), body: { appendChild() {} } },
    window: { addEventListener() {} }, history: { replaceState: (_, __, url) => navigations.push(url) },
    fetch: async (url, options) => {
      requests.push(String(url));
      options?.signal?.throwIfAborted();
      const file = path.join(root, String(url).split('?')[0]);
      return missing(String(url)) || !fs.existsSync(file) ? new Response('', { status: 404 }) : new Response(fs.readFileSync(file));
    },
  });
  for (const script of scripts) vm.runInContext(fs.readFileSync(path.join(root, script), 'utf8').replace(/\binit\(\);\s*$/, ''), context);
  const run = source => vm.runInContext(source, context);
  run('showToast=()=>{}');
  return { run, requests, downloads, navigations, nodes, context };
}

(async () => {
  for (const search of ['', '?organism=unknown', '?organism=__proto__', '?organism=ecoli']) {
    assert.equal(fixture(search).run('ORGANISM.id'), 'ecoli');
  }
  const ec = fixture('?organism=ecoli&aa=MET&q=AT1G01010');
  ec.run('initializeOrganismUI();readURLState()');
  assert.equal(ec.run('ORGANISM.id'), 'ecoli');
  await ec.run('loadData()');
  assert.equal(ec.run('state.rawByAA.size'), 39);
  assert.equal(ec.run('state.metadata.size'), 3085); // Three models have no successful scores.
  assert.equal(ec.run('state.rawByAA.get("ALA").length'), 6780);
  assert.ok(ec.run('state.rawByAA.get("ALA").every(row=>row.gene_id && row._compactDirectory === "Ec_results/L")'));
  assert.ok(ec.run('state.rawByAA.get("DALA").every(row=>row._compactDirectory === "Ec_results/D")'));
  ec.run('state.search="";state.aa="ALA";state.filtered=filterRows();updateURL()');
  assert.doesNotMatch(ec.navigations.at(-1), /organism=/);
  assert.ok(ec.run('state.filtered.length > 100'));
  assert.ok(ec.run('getRanking().every(row=>state.metadata.has(row.uniprot_id))'));
  assert.ok(ec.run('findProfileProteins("b0002").length > 0'));
  assert.ok(ec.run('findProfileProteins("thrA").length > 0'));
  assert.ok(ec.run('annotationSearchText("P00561").includes("thrA")'));
  assert.match(ec.run('rowsToDelimited(state.filtered.slice(0,2),"\t")'), /gene_id/);
  assert.doesNotMatch(ec.run('rowsToDelimited(state.filtered.slice(0,2),"\t")'), /tair_id/);
  await ec.run('loadGeneDescriptions().then(records=>globalThis.descriptions=records)');
  assert.match(ec.run('geneDescriptionSections(["b0001"],descriptions)'), /thr operon leader peptide/);
  assert.deepEqual(Array.from(ec.run('organismGeneIds("b0001 B0002 AT1G01010")')), ['b0001','b0002']);
  assert.match(ec.run('proteinTableIdentity("P00561")'), /uniprot.org/);
  assert.doesNotMatch(ec.run('proteinTableIdentity("P00561")'), /TAIR|Araport/);
  await ec.run('loadGOAnnotations().then(data=>{ globalThis.goData=data;globalThis.enrichment=computeGOEnrichment(data,goState); })');
  assert.ok(ec.run('enrichment.background.every(row=>state.metadata.has(row.uniprot_id))'));
  await ec.run('loadControlQCData().then(data=>{globalThis.controls=data;globalThis.qc=calculateControlQC(data,controlQCState)})');
  assert.ok(ec.run('controls.controls.every(row=>row.source.organism.includes("Escherichia coli"))'));
  await ec.run('prepareExplorerPoseElectrostatics()');
  assert.equal(ec.run('poseElectrostaticsTables.values().next().value.kind'), 'ready');
  await ec.run(`(async()=>{
    globalThis.row=state.rawByAA.get('ALA')[0];globalThis.request={row,aa:'ALA'};
    globalThis.positions=await loadLigandPositions(ligandPositionSource(request));
    globalThis.pointsTable=await loadPocketPointTable(pocketPointSource(request));
    globalThis.points=parsePocketPoints(pointsTable.get(row.pocket_id));
    globalThis.potentials=await loadValidatedPocketPotentials(request,points);
  })()`);
  assert.ok(ec.run('positions.has("1:1") && points.length > 0 && potentials.values.length === points.length'));
  assert.equal(ec.run('ligandPositionSource({aa:"ALA",row:{...row,_compactDirectory:"At_results/L"}})'), null);
  assert.equal(ec.run('pocketPointSource({aa:"ALA",row:{...row,_compactDirectory:"At_results/L"}})'), null);
  assert.equal(ec.run('poseElectrostaticsSource({...row,_compactDirectory:"At_results/L"},"ALA")'), null);
  ec.run('downloadText("scores.tsv","data")');
  assert.equal(ec.downloads.at(-1), 'ecoli_scores.tsv');
  assert.ok(ec.requests.every(url => url.startsWith('Ec_results/') || url.startsWith('annotations/ecoli/')));

  const at = fixture('?organism=arabidopsis');
  assert.equal(at.run('ORGANISM.id'), 'arabidopsis');
  at.run('updateURL()');
  assert.match(at.navigations.at(-1), /organism=arabidopsis/);
  await at.run(`(async()=>{const bundle=await loadCompactData();state.rawByAA.set('ALA',await loadCompactResultRows(AMINO_ACIDS[0],bundle));await loadUniProtAnnotations();})()`);
  assert.equal(at.run('state.rawByAA.get("ALA").length'), 33622);
  assert.equal(at.run('state.rawByAA.get("ALA")[0].pocket_id'), ec.run('row.pocket_id'));
  assert.notEqual(at.run('state.rawByAA.get("ALA")[0].uniprot_id'), ec.run('row.uniprot_id'));
  assert.deepEqual(Array.from(at.run('organismGeneIds("At1g01010.1 AT1G01010.2 b0001")')), ['AT1G01010']);
  await at.run('loadGeneDescriptions().then(records=>globalThis.descriptions=records)');
  assert.match(at.run('geneDescriptionSections(["AT1G01010"],descriptions)'), /NAC domain/);
  assert.ok(at.requests.every(url => url.startsWith('At_results/') || url.startsWith('annotations/arabidopsis/')));

  // Switching aborts old requests and blocks late exports. Rapid switches, including
  // back to the original organism before navigation commits, always navigate again.
  ec.run('switchOrganism("arabidopsis");switchOrganism("ecoli");switchOrganism("arabidopsis")');
  assert.match(ec.navigations.at(-1), /organism=arabidopsis/);
  assert.doesNotMatch(ec.navigations.at(-1), /[?&]q=/);
  assert.match(ec.navigations.at(-1), /aa=MET/);
  assert.equal(ec.run('organismAbortController.signal.aborted'), true);
  await assert.rejects(ec.run('atlasFetch("Ec_results/L/manifest.json")'), /abort/i);
  ec.run('downloadText("stale.tsv","stale")');
  assert.equal(ec.downloads.length, 1);
  assert.equal(at.run('state.annotations.has("P00561")'), false);
  const freshEc = fixture('?organism=ecoli');
  assert.ok(freshEc.run('state.rawByAA.size===0 && state.annotations.size===0 && state.rankingCache.size===0 && state.selectedProtein===null && state.selectedPocket===null && goResult===null && controlQCResult===null && ligandPositionCache.size===0 && pocketPotentialFiles.size===0 && poseElectrostaticsTables.size===0 && pocketPointTableCache===null'));

  const absent = fixture('?organism=ecoli', url => !url.startsWith('Ec_results/L/'));
  assert.equal(await absent.run('loadCompactData().then(data=>data.byAA.size)'), 20);
  await absent.run('loadUniProtAnnotations()');
  assert.equal(absent.run('state.annotations.size'), 0);
  await assert.rejects(absent.run('loadGeneDescriptions()'), /404/);
  await assert.rejects(absent.run('loadGOAnnotations()'), /404/);
  await assert.rejects(absent.run('loadControlQCData()'), /404/);
  await assert.rejects(absent.run('loadElectrostaticsManifest()'), /404/);
  await assert.rejects(absent.run('loadPocketPointTable("Ec_results/pocket_points.tsv.gz")'), /404/);
  absent.context.fetch = async url => { absent.requests.push(url); return new Response('', {status:404}); };
  await assert.rejects(absent.run('loadLigandPositions("Ec_results/L/aa_positions/positions_ala.tsv.gz")'), /404/);
  assert.ok(absent.requests.every(url => !url.includes('At_results') && !url.includes('arabidopsis')));
  console.log('Organisms passed: real L/D exports, annotations, GO/QC, poses, points/potentials, duplicate pocket IDs, missing files, downloads and rapid navigation isolation.');
})().catch(error => { console.error(error); process.exitCode = 1; });
