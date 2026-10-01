// Focused human integration: one real AA table, supplied annotations and optional provenance.
const assert = require('node:assert/strict');
const { fixture } = require('./organisms.test.cjs');

(async () => {
  const hs = fixture('?organism=human&aa=ALA');
  hs.run('initializeOrganismUI();readURLState()');
  assert.equal(hs.run('ORGANISM.scientificName'), 'Homo sapiens');
  assert.equal(fixture().run('ORGANISM.id'), 'ecoli');
  await hs.run(`(async()=>{
    globalThis.bundle=await loadCompactData();
    state.rawByAA.set('ALA',await loadCompactResultRows(AMINO_ACIDS[0],bundle));
    for(const row of state.rawByAA.get('ALA')) state.metadata.set(row.uniprot_id,row);
    await loadUniProtAnnotations();
  })()`);
  assert.equal(hs.run('bundle.byAA.size'), 20);
  assert.equal(hs.run('hasDControl()'), false);
  assert.equal(hs.run('getStereoControl("A2VEC9")'), null);
  assert.ok(hs.run('state.rawByAA.get("ALA").length > 30000'));
  assert.ok(hs.run('state.rawByAA.get("ALA").every(row=>row._compactDirectory==="Hs_results/L")'));
  assert.ok(hs.run('findProfileProteins("NQO1").some(row=>row.protein==="P15559")'));
  assert.match(hs.run('proteinTableIdentity("P15559")'), /NQO1/);
  assert.doesNotMatch(hs.run('proteinTableIdentity("P15559")'), /NQO1 · NQO1|TAIR|Araport/);

  await hs.run('loadGeneDescriptions().then(data=>globalThis.descriptions=data)');
  assert.deepEqual(Array.from(hs.run('matchingGeneDescriptions("P15559",[],descriptions).keys()')), ['ENSG00000181019']);
  assert.deepEqual(Array.from(hs.run('matchingGeneDescriptions("ABSENT",["NQO1"],descriptions).keys()')), ['ENSG00000181019']);
  assert.deepEqual(Array.from(hs.run('matchingGeneDescriptions("A2VEC9",[],descriptions).keys()')), ['ENSG00000197558']);
  const dialog = hs.nodes.get('#gene-description-dialog') || hs.context.document.querySelector('#gene-description-dialog');
  dialog.showModal = function () { this.open = true; };
  await hs.run('openGeneDescriptions("P15559",[],null)');
  assert.match(hs.nodes.get('#gene-description-body').innerHTML, /NAD\(P\)H|quinone/);
  assert.match(hs.nodes.get('#gene-description-ids').textContent, /ENSG00000181019/);
  assert.match(hs.nodes.get('#gene-description-body').innerHTML, /GO biological process/);

  // Identical symbols on two genes must not pick an arbitrary row; exact accessions win.
  hs.context.ambiguousCSV = 'gene_id,Symbol,UniProt,Protein_name\nENSG00000000001,DUP,P1,First\nENSG00000000002,DUP,P2,Second\n';
  hs.run('globalThis.ambiguous=parseGeneDescriptionCSV(ambiguousCSV)');
  assert.equal(hs.run('matchingGeneDescriptions("MISSING",["DUP"],ambiguous).size'), 0);
  assert.deepEqual(Array.from(hs.run('matchingGeneDescriptions("P2",["DUP"],ambiguous).keys()')), ['ENSG00000000002']);
  assert.equal(hs.run('matchingGeneDescriptions("P2-2",[],ambiguous).size'), 0, 'Do not strip isoform suffixes');
  assert.deepEqual(Array.from(hs.run('organismGeneIds("ENSG00000181019.2 NQO1")')), ['ENSG00000181019','NQO1']);

  hs.run(`globalThis.fragmentMetadata=[...bundle.byAA.get('ALA').pockets.values()].filter(row=>row.uniprot_id==='A2VEC9'&&row.pocket==='pocket3');
    globalThis.fragment=state.rawByAA.get('ALA').find(row=>row.protein==='AF-A2VEC9-F20-model_v6'&&row.pocket==='pocket3');
    state.metric='vina_affinity';state.p2rank=0;state.plddt=0;switchView=()=>{};
    selectProtein('A2VEC9',fragment);`);
  assert.equal(hs.run('fragmentMetadata.length'), 5);
  assert.equal(hs.run('new Set(fragmentMetadata.map(pocketKey)).size'), 5);
  assert.equal(hs.run('state.profilePocketAnchor.protein'), 'AF-A2VEC9-F20-model_v6');
  assert.equal(hs.run('getRanking().filter(row=>row.uniprot_id==="A2VEC9").length'), 1);
  const tsv = hs.run('rowsToDelimited([getRankedProtein("A2VEC9","ALA",state.metric)],"\t")');
  assert.match(tsv, /gene_id/); assert.match(tsv, /AF-A2VEC9-F/);
  hs.run('downloadText("scores.tsv","data");updateURL()');
  assert.equal(hs.downloads.at(-1), 'human_scores.tsv');
  assert.match(hs.navigations.at(-1), /organism=human/);

  await hs.run(`(async()=>{
    globalThis.goData=await loadGOAnnotations();
    globalThis.controls=await loadControlQCData();
    globalThis.row=state.rawByAA.get('ALA')[0];globalThis.request={row,aa:'ALA'};
    globalThis.poses=await loadLigandPositions(ligandPositionSource(request));
    globalThis.points=await loadPocketPointTable(pocketPointSource(request));
    globalThis.potentialMetadata=await loadPocketElectrostaticsMetadata();
    globalThis.poseBundle=await loadPoseElectrostaticsBundle();
  })()`);
  assert.ok(hs.run('goData.proteins.size > 0 && controls.controls.length > 0'));
  assert.ok(hs.run('poses.size > 0 && points.size > 0'));
  assert.equal(hs.run('potentialMetadata.pocketHash'), hs.run('row._compactPocketHash'));
  assert.equal(hs.run('poseBundle.pocketHash'), hs.run('row._compactPocketHash'));
  assert.equal(hs.run('potentialMetadata.summary.get("AF-Q6UB99-F1-model_v6|1").status'), 'points_error');
  hs.context.failedSummary = 'protein\tpocket_rank\tn_points\tpoints_sha256\tphi_mean\tphi_min\tphi_max\tfraction_positive\tfraction_negative\tstatus\terror\nAF-TEST-F1-model_v6\t1\t0\t\t\t\t\t\t\tpoints_error\tNo coordinates\n';
  assert.equal(hs.run('parsePocketPotentialSummary(failedSummary).size'), 1);
  assert.throws(()=>hs.run('parsePocketPotentialSummary(failedSummary.replace("points_error","success"))'), /Invalid/);
  assert.throws(()=>hs.run('parsePocketPotentialSummary(failedSummary.replace("points_error","fragment_only"))'), /Invalid/);
  assert.equal(hs.run('ligandPositionSource({aa:"ALA",row:{...row,_compactDirectory:"Ec_results/L"}})'), null);
  assert.ok(hs.requests.every(url=>url.startsWith('Hs_results/')||url.startsWith('annotations/human/')));

  // Check the two existing lookup modes without loading their full score datasets.
  for (const [organism,id,csv] of [
    ['ecoli','b0001','gene_id,product\nb0001,Leader peptide\n'],
    ['arabidopsis','AT1G01010','gene_id,Short_description\nAT1G01010,Transcription factor\n'],
  ]) {
    const other=fixture('?organism='+organism);
    other.context.testCSV=csv; other.context.testId=id;
    assert.equal(other.run('matchingGeneDescriptions("P1",[testId],parseGeneDescriptionCSV(testCSV)).get(testId).length'),1);
    other.run('initializeOrganismUI();switchOrganism("human")');
    assert.match(other.navigations.at(-1), /organism=human/);
  }
  hs.run('switchOrganism("ecoli");switchOrganism("human");switchOrganism("arabidopsis")');
  assert.equal(hs.run('organismAbortController.signal.aborted'), true);
  assert.match(hs.navigations.at(-1), /organism=arabidopsis/);
  const absent=fixture('?organism=human',()=>true);
  await assert.rejects(absent.run('loadGeneDescriptions()'), /404/);
  await assert.rejects(absent.run('loadElectrostaticsManifest()'), /404/);
  assert.ok(absent.requests.every(url=>url.startsWith('Hs_results/')||url.startsWith('annotations/human/')));
  console.log('Human integration passed: real L-only data, annotation joins, fragments, GO/QC, optional provenance, downloads and isolated organism navigation.');
})().catch(error=>{console.error(error);process.exitCode=1;});
