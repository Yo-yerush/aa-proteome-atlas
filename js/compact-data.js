// Optional pipeline export: shared pocket metadata plus small per-AA score tables.
const COMPACT_METADATA_FIELDS = [
  "uniprot_id", "tair_id", "protein", "pocket", "rank", "score", "probability",
  "center_x", "center_y", "center_z", "sas_points", "surf_atoms", "n_pocket_residues",
  "n_plddt_matched", "mean_pocket_plddt", "min_pocket_plddt", "fraction_plddt_ge70",
  "fraction_plddt_ge90", "residue_ids",
];
const COMPACT_SCORE_FIELDS = [
  "vina_affinity", "vina_status", "sfct_vina_score", "sfct_score", "vina_sfct_combined",
  "sfct_best_pose", "sfct_n_poses", "sfct_status",
];

function validateCompactFile(entry) {
  if (!entry || !/^[a-z0-9_]+\.tsv\.gz$/.test(entry.file)
    || !Number.isSafeInteger(entry.rows) || entry.rows < 0 || !/^[a-f0-9]{64}$/.test(entry.sha256)) {
    throw new Error("Invalid compact data manifest: malformed file entry. Regenerate the compact export");
  }
}

async function readCompactTable(entry, fields, directory) {
  validateCompactFile(entry);
  const path = `${directory}/${entry.file}?v=${entry.sha256}`;
  const text = await loadGzippedResultText(path);
  if (globalThis.crypto?.subtle) {
    const digest = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(text));
    const hash = [...new Uint8Array(digest)].map((value) => value.toString(16).padStart(2, "0")).join("");
    if (hash !== entry.sha256) {
      throw new Error(`Compact data checksum mismatch: ${entry.file}. Publish the complete regenerated bundle in ${directory}/ and refresh`);
    }
  }
  const headers = text.split("\n", 1)[0].replace(/\r$/, "").split("\t");
  // E. coli adds gene_id after the shared v1 metadata columns; score columns remain exact.
  const expected = ["pocket_id", ...fields];
  if (fields === COMPACT_METADATA_FIELDS && headers.at(-1) === "gene_id") expected.push("gene_id");
  if (headers.length !== expected.length || headers.some((field, i) => field !== expected[i])) {
    throw new Error(`Invalid compact TSV headers: ${entry.file}. Regenerate it with prepare_atlas_data.py`);
  }
  const rows = parseTSV(text);
  if (rows.length !== entry.rows) throw new Error(`Compact row-count mismatch: ${entry.file}`);
  return rows;
}

async function loadCompactBundle(directory, allowedLigands) {
  const path = `${directory}/manifest.json`;
  const response = await atlasFetch(path, { cache: "no-cache" });
  // Only absence triggers legacy loading. A broken compact bundle is never silently ignored.
  if (response.status === 404) return null;
  if (!response.ok) throw new Error(`Could not load ${path} (HTTP ${response.status})`);
  const manifest = await response.json();
  if (manifest?.format !== "aa-proteome-atlas-compact" || manifest.version !== 1
    || !Array.isArray(manifest.ligands) || !manifest.ligands.length) {
    throw new Error("Unsupported compact data manifest. Regenerate it with prepare_atlas_data.py");
  }
  const supported = new Set(allowedLigands.map((aa) => aa.code));
  const ligands = new Map();
  for (const entry of manifest.ligands) {
    validateCompactFile(entry);
    if (!supported.has(entry.code) || ligands.has(entry.code)) {
      throw new Error(`Invalid or duplicate AA in compact manifest: ${entry.code}`);
    }
    ligands.set(entry.code, entry);
  }
  const pockets = new Map();
  for (const row of await readCompactTable(manifest.pockets, COMPACT_METADATA_FIELDS, directory)) {
    if (!/^[1-9]\d*$/.test(row.pocket_id) || pockets.has(row.pocket_id)
      || !row.uniprot_id || !row.protein || !row.pocket) {
      throw new Error("Invalid or duplicate pocket metadata ID in compact data");
    }
    pockets.set(row.pocket_id, row);
  }
  return { directory, pockets, ligands, pocketHash: manifest.pockets.sha256 };
}

async function loadCompactData() {
  const [lBundle, dBundle] = await Promise.all([
    loadCompactBundle(`${RESULTS_DIRECTORY}/L`, AMINO_ACIDS),
    loadCompactBundle(`${RESULTS_DIRECTORY}/D`, D_AA_CONTROLS),
  ]);
  let bundles = [lBundle, dBundle].filter(Boolean);
  if (dBundle && !lBundle) throw new Error(`D-AA controls were found but ${RESULTS_DIRECTORY}/L/manifest.json is missing`);
  // Preserve older single-bundle sites, but never mix their stale data into the split layout.
  if (!bundles.length) {
    const legacy = await loadCompactBundle(RESULTS_DIRECTORY, DATA_LIGANDS);
    if (!legacy) return null;
    bundles = [legacy];
  }
  const byAA = new Map();
  for (const bundle of bundles) {
    for (const code of bundle.ligands.keys()) byAA.set(code, bundle);
  }
  return { byAA, ligandCount: byAA.size };
}

async function loadCompactResultRows(aa, data) {
  const bundle = data.byAA.get(aa.code);
  if (!bundle) return []; // Missing controls never remove successful L-AA results.
  const entry = bundle.ligands.get(aa.code);
  if (!entry) return []; // Unprovided AAs stay missing without removing other AAs' proteins.
  const result = [], seen = new Set();
  for (const scores of await readCompactTable(entry, COMPACT_SCORE_FIELDS, bundle.directory)) {
    const pocket = bundle.pockets.get(scores.pocket_id);
    if (!pocket) throw new Error(`Unknown compact pocket ID for ${aa.code}: ${scores.pocket_id}`);
    const key = `${pocket.uniprot_id}|${pocket.pocket}`;
    if (seen.has(key)) throw new Error(`Duplicate compact protein/pocket for ${aa.code}: ${key}`);
    seen.add(key);
    const row = { ...pocket, ...scores, _compactDirectory: bundle.directory,
      _compactAA: aa.code, _compactPocketHash: bundle.pocketHash };
    if (resultStatus(row, "vina") !== "success") row.vina_affinity = NaN;
    if (resultStatus(row, "sfct") !== "success") {
      for (const field of ["sfct_vina_score", "sfct_score", "vina_sfct_combined"]) row[field] = NaN;
    }
    if (isSuccessfulResult(row)) result.push(addDerivedScores(row));
  }
  return result;
}
