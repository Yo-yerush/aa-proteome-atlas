// Explorer-only annotation for the exported Vina MODEL 1, not a new ranking metric.
const POSE_ELECTROSTATICS_COLUMNS = ["pocket_id", "vina_pose", "n_atoms", "phi_mean", "phi_min", "phi_max", "qphi_kT", "status", "error"];
const POSE_QPHI_DESCRIPTION = "Vina pose 1 · qφ (kT). Negative is favorable; positive is unfavorable in the fixed receptor field. Not binding free energy or a desolvation-corrected score. May differ from the SFCT/Combined-selected pose.";
let poseElectrostaticsBundlePromise = null;
const poseElectrostaticsTables = new Map();
let explorerPoseObservation = null;

function poseElectrostaticsSource(row, aa) {
  // Explorer shows canonical L AAs. Never substitute L potentials for D/legacy results.
  if (!row || row._compactDirectory !== `${RESULTS_DIRECTORY}/L` || row._compactAA !== aa
    || !AMINO_ACIDS.some((entry) => entry.code === aa) || !isPotentialHash(row._compactPocketHash)
    || !/^[1-9]\d*$/.test(String(row.pocket_id))) return null;
  return { aa, directory: row._compactDirectory, pocketHash: row._compactPocketHash,
    key: `${row._compactDirectory}|${aa}|${row._compactPocketHash}` };
}

async function loadPoseElectrostaticsBundle() {
  if (poseElectrostaticsBundlePromise) return poseElectrostaticsBundlePromise;
  const pending = (async () => {
    const manifest = await loadElectrostaticsManifest();
    if (JSON.stringify(manifest.poses?.columns) !== JSON.stringify(POSE_ELECTROSTATICS_COLUMNS)
      || manifest.poses.pose !== "Vina MODEL 1" || manifest.poses.qphi_units !== "kT"
      || !Array.isArray(manifest.bundles)) throw new Error("Unsupported pose electrostatics format, units or pose identity");
    const text = await readPocketPotentialFile(`${RESULTS_DIRECTORY}/L/manifest.json`, null, false);
    const digest = await pocketPotentialSHA256(new TextEncoder().encode(text));
    const matches = manifest.bundles.filter((bundle) => bundle.manifest_sha256 === digest);
    const compact = JSON.parse(text);
    if (matches.length !== 1 || compact.format !== "aa-proteome-atlas-compact" || compact.version !== 1
      || !isPotentialHash(compact.pockets?.sha256)) throw new Error("Pose electrostatics does not match the current compact L bundle");
    return { bundle: matches[0], pocketHash: compact.pockets.sha256 };
  })();
  poseElectrostaticsBundlePromise = pending;
  try { return await pending; }
  catch (error) {
    if (poseElectrostaticsBundlePromise === pending) poseElectrostaticsBundlePromise = null;
    throw error;
  }
}

function parsePoseElectrostatics(text, expectedRows) {
  const table = new Map();
  for (const [id, pose, atoms, , , , rawValue, status, error] of pocketPotentialRows(text, POSE_ELECTROSTATICS_COLUMNS)) {
    if (!/^[1-9]\d*$/.test(id) || table.has(id) || !status || (pose !== "" && pose !== "1")) {
      throw new Error("Invalid/duplicate pocket ID or unexpected Vina pose in electrostatics");
    }
    if (status === "success" && (pose !== "1" || !/^[1-9]\d*$/.test(atoms))) {
      throw new Error("Successful pose electrostatics lacks Vina pose 1 or its atom count");
    }
    let value = NaN;
    if (status === "success" && rawValue.trim() !== "") {
      if (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(rawValue) || !Number.isFinite(Number(rawValue))) {
        throw new Error("Invalid qphi_kT value in pose electrostatics");
      }
      value = Number(rawValue);
    }
    // Unsuccessful/clashing rows never acquire a value, even if an exporter supplied one.
    const missingValue = status === "success" && !Number.isFinite(value);
    table.set(id, { value, pose: pose === "1" ? 1 : null, status: missingValue ? "missing_value" : status,
      error: missingValue ? "Successful calculation has no qphi_kT value" : error });
  }
  if (table.size !== expectedRows) throw new Error("Pose electrostatics row count does not match its manifest");
  return table;
}

function ensurePoseElectrostaticsTable(source, aa) {
  const key = source?.key || `unavailable|${aa}`;
  if (poseElectrostaticsTables.has(key)) {
    const entry = poseElectrostaticsTables.get(key);
    poseElectrostaticsTables.delete(key);
    poseElectrostaticsTables.set(key, entry);
    return entry;
  }
  const entry = { key, kind: source ? "loading" : "unavailable", table: null,
    error: source ? "" : "No matching compact L pocket provenance for this AA" };
  poseElectrostaticsTables.set(key, entry);
  while (poseElectrostaticsTables.size > 2) poseElectrostaticsTables.delete(poseElectrostaticsTables.keys().next().value);
  entry.promise = (async () => {
    if (!source) return;
    const { bundle, pocketHash } = await loadPoseElectrostaticsBundle();
    if (source.pocketHash !== pocketHash) throw new Error("Selected pocket export does not match pose electrostatics");
    const file = bundle.pose_files?.[source.aa];
    if (!bundle.codes?.includes(source.aa) || file?.file !== `pose_electrostatics/electrostatics_${source.aa.toLowerCase()}.tsv.gz`
      || !Number.isSafeInteger(file.rows) || file.rows < 0 || !isPotentialHash(file.sha256)) {
      throw new Error(`No valid ${source.aa} pose electrostatics file in the manifest`);
    }
    // The safe local bundle directory is authoritative; never fetch the manifest's server path.
    const text = await readPocketPotentialFile(`${source.directory}/${file.file}?v=${file.sha256}`, file.sha256);
    entry.table = parsePoseElectrostatics(text, file.rows);
    entry.kind = "ready";
  })().catch((error) => {
    entry.kind = "unavailable";
    entry.error = error.message;
    console.warn("Explorer pose electrostatics unavailable", error);
  }).then(() => entry);
  return entry;
}

function poseElectrostaticsForRow(row, aa = state.aa) {
  const source = poseElectrostaticsSource(row, aa);
  const entry = source && poseElectrostaticsTables.get(source.key);
  if (!source) return { value: NaN, pose: null, status: "unavailable", error: "No matching compact L pocket provenance" };
  if (!entry || entry.kind === "loading") return { value: NaN, pose: null, status: "loading", error: "Loading pose electrostatics" };
  if (entry.kind !== "ready") return { value: NaN, pose: null, status: "unavailable", error: entry.error };
  return entry.table.get(String(row.pocket_id)) || { value: NaN, pose: null, status: "missing", error: "Pocket missing from this AA's pose electrostatics table" };
}

function poseElectrostaticsCell(row) {
  const result = poseElectrostaticsForRow(row);
  const available = Number.isFinite(result.value);
  const poseId = savedSfctPoseId(row);
  const sfctPose = poseId === null ? "" : ` This row's saved SFCT/Combined pose: Vina MODEL ${poseId} (SFCT index ${row.sfct_best_pose}).`;
  const title = available ? POSE_QPHI_DESCRIPTION + sfctPose
    : `qφ (kT) unavailable (${result.status})${result.error ? `: ${result.error}` : ""}`;
  return `<td class="numeric qphi-cell" title="${escapeHTML(title)}">${available ? fmt(result.value, 3) : result.status === "loading" ? "…" : "—"}</td>`;
}

function updateExplorerPoseStatus(entry) {
  const loading = entry.kind === "loading";
  const status = $("#explorer-qphi-status");
  status.hidden = entry.kind === "ready";
  status.dataset.state = entry.kind;
  status.textContent = loading ? "Loading qφ (kT) for Vina pose 1…"
    : `qφ (kT) unavailable — ${entry.error}. Docking results are unchanged.`;
  $("#explorer-qphi-retry").hidden = entry.kind !== "unavailable";
  for (const selector of ["#download-tsv", "#download-csv"]) $(selector).disabled = loading;
}

function prepareExplorerPoseElectrostatics() {
  const aa = state.aa;
  const rows = state.rawByAA.get(aa);
  const source = poseElectrostaticsSource(rows?.[0], aa);
  const entry = ensurePoseElectrostaticsTable(source, aa);
  updateExplorerPoseStatus(entry);
  if (explorerPoseObservation?.entry === entry && explorerPoseObservation.rows === rows) return entry.promise;
  const observation = { aa, rows, entry };
  explorerPoseObservation = observation;
  // Update only the table after loading, preserving page, filters, scores, plots and viewer.
  entry.promise.then(() => {
    if (explorerPoseObservation !== observation || state.aa !== aa || state.rawByAA.get(aa) !== rows) return;
    state.filtered = sortRows(state.filtered);
    renderResults(state.filtered);
    updateExplorerPoseStatus(entry);
  });
  return entry.promise;
}

function retryExplorerPoseElectrostatics() {
  electrostaticsManifestPromise = null;
  poseElectrostaticsBundlePromise = null;
  poseElectrostaticsTables.clear();
  explorerPoseObservation = null;
  void prepareExplorerPoseElectrostatics();
  renderResults(state.filtered);
}
