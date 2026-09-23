// Optional coordinate overlays only: never used in scoring, ranking or QC statistics.
const LIGAND_POSITION_CACHE_LIMIT = 2;
const ligandPositionCache = new Map();
let molstarLigandData = null;
let molstarLigandSphere = null;
let ligandPoseGeneration = 0;

function setLigandViewerStatus(kind, message) {
  const status = $("#structure-ligand-status");
  status.dataset.state = kind;
  status.textContent = message;
  $("#structure-ligand-legend").hidden = kind !== "ready";
}

function parseLigandPositionTable(text) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
  if (lines.shift() !== "pocket_id\tpose_id\tpositions") throw new Error("Expected pocket_id, pose_id and positions columns");
  const positions = new Map();
  for (const line of lines) {
    if (!line.trim()) continue;
    const columns = line.split("\t");
    const [id, pose, coordinates] = columns;
    const key = `${id}:${pose}`;
    if (columns.length !== 3 || !/^[1-9]\d*$/.test(id) || !/^[1-9]\d*$/.test(pose)
      || !Number.isSafeInteger(Number(id)) || !Number.isSafeInteger(Number(pose)) || positions.has(key)) {
      throw new Error("Invalid or duplicate ligand pocket/pose ID");
    }
    // Keep compact strings in the cache; parse atoms only for the displayed pocket.
    positions.set(key, coordinates.trim());
  }
  return positions;
}

function ligandPoseSelection(request, mode = state.ligandPoseMode) {
  const metric = mode === "vina" ? "vina_affinity" : mode === "sfct" ? "sfct_score" : request.ligandMetric;
  if (!METRICS[metric]) return { error: "No score selected for the ligand pose" };
  if (!hasUsableScore(request.row, metric)) return { error: `No successful ${METRICS[metric].label} score for this pocket` };
  const vina = metric === "vina_affinity";
  const poseId = vina ? 1 : savedSfctPoseId(request.row);
  if (poseId === null) return { error: "Saved SFCT/Combined pose index is missing or invalid" };
  const manual = mode === "vina" || mode === "sfct";
  const label = manual ? `${vina ? "Vina-best" : "Saved SFCT/Combined"} (manual)`
    : `${METRICS[metric].label}${vina ? "" : " (saved pose)"}`;
  return { poseId, label: `Vina MODEL ${poseId} · ${label}` };
}

function syncLigandPoseControls(row = molstarLatestRequest?.row) {
  const select = $("#ligand-pose-select");
  select.value = state.ligandPoseMode;
  select.disabled = !row;
  $("#ligand-vina-pose-option").disabled = !hasUsableScore(row, "vina_affinity");
  const saved = $("#ligand-sfct-pose-option");
  const poseId = savedSfctPoseId(row);
  const available = poseId !== null && hasUsableScore(row, "sfct_score");
  saved.textContent = available ? `Saved SFCT/Combined - MODEL ${poseId}` : "Saved SFCT/Combined - unavailable";
  saved.disabled = !available;
}

function updateLigandPoseNote(request) {
  const manual = state.ligandPoseMode !== "auto";
  const selected = ligandPoseSelection(request);
  const automatic = ligandPoseSelection(request, "auto");
  const mismatch = manual && selected.poseId && automatic.poseId && selected.poseId !== automatic.poseId;
  const note = manual ? `Manual pose; scores and plots are unchanged.${mismatch ? ` Different pose from ${METRICS[request.ligandMetric].label} (Value type).` : ""}`
    : "Auto follows Value type; SFCT/Combined use the saved selection.";
  $("#structure-ligand-note").textContent = `${note} qφ (kT) remains Vina MODEL 1. Hydrogens hidden.`;
}

function parseLigandAtoms(positions) {
  if (!positions) return [];
  const atoms = positions.split(";").map((token) => {
    const parts = token.split(":");
    const coordinates = parts[1]?.split(",");
    if (parts.length !== 2 || !/^[CHNOS]$/.test(parts[0]) || coordinates?.length !== 3
      || coordinates.some((value) => !value.trim() || !Number.isFinite(Number(value)))) {
      throw new Error("Invalid ligand atom coordinates");
    }
    const [x, y, z] = coordinates.map(Number);
    // PDB has fixed-width coordinate fields. Never silently truncate or shift a position.
    if ([x, y, z].some((value) => value.toFixed(3).length > 8)) {
      throw new Error("Ligand coordinates exceed the supported PDB range");
    }
    return { element: parts[0], x, y, z };
  });
  if (atoms.length > 128 || !atoms.some((atom) => atom.element !== "H")) {
    throw new Error("Unexpected amino-acid atom count");
  }
  return atoms;
}

function ligandPDB(atoms) {
  // LIG and generated names avoid pretending that the export contains canonical atom names.
  // Mol* infers connectivity from these coordinates; bond orders/charges are not supplied.
  const lines = atoms.map((atom, i) => {
    const serial = String(i + 1).padStart(5);
    const name = `${atom.element}${i + 1}`.padEnd(4);
    const xyz = [atom.x, atom.y, atom.z].map((value) => value.toFixed(3).padStart(8)).join("");
    return `HETATM${serial} ${name} LIG Z   1    ${xyz}  1.00  0.00          ${atom.element.padStart(2)}  `;
  });
  return [...lines, "END", ""].join("\n");
}

function ligandPositionSource(request) {
  const aa = AMINO_ACIDS.find((item) => item.code === request.aa);
  const row = request.row;
  // IDs are local to a compact bundle. Never join L positions to D or legacy raw IDs.
  if (!aa || row._compactAA !== aa.code || row._compactDirectory !== `${RESULTS_DIRECTORY}/L`
    || !/^[1-9]\d*$/.test(String(row.pocket_id)) || !/^[a-f0-9]{64}$/.test(row._compactPocketHash || "")) return null;
  return `${row._compactDirectory}/aa_positions/positions_${aa.file}.tsv.gz?v=${row._compactPocketHash}&poses=2`;
}

async function loadLigandPositions(path) {
  if (ligandPositionCache.has(path)) {
    const cached = ligandPositionCache.get(path);
    ligandPositionCache.delete(path);
    ligandPositionCache.set(path, cached);
    return cached;
  }
  const pending = loadGzippedResultText(path).then(parseLigandPositionTable);
  ligandPositionCache.set(path, pending);
  while (ligandPositionCache.size > LIGAND_POSITION_CACHE_LIMIT) {
    ligandPositionCache.delete(ligandPositionCache.keys().next().value);
  }
  try {
    return await pending;
  } catch (error) {
    // An unavailable/incomplete file can be retried by toggling Show ligand.
    if (ligandPositionCache.get(path) === pending) ligandPositionCache.delete(path);
    throw error;
  }
}

function isCurrentStructureRequest(request) {
  return request === molstarLatestRequest && state.currentView === "protein"
    && request.aa === state.aa && request.protein === state.selectedProtein;
}

async function removeMolstarLigand(viewer) {
  molstarLigandSphere = null;
  if (!molstarLigandData) return;
  const ref = molstarLigandData.ref;
  // Remove only our ligand subtree, never the protein or its pLDDT representation.
  await viewer.plugin.build().delete(ref).commit();
  molstarLigandData = null;
}

function focusMolstarLigand(viewer, atoms, row, preserveCamera = false) {
  const heavy = atoms.filter((atom) => atom.element !== "H");
  const center = ["x", "y", "z"].map((axis) => heavy.reduce((sum, atom) => sum + atom[axis], 0) / heavy.length);
  const radius = Math.max(...heavy.map((atom) => Math.hypot(atom.x - center[0], atom.y - center[1], atom.z - center[2])));
  molstarLigandSphere = { center, radius };
  if (preserveCamera) return;
  const spheres = [molstarLigandSphere];
  if (typeof molstarPocketCloudSphere !== "undefined" && molstarPocketCloudSphere) spheres.push(molstarPocketCloudSphere);
  const pocketCenter = [row.center_x, row.center_y, row.center_z];
  if (pocketCenter.every(Number.isFinite)) spheres.push({ center: pocketCenter, radius: 5 });
  // Adjust only the camera: coordinates are never centered, aligned or minimized separately.
  viewer.plugin.managers.camera.focusSpheres(spheres, (sphere) => sphere, { minRadius: 8, extraRadius: 4 });
}

async function updateMolstarLigand(viewer, request, { preserveCamera = false } = {}) {
  if (!isCurrentStructureRequest(request)) return;
  const generation = ++ligandPoseGeneration;
  const isCurrent = () => isCurrentStructureRequest(request) && generation === ligandPoseGeneration;
  await removeMolstarLigand(viewer);
  if (!isCurrent()) return;
  syncLigandPoseControls(request.row);
  updateLigandPoseNote(request);
  if (!state.showLigand) {
    setLigandViewerStatus("hidden", "Ligand hidden");
    return;
  }
  const path = ligandPositionSource(request);
  if (!path) {
    setLigandViewerStatus("missing", "No matching L-AA coordinate bundle for this result");
    return;
  }
  const selection = ligandPoseSelection(request);
  if (selection.error) {
    setLigandViewerStatus("missing", selection.error);
    return;
  }
  setLigandViewerStatus("loading", `Loading ${request.aa} · ${selection.label}…`);
  try {
    const positions = await loadLigandPositions(path);
    if (!isCurrent() || !state.showLigand) return;
    const atoms = parseLigandAtoms(positions.get(`${request.row.pocket_id}:${selection.poseId}`));
    if (!atoms.length) {
      setLigandViewerStatus("missing", `${request.aa} · ${request.row.pocket}: Vina MODEL ${selection.poseId} coordinates unavailable; no other pose substituted`);
      return;
    }
    const label = `${request.aa} · ${request.row.pocket} · ${selection.label}`;
    const data = await viewer.plugin.builders.data.rawData({ data: ligandPDB(atoms), label });
    molstarLigandData = data;
    try {
      if (!isCurrent() || !state.showLigand) return;
      const trajectory = await viewer.plugin.builders.structure.parseTrajectory(data, "pdb");
      if (!isCurrent() || !state.showLigand) return;
      const model = await viewer.plugin.builders.structure.createModel(trajectory);
      const structure = await viewer.plugin.builders.structure.createStructure(model, { name: "model", params: {} });
      if (!isCurrent() || !state.showLigand) return;
      await viewer.plugin.builders.structure.representation.addRepresentation(structure, {
        type: "ball-and-stick", typeParams: { sizeFactor: 0.3, ignoreHydrogens: true }, color: "element-symbol",
        // Mol* defaults to chain-colored carbons and lightened elements; match our legend exactly.
        colorParams: { carbonColor: { name: "element-symbol", params: {} }, saturation: 0, lightness: 0 },
      });
      if (!isCurrent() || !state.showLigand) return;
      focusMolstarLigand(viewer, atoms, request.row, preserveCamera);
      setLigandViewerStatus("ready", `${label} · ${atoms.filter((atom) => atom.element !== "H").length} heavy atoms`);
    } finally {
      if (!isCurrent() || !state.showLigand) await removeMolstarLigand(viewer);
    }
  } catch (error) {
    await removeMolstarLigand(viewer);
    if (!isCurrent()) return;
    console.warn("Ligand overlay unavailable", error);
    setLigandViewerStatus("error", "Ligand unavailable — protein remains visible. Toggle Show ligand to retry.");
  }
}

function queueLigandOverlayUpdate() {
  const request = molstarLatestRequest;
  const generation = ++ligandPoseGeneration; // Invalidate an in-flight pose, not the protein/cloud request.
  syncLigandPoseControls(request?.row);
  if (!request || !isCurrentStructureRequest(request)) return;
  updateLigandPoseNote(request);
  setLigandViewerStatus(state.showLigand ? "loading" : "hidden", state.showLigand ? "Updating ligand pose…" : "Ligand hidden");
  molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
    if (generation !== ligandPoseGeneration || !isCurrentStructureRequest(request) || !molstarViewer
      || molstarLoadedModel !== dockingModelReference(request.row).id) return;
    await updateMolstarLigand(molstarViewer, request, { preserveCamera: true });
  });
}

function bindLigandViewerEvents() {
  syncLigandPoseControls();
  $("#show-pocket-ligand").checked = state.showLigand;
  $("#show-pocket-ligand").addEventListener("change", (event) => {
    state.showLigand = event.target.checked;
    queueLigandOverlayUpdate();
  });
  $("#ligand-pose-select").addEventListener("change", (event) => {
    if (!["auto", "vina", "sfct"].includes(event.target.value)) return;
    state.ligandPoseMode = event.target.value;
    queueLigandOverlayUpdate();
  });
}
