// P2Rank surface sample points are a visual overlay, never atoms used for scoring.
let pocketPointTableCache = null;
let molstarPocketCloudData = null;
let molstarPocketCloudSphere = null;
let pocketCloudGeneration = 0;

function setPocketCloudStatus(kind, message) {
  const status = $("#structure-cloud-status");
  status.dataset.state = kind;
  status.textContent = message;
  updatePocketPointColorUI();
}

function parsePocketPointTable(text) {
  const headerEnd = text.indexOf("\n");
  const header = (headerEnd < 0 ? text : text.slice(0, headerEnd)).replace(/^\uFEFF/, "").replace(/\r$/, "");
  if (header !== "pocket_id\tpoints") throw new Error("Unexpected pocket-point columns");
  const table = new Map();
  // Index compact strings without constructing millions of atom/point objects at startup.
  for (let start = headerEnd < 0 ? text.length : headerEnd + 1; start < text.length;) {
    const end = text.indexOf("\n", start);
    const line = text.slice(start, end < 0 ? text.length : end).replace(/\r$/, "");
    start = end < 0 ? text.length : end + 1;
    if (!line.trim()) continue;
    const tab = line.indexOf("\t");
    const id = line.slice(0, tab);
    if (tab < 0 || line.indexOf("\t", tab + 1) >= 0 || !/^[1-9]\d*$/.test(id) || table.has(id)) {
      throw new Error("Invalid or duplicate pocket-point ID");
    }
    table.set(id, line.slice(tab + 1).trim());
  }
  return table;
}

function parsePocketPoints(text) {
  if (!text) return [];
  const tokens = text.split(";");
  if (tokens.length > 9999) throw new Error("Too many points for the supported display format");
  return tokens.map((token) => {
    const fields = token.split(",");
    if (fields.length !== 3 || fields.some((field) => !field.trim() || !Number.isFinite(Number(field)))) {
      throw new Error("Invalid P2Rank point coordinates");
    }
    const point = fields.map(Number);
    if (point.some((value) => value.toFixed(3).length > 8)) throw new Error("Point exceeds the PDB coordinate range");
    return point;
  });
}

function pocketPointsPDB(points) {
  // Pseudoatoms are solely a Mol* rendering carrier, following P2Rank's points-PDB approach.
  // Uniform small spheres only: no bonds, inferred surface, recentering or resampling.
  const lines = points.map((point, index) => {
    const serial = String(index + 1).padStart(5);
    const name = `P${index.toString(36).toUpperCase().padStart(3, "0")}`;
    const xyz = point.map((value) => value.toFixed(3).padStart(8)).join("");
    return `HETATM${serial} ${name} PNT Q   1    ${xyz}  1.00  0.00           H  `;
  });
  return [...lines, "END", ""].join("\n");
}

function pocketPointSource(request) {
  const row = request.row;
  // This shared table uses L/pockets.tsv.gz IDs, not D-bundle or legacy IDs.
  if (row._compactDirectory !== `${RESULTS_DIRECTORY}/L` || row._compactAA !== request.aa
    || !AMINO_ACIDS.some((aa) => aa.code === request.aa)
    || !/^[1-9]\d*$/.test(String(row.pocket_id)) || !/^[a-f0-9]{64}$/.test(row._compactPocketHash || "")) return null;
  return `${RESULTS_DIRECTORY}/pocket_points.tsv.gz?v=${row._compactPocketHash}`;
}

async function loadPocketPointTable(path) {
  if (pocketPointTableCache?.path === path) return pocketPointTableCache.promise;
  const entry = { path, promise: loadGzippedResultText(path).then(parsePocketPointTable) };
  pocketPointTableCache = entry;
  try {
    return await entry.promise;
  } catch (error) {
    if (pocketPointTableCache === entry) pocketPointTableCache = null;
    throw error;
  }
}

function isCurrentPocketCloud(request, generation) {
  return generation === pocketCloudGeneration && isCurrentStructureRequest(request) && state.showPocketCloud;
}

async function removeMolstarPocketCloud(viewer) {
  resetPocketPointColorCloud();
  if (molstarPocketCloudData) await viewer.plugin.build().delete(molstarPocketCloudData.ref).commit();
  molstarPocketCloudData = null;
  molstarPocketCloudSphere = null;
}

function pocketPointSphere(points) {
  const center = [0, 1, 2].map((axis) => points.reduce((sum, point) => sum + point[axis], 0) / points.length);
  const radius = points.reduce((max, point) => Math.max(max, Math.hypot(...point.map((value, axis) => value - center[axis]))), 0);
  return { center, radius: radius + 0.55 };
}

async function drawMolstarPocketCloud(viewer, request, generation, points) {
  if (!isCurrentPocketCloud(request, generation)) return;
  await removeMolstarPocketCloud(viewer);
  if (!isCurrentPocketCloud(request, generation)) return;
  if (!points.length) {
    setPocketCloudStatus("missing", "No exported P2Rank points for this pocket");
    return;
  }
  try {
    molstarPocketCloudData = await viewer.plugin.builders.data.rawData({
      data: pocketPointsPDB(points), label: `${request.row.pocket} · P2Rank sample points (not atoms)`,
    });
    if (!isCurrentPocketCloud(request, generation)) return;
    const trajectory = await viewer.plugin.builders.structure.parseTrajectory(molstarPocketCloudData, "pdb");
    if (!isCurrentPocketCloud(request, generation)) return;
    const model = await viewer.plugin.builders.structure.createModel(trajectory);
    const structure = await viewer.plugin.builders.structure.createStructure(model, { name: "model", params: {} });
    if (!isCurrentPocketCloud(request, generation)) return;
    if (structure.obj?.data?.elementCount !== points.length) throw new Error("Pocket point count changed during rendering");
    const representation = await viewer.plugin.builders.structure.representation.addRepresentation(structure, {
      type: "spacefill", color: "uniform", colorParams: { value: state.pocketPointColor === "potential" ? POCKET_POTENTIAL_MISSING_COLOR : 0xffd700 },
      size: "uniform", sizeParams: { value: 0.55 },
      typeParams: { alpha: 0.3, sizeFactor: 1, ignoreHydrogens: false, visuals: ["element-sphere"], bumpFrequency: 0 },
    });
    if (!isCurrentPocketCloud(request, generation)) return;
    molstarPocketCloudSphere = pocketPointSphere(points);
    // Use data bounds, not the previous camera radius (which would grow on every toggle).
    const spheres = [molstarPocketCloudSphere];
    if (molstarLigandSphere) spheres.push(molstarLigandSphere);
    const pocketCenter = [request.row.center_x, request.row.center_y, request.row.center_z];
    if (pocketCenter.every(Number.isFinite)) spheres.push({ center: pocketCenter, radius: 5 });
    viewer.plugin.managers.camera.focusSpheres(spheres, (sphere) => sphere, { minRadius: 8, extraRadius: 2 });
    const expected = request.row.sas_points;
    const countNote = Number.isFinite(expected) && expected !== points.length ? ` (metadata: ${expected})` : "";
    setPocketCloudStatus("ready", `${request.row.pocket} · ${points.length.toLocaleString()} exported points${countNote}`);
    registerPocketPointColorCloud(viewer, request, generation, points, structure, representation);
  } catch (error) {
    await removeMolstarPocketCloud(viewer);
    throw error;
  } finally {
    if (!isCurrentPocketCloud(request, generation)) await removeMolstarPocketCloud(viewer);
  }
}

function pocketCloudError(request, generation, error) {
  if (!isCurrentPocketCloud(request, generation)) return;
  console.warn("P2Rank point cloud unavailable", error);
  setPocketCloudStatus("error", "Cloud unavailable — other layers remain visible. Toggle to retry.");
}

async function prepareMolstarPocketCloud(viewer, request) {
  if (!isCurrentStructureRequest(request)) return;
  const generation = ++pocketCloudGeneration;
  if (!state.showPocketCloud) {
    setPocketCloudStatus("hidden", "Pocket cloud hidden");
    return;
  }
  const path = pocketPointSource(request);
  if (!path) {
    setPocketCloudStatus("missing", "No matching P2Rank point bundle for this result");
    return;
  }
  setPocketCloudStatus("loading", "Loading P2Rank points (shared download, about 21 MB)…");
  try {
    // The large download must not block the model/pocket/ligand update queue.
    const table = await loadPocketPointTable(path);
    if (!isCurrentPocketCloud(request, generation)) return;
    const points = parsePocketPoints(table.get(String(request.row.pocket_id)));
    molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
      try {
        await drawMolstarPocketCloud(viewer, request, generation, points);
      } catch (error) {
        pocketCloudError(request, generation, error);
      }
    });
  } catch (error) {
    pocketCloudError(request, generation, error);
  }
}

function bindPocketCloudEvents() {
  bindPocketPointColorEvents();
  $("#show-pocket-cloud").checked = state.showPocketCloud;
  $("#show-pocket-cloud").addEventListener("change", (event) => {
    state.showPocketCloud = event.target.checked;
    ++pocketCloudGeneration;
    if (!state.showPocketCloud) {
      setPocketCloudStatus("hidden", "Pocket cloud hidden");
      molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
        if (!state.showPocketCloud && molstarViewer) await removeMolstarPocketCloud(molstarViewer);
      });
      return;
    }
    const request = molstarLatestRequest;
    if (request && isCurrentStructureRequest(request) && molstarViewer
      && molstarLoadedModel === dockingModelReference(request.row).id) {
      void prepareMolstarPocketCloud(molstarViewer, request);
    }
  });
}
