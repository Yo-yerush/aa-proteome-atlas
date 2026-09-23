// Optional, precomputed receptor potentials sampled at the original P2Rank points.
// No scoring changes, client-side electrostatics calculation, or protein recoloring.
const POCKET_POTENTIAL_MISSING_COLOR = 0x9aa0a6;
const POCKET_POTENTIAL_THEME = "atlas-pocket-electrostatics";
const POCKET_POTENTIAL_COLUMNS = ["pocket_rank", "point_index", "phi_kT_per_e"];
const POCKET_POTENTIAL_SUMMARY_COLUMNS = ["protein", "pocket_rank", "n_points", "points_sha256",
  "phi_mean", "phi_min", "phi_max", "fraction_positive", "fraction_negative", "status", "error"];
let pocketElectrostaticsMetadata = null;
let electrostaticsManifestPromise = null;
const pocketPotentialFiles = new Map();
const pocketPotentialThemes = new WeakMap();
const pocketPotentialRegistries = new WeakSet();
let pocketPotentialCloud = null;
let pocketPotentialGeneration = 0;
let pocketPotentialStatus = { kind: "idle", message: "" };

function isPotentialHash(value) { return /^[a-f0-9]{64}$/.test(value || ""); }

async function pocketPotentialSHA256(bytes) {
  if (!globalThis.crypto?.subtle) throw new Error("Checksum validation requires HTTPS or localhost");
  const digest = await globalThis.crypto.subtle.digest("SHA-256", bytes);
  return [...new Uint8Array(digest)].map((byte) => byte.toString(16).padStart(2, "0")).join("");
}

async function readPocketPotentialFile(path, checksum, gzip = true) {
  const response = await fetch(path, { cache: checksum ? "default" : "no-cache" });
  if (!response.ok) throw new Error(`Missing electrostatics file (HTTP ${response.status}): ${path.split("?")[0]}`);
  const bytes = new Uint8Array(await response.arrayBuffer());
  // This pipeline hashes compressed bytes, unlike the compact score export's text hashes.
  if (checksum && await pocketPotentialSHA256(bytes) !== checksum) {
    throw new Error("File checksum mismatch; publish matching files and serve .gz files without HTTP Content-Encoding");
  }
  if (gzip && bytes[0] === 0x1f && bytes[1] === 0x8b) {
    if (typeof DecompressionStream === "undefined") throw new Error("This browser cannot decompress electrostatics files");
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
    return new TextDecoder("utf-8", { fatal: true }).decode(await new Response(stream).arrayBuffer());
  }
  if (gzip && checksum) throw new Error("Expected a gzip electrostatics file");
  return new TextDecoder("utf-8", { fatal: true }).decode(bytes);
}

function validatePocketPotentialManifest(manifest) {
  const points = manifest?.points;
  const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
  if (manifest?.format !== "aa-pocket-electrostatics" || manifest.version !== 1
    || !same(points?.columns, POCKET_POTENTIAL_COLUMNS) || points.point_index_base !== 0
    || !same(points.join, ["protein", "pocket_rank"]) || points.potential_units !== "kT/e"
    || points.order !== "ATOM/HETATM encounter order within the original P2Rank pocket rank"
    || points.coordinate_fingerprint !== "SHA256 of ASCII x,y,z;x,y,z;...; fixed 3 decimals, -0.000 normalized to 0.000; no newline"
    || points.file_template !== "points_{protein}.tsv.gz"
    || !same(manifest.summary?.columns, POCKET_POTENTIAL_SUMMARY_COLUMNS)
    || manifest.summary.file !== "pocket_summary.tsv.gz" || !isPotentialHash(manifest.summary.sha256)
    || !Array.isArray(manifest.bundles) || !manifest.bundles.length || !manifest.proteins) {
    throw new Error("Unsupported electrostatics manifest or point-order convention");
  }
  const scale = points.color_scale;
  if (scale?.minimum !== -5 || scale.midpoint !== 0 || scale.maximum !== 5
    || !same(scale.colors?.map((color) => String(color).toLowerCase()), ["#d73027", "#ffffff", "#4575b4"])) {
    throw new Error("Unsupported electrostatic color scale; expected fixed -5 / 0 / +5 kT/e");
  }
  // Discard pose tables/ligand charge models: this feature needs only point metadata.
  return { proteins: manifest.proteins, bundles: manifest.bundles.map(({ manifest_sha256 }) => manifest_sha256),
    summary: manifest.summary, scale };
}

async function loadElectrostaticsManifest() {
  if (electrostaticsManifestPromise) return electrostaticsManifestPromise;
  const pending = (async () => {
    await pocketPotentialSHA256(new Uint8Array());
    const manifest = JSON.parse(await readPocketPotentialFile(`${RESULTS_DIRECTORY}/electrostatics/manifest.json.gz`));
    if (manifest?.format !== "aa-pocket-electrostatics" || manifest.version !== 1) {
      throw new Error("Unsupported electrostatics manifest");
    }
    // Shared by the point viewer and Explorer. Do not retain embedded ligand molblocks.
    const { format, version, bundles, proteins, points, summary, poses } = manifest;
    return { format, version, bundles, proteins, points, summary, poses };
  })();
  electrostaticsManifestPromise = pending;
  try { return await pending; }
  catch (error) {
    if (electrostaticsManifestPromise === pending) electrostaticsManifestPromise = null;
    throw error;
  }
}

function pocketPotentialRows(text, columns) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
  if (lines.shift() !== columns.join("\t")) throw new Error("Unexpected electrostatics table columns");
  return lines.filter((line) => line !== "").map((line) => {
    const fields = line.split("\t");
    if (fields.length !== columns.length) throw new Error("Malformed electrostatics table row");
    return fields;
  });
}

function parsePocketPotentialSummary(text) {
  const summary = new Map();
  for (const fields of pocketPotentialRows(text, POCKET_POTENTIAL_SUMMARY_COLUMNS)) {
    const [model, rank, count, hash] = fields;
    const key = `${model}|${rank}`;
    if (!/^[1-9]\d*$/.test(rank) || !/^\d+$/.test(count) || !Number.isSafeInteger(Number(count))
      || !isPotentialHash(hash) || summary.has(key)) throw new Error("Invalid or duplicate electrostatics pocket metadata");
    summary.set(key, { count: Number(count), hash, status: fields[9], error: fields[10] });
  }
  return summary;
}

async function loadPocketElectrostaticsMetadata() {
  if (pocketElectrostaticsMetadata) return pocketElectrostaticsMetadata;
  const directory = `${RESULTS_DIRECTORY}/electrostatics`;
  const pending = (async () => {
    const manifest = validatePocketPotentialManifest(await loadElectrostaticsManifest());
    const [summaryText, compactText] = await Promise.all([
      readPocketPotentialFile(`${directory}/${manifest.summary.file}?v=${manifest.summary.sha256}`, manifest.summary.sha256),
      readPocketPotentialFile(`${RESULTS_DIRECTORY}/L/manifest.json`, null, false),
    ]);
    // Bridge the two hash conventions through the exact compact manifest used by the pipeline.
    const compactHash = await pocketPotentialSHA256(new TextEncoder().encode(compactText));
    const compact = JSON.parse(compactText);
    if (!manifest.bundles.includes(compactHash) || compact.format !== "aa-proteome-atlas-compact"
      || compact.version !== 1 || !isPotentialHash(compact.pockets?.sha256)) {
      throw new Error("Electrostatics was generated from a different compact L bundle");
    }
    return { ...manifest, pocketHash: compact.pockets.sha256, summary: parsePocketPotentialSummary(summaryText) };
  })();
  pocketElectrostaticsMetadata = pending;
  try { return await pending; }
  catch (error) {
    if (pocketElectrostaticsMetadata === pending) pocketElectrostaticsMetadata = null;
    throw error;
  }
}

function parsePocketPotentials(text) {
  const pockets = new Map();
  for (const [rank, index, field] of pocketPotentialRows(text, POCKET_POTENTIAL_COLUMNS)) {
    if (!/^[1-9]\d*$/.test(rank) || !/^(0|[1-9]\d*)$/.test(index)) throw new Error("Invalid potential point index");
    if (!pockets.has(rank)) pockets.set(rank, []);
    const values = pockets.get(rank);
    if (Number(index) !== values.length) throw new Error("Potential point indices are duplicated, missing or out of order");
    // Number('') is zero! Explicitly preserve missing samples, never silently whiten them.
    const value = field.trim() === "" ? null : Number(field);
    if (value !== null && (!/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?$/i.test(field) || !Number.isFinite(value))) {
      throw new Error("Invalid electrostatic potential sample");
    }
    values.push(value);
  }
  return pockets;
}

async function loadProteinPocketPotentials(model, entry) {
  if (entry?.points_file !== `points_${model}.tsv.gz` || !isPotentialHash(entry.points_file_sha256)) {
    throw new Error("Missing or invalid potential file for the recorded model");
  }
  const path = `${RESULTS_DIRECTORY}/electrostatics/${entry.points_file}?v=${entry.points_file_sha256}`;
  let pending = pocketPotentialFiles.get(path);
  if (pending) pocketPotentialFiles.delete(path);
  else pending = readPocketPotentialFile(path, entry.points_file_sha256).then(parsePocketPotentials);
  pocketPotentialFiles.set(path, pending);
  // Bounded demand-only cache: never preload the proteome's potential files.
  while (pocketPotentialFiles.size > 2) pocketPotentialFiles.delete(pocketPotentialFiles.keys().next().value);
  try { return await pending; }
  catch (error) {
    if (pocketPotentialFiles.get(path) === pending) pocketPotentialFiles.delete(path);
    throw error;
  }
}

function pocketCoordinateFingerprintText(points) {
  return points.map((point) => {
    if (point.length !== 3 || !point.every(Number.isFinite)) throw new Error("Invalid point coordinates");
    return point.map((value) => {
      const fixed = value.toFixed(3);
      return fixed === "-0.000" ? "0.000" : fixed;
    }).join(",");
  }).join(";");
}

async function loadValidatedPocketPotentials(request, points, stillCurrent = () => true) {
  const model = dockingModelReference(request.row).id;
  const rank = request.row.rank;
  if (!pocketPointSource(request) || !Number.isSafeInteger(rank) || rank < 1
    || request.row.pocket !== `pocket${rank}`) throw new Error("The selected pocket lacks matching model/rank provenance");
  const metadata = await loadPocketElectrostaticsMetadata();
  if (!stillCurrent()) return null;
  if (metadata.pocketHash !== request.row._compactPocketHash) throw new Error("Pocket export checksum does not match electrostatics");
  const pocket = metadata.summary.get(`${model}|${rank}`);
  if (!pocket) throw new Error("No electrostatics result for this model and pocket");
  if (pocket.status !== "success") throw new Error(`Calculation unavailable (${pocket.status || "missing status"})${pocket.error ? `: ${pocket.error}` : ""}`);
  if (!points.length || pocket.count !== points.length) throw new Error("P2Rank point count does not match electrostatics");
  const fingerprint = await pocketPotentialSHA256(new TextEncoder().encode(pocketCoordinateFingerprintText(points)));
  if (!stillCurrent()) return null;
  if (fingerprint !== pocket.hash) throw new Error("P2Rank coordinate/order checksum does not match electrostatics");
  const entry = metadata.proteins[model];
  const pockets = await loadProteinPocketPotentials(model, entry);
  if (!stillCurrent()) return null;
  const values = pockets.get(String(rank));
  if (!values || values.length !== points.length) throw new Error("Potential sample count does not match the selected pocket");
  const missing = values.filter((value) => value === null).length;
  if (missing === values.length) throw new Error("All potential samples are missing for this pocket");
  return { values, missing, scale: metadata.scale, fragment: entry.fragment_only === true };
}

function pocketPotentialColor(value, scale) {
  if (!Number.isFinite(value)) return POCKET_POTENTIAL_MISSING_COLOR;
  const lower = value <= scale.midpoint;
  const start = lower ? scale.minimum : scale.midpoint;
  const end = lower ? scale.midpoint : scale.maximum;
  const fraction = Math.max(0, Math.min(1, (value - start) / (end - start)));
  const a = parseInt(scale.colors[lower ? 0 : 1].slice(1), 16);
  const b = parseInt(scale.colors[lower ? 1 : 2].slice(1), 16);
  return [16, 8, 0].reduce((color, shift) => color | Math.round(((a >> shift) & 255) * (1 - fraction)
    + ((b >> shift) & 255) * fraction) << shift, 0);
}

function pocketPotentialTheme(ctx, props) {
  const data = pocketPotentialThemes.get(ctx.structure?.root);
  return { factory: pocketPotentialTheme, granularity: "group", props, contextHash: data?.revision || 0,
    description: "Receptor potential at P2Rank points (kT/e); gray indicates unavailable samples.",
    color: (location) => {
      if (location.kind !== "element-location") return POCKET_POTENTIAL_MISSING_COLOR;
      const name = location.unit.model.atomicHierarchy.atoms.label_atom_id.value(location.element);
      const index = /^P[0-9A-Z]{3}$/.test(name) ? parseInt(name.slice(1), 36) : -1;
      return data?.colors[index] ?? POCKET_POTENTIAL_MISSING_COLOR;
    } };
}

function registerPocketPotentialTheme(viewer, cloud, result, revision) {
  const root = cloud.structure.obj?.data?.root;
  const registry = viewer.plugin.representation?.structure?.themes?.colorThemeRegistry;
  if (!root || !registry) throw new Error("Mol* point-color theme is unavailable");
  // PDB parsing may reorder atoms. Use our unique point IDs rather than Mol* atom offsets.
  const indices = new Set();
  for (const unit of root.units) {
    for (const element of unit.elements) {
      const name = unit.model.atomicHierarchy.atoms.label_atom_id.value(element);
      const index = /^P[0-9A-Z]{3}$/.test(name) ? parseInt(name.slice(1), 36) : -1;
      if (index < 0 || index >= result.values.length || indices.has(index)) throw new Error("Rendered point IDs do not match potential indices");
      indices.add(index);
    }
  }
  if (indices.size !== result.values.length) throw new Error("Rendered point IDs are incomplete");
  pocketPotentialThemes.set(root, { revision, colors: result.values.map((value) => pocketPotentialColor(value, result.scale)) });
  if (!pocketPotentialRegistries.has(registry)) {
    registry.add({ name: POCKET_POTENTIAL_THEME, label: "Pocket electrostatic potential", category: "Miscellaneous",
      factory: pocketPotentialTheme, getParams: () => ({}), defaultValues: {},
      isApplicable: (ctx) => !!ctx.structure && pocketPotentialThemes.has(ctx.structure.root) });
    pocketPotentialRegistries.add(registry);
  }
}

function updatePocketPointColorUI() {
  const potential = state.pocketPointColor === "potential";
  const ready = $("#structure-cloud-status").dataset.state === "ready";
  const colored = ["ready", "partial"].includes(pocketPotentialStatus.kind);
  $("#pocket-point-color").value = state.pocketPointColor;
  $("#structure-cloud-legend").hidden = !ready || potential;
  $("#pocket-potential-legend").hidden = !ready || !potential || !colored;
  $("#pocket-potential-missing-key").hidden = !ready || !potential || (colored && !pocketPotentialStatus.missing);
  const status = $("#pocket-potential-status");
  status.hidden = !ready || !potential;
  status.dataset.state = pocketPotentialStatus.kind;
  status.textContent = pocketPotentialStatus.message;
  $("#pocket-potential-retry").hidden = !ready || !potential || pocketPotentialStatus.kind !== "unavailable";
}

function resetPocketPointColorCloud() {
  ++pocketPotentialGeneration;
  pocketPotentialCloud = null;
  pocketPotentialStatus = { kind: "idle", message: "" };
  updatePocketPointColorUI();
}

function registerPocketPointColorCloud(viewer, request, generation, points, structure, representation) {
  pocketPotentialCloud = { viewer, request, generation, points, structure, representation };
  if (state.pocketPointColor === "potential") void refreshPocketPointColor();
  else updatePocketPointColorUI();
}

async function refreshPocketPointColor() {
  const revision = ++pocketPotentialGeneration;
  const cloud = pocketPotentialCloud;
  const mode = state.pocketPointColor;
  const current = () => cloud && cloud === pocketPotentialCloud && revision === pocketPotentialGeneration
    && mode === state.pocketPointColor && isCurrentPocketCloud(cloud.request, cloud.generation);
  if (!current()) { updatePocketPointColorUI(); return; }
  const paint = (result) => {
    const work = molstarUpdateQueue.catch(() => {}).then(async () => {
      if (!current()) return;
      if (result) registerPocketPotentialTheme(cloud.viewer, cloud, result, revision);
      const theme = result ? { name: POCKET_POTENTIAL_THEME, params: {} }
        : { name: "uniform", params: { value: mode === "gold" ? 0xffd700 : POCKET_POTENTIAL_MISSING_COLOR, saturation: 0, lightness: 0 } };
      const update = cloud.viewer.plugin.build();
      // Only change this existing representation's color theme. No camera/geometry updates.
      update.to(cloud.representation).update((params) => { params.colorTheme = theme; });
      await update.commit({ doNotUpdateCurrent: true });
    });
    molstarUpdateQueue = work.catch(() => {});
    return work;
  };
  pocketPotentialStatus = { kind: mode === "gold" ? "idle" : "loading", message: "Loading and validating electrostatic potential…" };
  updatePocketPointColorUI();
  try {
    await paint(null); // Gray while pending; white is reserved for a measured zero potential.
    if (!current() || mode === "gold") return;
    // Network and hashing stay OUTSIDE the Mol* queue, keeping all other controls responsive.
    const result = await loadValidatedPocketPotentials(cloud.request, cloud.points, current);
    if (!current()) return;
    await paint(result);
    if (!current()) return;
    pocketPotentialStatus = { kind: result.missing ? "partial" : "ready", missing: result.missing,
      message: `Validated receptor potential · ${result.values.length - result.missing}/${result.values.length} points`
        + (result.missing ? ` · ${result.missing} unavailable (gray)` : "")
        + " · Colors saturate at ±5 kT/e" + (result.fragment ? " · Fragment-only model" : "") };
  } catch (error) {
    if (!current()) return;
    console.warn("Pocket electrostatic potential unavailable", error);
    pocketPotentialStatus = { kind: "unavailable", message: `Electrostatic potential unavailable — ${error.message}. Choose Gold or retry.` };
    try { await paint(null); }
    catch (paintError) {
      console.warn("Point cloud color update failed", paintError);
      pocketPotentialStatus.message += " Point recoloring failed; reload the pocket.";
    }
  }
  if (current()) updatePocketPointColorUI();
}

function bindPocketPointColorEvents() {
  updatePocketPointColorUI();
  $("#pocket-point-color").addEventListener("change", (event) => {
    if (!["gold", "potential"].includes(event.target.value)) return;
    state.pocketPointColor = event.target.value;
    void refreshPocketPointColor();
  });
  $("#pocket-potential-retry").addEventListener("click", () => {
    // Explicit retry also rechecks newly published metadata after a failed calculation/export.
    electrostaticsManifestPromise = null;
    pocketElectrostaticsMetadata = null;
    pocketPotentialFiles.clear();
    void refreshPocketPointColor();
  });
}
