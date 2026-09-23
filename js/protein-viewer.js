// Switch only the receptor representation. Pocket, ligand and cloud layers are independent.
let molstarProteinComponent = null;
let molstarProteinRepresentations = {};
let molstarProteinStyle = "cartoon";
let molstarPocketColorSchema = null;
let molstarAppliedProteinColors = { color: "plddt", highlight: true };
const PROTEIN_RESIDUE_COLORS = {
  negative: 0xd1495b, positive: 0x3675c5, histidine: 0xb88732, neutral: 0xb8bdc5, pocket: 0x39ff14,
  polar: 0x4d86b8, nonpolar: 0xc49a6c,
};
const RESIDUE_CHARGE_GROUPS = [
  { names: ["ASP", "GLU"], color: PROTEIN_RESIDUE_COLORS.negative },
  { names: ["LYS", "ARG"], color: PROTEIN_RESIDUE_COLORS.positive },
  { names: ["HIS"], color: PROTEIN_RESIDUE_COLORS.histidine },
];
// Canonical-residue grouping following ProDy's polar/hydrophobic selections.
// This convention includes GLY, CYS and TYR as polar; it is not a hydropathy scale.
// https://www.bahargroup.org/prody/manual/reference/atomic/flags.html#protein
const RESIDUE_POLARITY_GROUPS = [
  { label: "Polar", names: ["ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "LYS", "SER", "THR", "TYR"], color: PROTEIN_RESIDUE_COLORS.polar },
  { label: "Nonpolar", names: ["ALA", "ILE", "LEU", "MET", "PHE", "PRO", "TRP", "VAL"], color: PROTEIN_RESIDUE_COLORS.nonpolar },
];

function resetProteinRepresentations() {
  molstarProteinComponent = null;
  molstarProteinRepresentations = {};
  molstarProteinStyle = "cartoon";
  molstarPocketColorSchema = null;
  molstarAppliedProteinColors = { color: "plddt", highlight: true };
}

function addProteinResiduePaint(update, representation, schema, colorMode, highlight) {
  const lib = window.molstar?.lib;
  const transform = lib?.plugin?.StateTransforms?.Representation?.OverpaintStructureRepresentation3DFromBundle;
  if (!transform) throw new Error("Protein pocket coloring is unavailable");
  // Update only color, preserving geometry, opacity, size, visibility and camera.
  update.to(representation).update((params) => {
    if (colorMode === "plddt") {
      if (params.colorTheme?.name !== "plddt-confidence") {
        params.colorTheme = { name: "plddt-confidence", params: {} };
      }
    } else {
      params.colorTheme = { name: "uniform", params: { value: PROTEIN_RESIDUE_COLORS.neutral, saturation: 0, lightness: 0 } };
    }
  });
  const layers = [];
  const addLayer = (selection, color) => {
    // Overpaint bundles refer to the representation's own root, not the ligand/point models.
    const root = representation.obj?.data?.sourceData?.root;
    if (!root) throw new Error("Protein representation has no source structure");
    const loci = lib.structure.StructureElement.Loci.fromSchema(root, selection);
    if (!lib.structure.StructureElement.Loci.isEmpty(loci)) {
      layers.push({ bundle: lib.structure.StructureElement.Bundle.fromLoci(loci), color, clear: false });
    }
  };
  // Qualitative residue-name classes, not calculated charges, hydropathy or electrostatics.
  const groups = colorMode === "charge" ? RESIDUE_CHARGE_GROUPS
    : colorMode === "polarity" ? RESIDUE_POLARITY_GROUPS : [];
  for (const group of groups) {
    addLayer({ items: group.names.map((name) => ({ label_comp_id: name })) }, group.color);
  }
  // Last layer wins: disabling the pocket override reveals the selected underlying color.
  if (highlight && schema) addLayer(schema, PROTEIN_RESIDUE_COLORS.pocket);
  // Replace all layers so previous pockets and residue-color classes never accumulate.
  update.to(representation).applyOrUpdateTagged("atlas-pocket-color", transform, { layers });
}

async function refreshProteinColors(viewer) {
  const representations = Object.values(molstarProteinRepresentations);
  if (!representations.length) return;
  const color = state.proteinColor;
  const highlight = state.showPocketHighlight;
  const update = viewer.plugin.build();
  for (const representation of representations) {
    addProteinResiduePaint(update, representation, molstarPocketColorSchema, color, highlight);
  }
  await update.commit({ doNotUpdateCurrent: true });
  molstarAppliedProteinColors = { color, highlight };
}

async function setProteinPocketHighlight(viewer, schema) {
  // Keep the current pocket even while its override is off, so it can be restored instantly.
  molstarPocketColorSchema = schema;
  await refreshProteinColors(viewer);
}

function registerProteinRepresentations(loaded) {
  molstarProteinComponent = loaded?.representation?.components?.polymer;
  molstarProteinRepresentations = { cartoon: loaded?.representation?.representations?.polymer };
  if (!molstarProteinComponent || !molstarProteinRepresentations.cartoon) {
    throw new Error("Protein cartoon representation is unavailable");
  }
  molstarProteinStyle = "cartoon";
}

function updateProteinStyleControls(message = "") {
  $$("[data-protein-style]").forEach((button) => {
    button.setAttribute("aria-pressed", String(button.dataset.proteinStyle === state.proteinStyle));
  });
  const status = $("#structure-style-status");
  status.textContent = message;
  status.hidden = !message;
}

async function applyProteinStyle(viewer, request) {
  if (request !== molstarLatestRequest || !molstarProteinComponent) return;
  const style = state.proteinStyle === "surface" ? "surface" : "cartoon";
  const component = molstarProteinComponent;
  try {
    if (style === "surface" && !molstarProteinRepresentations.surface) {
      updateProteinStyleControls("Preparing protein surface…");
      // Build hidden first; keep the cartoon visible until the surface is ready.
      const surface = await viewer.plugin.builders.structure.representation.addRepresentation(component, {
        type: "molecular-surface", color: "plddt-confidence",
        typeParams: { alpha: 0.45, quality: "medium", ignoreHydrogens: true },
      }, { initialState: { isHidden: true }, tag: "atlas-protein-surface" });
      if (!surface) throw new Error("Protein surface could not be created");
      if (component !== molstarProteinComponent) {
        await viewer.plugin.build().delete(surface.ref).commit();
        return;
      }
      molstarProteinRepresentations.surface = surface;
    }
    if (request !== molstarLatestRequest || style !== state.proteinStyle) return;
    // Also color a newly created surface when no pocket override is enabled.
    await refreshProteinColors(viewer);
    if (request !== molstarLatestRequest || style !== state.proteinStyle) return;
    for (const [name, representation] of Object.entries(molstarProteinRepresentations)) {
      viewer.plugin.state.data.updateCellState(representation.ref, { isHidden: name !== style });
    }
    molstarProteinStyle = style;
    updateProteinStyleControls();
  } catch (error) {
    if (request !== molstarLatestRequest || style !== state.proteinStyle) return;
    console.warn("Protein representation change failed", error);
    state.proteinStyle = molstarProteinStyle;
    updateProteinStyleControls("Could not switch protein style. Previous view kept; click to retry.");
  }
}

function setPocketSticksVisibility(viewer) {
  if (!molstarPocketSticks) return;
  viewer.plugin.state.data.updateCellState(molstarPocketSticks.ref, { isHidden: !state.showPocketSticks });
}

function updateProteinColorControls(message = "") {
  $("#protein-color").value = state.proteinColor;
  $("#highlight-pocket-residues").checked = state.showPocketHighlight;
  $("#structure-charge-note").hidden = state.proteinColor !== "charge";
  $("#structure-polarity-note").hidden = state.proteinColor !== "polarity";
  $("#protein-color").setAttribute("aria-describedby", state.proteinColor === "charge" ? "structure-charge-note"
    : state.proteinColor === "polarity" ? "structure-polarity-note" : "");
  const status = $("#structure-color-status");
  status.textContent = message;
  status.hidden = !message;
  const key = (color, label, description = "") => `<span${description ? ` title="${escapeHTML(description)}"` : ""}><i style="background:#${color.toString(16).padStart(6, "0")}"></i>${label}</span>`;
  let legend;
  if (state.proteinColor === "charge") {
    legend = `<strong>Charge</strong>${key(PROTEIN_RESIDUE_COLORS.negative, "ASP/GLU −")}${key(PROTEIN_RESIDUE_COLORS.positive, "LYS/ARG +")}${key(PROTEIN_RESIDUE_COLORS.histidine, "HIS (variable)")}${key(PROTEIN_RESIDUE_COLORS.neutral, "Other")}`;
  } else if (state.proteinColor === "polarity") {
    legend = `<strong>Polarity</strong>${RESIDUE_POLARITY_GROUPS.map((group) => key(group.color, group.label, group.names.join(", "))).join("")}${key(PROTEIN_RESIDUE_COLORS.neutral, "Other/unknown", "Unrecognized or modified residues are not assigned a polarity class.")}`;
  } else if (state.proteinColor === "none") {
    legend = `<strong>Uniform</strong>${key(PROTEIN_RESIDUE_COLORS.neutral, "Protein")}`;
  } else {
    legend = '<strong>pLDDT</strong><span><i class="confidence-very-high"></i>≥90</span><span><i class="confidence-high"></i>70–90</span><span><i class="confidence-low"></i>50–70</span><span><i class="confidence-very-low"></i>&lt;50</span>';
  }
  if (state.showPocketHighlight || state.showPocketSticks) {
    const label = state.showPocketHighlight ? "Selected pocket" : "Pocket sticks";
    legend += `<span class="pocket-highlight-key"><i class="confidence-pocket"></i>${label}</span>`;
  }
  $("#structure-protein-legend").innerHTML = legend;
}

function bindProteinColorEvents() {
  updateProteinColorControls();
  const applyColors = () => {
    updateProteinColorControls();
    const request = molstarLatestRequest;
    if (!request || !isCurrentStructureRequest(request) || !molstarViewer) return;
    molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
      if (!isCurrentStructureRequest(request)) return;
      const color = state.proteinColor, highlight = state.showPocketHighlight;
      try {
        await refreshProteinColors(molstarViewer);
      } catch (error) {
        console.warn("Protein coloring failed", error);
        if (!isCurrentStructureRequest(request) || color !== state.proteinColor || highlight !== state.showPocketHighlight) return;
        state.proteinColor = molstarAppliedProteinColors.color;
        state.showPocketHighlight = molstarAppliedProteinColors.highlight;
        // Restore the previous colors as well if a failed commit partially updated a layer.
        try { await refreshProteinColors(molstarViewer); } catch (restoreError) { console.warn("Protein color restore failed", restoreError); }
        updateProteinColorControls("Could not apply protein colors. Please try again.");
      }
    });
  };
  $("#protein-color").addEventListener("change", (event) => {
    if (!["plddt", "charge", "polarity", "none"].includes(event.target.value)) return;
    state.proteinColor = event.target.value;
    applyColors();
  });
  $("#highlight-pocket-residues").addEventListener("change", (event) => {
    state.showPocketHighlight = event.target.checked;
    applyColors();
  });
}

function bindPocketSticksEvents() {
  $("#show-pocket-sticks").checked = state.showPocketSticks;
  $("#show-pocket-sticks").addEventListener("change", (event) => {
    state.showPocketSticks = event.target.checked;
    updateProteinColorControls();
    const request = molstarLatestRequest;
    if (!request || !isCurrentStructureRequest(request) || !molstarViewer) return;
    molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(() => {
      if (isCurrentStructureRequest(request)) setPocketSticksVisibility(molstarViewer);
    });
  });
}

function bindProteinStyleEvents() {
  updateProteinStyleControls();
  $$("[data-protein-style]").forEach((button) => button.addEventListener("click", () => {
    const style = button.dataset.proteinStyle;
    if (!["cartoon", "surface"].includes(style)) return;
    state.proteinStyle = style;
    updateProteinStyleControls();
    const request = molstarLatestRequest;
    // During initial loading, performMolstarUpdate will apply the requested style when ready.
    if (!request || !isCurrentStructureRequest(request) || !molstarViewer) return;
    molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
      if (!isCurrentStructureRequest(request)) return;
      await applyProteinStyle(molstarViewer, request);
    });
  }));
}
