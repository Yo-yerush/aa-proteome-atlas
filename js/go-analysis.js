// Local GO overrepresentation analysis. Source annotations and docking rows are never modified.
const GO_ANNOTATION_PATH = ORGANISM.go;
const GO_ASPECTS = [
  { code: "BP", name: "Biological process", column: "Gene Ontology (biological process)" },
  { code: "MF", name: "Molecular function", column: "Gene Ontology (molecular function)" },
  { code: "CC", name: "Cellular component", column: "Gene Ontology (cellular component)" },
];
const GO_TIERS = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
const GO_BACKGROUNDS = { filtered: "Parameter-filtered", total: "Total (target AA)" };
const goState = { aa: "ALA", top: 5, metric: "vina_sfct_combined_50", p2rank: 0.7, plddt: 90, maxCompetitors: 19, requireLPreference: false,
  background: "filtered", aspect: "BP", fdr: 0.05, significant: true, search: "", sort: "fdr", direction: 1, page: 1, proteinPage: 1, term: null };
const GO_PAGE_SIZE = 10;
let goAnnotationPromise = null;
let goResult = null;
let goRequest = 0;
let goInputTimer;

function parseGOAnnotations(text) {
  const lines = text.replace(/^\uFEFF/, "").split(/\r?\n/);
  const headers = (lines.shift() || "").split("\t").map((value) => value.trim());
  const required = ["Entry", ...GO_ASPECTS.map((aspect) => aspect.column)];
  const missing = required.filter((column) => !headers.includes(column));
  if (missing.length) throw new Error(`Invalid GO TSV: missing columns ${missing.join(", ")}.`);
  const entryIndex = headers.indexOf("Entry");
  const columns = GO_ASPECTS.map((aspect) => ({ ...aspect, index: headers.indexOf(aspect.column) }));
  const proteins = new Map(), terms = new Map();
  let sourceRows = 0;
  for (let index = 0; index < lines.length; index++) {
    if (!lines[index].trim()) continue;
    const cells = lines[index].split("\t");
    if (cells.length > headers.length) throw new Error(`Unexpected column count in GO TSV row ${index + 2}.`);
    const protein = (cells[entryIndex] || "").trim();
    if (!protein) continue;
    sourceRows++;
    if (!proteins.has(protein)) proteins.set(protein, { BP: new Set(), MF: new Set(), CC: new Set() });
    const annotations = proteins.get(protein);
    for (const aspect of columns) {
      for (const entry of (cells[aspect.index] || "").split(";")) {
        const match = entry.match(/\bGO:\d{7}\b/);
        if (!match) continue;
        const id = match[0];
        const name = entry.slice(0, match.index).replace(/[\s\[]+$/, "").trim() || id;
        const previous = terms.get(id);
        if (previous && previous.aspect !== aspect.code) throw new Error(`GO term ${id} occurs in conflicting GO aspects.`);
        if (!previous) terms.set(id, { id, name, aspect: aspect.code });
        annotations[aspect.code].add(id);
      }
    }
  }
  if (!proteins.size || !terms.size) throw new Error("The GO file contains no usable protein-to-GO mappings.");
  return { proteins, terms, sourceRows };
}

async function loadGOAnnotations() {
  if (!goAnnotationPromise) {
    goAnnotationPromise = (async () => {
      const response = await atlasFetch(GO_ANNOTATION_PATH, { cache: "no-cache" });
      if (!response.ok) throw new Error(`Could not load ${GO_ANNOTATION_PATH} (HTTP ${response.status}).`);
      const bytes = new Uint8Array(await response.arrayBuffer());
      let text;
      if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
        if (typeof DecompressionStream === "undefined") throw new Error("Please use an up-to-date browser to read the compressed GO annotations.");
        try {
          text = await new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"))).text();
        } catch { throw new Error(`Could not decompress ${GO_ANNOTATION_PATH}. Check that it is a valid gzip file.`); }
      } else text = new TextDecoder().decode(bytes);
      return parseGOAnnotations(text);
    })().catch((error) => { goAnnotationPromise = null; throw error; });
  }
  return goAnnotationPromise;
}

function goRequiresLPreference(options) {
  return Boolean(options.requireLPreference && hasDControl(options.aa));
}

function goProteinSets(options) {
  const qualityEligible = getRanking(options.aa, options.metric, options);
  const background = options.background === "total" ? getRanking(options.aa, options.metric, null) : qualityEligible;
  const tier = GO_TIERS.includes(options.top) ? options.top : 5;
  // Background choice never changes the selected hits or their QC-proteome ranks.
  const selected = qualityEligible.filter((row) => {
    if (row.proteome_percentile > tier) return false;
    if (goRequiresLPreference(options) && !getStereoControl(row.uniprot_id, options.aa, options.metric, null, options)?.lPreferred) return false;
    if (options.maxCompetitors === 19) return true;
    const count = getComparison(row.uniprot_id, options.aa, options.metric, null, options).nearCompetitors;
    return Number.isFinite(count) && count <= options.maxCompetitors;
  });
  return { background, selected };
}

function goLogFactorials(size) {
  const values = new Float64Array(size + 1);
  for (let i = 2; i <= size; i++) values[i] = values[i - 1] + Math.log(i);
  return values;
}

// Log of P(X >= k), X ~ Hypergeometric(N, K, n): the one-sided Fisher exact test.
function goHypergeometricLogP(N, K, n, k, factorials = goLogFactorials(N)) {
  const low = Math.max(0, n + K - N), high = Math.min(n, K);
  if (k <= low) return 0;
  if (k > high) return -Infinity;
  const choose = (a, b) => factorials[a] - factorials[b] - factorials[a - b];
  let logMass = choose(K, k) + choose(N - K, n - k) - choose(N, n);
  let logTail = logMass;
  for (let x = k; x < high; x++) {
    logMass += Math.log(K - x) + Math.log(n - x) - Math.log(x + 1) - Math.log(N - K - n + x + 1);
    const largest = Math.max(logTail, logMass);
    logTail = largest + Math.log1p(Math.exp(Math.min(logTail, logMass) - largest));
  }
  return Math.min(0, logTail);
}

function goAdjustFDR(rows) {
  const ordered = [...rows].sort((a, b) => a.logP - b.logP || a.id.localeCompare(b.id));
  let adjusted = 0;
  for (let index = ordered.length - 1; index >= 0; index--) {
    adjusted = Math.min(adjusted, ordered[index].logP + Math.log(ordered.length / (index + 1)));
    ordered[index].logFDR = adjusted;
  }
  return rows;
}

function computeGOEnrichment(data, options) {
  const { background, selected } = goProteinSets(options);
  const selectedIDs = new Set(selected.map((row) => row.uniprot_id));
  const coverage = Object.fromEntries(GO_ASPECTS.map(({ code }) => [code, { N: 0, n: 0 }]));
  const counts = new Map();
  let annotatedBackground = 0, annotatedSelected = 0;
  for (const row of background) {
    const annotations = data.proteins.get(row.uniprot_id);
    if (!annotations) continue;
    const inSelected = selectedIDs.has(row.uniprot_id);
    if (GO_ASPECTS.some(({ code }) => annotations[code].size)) {
      annotatedBackground++;
      if (inSelected) annotatedSelected++;
    }
    for (const { code } of GO_ASPECTS) {
      if (!annotations[code].size) continue;
      coverage[code].N++;
      if (inSelected) coverage[code].n++;
      for (const id of annotations[code]) {
        if (!counts.has(id)) counts.set(id, { ...data.terms.get(id), K: 0, k: 0, proteins: [] });
        const term = counts.get(id);
        term.K++;
        if (inSelected) { term.k++; term.proteins.push(row.uniprot_id); }
      }
    }
  }
  const factorials = goLogFactorials(background.length);
  const rows = [...counts.values()].map((term) => {
    const { N, n } = coverage[term.aspect];
    return { ...term, N, n, expected: n * term.K / N, fold: n ? (term.k / n) / (term.K / N) : NaN,
      logP: goHypergeometricLogP(N, term.K, n, term.k, factorials) };
  });
  // All background terms, including zero-hit terms and all three aspects, belong to one testing family.
  goAdjustFDR(rows);
  rows.sort((a, b) => a.logFDR - b.logFDR || a.logP - b.logP || a.id.localeCompare(b.id));
  return { options: { ...options }, data, background, selected, coverage, annotatedBackground, annotatedSelected,
    rows, terms: new Map(rows.map((row) => [row.id, row])), tested: rows.length };
}

function goProbability(logValue, fullPrecision = false) {
  if (!Number.isFinite(logValue)) return "—";
  if (logValue > Math.log(0.001)) return fullPrecision ? Math.exp(logValue).toPrecision(15) : Math.exp(logValue).toFixed(3);
  const power = logValue / Math.LN10;
  let exponent = Math.floor(power), mantissa = Number((10 ** (power - exponent)).toPrecision(fullPrecision ? 15 : 3));
  if (mantissa >= 10) { mantissa /= 10; exponent++; }
  return `${mantissa}e${exponent}`;
}

function goTermURL(id) { return `https://amigo.geneontology.org/amigo/term/${encodeURIComponent(id)}`; }
function goIsEnriched(row) { return row.fold > 1 && row.logFDR <= Math.log(goState.fdr); }

function filteredGOTerms() {
  const search = goState.search.trim().toLowerCase();
  const value = (row) => ({ fdr: row.logFDR, p: row.logP, fold: row.fold, hits: row.k, term: row.id }[goState.sort]);
  return (goResult?.rows || []).filter((row) => row.k > 0
    && (goState.aspect === "ALL" || row.aspect === goState.aspect)
    && (!goState.significant || goIsEnriched(row))
    && (!search || `${row.id} ${row.name}`.toLowerCase().includes(search)))
    .sort((a, b) => (typeof value(a) === "string" ? value(a).localeCompare(value(b)) : value(a) - value(b)) * goState.direction
      || a.logP - b.logP || a.id.localeCompare(b.id));
}

function goPagination(id, pageKey, total, noun) {
  const pages = Math.max(1, Math.ceil(total / GO_PAGE_SIZE));
  goState[pageKey] = Math.max(1, Math.min(pages, goState[pageKey]));
  const page = goState[pageKey];
  $(id).innerHTML = `<span>${total.toLocaleString()} ${noun} · ${GO_PAGE_SIZE} per page</span><div class="pagination-controls"><button type="button" data-go-page="${pageKey}" data-step="-1" ${page === 1 ? "disabled" : ""}>← Previous</button><span>Page ${page} of ${pages}</span><button type="button" data-go-page="${pageKey}" data-step="1" ${page === pages ? "disabled" : ""}>Next →</button></div>`;
  return (page - 1) * GO_PAGE_SIZE;
}

function goEnrichmentChart(rows) {
  const shown = rows.filter(goIsEnriched).sort((a, b) => a.logFDR - b.logFDR || a.logP - b.logP).slice(0, 12);
  if (!shown.length) return `<p class="analysis-empty">No enriched terms at FDR ≤ ${goState.fdr} match the current view. The table can still show annotation counts.</p>`;
  const labels = shown.map((row) => row.name.length > 43 ? `${row.name.slice(0, 40)}…` : row.name);
  // Reserve only the width these labels need; give the rest to the bars.
  const labelContext = document.createElement("canvas").getContext("2d");
  if (labelContext) labelContext.font = '11px "Segoe UI", sans-serif';
  const left = Math.ceil(Math.max(...labels.map((label) => labelContext ? labelContext.measureText(label).width : label.length * 6))) + 20;
  const plotBottom = 16 + 27 * shown.length, height = Math.max(126, plotBottom + 48), width = 850, plotWidth = width - left - 140;
  const limit = Math.ceil(Math.max(2, ...shown.map((row) => row.fold)));
  const x = (value) => left + value / limit * plotWidth;
  const colorMin = -Math.log(goState.fdr) / Math.LN10;
  const colorMax = Math.max(colorMin + 0.1, ...shown.map((row) => -row.logFDR / Math.LN10));
  const pale = [205, 189, 226], dark = [78, 39, 123];
  const color = (value) => {
    const fraction = Math.max(0, Math.min(1, (value - colorMin) / (colorMax - colorMin)));
    return `rgb(${pale.map((channel, index) => Math.round(channel + fraction * (dark[index] - channel))).join(", ")})`;
  };
  const bars = shown.map((row, index) => {
    const y = 16 + index * 27, label = labels[index];
    const significance = -row.logFDR / Math.LN10;
    return `<g><title>${escapeHTML(`${row.id} · ${row.name}: ${row.k}/${row.n} selected, ${row.K}/${row.N} background; fold ${fmt(row.fold, 2)}, FDR ${goProbability(row.logFDR)}, −log10(FDR) ${fmt(significance, 2)}`)}</title><text x="${left - 10}" y="${y + 13}" text-anchor="end">${escapeHTML(label)}</text><rect class="go-chart-bar" x="${left}" y="${y}" width="${x(row.fold) - left}" height="17" rx="3" fill="${color(significance)}"/></g>`;
  }).join("");
  const ticks = [0, .25, .5, .75, 1].map((fraction) => `<line class="go-chart-grid" x1="${x(limit * fraction)}" x2="${x(limit * fraction)}" y1="10" y2="${plotBottom}"/><text x="${x(limit * fraction)}" y="${plotBottom + 17}" text-anchor="middle">${fmt(limit * fraction, 1)}</text>`).join("");
  const legendX = left + plotWidth + 28, legendWidth = 14;
  const legendHeight = Math.min(180, Math.max(64, plotBottom - 40));
  const legendY = Math.max(30, (plotBottom - legendHeight) / 2);
  const legendTicks = [0, .5, 1].map((fraction) => `<text x="${legendX + legendWidth + 8}" y="${legendY + (1 - fraction) * legendHeight + 4}">${fmt(colorMin + fraction * (colorMax - colorMin), 2)}</text>`).join("");
  const legend = `<defs><linearGradient id="go-fdr-gradient" x1="0%" y1="100%" x2="0%" y2="0%" color-interpolation="sRGB"><stop offset="0%" stop-color="${color(colorMin)}"/><stop offset="100%" stop-color="${color(colorMax)}"/></linearGradient></defs><g class="go-chart-legend"><title>Darker purple indicates greater −log10(FDR), or a smaller FDR. Color scale starts at FDR ${goState.fdr}.</title><text x="${legendX - 4}" y="${legendY - 12}">−log10(FDR)</text><rect x="${legendX}" y="${legendY}" width="${legendWidth}" height="${legendHeight}" rx="3" fill="url(#go-fdr-gradient)"/>${legendTicks}</g>`;
  return `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="Up to 12 enriched GO terms ranked by adjusted significance; bar length shows fold enrichment and darker purple shows higher minus log10 FDR">${ticks}${bars}<line class="go-chart-axis" x1="${left}" x2="${left + plotWidth}" y1="${plotBottom}" y2="${plotBottom}"/><text x="${left + plotWidth / 2}" y="${plotBottom + 36}" text-anchor="middle">Fold enrichment</text>${legend}</svg>`;
}

function renderGOTerms() {
  const rows = filteredGOTerms(), start = goPagination("#go-pagination", "page", rows.length, "terms");
  const headings = [["term", "GO term"], [null, "Aspect"], ["hits", "Selected k/n"], [null, "Background K/N"], ["fold", "Fold enrichment"], ["p", "P-value"], ["fdr", "FDR (BH)"]];
  $("#go-term-head").innerHTML = `<tr>${headings.map(([key, label]) => `<th${key ? ` aria-sort="${goState.sort === key ? goState.direction === 1 ? "ascending" : "descending" : "none"}"` : ""}>${key ? `<button type="button" class="sort-button" data-go-sort="${key}">${label} <span>${goState.sort === key ? goState.direction === 1 ? "↑" : "↓" : "↕"}</span></button>` : label}</th>`).join("")}</tr>`;
  $("#go-term-body").innerHTML = rows.slice(start, start + GO_PAGE_SIZE).map((row) => `<tr class="${goState.term === row.id ? "go-active-term" : ""}"><td class="go-term-cell"><a href="${goTermURL(row.id)}" target="_blank" rel="noopener noreferrer">${escapeHTML(row.id)}</a><span>${escapeHTML(row.name)}</span></td><td><abbr title="${GO_ASPECTS.find((aspect) => aspect.code === row.aspect).name}">${row.aspect}</abbr></td><td><button type="button" class="go-hit-button" data-go-term="${row.id}" title="Show selected proteins annotated to ${row.id}">${row.k}/${row.n} ↘</button></td><td>${row.K}/${row.N}</td><td>${fmt(row.fold, 2)}×</td><td>${goProbability(row.logP)}</td><td class="${goIsEnriched(row) ? "go-significant" : ""}">${goProbability(row.logFDR)}</td></tr>`).join("") || '<tr><td colspan="7" class="analysis-empty">No GO terms match this view. Try another tier or relax the term filters.</td></tr>';
  $("#go-term-note").textContent = `${rows.length.toLocaleString()} terms with selected hits shown before pagination. FDR was calculated across all ${goResult.tested.toLocaleString()} background terms, including zero-hit terms, before these display filters. Click k/n to inspect matching proteins.`;
  $("#go-chart").innerHTML = goEnrichmentChart(rows);
  $("#go-plot-download-menu").hidden = !rows.some(goIsEnriched);
  if ($("#go-plot-download-menu").hidden) $("#go-plot-download-menu").open = false;
}

function buildGOPlotExport() {
  const source = $("#go-chart svg");
  if (!goResult || !source) return null;
  const [, , width, chartHeight] = source.getAttribute("viewBox").split(/\s+/).map(Number);
  if (!(width > 0 && chartHeight > 0)) return null;
  const height = chartHeight + 90, options = goResult.options;
  const aspect = GO_ASPECTS.find(({ code }) => code === goState.aspect)?.name || "All three aspects";
  const tier = `Top ${options.top}%`;
  const subtitle = `${aspect} · ${METRICS[options.metric].label} · ${tier} · FDR ≤ ${goState.fdr} · Background: ${GO_BACKGROUNDS[options.background]}`;
  const quality = `P2Rank ≥ ${options.p2rank} · Pocket mean pLDDT ≥ ${options.plddt} · Competitive AAs ≤ ${options.maxCompetitors}${goRequiresLPreference(options) ? " · L>D" : ""}`;
  const description = `${subtitle}. ${quality}. Bar length: fold enrichment. Darker purple: higher −log10(FDR). Term search: ${goState.search || "none"}. Source: ${GO_ANNOTATION_PATH}.`;
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="go-export-title go-export-description">
    <title id="go-export-title">Enriched GO terms · ${escapeHTML(options.aa)}</title>
    <desc id="go-export-description">${escapeHTML(description)}</desc>
    <rect width="100%" height="100%" fill="#fffdf9"/>
    <style>
      text { fill: #746d7e; font-family: "Segoe UI", sans-serif; font-size: 11px; }
      .go-export-title { fill: #241d2e; font-family: Georgia, serif; font-size: 20px; }
      .go-chart-grid { stroke: #e5dfd8; stroke-width: 1; }
      .go-chart-axis { stroke: #746d7e; stroke-width: 1; }
    </style>
    <text class="go-export-title" x="20" y="26">Enriched GO terms · ${escapeHTML(options.aa)}</text>
    <text x="20" y="45">${escapeHTML(subtitle)}</text>
    <text x="20" y="62">${escapeHTML(quality)}</text>
    <g transform="translate(0 76)">${source.innerHTML}</g>
  </svg>`;
  const filename = `go_${options.aa.toLowerCase()}_top${options.top}_${METRICS[options.metric].short}_${goState.aspect}_fdr${goState.fdr}_bg_${options.background}${goRequiresLPreference(options) ? "_ld" : ""}_enrichment`.replace(/[^a-zA-Z0-9_.-]+/g, "_");
  return { svg, width, height, filename };
}

async function downloadGOPlot(format) {
  if (!["png", "svg"].includes(format)) return;
  let svgURL;
  try {
    const plot = buildGOPlotExport();
    if (!plot) { showToast("No enriched GO plot is available to download."); return; }
    const svgBlob = new Blob([plot.svg], { type: "image/svg+xml;charset=utf-8" });
    if (format === "svg") { downloadBlob(`${plot.filename}.svg`, svgBlob); return; }
    svgURL = URL.createObjectURL(svgBlob);
    const image = new Image();
    await new Promise((resolve, reject) => {
      image.onload = resolve;
      image.onerror = reject;
      image.src = svgURL;
    });
    const canvas = document.createElement("canvas"), scale = 2;
    canvas.width = plot.width * scale;
    canvas.height = plot.height * scale;
    const context = canvas.getContext("2d");
    if (!context) throw new Error("Canvas is unavailable");
    context.scale(scale, scale);
    context.drawImage(image, 0, 0, plot.width, plot.height);
    const pngBlob = await new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
    if (!pngBlob) throw new Error("PNG export failed");
    downloadBlob(`${plot.filename}.png`, pngBlob);
  } catch (error) {
    console.error(error);
    showToast("Could not export the GO plot. Please try SVG instead.");
  } finally {
    if (svgURL) URL.revokeObjectURL(svgURL);
  }
}

function goDisplayedProteins() {
  if (!goResult) return [];
  const term = goResult.terms.get(goState.term);
  const ids = term ? new Set(term.proteins) : null;
  return goResult.selected.filter((row) => !ids || ids.has(row.uniprot_id));
}

function renderGOProteins() {
  const rows = goDisplayedProteins(), start = goPagination("#go-protein-pagination", "proteinPage", rows.length, "proteins");
  const options = goResult.options;
  const term = goResult.terms.get(goState.term);
  $("#go-protein-heading").textContent = term ? `${term.id} · Selected proteins` : "Selected proteins";
  $("#go-protein-note").textContent = term ? term.name : "One best pocket per protein, using the protein-selection settings above, including L>D when enabled. Background choice does not change these hits. Proteins without GO annotations remain visible but do not enter enrichment tests.";
  $("#go-clear-term").hidden = !term;
  $("#go-protein-body").innerHTML = rows.slice(start, start + GO_PAGE_SIZE).map((row) => {
    const annotation = goResult.data.proteins.get(row.uniprot_id);
    const count = annotation ? GO_ASPECTS.reduce((sum, { code }) => sum + annotation[code].size, 0) : 0;
    const comparison = getComparison(row.uniprot_id, options.aa, options.metric, null, options);
    const competitors = Number.isFinite(comparison.nearCompetitors) ? `${comparison.nearCompetitors} / ${comparison.otherCount}` : "—";
    const competitorNote = Number.isFinite(comparison.nearCompetitors)
      ? `Other AA best successful score ≤ ${options.aa} score + ${fmt(METRICS[options.metric].nearWindow, 2)} (${METRICS[options.metric].label}). Count / available other canonical AAs. Target AA and D-AA controls excluded. Like Explorer, other-AA scores pass the same QC thresholds before best-pocket selection; they are not tier-filtered.`
      : "Missing comparison: no successful target score or no successful scores for other canonical AAs. Missing results are not counted as noncompetitive.";
    return `<tr><td>${proteinTableIdentity(row.uniprot_id, row)}</td><td>${fmt(row[options.metric])}</td><td>${fmt(row.proteome_percentile, 2)}%</td><td class="numeric" title="${escapeHTML(competitorNote)}">${competitors}</td><td>${escapeHTML(row.pocket)}</td><td>${fmt(row.probability)}</td><td>${fmt(row.mean_pocket_plddt, 1)}</td><td>${count || "Missing"}</td><td><button type="button" class="open-row" data-go-protein="${escapeHTML(row.uniprot_id)}" aria-label="Open ${escapeHTML(row.uniprot_id)} protein profile" title="Open protein profile">→</button></td></tr>`;
  }).join("") || '<tr><td colspan="9" class="analysis-empty">No proteins meet this selection.</td></tr>';
}

function renderGOResults() {
  const result = goResult;
  const enriched = result.rows.filter(goIsEnriched).length;
  $("#go-summary").innerHTML = analysisCards([
    ["Selected proteins", result.selected.length.toLocaleString(), `Top ${result.options.top}% + QC · Competitive AAs ≤ ${result.options.maxCompetitors}${goRequiresLPreference(result.options) ? " · L>D" : ""}`],
    ["GO-annotated selected", result.annotatedSelected.toLocaleString(), `${result.selected.length - result.annotatedSelected} without GO annotations`],
    ["Background proteins", result.background.length.toLocaleString(), `${GO_BACKGROUNDS[result.options.background]} · ${result.annotatedBackground.toLocaleString()} with GO annotations`],
    ["Enriched terms", enriched.toLocaleString(), `All aspects · FDR ≤ ${goState.fdr}`],
  ]);
  $("#go-coverage-body").innerHTML = GO_ASPECTS.map(({ code, name }) => {
    const { n, N } = result.coverage[code];
    return `<tr><th>${name} (${code})</th><td>${n}/${result.selected.length}</td><td>${N}/${result.background.length}</td></tr>`;
  }).join("");
  const warning = !result.selected.length ? `No proteins pass these settings. Try a broader tier, lower quality thresholds or a higher Competitive AAs limit.${goRequiresLPreference(result.options) ? " You can also turn off L>D." : ""}`
    : !result.annotatedSelected ? "None of the selected proteins has GO annotations in this file. No enrichment can be assessed."
    : result.annotatedSelected === result.annotatedBackground ? "All GO-annotated background proteins are selected. There is no annotated comparison population outside the selection."
    : "";
  $("#go-warning").hidden = !warning;
  $("#go-warning").textContent = warning;
  $("#go-source-note").textContent = `Source: ${GO_ANNOTATION_PATH} · ${result.data.proteins.size.toLocaleString()} UniProt entries · ${result.data.terms.size.toLocaleString()} GO terms. Local annotations only; no external enrichment API.`;
  renderGOTerms();
  renderGOProteins();
}

async function renderGOAnalysis() {
  clearTimeout(goInputTimer);
  const request = ++goRequest, options = { ...goState };
  goResult = null;
  $("#go-results").hidden = true;
  $("#go-error").hidden = true;
  $("#go-status").hidden = false;
  $("#go-status").textContent = "Loading GO annotations and calculating enrichment…";
  $("#go-workspace").setAttribute("aria-busy", "true");
  try {
    const data = await loadGOAnnotations();
    await new Promise((resolve) => setTimeout(resolve, 0));
    if (request !== goRequest) return;
    goResult = computeGOEnrichment(data, options);
    goState.page = 1; goState.proteinPage = 1; goState.term = null;
    renderGOResults();
    $("#go-results").hidden = false;
  } catch (error) {
    if (request !== goRequest) return;
    goResult = null;
    $("#go-error-message").textContent = `Unavailable — ${error.message}`;
    $("#go-error").hidden = false;
  } finally {
    if (request === goRequest) {
      $("#go-status").hidden = true;
      $("#go-workspace").setAttribute("aria-busy", "false");
    }
  }
}

function goDelimited(headers, rows, delimiter) {
  const cell = (value) => {
    if (typeof value === "number") return Number.isFinite(value) ? String(value) : "";
    const text = String(value ?? ""), safe = /^[=+@-]/.test(text) ? `'${text}` : text;
    return safe.includes(delimiter) || /["\r\n]/.test(safe) ? `"${safe.replaceAll('"', '""')}"` : safe;
  };
  return [headers, ...rows].map((row) => row.map(cell).join(delimiter)).join("\n");
}

function downloadGOTerms(delimiter) {
  if (!goResult) return;
  const options = goResult.options;
  const headers = ["go_id", "go_name", "aspect", "selected_hits_k", "selected_annotated_n", "background_hits_K", "background_annotated_N", "expected_hits", "fold_enrichment", "p_value", "fdr_bh", "minus_log10_p", "minus_log10_fdr", "enriched_at_fdr_threshold", "fdr_threshold", "tested_terms_all_aspects", "target_aa", "proteome_tier_percent", "ranking_score", "p2rank_min", "plddt_min", "max_competitive_aas", "require_l_preference", "background_mode", "selected_uniprot_ids", "annotation_source"];
  const rows = filteredGOTerms().map((row) => [row.id, row.name, row.aspect, row.k, row.n, row.K, row.N, row.expected, row.fold, goProbability(row.logP, true), goProbability(row.logFDR, true), -row.logP / Math.LN10, -row.logFDR / Math.LN10, goIsEnriched(row), goState.fdr, goResult.tested, options.aa, options.top, options.metric, options.p2rank, options.plddt, options.maxCompetitors, goRequiresLPreference(options), options.background, row.proteins.join(";"), GO_ANNOTATION_PATH]);
  const extension = delimiter === "\t" ? "tsv" : "csv";
  downloadText(`go_${options.aa.toLowerCase()}_top${options.top}_${METRICS[options.metric].short}_${goState.aspect}_bg_${options.background}${goRequiresLPreference(options) ? "_ld" : ""}.${extension}`, goDelimited(headers, rows, delimiter), delimiter === "\t" ? "text/tab-separated-values" : "text/csv");
}

function downloadGOProteins(background = false) {
  if (!goResult) return;
  const options = goResult.options;
  const headers = ["uniprot_id", "gene_symbol", ORGANISM.rowIdentifier, "target_aa", "tier_percent", "ranking_score", "score", "pocket", "proteome_percentile", "competitive_aas", "available_other_aa_count", "competitive_aa_tolerance", "p2rank_probability", "pocket_mean_plddt", "p2rank_min", "plddt_min", "max_competitive_aas", "require_l_preference", "background_mode", "bp_terms", "mf_terms", "cc_terms", "inspected_go_term", "annotation_source"];
  const rows = (background ? goResult.background : goDisplayedProteins()).map((row) => {
    const annotation = goResult.data.proteins.get(row.uniprot_id);
    const comparison = getComparison(row.uniprot_id, options.aa, options.metric, null, options);
    return [row.uniprot_id, geneSymbol(state.annotations.get(row.uniprot_id)), row[ORGANISM.rowIdentifier], options.aa, options.top, options.metric, row[options.metric], row.pocket, row.proteome_percentile, comparison.nearCompetitors, comparison.otherCount, METRICS[options.metric].nearWindow, row.probability, row.mean_pocket_plddt, options.p2rank, options.plddt, options.maxCompetitors, goRequiresLPreference(options), options.background,
      ...GO_ASPECTS.map(({ code }) => [...(annotation?.[code] || [])].join(";")), background ? "" : goState.term || "", GO_ANNOTATION_PATH];
  });
  downloadText(`go_${options.aa.toLowerCase()}_${background ? "background" : goState.term?.replace(":", "_") || `top${options.top}_selected`}_bg_${options.background}${goRequiresLPreference(options) ? "_ld" : ""}_proteins.tsv`, goDelimited(headers, rows, "\t"));
}

function syncGOControls() {
  if (!GO_TIERS.includes(goState.top)) goState.top = 5;
  if (!Object.hasOwn(GO_BACKGROUNDS, goState.background)) goState.background = "filtered";
  for (const key of ["aa", "top", "metric", "p2rank", "plddt", "maxCompetitors", "background", "aspect", "fdr"]) $(`#go-${key}`).value = goState[key];
  $("#go-stereo-control").hidden = !hasDControl(goState.aa);
  $("#go-requireLPreference").checked = goState.requireLPreference;
  const stereoDescription = `Require predicted L-${goState.aa} preference over D-${goState.aa}; filters selected proteins, not the GO background. See Methods.`;
  $("#go-stereo-control").title = stereoDescription;
  $("#go-requireLPreference").setAttribute("aria-label", stereoDescription);
  $("#go-significant").checked = goState.significant;
  $("#go-significant-label").textContent = `FDR ≤ ${goState.fdr} only`;
  $("#go-chart-fdr").textContent = `FDR ≤ ${goState.fdr}`;
  $("#go-fdr").setCustomValidity("");
}

function bindGOEvents() {
  $("#go-aa").innerHTML = analysisOptions(AMINO_ACIDS.map(({ code, name }) => [code, `${code} · ${name}`]), goState.aa);
  $("#go-metric").innerHTML = analysisOptions(Object.entries(METRICS).map(([key, metric]) => [key, metric.label]), goState.metric);
  syncGOControls();
  const recalculate = () => {
    clearTimeout(goInputTimer); goRequest++; goResult = null;
    $("#go-results").hidden = true;
    $("#go-status").hidden = false;
    $("#go-status").textContent = "Updating GO analysis…";
    $("#go-workspace").setAttribute("aria-busy", "true");
    goInputTimer = setTimeout(renderGOAnalysis, 160);
  };
  for (const key of ["aa", "top", "metric", "p2rank", "plddt", "maxCompetitors", "background", "requireLPreference"]) {
    $(`#go-${key}`).addEventListener("change", (event) => {
      const value = event.target.value;
      goState[key] = key === "top" ? Number(value) : key === "p2rank" || key === "plddt" ? Math.min(key === "p2rank" ? 1 : 100, Math.max(0, Number(value) || 0)) : value;
      if (key === "maxCompetitors") goState[key] = Math.min(19, Math.max(0, Math.floor(Number(value) || 0)));
      if (key === "requireLPreference") goState[key] = event.target.checked;
      syncGOControls();
      recalculate();
    });
  }
  $("#go-use-explorer").addEventListener("click", () => {
    for (const key of ["aa", "metric", "p2rank", "plddt", "maxCompetitors", "requireLPreference"]) goState[key] = state[key];
    if (GO_TIERS.includes(state.top)) goState.top = state.top;
    syncGOControls(); recalculate();
  });
  $("#go-retry").addEventListener("click", renderGOAnalysis);
  $("#go-fdr").addEventListener("input", (event) => { event.target.setCustomValidity(""); });
  $("#go-fdr").addEventListener("change", (event) => {
    const value = Number(event.target.value);
    if (!Number.isFinite(value) || value <= 0 || value > 1) {
      event.target.setCustomValidity("Enter an FDR threshold greater than 0 and at most 1.");
      event.target.reportValidity();
      return;
    }
    goState.fdr = value; goState.page = 1;
    syncGOControls();
    // A display threshold must not change the protein sets or the BH testing family.
    if (goResult) renderGOResults();
  });
  $("#go-aspect").addEventListener("change", (event) => { goState.aspect = event.target.value; goState.page = 1; if (goResult) renderGOTerms(); });
  $("#go-significant").addEventListener("change", (event) => { goState.significant = event.target.checked; goState.page = 1; if (goResult) renderGOTerms(); });
  $("#go-term-search").addEventListener("input", (event) => { goState.search = event.target.value; goState.page = 1; if (goResult) renderGOTerms(); });
  $("#go-clear-term").addEventListener("click", () => { goState.term = null; goState.proteinPage = 1; if (goResult) { renderGOProteins(); renderGOTerms(); } });
  $("#go-workspace").addEventListener("click", (event) => {
    if (!goResult) return;
    const sort = event.target.closest("[data-go-sort]");
    if (sort) {
      const key = sort.dataset.goSort;
      goState.direction = goState.sort === key ? -goState.direction : ["hits", "fold"].includes(key) ? -1 : 1;
      goState.sort = key; goState.page = 1; renderGOTerms();
    }
    const page = event.target.closest("[data-go-page]");
    if (page) { goState[page.dataset.goPage] += Number(page.dataset.step); if (page.dataset.goPage === "page") renderGOTerms(); else renderGOProteins(); }
    const term = event.target.closest("[data-go-term]");
    if (term && goResult.terms.has(term.dataset.goTerm)) {
      goState.term = term.dataset.goTerm; goState.proteinPage = 1; renderGOProteins(); renderGOTerms();
      $("#go-proteins-panel").scrollIntoView({ behavior: "smooth", block: "nearest" });
    }
    const protein = event.target.closest("[data-go-protein]");
    if (protein) {
      state.aa = goResult.options.aa; state.metric = goResult.options.metric;
      state.p2rank = goResult.options.p2rank; state.plddt = goResult.options.plddt;
      state.profileValue = Object.keys(PROFILE_VALUES).find((key) => key !== "percentile" && PROFILE_VALUES[key].metric === state.metric);
      if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
      syncControls(); renderExplorer(); selectProtein(protein.dataset.goProtein, goResult.selected.find((row) => row.uniprot_id === protein.dataset.goProtein));
    }
  });
  for (const [format, delimiter] of [["tsv", "\t"], ["csv", ","]]) {
    $(`#go-download-${format}`).addEventListener("click", () => { downloadGOTerms(delimiter); $("#go-download-menu").open = false; });
  }
  $("#go-download-proteins").addEventListener("click", () => downloadGOProteins());
  $("#go-download-background").addEventListener("click", () => downloadGOProteins(true));
  for (const format of ["png", "svg"]) {
    $(`#go-download-plot-${format}`).addEventListener("click", () => {
      $("#go-plot-download-menu").open = false;
      downloadGOPlot(format);
    });
  }
}
