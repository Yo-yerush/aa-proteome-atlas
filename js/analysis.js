// Independent proteome analyses. Source rows remain untouched.
const analysisState = {
  compare: { x: "ALA", y: "LEU", metric: "vina_sfct_combined_50", p2: 0.7, plddt: 90, preference: "all", minDelta: "", search: "", page: 1, sort: "delta", direction: -1 },
  overlap: { aas: ["ALA", "LEU"], metric: "vina_sfct_combined_50", p2: 0.7, plddt: 90, tier: 5, group: "shared", search: "", page: 1 },
  statistics: { aa: "ALA", metric: "vina_sfct_combined_50", p2: 0.7, plddt: 90, unit: "best", scatter: "vina_sfct" },
};
const analysisResults = {};
const ANALYSIS_PAGE_SIZE = 10;

function qualityPasses(row, options) {
  return isSuccessfulResult(row) && passesPocketQuality(row, options);
}

function qualityRanking(aa, options) {
  return rankedPopulation(aa, options.metric, options);
}

function analysisCutoffs(ranking, metric) {
  return [1, 5, 10].map((tier) => {
    const count = ranking.rows.length ? Math.floor(tier / 100 * (ranking.rows.length - 1)) + 1 : 0;
    return { tier, count, value: count ? ranking.rows[count - 1][metric] : NaN };
  });
}

function pairedAAResults(options) {
  const x = qualityRanking(options.x, options);
  const y = qualityRanking(options.y, options);
  const rows = x.rows.filter((row) => y.byProtein.has(row.uniprot_id)).map((row) => {
    const other = y.byProtein.get(row.uniprot_id);
    return { protein: row.uniprot_id, x: row[options.metric], y: other[options.metric], delta: other[options.metric] - row[options.metric],
      pocketX: row.pocket, pocketY: other.pocket, percentileX: x.percentile.get(row.uniprot_id), percentileY: y.percentile.get(row.uniprot_id) };
  });
  return { rows, missingX: y.rows.length - rows.length, missingY: x.rows.length - rows.length };
}

function topHitOverlap(options) {
  const rankings = options.aas.map((aa) => qualityRanking(aa, options));
  const sets = rankings.map((ranking) => new Set(ranking.rows.filter((row) => ranking.percentile.get(row.uniprot_id) <= options.tier).map((row) => row.uniprot_id)));
  const union = new Set(sets.flatMap((set) => [...set]));
  // Count breadth across all canonical AAs, independently of the selected overlap sets.
  const allRankings = union.size ? AMINO_ACIDS.map(({ code }) => ({ code, ranking: qualityRanking(code, options) })) : [];
  const groups = new Map();
  const rows = [...union].sort().map((protein) => {
    const membership = sets.map((set) => set.has(protein));
    const key = membership.map((member) => member ? "1" : "0").join("");
    const statuses = rankings.map((ranking, index) => !ranking.byProtein.has(protein) ? "Missing / fails QC" : membership[index] ? `Top ${options.tier}%` : "Outside tier");
    const topHitAAs = allRankings.filter(({ ranking }) => {
      const percentile = ranking.percentile.get(protein);
      return Number.isFinite(percentile) && percentile <= options.tier;
    }).map(({ code }) => code);
    const topHitAACount = topHitAAs.length;
    const availableAACount = allRankings.filter(({ ranking }) => ranking.byProtein.has(protein)).length;
    const availableSelectedAACount = rankings.filter((ranking) => ranking.byProtein.has(protein)).length;
    groups.set(key, (groups.get(key) || 0) + 1);
    return { protein, membership, key, statuses, percentiles: rankings.map((ranking) => ranking.percentile.get(protein)),
      topHitAAs, topHitAACount, availableAACount, availableSelectedAACount };
  });
  return { rows, sets, rankings, groups: [...groups].map(([key, count]) => ({ key, count })).sort((a, b) => b.count - a.count || a.key.localeCompare(b.key)) };
}

function pairStatistics(pairs, retainPairs = true) {
  const valid = pairs.filter((pair) => Number.isFinite(pair.x) && Number.isFinite(pair.y));
  const xs = valid.map((pair) => pair.x), ys = valid.map((pair) => pair.y);
  return { ...(retainPairs ? { pairs: valid } : {}), n: valid.length, r: pearson(xs, ys), rho: pearson(rankValues(xs), rankValues(ys)) };
}

function numericSummary(values) {
  const valid = values.filter(Number.isFinite).sort((a, b) => a - b);
  const atQuantile = (q) => {
    if (!valid.length) return NaN;
    const position = (valid.length - 1) * q, base = Math.floor(position);
    return valid[base] + ((valid[base + 1] ?? valid[base]) - valid[base]) * (position - base);
  };
  return { n: valid.length, mean: mean(valid), sd: valid.length ? standardDeviation(valid) : NaN, median: atQuantile(.5), q1: atQuantile(.25), q3: atQuantile(.75) };
}

function macroStatistics(options) {
  const aas = options.aa === "ALL" ? AMINO_ACIDS.map((aa) => aa.code) : [options.aa];
  const records = [];
  const summary = [];
  const retainedPockets = new Set();
  const proteins = new Set();
  const pocketsPerPair = [];
  for (const aa of aas) {
    const ranking = qualityRanking(aa, options);
    const eligible = (state.rawByAA.get(aa) || []).filter((row) => qualityPasses(row, options) && hasUsableScore(row, options.metric));
    const counts = new Map();
    for (const row of eligible) {
      retainedPockets.add(`${row.uniprot_id}|${row.pocket}`);
      proteins.add(row.uniprot_id);
      if (!counts.has(row.uniprot_id)) counts.set(row.uniprot_id, new Set());
      counts.get(row.uniprot_id).add(row.pocket);
    }
    pocketsPerPair.push(...[...counts.values()].map((set) => set.size));
    const chosen = options.unit === "all" ? eligible : ranking.rows;
    for (const row of chosen) records.push({ aa, row });
    summary.push({ aa, proteins: ranking.rows.length, ...numericSummary(chosen.map((row) => row[options.metric])), cutoffs: analysisCutoffs(ranking, options.metric) });
  }
  // Selectivity always compares independent best quality-passing pockets, never missing scores as zeros.
  const profiles = new Map();
  for (const { code } of AMINO_ACIDS) {
    for (const row of qualityRanking(code, options).rows) {
      if (!profiles.has(row.uniprot_id)) profiles.set(row.uniprot_id, []);
      profiles.get(row.uniprot_id).push({ aa: code, value: row[options.metric] });
    }
  }
  for (const profile of profiles.values()) profile.sort((a, b) => a.value - b.value);
  const deltas = [], percentiles = [];
  for (const aa of aas) {
    const ranking = qualityRanking(aa, options);
    percentiles.push(...ranking.percentile.values());
    if (!AMINO_ACIDS.some((entry) => entry.code === aa)) continue;
    for (const row of ranking.rows) {
      const profile = profiles.get(row.uniprot_id) || [];
      const otherCount = profile.length - 1;
      if (otherCount < 1) continue;
      const omitted = profile.findIndex((entry) => entry.aa === aa);
      const valueAt = (index) => profile[index >= omitted ? index + 1 : index].value;
      const middle = Math.floor(otherCount / 2);
      const otherMedian = otherCount % 2 ? valueAt(middle) : (valueAt(middle - 1) + valueAt(middle)) / 2;
      deltas.push(otherMedian - row[options.metric]);
    }
  }
  const specs = [
    ["vina_sfct", "Vina", "SFCT", "vina_affinity", "sfct_score"],
    ["vina_combined", "Vina", "Combined 80%", "vina_affinity", "vina_sfct_combined"],
    ["vina_combined_50", "Vina", "Combined 50%", "vina_affinity", "vina_sfct_combined_50"],
    ["combined_weights", "Combined 80%", "Combined 50%", "vina_sfct_combined", "vina_sfct_combined_50"],
    ["p2_score", "P2Rank probability", METRICS[options.metric].label, "probability", options.metric],
    ["plddt_score", "Pocket mean pLDDT", METRICS[options.metric].label, "mean_pocket_plddt", options.metric],
  ];
  // Keep the source records once; materialize all plotted pairs only for the selected variables.
  const correlations = specs.map(([id, xLabel, yLabel, xKey, yKey]) => ({ id, xLabel, yLabel, xKey, yKey,
    ...pairStatistics(records.map(({ row }) => ({ x: row[xKey], y: row[yKey] })), false) }));
  let stereo = [];
  if (hasDControl(options.aa)) {
    const aa = canonicalAACode(options.aa);
    stereo = pairedAAResults({ ...options, x: aa, y: dControlCode(aa) }).rows;
    correlations.push({ id: "l_d", xLabel: `L-${aa}`, yLabel: `D-${aa}`, ...pairStatistics(stereo, false) });
  }
  return { records, summary, proteins: proteins.size, pockets: retainedPockets.size, pocketsPerPair, deltas, percentiles, correlations, stereo };
}

function analysisOptions(items, selected) {
  return items.map(([value, label]) => `<option value="${escapeHTML(value)}"${String(value) === String(selected) ? " selected" : ""}>${escapeHTML(label)}</option>`).join("");
}

function analysisSelect(id, label, options, selected) {
  return `<label class="field"><span>${escapeHTML(label)}</span><select id="${id}">${analysisOptions(options, selected)}</select></label>`;
}

function analysisNumber(id, label, value, min, max, step = 1, placeholder = "") {
  return `<label class="field"><span>${escapeHTML(label)}</span><input id="${id}" type="number" value="${value}" ${min !== null ? `min="${min}"` : ""} ${max !== null ? `max="${max}"` : ""} step="${step}" placeholder="${placeholder}" /></label>`;
}

function commonAnalysisControls(view) {
  const options = analysisState[view];
  return analysisSelect(`${view}-metric`, "Ranking score", Object.entries(METRICS).map(([key, metric]) => [key, metric.label]), options.metric)
    + analysisNumber(`${view}-p2`, "P2Rank probability ≥", options.p2, 0, 1, .05)
    + analysisNumber(`${view}-plddt`, "Pocket mean pLDDT ≥", options.plddt, 0, 100, 5);
}

function analysisCards(items) {
  return items.map(([label, value, note]) => `<article><span>${escapeHTML(label)}</span><strong>${escapeHTML(value)}</strong><small>${escapeHTML(note || "")}</small></article>`).join("");
}

function analysisProteinCell(protein) {
  return proteinTableIdentity(protein);
}

function analysisProteinProfileButton(protein) {
  return `<button class="open-row" type="button" data-analysis-protein="${escapeHTML(protein)}" aria-label="Open ${escapeHTML(protein)} protein profile" title="Open protein profile">→</button>`;
}

function analysisMatches(protein, search) {
  return `${protein} ${annotationSearchText(protein)}`.toLowerCase().includes(search.trim().toLowerCase());
}

function analysisPagination(view, count) {
  const options = analysisState[view];
  const pages = Math.max(1, Math.ceil(count / ANALYSIS_PAGE_SIZE));
  options.page = Math.min(pages, Math.max(1, options.page));
  $(`#${view}-pagination`).innerHTML = `<span>${count.toLocaleString()} proteins · ${ANALYSIS_PAGE_SIZE} per page</span><div class="pagination-controls"><button type="button" data-analysis-page="${view}" data-step="-1" ${options.page === 1 ? "disabled" : ""}>← Previous</button><span>Page ${options.page} of ${pages}</span><button type="button" data-analysis-page="${view}" data-step="1" ${options.page === pages ? "disabled" : ""}>Next →</button></div>`;
  return (options.page - 1) * ANALYSIS_PAGE_SIZE;
}

function analysisTSV(headers, rows) {
  const cell = (value) => {
    if (value === null || value === undefined || (typeof value === "number" && !Number.isFinite(value))) return "";
    const text = String(value);
    // Spreadsheet-formula protection for identifiers/annotations, preserving numeric negatives.
    const safe = typeof value === "string" && /^[=+@-]/.test(text) ? `'${text}` : text;
    return /[\t\r\n"]/.test(safe) ? `"${safe.replaceAll('"', '""')}"` : safe;
  };
  return [headers, ...rows].map((row) => row.map(cell).join("\t")).join("\n");
}

function analysisHistogram(values, label, { integer = false, zero = false } = {}) {
  const valid = values.filter(Number.isFinite);
  if (!valid.length) return `<p class="analysis-empty">No eligible observations for this distribution.</p>`;
  let low = Infinity, high = -Infinity;
  for (const value of valid) { low = Math.min(low, value); high = Math.max(high, value); }
  if (integer) { low -= .5; high += .5; }
  else if (low === high) { low -= .5; high += .5; }
  const count = integer ? Math.min(60, Math.max(1, Math.round(high - low))) : 40;
  const bins = Array(count).fill(0), step = (high - low) / count;
  for (const value of valid) bins[Math.min(count - 1, Math.max(0, Math.floor((value - low) / step)))]++;
  const peak = Math.max(1, ...bins), x = (value) => 62 + (value - low) / (high - low) * 506;
  const y = (value) => 254 - value / peak * 210;
  const ticks = [0, .5, 1].map((fraction) => `<line class="analysis-grid-line" x1="62" x2="568" y1="${y(peak * fraction)}" y2="${y(peak * fraction)}"/><text x="54" y="${y(peak * fraction) + 4}" text-anchor="end">${Math.round(peak * fraction).toLocaleString()}</text>`).join("");
  const bars = bins.map((n, index) => `<rect class="analysis-bar" x="${x(low + index * step) + .5}" y="${y(n)}" width="${Math.max(1, 506 / count - 1)}" height="${254 - y(n)}"><title>${fmt(low + index * step, 2)} to ${fmt(low + (index + 1) * step, 2)}: ${n.toLocaleString()} observations</title></rect>`).join("");
  const axis = [0, .25, .5, .75, 1].map((fraction) => `<text x="${62 + fraction * 506}" y="275" text-anchor="middle">${fmt(low + fraction * (high - low), integer ? 1 : 2)}</text>`).join("");
  return `<svg viewBox="0 0 600 320" role="img" aria-label="${escapeHTML(label)} histogram, ${valid.length} observations">${ticks}${bars}${zero && low <= 0 && high >= 0 ? `<line class="analysis-reference" x1="${x(0)}" x2="${x(0)}" y1="35" y2="254"><title>Zero difference</title></line>` : ""}${axis}<text x="315" y="305" text-anchor="middle">${escapeHTML(label)}</text><text transform="translate(16 145) rotate(-90)" text-anchor="middle">Observations</text></svg>`;
}

function sampleAnalysisPoints(pairs) {
  const limit = Math.min(2000, pairs.length);
  return Array.from({ length: limit }, (_, i) => pairs[Math.floor(i * pairs.length / limit)]);
}

function analysisScatterGeometry(valid, identity) {
  let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
  for (const pair of valid) { xMin = Math.min(xMin, pair.x); xMax = Math.max(xMax, pair.x); yMin = Math.min(yMin, pair.y); yMax = Math.max(yMax, pair.y); }
  if (identity) { xMin = yMin = Math.min(xMin, yMin); xMax = yMax = Math.max(xMax, yMax); }
  const xPad = (xMax - xMin || 1) * .05, yPad = (yMax - yMin || 1) * .05;
  xMin -= xPad; xMax += xPad; yMin -= yPad; yMax += yPad;
  const x = (value) => 74 + (value - xMin) / (xMax - xMin) * 482;
  const y = (value) => 294 - (value - yMin) / (yMax - yMin) * 254;
  return { x, y, xMin, xMax, yMin, yMax };
}

function analysisScatterFrame(geometry, xLabel, yLabel, identity, marks, shown, total) {
  const { x, y, xMin, xMax, yMin, yMax } = geometry;
  const grid = [0, .25, .5, .75, 1].map((f) => `<line class="analysis-grid-line" x1="74" x2="556" y1="${294 - f * 254}" y2="${294 - f * 254}"/><text x="${74 + f * 482}" y="314" text-anchor="middle">${fmt(xMin + f * (xMax - xMin), 2)}</text><text x="66" y="${298 - f * 254}" text-anchor="end">${fmt(yMin + f * (yMax - yMin), 2)}</text>`).join("");
  return `<svg viewBox="0 0 600 360" role="img" aria-label="${escapeHTML(xLabel)} versus ${escapeHTML(yLabel)}; ${shown.toLocaleString()} of ${total.toLocaleString()} points shown">${grid}${identity ? `<line class="analysis-reference" x1="${x(xMin)}" y1="${y(xMin)}" x2="${x(xMax)}" y2="${y(xMax)}"><title>Equal scores</title></line>` : ""}${marks}<text x="315" y="347" text-anchor="middle">${escapeHTML(xLabel)}</text><text transform="translate(16 165) rotate(-90)" text-anchor="middle">${escapeHTML(yLabel)}</text><text x="556" y="20" text-anchor="end">${shown.toLocaleString()} of ${total.toLocaleString()} points shown</text></svg>`;
}

function analysisScatter(pairs, xLabel, yLabel, identity = false, total = pairs.length) {
  const valid = pairs.filter((pair) => Number.isFinite(pair.x) && Number.isFinite(pair.y));
  if (!valid.length) return `<p class="analysis-empty">No paired observations pass these filters.</p>`;
  const geometry = analysisScatterGeometry(valid, identity);
  const { x, y } = geometry;
  // Keep interaction responsive; every point still contributes to correlations and downloads.
  const shown = sampleAnalysisPoints(valid);
  const marks = shown.map((pair) => `<circle class="analysis-dot" cx="${x(pair.x)}" cy="${y(pair.y)}" r="2.4"><title>${escapeHTML(pair.label || [pair.protein, pair.aa, pair.pocket].filter(Boolean).join(" · ") || "Observation")} · ${escapeHTML(xLabel)} ${fmt(pair.x)} · ${escapeHTML(yLabel)} ${fmt(pair.y)}</title></circle>`).join("");
  return analysisScatterFrame(geometry, xLabel, yLabel, identity, marks, shown.length, total);
}

function renderCompare() {
  const options = analysisState.compare;
  const result = pairedAAResults(options);
  analysisResults.compare = result;
  const positive = result.rows.filter((row) => row.delta > 0).length;
  const negative = result.rows.filter((row) => row.delta < 0).length;
  $("#compare-summary").innerHTML = analysisCards([
    ["Paired proteins", result.rows.length.toLocaleString(), "Both AAs pass QC"],
    [`${options.x} scores better`, positive.toLocaleString(), "Δ > 0"],
    [`${options.y} scores better`, negative.toLocaleString(), "Δ < 0"],
    ["Median Δ", fmt(median(result.rows.map((row) => row.delta))), `${result.rows.length - positive - negative} tied scores`],
  ]);
  const correlations = pairStatistics(result.rows);
  $("#compare-correlations").textContent = `Pearson r ${fmt(correlations.r)} · Spearman ρ ${fmt(correlations.rho)} · n = ${correlations.n.toLocaleString()}`;
  $("#compare-scatter").innerHTML = analysisScatter(result.rows, `${options.x} ${METRICS[options.metric].short}`, `${options.y} ${METRICS[options.metric].short}`, true);
  $("#compare-list-note").textContent = `${result.missingX.toLocaleString()} proteins have only ${options.y} eligible; ${result.missingY.toLocaleString()} have only ${options.x} eligible. These unpaired proteins are excluded from the comparison, not from the atlas. Table filters below do not change the paired-score overview.`;
  renderCompareTable();
}

function filteredComparison() {
  const options = analysisState.compare;
  return (analysisResults.compare?.rows || []).filter((row) => analysisMatches(row.protein, options.search)
    && (options.preference === "all" || (options.preference === "x" ? row.delta > 0 : row.delta < 0))
    && (options.minDelta === "" || row.delta >= Number(options.minDelta)))
    .sort((a, b) => options.direction * (options.sort === "protein" ? a.protein.localeCompare(b.protein) : a[options.sort] - b[options.sort]) || a.protein.localeCompare(b.protein));
}

function renderCompareTable() {
  const options = analysisState.compare, rows = filteredComparison();
  const start = analysisPagination("compare", rows.length);
  const heading = (key, label) => `<th aria-sort="${options.sort === key ? options.direction === 1 ? "ascending" : "descending" : "none"}"><button type="button" class="sort-button" data-compare-sort="${key}">${escapeHTML(label)} ${options.sort === key ? options.direction === 1 ? "↑" : "↓" : "↕"}</button></th>`;
  $("#compare-head").innerHTML = `<tr>${heading("protein", "Protein")}${heading("x", `${options.x} score`)}${heading("y", `${options.y} score`)}${heading("delta", "Δ (Y − X)")}<th>${options.x} pocket</th><th>${options.y} pocket</th>${heading("percentileX", `${options.x} top %`)}${heading("percentileY", `${options.y} top %`)}<th><span class="sr-only">Protein profile</span></th></tr>`;
  $("#compare-body").innerHTML = rows.slice(start, start + ANALYSIS_PAGE_SIZE).map((row) => `<tr><td>${analysisProteinCell(row.protein)}</td><td>${fmt(row.x)}</td><td>${fmt(row.y)}</td><td class="${row.delta > 0 ? "positive" : row.delta < 0 ? "negative" : ""}">${fmt(row.delta)}</td><td>${escapeHTML(row.pocketX)}</td><td>${escapeHTML(row.pocketY)}</td><td>${fmt(row.percentileX, 2)}%</td><td>${fmt(row.percentileY, 2)}%</td><td>${analysisProteinProfileButton(row.protein)}</td></tr>`).join("") || `<tr><td colspan="9" class="analysis-empty">No paired proteins match. Try clearing the table filters or lowering QC thresholds.</td></tr>`;
}

function overlapGroupLabel(key, aas) {
  const included = aas.filter((_, index) => key[index] === "1");
  return included.length === aas.length ? "Shared by all selected AAs" : `${included.join(" + ")} only`;
}

function overlapPlot(result, aas) {
  if (!result.rows.length) return `<p class="analysis-empty">No top hits pass these filters.</p>`;
  if (aas.length === 2) {
    const shared = result.groups.find((group) => group.key === "11")?.count || 0;
    const left = result.sets[0].size - shared, right = result.sets[1].size - shared;
    return `<svg viewBox="0 0 600 310" role="img" aria-label="Top-hit Venn diagram: ${aas[0]} only ${left}, shared ${shared}, ${aas[1]} only ${right}"><circle cx="235" cy="154" r="105" fill="#7544a9" fill-opacity=".2" stroke="#7544a9"/><circle cx="365" cy="154" r="105" fill="#29869f" fill-opacity=".2" stroke="#29869f"/><text x="201" y="140" text-anchor="middle">${aas[0]} only</text><text class="overlap-count" x="201" y="172" text-anchor="middle">${left.toLocaleString()}</text><text x="300" y="140" text-anchor="middle">Shared</text><text class="overlap-count" x="300" y="172" text-anchor="middle">${shared.toLocaleString()}</text><text x="399" y="140" text-anchor="middle">${aas[1]} only</text><text class="overlap-count" x="399" y="172" text-anchor="middle">${right.toLocaleString()}</text></svg>`;
  }
  const groups = result.groups.slice(0, 20), width = Math.max(580, 105 + groups.length * 43), height = 230 + aas.length * 24;
  const maximum = Math.max(1, ...groups.map((group) => group.count));
  return `<svg style="min-width:${width}px" viewBox="0 0 ${width} ${height}" role="img" aria-label="UpSet plot of exact top-hit intersections"><text x="90" y="18">Proteins per exact intersection</text>${aas.map((aa, i) => `<text x="64" y="${220 + i * 24}" text-anchor="end">${aa}</text>`).join("")}${groups.map((group, index) => {
    const x = 103 + index * 43, h = group.count / maximum * 135;
    const active = aas.map((_, i) => i).filter((i) => group.key[i] === "1");
    return `<g><title>${escapeHTML(overlapGroupLabel(group.key, aas))}: ${group.count} proteins</title><rect class="analysis-bar" x="${x - 12}" y="${175 - h}" width="24" height="${h}"/><text x="${x}" y="${168 - h}" text-anchor="middle">${group.count}</text><line x1="${x}" x2="${x}" y1="${216 + active[0] * 24}" y2="${216 + active[active.length - 1] * 24}" stroke="#7544a9" stroke-width="2"/>${aas.map((_, i) => `<circle cx="${x}" cy="${216 + i * 24}" r="5" fill="${group.key[i] === "1" ? "#7544a9" : "#e7e1eb"}"/>`).join("")}</g>`;
  }).join("")}</svg>`;
}

function renderOverlap() {
  const options = analysisState.overlap, result = topHitOverlap(options);
  analysisResults.overlap = result;
  const shared = result.rows.filter((row) => row.membership.every(Boolean)).length;
  $("#overlap-summary").innerHTML = analysisCards([
    ["Top-hit union", result.rows.length.toLocaleString(), `Top ${options.tier}% in at least one selected AA`],
    ["Shared by all", shared.toLocaleString(), `${options.aas.length} selected AAs`],
    ["In one list only", result.rows.filter((row) => row.membership.filter(Boolean).length === 1).length.toLocaleString(), "Not proof of exclusive binding"],
    ["Incomplete coverage", result.rows.filter((row) => row.statuses.includes("Missing / fails QC")).length.toLocaleString(), "Top hits missing ≥1 eligible selected AA"],
  ]);
  $("#overlap-plot-title").textContent = options.aas.length === 2 ? "Two-AA overlap" : "Exact intersections (UpSet)";
  $("#overlap-plot").innerHTML = overlapPlot(result, options.aas);
  $("#overlap-plot-note").textContent = options.aas.length === 2 ? "Schematic circles; counts are exact, areas are not proportional." : `Showing the ${Math.min(20, result.groups.length)} largest of ${result.groups.length} exact intersections. Select any intersection in the table, including those not plotted.`;
  $("#overlap-matrix").innerHTML = `<table class="data-table"><thead><tr><th>Jaccard</th>${options.aas.map((aa) => `<th>${aa}</th>`).join("")}</tr></thead><tbody>${options.aas.map((aa, i) => `<tr><th>${aa}<small class="analysis-gene">${result.sets[i].size} hits / ${result.rankings[i].rows.length} eligible</small></th>${options.aas.map((other, j) => {
    const sharedCount = [...result.sets[i]].filter((protein) => result.sets[j].has(protein)).length;
    const unionCount = result.sets[i].size + result.sets[j].size - sharedCount;
    const value = unionCount ? sharedCount / unionCount : NaN;
    return `<td style="background:${Number.isFinite(value) ? `rgba(117,68,169,${.05 + value * .3})` : "#eee"}" title="${aa} × ${other}: ${sharedCount} shared / ${unionCount} union">${Number.isFinite(value) ? `${fmt(value * 100, 1)}%` : "—"}</td>`;
  }).join("")}</tr>`).join("")}</tbody></table>`;
  const groups = [["all", "All top-hit proteins"], ["shared", "Shared by all selected AAs"], ["unique", "In exactly one top-hit list"], ...result.groups.map(({ key, count }) => [key, `${overlapGroupLabel(key, options.aas)} (${count})`])];
  if (!groups.some(([key]) => key === options.group)) options.group = "shared";
  $("#overlap-group").innerHTML = analysisOptions(groups, options.group);
  renderOverlapTable();
}

function filteredOverlap() {
  const options = analysisState.overlap;
  return (analysisResults.overlap?.rows || []).filter((row) => analysisMatches(row.protein, options.search)
    && (options.group === "all" || (options.group === "shared" ? row.membership.every(Boolean) : options.group === "unique" ? row.membership.filter(Boolean).length === 1 : row.key === options.group)));
}

function renderOverlapTable() {
  const options = analysisState.overlap, rows = filteredOverlap(), start = analysisPagination("overlap", rows.length);
  $("#overlap-list-note").textContent = `Competitive AAs counts top-${options.tier}% membership across all ${AMINO_ACIDS.length} canonical AAs using the same score and P2Rank/pLDDT filters, not just the selected overlap AAs. Hover over a count for AA names and eligible coverage. Missing / fails QC remains unknown. Downloads include the entire matching list.`;
  $("#overlap-head").innerHTML = `<tr><th>Protein</th><th>Top-hit membership</th><th>Competitive AAs</th>${options.aas.map((aa) => `<th>${aa}</th>`).join("")}<th><span class="sr-only">Protein profile</span></th></tr>`;
  $("#overlap-body").innerHTML = rows.slice(start, start + ANALYSIS_PAGE_SIZE).map((row) => `<tr><td>${analysisProteinCell(row.protein)}</td><td>${escapeHTML(options.aas.filter((_, i) => row.membership[i]).join(" + "))}</td><td title="Top ${options.tier}% AAs: ${escapeHTML(row.topHitAAs.join(", "))}. Eligible results: ${row.availableAACount}/${AMINO_ACIDS.length}; missing / fails QC is unknown.">${row.topHitAACount}/${AMINO_ACIDS.length}${row.availableAACount < AMINO_ACIDS.length ? `<small class="analysis-gene">${row.availableAACount}/${AMINO_ACIDS.length} eligible</small>` : ""}</td>${row.statuses.map((status, i) => `<td><span class="${row.membership[i] ? "positive" : ""}">${escapeHTML(status)}</span>${Number.isFinite(row.percentiles[i]) ? `<small class="analysis-gene">${fmt(row.percentiles[i], 2)}%</small>` : ""}</td>`).join("")}<td>${analysisProteinProfileButton(row.protein)}</td></tr>`).join("") || `<tr><td colspan="${options.aas.length + 4}" class="analysis-empty">No proteins in this intersection match your search.</td></tr>`;
}

function renderStatistics() {
  const options = analysisState.statistics;
  const key = [options.aa, options.metric, options.p2, options.plddt, options.unit].join("|");
  const result = analysisResults.statisticsKey === key ? analysisResults.statistics : macroStatistics(options);
  analysisResults.statistics = result;
  analysisResults.statisticsKey = key;
  const pooled = result.pooled || (result.pooled = numericSummary(result.records.map(({ row }) => row[options.metric])));
  const allAA = options.aa === "ALL";
  $("#statistics-population-note").textContent = `Successful finite ${METRICS[options.metric].short} scores, after P2Rank ≥ ${options.p2} and pLDDT ≥ ${options.plddt}. ${options.unit === "best" ? "One best pocket per protein × AA." : "All passing pocket rows; proteins with more pockets contribute more observations."} ${allAA ? "All 20 canonical AAs are pooled (D-AA controls excluded); the same protein contributes separately for each available AA. These observations are not independent proteins." : "The selected AA is analyzed across the proteome."} Retained-pocket counts are unique protein × pocket IDs, not duplicated across AAs.`;
  $("#statistics-summary").innerHTML = analysisCards([
    ["Proteins", result.proteins.toLocaleString(), "At least one eligible result"],
    ["Retained pockets", result.pockets.toLocaleString(), "Unique QC-passing pockets"],
    ["Score observations", pooled.n.toLocaleString(), options.unit === "best" ? "Protein × AA best scores" : "Pocket × AA scores"],
    ["Mean ± SD", `${fmt(pooled.mean)} ± ${fmt(pooled.sd)}`, `Median ${fmt(pooled.median)}`],
  ]);
  $("#statistics-aa-body").innerHTML = result.summary.map((row) => `<tr><th>${row.aa}</th><td>${row.proteins.toLocaleString()}</td><td>${row.n.toLocaleString()}</td>${[row.mean, row.sd, row.median, row.q1, row.q3, ...row.cutoffs.map((cutoff) => cutoff.value)].map((value) => `<td>${fmt(value)}</td>`).join("")}</tr>`).join("");
  const distributions = [
    ["Selected score", METRICS[options.metric].label, result.records.map(({ row }) => row[options.metric]), "Uses the selected observation unit."],
    ["P2Rank probability", "P2Rank probability", result.records.map(({ row }) => row.probability), "Pocket quality for the selected observations after QC."],
    ["Pocket confidence", "Pocket mean pLDDT", result.records.map(({ row }) => row.mean_pocket_plddt), "Pocket quality for the selected observations after QC."],
    ["Retained pockets per protein × AA", "Number of retained pockets", result.pocketsPerPair, "Each eligible protein × AA contributes once; no zero is assigned to failed or missing results.", { integer: true }],
    ["Median selectivity", "Median(other available AAs) − target score", result.deltas, "Independent best QC-passing pockets; canonical AAs only. Requires at least one other available AA. Positive values favor the target.", { zero: true }],
  ];
  if (hasDControl(options.aa)) {
    const aa = canonicalAACode(options.aa);
    distributions.push([`L-${aa} / D-${aa} stereoselectivity`, `D-${aa} score − L-${aa} score`, result.stereo.map((row) => row.delta), "Requires successful scores and QC-passing pockets for both configurations. Positive values favor L; this is a docking comparison, not proof of binding or nonbinding.", { zero: true }]);
  }
  $("#statistics-distributions").innerHTML = distributions.map(([title, label, values, note, config]) => `<section class="panel analysis-chart-panel"><header class="panel-header"><h2>${escapeHTML(title)}</h2><span class="analysis-note">n = ${values.length.toLocaleString()}</span></header><div class="analysis-chart">${analysisHistogram(values, label, config)}</div><p class="analysis-note analysis-inset">${escapeHTML(note)}</p></section>`).join("");
  $("#statistics-correlation-body").innerHTML = result.correlations.map((row) => `<tr><th>${escapeHTML(row.xLabel)} vs ${escapeHTML(row.yLabel)}</th><td>${row.n.toLocaleString()}</td><td>${fmt(row.r)}</td><td>${fmt(row.rho)}</td></tr>`).join("");
  if (!result.correlations.some((row) => row.id === options.scatter)) options.scatter = result.correlations[0].id;
  $("#statistics-scatter-select").innerHTML = analysisOptions(result.correlations.map((row) => [row.id, `${row.xLabel} vs ${row.yLabel}`]), options.scatter);
  renderStatisticsScatter();
}

function statisticsCorrelationPairs(result, correlation) {
  if (correlation.id === "l_d") return result.stereo.filter((pair) => Number.isFinite(pair.x) && Number.isFinite(pair.y));
  const pairs = [];
  for (const { aa, row } of result.records) {
    const x = row[correlation.xKey], y = row[correlation.yKey];
    if (Number.isFinite(x) && Number.isFinite(y)) pairs.push({ x, y, protein: row.uniprot_id, aa, pocket: row.pocket });
  }
  return pairs;
}

function renderStatisticsScatter() {
  const result = analysisResults.statistics;
  const row = result?.correlations.find((entry) => entry.id === analysisState.statistics.scatter);
  if (!row) return;
  const container = $("#statistics-scatter");
  const pairs = statisticsCorrelationPairs(result, row);
  $("#statistics-scatter-note").textContent = `Pearson r ${fmt(row.r)} · Spearman ρ ${fmt(row.rho)} · ${row.n.toLocaleString()} paired observations. All eligible paired observations are drawn and used for correlations. Overlapping points may appear as a single darker mark. Hover over points for values. Undefined correlations (too few pairs or no variation) are shown as —.`;
  if (!pairs.length) {
    container.innerHTML = `<p class="analysis-empty">No paired observations pass these filters.</p>`;
    return;
  }
  const identity = row.id === "l_d";
  const geometry = analysisScatterGeometry(pairs, identity);
  // Canvas renders every observation without creating hundreds of thousands of SVG nodes.
  container.innerHTML = `<div class="statistics-scatter-layer">${analysisScatterFrame(geometry, row.xLabel, row.yLabel, identity, "", pairs.length, row.n)}<canvas aria-hidden="true"></canvas><div class="statistics-scatter-tooltip" hidden></div></div>`;
  const canvas = container.querySelector("canvas"), tooltip = container.querySelector(".statistics-scatter-tooltip");
  const scale = Math.min(3, Math.max(1, globalThis.devicePixelRatio || 1));
  canvas.width = Math.round(600 * scale);
  canvas.height = Math.round(360 * scale);
  const context = canvas.getContext("2d");
  if (!context) {
    container.innerHTML = `<p class="analysis-empty">This browser cannot draw the correlation plot (Canvas 2D is unavailable).</p>`;
    return;
  }
  context.setTransform(canvas.width / 600, 0, 0, canvas.height / 360, 0, 0);
  context.fillStyle = "rgba(117, 68, 169, 0.35)";
  // Index all points spatially for hover; this groups lookup work, not the plotted data.
  const positions = new Float32Array(pairs.length * 2), cells = new Map();
  const cellSize = 8, columns = 75;
  for (let index = 0; index < pairs.length; index++) {
    const x = geometry.x(pairs[index].x), y = geometry.y(pairs[index].y);
    positions[index * 2] = x;
    positions[index * 2 + 1] = y;
    const key = Math.floor(x / cellSize) + columns * Math.floor(y / cellSize);
    if (!cells.has(key)) cells.set(key, []);
    cells.get(key).push(index);
    context.beginPath();
    context.arc(x, y, 2.4, 0, Math.PI * 2);
    context.fill();
  }
  canvas.addEventListener("pointermove", (event) => {
    const bounds = canvas.getBoundingClientRect();
    if (!bounds.width || !bounds.height) return;
    const x = (event.clientX - bounds.left) * 600 / bounds.width;
    const y = (event.clientY - bounds.top) * 360 / bounds.height;
    const column = Math.floor(x / cellSize), cellRow = Math.floor(y / cellSize);
    let closest = -1, distance = 36;
    for (let dy = -1; dy <= 1; dy++) {
      for (let dx = -1; dx <= 1; dx++) {
        for (const index of cells.get(column + dx + columns * (cellRow + dy)) || []) {
          const d = (positions[index * 2] - x) ** 2 + (positions[index * 2 + 1] - y) ** 2;
          if (d < distance) { closest = index; distance = d; }
        }
      }
    }
    tooltip.hidden = closest < 0;
    if (closest >= 0) {
      const pair = pairs[closest];
      const label = pair.label || [pair.protein, pair.aa, pair.pocket].filter(Boolean).join(" · ") || "Observation";
      tooltip.textContent = `${label} · ${row.xLabel} ${fmt(pair.x)} · ${row.yLabel} ${fmt(pair.y)}`;
    }
  });
  canvas.addEventListener("pointerleave", () => { tooltip.hidden = true; });
}

function renderAnalysisView(view) {
  if (view === "compare") renderCompare();
  if (view === "overlap") renderOverlap();
  if (view === "statistics") scheduleStatistics();
}

let statisticsRenderTimer;
function scheduleStatistics() {
  const workspace = $("#statistics-workspace");
  const status = $("#statistics-status");
  clearTimeout(statisticsRenderTimer);
  workspace.setAttribute("aria-busy", "true");
  status.textContent = "Calculating proteome statistics… Pooled all-pocket summaries may take a few moments.";
  status.hidden = false;
  $("#statistics-download").disabled = true;
  // Allow the status to paint before the potentially large pooled calculation.
  statisticsRenderTimer = setTimeout(() => {
    let complete = false;
    try { renderStatistics(); complete = true; }
    catch (error) { console.error(error); status.textContent = "Statistics could not be calculated. Try a single AA or best-pocket mode; any results below belong to the previous selection."; }
    finally { workspace.setAttribute("aria-busy", "false"); status.hidden = complete; $("#statistics-download").disabled = !complete; }
  }, 50);
}

function bindAnalysisEvents() {
  const aaOptions = DATA_LIGANDS.filter(({ code }) => state.rawByAA.has(code)).map(({ code, name }) => [code, `${code} · ${name}`]);
  $("#compare-controls").innerHTML = analysisSelect("compare-x", "Target AA (X)", aaOptions, analysisState.compare.x)
    + analysisSelect("compare-y", "Comparison AA (Y)", aaOptions, analysisState.compare.y) + commonAnalysisControls("compare");
  $("#compare-table-controls").innerHTML = analysisSelect("compare-preference", "Preference", [["all", "Either AA"], ["x", "X scores better"], ["y", "Y scores better"]], "all")
    + analysisNumber("compare-minDelta", "Δ (Y − X) ≥", "", null, null, .05, "Any")
    + `<label class="field"><span>Protein or gene</span><input id="compare-search" type="search" placeholder="Protein or gene…" /></label>`;
  $("#overlap-controls").innerHTML = analysisSelect("overlap-tier", "Top-hit tier", [[1, "Top 1%"], [5, "Top 5%"], [10, "Top 10%"]], 5) + commonAnalysisControls("overlap");
  $("#overlap-aas").innerHTML = AMINO_ACIDS.map(({ code, name }) => `<label title="${name}"><input type="checkbox" value="${code}" ${analysisState.overlap.aas.includes(code) ? "checked" : ""}/><span>${code}</span></label>`).join("");
  $("#overlap-table-controls").innerHTML = analysisSelect("overlap-group", "Protein set", [["shared", "Shared by all selected AAs"], ["all", "All top-hit proteins"]], analysisState.overlap.group)
    + `<label class="field"><span>Protein or gene</span><input id="overlap-search" type="search" placeholder="Protein or gene…" /></label>`;
  $("#statistics-controls").innerHTML = analysisSelect("statistics-aa", "Amino acid", [["ALL", "All 20 canonical AAs"], ...aaOptions], analysisState.statistics.aa)
    + commonAnalysisControls("statistics") + analysisSelect("statistics-unit", "Observation unit", [["best", "Best pocket per protein × AA"], ["all", "All pockets"]], "best");
  for (const view of ["compare", "overlap", "statistics"]) {
    $(`#${view}-workspace`).addEventListener("change", (event) => {
      const input = event.target, prefix = `${view}-`;
      if (!input.id.startsWith(prefix)) return;
      const key = input.id.slice(prefix.length), options = analysisState[view];
      if (!(key in options)) return;
      let value = input.value;
      if (key === "p2" || key === "plddt") { value = Math.min(key === "p2" ? 1 : 100, Math.max(0, Number(value) || 0)); input.value = value; }
      if (key === "tier") value = Number(value);
      if (key === "minDelta" && value !== "" && !Number.isFinite(Number(value))) value = "";
      const previous = options[key];
      options[key] = value;
      if (view === "compare" && options.x === options.y) {
        options[key === "x" ? "y" : "x"] = previous;
        $("#compare-x").value = options.x; $("#compare-y").value = options.y;
      }
      options.page = 1;
      if (["search", "preference", "minDelta", "group"].includes(key)) {
        if (view === "compare") renderCompareTable();
        if (view === "overlap") renderOverlapTable();
      } else if (key === "scatter") renderStatisticsScatter();
      else renderAnalysisView(view);
    });
    $(`#${view}-workspace`).addEventListener("click", (event) => {
      const page = event.target.closest("[data-analysis-page]");
      if (page) {
        analysisState[view].page += Number(page.dataset.step);
        if (view === "compare") renderCompareTable(); else renderOverlapTable();
      }
      const protein = event.target.closest("[data-analysis-protein]");
      if (protein) {
        const aa = view === "compare" ? analysisState.compare.x : analysisState.overlap.aas.find((code) => qualityRanking(code, analysisState.overlap).byProtein.has(protein.dataset.analysisProtein));
        if (AMINO_ACIDS.some((entry) => entry.code === aa)) state.aa = aa;
        else if (dControlCode(aa)) state.aa = canonicalAACode(aa);
        state.metric = analysisState[view].metric;
        state.p2rank = analysisState[view].p2;
        state.plddt = analysisState[view].plddt;
        state.profileValue = Object.keys(PROFILE_VALUES).find((key) => key !== "percentile" && PROFILE_VALUES[key].metric === state.metric);
        if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
        syncControls();
        renderExplorer();
        selectProtein(protein.dataset.analysisProtein, getRankedProtein(protein.dataset.analysisProtein, state.aa, state.metric));
      }
    });
  }
  for (const view of ["compare", "overlap"]) {
    $(`#${view}-search`).addEventListener("input", (event) => {
      analysisState[view].search = event.target.value; analysisState[view].page = 1;
      if (view === "compare") renderCompareTable(); else renderOverlapTable();
    });
  }
  $("#compare-head").addEventListener("click", (event) => {
    const button = event.target.closest("[data-compare-sort]");
    if (!button) return;
    const options = analysisState.compare, key = button.dataset.compareSort;
    options.direction = options.sort === key ? -options.direction : key === "delta" ? -1 : 1;
    options.sort = key; options.page = 1;
    renderCompareTable();
  });
  $("#overlap-aas").addEventListener("change", (event) => {
    const selected = $$("#overlap-aas input:checked").map((input) => input.value);
    if (selected.length < 2) { event.target.checked = true; showToast("Choose at least two amino acids for overlap analysis"); return; }
    Object.assign(analysisState.overlap, { aas: selected, page: 1, group: "shared" });
    renderOverlap();
  });
  $("#statistics-scatter-select").addEventListener("change", (event) => { analysisState.statistics.scatter = event.target.value; renderStatisticsScatter(); });
  $("#compare-download").addEventListener("click", () => {
    const options = analysisState.compare;
    downloadText(`aa_comparison_${options.x}_${options.y}_${METRICS[options.metric].short}.tsv`, analysisTSV(
      ["uniprot_id", "gene_symbol", "aa_x", "aa_y", "metric", "score_x", "score_y", "delta_y_minus_x", "pocket_x", "pocket_y", "qc_percentile_x", "qc_percentile_y", "p2rank_min", "plddt_min"],
      filteredComparison().map((row) => [row.protein, geneSymbol(state.annotations.get(row.protein)), options.x, options.y, options.metric, row.x, row.y, row.delta, row.pocketX, row.pocketY, row.percentileX, row.percentileY, options.p2, options.plddt])));
  });
  $("#overlap-download").addEventListener("click", () => {
    const options = analysisState.overlap;
    downloadText(`aa_overlap_top${options.tier}_${options.aas.join("_")}.tsv`, analysisTSV(
      ["uniprot_id", "gene_symbol", "metric", "tier_percent", "p2rank_min", "plddt_min", "top_hit_membership", "top_hit_aa_count", "canonical_aa_count", "available_canonical_aa_count", "all_top_hit_aas", "selected_top_hit_aa_count", "selected_aa_count", "available_selected_aa_count", ...options.aas.flatMap((aa) => [`${aa}_status`, `${aa}_qc_percentile`])],
      filteredOverlap().map((row) => [row.protein, geneSymbol(state.annotations.get(row.protein)), options.metric, options.tier, options.p2, options.plddt, options.aas.filter((_, i) => row.membership[i]).join(";"), row.topHitAACount, AMINO_ACIDS.length, row.availableAACount, row.topHitAAs.join(";"), row.membership.filter(Boolean).length, options.aas.length, row.availableSelectedAACount, ...row.statuses.flatMap((status, i) => [status, row.percentiles[i]])])));
  });
  $("#statistics-download").addEventListener("click", () => {
    const options = analysisState.statistics;
    downloadText(`aa_statistics_${options.aa}_${METRICS[options.metric].short}.tsv`, analysisTSV(
      ["aa", "metric", "unit", "p2rank_min", "plddt_min", "proteins", "observations", "mean", "population_sd", "median", "q1", "q3", "top1_best_pocket_cutoff", "top5_best_pocket_cutoff", "top10_best_pocket_cutoff"],
      (analysisResults.statistics?.summary || []).map((row) => [row.aa, options.metric, options.unit, options.p2, options.plddt, row.proteins, row.n, row.mean, row.sd, row.median, row.q1, row.q3, ...row.cutoffs.map((cutoff) => cutoff.value)])));
  });
}
