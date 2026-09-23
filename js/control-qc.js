// Local experimental controls: descriptive ranking QC and a separate residue-overlap diagnostic.
const CONTROL_QC_PATH = "annotations/arabidopsis/strict_WT_single_protein_AA_controls.tsv";
const CONTROL_QC_METRICS = [
  { id: "vina", key: "vina_affinity", label: "Vina", color: "#7650a1", dash: "" },
  { id: "sfct", key: "sfct_score", label: "SFCT", color: "#437884", dash: "7 3" },
  { id: "combined", key: "vina_sfct_combined", label: "Combined 80%", color: "#303847", dash: "2 3" },
  { id: "combined50", key: "vina_sfct_combined_50", label: "Combined 50%", color: "#92795d", dash: "9 3 2 3" },
];
const CONTROL_QC_POCKET_MODES = { best: "Best-scoring pocket", matched: "Control-matched pocket" };
const controlQCState = { aa: "ALL", denominator: "eligible", pocketMode: "best", overlap: "50", p2: 0.7, plddt: 90 };
const controlQCTableState = { sort: null, direction: 1 };
const CONTROL_QC_TABLE_SORTS = {
  protein: { label: "Protein", value: (row) => row.protein },
  aa: { label: "AA", value: (row) => row.aa },
  pdb: { label: "PDB ID", value: (row) => row.pdb },
  ...Object.fromEntries(CONTROL_QC_METRICS.flatMap(({ id, label }, i) => [
    [id, { label: `${label} percentile`, value: (row) => row.percentiles[i] }],
    [`stereo_${id}`, { label: `${label} L/D delta`, value: (row) => row.stereo[i]?.delta, descending: true }],
    [`site_${id}`, { label: `${label} site coverage`, value: (row) => row.scoreOverlaps[i]?.fraction, descending: true }],
  ])),
  docking: { label: "Docking coverage", value: (row) => row.dockingStatus },
  match: { label: "Control-site match", value: (row) => row.matched?.fraction, descending: true },
};
let controlQCDataPromise = null, controlQCRequest = 0, controlQCResult = null, controlQCIndex = null;
const controlQCCache = new Map();

function parseControlTSV(text) {
  let headers = null, values = [], field = "", quoted = false, closed = false;
  const records = [];
  const finishRow = () => {
    values.push(field); field = ""; closed = false;
    if (values.some((value) => value.trim())) {
      if (!headers) headers = values.map((value) => value.trim());
      else {
        if (values.length !== headers.length) throw new Error("Unexpected column count in the control TSV.");
        records.push(Object.fromEntries(headers.map((header, i) => [header, values[i].trim()])));
      }
    }
    values = [];
  };
  const source = text.replace(/^\uFEFF/, "");
  for (let i = 0; i < source.length; i++) {
    const c = source[i];
    if (quoted) {
      if (c === '"' && source[i + 1] === '"') { field += '"'; i++; }
      else if (c === '"') { quoted = false; closed = true; }
      else field += c;
    } else if (c === "\t") { values.push(field); field = ""; closed = false; }
    else if (c === "\r" || c === "\n") {
      if (c === "\r" && source[i + 1] === "\n") i++;
      finishRow();
    } else if (c === '"' && !field && !closed) quoted = true;
    else {
      if (closed || c === '"') throw new Error("Invalid quoting in the control TSV.");
      field += c;
    }
  }
  if (quoted) throw new Error("Unterminated quoted field in the control TSV.");
  if (field || values.length || closed) finishRow();
  const required = ["AA", "UniProt_accession", "PDB_ID", "control_eligible", "UniProt_mapping_status", "UniProt_contact_positions"];
  if (!headers || new Set(headers).size !== headers.length || required.some((key) => !headers.includes(key))) {
    throw new Error("Control TSV is missing required columns or has duplicate headers.");
  }
  const canonical = new Set(AMINO_ACIDS.map(({ code }) => code)), seen = new Set(), controls = [];
  for (const row of records) {
    const aa = row.AA.toUpperCase(), protein = row.UniProt_accession.toUpperCase();
    if (!/^(1|true|yes)$/i.test(row.control_eligible) || !canonical.has(aa)) continue;
    if (!protein) throw new Error("An eligible control has no UniProt accession.");
    const key = `${aa}|${protein}`;
    if (seen.has(key)) throw new Error(`Duplicate AA–protein control: ${aa} / ${protein}. Use one representative per pair.`);
    seen.add(key);
    const tokens = row.UniProt_contact_positions.split(";").map((value) => value.trim());
    const validPositions = tokens.length > 0 && tokens.every((value) => /^[1-9]\d*$/.test(value) && Number.isSafeInteger(Number(value)));
    const positions = row.UniProt_mapping_status.toLowerCase() === "contact_overlap" && validPositions
      ? [...new Set(tokens.map(Number))] : [];
    controls.push({ aa, protein, pdb: row.PDB_ID, positions, source: row });
  }
  if (!controls.length) throw new Error("No eligible canonical AA controls were found.");
  return { controls, total: records.length, excluded: records.length - controls.length };
}

async function loadControlQCData() {
  if (!controlQCDataPromise) {
    controlQCDataPromise = (async () => {
      const response = await fetch(CONTROL_QC_PATH);
      if (!response.ok) throw new Error(`Control file could not be loaded (HTTP ${response.status}).`);
      return parseControlTSV(await response.text());
    })().catch((error) => { controlQCDataPromise = null; throw error; });
  }
  return controlQCDataPromise;
}

function controlQCPocket(row) {
  if (!row || !String(row.protein || "").startsWith(`AF-${row.uniprot_id}-F1-model_`)) return null;
  const tokens = String(row.residue_ids || "").trim().split(/\s+/);
  if (!tokens.length || !tokens.every((token) => /^A_[1-9]\d*$/.test(token))) return null;
  const residues = new Set(tokens.map((token) => Number(token.slice(2))));
  const identity = [row.protein, row.pocket, [...residues].sort((a, b) => a - b).join(","), row.center_x, row.center_y, row.center_z].join("|");
  return { row, residues, identity };
}

function controlQCOverlap(control, pocket) {
  if (!control.positions.length || !pocket) return null;
  const hits = control.positions.filter((position) => pocket.residues.has(position));
  return { ...pocket, hits, fraction: hits.length / control.positions.length };
}

function controlQCSourceIndex(data) {
  const sources = [...state.rawByAA.entries()];
  if (controlQCIndex?.data === data && sources.length === controlQCIndex.sources.length
    && sources.every(([aa, rows], i) => aa === controlQCIndex.sources[i][0] && rows === controlQCIndex.sources[i][1])) return controlQCIndex;
  controlQCCache.clear();
  const proteins = new Set(data.controls.map((control) => control.protein));
  const pairs = new Map(), pockets = new Map();
  for (const [aa, rows] of sources) {
    for (const row of rows) {
      if (!proteins.has(row.uniprot_id) || !isSuccessfulResult(row)) continue;
      const pairKey = `${aa}|${row.uniprot_id}`;
      if (!pairs.has(pairKey)) pairs.set(pairKey, []);
      pairs.get(pairKey).push(row);
      // D controls must not change the existing L-AA experimental-site reference.
      if (!AMINO_ACIDS.some((entry) => entry.code === aa)) continue;
      const pocket = controlQCPocket(row);
      if (!pocket) continue;
      if (!pockets.has(row.uniprot_id)) pockets.set(row.uniprot_id, new Map());
      const previous = pockets.get(row.uniprot_id).get(pocket.identity);
      if (!previous || (Number.isFinite(row.rank) ? row.rank : Infinity) < (Number.isFinite(previous.row.rank) ? previous.row.rank : Infinity)) {
        pockets.get(row.uniprot_id).set(pocket.identity, pocket);
      }
    }
  }
  controlQCIndex = { data, sources, pairs, pockets };
  return controlQCIndex;
}

// Locate this pocket's score among OTHER proteins' best QC-passing scores.
// This preserves the full-proteome reference and avoids counting the control twice.
function controlQCPocketPercentile(ranking, row, metric) {
  return pocketProteomePosition(ranking, row, metric)?.proteome_percentile ?? NaN;
}

function calculateControlQC(data, options) {
  const mode = options.pocketMode || "best", threshold = options.overlap || "50";
  const index = controlQCSourceIndex(data), key = [options.aa, options.p2, options.plddt, mode, threshold].join("|");
  if (controlQCCache.has(key)) return controlQCCache.get(key);
  const controls = data.controls.filter((control) => options.aa === "ALL" || control.aa === options.aa);
  const rankings = new Map(), dRankings = new Map();
  for (const aa of new Set(controls.map((control) => control.aa))) {
    rankings.set(aa, CONTROL_QC_METRICS.map(({ key: metric }) => qualityRanking(aa, { ...options, metric })));
    if (dControlCode(aa)) dRankings.set(aa, CONTROL_QC_METRICS.map(({ key: metric }) => qualityRanking(dControlCode(aa), { ...options, metric })));
  }
  const rows = controls.map((control) => {
    const aaRankings = rankings.get(control.aa);
    const raw = index.pairs.get(`${control.aa}|${control.protein}`) || [];
    const passing = raw.filter((row) => qualityPasses(row, options));
    const allPockets = [...(index.pockets.get(control.protein)?.values() || [])];
    const mapped = control.positions.length > 0;
    // Select the experimental-site match using geometry only, BEFORE score/QC availability.
    const overlaps = allPockets.map((pocket) => controlQCOverlap(control, pocket)).filter(Boolean).sort((a, b) =>
      b.hits.length - a.hits.length
      || (Number.isFinite(a.row.rank) ? a.row.rank : Infinity) - (Number.isFinite(b.row.rank) ? b.row.rank : Infinity)
      || a.identity.localeCompare(b.identity));
    const candidate = overlaps[0] || null;
    const matched = controlQCCaptured(candidate, threshold) ? candidate : null;
    const matchedRows = matched ? raw.filter((row) => samePocketGeometry(matched.row, row)) : [];
    const matchedPassing = matchedRows.filter((row) => qualityPasses(row, options));
    const scoreRows = CONTROL_QC_METRICS.map(({ key: metric }, i) => mode === "best" ? aaRankings[i].byProtein.get(control.protein)
      : matchedPassing.filter((row) => hasUsableScore(row, metric)).sort((a, b) => a[metric] - b[metric])[0]);
    const percentiles = scoreRows.map((row, i) => row ? controlQCPocketPercentile(aaRankings[i], row, CONTROL_QC_METRICS[i].key) : NaN);
    const eligible = percentiles.every(Number.isFinite);
    let dockingStatus = !raw.length ? "No successful result (missing / failed)" : !passing.length ? "No pocket passes QC"
      : !eligible ? "Missing finite score(s)" : "Eligible for all scores";
    const matchStatus = !mapped ? "Unknown: no confirmed residue mapping" : !allPockets.length ? "Unknown: no usable retained pocket"
      : !matched ? "No pocket meets site-overlap threshold" : !matchedRows.length ? "Matched pocket has no successful target-AA result"
      : !matchedPassing.length ? "Matched pocket fails QC" : "Control site matched";
    if (mode === "matched" && !eligible && matchStatus !== "Control site matched") dockingStatus = matchStatus;
    const scoreOverlaps = scoreRows.map((row) => controlQCOverlap(control, controlQCPocket(row)));
    const dRaw = index.pairs.get(`${dControlCode(control.aa)}|${control.protein}`) || [];
    const stereo = scoreRows.map((lRow, i) => {
      if (!lRow || !dControlCode(control.aa)) return null;
      const metric = CONTROL_QC_METRICS[i].key;
      const dRow = mode === "best" ? dRankings.get(control.aa)?.[i].byProtein.get(control.protein)
        : dRaw.find((row) => qualityPasses(row, options) && hasUsableScore(row, metric) && samePocketGeometry(lRow, row));
      if (!dRow) return null;
      return { lRow, dRow, delta: dRow[metric] - lRow[metric] };
    });
    return { ...control, percentiles, scoreRows, populations: aaRankings.map((ranking) => ranking.rows.length), eligible,
      dockingStatus, mapped, candidate, matched, matchStatus, scoreOverlaps, stereo };
  });
  const result = { rows, options: { aa: options.aa, p2: options.p2, plddt: options.plddt, pocketMode: mode, overlap: threshold }, excluded: data.excluded };
  if (controlQCCache.size >= 4) controlQCCache.delete(controlQCCache.keys().next().value);
  controlQCCache.set(key, result);
  return result;
}

function controlQCRecoveryRows(result, denominator) {
  return denominator === "all" ? result.rows : result.rows.filter((row) => row.eligible);
}

function controlQCRecovery(result, denominator) {
  const rows = controlQCRecoveryRows(result, denominator);
  return { n: rows.length, curves: CONTROL_QC_METRICS.map((metric, i) => ({ ...metric,
    percentiles: rows.map((row) => row.percentiles[i]).filter(Number.isFinite).sort((a, b) => a - b),
    counts: [1, 5, 10].map((tier) => rows.filter((row) => Number.isFinite(row.percentiles[i]) && row.percentiles[i] <= tier).length),
  })) };
}

function controlQCFraction(count, total) {
  return total ? `${fmt(100 * count / total, 1)}% (${count}/${total})` : "— (0 evaluable)";
}

function controlQCStereoSummary(result, options) {
  const selected = controlQCRecoveryRows(result, options.denominator);
  const chiral = selected.filter((row) => dControlCode(row.aa));
  return { selected: selected.length, glycine: selected.length - chiral.length, summaries: CONTROL_QC_METRICS.map((metric, i) => {
    const values = chiral.map((row) => row.stereo[i]?.delta).filter(Number.isFinite);
    return { ...metric, n: values.length, lPreferred: values.filter((value) => value > 0).length,
      dPreferred: values.filter((value) => value < 0).length, tied: values.filter((value) => value === 0).length,
      missing: chiral.length - values.length, median: median(values) };
  }) };
}

function renderControlQCStereo(result, options) {
  const summary = controlQCStereoSummary(result, options);
  $("#control-qc-stereo-summary").innerHTML = summary.summaries.map((row) => `<tr><th>${row.label}</th><td>${controlQCFraction(row.lPreferred, row.n)}</td><td>${controlQCFraction(row.dPreferred, row.n)}</td><td>${row.tied}</td><td>${row.missing}</td><td>${fmt(row.median)}</td></tr>`).join("");
  $("#control-qc-stereo-note").textContent = `Within the ${summary.selected} controls in the selected recovery denominator. ${summary.glycine} glycine controls excluded (no L/D pair). Δ = D score − L score; positive favors L. Percentages use successful L/D pairs for each metric, not missing scores. Both configurations must pass P2Rank ≥ ${options.p2} and pLDDT ≥ ${options.plddt}. ${result.options.pocketMode === "matched" ? "L and D use the same control-matched pocket and receptor geometry; no best-pocket fallback." : "L and D use independently selected best QC-passing pockets, which may differ."} This summary does not filter or change the recovery curves. D-AA docking is a stereochemical comparison, not a proven nonbinding control or an FDR estimate.`;
}

function controlQCCaptured(overlap, threshold) {
  return Boolean(overlap && overlap.hits.length > 0 && (threshold === "any" || overlap.fraction >= Number(threshold) / 100));
}

function controlQCChart(recovery) {
  if (!recovery.n) return '<p class="analysis-empty">No controls in this denominator. Lower the quality thresholds or choose All known controls.</p>';
  const width = 610, height = 412, left = 66, top = 24, plotWidth = 516, plotHeight = 264;
  const x = (value) => left + plotWidth * value / 100, y = (value) => top + plotHeight * (1 - value / 100);
  const ticks = [0, 20, 40, 60, 80, 100];
  const grid = ticks.map((v) => `<path d="M${x(v)} ${top}V${y(0)} M${left} ${y(v)}H${x(100)}" stroke="#e8e1da" fill="none"/><text x="${x(v)}" y="${y(0) + 20}" text-anchor="middle">${v}</text><text x="${left - 10}" y="${y(v) + 4}" text-anchor="end">${v}</text>`).join("");
  const curves = recovery.curves.map((curve) => {
    let path = `M${x(0)} ${y(0)}`, count = 0;
    for (let i = 0; i < curve.percentiles.length;) {
      const percentile = curve.percentiles[i];
      while (i < curve.percentiles.length && curve.percentiles[i] === percentile) { count++; i++; }
      path += `H${x(percentile).toFixed(3)}V${y(100 * count / recovery.n).toFixed(3)}`;
    }
    path += `H${x(100)}`;
    return `<path d="${path}" fill="none" stroke="${curve.color}" stroke-width="2.5" stroke-dasharray="${curve.dash}"><title>${curve.label}: ${curve.percentiles.length}/${recovery.n} recovered by Top 100%</title></path>`;
  }).join("");
  const legend = [...recovery.curves, { label: "Random (y = x)", color: "#999099", dash: "5 5" }].map((item, i) => {
    const lx = 40 + (i % 3) * 190, ly = 365 + Math.floor(i / 3) * 24;
    return `<path d="M${lx} ${ly}h22" stroke="${item.color}" stroke-width="2.5" stroke-dasharray="${item.dash}"/><text x="${lx + 28}" y="${ly + 4}">${item.label}</text>`;
  }).join("");
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Known AA–protein control recovery versus top proteome percentile for Vina, SFCT, Combined 80% and Combined 50%"><title>Positive-control recovery</title><desc>Each control is one AA–protein pair. Curves show cumulative recovery in percent. The diagonal is the ideal random-ranking reference. Exact Top 1, 5 and 10 percent recoveries are in the adjacent table.</desc>${grid}<path d="M${x(0)} ${y(0)}L${x(100)} ${y(100)}" stroke="#999099" stroke-width="1.5" stroke-dasharray="5 5"/>${curves}<text x="${left + plotWidth / 2}" y="337" text-anchor="middle">Top proteome percentile (%)</text><text transform="translate(17 ${top + plotHeight / 2}) rotate(-90)" text-anchor="middle">Known AA–protein controls recovered (%)</text>${legend}</svg>`;
}

function controlQCProfileButton(row) {
  const available = state.metadata.has(row.protein);
  return `<button type="button" class="open-row" data-control-protein="${escapeHTML(row.protein)}" data-control-aa="${row.aa}" ${available ? "" : "disabled"} title="${available ? "Open protein profile" : "No loaded results for this control protein"}" aria-label="Open ${escapeHTML(row.protein)} profile">→</button>`;
}

function controlQCOrderedRows(result, options = controlQCState) {
  const rows = [...controlQCRecoveryRows(result, options.denominator)], spec = CONTROL_QC_TABLE_SORTS[controlQCTableState.sort];
  if (!spec) return rows; // Preserve file order until a column is chosen.
  const missing = (value) => value === null || value === undefined || value === ""
    || (typeof value === "number" && !Number.isFinite(value));
  return rows.sort((a, b) => {
    const av = spec.value(a), bv = spec.value(b), am = missing(av), bm = missing(bv);
    // Unknown scores/overlaps stay last in either direction; never treat them as zero.
    if (am !== bm) return am ? 1 : -1;
    const difference = am ? 0 : typeof av === "number" ? av - bv
      : String(av).localeCompare(String(bv), undefined, { numeric: true, sensitivity: "base" });
    return difference * controlQCTableState.direction || a.aa.localeCompare(b.aa)
      || a.protein.localeCompare(b.protein) || String(a.pdb || "").localeCompare(String(b.pdb || ""));
  });
}

function controlQCPDBLink(pdb) {
  const id = String(pdb || "").trim();
  if (!id) return "—";
  return `<a class="control-qc-pdb-link" href="https://www.rcsb.org/structure/${encodeURIComponent(id)}" target="_blank" rel="noopener noreferrer" title="Open ${escapeHTML(id)} on RCSB PDB">${escapeHTML(id)}</a>`;
}

function controlQCStereoCell(row) {
  if (!dControlCode(row.aa)) return '<td><small class="analysis-gene">Not applicable (glycine)</small></td>';
  return `<td>${row.stereo.map((pair, i) => `<small class="analysis-gene" title="${pair ? `L: ${fmt(pair.lRow[CONTROL_QC_METRICS[i].key])} (${escapeHTML(pair.lRow.pocket)}); D: ${fmt(pair.dRow[CONTROL_QC_METRICS[i].key])} (${escapeHTML(pair.dRow.pocket)})` : "Successful QC-passing L/D pair required"}">${CONTROL_QC_METRICS[i].label}: ${pair ? fmt(pair.delta) : "Missing"}</small>`).join("")}</td>`;
}

function renderControlQCTable(result, options = controlQCState) {
  const selected = controlQCTableState.sort;
  const sortDirection = controlQCTableState.direction === 1 ? "ascending" : "descending";
  const sortButton = (key, label = CONTROL_QC_TABLE_SORTS[key].label) => `<button type="button" class="sort-button" data-control-sort="${key}" aria-label="Sort by ${CONTROL_QC_TABLE_SORTS[key].label}" title="${key === "match" ? "Sort by matched experimental-residue overlap" : `Sort by ${CONTROL_QC_TABLE_SORTS[key].label}`}">${label} <span aria-hidden="true">${selected === key ? controlQCTableState.direction === 1 ? "↑" : "↓" : "↕"}</span></button>`;
  const header = (key) => `<th scope="col" aria-sort="${selected === key ? sortDirection : "none"}">${sortButton(key)}</th>`;
  const stereoHeader = `<th scope="col" aria-sort="${selected?.startsWith("stereo_") ? sortDirection : "none"}">L/D Δ (D − L)<div class="control-qc-site-sorts">${CONTROL_QC_METRICS.map(({ id, label }) => sortButton(`stereo_${id}`, label)).join("")}</div></th>`;
  $("#control-qc-head").innerHTML = `<tr>${["protein", "pdb", "aa", ...CONTROL_QC_METRICS.map(({ id }) => id), "docking"].map(header).join("")}${stereoHeader}<th scope="col" aria-sort="${selected?.startsWith("site_") ? sortDirection : "none"}">Analyzed-site coverage<div class="control-qc-site-sorts">${CONTROL_QC_METRICS.map(({ id, label }) => sortButton(`site_${id}`, label)).join("")}</div></th>${header("match")}<th scope="col"><span class="sr-only">Protein profile</span></th></tr>`;
  const rows = controlQCOrderedRows(result, options);
  $("#control-qc-table-note").textContent = `Showing ${rows.length} of ${result.rows.length} controls for the selected AA set, using the same denominator as the recovery plot. ${options.denominator === "all" ? "All known controls includes unmatched, QC-excluded and missing results as shown in the plot denominator." : "Only controls with usable scores for all four metrics after the selected pocket mode and quality filters are shown."}`;
  $("#control-qc-body").innerHTML = rows.map((row) => `<tr><td>${proteinTableIdentity(row.protein)}</td><td>${controlQCPDBLink(row.pdb)}</td><td>${escapeHTML(row.aa)}</td>${row.percentiles.map((p, i) => `<td>${Number.isFinite(p) ? `${fmt(p, 2)}%<small class="analysis-gene">${escapeHTML(row.scoreRows[i].pocket)} · score ${fmt(row.scoreRows[i][CONTROL_QC_METRICS[i].key])}</small>` : "—"}</td>`).join("")}<td>${escapeHTML(row.dockingStatus)}</td>${controlQCStereoCell(row)}<td>${row.scoreOverlaps.map((overlap, i) => `<small class="analysis-gene">${CONTROL_QC_METRICS[i].label}: ${overlap ? `${fmt(100 * overlap.fraction, 1)}% (${overlap.hits.length}/${row.positions.length})` : "Missing"}</small>`).join("")}</td><td>${row.matched ? `${escapeHTML(row.matched.row.pocket)} · ${fmt(100 * row.matched.fraction, 1)}%` : "—"}<small class="analysis-gene">${escapeHTML(row.matchStatus)}</small></td><td>${controlQCProfileButton(row)}</td></tr>`).join("") || '<tr><td colspan="12" class="analysis-empty">No controls match the current QC settings and recovery denominator.</td></tr>';
}

function renderControlQCResults(result, options) {
  renderControlQCStereo(result, options);
  const recovery = controlQCRecovery(result, options.denominator), eligible = result.rows.filter((row) => row.eligible).length;
  $("#control-qc-chart").innerHTML = controlQCChart(recovery);
  $("#control-qc-summary").innerHTML = recovery.curves.map((curve) => `<tr><th><span class="control-qc-key" style="--qc-color:${curve.color}"></span>${curve.label}</th>${curve.counts.map((count) => `<td>${controlQCFraction(count, recovery.n)}</td>`).join("")}</tr>`).join("");
  $("#control-qc-coverage").textContent = `${CONTROL_QC_POCKET_MODES[result.options.pocketMode]} · ${result.rows.length} known controls · ${eligible} eligible for all four scores · ${result.rows.length - eligible} unmatched, missing a score or excluded by QC. Recovery denominator: ${recovery.n}. Percentages count AA–protein pairs, not unique proteins.${result.excluded ? ` ${result.excluded} source rows excluded as ineligible or noncanonical.` : ""}`;
  const modeNote = result.options.pocketMode === "matched"
    ? "Each control uses the pocket most closely matching its experimental residues, chosen without docking scores. Its score is placed among other proteins’ best-pocket scores; the control’s own best score is excluded. This is known-site-conditioned QC, not blind pocket selection."
    : "Each control uses its best quality-passing pocket independently for Vina, SFCT, Combined 80% and Combined 50%.";
  $("#control-qc-curve-note").textContent = `${modeNote} Reference: each AA’s full best-pocket proteome after P2Rank ≥ ${result.options.p2} and pLDDT ≥ ${result.options.plddt}. Lower scores are better; ties use UniProt order. ${options.denominator === "all" ? "Unmatched/missing results remain in the denominator without an invented rank." : "All curves use the same controls with available scores for all four metrics; unmatched/missing controls are excluded and reported above."} Combined 50% reweights the saved 80%-selected pose; other poses are unavailable. The y = x line is an ideal random-ranking reference, not adjusted for missing coverage or control-site matching.`;
  const mapped = result.rows.filter((row) => row.mapped).length;
  const siteSummaries = CONTROL_QC_METRICS.map((metric, i) => {
    const sites = result.rows.filter((row) => row.scoreOverlaps[i]);
    return { label: metric.label, n: sites.length, captured: sites.filter((row) => controlQCCaptured(row.scoreOverlaps[i], options.overlap)).length };
  });
  $("#control-qc-pocket-summary").innerHTML = `<div class="table-wrap"><table class="data-table control-qc-summary-table"><thead><tr><th>Score</th><th>Captured / evaluable</th><th>Evaluable / known</th></tr></thead><tbody>${siteSummaries.map((site) => `<tr><th>${site.label}</th><td>${controlQCFraction(site.captured, site.n)}</td><td>${site.n}/${result.rows.length}</td></tr>`).join("")}</tbody></table></div><p class="analysis-note">${mapped}/${result.rows.length} controls have mapped experimental residues; ${result.rows.filter((row) => row.matched).length}/${mapped} mapped sites have a pocket matching the overlap threshold before target-AA score/QC checks. Evaluable requires a usable score and mapped residues for the analyzed pocket; the remaining controls are unknown/excluded from that metric’s site fraction.${result.options.pocketMode === "matched" ? " In control-matched mode, overlap is enforced during pocket selection, so 100% among evaluable matches is expected—not an independent recovery test." : ""}</p>`;
  renderControlQCTable(result, options);
}

async function renderControlQC() {
  const request = ++controlQCRequest;
  const options = { ...controlQCState };
  const panel = $("#control-qc-panel"), status = $("#control-qc-status");
  panel.setAttribute("aria-busy", "true");
  status.hidden = false; status.textContent = "Loading and evaluating experimental controls…";
  $("#control-qc-results").hidden = true;
  $("#control-qc-download").disabled = true;
  $("#control-qc-retry").hidden = true;
  controlQCResult = null;
  try {
    const data = await loadControlQCData();
    if (request !== controlQCRequest) return;
    $("#control-qc-aa").innerHTML = analysisOptions([["ALL", `All control AAs (${data.controls.length})`], ...AMINO_ACIDS.map(({ code, name }) => [code, `${code} · ${name} (${data.controls.filter((row) => row.aa === code).length})`])], options.aa);
    const result = calculateControlQC(data, options);
    renderControlQCResults(result, options);
    controlQCResult = { result, options };
    status.hidden = true;
    $("#control-qc-results").hidden = false;
    $("#control-qc-download").disabled = !controlQCRecoveryRows(result, options.denominator).length;
  } catch (error) {
    if (request !== controlQCRequest) return;
    console.error(error);
    status.textContent = `Positive-control QC unavailable: ${error.message} Other tabs remain available.`;
    $("#control-qc-retry").hidden = false;
  } finally {
    if (request === controlQCRequest) panel.setAttribute("aria-busy", "false");
  }
}

function downloadControlQC() {
  if (!controlQCResult) return;
  const { result, options } = controlQCResult;
  const recovery = controlQCRecovery(result, options.denominator);
  const siteCounts = CONTROL_QC_METRICS.map((_, i) => result.rows.filter((row) => row.scoreOverlaps[i]).length);
  const headers = ["aa", "uniprot_id", "pdb_id", "pocket_selection", "docking_status", "eligible_all_scores", "recovery_denominator", "recovery_n", "included_in_recovery_denominator",
    ...CONTROL_QC_METRICS.flatMap(({ key }) => [`${key}_score`, `${key}_pocket`, `${key}_proteome_n`, `${key}_percentile`, `${key}_top_1`, `${key}_top_5`, `${key}_top_10`, `${key}_matched_positions`, `${key}_site_coverage_fraction`, `${key}_site_captured`, `${key}_site_evaluable_n`, `${key}_d_score`, `${key}_d_pocket`, `${key}_d_minus_l_delta`, `${key}_l_preferred`]),
    "control_match_status", "mapped_experimental_positions", "site_threshold", "control_matched_pocket", "control_matched_model", "control_matched_positions", "control_matched_fraction", "experimental_residue_count", "percentile_reference", "p2rank_min", "plddt_min", "source"];
  const rows = controlQCOrderedRows(result, options).map((row) => [row.aa, row.protein, row.pdb, result.options.pocketMode, row.dockingStatus, row.eligible, options.denominator, recovery.n, options.denominator === "all" || row.eligible,
    ...CONTROL_QC_METRICS.flatMap(({ key }, i) => [row.scoreRows[i]?.[key], row.scoreRows[i]?.pocket, row.populations[i], row.percentiles[i], ...[1, 5, 10].map((tier) => Number.isFinite(row.percentiles[i]) ? row.percentiles[i] <= tier : null), row.scoreOverlaps[i]?.hits.join(";"), row.scoreOverlaps[i]?.fraction, row.scoreOverlaps[i] ? controlQCCaptured(row.scoreOverlaps[i], options.overlap) : null, siteCounts[i], row.stereo[i]?.dRow[key], row.stereo[i]?.dRow.pocket, row.stereo[i]?.delta, row.stereo[i] ? row.stereo[i].delta > 0 : null]),
    row.matchStatus, row.positions.join(";"), options.overlap === "any" ? "at_least_one_residue" : `at_least_${options.overlap}_percent`, row.matched?.row.pocket, row.matched?.row.protein, row.matched?.hits.join(";"), row.matched?.fraction, row.mapped ? row.positions.length : null, "other_proteins_best_qc_pocket_scores", result.options.p2, result.options.plddt, CONTROL_QC_PATH]);
  downloadText(`positive_control_qc_${options.aa}_${result.options.pocketMode}_overlap_${options.overlap}_p2_${options.p2}_plddt_${options.plddt}.tsv`, analysisTSV(headers, rows));
}

function bindControlQCEvents() {
  $("#control-qc-head").addEventListener("click", (event) => {
    const button = event.target.closest("[data-control-sort]");
    if (!button || !controlQCResult) return;
    const key = button.dataset.controlSort;
    if (!Object.hasOwn(CONTROL_QC_TABLE_SORTS, key)) return;
    controlQCTableState.direction = controlQCTableState.sort === key ? -controlQCTableState.direction
      : CONTROL_QC_TABLE_SORTS[key].descending ? -1 : 1;
    controlQCTableState.sort = key;
    renderControlQCTable(controlQCResult.result, controlQCResult.options);
    $(`#control-qc-head [data-control-sort="${key}"]`)?.focus({ preventScroll: true });
  });
  for (const key of ["aa", "denominator", "pocketMode", "overlap", "p2", "plddt"]) {
    $(`#control-qc-${key}`).addEventListener("change", (event) => {
      let value = event.target.value;
      if (key === "p2" || key === "plddt") {
        value = Math.min(key === "p2" ? 1 : 100, Math.max(0, Number(value) || 0));
        event.target.value = value;
      }
      controlQCState[key] = value;
      renderControlQC();
    });
  }
  $("#control-qc-retry").addEventListener("click", () => renderControlQC());
  $("#control-qc-download").addEventListener("click", downloadControlQC);
  $("#control-qc-body").addEventListener("click", (event) => {
    const button = event.target.closest("[data-control-protein]");
    if (!button) return;
    if (!state.metadata.has(button.dataset.controlProtein)) { showToast("No loaded results for this control protein"); return; }
    state.aa = button.dataset.controlAa;
    const options = controlQCResult?.options || controlQCState;
    state.p2rank = options.p2; state.plddt = options.plddt;
    state.profileValue = Object.keys(PROFILE_VALUES).find((key) => !PROFILE_VALUES[key].percentile && PROFILE_VALUES[key].metric === state.metric);
    const control = controlQCResult?.result.rows.find((row) => row.protein === button.dataset.controlProtein && row.aa === state.aa);
    const metricIndex = CONTROL_QC_METRICS.findIndex((metric) => metric.key === state.metric);
    const pocket = control?.scoreRows[metricIndex] || (options.pocketMode === "matched" ? control?.matched?.row : null);
    if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
    syncControls(); renderExplorer(); selectProtein(button.dataset.controlProtein, pocket);
  });
}
