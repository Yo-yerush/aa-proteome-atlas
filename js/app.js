const AMINO_ACIDS = [
  { code: "ALA", file: "ala", name: "Alanine" },
  { code: "ARG", file: "arg", name: "Arginine" },
  { code: "ASN", file: "asn", name: "Asparagine" },
  { code: "ASP", file: "asp", name: "Aspartate" },
  { code: "CYS", file: "cys", name: "Cysteine" },
  { code: "GLN", file: "gln", name: "Glutamine" },
  { code: "GLU", file: "glu", name: "Glutamate" },
  { code: "GLY", file: "gly", name: "Glycine" },
  { code: "HIS", file: "his", name: "Histidine" },
  { code: "ILE", file: "ile", name: "Isoleucine" },
  { code: "LEU", file: "leu", name: "Leucine" },
  { code: "LYS", file: "lys", name: "Lysine" },
  { code: "MET", file: "met", name: "Methionine" },
  { code: "PHE", file: "phe", name: "Phenylalanine" },
  { code: "PRO", file: "pro", name: "Proline" },
  { code: "SER", file: "ser", name: "Serine" },
  { code: "THR", file: "thr", name: "Threonine" },
  { code: "TRP", file: "trp", name: "Tryptophan" },
  { code: "TYR", file: "tyr", name: "Tyrosine" },
  { code: "VAL", file: "val", name: "Valine" },
];

const D_AA_CONTROLS = AMINO_ACIDS.filter(({ code }) => code !== "GLY").map((aa) => ({
  code: `D${aa.code}`, file: `d${aa.file}`, name: `D-${aa.name}`, canonicalCode: aa.code,
}));
const DATA_LIGANDS = [...AMINO_ACIDS, ...D_AA_CONTROLS];
// Retained only for legacy Vina/SFCT sites predating the split compact layout.
const LEGACY_DATA_LIGANDS = [...AMINO_ACIDS, D_AA_CONTROLS.find(({ code }) => code === "DMET")];

function canonicalAACode(code) {
  return D_AA_CONTROLS.find((aa) => aa.code === code)?.canonicalCode || code;
}

function dControlCode(code) {
  return D_AA_CONTROLS.find((aa) => aa.canonicalCode === canonicalAACode(code))?.code || null;
}

function hasDControl(code = state.aa) {
  return Boolean(dControlCode(code) && state.rawByAA.has(dControlCode(code)));
}
const UNIPROT_ANNOTATION_PATH = ORGANISM.annotations;

const METRICS = {
  vina_sfct_combined: { label: "Combined 80%", short: "combined80", digits: 3, nearWindow: 0.10 },
  vina_sfct_combined_50: { label: "Combined 50%", short: "combined50", digits: 3, nearWindow: 0.10 },
  vina_affinity: { label: "AutoDock Vina", short: "Vina", digits: 3, nearWindow: 0.25 },
  sfct_score: { label: "OnionNet-SFCT", short: "SFCT", digits: 3, nearWindow: 0.15 },
};

const SORT_KEYS = new Set(["protein", "pocket", "score", "qphi_kT", "aa_rank", "competitors", "percentile", "p2rank", "plddt", "delta", "z", "stereo"]);

const PROFILE_VALUES = {
  raw_combined: { label: "Combined 80%", digits: 3, metric: "vina_sfct_combined" },
  raw_combined_50: { label: "Combined 50%", digits: 3, metric: "vina_sfct_combined_50" },
  raw_vina: { label: "Vina score", axisLabel: "Vina score (predicted affinity, kcal/mol)", digits: 3, metric: "vina_affinity" },
  raw_sfct: { label: "SFCT score", digits: 3, metric: "sfct_score" },
  percentile: { label: "Percentile · Combined 80%", digits: 1, metric: "vina_sfct_combined", percentile: true },
  percentile_50: { label: "Percentile · Combined 50%", digits: 1, metric: "vina_sfct_combined_50", percentile: true },
};

const MOLSTAR_VERSION = "5.11.0";
const MOLSTAR_ASSET_ROOT = `https://cdn.jsdelivr.net/npm/molstar@${MOLSTAR_VERSION}/build/viewer`;

const NUMERIC_FIELDS = new Set([
  "rank", "score", "probability", "center_x", "center_y", "center_z",
  "sas_points", "surf_atoms", "n_pocket_residues", "n_plddt_matched",
  "mean_pocket_plddt", "min_pocket_plddt", "fraction_plddt_ge70",
  "fraction_plddt_ge90", "sfct_best_pose", "sfct_vina_score", "sfct_score", "vina_affinity",
  "vina_sfct_combined", "vina_sfct_combined_50", "sfct_n_poses",
]);

const state = {
  unavailableDatasets: [],
  rawByAA: new Map(),
  rankingCache: new Map(),
  metadata: new Map(),
  annotations: new Map(),
  aa: "ALA",
  metric: "vina_sfct_combined_50",
  top: 100,
  p2rank: 0.7,
  plddt: 90,
  aaRank: 20,
  minDelta: null,
  maxCompetitors: 19,
  requireLPreference: false,
  sortKey: "percentile",
  sortDirection: "asc",
  pocketMode: "best",
  profileOrder: "score_asc",
  profileValue: "raw_combined_50",
  profileNormalize: false,
  search: "",
  selectedProtein: null,
  selectedPocket: null,
  profilePocket: null,
  profilePocketAnchor: null,
  showLigand: true,
  ligandPoseMode: "auto",
  showPocketCloud: true,
  pocketPointColor: "gold",
  proteinStyle: "cartoon",
  proteinColor: "plddt",
  showPocketHighlight: true,
  showPocketSticks: false,
  pageSize: 10,
  page: 1,
  filtered: [],
  distributionSource: "all",
  currentView: "explorer",
};

const $ = (selector) => document.querySelector(selector);
const $$ = (selector) => [...document.querySelectorAll(selector)];

let molstarAssetsPromise = null;
let molstarViewerPromise = null;
let molstarViewer = null;
let molstarLoadedModel = null;
let molstarProteinStructure = null;
let molstarPocketComponent = null;
let molstarPocketSticks = null;
let molstarLatestRequest = null;
let molstarUpdateQueue = Promise.resolve();

function parseTSV(text) {
  const lines = text.trim().split(/\r?\n/);
  const headers = lines.shift().split("\t").map((header) => header.replace(/^\uFEFF/, ""));
  return lines.filter(Boolean).map((line) => {
    const values = line.split("\t");
    const row = {};
    headers.forEach((header, index) => {
      const value = values[index] ?? "";
      row[header] = NUMERIC_FIELDS.has(header) || /^vina_[a-z]+_affinity$/.test(header)
        ? (value.trim() === "" ? NaN : Number(value)) : value;
    });
    return row;
  });
}

function escapeHTML(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function resultStatus(row, source = "sfct") {
  return String(row?.status ?? row?.[`${source}_status`] ?? "").trim().toLowerCase();
}

function isSuccessfulResult(row) {
  return resultStatus(row, "vina") === "success" || resultStatus(row, "sfct") === "success";
}

function hasUsableScore(row, metric) {
  return resultStatus(row, metric === "vina_affinity" ? "vina" : "sfct") === "success"
    && Number.isFinite(row?.[metric]);
}

function savedSfctPoseId(row) {
  // The saved SFCT index is zero-based; exported pose_id is the actual Vina MODEL number.
  const index = row?.sfct_best_pose;
  return Number.isSafeInteger(index) && index >= 0 && Number.isSafeInteger(index + 1) ? index + 1 : null;
}

function addDerivedScores(row) {
  // origin_score is exported as sfct_vina_score: both inputs describe the saved SFCT pose.
  // Never replace it with the independent Vina minimum, or overwrite the supplied 80% score.
  const valid = resultStatus(row, "sfct") === "success"
    && Number.isFinite(row.sfct_vina_score) && Number.isFinite(row.sfct_score);
  return { ...row, vina_sfct_combined_50: valid ? 0.5 * row.sfct_vina_score + 0.5 * row.sfct_score : NaN };
}

function fmt(value, digits = 3) {
  return Number.isFinite(value) ? value.toFixed(digits) : "—";
}

function geneLabel(row) {
  return (row?.[ORGANISM.rowIdentifier] || `No ${ORGANISM.identifierLabel}`).split(";").filter(Boolean).join(" · ");
}

function annotationValue(annotation, field, fallback = "Not available") {
  return annotation?.[field]?.trim() || fallback;
}

function locusLabel(annotation, fallback = `No ${ORGANISM.identifierLabel}`) {
  const ids = annotationValue(annotation, ORGANISM.annotationIdentifier, "").split(/[;\s]+/).filter(Boolean);
  return ids.length ? ids.join(" · ") : fallback;
}

function geneSymbol(annotation) {
  return annotationValue(annotation, "Gene Names", "").split(/\s+/).filter(Boolean)[0] || "—";
}

function annotationSearchText(protein) {
  const annotation = state.annotations.get(protein);
  return annotation ? [annotation.Entry, annotation["Entry Name"], annotation["Gene Names"], annotation[ORGANISM.annotationIdentifier]].join(" ") : "";
}

function median(values) {
  if (!values.length) return NaN;
  const sorted = [...values].sort((a, b) => a - b);
  const midpoint = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[midpoint] : (sorted[midpoint - 1] + sorted[midpoint]) / 2;
}

function mean(values) {
  return values.length ? values.reduce((sum, value) => sum + value, 0) / values.length : NaN;
}

function standardDeviation(values) {
  if (values.length < 2) return 0;
  const average = mean(values);
  return Math.sqrt(values.reduce((sum, value) => sum + (value - average) ** 2, 0) / values.length);
}

function quantile(values, q) {
  if (!values.length) return NaN;
  const sorted = [...values].sort((a, b) => a - b);
  const position = (sorted.length - 1) * q;
  const base = Math.floor(position);
  const rest = position - base;
  return sorted[base + 1] === undefined ? sorted[base] : sorted[base] + rest * (sorted[base + 1] - sorted[base]);
}

function passesPocketQuality(row, options = state) {
  if (options === null) return true; // Explicit unfiltered GO background only.
  const p2 = options.p2 ?? options.p2rank ?? state.p2rank;
  const plddt = options.plddt ?? state.plddt;
  return Number.isFinite(row.probability) && Number.isFinite(row.mean_pocket_plddt)
    && row.probability >= p2 && row.mean_pocket_plddt >= plddt;
}

function compareProteinScores(a, b, metric) {
  return a[metric] - b[metric] || a.uniprot_id.localeCompare(b.uniprot_id);
}

function proteomePercentile(zeroBasedRank, populationSize) {
  return populationSize <= 1 ? 0 : 100 * zeroBasedRank / (populationSize - 1);
}

function comparePocketScores(a, b, metric) {
  return a[metric] - b[metric]
    || (Number.isFinite(a.rank) ? a.rank : Infinity) - (Number.isFinite(b.rank) ? b.rank : Infinity)
    || String(a.pocket).localeCompare(String(b.pocket), undefined, { numeric: true })
    || String(a.protein).localeCompare(String(b.protein));
}

function rankedPopulation(aa, metric, options = state) {
  const qcKey = options === null ? "unfiltered" : `${options.p2 ?? options.p2rank ?? state.p2rank}|${options.plddt ?? state.plddt}`;
  const cacheKey = `${aa}|${metric}|${qcKey}`;
  const source = state.rawByAA.get(aa);
  const cached = state.rankingCache.get(cacheKey);
  if (cached?.source === source && cached) return cached.result;
  const bestByProtein = new Map();
  for (const row of source || []) {
    if (!hasUsableScore(row, metric) || !passesPocketQuality(row, options)) continue;
    const current = bestByProtein.get(row.uniprot_id);
    if (!current || comparePocketScores(row, current, metric) < 0) bestByProtein.set(row.uniprot_id, row);
  }
  const ranked = [...bestByProtein.values()]
    .sort((a, b) => compareProteinScores(a, b, metric))
    .map((row, index, all) => ({
      ...row,
      aa,
      proteome_rank: index + 1,
      proteome_percentile: proteomePercentile(index, all.length),
    }));
  const result = { rows: ranked, byProtein: new Map(ranked.map((row) => [row.uniprot_id, row])),
    percentile: new Map(ranked.map((row) => [row.uniprot_id, row.proteome_percentile])) };
  if (state.rankingCache.size >= 64) state.rankingCache.delete(state.rankingCache.keys().next().value);
  state.rankingCache.set(cacheKey, { source, result });
  return result;
}

function getRanking(aa = state.aa, metric = state.metric, options = state) {
  return rankedPopulation(aa, metric, options).rows;
}

function getRankedProtein(protein, aa, metric, options = state) {
  return rankedPopulation(aa, metric, options).byProtein.get(protein);
}

function pocketProteomePosition(population, row, metric) {
  if (!row || !hasUsableScore(row, metric) || !population.byProtein.has(row.uniprot_id)) return null;
  let low = 0, high = population.rows.length;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (compareProteinScores(population.rows[middle], row, metric) < 0) low = middle + 1;
    else high = middle;
  }
  // Replace this protein's best score, rather than counting it as another protein.
  const ownBest = population.byProtein.get(row.uniprot_id);
  const betterCount = low - (compareProteinScores(ownBest, row, metric) < 0 ? 1 : 0);
  return { proteome_rank: betterCount + 1,
    proteome_percentile: proteomePercentile(betterCount, population.rows.length) };
}

function addPocketProteomePosition(row, aa, metric, options = state) {
  if (!hasUsableScore(row, metric) || !passesPocketQuality(row, options)) return null;
  const position = pocketProteomePosition(rankedPopulation(aa, metric, options), row, metric);
  if (!position) return null;
  return {
    ...row, aa, ...position,
  };
}

function profilePocketAnchor(protein, pocket, targetAA = state.aa) {
  if (!pocket) return null;
  if (typeof pocket === "object") return pocket.uniprot_id === protein ? pocket : null;
  const pinned = state.profilePocketAnchor;
  if (pinned?.uniprot_id === protein && pinned.pocket === pocket) return pinned;
  return (state.rawByAA.get(targetAA) || []).find((row) => row.uniprot_id === protein && row.pocket === pocket) || null;
}

function getProteinProfile(protein, metric = state.metric, pocket = null, options = state, targetAA = state.aa) {
  const anchor = profilePocketAnchor(protein, pocket, targetAA);
  return AMINO_ACIDS.map(({ code, name }) => {
    const row = pocket
      ? anchor && (state.rawByAA.get(code) || []).find((entry) => samePocketGeometry(anchor, entry)
        && hasUsableScore(entry, metric) && passesPocketQuality(entry, options))
      : getRankedProtein(protein, code, metric, options);
    const positioned = row && pocket ? addPocketProteomePosition(row, code, metric, options) : row;
    return positioned ? { ...positioned, code, name, value: positioned[metric] } : null;
  }).filter(Boolean).sort((a, b) => a.value - b.value).map((entry, index) => ({ ...entry, aa_rank: index + 1 }));
}

function getProfileValueConfig() {
  const config = PROFILE_VALUES[state.profileValue];
  return state.profileNormalize && !config.percentile
    ? { ...config, normalized: true, digits: 2, label: `AA-normalized Z-score · ${METRICS[config.metric].label}`, axisLabel: "AA-normalized Z-score" }
    : config;
}

function syncProfileValueControls() {
  const percentile = Boolean(PROFILE_VALUES[state.profileValue].percentile);
  if (percentile) state.profileNormalize = false;
  $("#profile-value-select").value = state.profileValue;
  $("#profile-normalize").checked = state.profileNormalize;
  $("#profile-normalize").disabled = percentile;
  $("#profile-normalize-control").classList.toggle("is-disabled", percentile);
  $("#profile-normalize-control").title = percentile
    ? "Percentiles already express relative standing; Z-score normalization is unavailable."
    : "Normalize the selected Value type separately for each AA, using its QC-passing proteome distribution.";
}

function getProfileScoreSummary(aa, metric) {
  const population = rankedPopulation(aa, metric);
  // Cache only with this AA/metric/QC population, never with the displayed protein or table filters.
  if (!population.profileScoreSummary) {
    const values = population.rows.map((row) => row[metric]);
    population.profileScoreSummary = {
      count: values.length,
      mean: mean(values),
      // Explicitly detect constant scores, including floating-point rounding of their mean.
      sd: values.length < 2 || values.every((value) => value === values[0]) ? 0 : standardDeviation(values),
      // Last score in the existing rank-based Top 5% tier, not a normal-distribution estimate.
      top5Score: values[Math.floor((values.length - 1) * 5 / 100)] ?? NaN,
    };
  }
  return population.profileScoreSummary;
}

function getProfilePlotData(protein, pocket = null) {
  const config = getProfileValueConfig();
  return getProteinProfile(protein, config.metric, pocket).map((entry) => {
    if (!config.normalized) return { ...entry, plotValue: config.percentile ? entry.proteome_percentile : entry[config.metric] };
    const normalization = getProfileScoreSummary(entry.code, config.metric);
    const usable = normalization.count >= 2 && Number.isFinite(normalization.mean)
      && Number.isFinite(normalization.sd) && normalization.sd > 0;
    return { ...entry, normalization,
      plotValue: usable ? (entry[config.metric] - normalization.mean) / normalization.sd : NaN,
      top5Z: usable ? (normalization.top5Score - normalization.mean) / normalization.sd : NaN };
  });
}

function getProfileNormalizationNote(profile) {
  const config = getProfileValueConfig();
  if (!config.normalized) return "";
  const missing = AMINO_ACIDS.filter(({ code }) => !profile.some((entry) => entry.code === code && Number.isFinite(entry.plotValue)));
  return `${METRICS[config.metric].label} (selected Value type). Z = (score − AA mean) / AA SD, using each AA’s best successful pocket per protein after P2Rank ≥ ${state.p2rank} and pLDDT ≥ ${state.plddt}. Negative is favorable; 0 is typical. Dashed purple segments mark each AA’s actual Top 5% score boundary, not statistical significance. Not a binding probability.${missing.length ? ` Missing/unavailable Z: ${missing.map(({ code }) => code).join(", ")}.` : ""}`;
}

function getProfileReference(profile) {
  if (getProfileValueConfig().normalized) return { value: 0, label: "Z = 0" };
  if (PROFILE_VALUES[state.profileValue].percentile) return { value: 10, label: "Top 10% cutoff" };
  const metric = PROFILE_VALUES[state.profileValue].metric;
  const target = profile.find((entry) => entry.code === state.aa);
  if (!target) return { value: NaN, label: "Target score missing" };
  const value = target.plotValue + METRICS[metric].nearWindow;
  return { value, label: `Competitive-AA threshold ${fmt(value)}` };
}

function getProfileTier(percentile) {
  if (percentile <= 1) return { className: "tier-top1", label: "Top 1%" };
  if (percentile <= 5) return { className: "tier-top5", label: "Top 5%" };
  if (percentile <= 10) return { className: "tier-top10", label: "Top 10%" };
  return { className: "tier-outside", label: "Outside top 10%" };
}

function samePocketGeometry(left, right) {
  if (!left || !right || !left.protein || !right.protein || !String(left.residue_ids || "").trim() || !String(right.residue_ids || "").trim()) return false;
  const residues = (row) => [...new Set(String(row.residue_ids).trim().split(/\s+/))].sort().join(" ");
  return left.uniprot_id === right.uniprot_id && left.protein === right.protein && left.pocket === right.pocket
    && residues(left) === residues(right)
    && ["center_x", "center_y", "center_z"].every((field) =>
      Number.isFinite(left[field]) && Number.isFinite(right[field]) && Math.abs(left[field] - right[field]) <= 1e-9);
}

function getStereoControl(protein, targetAA = state.aa, metric = state.metric, pocket = null, options = state) {
  const dAA = dControlCode(targetAA);
  if (!dAA) return null;
  const lAA = canonicalAACode(targetAA);
  const anchor = profilePocketAnchor(protein, pocket, lAA);
  const lRow = pocket
    ? anchor && (state.rawByAA.get(lAA) || []).find((entry) => samePocketGeometry(anchor, entry) && hasUsableScore(entry, metric) && passesPocketQuality(entry, options))
    : getRankedProtein(protein, lAA, metric, options);
  const dRow = pocket
    ? anchor && (state.rawByAA.get(dAA) || []).find((entry) => samePocketGeometry(anchor, entry) && hasUsableScore(entry, metric) && passesPocketQuality(entry, options))
    : getRankedProtein(protein, dAA, metric, options);
  if (!lRow || !dRow || (pocket && !samePocketGeometry(lRow, dRow))) return null;
  const delta = dRow[metric] - lRow[metric];
  return { lRow, dRow, delta, lPreferred: delta > 0 };
}

function getComparison(protein, targetAA = state.aa, metric = state.metric, pocket = null, options = state) {
  const profile = getProteinProfile(protein, metric, pocket, options, targetAA);
  const target = profile.find((entry) => entry.code === targetAA);
  if (!target) return { target: null, profile, otherCount: profile.length, delta: NaN, z: NaN, otherMedian: NaN, nearCompetitors: NaN, stereoControl: null };
  const others = profile.filter((entry) => entry.code !== targetAA).map((entry) => entry.value);
  const otherMedian = median(others);
  const sd = standardDeviation(others);
  const nearWindow = METRICS[metric].nearWindow;
  return {
    target,
    profile,
    otherCount: others.length,
    otherMedian,
    delta: otherMedian - target.value,
    z: others.length >= 2 && sd > 0 ? (mean(others) - target.value) / sd : NaN,
    nearCompetitors: others.length ? others.filter((value) => value <= target.value + nearWindow).length : NaN,
    stereoControl: getStereoControl(protein, targetAA, metric, pocket, options),
  };
}

function getPockets(protein, aa = state.aa, metric = state.metric, options = state) {
  return (state.rawByAA.get(aa) || [])
    .filter((row) => row.uniprot_id === protein && hasUsableScore(row, metric) && passesPocketQuality(row, options))
    .sort((a, b) => comparePocketScores(a, b, metric));
}

function sortValue(row, key) {
  if (key === "qphi_kT") return poseElectrostaticsForRow(row).value;
  const values = {
    protein: `${row.uniprot_id} ${geneLabel(row)}`,
    pocket: row.pocket,
    score: row[state.metric],
    aa_rank: row.comparison.target?.aa_rank,
    competitors: row.comparison.nearCompetitors,
    percentile: row.proteome_percentile,
    p2rank: row.probability,
    plddt: row.mean_pocket_plddt,
    delta: row.comparison.delta,
    z: row.comparison.z,
    stereo: row.comparison.stereoControl?.delta,
  };
  return values[key];
}

function sortRows(rows) {
  const direction = state.sortDirection === "asc" ? 1 : -1;
  return [...rows].sort((a, b) => {
    const left = sortValue(a, state.sortKey);
    const right = sortValue(b, state.sortKey);
    const leftMissing = left === null || left === undefined || Number.isNaN(left);
    const rightMissing = right === null || right === undefined || Number.isNaN(right);
    if (leftMissing || rightMissing) {
      return leftMissing === rightMissing ? a.proteome_rank - b.proteome_rank : leftMissing ? 1 : -1;
    }
    const comparison = typeof left === "string"
      ? left.localeCompare(String(right), undefined, { numeric: true, sensitivity: "base" })
      : left - right;
    return comparison === 0 ? a.proteome_rank - b.proteome_rank : comparison * direction;
  });
}

function readURLState() {
  const params = new URLSearchParams(location.search);
  const aa = params.get("aa")?.toUpperCase();
  // Preserve old shared links to the Vina option while using its independent source now.
  const metric = params.get("metric") === "sfct_vina_score" ? "vina_affinity" : params.get("metric");
  if (AMINO_ACIDS.some((entry) => entry.code === aa)) state.aa = aa;
  if (METRICS[metric]) state.metric = metric;
  if ([1, 5, 10, 100].includes(Number(params.get("top")))) state.top = Number(params.get("top"));
  if (params.has("p2rank")) state.p2rank = Math.min(1, Math.max(0, Number(params.get("p2rank")) || 0));
  if (params.has("plddt")) state.plddt = Math.min(100, Math.max(0, Number(params.get("plddt")) || 0));
  state.aaRank = Math.min(20, Math.max(1, Number(params.get("aarank")) || 20));
  const requestedDelta = params.has("delta") && params.get("delta") !== "" ? Number(params.get("delta")) : null;
  state.minDelta = Number.isFinite(requestedDelta) ? requestedDelta : null;
  state.maxCompetitors = params.has("competitors")
    ? Math.min(19, Math.max(0, Number(params.get("competitors")) || 0))
    : 19;
  state.requireLPreference = params.get("l_preferred") === "1" || (state.aa === "MET" && params.get("lmet_preferred") === "1");
  if (SORT_KEYS.has(params.get("sort"))) state.sortKey = params.get("sort");
  if (["asc", "desc"].includes(params.get("dir"))) state.sortDirection = params.get("dir");
  if (["best", "all"].includes(params.get("pockets"))) state.pocketMode = params.get("pockets");
  if (params.get("distribution") === "filtered") state.distributionSource = "filtered";
  if (["aa", "score_asc"].includes(params.get("profile_order"))) state.profileOrder = params.get("profile_order");
  const profileValue = params.get("profile_value");
  if (profileValue === "aa_zscore") {
    // Migrate old normalized-profile links using their original Explorer score as the base value.
    state.profileValue = Object.keys(PROFILE_VALUES).find((key) => !PROFILE_VALUES[key].percentile && PROFILE_VALUES[key].metric === state.metric);
    state.profileNormalize = true;
  } else {
    if (PROFILE_VALUES[profileValue]) state.profileValue = profileValue;
    state.profileNormalize = params.get("profile_z") === "1" && !PROFILE_VALUES[state.profileValue].percentile;
  }
  state.search = params.get("q") || "";
}

function updateURL() {
  if (organismLeaving) return;
  const params = new URLSearchParams();
  if (ORGANISM.id !== DEFAULT_ORGANISM_ID) params.set("organism", ORGANISM.id);
  if (state.aa !== "ALA") params.set("aa", state.aa);
  if (state.metric !== "vina_sfct_combined_50") params.set("metric", state.metric);
  if (state.top !== 100) params.set("top", state.top);
  if (state.p2rank !== 0.7) params.set("p2rank", state.p2rank);
  if (state.plddt !== 90) params.set("plddt", state.plddt);
  if (state.aaRank !== 20) params.set("aarank", state.aaRank);
  if (state.minDelta !== null && Number.isFinite(state.minDelta)) params.set("delta", state.minDelta);
  if (state.maxCompetitors !== 19) params.set("competitors", state.maxCompetitors);
  if (hasDControl() && state.requireLPreference) params.set("l_preferred", "1");
  if (state.sortKey !== "percentile") params.set("sort", state.sortKey);
  if (state.sortDirection !== "asc") params.set("dir", state.sortDirection);
  if (state.pocketMode !== "best") params.set("pockets", state.pocketMode);
  if (state.distributionSource === "filtered") params.set("distribution", "filtered");
  if (state.profileOrder !== "score_asc") params.set("profile_order", state.profileOrder);
  if (state.profileValue !== "raw_combined_50") params.set("profile_value", state.profileValue);
  if (getProfileValueConfig().normalized) params.set("profile_z", "1");
  if (state.search) params.set("q", state.search);
  const query = params.toString();
  history.replaceState(null, "", `${location.pathname}${query ? `?${query}` : ""}${location.hash}`);
}

function filterRows() {
  const query = state.search.trim().toLowerCase();
  const rankedProteins = getRanking().map((row) => ({ ...row, comparison: getComparison(row.uniprot_id), isBestPocket: true }));
  const proteinEligible = (row) => (state.top === 100 || row.proteome_percentile <= state.top)
      && (row.comparison.target?.aa_rank ?? 99) <= state.aaRank
      && (state.minDelta === null || row.comparison.delta >= state.minDelta)
      && (state.maxCompetitors === 19 || row.comparison.nearCompetitors <= state.maxCompetitors)
      && (!hasDControl() || !state.requireLPreference || row.comparison.stereoControl?.lPreferred);

  let filtered;
  if (state.pocketMode === "best") {
    filtered = rankedProteins.filter((row) => {
      const searchable = [row.uniprot_id, row[ORGANISM.rowIdentifier], row.protein, row.pocket, annotationSearchText(row.uniprot_id)].join(" ").toLowerCase();
      return proteinEligible(row) && (!query || searchable.includes(query));
    });
  } else {
    const eligibleByProtein = new Map(rankedProteins.filter(proteinEligible).map((row) => [row.uniprot_id, row]));
    filtered = (state.rawByAA.get(state.aa) || []).filter((row) => {
      const searchable = [row.uniprot_id, row[ORGANISM.rowIdentifier], row.protein, row.pocket, annotationSearchText(row.uniprot_id)].join(" ").toLowerCase();
      return eligibleByProtein.has(row.uniprot_id)
        && hasUsableScore(row, state.metric)
        && (!query || searchable.includes(query))
        && passesPocketQuality(row);
    }).map((row) => {
      const best = eligibleByProtein.get(row.uniprot_id);
      return {
        ...row,
        aa: state.aa,
        proteome_rank: best.proteome_rank,
        proteome_percentile: best.proteome_percentile,
        comparison: best.comparison,
        isBestPocket: samePocketGeometry(row, best),
      };
    });
  }
  state.filtered = sortRows(filtered);
  return state.filtered;
}

function populateAASelects() {
  const options = AMINO_ACIDS.map(({ code, name }) => `<option value="${code}">${code} · ${name}</option>`).join("");
  $("#aa-select").innerHTML = options;
  $("#profile-aa-select").innerHTML = options;
  $("#aa-select").value = state.aa;
  $("#profile-aa-select").value = state.aa;
}

function syncStereoControls() {
  $("#stereo-control-filter").hidden = !hasDControl();
  $("#stereo-filter-title").textContent = "L>D";
  const description = `Require predicted L-${state.aa} preference over D-${state.aa}; both scores required. See Methods for the comparison rule.`;
  $("#stereo-control-filter").title = description;
  $("#stereo-preference-filter").setAttribute("aria-label", description);
}

function syncControls() {
  $("#aa-select").value = state.aa;
  $("#profile-aa-select").value = state.aa;
  $("#metric-select").value = state.metric;
  $("#protein-search").value = state.search;
  $("#p2rank-filter").value = state.p2rank;
  $("#plddt-filter").value = state.plddt;
  $("#aa-rank-filter").value = state.aaRank;
  $("#delta-filter").value = state.minDelta ?? "";
  $("#competitor-filter").value = state.maxCompetitors;
  $("#stereo-preference-filter").checked = state.requireLPreference;
  $("#profile-order-select").value = state.profileOrder;
  syncProfileValueControls();
  $("#rows-per-page").value = String(state.pageSize);
  $("#distribution-source").value = state.distributionSource;
  syncStereoControls();
  $$("#top-filter button").forEach((button) => button.classList.toggle("active", Number(button.dataset.top) === state.top));
  $$("#pocket-mode-toggle button").forEach((button) => button.classList.toggle("active", button.dataset.pocketMode === state.pocketMode));
}

function renderActiveQuery() {
  const aaInfo = AMINO_ACIDS.find((entry) => entry.code === state.aa);
  const chips = [
    `${aaInfo.code} · ${aaInfo.name}`,
    METRICS[state.metric].label,
    state.top === 100 ? "All proteome ranks" : `Top ${state.top}%`,
    state.pocketMode === "best" ? "Best pocket per protein" : "All retained pockets",
  ];
  if (state.p2rank) chips.push(`P2Rank ≥ ${state.p2rank.toFixed(2)}`);
  if (state.plddt) chips.push(`pLDDT ≥ ${state.plddt}`);
  if (state.aaRank !== 20) chips.push(`AA rank ≤ ${state.aaRank}`);
  if (state.minDelta !== null) chips.push(`Median Δ ≥ ${state.minDelta}`);
  if (state.maxCompetitors !== 19) chips.push(`Competitive AAs ≤ ${state.maxCompetitors}`);
  if (hasDControl() && state.requireLPreference) chips.push("L>D");
  if (state.search) chips.push(`Search: ${state.search}`);
  $("#active-query").innerHTML = chips.map((chip) => `<span class="query-chip">${escapeHTML(chip)}</span>`).join("");
  syncStereoControls();
}

function renderMetrics(rows) {
  const values = rows.map((row) => row[state.metric]);
  const proteinSet = new Set(rows.map((row) => row.uniprot_id));
  const pockets = state.pocketMode === "all"
    ? rows.length
    : (state.rawByAA.get(state.aa) || []).filter((row) => proteinSet.has(row.uniprot_id) && hasUsableScore(row, state.metric) && passesPocketQuality(row)).length;
  const deltas = [...proteinSet].map((protein) => getComparison(protein).delta).filter(Number.isFinite);
  $("#metric-proteins").textContent = proteinSet.size.toLocaleString();
  $("#metric-pockets").textContent = pockets.toLocaleString();
  $("#metric-median").textContent = fmt(median(values));
  $("#metric-selectivity").textContent = fmt(median(deltas));
  $("#metric-label").textContent = METRICS[state.metric].short;
  $("#metric-protein-note").textContent = state.top === 100 ? "matching filters" : `within top ${state.top}%`;
  $("#metric-pocket-note").textContent = state.pocketMode === "all" ? "visible after pocket filters" : "QC-passing for matching proteins";
}

function renderResults(rows) {
  const body = $("#results-body");
  const totalRows = rows.length;
  const pageCount = Math.max(1, Math.ceil(totalRows / state.pageSize));
  state.page = Math.min(Math.max(1, state.page), pageCount);
  const start = (state.page - 1) * state.pageSize;
  const end = Math.min(start + state.pageSize, totalRows);
  const visibleRows = rows.slice(start, end);
  body.innerHTML = visibleRows.map((row) => {
    const comparison = row.comparison;
    const aaRank = comparison.target?.aa_rank ?? "—";
    const percentile = row.proteome_percentile;
    return `<tr>
      <td>${proteinTableIdentity(row.uniprot_id, row)}</td>
      <td><span class="pocket-tag ${row.isBestPocket ? "" : "secondary-pocket"}">${escapeHTML(row.pocket)}</span>${row.isBestPocket && state.pocketMode === "all" ? '<span class="best-pocket-label">best</span>' : ""}</td>
      <td class="numeric score-value">${fmt(row[state.metric])}</td>
      ${poseElectrostaticsCell(row)}
      <td class="numeric">${aaRank} / ${comparison.profile.length}</td>
      <td class="numeric">${Number.isFinite(comparison.nearCompetitors) ? `${comparison.nearCompetitors} / ${comparison.otherCount}` : "—"}</td>
      <td class="numeric"><span class="percentile-cell"><span>${fmt(percentile, 1)}%</span><i class="percentile-track"><i style="width:${Math.max(4, 100 - percentile)}%"></i></i></span></td>
      <td class="numeric">${fmt(row.probability)}</td>
      <td class="numeric">${fmt(row.mean_pocket_plddt, 1)}</td>
      <td class="numeric ${comparison.delta >= 0 ? "positive" : "negative"}">${fmt(comparison.delta)}</td>
      <td class="numeric ${comparison.z >= 0 ? "positive" : "negative"}">${fmt(comparison.z, 2)}</td>
      <td class="numeric stereo-cell ${comparison.stereoControl?.delta >= 0 ? "positive" : "negative"}" ${hasDControl() ? "" : "hidden"}>${comparison.stereoControl ? fmt(comparison.stereoControl.delta) : "—"}</td>
      <td><button class="open-row" type="button" data-protein="${escapeHTML(row.uniprot_id)}" data-pocket="${escapeHTML(row.pocket)}" aria-label="Open ${escapeHTML(row.uniprot_id)} ${escapeHTML(row.pocket)} profile">→</button></td>
    </tr>`;
  }).join("");
  const rowType = state.pocketMode === "best" ? (totalRows === 1 ? "protein" : "proteins") : (totalRows === 1 ? "pocket row" : "pocket rows");
  $("#row-count").textContent = totalRows ? `${start + 1}–${end} of ${totalRows} ${rowType}` : `0 ${rowType}`;
  $("#results-page-status").textContent = `Page ${state.page} of ${pageCount}`;
  $("#previous-results-page").disabled = state.page === 1;
  $("#next-results-page").disabled = state.page === pageCount;
  $("#results-page-controls").hidden = totalRows <= state.pageSize;
  $("#results-heading").textContent = state.pocketMode === "best" ? "Protein hits" : "Pocket hits";
  $("#stereo-header").hidden = !hasDControl();
  $$('[data-sort-header]').forEach((header) => {
    const active = header.dataset.sortHeader === state.sortKey;
    header.setAttribute("aria-sort", active ? (state.sortDirection === "asc" ? "ascending" : "descending") : "none");
    const indicator = header.querySelector(".sort-button span");
    if (indicator) indicator.textContent = active ? (state.sortDirection === "asc" ? "↑" : "↓") : "↕";
  });
  $("#empty-state").hidden = rows.length !== 0;
  $("#results-table").hidden = rows.length === 0;
}

function getDistributionData(rows = state.filtered) {
  const reference = rankedPopulation(state.aa, state.metric);
  const bestPassingPockets = (candidates) => {
    const bestByProtein = new Map();
    for (const row of candidates) {
      if (!hasUsableScore(row, state.metric) || !passesPocketQuality(row)) continue;
      const current = bestByProtein.get(row.uniprot_id);
      if (!current || comparePocketScores(row, current, state.metric) < 0) bestByProtein.set(row.uniprot_id, row);
    }
    return [...bestByProtein.values()].sort((a, b) => compareProteinScores(a, b, state.metric));
  };
  const ranking = reference.rows;
  const population = state.distributionSource === "filtered" ? bestPassingPockets(rows) : ranking;
  const cutoffs = [1, 5, 10].map((tier) => {
    // Apply the rank-percentile formula to the confidence-qualified population.
    const lastIndex = Math.floor((ranking.length - 1) * tier / 100);
    return { tier, value: ranking[lastIndex]?.[state.metric] ?? NaN };
  });
  return {
    values: population.map((row) => row[state.metric]),
    referenceValues: ranking.map((row) => row[state.metric]),
    percentiles: population.map((row) => reference.percentile.get(row.uniprot_id)),
    cutoffs,
  };
}

function histogramSVG(values, cutoffs = [], referenceValues = [], percentiles = null) {
  if (!values.length) return `<div class="empty-state"><p>No scores to plot.</p></div>`;
  const width = 360;
  const height = 235;
  const margin = { top: 18, right: 12, bottom: 38, left: 44 };
  const extent = [...values, ...referenceValues, ...cutoffs.map((cutoff) => cutoff.value)].filter(Number.isFinite);
  let min = extent.reduce((lowest, value) => Math.min(lowest, value), Infinity);
  let max = extent.reduce((highest, value) => Math.max(highest, value), -Infinity);
  if (min === max) { min -= 0.5; max += 0.5; }
  const binCount = 30;
  const binWidth = (max - min) / binCount;
  // Main call supplies the shared QC-protein ranks, including deterministic tie order.
  const positions = percentiles || (() => {
    const ranks = Array(values.length);
    values.map((value, index) => ({ value, index })).sort((a, b) => a.value - b.value || a.index - b.index)
      .forEach((entry, index) => { ranks[entry.index] = proteomePercentile(index, values.length); });
    return ranks;
  })();
  const tierLabels = { 1: "Top 1%", 5: "Top 1–5%", 10: "Top 5–10%", other: "Outside top 10%" };
  const bins = Array.from({ length: binCount }, (_, index) => ({
    start: min + index * binWidth,
    end: min + (index + 1) * binWidth,
    count: 0,
    tierCounts: { 1: 0, 5: 0, 10: 0, other: 0 },
  }));
  values.forEach((value, index) => {
    const bin = bins[Math.min(binCount - 1, Math.floor((value - min) / binWidth))];
    const tier = [1, 5, 10].find((cutoff) => Number.isFinite(positions[index]) && positions[index] <= cutoff) ?? "other";
    bin.count++;
    bin.tierCounts[tier]++;
  });
  const maxCount = Math.max(...bins.map((bin) => bin.count), 1);
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const x = (value) => margin.left + ((value - min) / (max - min)) * plotWidth;
  const y = (count) => margin.top + plotHeight - (count / maxCount) * plotHeight;
  const ticks = Array.from({ length: 4 }, (_, index) => min + (index / 3) * (max - min));
  const medianValue = median(values);
  return `<svg viewBox="0 0 ${width} ${height}" aria-hidden="true">
    ${[0, .5, 1].map((fraction) => `<line class="grid-line" x1="${margin.left}" x2="${width - margin.right}" y1="${y(fraction * maxCount)}" y2="${y(fraction * maxCount)}"/><text class="axis-label" text-anchor="end" x="${margin.left - 6}" y="${y(fraction * maxCount) + 3}">${fmt(fraction * maxCount, fraction * maxCount % 1 ? 1 : 0)}</text>`).join("")}
    ${bins.map((bin) => {
      const barX = x(bin.start) + .7;
      let cumulative = 0;
      const segments = [1, 5, 10, "other"].map((tier) => {
        const count = bin.tierCounts[tier];
        if (!count) return "";
        const bottom = y(cumulative);
        cumulative += count;
        const top = y(cumulative);
        return `<rect class="hist-bar hist-tier-${tier}" x="${barX}" y="${top}" width="${Math.max(1, plotWidth / binCount - 1.4)}" height="${bottom - top}"><title>${tierLabels[tier]}: ${count} protein${count === 1 ? "" : "s"} · ${fmt(bin.start)} to ${fmt(bin.end)}</title></rect>`;
      }).join("");
      return `<g class="hist-bin"><title>${bin.count} protein${bin.count === 1 ? "" : "s"}: ${fmt(bin.start)} to ${fmt(bin.end)}</title>${segments}</g>`;
    }).join("")}
    <line class="median-line" x1="${x(medianValue)}" x2="${x(medianValue)}" y1="${margin.top}" y2="${margin.top + plotHeight}"/>
    <text class="axis-label" text-anchor="${x(medianValue) > width - 60 ? "end" : "start"}" x="${x(medianValue) + (x(medianValue) > width - 60 ? -4 : 4)}" y="${margin.top + 10}">median</text>
    ${ticks.map((tick) => `<text class="axis-label" text-anchor="middle" x="${x(tick)}" y="${height - 21}">${fmt(tick, 2)}</text>`).join("")}
    <text class="axis-label" text-anchor="middle" x="${width / 2}" y="${height - 3}">${escapeHTML(METRICS[state.metric].label)} score</text>
    <text class="axis-label" text-anchor="middle" transform="translate(10 ${margin.top + plotHeight / 2}) rotate(-90)">Proteins</text>
  </svg>`;
}

function renderDistribution(rows) {
  const { values, referenceValues, cutoffs, percentiles } = getDistributionData(rows);
  $("#distribution-aa").textContent = state.aa;
  $("#distribution-source").value = state.distributionSource;
  $("#distribution-summary").textContent = `${values.length.toLocaleString()} ${values.length === 1 ? "protein" : "proteins"}`;
  $("#histogram").innerHTML = histogramSVG(values, cutoffs, referenceValues, percentiles);
  $("#distribution-cutoffs").innerHTML = cutoffs.map(({ tier, value }) => `<div><span><i class="hist-tier-swatch hist-tier-${tier}"></i>Top ${tier}%</span><strong>${fmt(value)}</strong></div>`).join("");
  $("#distribution-stats").innerHTML = [
    ["Q1", quantile(values, .25)],
    ["Median", median(values)],
    ["Q3", quantile(values, .75)],
  ].map(([label, value]) => `<div><strong>${fmt(value)}</strong><span>${label}</span></div>`).join("");
}

function renderExplorer() {
  void prepareExplorerPoseElectrostatics();
  state.page = 1;
  const rows = filterRows();
  renderActiveQuery();
  renderMetrics(rows);
  renderResults(rows);
  renderDistribution(rows);
  updateURL();
}

function selectProtein(protein, pocket = null) {
  state.selectedProtein = protein;
  const anchor = typeof pocket === "object" ? pocket
    : pocket && (state.rawByAA.get(state.aa) || []).find((row) => row.uniprot_id === protein && row.pocket === pocket);
  state.selectedPocket = anchor?.pocket || null;
  state.profilePocket = anchor?.pocket || null;
  state.profilePocketAnchor = anchor || null;
  switchView("protein");
}

let profileProteinSearchIndex = null;

function findProfileProteins(query) {
  const search = query.trim().toLowerCase();
  if (!search) return [];
  if (!profileProteinSearchIndex) {
    profileProteinSearchIndex = [...state.metadata].map(([protein, row]) => {
      const annotation = state.annotations.get(protein);
      const details = [row?.[ORGANISM.rowIdentifier], annotation?.[ORGANISM.annotationIdentifier], annotation?.["Gene Names"], annotation?.["Entry Name"]].filter(Boolean);
      const aliases = details.flatMap((value) => String(value).toLowerCase().split(/[\s;,|]+/)).filter(Boolean);
      return { protein, label: [...new Set(details)].join(" · "), aliases, searchText: [protein, ...details].join(" ").toLowerCase() };
    });
  }
  return profileProteinSearchIndex.filter((entry) => entry.searchText.includes(search)).map((entry) => ({
    ...entry, priority: entry.protein.toLowerCase() === search ? 0 : entry.aliases.includes(search) ? 1 : entry.protein.toLowerCase().startsWith(search) ? 2 : 3,
  })).sort((a, b) => a.priority - b.priority || a.protein.localeCompare(b.protein));
}

function renderProfileProteinSearch() {
  const query = $("#profile-protein-search").value;
  const matches = findProfileProteins(query);
  $("#profile-protein-options").innerHTML = matches.slice(0, 20).map((entry) => `<option value="${escapeHTML(entry.protein)}" label="${escapeHTML(entry.label)}"></option>`).join("");
  $("#profile-protein-search-status").textContent = !query.trim() ? `Type a UniProt ID, gene symbol or ${ORGANISM.identifierLabel}.`
    : !matches.length ? "No loaded proteins match this search."
    : `${matches.length} matching proteins. ${matches.length > 20 ? "Showing the first 20 suggestions; refine your search for more specific matches." : "Choose a suggestion or press Enter for a unique match."}`;
  return matches;
}

function openProfileProteinSearchResult(entry) {
  if (!entry || !state.metadata.has(entry.protein)) return;
  $("#profile-protein-search").value = "";
  $("#profile-protein-options").innerHTML = "";
  $("#profile-protein-search-status").textContent = `Showing profile for ${entry.protein}.`;
  selectProtein(entry.protein);
}

function bindProfileProteinSearchEvents() {
  $("#profile-protein-search").addEventListener("input", renderProfileProteinSearch);
  $("#profile-protein-search").addEventListener("change", () => {
    // A native suggestion inserts a UniProt accession. Free-text aliases wait for submission.
    const exactID = findProfileProteins($("#profile-protein-search").value).find((entry) => entry.priority === 0);
    if (exactID) openProfileProteinSearchResult(exactID);
  });
  $("#profile-protein-search-form").addEventListener("submit", (event) => {
    event.preventDefault();
    if (!$("#profile-protein-search").value.trim()) return;
    const matches = renderProfileProteinSearch();
    const exactID = matches.find((entry) => entry.priority === 0);
    const exactAliases = matches.filter((entry) => entry.priority === 1);
    const match = exactID || (exactAliases.length === 1 ? exactAliases[0] : matches.length === 1 ? matches[0] : null);
    if (match) openProfileProteinSearchResult(match);
    else {
      showToast(matches.length ? "Several proteins match. Choose a UniProt ID from the suggestions or refine your search." : "No loaded proteins match this search.");
      $("#profile-protein-search").focus();
    }
  });
}

function profileChartSVG(profile) {
  if (!profile.length) return '<div class="empty-state"><p>No successful scores available for this profile.</p></div>';
  const config = getProfileValueConfig();
  const available = profile.filter((entry) => Number.isFinite(entry.plotValue));
  if (!available.length) return '<div class="empty-state"><p>No AA-normalized Z-scores available. Each AA needs at least two proteins with successful QC-passing scores and nonzero score variation.</p></div>';
  const ordered = state.profileOrder === "score_asc"
    ? [...available].sort((a, b) => a.plotValue - b.plotValue)
    : AMINO_ACIDS.map(({ code }) => available.find((entry) => entry.code === code)).filter(Boolean);
  const values = ordered.map((entry) => entry.plotValue);
  const reference = getProfileReference(profile);
  const extent = [...values, reference.value, ...(config.normalized ? ordered.map((entry) => entry.top5Z) : [])].filter(Number.isFinite);
  const width = 820;
  const height = 285;
  const margin = { top: 18, right: 12, bottom: 42, left: 48 };
  const min = Math.min(...extent);
  const max = Math.max(...extent);
  const pad = Math.max((max - min) * .12, .25);
  const yMin = min - pad;
  const yMax = max + pad;
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const x = (index) => margin.left + (index / Math.max(1, ordered.length - 1)) * plotWidth;
  const cutoffHalfWidth = Math.min(14, plotWidth / Math.max(1, ordered.length - 1) * .35);
  const y = (value) => margin.top + ((yMax - value) / (yMax - yMin)) * plotHeight;
  const baseline = margin.top + plotHeight;
  const gridValues = Array.from({ length: 4 }, (_, index) => yMin + (index / 3) * (yMax - yMin));
  const path = ordered.map((entry, index) => `${index ? "L" : "M"}${x(index)},${y(entry.plotValue)}`).join(" ");
  const cutoffLabelIndex = config.normalized ? ordered.map((entry) => Number.isFinite(entry.top5Z)).lastIndexOf(true) : -1;
  const cutoffLabelY = cutoffLabelIndex >= 0 ? y(ordered[cutoffLabelIndex].top5Z) - 6 : NaN;
  // Keep both labels readable when the rightmost Top 5% cutoff is close to Z = 0.
  const referenceLabelX = cutoffLabelIndex > 0 && Math.abs(cutoffLabelY - (y(reference.value) - 6)) < 12
    ? width - margin.right - 52 : width - margin.right;
  return `<svg viewBox="0 0 ${width} ${height}" aria-hidden="true">
    ${config.normalized ? `<desc>${escapeHTML(getProfileNormalizationNote(profile))}</desc>` : ""}
    ${gridValues.map((tick) => `<line class="grid-line" x1="${margin.left}" x2="${width - margin.right}" y1="${y(tick)}" y2="${y(tick)}"/><text class="axis-label" text-anchor="end" x="${margin.left - 8}" y="${y(tick) + 3}">${fmt(tick, config.digits)}</text>`).join("")}
    <path d="${path}" fill="none" stroke="#cbb8da" stroke-width="1.5"/>
    ${config.normalized ? ordered.map((entry, index) => Number.isFinite(entry.top5Z)
      ? `<line class="profile-top5-cutoff" data-aa="${entry.code}" x1="${Math.max(margin.left, x(index) - cutoffHalfWidth)}" x2="${Math.min(width - margin.right, x(index) + cutoffHalfWidth)}" y1="${y(entry.top5Z)}" y2="${y(entry.top5Z)}"><title>${entry.code} actual Top 5% boundary: Z = ${fmt(entry.top5Z, config.digits)} · raw score ${fmt(entry.normalization.top5Score)} · ${entry.normalization.count} reference proteins. Boundary ties follow rank-based tier colors.</title></line>` : "").join("") : ""}
    ${ordered.map((entry, index) => {
      const isTarget = entry.code === state.aa;
      const meetsReference = !isTarget && entry.plotValue <= reference.value;
      const tier = getProfileTier(entry.proteome_percentile);
      const pointClass = `${tier.className}${isTarget ? " target" : ""}`;
      const referenceNote = config.normalized ? (entry.plotValue < 0 ? " · better than AA mean" : "") : (meetsReference ? " · meets cutoff" : "");
      const normalizationNote = entry.normalization ? ` · raw score ${fmt(entry.value)} · AA mean ${fmt(entry.normalization.mean)} · population SD ${fmt(entry.normalization.sd)} · ${entry.normalization.count} reference proteins · actual Top 5% boundary Z = ${fmt(entry.top5Z, config.digits)}` : "";
      return `<line class="profile-stem ${pointClass}" x1="${x(index)}" x2="${x(index)}" y1="${baseline}" y2="${y(entry.plotValue)}"/>
        <circle class="profile-dot ${pointClass}" cx="${x(index)}" cy="${y(entry.plotValue)}" r="${isTarget ? 5 : 3.5}"><title>${entry.code}: ${fmt(entry.plotValue, config.digits)} · ${escapeHTML(config.label)}${normalizationNote} · ${tier.label} (percentile ${fmt(entry.proteome_percentile, 1)}%)${referenceNote}${isTarget ? " · target AA" : ""}</title></circle>
        <text class="profile-label" text-anchor="middle" x="${x(index)}" y="${height - 13}">${entry.code}</text>`;
    }).join("")}
    ${Number.isFinite(reference.value) ? `<line class="threshold-line" x1="${margin.left}" x2="${width - margin.right}" y1="${y(reference.value)}" y2="${y(reference.value)}"><title>${escapeHTML(reference.label)}</title></line>
    <text class="threshold-label" text-anchor="end" x="${referenceLabelX}" y="${y(reference.value) - 6}">${escapeHTML(reference.label)}</text>` : ""}
    ${cutoffLabelIndex >= 0 ? `<text class="profile-top5-label" data-aa="${ordered[cutoffLabelIndex].code}" text-anchor="${cutoffLabelIndex === 0 ? "start" : "end"}" x="${x(cutoffLabelIndex)}" y="${cutoffLabelY}">Top 5%</text>` : ""}
    <text class="axis-label" text-anchor="middle" transform="translate(10 ${height / 2}) rotate(-90)">${escapeHTML(config.axisLabel || config.label)}</text>
  </svg>`;
}

function renderProfileTable(profile) {
  $("#profile-aa-body").innerHTML = profile.map((entry) => `<tr class="${entry.code === state.aa ? "target-row" : ""}">
    <td>${entry.aa_rank}</td>
    <td><span class="aa-code"><i></i>${entry.code} <small>${escapeHTML(entry.name)}</small></span></td>
    <td class="numeric">${fmt(entry.vina_affinity)}</td>
    <td class="numeric">${fmt(entry.sfct_score)}</td>
    <td class="numeric">${fmt(entry.vina_sfct_combined)}</td>
    <td class="numeric">${fmt(entry.vina_sfct_combined_50)}</td>
    <td><span class="pocket-tag">${escapeHTML(entry.pocket)}</span></td>
  </tr>`).join("") + AMINO_ACIDS.filter(({ code }) => !profile.some((entry) => entry.code === code)).map(({ code, name }) => `<tr class="${code === state.aa ? "target-row" : ""}">
    <td>—</td><td><span class="aa-code"><i></i>${code} <small>${escapeHTML(name)}</small></span></td>
    <td colspan="5">Missing - no successful QC-passing score at the required pocket/model</td>
  </tr>`).join("");
}

function loadMolstarAssets() {
  if (window.molstar?.Viewer) return Promise.resolve();
  if (molstarAssetsPromise) return molstarAssetsPromise;
  molstarAssetsPromise = new Promise((resolve, reject) => {
    let stylesheet = document.querySelector("link[data-molstar-styles]");
    if (!stylesheet) {
      stylesheet = document.createElement("link");
      stylesheet.rel = "stylesheet";
      stylesheet.href = `${MOLSTAR_ASSET_ROOT}/molstar.css`;
      stylesheet.dataset.molstarStyles = "true";
      document.head.append(stylesheet);
    }

    const existing = document.querySelector("script[data-molstar-script]");
    if (existing) {
      existing.addEventListener("load", resolve, { once: true });
      existing.addEventListener("error", () => reject(new Error("Mol* could not be loaded")), { once: true });
      return;
    }

    const script = document.createElement("script");
    script.src = `${MOLSTAR_ASSET_ROOT}/molstar.js`;
    script.dataset.molstarScript = "true";
    script.onload = resolve;
    script.onerror = () => {
      script.remove();
      reject(new Error("Mol* could not be loaded"));
    };
    document.head.append(script);
  }).catch((error) => {
    molstarAssetsPromise = null;
    throw error;
  });
  return molstarAssetsPromise;
}

async function getMolstarViewer() {
  if (molstarViewer) return molstarViewer;
  if (!molstarViewerPromise) {
    molstarViewerPromise = loadMolstarAssets().then(() => {
      if (!window.molstar?.Viewer) throw new Error("Mol* viewer is unavailable");
      return window.molstar.Viewer.create("molstar-viewer", {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        viewportBackgroundColor: "#f7f4f8",
      });
    }).then((viewer) => {
      molstarViewer = viewer;
      return viewer;
    }).catch((error) => {
      molstarViewerPromise = null;
      throw error;
    });
  }
  return molstarViewerPromise;
}

function setStructureViewerState(kind, message, detail) {
  const overlay = $("#structure-viewer-overlay");
  overlay.dataset.state = kind;
  overlay.hidden = kind === "ready";
  $("#structure-viewer-message").textContent = message;
  $("#structure-viewer-detail").textContent = detail;
}

function parsePocketResidues(residueIds) {
  return String(residueIds || "").trim().split(/\s+/).map((token) => {
    const match = token.match(/^(.+)_(-?\d+)$/);
    return match ? { chain: match[1], number: Number(match[2]) } : null;
  }).filter((residue) => residue && Number.isFinite(residue.number));
}

async function removeMolstarPocketSticks(viewer) {
  molstarPocketSticks = null;
  if (!molstarPocketComponent) return;
  await viewer.plugin.build().delete(molstarPocketComponent.ref).commit();
  molstarPocketComponent = null;
}

async function focusMolstarPocket(viewer, row, request) {
  await removeMolstarPocketSticks(viewer);
  if (request !== molstarLatestRequest) return;
  await setProteinPocketHighlight(viewer, null);
  if (request !== molstarLatestRequest) return;
  const residues = [...new Map(parsePocketResidues(row.residue_ids)
    .map((residue) => [`${residue.chain}_${residue.number}`, residue])).values()];
  // Pocket color is a persistent representation overpaint, not a transient selection marker.
  await viewer.structureInteractivity({ action: "select" });
  if (request !== molstarLatestRequest) return;
  if (!residues.length) {
    $("#structure-pocket-status").textContent = `${row.pocket}: no residue list for pocket coloring or sticks`;
    return;
  }
  const protein = molstarProteinStructure;
  const structure = protein?.obj?.data;
  const StructureElement = window.molstar?.lib?.structure?.StructureElement;
  if (!structure || !StructureElement) throw new Error("Protein structure selection is unavailable");
  // Select on the recorded receptor alone, never on a ligand or another loaded structure.
  const schema = {
    items: residues.map((residue) => ({ auth_asym_id: residue.chain, auth_seq_id: residue.number })),
  };
  const loci = StructureElement.Loci.fromSchema(structure, schema);
  if (StructureElement.Loci.isEmpty(loci)) {
    $("#structure-pocket-status").textContent = `${row.pocket}: listed residues not found in this model`;
    return;
  }
  await setProteinPocketHighlight(viewer, schema);
  if (request !== molstarLatestRequest) return;
  try {
    molstarPocketComponent = await viewer.plugin.builders.structure.tryCreateComponent(protein, {
      type: { name: "bundle", params: StructureElement.Bundle.fromLoci(loci) },
      label: `${row.pocket} · P2Rank residues`, nullIfEmpty: true,
    }, "atlas-selected-pocket");
    if (request !== molstarLatestRequest) return;
    if (!molstarPocketComponent) throw new Error("Pocket residue component is empty");
    molstarPocketSticks = await viewer.plugin.builders.structure.representation.addRepresentation(molstarPocketComponent, {
      type: "ball-and-stick",
      // Bond cylinders only: pocket sticks remain distinct from ligand atom spheres.
      typeParams: { visuals: ["intra-bond", "inter-bond"], sizeFactor: 0.18, sizeAspectRatio: 1, ignoreHydrogens: true },
      color: "uniform", colorParams: { value: 0x39ff14 },
    }, { initialState: { isHidden: true } });
    if (request !== molstarLatestRequest) return;
    if (!molstarPocketSticks) throw new Error("Pocket sticks representation is unavailable");
    setPocketSticksVisibility(viewer);
    viewer.plugin.managers.camera.focusLoci(loci, { minRadius: 8, extraRadius: 4 });
    $("#structure-pocket-status").textContent = `${row.pocket}: P2Rank pocket residues mapped`;
  } catch (error) {
    await removeMolstarPocketSticks(viewer);
    throw error;
  } finally {
    if (request !== molstarLatestRequest) await removeMolstarPocketSticks(viewer);
  }
}

function dockingModelReference(row) {
  const id = String(row?.protein || "").trim().replace(/\.(?:pdb|cif|bcif)(?:\.gz)?$/i, "");
  const match = id.match(/^AF-([A-Za-z0-9]+(?:-\d+)?)-F([1-9]\d*)-model_v([1-9]\d*)$/);
  if (!match || match[1] !== row.uniprot_id) throw new Error("The result does not identify a supported, versioned AlphaFold docking model.");
  return { id, url: `https://alphafold.ebi.ac.uk/files/${id}.cif` };
}

async function performMolstarUpdate(request) {
  if (request !== molstarLatestRequest) return;
  const model = dockingModelReference(request.row);
  setStructureViewerState("loading", "Loading docking model", `Fetching ${model.id} on demand…`);
  const viewer = await getMolstarViewer();
  if (request !== molstarLatestRequest) return;
  if (molstarLoadedModel !== model.id) {
    await viewer.plugin.clear();
    molstarLigandData = null;
    molstarLigandSphere = null;
    molstarPocketCloudData = null;
    molstarPocketCloudSphere = null;
    molstarPocketComponent = null;
    molstarPocketSticks = null;
    molstarProteinStructure = null;
    resetProteinRepresentations();
    molstarLoadedModel = null;
    const data = await viewer.plugin.builders.data.download({ url: model.url, isBinary: false }, { state: { isGhost: true } });
    const trajectory = await viewer.plugin.builders.structure.parseTrajectory(data, "mmcif");
    const loaded = await viewer.plugin.builders.structure.hierarchy.applyPreset(trajectory, "default", {
      structure: { name: "model", params: {} },
      // Explicit cartoon avoids the size-dependent auto preset's all-atom/surface styles.
      representationPreset: "polymer-cartoon",
      representationPresetParams: { theme: { globalName: "plddt-confidence" }, ignoreHydrogens: true },
    });
    molstarProteinStructure = loaded?.structureProperties || loaded?.structure;
    if (!molstarProteinStructure?.obj?.data) throw new Error("The recorded protein model could not be built");
    registerProteinRepresentations(loaded);
    molstarLoadedModel = model.id;
  }
  if (request !== molstarLatestRequest) return;
  await applyProteinStyle(viewer, request);
  if (request !== molstarLatestRequest) return;
  await removeMolstarLigand(viewer);
  await removeMolstarPocketCloud(viewer);
  if (request !== molstarLatestRequest) return;
  setStructureViewerState("ready", "Structure loaded", "");
  try {
    await focusMolstarPocket(viewer, request.row, request);
    if (request === molstarLatestRequest) $("#structure-pocket-status").textContent += ` · ${model.id}`;
  } catch (error) {
    console.warn("Pocket sticks unavailable", error);
    if (request === molstarLatestRequest) $("#structure-pocket-status").textContent = `${request.row.pocket}: protein loaded; pocket sticks unavailable`;
  }
  if (request === molstarLatestRequest) {
    void prepareMolstarPocketCloud(viewer, request);
    await updateMolstarLigand(viewer, request);
  }
}

function updateMolstarPocket(row) {
  if (!row || state.currentView !== "protein") return;
  syncLigandPoseControls(row);
  const protein = state.selectedProtein;
  $("#alphafold-entry-link").href = `https://alphafold.ebi.ac.uk/entry/${encodeURIComponent(protein)}`;
  $("#structure-pocket-status").textContent = `${row.pocket}: preparing pocket residues`;
  setLigandViewerStatus(state.showLigand ? "loading" : "hidden", state.showLigand ? "Preparing docked ligand…" : "Ligand hidden");
  setPocketCloudStatus(state.showPocketCloud ? "loading" : "hidden", state.showPocketCloud ? "Preparing P2Rank cloud…" : "Pocket cloud hidden");
  // Snapshot the profile's underlying score, including normalized/percentile views.
  const request = { protein, row, aa: state.aa, ligandMetric: PROFILE_VALUES[state.profileValue].metric };
  molstarLatestRequest = request;
  molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
    try {
      await performMolstarUpdate(request);
    } catch (error) {
      if (request !== molstarLatestRequest) return;
      console.warn("AlphaFold structure loading failed", error);
      molstarLoadedModel = null;
      molstarProteinStructure = null;
      resetProteinRepresentations();
      setLigandViewerStatus("missing", "Ligand requires the recorded receptor model");
      setPocketCloudStatus("missing", "Pocket cloud requires the recorded receptor model");
      setStructureViewerState("error", "Exact docking model unavailable", `Could not load ${row.protein || "the model recorded in this result"}. No alternate version was substituted.`);
      $("#structure-pocket-status").textContent = `${row.pocket}: the recorded model is required for pocket sticks`;
    }
  });
}

function renderPocketDetail(row) {
  $(".structure-viewer-block").hidden = !row;
  if (!row) {
    molstarLatestRequest = null;
    syncLigandPoseControls(null);
    setLigandViewerStatus("missing", "No pocket selected");
    setPocketCloudStatus("missing", "No pocket selected");
    molstarUpdateQueue = molstarUpdateQueue.catch(() => {}).then(async () => {
      if (molstarViewer) {
        await removeMolstarLigand(molstarViewer);
        await removeMolstarPocketSticks(molstarViewer);
        await setProteinPocketHighlight(molstarViewer, null);
        await removeMolstarPocketCloud(molstarViewer);
      }
    });
    state.selectedPocket = null;
    $("#pocket-detail-title").textContent = "No successful pocket";
    $("#pocket-detail-rank").textContent = "—";
    $("#pocket-detail-list").innerHTML = "";
    $("#pocket-residues").textContent = "—";
    return;
  }
  state.selectedPocket = row.pocket;
  $("#pocket-detail-title").textContent = row.pocket;
  $("#pocket-detail-rank").textContent = `P2Rank #${row.rank}`;
  $("#pocket-detail-list").innerHTML = [
    ["Probability", fmt(row.probability)],
    ["P2Rank score", fmt(row.score, 2)],
    ["Mean pLDDT", fmt(row.mean_pocket_plddt, 1)],
    ["Minimum pLDDT", fmt(row.min_pocket_plddt, 1)],
    ["Residues", row.n_pocket_residues],
    ["Surface atoms", row.surf_atoms],
    ["Saved SFCT/Combined pose", savedSfctPoseId(row) === null ? "—" : `Vina MODEL ${savedSfctPoseId(row)}`],
    ["Status", resultStatus(row)],
  ].map(([label, value]) => `<div><dt>${escapeHTML(label)}</dt><dd>${escapeHTML(value)}</dd></div>`).join("");
  $("#pocket-residues").textContent = row.residue_ids || "No residue list";
  $$("#pockets-body tr").forEach((tr) => tr.classList.toggle("target-row", tr.dataset.pocket === row.pocket));
  updateMolstarPocket(row);
}

function renderPockets(protein, target) {
  const pockets = getPockets(protein);
  const bestPocket = getRanking(state.aa, state.metric).find((row) => row.uniprot_id === protein)?.pocket;
  $("#pockets-heading").textContent = `${state.aa} retained pockets`;
  $("#pocket-count").textContent = `${pockets.length} ${pockets.length === 1 ? "pocket" : "pockets"}`;
  $("#pockets-body").innerHTML = pockets.map((row) => `<tr data-pocket="${escapeHTML(row.pocket)}">
    <td><span class="pocket-name-cell"><span class="pocket-tag">${escapeHTML(row.pocket)}</span>${row.pocket === bestPocket ? '<span class="pocket-best-badge">Best</span>' : ""}</span></td>
    <td class="numeric">#${row.rank}</td>
    <td class="numeric">${fmt(row.probability)}</td>
    <td class="numeric">${fmt(row.mean_pocket_plddt, 1)}</td>
    <td class="numeric">${fmt(row.vina_affinity)}</td>
    <td class="numeric">${fmt(row.sfct_score)}</td>
    <td class="numeric score-value">${fmt(row.vina_sfct_combined)}</td>
    <td class="numeric score-value">${fmt(row.vina_sfct_combined_50)}</td>
    <td>${fmt(row.center_x, 1)}, ${fmt(row.center_y, 1)}, ${fmt(row.center_z, 1)}</td>
    <td><button type="button" class="pocket-select-button" data-pocket="${escapeHTML(row.pocket)}">Inspect</button></td>
  </tr>`).join("");
  if (!pockets.length) $("#pockets-body").innerHTML = '<tr><td colspan="10">Missing — no successful QC-passing pocket results for this AA and score.</td></tr>';
  const selected = state.profilePocket
    ? pockets.find((row) => samePocketGeometry(profilePocketAnchor(protein, state.profilePocket), row))
    : target || pockets[0];
  renderPocketDetail(selected);
}

function renderProtein() {
  if (!state.selectedProtein) {
    const fallback = state.filtered[0] || getRanking()[0];
    if (fallback) state.selectedProtein = fallback.uniprot_id;
  }
  const protein = state.selectedProtein;
  if (!protein) return;
  const profilePocket = state.profilePocket;
  $("#profile-qc-note").textContent = `P2Rank ≥ ${state.p2rank} · pocket mean pLDDT ≥ ${state.plddt} (shared with Explorer). QC applies before selecting and ranking pockets for every AA. ${profilePocket ? "Inspected comparisons require the same model, residues and pocket center; mismatched, failed or QC-excluded results are missing." : "Each AA uses its best QC-passing pocket."}`;
  const comparison = getComparison(protein, state.aa, state.metric, profilePocket);
  const plotProfile = getProfilePlotData(protein, profilePocket);
  const plotReference = getProfileReference(plotProfile);
  const target = comparison.target;
  const meta = target || state.profilePocketAnchor || state.metadata.get(protein);
  const annotation = state.annotations.get(protein);
  $("#protein-title").innerHTML = `<a class="protein-uniprot-link" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(protein)}/entry" target="_blank" rel="noopener noreferrer" title="Open ${escapeHTML(protein)} on UniProt">${escapeHTML(protein)}</a>`;
  $("#protein-subtitle").textContent = `${geneLabel(meta)} · ${meta?.protein || "AlphaFold model"}`;
  $("#profile-entry-name").textContent = annotationValue(annotation, "Entry Name", "No UniProt annotation");
  $("#profile-gene-names").textContent = annotationValue(annotation, "Gene Names", "Gene names not available");
  $("#profile-araport").textContent = locusLabel(annotation, geneLabel(meta));
  const geneInfoButton = $("#profile-gene-info");
  geneInfoButton.dataset.geneInfo = protein;
  geneInfoButton.dataset.geneIds = proteinGeneIds(protein, meta).join(";");
  geneInfoButton.setAttribute("aria-label", `Gene descriptions for ${protein}`);
  geneInfoButton.hidden = false;
  $("#profile-aa-select").value = state.aa;
  $("#profile-score").textContent = fmt(target?.value);
  $("#profile-score-label").textContent = METRICS[state.metric].label;
  $("#profile-rank").textContent = target ? `#${target.aa_rank}` : "—";
  $("#profile-rank-note").textContent = target ? `${Number.isFinite(comparison.nearCompetitors) ? comparison.nearCompetitors : "—"} competitive AAs · ${comparison.profile.length} successful AAs` : `Target missing · ${comparison.profile.length} successful AAs`;
  $("#profile-delta").textContent = fmt(comparison.delta);
  $("#profile-delta-note").textContent = `${comparison.otherCount} successful other AAs`;
  $("#profile-z").textContent = fmt(comparison.z, 2);
  $("#profile-pocket").textContent = target?.pocket || "—";
  $("#profile-pocket-label").textContent = profilePocket ? "Inspected pocket" : "Best pocket";
  $("#profile-confidence").textContent = target ? `P2Rank ${fmt(target.probability)}` : "P2Rank —";
  $("#profile-comparison-heading").textContent = profilePocket ? `20-AA profile · ${profilePocket}` : "Full 20-AA profile";
  $("#profile-pocket-reset").hidden = !profilePocket;
  syncProfileValueControls();
  $("#profile-threshold-key").innerHTML = `<i></i>${escapeHTML(plotReference.label)}`;
  $("#profile-normalization-note").hidden = !getProfileValueConfig().normalized;
  $("#profile-normalization-note").textContent = getProfileNormalizationNote(plotProfile);
  $("#profile-top5-key").hidden = !getProfileValueConfig().normalized;
  const showStereo = Boolean(dControlCode(state.aa));
  const stereo = comparison.stereoControl;
  $("#profile-stereo-card").hidden = !showStereo;
  $(".profile-summary").classList.toggle("has-stereo", showStereo);
  if (showStereo) {
    $("#profile-stereo-label").textContent = `D-${state.aa} control Δ`;
    $("#profile-stereo-delta").textContent = stereo ? fmt(stereo.delta) : "Missing";
    $("#profile-stereo-delta").className = stereo ? (stereo.delta > 0 ? "positive" : stereo.delta < 0 ? "negative" : "") : "";
    $("#profile-stereo-note").textContent = !stereo ? "Successful comparable L and D scores required"
      : stereo.delta === 0 ? "Equal L/D scores" : `${stereo.lPreferred ? "L" : "D"}-${state.aa} predicted better`;
    $("#profile-stereo-card").title = stereo
      ? `L-${state.aa}: ${fmt(stereo.lRow[state.metric])} (${stereo.lRow.pocket}); D-${state.aa}: ${fmt(stereo.dRow[state.metric])} (${stereo.dRow.pocket}). ${profilePocket ? "Same inspected pocket." : "Independent best successful pockets."}`
      : "Missing scores are not evidence of L preference.";
  }
  $("#profile-chart").innerHTML = profileChartSVG(plotProfile);
  renderProfileTable(comparison.profile);
  renderPockets(protein, target);
}

function rankValues(values) {
  const indexed = values.map((value, index) => ({ value, index })).sort((a, b) => a.value - b.value);
  const ranks = Array(values.length);
  let start = 0;
  while (start < indexed.length) {
    let end = start;
    while (end + 1 < indexed.length && indexed[end + 1].value === indexed[start].value) end++;
    const averageRank = (start + end + 2) / 2;
    for (let index = start; index <= end; index++) ranks[indexed[index].index] = averageRank;
    start = end + 1;
  }
  return ranks;
}

function pearson(xs, ys) {
  if (xs.length < 2) return NaN;
  const mx = mean(xs);
  const my = mean(ys);
  let numerator = 0;
  let dx = 0;
  let dy = 0;
  for (let index = 0; index < xs.length; index++) {
    const a = xs[index] - mx;
    const b = ys[index] - my;
    numerator += a * b;
    dx += a * a;
    dy += b * b;
  }
  return dx && dy ? numerator / Math.sqrt(dx * dy) : NaN;
}

function spearman(aaOne, aaTwo, metric) {
  const mapOne = new Map(getRanking(aaOne, metric).map((row) => [row.uniprot_id, row[metric]]));
  const pairs = getRanking(aaTwo, metric)
    .filter((row) => mapOne.has(row.uniprot_id))
    .map((row) => [mapOne.get(row.uniprot_id), row[metric]]);
  return pearson(rankValues(pairs.map((pair) => pair[0])), rankValues(pairs.map((pair) => pair[1])));
}

function mixColor(a, b, amount) {
  return a.map((channel, index) => Math.round(channel + (b[index] - channel) * amount));
}

function correlationColor(value) {
  if (!Number.isFinite(value)) return "#ece9e5";
  const low = [215, 101, 77];
  const middle = [240, 233, 222];
  const high = [101, 53, 148];
  const rgb = value < 0 ? mixColor(low, middle, value + 1) : mixColor(middle, high, value);
  return `rgb(${rgb.join(",")})`;
}

function renderMatrix() {
  const metric = $("#matrix-metric-select").value;
  $("#matrix-qc-note").textContent = `P2Rank ≥ ${state.p2rank} · pocket mean pLDDT ≥ ${state.plddt} (shared with Explorer). Each AA uses its best QC-passing pocket per protein.`;
  $("#matrix-title").textContent = `${METRICS[metric].label}-score correlations`;
  const html = [`<span class="matrix-label"></span>`, ...AMINO_ACIDS.map((aa) => `<span class="matrix-label">${aa.code}</span>`)];
  AMINO_ACIDS.forEach((rowAA) => {
    html.push(`<span class="matrix-label row">${rowAA.code}</span>`);
    AMINO_ACIDS.forEach((columnAA) => {
      const correlation = spearman(rowAA.code, columnAA.code, metric);
      const light = !Number.isFinite(correlation) || Math.abs(correlation) < .48 ? "light-text" : "";
      const label = Number.isFinite(correlation) ? `ρ = ${fmt(correlation, 2)}` : "Missing: fewer than two shared successful scores, or no score variation";
      html.push(`<button type="button" class="matrix-cell ${light}" style="background:${correlationColor(correlation)}" data-matrix-aa="${rowAA.code}" title="${rowAA.code} × ${columnAA.code}: ${label}">${fmt(correlation, 2)}</button>`);
    });
  });
  $("#correlation-matrix").innerHTML = html.join("");
}

function switchView(view) {
  state.currentView = view;
  $$("[data-view-panel]").forEach((panel) => {
    const active = panel.dataset.viewPanel === view;
    panel.hidden = !active;
    panel.classList.toggle("active", active);
  });
  $$("[data-view]").forEach((button) => button.classList.toggle("active", button.dataset.view === view));
  if (view === "protein") renderProtein();
  if (view === "matrix") renderMatrix();
  if (["compare", "overlap", "statistics"].includes(view)) renderAnalysisView(view);
  if (view === "go") renderGOAnalysis();
  if (view === "control-qc") renderControlQC();
  location.hash = view;
  window.scrollTo({ top: 0, behavior: "smooth" });
}

function downloadText(filename, text, mime = "text/tab-separated-values") {
  const blob = new Blob([text], { type: `${mime};charset=utf-8` });
  downloadBlob(filename, blob);
}

function downloadBlob(filename, blob) {
  if (organismLeaving) return;
  filename = `${ORGANISM.id}_${filename}`;
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  URL.revokeObjectURL(url);
  showToast(`Prepared ${filename}`);
}

function rowsToDelimited(rows, delimiter) {
  const headers = ["uniprot_id", "gene_symbol", ORGANISM.rowIdentifier, "amino_acid", "pocket", "is_best_pocket", "score_metric", "score", "proteome_rank", "proteome_percentile", "aa_rank_within_protein", "near_competing_aas", "p2rank_probability", "mean_pocket_plddt", "delta_to_median_other_available_aas", "selectivity_z", "d_minus_l_delta", "vina_score", "sfct_score", "combined_score", "combined_50_score", "available_aa_count", "available_other_aa_count", "d_control_aa", "d_control_score", "d_control_pocket", "qphi_kT", "qphi_status", "qphi_vina_pose", "qphi_error"];
  const quote = (value) => {
    if (typeof value === "number" && !Number.isFinite(value)) return "";
    const string = String(value ?? "");
    if (string.includes(delimiter) || /["\r\n]/.test(string)) return `"${string.replaceAll('"', '""')}"`;
    return string;
  };
  const data = rows.map((row) => {
    const comparison = row.comparison || getComparison(row.uniprot_id);
    const qphi = poseElectrostaticsForRow(row);
    return [row.uniprot_id, geneSymbol(state.annotations.get(row.uniprot_id)), row[ORGANISM.rowIdentifier], state.aa, row.pocket, row.isBestPocket, state.metric, row[state.metric], row.proteome_rank, row.proteome_percentile, comparison.target?.aa_rank, comparison.nearCompetitors, row.probability, row.mean_pocket_plddt, comparison.delta, comparison.z, comparison.stereoControl?.delta ?? "", row.vina_affinity, row.sfct_score, row.vina_sfct_combined, row.vina_sfct_combined_50, comparison.profile.length, comparison.otherCount, dControlCode(state.aa) || "", comparison.stereoControl?.dRow[state.metric] ?? "", comparison.stereoControl?.dRow.pocket ?? "", qphi.value, qphi.status, qphi.pose, qphi.error].map(quote).join(delimiter);
  });
  return [headers.join(delimiter), ...data].join("\n");
}

async function downloadFiltered(delimiter) {
  const aa = state.aa, sourceRows = state.rawByAA.get(aa);
  const poseEntry = await prepareExplorerPoseElectrostatics();
  if (state.aa !== aa || state.rawByAA.get(aa) !== sourceRows || poseElectrostaticsTables.get(poseEntry.key) !== poseEntry) {
    showToast("The target AA changed while preparing the download. Please download again.");
    return;
  }
  const extension = delimiter === "\t" ? "tsv" : "csv";
  downloadText(`aa_atlas_${state.aa.toLowerCase()}_${state.pocketMode}_pockets_filtered.${extension}`, rowsToDelimited(state.filtered, delimiter), delimiter === "\t" ? "text/tab-separated-values" : "text/csv");
}

function profileTableToDelimited(profile, delimiter) {
  const headers = ["Rank", "AA", "Vina", "SFCT", "Combined 80%", "Combined 50%", "Pocket"];
  const rows = profile.map((entry) => [entry.aa_rank, `${entry.code} ${entry.name}`, entry.vina_affinity, entry.sfct_score, entry.vina_sfct_combined, entry.vina_sfct_combined_50, entry.pocket]);
  // Match the table: ranked successful results first, then missing canonical AAs.
  for (const { code, name } of AMINO_ACIDS) {
    if (!profile.some((entry) => entry.code === code)) rows.push([null, `${code} ${name}`, null, null, null, null, null]);
  }
  const quote = (value) => {
    if (typeof value === "number") return Number.isFinite(value) ? String(value) : "";
    const text = String(value ?? "");
    const safe = /^[=+@-]/.test(text) ? `'${text}` : text;
    return safe.includes(delimiter) || /["\r\n]/.test(safe) ? `"${safe.replaceAll('"', '""')}"` : safe;
  };
  return [headers, ...rows].map((row) => row.map(quote).join(delimiter)).join("\n");
}

function downloadProfileTable(delimiter) {
  if (!state.selectedProtein) return;
  const profile = getProteinProfile(state.selectedProtein, state.metric, state.profilePocket);
  const extension = delimiter === "\t" ? "tsv" : "csv";
  const basis = state.profilePocket || "best_pockets";
  const filename = `${state.selectedProtein}_${basis}_${METRICS[state.metric].short.toLowerCase()}_20aa_table.${extension}`;
  downloadText(filename, profileTableToDelimited(profile, delimiter), delimiter === "\t" ? "text/tab-separated-values" : "text/csv");
}

function downloadProfile() {
  if (!state.selectedProtein) return;
  const profile = getProteinProfile(state.selectedProtein, state.metric, state.profilePocket);
  const headers = ["uniprot_id", ORGANISM.rowIdentifier, "amino_acid", "aa_name", "aa_rank_within_protein", "pocket", "vina_score", "sfct_score", "combined_score", "combined_50_score"];
  const rows = profile.map((row) => [state.selectedProtein, row[ORGANISM.rowIdentifier], row.code, row.name, row.aa_rank, row.pocket, row.vina_affinity, row.sfct_score, row.vina_sfct_combined, row.vina_sfct_combined_50].map((value) => typeof value === "number" && !Number.isFinite(value) ? "" : value).join("\t"));
  const basis = state.profilePocket ? `_${state.profilePocket}` : "_best_pockets";
  downloadText(`${state.selectedProtein}${basis}_20aa_profile.tsv`, [headers.join("\t"), ...rows].join("\n"));
}

function buildProfilePlotExportSVG() {
  const source = $("#profile-chart svg");
  if (!source || !state.selectedProtein) return null;

  const namespace = "http://www.w3.org/2000/svg";
  const exportSVG = document.createElementNS(namespace, "svg");
  exportSVG.setAttribute("xmlns", namespace);
  exportSVG.setAttribute("width", "820");
  exportSVG.setAttribute("height", "390");
  exportSVG.setAttribute("viewBox", "0 0 820 390");

  const background = document.createElementNS(namespace, "rect");
  background.setAttribute("width", "820");
  background.setAttribute("height", "390");
  background.setAttribute("fill", "#fffdf9");
  exportSVG.append(background);

  const style = document.createElementNS(namespace, "style");
  style.textContent = `
    svg { font-family: "Segoe UI", Arial, sans-serif; }
    .export-title { fill: #241d2e; font-family: Georgia, serif; font-size: 18px; font-weight: 600; }
    .export-subtitle { fill: #746d7e; font-size: 11px; }
    .export-legend { fill: #746d7e; font-size: 10px; }
    .axis-label { fill: #746d7e; font-size: 9px; }
    .grid-line { stroke: #e9e4de; stroke-width: 1; }
    .threshold-line { stroke: #257b96; stroke-width: 1.5; stroke-dasharray: 6 5; }
    .threshold-label { fill: #257b96; font-size: 8px; font-weight: 800; }
    .profile-top5-cutoff { stroke: #8e65a5; stroke-width: 2; stroke-dasharray: 4 3; }
    .profile-top5-label { fill: #8e65a5; font-size: 8px; font-weight: 800; paint-order: stroke; stroke: #fffdf9; stroke-width: 3; stroke-linejoin: round; }
    .profile-label { fill: #746d7e; font-size: 8px; font-weight: 800; }
    .profile-dot { stroke-width: 2; }
    .profile-dot.tier-top1 { fill: #5f2d7a; stroke: #4e2168; }
    .profile-dot.tier-top5 { fill: #9b63b7; stroke: #7e469c; }
    .profile-dot.tier-top10 { fill: #d8c2e4; stroke: #aa7fc0; }
    .profile-dot.tier-outside { fill: #fffdf9; stroke: #bdb4c1; }
    .profile-dot.target { stroke: #ee6b3b; stroke-width: 3; }
    .profile-stem { stroke-width: 2; }
    .profile-stem.tier-top1 { stroke: #5f2d7a; }
    .profile-stem.tier-top5 { stroke: #9b63b7; }
    .profile-stem.tier-top10 { stroke: #c4a7d3; }
    .profile-stem.tier-outside { stroke: #ded8e0; }
    .profile-stem.target { stroke: #ee6b3b; }
  `;
  exportSVG.append(style);

  const title = document.createElementNS(namespace, "text");
  title.setAttribute("x", "24");
  title.setAttribute("y", "27");
  title.setAttribute("class", "export-title");
  title.textContent = `${state.selectedProtein} · ${state.aa} 20-AA profile`;
  exportSVG.append(title);

  const subtitle = document.createElementNS(namespace, "text");
  subtitle.setAttribute("x", "24");
  subtitle.setAttribute("y", "45");
  subtitle.setAttribute("class", "export-subtitle");
  const orderLabel = state.profileOrder === "score_asc" ? "best → worst" : "canonical AA order";
  const pocketLabel = state.profilePocket ? `${state.profilePocket} across all AAs` : "each AA’s best pocket";
  subtitle.textContent = `${getProfileValueConfig().label} · ${orderLabel} · ${pocketLabel}`;
  exportSVG.append(subtitle);

  const chart = source.cloneNode(true);
  chart.setAttribute("x", "0");
  chart.setAttribute("y", "54");
  chart.setAttribute("width", "820");
  chart.setAttribute("height", "285");
  exportSVG.append(chart);

  const legend = [
    [30, "#5f2d7a", "#5f2d7a", "Top 1%"],
    [120, "#9b63b7", "#9b63b7", "Top 5%"],
    [210, "#d8c2e4", "#aa7fc0", "Top 10%"],
    [310, "#fffdf9", "#bdb4c1", "Outside top 10%"],
    [470, "#fffdf9", "#ee6b3b", "Target AA"],
  ];
  legend.forEach(([x, fill, stroke, label], index) => {
    const circle = document.createElementNS(namespace, "circle");
    circle.setAttribute("cx", String(x + 5));
    circle.setAttribute("cy", "363");
    circle.setAttribute("r", "5");
    circle.setAttribute("fill", fill);
    circle.setAttribute("stroke", stroke);
    circle.setAttribute("stroke-width", index === legend.length - 1 ? "2.5" : "1.5");
    exportSVG.append(circle);

    const labelText = document.createElementNS(namespace, "text");
    labelText.setAttribute("x", String(x + 15));
    labelText.setAttribute("y", "366");
    labelText.setAttribute("class", "export-legend");
    labelText.textContent = label;
    exportSVG.append(labelText);
  });

  if (getProfileValueConfig().normalized) {
    const cutoffKey = document.createElementNS(namespace, "line");
    cutoffKey.setAttribute("x1", "590");
    cutoffKey.setAttribute("x2", "612");
    cutoffKey.setAttribute("y1", "363");
    cutoffKey.setAttribute("y2", "363");
    cutoffKey.setAttribute("class", "profile-top5-cutoff");
    exportSVG.append(cutoffKey);
    const cutoffLabel = document.createElementNS(namespace, "text");
    cutoffLabel.setAttribute("x", "620");
    cutoffLabel.setAttribute("y", "366");
    cutoffLabel.setAttribute("class", "export-legend");
    cutoffLabel.textContent = "AA-specific Top 5% cutoff";
    exportSVG.append(cutoffLabel);
  }

  return new XMLSerializer().serializeToString(exportSVG);
}

async function downloadProfilePlot(format) {
  const svgText = buildProfilePlotExportSVG();
  if (!svgText || !state.selectedProtein) return;
  const pocketFilename = state.profilePocket || "best_pockets";
  const config = getProfileValueConfig();
  const valueFilename = config.normalized ? `aa_zscore_${METRICS[config.metric].short}` : state.profileValue;
  const baseFilename = `${state.selectedProtein}_${state.aa}_${valueFilename}_${pocketFilename}_20aa_profile`.replace(/[^a-zA-Z0-9_.-]+/g, "_");

  if (format === "svg") {
    downloadBlob(`${baseFilename}.svg`, new Blob([svgText], { type: "image/svg+xml;charset=utf-8" }));
    return;
  }

  const svgURL = URL.createObjectURL(new Blob([svgText], { type: "image/svg+xml;charset=utf-8" }));
  try {
    const image = new Image();
    await new Promise((resolve, reject) => {
      image.onload = resolve;
      image.onerror = reject;
      image.src = svgURL;
    });

    const width = 820;
    const height = 390;
    const scale = 2;
    const canvas = document.createElement("canvas");
    canvas.width = width * scale;
    canvas.height = height * scale;
    const context = canvas.getContext("2d");
    context.scale(scale, scale);
    context.drawImage(image, 0, 0, width, height);

    const pngBlob = await new Promise((resolve) => canvas.toBlob(resolve, "image/png"));
    if (!pngBlob) throw new Error("PNG export failed");
    downloadBlob(`${baseFilename}.png`, pngBlob);
  } catch (error) {
    console.error(error);
    showToast("Could not export the profile plot");
  } finally {
    URL.revokeObjectURL(svgURL);
  }
}

let toastTimer;
function showToast(message) {
  const toast = $("#toast");
  toast.textContent = message;
  toast.classList.add("show");
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => toast.classList.remove("show"), 2400);
}

function bindEvents() {
  bindProteinStyleEvents();
  bindProteinColorEvents();
  bindPocketSticksEvents();
  bindLigandViewerEvents();
  bindPocketCloudEvents();
  bindGeneDescriptionEvents();
  bindAnalysisEvents();
  bindControlQCEvents();
  bindGOEvents();
  bindProfileProteinSearchEvents();
  $$("[data-view]").forEach((button) => button.addEventListener("click", () => switchView(button.dataset.view)));
  $$("[data-go-view]").forEach((button) => button.addEventListener("click", () => switchView(button.dataset.goView)));
  $("#aa-select").addEventListener("change", (event) => {
    state.aa = event.target.value;
    if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
    state.selectedPocket = null;
    state.profilePocket = null;
    state.profilePocketAnchor = null;
    $("#profile-aa-select").value = state.aa;
    renderExplorer();
    if (state.selectedProtein) renderProtein();
  });
  $("#profile-aa-select").addEventListener("change", (event) => {
    state.aa = event.target.value;
    if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
    state.selectedPocket = null;
    state.profilePocket = null;
    state.profilePocketAnchor = null;
    $("#aa-select").value = state.aa;
    renderExplorer();
    renderProtein();
  });
  $("#profile-value-select").addEventListener("change", (event) => {
    state.profileValue = event.target.value;
    syncProfileValueControls();
    renderProtein();
    updateURL();
  });
  $("#profile-normalize").addEventListener("change", (event) => {
    state.profileNormalize = event.target.checked && !PROFILE_VALUES[state.profileValue].percentile;
    syncProfileValueControls();
    renderProtein();
    updateURL();
  });
  $("#profile-order-select").addEventListener("change", (event) => {
    state.profileOrder = event.target.value;
    renderProtein();
    updateURL();
  });
  $("#metric-select").addEventListener("change", (event) => {
    state.metric = event.target.value;
    state.selectedPocket = null;
    state.profilePocket = null;
    state.profilePocketAnchor = null;
    renderExplorer();
    if (state.selectedProtein) renderProtein();
  });
  $("#matrix-metric-select").addEventListener("change", renderMatrix);
  $("#distribution-source").addEventListener("change", (event) => {
    state.distributionSource = event.target.value === "filtered" ? "filtered" : "all";
    renderDistribution(state.filtered);
    updateURL();
  });
  $("#protein-search").addEventListener("input", (event) => {
    state.search = event.target.value;
    renderExplorer();
  });
  $("#rows-per-page").addEventListener("change", (event) => {
    state.pageSize = Number(event.target.value) || 10;
    state.page = 1;
    renderResults(state.filtered);
  });
  $("#previous-results-page").addEventListener("click", () => {
    if (state.page <= 1) return;
    state.page -= 1;
    renderResults(state.filtered);
  });
  $("#next-results-page").addEventListener("click", () => {
    const pageCount = Math.max(1, Math.ceil(state.filtered.length / state.pageSize));
    if (state.page >= pageCount) return;
    state.page += 1;
    renderResults(state.filtered);
  });
  $("#p2rank-filter").addEventListener("input", (event) => {
    state.p2rank = Math.min(1, Math.max(0, Number(event.target.value) || 0));
    renderExplorer();
  });
  $("#plddt-filter").addEventListener("input", (event) => {
    state.plddt = Math.min(100, Math.max(0, Number(event.target.value) || 0));
    renderExplorer();
  });
  $("#aa-rank-filter").addEventListener("input", (event) => {
    state.aaRank = Math.min(20, Math.max(1, Number(event.target.value) || 20));
    renderExplorer();
  });
  $("#delta-filter").addEventListener("input", (event) => {
    state.minDelta = event.target.value === "" ? null : Number(event.target.value);
    renderExplorer();
  });
  $("#competitor-filter").addEventListener("input", (event) => {
    state.maxCompetitors = event.target.value === "" ? 19 : Math.min(19, Math.max(0, Number(event.target.value) || 0));
    renderExplorer();
  });
  $("#stereo-preference-filter").addEventListener("change", (event) => {
    state.requireLPreference = event.target.checked;
    renderExplorer();
  });
  $("#top-filter").addEventListener("click", (event) => {
    const button = event.target.closest("button[data-top]");
    if (!button) return;
    state.top = Number(button.dataset.top);
    syncControls();
    renderExplorer();
  });
  $("#reset-filters").addEventListener("click", () => {
    Object.assign(state, { aa: "ALA", metric: "vina_sfct_combined_50", top: 100, p2rank: 0.7, plddt: 90, aaRank: 20, minDelta: null, maxCompetitors: 19, requireLPreference: false, search: "", sortKey: "percentile", sortDirection: "asc", pocketMode: "best", profileOrder: "score_asc", profileValue: "raw_combined_50", profileNormalize: false, selectedPocket: null, profilePocket: null, profilePocketAnchor: null });
    syncControls();
    renderExplorer();
  });
  $("#results-body").addEventListener("click", (event) => {
    const button = event.target.closest("[data-protein]");
    if (button) selectProtein(button.dataset.protein, button.dataset.pocket);
  });
  $("#pocket-mode-toggle").addEventListener("click", (event) => {
    const button = event.target.closest("button[data-pocket-mode]");
    if (!button) return;
    state.pocketMode = button.dataset.pocketMode;
    $$("#pocket-mode-toggle button").forEach((option) => option.classList.toggle("active", option === button));
    renderExplorer();
  });
  $("#results-table thead").addEventListener("click", (event) => {
    const button = event.target.closest("button[data-sort]");
    if (!button) return;
    const key = button.dataset.sort;
    if (state.sortKey === key) state.sortDirection = state.sortDirection === "asc" ? "desc" : "asc";
    else {
      state.sortKey = key;
      state.sortDirection = key === "protein" || key === "pocket" || key === "qphi_kT" ? "asc" : "desc";
    }
    renderExplorer();
  });
  $("#pockets-body").addEventListener("click", (event) => {
    const button = event.target.closest("button[data-pocket]");
    if (!button) return;
    const row = getPockets(state.selectedProtein).find((pocket) => pocket.pocket === button.dataset.pocket);
    if (!row) return;
    state.selectedPocket = row.pocket;
    state.profilePocket = row.pocket;
    state.profilePocketAnchor = row;
    renderProtein();
  });
  $("#profile-pocket-reset").addEventListener("click", () => {
    state.selectedPocket = null;
    state.profilePocket = null;
    state.profilePocketAnchor = null;
    renderProtein();
  });
  $("#explorer-qphi-retry").addEventListener("click", () => retryExplorerPoseElectrostatics());
  $("#download-tsv").addEventListener("click", () => {
    downloadFiltered("\t");
    $("#table-download-menu").open = false;
  });
  $("#download-csv").addEventListener("click", () => {
    downloadFiltered(",");
    $("#table-download-menu").open = false;
  });
  document.addEventListener("click", (event) => {
    $$(".download-menu[open]").forEach((menu) => {
      if (!menu.contains(event.target)) menu.open = false;
    });
  });
  $("#download-profile").addEventListener("click", downloadProfile);
  for (const [format, delimiter] of [["tsv", "\t"], ["csv", ","]]) {
    $(`#download-profile-table-${format}`).addEventListener("click", () => {
      downloadProfileTable(delimiter);
      $("#profile-table-download-menu").open = false;
    });
  }
  $("#download-profile-png").addEventListener("click", () => {
    $("#plot-download-menu").open = false;
    downloadProfilePlot("png");
  });
  $("#download-profile-svg").addEventListener("click", () => {
    $("#plot-download-menu").open = false;
    downloadProfilePlot("svg");
  });
  $("#correlation-matrix").addEventListener("click", (event) => {
    const cell = event.target.closest("[data-matrix-aa]");
    if (!cell) return;
    state.aa = cell.dataset.matrixAa;
    state.selectedPocket = null;
    state.profilePocket = null;
    state.profilePocketAnchor = null;
    if (!hasDControl() && state.sortKey === "stereo") Object.assign(state, { sortKey: "percentile", sortDirection: "asc" });
    syncControls();
    renderExplorer();
    showToast(`${state.aa} selected as the target amino acid`);
  });
}

async function loadGzippedResultText(path) {
  const response = await atlasFetch(path, { cache: "no-cache" });
  if (!response.ok) throw new Error(`Could not load ${path} (HTTP ${response.status}). Check that this file is included in the published site`);
  const bytes = new Uint8Array(await response.arrayBuffer());
  let text;
  // Some servers decode gzip as HTTP content encoding; check the bytes to avoid decoding twice.
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    if (typeof DecompressionStream === "undefined") {
      throw new Error("This browser cannot decompress the AA results. Please open the app in an up-to-date browser");
    }
    try {
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      text = await new Response(stream).text();
    } catch {
      throw new Error(`Could not decompress ${path}. Check that it is a complete, valid gzip file`);
    }
  } else {
    text = new TextDecoder().decode(bytes);
  }

  return text;
}

async function loadResultFile(aa, source) {
  const path = `${RESULTS_DIRECTORY}/${source}/${source}_${aa.file}_all_pockets.tsv.gz`;
  const text = await loadGzippedResultText(path);
  const headers = text.slice(0, text.indexOf("\n") === -1 ? text.length : text.indexOf("\n")).replace(/^\uFEFF/, "").trimEnd().split("\t");
  const scoreColumns = source === "vina" ? [`vina_${aa.file}_affinity`]
    : ["sfct_vina_score", "sfct_score", "vina_sfct_combined"];
  const required = ["uniprot_id", "protein", "pocket", "probability", "mean_pocket_plddt", ...scoreColumns];
  const missing = required.filter((column) => !headers.includes(column));
  if (!headers.includes("status") && !headers.includes(`${source}_status`)) missing.push(`status or ${source}_status`);
  if (missing.length) throw new Error(`Invalid TSV in ${path}: missing columns ${missing.join(", ")}`);
  if (new Set(headers).size !== headers.length) throw new Error(`Invalid TSV in ${path}: duplicate column names`);
  return parseTSV(text);
}

function mergeResultRows(aa, sfctRows, vinaRows) {
  const pockets = new Map();
  for (const [source, rows] of [["sfct", sfctRows], ["vina", vinaRows]]) {
    const seen = new Set();
    for (const row of rows) {
      if (![row.uniprot_id, row.protein, row.pocket].every((value) => String(value || "").trim())) {
        throw new Error(`Invalid ${source} result for ${aa.code}: missing UniProt, protein model or pocket ID`);
      }
      const key = `${row.uniprot_id}|${row.pocket}`;
      if (seen.has(key)) throw new Error(`Duplicate ${source} pocket for ${aa.code}: ${key}`);
      seen.add(key);
      let merged = pockets.get(key);
      if (merged) {
        // Never combine scores from different receptor models or different pocket geometries.
        const sameNumber = (a, b) => a === b || (Number.isNaN(a) && Number.isNaN(b))
          || (Number.isFinite(a) && Number.isFinite(b) && Math.abs(a - b) <= 1e-9);
        const numericFields = ["rank", "probability", "mean_pocket_plddt", "center_x", "center_y", "center_z"];
        const residues = (value) => String(value || "").trim().split(/\s+/).sort().join(" ");
        if (merged.protein !== row.protein || numericFields.some((field) => !sameNumber(merged[field], row[field]))
          || residues(merged.residue_ids) !== residues(row.residue_ids)) {
          throw new Error(`Vina/SFCT pocket metadata mismatch for ${aa.code}: ${key}. Use results from the same receptor model and pockets`);
        }
      } else {
        // This is a processed copy. Source files and parsed source measurements stay unchanged.
        const { status, ...metadata } = row;
        merged = { ...metadata, vina_status: "missing", sfct_status: "missing", vina_affinity: NaN,
          sfct_vina_score: NaN, sfct_score: NaN, vina_sfct_combined: NaN };
        pockets.set(key, merged);
      }
      const status = resultStatus(row, source);
      merged[`${source}_status`] = status;
      if (source === "vina") {
        // No fallback to sfct_vina_score: that is the Vina value of a different, combined-selected pose.
        merged.vina_affinity = status === "success" ? row[`vina_${aa.file}_affinity`] : NaN;
      } else {
        for (const field of ["sfct_vina_score", "sfct_score", "vina_sfct_combined"]) {
          merged[field] = status === "success" ? row[field] : NaN;
        }
      }
    }
  }
  // A failure in one scorer must not remove a successful score from the other source.
  return [...pockets.values()].filter(isSuccessfulResult).map(addDerivedScores);
}

async function loadResultRows(aa) {
  const [sfctRows, vinaRows] = await Promise.all([loadResultFile(aa, "sfct"), loadResultFile(aa, "vina")]);
  return mergeResultRows(aa, sfctRows, vinaRows);
}

async function loadUniProtAnnotations() {
  state.annotations.clear();
  const response = await atlasFetch(UNIPROT_ANNOTATION_PATH, { cache: "no-cache" });
  if (!response.ok) {
    state.unavailableDatasets.push("UniProt annotations");
    return;
  }
  const bytes = new Uint8Array(await response.arrayBuffer());
  let text;
  // A server may already decode gzip via Content-Encoding; inspect the actual bytes.
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    if (typeof DecompressionStream === "undefined") throw new Error("Please use an up-to-date browser to read the compressed UniProt annotations.");
    try {
      const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
      text = await new Response(stream).text();
    } catch {
      throw new Error(`Could not decompress ${UNIPROT_ANNOTATION_PATH}. Check that it is a complete, valid gzip file`);
    }
  } else text = new TextDecoder().decode(bytes);
  const rows = parseTSV(text);
  if (!rows.length || !("Entry" in rows[0])) throw new Error(`Invalid TSV in ${UNIPROT_ANNOTATION_PATH}: missing Entry column or annotation rows`);
  rows.forEach((row) => state.annotations.set(row.Entry, row));
}

function updateAtlasLoadingProgress(percent) {
  const progress = $("#atlas-loading-progress");
  const label = $("#atlas-loading-percent");
  // A failed sibling download may finish after the loader was replaced by an error message.
  if (!progress || !label) return;
  const value = Math.round(Math.max(0, Math.min(100, Number.isFinite(percent) ? percent : 0)));
  progress.value = value;
  label.textContent = `${value}%`;
}

async function loadData() {
  updateAtlasLoadingProgress(0);
  if (location.protocol === "file:") {
    throw new Error("Open this app on GitHub Pages or through a local web server. To run locally, start python -m http.server 8765 in the app folder and open http://localhost:8765/");
  }
  state.rawByAA.clear();
  state.rankingCache.clear();
  state.metadata.clear();
  profileProteinSearchIndex = null;
  const compact = await loadCompactData();
  const ligands = compact ? DATA_LIGANDS.filter(({ code }) => compact.byAA.has(code)) : LEGACY_DATA_LIGANDS;
  // Count completed work: bundle discovery, each AA, annotations, and the initial view.
  // This is not byte/time progress. The last step completes only after init renders the app.
  const totalSteps = ligands.length + 3;
  let completedSteps = 1;
  updateAtlasLoadingProgress(100 * completedSteps / totalSteps);
  // Limit concurrent decompression to keep peak memory lower with full-proteome files.
  const batchSize = 2;
  for (let start = 0; start < ligands.length; start += batchSize) {
    const batch = ligands.slice(start, start + batchSize);
    const results = await Promise.all(batch.map(async (aa) => {
      const rows = await (compact ? loadCompactResultRows(aa, compact) : loadResultRows(aa));
      updateAtlasLoadingProgress(100 * ++completedSteps / totalSteps);
      return [aa.code, rows];
    }));
    results.forEach(([aa, rows]) => {
      state.rawByAA.set(aa, rows);
      rows.forEach((row) => {
        if (!state.metadata.has(row.uniprot_id)) state.metadata.set(row.uniprot_id, row);
      });
    });
  }

  await loadUniProtAnnotations();
  updateAtlasLoadingProgress(100 * ++completedSteps / totalSteps);
}

async function init() {
  try {
    initializeOrganismUI();
    readURLState();
    await loadData();
    if (organismLeaving) return;
    populateAASelects();
    syncControls();
    bindEvents();
    renderExplorer();
    const proteinCount = state.metadata.size;
    $("#release-summary").textContent = `${AMINO_ACIDS.length} amino acids · ${proteinCount} proteins`;
    const missing = [...state.unavailableDatasets, ...DATA_LIGANDS.filter(({ code }) => !state.rawByAA.has(code)).map(({ code }) => `${code} scores`)];
    $("#dataset-availability").hidden = !missing.length;
    $("#dataset-availability").textContent = missing.length ? `Unavailable for ${ORGANISM.name}: ${missing.join(", ")}.` : "";
    const requestedView = location.hash.slice(1);
    if (["explorer", "protein", "matrix", "compare", "overlap", "go", "control-qc", "statistics", "methods"].includes(requestedView) && requestedView !== "explorer") switchView(requestedView);
    updateAtlasLoadingProgress(100);
  } catch (error) {
    if (organismLeaving) return;
    console.error(error);
    $("#loading-screen").innerHTML = `<div class="empty-state"><div>!</div><h3>Atlas data could not load</h3><p>${escapeHTML(error.message)}</p></div>`;
    return;
  }
  requestAnimationFrame(() => { if (!organismLeaving) $("#loading-screen").classList.add("hidden"); });
}

init();
