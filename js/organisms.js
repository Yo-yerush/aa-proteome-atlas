// One immutable organism per document. Navigation disposes all analysis and viewer state.
const ORGANISMS = Object.freeze({
  arabidopsis: Object.freeze({
    id: "arabidopsis", name: "Arabidopsis", scientificName: "Arabidopsis thaliana",
    resultsDirectory: "At_results",
    annotations: "annotations/arabidopsis/arabidopsis_uniprot.tsv.gz",
    descriptions: "annotations/arabidopsis/At_custom_description_file.csv.gz",
    go: "annotations/arabidopsis/arabidopsis_uniprot_go.tsv.gz",
    controls: "annotations/arabidopsis/strict_WT_single_protein_AA_controls.tsv",
    rowIdentifier: "tair_id", annotationIdentifier: "Araport", identifierLabel: "TAIR/Araport ID",
    genePattern: /\bAT(?:[1-5]|C|M)G\d{5}(?:\.\d+)?\b/gi,
    normalizeGene: (id) => id.toUpperCase().replace(/\.\d+$/, ""),
    descriptionFields: [["Short_description", "Short description"], ["Gene_description", "Gene description"], ["Computational_description", "Computational description"]],
  }),
  ecoli: Object.freeze({
    id: "ecoli", name: "E. coli", scientificName: "Escherichia coli",
    resultsDirectory: "Ec_results",
    annotations: "annotations/ecoli/ecoli_uniprot.tsv.gz",
    descriptions: "annotations/ecoli/Ec_custom_description_file.csv.gz",
    go: "annotations/ecoli/ecoli_uniprot_go.tsv.gz",
    controls: "annotations/ecoli/strict_WT_single_protein_AA_controls.tsv",
    rowIdentifier: "gene_id", annotationIdentifier: "Gene Names (ordered locus)", identifierLabel: "Locus ID",
    genePattern: /\bb\d{4}\b/gi,
    normalizeGene: (id) => id.toLowerCase(),
    descriptionFields: [["product", "Product"]],
  }),
});

const DEFAULT_ORGANISM_ID = "ecoli";

function resolveOrganism(search = globalThis.location?.search || "") {
  const id = search ? new URLSearchParams(search).get("organism") : null;
  return Object.hasOwn(ORGANISMS, id) ? ORGANISMS[id] : ORGANISMS[DEFAULT_ORGANISM_ID];
}

const ORGANISM = resolveOrganism();
const RESULTS_DIRECTORY = ORGANISM.resultsDirectory;
let organismAbortController = null;
let organismLeaving = false;

function atlasFetch(path, options = {}) {
  return fetch(path, { ...options, signal: organismAbortController?.signal });
}

function organismGeneIds(value) {
  return [...new Set((String(value || "").match(ORGANISM.genePattern) || []).map(ORGANISM.normalizeGene))];
}

function organismURL(id) {
  if (!Object.hasOwn(ORGANISMS, id)) throw new Error("Unknown organism");
  const url = new URL(location.href);
  url.searchParams.set("organism", id);
  // These identify entries in the old organism, unlike shared AA/quality filters.
  for (const key of ["q", "protein", "pocket"]) url.searchParams.delete(key);
  return url.href;
}

function switchOrganism(id) {
  if (id === ORGANISM.id && !organismLeaving) return;
  const url = organismURL(id);
  organismLeaving = true;
  organismAbortController?.abort();
  // Keep the selector usable during loading, including a second rapid switch.
  document.querySelector("main").inert = true;
  document.querySelector(".primary-nav").inert = true;
  document.querySelector("#loading-screen").classList.remove("hidden");
  document.querySelector("#loading-screen").textContent = `Loading ${ORGANISMS[id].name} atlas…`;
  location.assign(url);
}

function initializeOrganismUI() {
  organismAbortController = new AbortController();
  const selector = document.querySelector("#organism-select");
  selector.value = ORGANISM.id;
  selector.addEventListener("change", (event) => switchOrganism(event.target.value));
  for (const element of document.querySelectorAll("[data-organism-name]")) element.textContent = ORGANISM.name;
  for (const element of document.querySelectorAll("[data-organism-scientific-name]")) element.textContent = ORGANISM.scientificName;
  for (const element of document.querySelectorAll("[data-organism-path]")) {
    const key = element.dataset.organismPath;
    element.textContent = key === "results" ? RESULTS_DIRECTORY + (element.dataset.pathSuffix || "") : ORGANISM[key];
  }
  document.querySelector("#protein-search").placeholder = `UniProt, ${ORGANISM.identifierLabel}, gene symbol…`;
  document.querySelector("#profile-protein-search").placeholder = `Protein / gene / ${ORGANISM.identifierLabel}`;
  document.querySelector("#profile-protein-search-status").textContent = `Type a UniProt ID, gene symbol or ${ORGANISM.identifierLabel}.`;
  document.title = `AA Proteome Atlas · ${ORGANISM.name}`;
  document.querySelector(".loading-title").textContent = `Loading ${ORGANISM.name} atlas…`;
  // A restored page must not retain an aborted session from a previous switch.
  window.addEventListener("pageshow", (event) => { if (event.persisted && organismLeaving) location.reload(); });
}
