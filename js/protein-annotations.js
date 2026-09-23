// On-demand local descriptions, shared by every protein-row table.
const GENE_DESCRIPTION_PATH = ORGANISM.descriptions;
const GENE_DESCRIPTION_FIELDS = ORGANISM.descriptionFields;
const GENE_DETAIL_FIELDS = [
  ["note", "Notes"], ["Protein.families", "Protein families"], ["EC_number", "EC number"],
  ["GO.biological.process", "GO biological process"], ["GO.cellular.component", "GO cellular component"],
  ["GO.molecular.function", "GO molecular function"], ["AraCyc.Name", "AraCyc pathway"],
  ["AraCyc.Db", "AraCyc database"], ["EC", "EC number"], ["KEGG_pathway", "KEGG pathway"],
  ["gene_model_type", "Gene model type"],
  ["refseq_id", "RefSeq IDs"], ["PMID", "PubMed IDs"], ["Derives_from", "Derived from"],
];
let geneDescriptionsPromise = null;
let geneDescriptionRequest = 0;
let geneDescriptionOpener = null;

function proteinGeneIds(protein, row = null) {
  const annotation = state.annotations.get(protein);
  return organismGeneIds([row?.[ORGANISM.rowIdentifier], state.metadata.get(protein)?.[ORGANISM.rowIdentifier],
    annotation?.[ORGANISM.annotationIdentifier], annotation?.["Gene Names"]].filter(Boolean).join(";"));
}

function proteinTableIdentity(protein, row = null) {
  const annotation = state.annotations.get(protein);
  const ids = proteinGeneIds(protein, row);
  const symbol = geneSymbol(annotation);
  const annotationLine = [locusLabel(annotation, "") || ids.join(" · "), symbol === "—" ? "" : symbol].filter(Boolean).join(" · ");
  const uniprotURL = `https://www.uniprot.org/uniprotkb/${encodeURIComponent(protein)}/entry`;
  return `<div class="protein-cell"><div class="protein-id-line"><a class="protein-uniprot-link" href="${escapeHTML(uniprotURL)}" target="_blank" rel="noopener noreferrer" title="Open ${escapeHTML(protein)} on UniProt">${escapeHTML(protein)}</a><button class="gene-info-button" type="button" data-gene-info="${escapeHTML(protein)}" data-gene-ids="${escapeHTML(ids.join(";"))}" aria-haspopup="dialog" aria-controls="gene-description-dialog" aria-label="Gene descriptions for ${escapeHTML(protein)}" title="View gene descriptions">i</button></div>${annotationLine ? `<small>${escapeHTML(annotationLine)}</small>` : ""}</div>`;
}

function parseGeneDescriptionCSV(text) {
  const records = new Map();
  let headers = null, values = [], field = "", quoted = false, rowNumber = 0;
  const finishRow = () => {
    values.push(field); field = "";
    if (values.every((value) => !value.trim())) { values = []; return; }
    rowNumber++;
    if (!headers) {
      headers = values.map((value) => value.trim().replace(/^\uFEFF/, ""));
      if (!headers.includes("gene_id") || !GENE_DESCRIPTION_FIELDS.some(([name]) => headers.includes(name))) {
        throw new Error("The annotation CSV must contain gene_id and a description column.");
      }
    } else {
      if (values.length !== headers.length) throw new Error(`Unexpected column count in annotation CSV record ${rowNumber}.`);
      const row = Object.fromEntries(headers.map((header, index) => {
        const value = values[index].trim();
        return [header, /^(?:NA|N\/A|null|nan)$/i.test(value) ? "" : value];
      }));
      for (const id of organismGeneIds(row.gene_id)) {
        if (!records.has(id)) records.set(id, []);
        records.get(id).push(row);
      }
    }
    values = [];
  };
  const csv = text.replace(/^\uFEFF/, "");
  for (let index = 0; index < csv.length; index++) {
    const character = csv[index];
    if (character === '"') {
      if (quoted && csv[index + 1] === '"') { field += '"'; index++; }
      else quoted = !quoted;
    } else if (character === "," && !quoted) { values.push(field); field = ""; }
    else if ((character === "\n" || character === "\r") && !quoted) {
      if (character === "\r" && csv[index + 1] === "\n") index++;
      finishRow();
    } else field += character;
  }
  if (quoted) throw new Error("The annotation CSV contains an unterminated quoted field.");
  if (field.length || values.length) finishRow();
  if (!headers || !records.size) throw new Error(`No ${ORGANISM.name} gene descriptions were found in the annotation CSV.`);
  return records;
}

async function loadGeneDescriptions() {
  if (!geneDescriptionsPromise) {
    geneDescriptionsPromise = (async () => {
      const response = await atlasFetch(GENE_DESCRIPTION_PATH, { cache: "no-cache" });
      if (!response.ok) throw new Error(`Could not load ${GENE_DESCRIPTION_PATH} (HTTP ${response.status}).`);
      const bytes = new Uint8Array(await response.arrayBuffer());
      let text;
      if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
        if (typeof DecompressionStream === "undefined") throw new Error("Please use an up-to-date browser to read the compressed annotations.");
        const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream("gzip"));
        text = await new Response(stream).text();
      } else text = new TextDecoder().decode(bytes);
      return parseGeneDescriptionCSV(text);
    })().catch((error) => { geneDescriptionsPromise = null; throw error; });
  }
  return geneDescriptionsPromise;
}

function geneDescriptionSections(ids, records) {
  return ids.map((id) => {
    const matches = records.get(id);
    if (!matches?.length) return `<section class="gene-description-record"><h3>${escapeHTML(id)}</h3><p class="gene-description-missing">No matching description is available in the annotation file.</p></section>`;
    return matches.map((record) => {
      const descriptions = GENE_DESCRIPTION_FIELDS.filter(([key]) => record[key]);
      const details = GENE_DETAIL_FIELDS.filter(([key]) => record[key]);
      const fields = (items) => items.map(([key, label]) => `<div><dt><strong>${escapeHTML(label)}</strong></dt><dd>${escapeHTML(record[key])}</dd></div>`).join("");
      const seenSymbols = new Set([String(record.Symbol || "").trim().toLowerCase()]);
      const aliases = String(record.old_symbols || record.gene_synonym || "").split(/[\s,;|]+/).filter((symbol) => {
        const key = symbol.toLowerCase();
        if (!key || seenSymbols.has(key)) return false;
        seenSymbols.add(key);
        return true;
      });
      const otherSymbols = aliases.length ? `<span class="gene-description-aliases"><strong>Other symbols:</strong> ${escapeHTML(aliases.join(" "))}</span>` : "";
      return `<section class="gene-description-record"><h3>${escapeHTML(id)}${record.Symbol ? `<span>${escapeHTML(record.Symbol)}</span>` : ""}${otherSymbols}</h3>${descriptions.length ? "" : '<p class="gene-description-missing">No description text is available for this gene.</p>'}${descriptions.length || details.length ? `<dl class="gene-description-fields">${fields([...descriptions, ...details])}</dl>` : ""}</section>`;
    }).join("");
  }).join("");
}

async function openGeneDescriptions(protein, ids, opener) {
  const dialog = $("#gene-description-dialog"), body = $("#gene-description-body");
  const request = ++geneDescriptionRequest;
  geneDescriptionOpener = opener;
  $("#gene-description-title").textContent = `${protein} · ${ORGANISM.name} gene descriptions`;
  $("#gene-description-ids").textContent = ids.length ? `${ORGANISM.identifierLabel}: ${ids.join(" · ")}` : `No mapped ${ORGANISM.identifierLabel}`;
  body.setAttribute("aria-busy", ids.length ? "true" : "false");
  body.innerHTML = `<p class="gene-description-missing">${ids.length ? "Loading gene descriptions…" : `Unavailable — no ${ORGANISM.identifierLabel} is mapped to this protein in the result or UniProt annotation files.`}</p>`;
  if (!dialog.open) dialog.showModal();
  if (!ids.length) return;
  try {
    const records = await loadGeneDescriptions();
    if (request === geneDescriptionRequest && dialog.open) body.innerHTML = geneDescriptionSections(ids, records);
  } catch (error) {
    if (request === geneDescriptionRequest && dialog.open) body.innerHTML = `<p class="gene-description-error">Unavailable — ${escapeHTML(error.message)} Close and reopen this window to retry. Protein results are still available.</p>`;
  } finally {
    if (request === geneDescriptionRequest) body.setAttribute("aria-busy", "false");
  }
}

function bindGeneDescriptionEvents() {
  const dialog = $("#gene-description-dialog");
  document.addEventListener("click", (event) => {
    const button = event.target.closest("[data-gene-info]");
    if (!button) return;
    const ids = organismGeneIds(button.dataset.geneIds);
    openGeneDescriptions(button.dataset.geneInfo, ids.length ? ids : proteinGeneIds(button.dataset.geneInfo), button);
  });
  $("#gene-description-close").addEventListener("click", () => dialog.close());
  dialog.addEventListener("click", (event) => {
    if (event.target !== dialog) return;
    const rect = dialog.getBoundingClientRect();
    if (event.clientX < rect.left || event.clientX > rect.right || event.clientY < rect.top || event.clientY > rect.bottom) dialog.close();
  });
  dialog.addEventListener("close", () => {
    geneDescriptionRequest++;
    if (geneDescriptionOpener?.isConnected) geneDescriptionOpener.focus();
    geneDescriptionOpener = null;
  });
}
