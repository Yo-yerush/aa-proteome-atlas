#!/usr/bin/env python3
"""Build auditable canonical free-amino-acid datasets using official RCSB APIs.

Python >=3.10; standard library only. See README.md for scope and interpretation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import sys
import time
import urllib.error
import urllib.request
import uuid
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

AMINO_ACIDS = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
ORGANISMS = {"Arabidopsis": "3701", "E. coli": "562", "Human": "9606", "Mouse": "10090", "Yeast": "4932"}
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/graphql"
FREE = "HAS_NO_COVALENT_LINKAGE"
VERSION = "1.0.0"

# The raw response retains all requested metadata, all nonpolymer instances,
# and all API neighbors, including ones subsequently excluded from the tables.
ENTRY_QUERY = """query Entries($ids: [String!]!) {
 entries(entry_ids: $ids) {
  rcsb_id exptl { method } rcsb_entry_info { resolution_combined }
  rcsb_accession_info { initial_release_date }
  polymer_entities {
   rcsb_id entity_poly { rcsb_entity_polymer_type }
   rcsb_polymer_entity { pdbx_description pdbx_mutation }
   rcsb_polymer_entity_container_identifiers { entity_id uniprot_ids }
   rcsb_entity_source_organism {
    ncbi_taxonomy_id ncbi_scientific_name beg_seq_num end_seq_num
    taxonomy_lineage { id }
   }
   rcsb_polymer_entity_align {
    reference_database_name reference_database_accession provenance_source
    aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
   }
   polymer_entity_instances {
    rcsb_id rcsb_polymer_entity_instance_container_identifiers {
     asym_id auth_asym_id auth_to_entity_poly_seq_mapping
    }
   }
  }
  nonpolymer_entities {
   rcsb_id pdbx_entity_nonpoly { comp_id }
   nonpolymer_entity_instances {
    rcsb_id rcsb_nonpolymer_entity_instance_container_identifiers {
     asym_id auth_asym_id auth_seq_id comp_id entity_id
    }
    rcsb_nonpolymer_instance_annotation { comp_id type }
    rcsb_target_neighbors {
     alt_id atom_id comp_id distance target_asym_id target_atom_id
     target_auth_seq_id target_comp_id target_entity_id target_is_bound
     target_model_id target_seq_id
    }
   }
  }
 }
}"""

FIELDS = [
    "AA", "PDB_ID", "ligand_instance", "ligand_label_asym_id",
    "ligand_auth_asym_id", "ligand_auth_seq_id", "ligand_entity_id",
    "protein_chain", "protein_auth_chain", "protein_entity_id",
    "UniProt_accession", "protein_name", "organism", "organism_taxonomy_ids",
    "organism_group", "entry_organism_groups", "experimental_method",
    "resolution_A", "all_resolutions_A", "initial_release_date",
    "contact_residue_count", "contact_residues", "contact_residues_json",
    "minimum_contact_distance_A", "contact_model_ids", "contact_source",
    "contact_cutoff_A", "UniProt_mapping_status", "UniProt_contact_residue_count",
    "UniProt_contact_positions", "UniProt_organism", "UniProt_taxonomy_id",
    "UniProt_organism_group", "UniProt_source_status", "protein_mutations", "has_metal_coordination",
    "status", "control_eligible", "raw_response_file",
]
AUDIT_FIELDS = ["AA", "PDB_ID", "search_hit", "ligand_instance", "status", "output_rows", "raw_response_file"]


def now():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def write_tsv(path, rows, fields):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


class ApiClient:
    def __init__(self, raw_dir, timeout=120, retries=5):
        self.raw_dir = Path(raw_dir)
        self.timeout = timeout
        self.retries = retries

    def post(self, url, payload, category):
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        # 128-bit cache key keeps paths usable in deeply nested Windows workspaces.
        # Cache reads also verify the full request, so collisions cannot mix data.
        digest = hashlib.sha256(url.encode() + encoded).hexdigest()[:32]
        path = self.raw_dir / category / (digest + ".json")
        if path.exists():
            saved = json.loads(path.read_text(encoding="utf-8"))
            if saved["url"] != url or saved["request"] != payload:
                raise RuntimeError(f"Invalid cache: {path}")
            return saved["response"], path
        for attempt in range(self.retries):
            request = urllib.request.Request(url, data=encoded, headers={
                "Content-Type": "application/json", "Accept": "application/json",
                "User-Agent": f"canonical-aa-pdb-dataset/{VERSION}",
            })
            try:
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    body = response.read().decode("utf-8")
                    data = json.loads(body) if body.strip() else None
                    status = response.status
                if isinstance(data, dict) and data.get("errors"):
                    # Never cache or silently accept GraphQL partial success.
                    raise RuntimeError("GraphQL errors: " + json.dumps(data["errors"]))
                atomic_json(path, {"url": url, "request": payload, "retrieved_at_utc": now(),
                                   "http_status": status, "response": data})
                return data, path
            except urllib.error.HTTPError as exc:
                detail = exc.read().decode("utf-8", errors="replace")
                if exc.code not in {408, 429, 500, 502, 503, 504} or attempt + 1 == self.retries:
                    raise RuntimeError(f"HTTP {exc.code} from {url}: {detail[:2000]}") from exc
                delay = min(30.0, 2 ** attempt + random.random())
                try:
                    delay = min(60.0, max(delay, float(exc.headers.get("Retry-After", "0"))))
                except ValueError:
                    pass
            except (urllib.error.URLError, TimeoutError, ConnectionError, json.JSONDecodeError):
                if attempt + 1 == self.retries:
                    raise
                delay = min(30.0, 2 ** attempt + random.random())
            print(f"Retrying {category} in {delay:.1f}s", file=sys.stderr, flush=True)
            time.sleep(delay)
        raise AssertionError("unreachable")


def terminal(attribute, value, operator="exact_match"):
    return {"type": "terminal", "service": "text", "parameters": {
        "attribute": attribute, "operator": operator, "value": value}}


def search_query(aa, organisms, start, page_size):
    return {
        "query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "group", "logical_operator": "and", "label": "nested-attribute", "nodes": [
                terminal("rcsb_nonpolymer_instance_annotation.comp_id", aa),
                terminal("rcsb_nonpolymer_instance_annotation.type", FREE),
            ]},
            terminal("rcsb_entity_source_organism.taxonomy_lineage.id", list(organisms.values()), "in"),
            terminal("rcsb_entry_info.polymer_entity_count_protein", 0, "greater"),
        ]},
        "return_type": "non_polymer_entity",
        "request_options": {"paginate": {"start": start, "rows": page_size},
                            "results_content_type": ["experimental"],
                            "sort": [{"sort_by": "rcsb_id", "direction": "asc"}]},
    }


def search_aa(client, aa, organisms, page_size):
    hits, seen, expected = [], set(), None
    while True:
        data, path = client.post(SEARCH_URL, search_query(aa, organisms, len(hits), page_size), "search")
        count = data.get("total_count", 0) if data else 0
        if expected is None:
            expected = count
        elif count != expected:
            raise RuntimeError(f"{aa}: archive changed during pagination ({expected} -> {count}); start a fresh run")
        page = (data or {}).get("result_set", [])
        for hit in page:
            identifier = hit["identifier"]
            if identifier in seen:
                raise RuntimeError(f"{aa}: duplicate paginated hit {identifier}; start a fresh run")
            seen.add(identifier)
            hits.append({"AA": aa, "search_hit": identifier, "score": hit.get("score", ""),
                         "raw_response_file": path.relative_to(client.raw_dir).as_posix()})
        if len(hits) == expected:
            break
        if not page or len(hits) > expected:
            raise RuntimeError(f"{aa}: incomplete pagination: {len(hits)}/{expected}")
    print(f"{aa}: {len(hits)} experimental nonpolymer entity hits", flush=True)
    return hits


def fetch_entries(client, ids):
    data, path = client.post(DATA_URL, {"query": ENTRY_QUERY, "variables": {"ids": ids}}, "data")
    entries = (data or {}).get("data", {}).get("entries") or []
    found = {entry["rcsb_id"] for entry in entries if entry}
    if found != set(ids) or len(entries) != len(ids):
        raise RuntimeError(f"Data API did not return every entry: missing {sorted(set(ids) - found)}")
    return entries, path


def source_groups(source, organisms):
    lineage = {str(x["id"]) for x in source.get("taxonomy_lineage") or []}
    lineage.add(str(source.get("ncbi_taxonomy_id", "")))
    return [name for name, taxid in organisms.items() if taxid in lineage]


def source_overlaps(source, positions):
    begin, end = source.get("beg_seq_num"), source.get("end_seq_num")
    if begin is None and end is None:
        return True
    return any((begin is None or p >= begin) and (end is None or p <= end) for p in positions)


def finite_resolution(values):
    return [float(x) for x in values or [] if isinstance(x, (int, float)) and math.isfinite(x) and x > 0]


def summarize_contacts(neighbors, chain_ids):
    residues = {}
    author_mapping = chain_ids.get("auth_to_entity_poly_seq_mapping") or []
    for n in neighbors:
        pos = n.get("target_seq_id")
        author_id = str(n.get("target_auth_seq_id", ""))
        if isinstance(pos, int) and 1 <= pos <= len(author_mapping) and author_mapping[pos - 1] not in (None, "?", "."):
            author_id = str(author_mapping[pos - 1])
        key = (pos, author_id, n.get("target_comp_id"))
        if key not in residues:
            residues[key] = {"label_seq_id": pos, "auth_seq_id": author_id,
                             "comp_id": n.get("target_comp_id"), "minimum_distance_A": n["distance"],
                             "model_ids": set()}
        residues[key]["minimum_distance_A"] = min(residues[key]["minimum_distance_A"], n["distance"])
        if n.get("target_model_id") is not None:
            residues[key]["model_ids"].add(n["target_model_id"])
    result = sorted(residues.values(), key=lambda r: (r["label_seq_id"] is None, r["label_seq_id"] or 0, r["auth_seq_id"], r["comp_id"] or ""))
    for residue in result:
        residue["model_ids"] = sorted(residue["model_ids"])
    return result


def uniprot_mappings(polymer, positions):
    identifiers = polymer.get("rcsb_polymer_entity_container_identifiers") or {}
    alignments = [a for a in polymer.get("rcsb_polymer_entity_align") or []
                  if a.get("reference_database_name") == "UniProt"]
    accessions = set(identifiers.get("uniprot_ids") or [])
    accessions.update(a["reference_database_accession"] for a in alignments if a.get("reference_database_accession"))
    for accession in sorted(accessions):
        mapped, covered, has_regions = set(), set(), False
        for alignment in alignments:
            if alignment.get("reference_database_accession") != accession:
                continue
            for region in alignment.get("aligned_regions") or []:
                begin, ref, length = (region.get(k) for k in ("entity_beg_seq_id", "ref_beg_seq_id", "length"))
                if not all(isinstance(v, int) for v in (begin, ref, length)):
                    continue
                has_regions = True
                for pos in positions:
                    if begin <= pos < begin + length:
                        covered.add(pos)
                        mapped.add(ref + pos - begin)
        status = ("contact_overlap" if covered else "no_contact_overlap" if has_regions else
                  "single_accession_no_alignment" if len(accessions) == 1 else "multiple_accessions_unresolved")
        yield accession, status, len(covered), ";".join(map(str, sorted(mapped))), covered
    if not accessions:
        yield "", "missing", 0, "", set()


def process_entry(entry, hits, organisms, cutoff, raw_file):
    rows, audit = [], []
    pdb = entry["rcsb_id"]
    methods = sorted({m["method"] for m in entry.get("exptl") or [] if m.get("method")})
    if not methods:
        raise RuntimeError(f"{pdb}: experimental search hit lacks experimental method")
    polymers = [p for p in entry.get("polymer_entities") or []
                if (p.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
    chains = {}
    entry_groups = set()
    for polymer in polymers:
        for source in polymer.get("rcsb_entity_source_organism") or []:
            entry_groups.update(source_groups(source, organisms))
        for chain in polymer.get("polymer_entity_instances") or []:
            ids = chain["rcsb_polymer_entity_instance_container_identifiers"]
            chains[ids["asym_id"]] = (polymer, ids)
    resolutions = finite_resolution((entry.get("rcsb_entry_info") or {}).get("resolution_combined"))
    nonpolymers = {n["rcsb_id"]: n for n in entry.get("nonpolymer_entities") or []}
    for hit in hits:
        aa, hit_id = hit["AA"], hit["search_hit"]
        if hit_id not in nonpolymers:
            raise RuntimeError(f"{hit_id}: Search/Data nonpolymer entity mismatch")
        nonpolymer = nonpolymers[hit_id]
        if (nonpolymer.get("pdbx_entity_nonpoly") or {}).get("comp_id") != aa:
            audit.append(dict(AA=aa, PDB_ID=pdb, search_hit=hit_id, ligand_instance="",
                              status="entity_component_mismatch", output_rows=0, raw_response_file=raw_file))
            continue
        instances = nonpolymer.get("nonpolymer_entity_instances") or []
        if not instances:
            raise RuntimeError(f"{hit_id}: nonpolymer entity has no instances")
        for ligand in instances:
            lid = ligand["rcsb_nonpolymer_entity_instance_container_identifiers"]
            annotations = ligand.get("rcsb_nonpolymer_instance_annotation") or []
            item = dict(AA=aa, PDB_ID=pdb, search_hit=hit_id, ligand_instance=ligand["rcsb_id"],
                        status="", output_rows=0, raw_response_file=raw_file)
            if lid.get("comp_id") != aa or not any(a.get("comp_id") == aa and a.get("type") == FREE for a in annotations):
                item["status"] = "instance_not_annotated_free"
                audit.append(item)
                continue
            # The entry-level search can match the source of RNA in a mixed
            # species complex. Require a target-organism protein in the Data API.
            if not entry_groups:
                item["status"] = "no_target_organism_protein_in_entry"
                audit.append(item)
                continue
            base = {key: "" for key in FIELDS}
            base.update(AA=aa, PDB_ID=pdb, ligand_instance=ligand["rcsb_id"],
                        ligand_label_asym_id=lid["asym_id"], ligand_auth_asym_id=lid.get("auth_asym_id") or "",
                        ligand_auth_seq_id=lid.get("auth_seq_id") or "", ligand_entity_id=lid["entity_id"],
                        entry_organism_groups=";".join(sorted(entry_groups)), experimental_method=";".join(methods),
                        resolution_A=min(resolutions) if resolutions else "",
                        all_resolutions_A=";".join(map(str, resolutions)),
                        initial_release_date=(entry.get("rcsb_accession_info") or {}).get("initial_release_date") or "",
                        contact_source="RCSB Data API rcsb_target_neighbors", contact_cutoff_A=cutoff,
                        has_metal_coordination=int(any(a.get("type") == "HAS_METAL_COORDINATION_LINKAGE" for a in annotations)),
                        control_eligible=0, raw_response_file=raw_file)
            neighbors = ligand.get("rcsb_target_neighbors") or []
            grouped = defaultdict(list)
            for neighbor in neighbors:
                distance = neighbor.get("distance")
                if (neighbor.get("comp_id") == aa and isinstance(distance, (int, float))
                        and 0 <= distance <= cutoff and neighbor.get("target_asym_id") in chains):
                    grouped[neighbor["target_asym_id"]].append(neighbor)
            emitted = 0
            for chain_id, contacts in sorted(grouped.items()):
                polymer, chain_ids = chains[chain_id]
                positions = {n["target_seq_id"] for n in contacts if isinstance(n.get("target_seq_id"), int)}
                sources = polymer.get("rcsb_entity_source_organism") or []
                groups = {g for s in sources if source_overlaps(s, positions) for g in source_groups(s, organisms)}
                if not groups:
                    continue
                residues = summarize_contacts(contacts, chain_ids)
                identifiers = polymer["rcsb_polymer_entity_container_identifiers"]
                description = polymer.get("rcsb_polymer_entity") or {}
                chain_row = dict(base, protein_chain=chain_id, protein_auth_chain=chain_ids.get("auth_asym_id") or "",
                                 protein_entity_id=identifiers["entity_id"], protein_name=description.get("pdbx_description") or "",
                                 organism=";".join(sorted({s["ncbi_scientific_name"] for s in sources if s.get("ncbi_scientific_name")})),
                                 organism_taxonomy_ids=";".join(sorted({str(s["ncbi_taxonomy_id"]) for s in sources if s.get("ncbi_taxonomy_id")})),
                                 organism_group=";".join(sorted(groups)), contact_residue_count=len(residues),
                                 contact_residues=";".join(f'{r["comp_id"]}:{r["auth_seq_id"]}[label={r["label_seq_id"]}]' for r in residues),
                                 contact_residues_json=json.dumps(residues, separators=(",", ":")),
                                 minimum_contact_distance_A=min(n["distance"] for n in contacts),
                                 contact_model_ids=";".join(map(str, sorted({n["target_model_id"] for n in contacts if n.get("target_model_id") is not None}))),
                                 protein_mutations=description.get("pdbx_mutation") or "", status="target_protein_contact")
                for accession, mapping_status, mapped_count, mapped_positions, covered in uniprot_mappings(polymer, positions):
                    # An accession and a source must overlap the SAME contacted
                    # sequence positions, not opposite ends of a fusion protein.
                    accession_groups = {g for s in sources if source_overlaps(s, covered or positions)
                                        for g in source_groups(s, organisms)}
                    row = dict(chain_row, UniProt_accession=accession, UniProt_mapping_status=mapping_status,
                               UniProt_contact_residue_count=mapped_count, UniProt_contact_positions=mapped_positions,
                               UniProt_organism_group=";".join(sorted(accession_groups)),
                               control_eligible=int(bool(accession_groups) and mapping_status in {"contact_overlap", "single_accession_no_alignment"}))
                    rows.append(row)
                    emitted += 1
            if not emitted:
                status = ("no_target_organism_protein_contact" if grouped else
                          "no_protein_contact_within_cutoff" if neighbors else "no_api_neighbors")
                rows.append(dict(base, status=status, contact_residue_count=0, contact_residues_json="[]"))
                emitted = 1
                item["status"] = status
            else:
                item["status"] = "target_protein_contact"
            item["output_rows"] = emitted
            audit.append(item)
    return rows, audit


def fetch_uniprot_sources(client, accessions):
    # Aliases allow batching the Data API's singular uniprot field.
    query = "{ " + " ".join(
        f'u{i}: uniprot(uniprot_id: {json.dumps(accession)}) {{ rcsb_uniprot_protein {{ source_organism {{ scientific_name taxonomy_id }} }} }}'
        for i, accession in enumerate(accessions)) + " }"
    response, _ = client.post(DATA_URL, {"query": query}, "uniprot")
    data = (response or {}).get("data") or {}
    if set(data) != {f"u{i}" for i in range(len(accessions))}:
        raise RuntimeError("Incomplete UniProt source response")
    return {accession: ((data[f"u{i}"] or {}).get("rcsb_uniprot_protein") or {}).get("source_organism")
            for i, accession in enumerate(accessions)}


def validate_uniprot_sources(rows, sources, organisms):
    for row in rows:
        accession = row["UniProt_accession"]
        if not accession:
            continue
        source = sources.get(accession)
        if not source:
            row["UniProt_source_status"] = "unavailable"
            # Do not claim a species-specific control for an unresolved chimera.
            if ";" in row["organism_taxonomy_ids"]:
                row["control_eligible"] = 0
            continue
        name, taxid = source["scientific_name"], str(source.get("taxonomy_id") or "")
        row["UniProt_organism"] = name
        row["UniProt_taxonomy_id"] = taxid
        # UniProt exposes species/strain names but not lineage in this endpoint.
        # Exact species plus a suffix includes strain names; genus scope is explicit.
        roots = {"Arabidopsis": "Arabidopsis" if organisms["Arabidopsis"] == "3701" else "Arabidopsis thaliana",
                 "E. coli": "Escherichia coli", "Human": "Homo sapiens", "Mouse": "Mus musculus",
                 "Yeast": "Saccharomyces cerevisiae"}
        groups = {g for g, root in roots.items() if taxid == organisms[g] or name == root or name.startswith(root + " ")}
        supported = groups.intersection(row["UniProt_organism_group"].split(";"))
        row["UniProt_organism_group"] = ";".join(sorted(supported))
        row["UniProt_source_status"] = "confirmed_target_source" if supported else "outside_contacting_target_source"
        if not supported:
            row["control_eligible"] = 0


def representative_rank(row):
    resolution = float(row["resolution_A"]) if row["resolution_A"] != "" else math.inf
    # Resolution first; prefer more contacting residues only on a resolution tie.
    return (resolution, -int(row["contact_residue_count"]), row["PDB_ID"],
            row["ligand_instance"], row["protein_chain"])


def select_representatives(rows):
    groups = defaultdict(list)
    for row in rows:
        if row["control_eligible"] and row["UniProt_accession"]:
            groups[row["AA"], row["UniProt_accession"]].append(row)
    result = []
    for key, candidates in sorted(groups.items()):
        representative = min(candidates, key=representative_rank)
        result.append(dict(representative, candidate_structure_count=len({r["PDB_ID"] for r in candidates}),
                           candidate_row_count=len(candidates),
                           representative_selection="lowest_available_resolution_then_most_contacts_then_identifiers"))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("."), help="Output directory (default: current directory)")
    parser.add_argument("--resume", type=Path, help="Reuse an existing raw run directory; preserves its original snapshot")
    parser.add_argument("--organisms", nargs="+", choices=["ecoli", "arabidopsis", "human", "mouse", "yeast"], default=["ecoli", "arabidopsis", "human", "mouse"])
    parser.add_argument("--aa", nargs="+", choices=AMINO_ACIDS, default=AMINO_ACIDS)
    parser.add_argument("--arabidopsis-taxid", choices=["3701", "3702"], default="3701",
                        help="3701: genus (default); 3702: A. thaliana only")
    parser.add_argument("--contact-cutoff", type=float, default=5.0,
                        help="Keep API neighbors at or below this distance in angstroms, maximum 5 (default: 5)")
    parser.add_argument("--workers", type=int, default=4, help="Concurrent Data API batches (default: 4)")
    parser.add_argument("--batch-size", type=int, default=10)
    parser.add_argument("--page-size", type=int, default=1000)
    parser.add_argument("--timeout", type=float, default=120)
    args = parser.parse_args(argv)
    if not 0 < args.contact_cutoff <= 5 or not 1 <= args.workers <= 8 or not 1 <= args.batch_size <= 50 or not 1 <= args.page_size <= 10000 or args.timeout <= 0:
        parser.error("Require 0 < cutoff <= 5, 1..8 workers, 1..50 batch size, 1..10000 page size, positive timeout")
    args.out.mkdir(parents=True, exist_ok=True)
    args.out = args.out.resolve()
    organisms = dict(ORGANISMS, Arabidopsis=args.arabidopsis_taxid)
    aliases = {"ecoli": "E. coli", "arabidopsis": "Arabidopsis", "human": "Human", "mouse": "Mouse", "yeast": "Yeast"}
    search_organisms = {aliases[name]: organisms[aliases[name]] for name in args.organisms}
    aa_list = sorted(set(args.aa), key=AMINO_ACIDS.index)
    config = dict(version=VERSION, amino_acids=aa_list, organisms=search_organisms,
                  contact_cutoff_A=args.contact_cutoff, page_size=args.page_size, batch_size=args.batch_size,
                  data_query_sha256=hashlib.sha256(ENTRY_QUERY.encode()).hexdigest())
    raw_dir = (args.resume.resolve() if args.resume else args.out / "raw" /
               (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]))
    manifest_path = raw_dir / "manifest.json"
    if args.resume:
        if not manifest_path.exists():
            parser.error("Resume directory must contain manifest.json")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["configuration"] != config:
            parser.error("Resume configuration differs; use the original AA, taxonomy, cutoff, batch and page options")
    else:
        manifest = {"created_at_utc": now(), "configuration": config, "search_url": SEARCH_URL,
                    "data_url": DATA_URL, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    manifest.update(status="running", last_started_at_utc=now())
    atomic_json(manifest_path, manifest)
    client = ApiClient(raw_dir, timeout=args.timeout)
    print(f"Raw snapshot: {raw_dir}", flush=True)
    try:
        hits = []
        for aa in aa_list:
            hits.extend(search_aa(client, aa, search_organisms, args.page_size))
            write_tsv(raw_dir / "raw_search_hits.tsv", hits, ["AA", "search_hit", "score", "raw_response_file"])
        by_entry = defaultdict(list)
        for hit in hits:
            by_entry[hit["search_hit"].rsplit("_", 1)[0]].append(hit)
        ids = sorted(by_entry)
        print(f"Retrieving {len(ids)} unique experimental structures", flush=True)
        rows, audit, done = [], [], 0
        batches = [ids[i:i + args.batch_size] for i in range(0, len(ids), args.batch_size)]
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(fetch_entries, client, batch) for batch in batches]
            for future in as_completed(futures):
                entries, path = future.result()
                relative = os.path.relpath(path, args.out).replace(os.sep, "/")
                for entry in entries:
                    entry_rows, entry_audit = process_entry(entry, by_entry[entry["rcsb_id"]], organisms,
                                                          args.contact_cutoff, relative)
                    rows.extend(entry_rows)
                    audit.extend(entry_audit)
                done += len(entries)
                print(f"Processed {done}/{len(ids)} structures; {len(rows)} rows", flush=True)
        accessions = sorted({r["UniProt_accession"] for r in rows if r["UniProt_accession"]})
        print(f"Verifying source organisms for {len(accessions)} UniProt accessions", flush=True)
        uniprot_sources = {}
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(fetch_uniprot_sources, client, accessions[i:i + 25])
                       for i in range(0, len(accessions), 25)]
            for future in as_completed(futures):
                uniprot_sources.update(future.result())
        validate_uniprot_sources(rows, uniprot_sources, organisms)
        rows.sort(key=lambda r: (r["AA"], r["PDB_ID"], r["ligand_instance"], r["protein_chain"], r["UniProt_accession"]))
        audit.sort(key=lambda r: (r["AA"], r["search_hit"], r["ligand_instance"]))
        if {h["search_hit"] for h in hits} != {a["search_hit"] for a in audit}:
            raise RuntimeError("Not all search hits were accounted for in the instance audit")
        controls = select_representatives(rows)
        files = {
            "all_AA_PDB_complexes.tsv": (rows, FIELDS),
            "nonredundant_AA_protein_controls.tsv": (controls, FIELDS + ["candidate_structure_count", "candidate_row_count", "representative_selection"]),
            "ligand_instance_audit.tsv": (audit, AUDIT_FIELDS),
        }
        # Store immutable run tables alongside raw responses, then publish copies.
        for filename, (table, fields) in files.items():
            write_tsv(raw_dir / filename, table, fields)
            write_tsv(args.out / filename, table, fields)
        summaries = []
        for group in organisms:
            for aa in aa_list:
                subset = [r for r in rows if r["AA"] == aa and group in r["organism_group"].split(";")]
                reps = [r for r in controls if r["AA"] == aa and group in r["UniProt_organism_group"].split(";")]
                summaries.append(dict(organism_group=group, AA=aa, structures=len({r["PDB_ID"] for r in subset}),
                                      ligand_instances=len({r["ligand_instance"] for r in subset}),
                                      complex_rows=len(subset), nonredundant_controls=len(reps)))
        summary_fields = ["organism_group", "AA", "structures", "ligand_instances", "complex_rows", "nonredundant_controls"]
        write_tsv(args.out / "AA_organism_summary.tsv", summaries, summary_fields)
        write_tsv(raw_dir / "AA_organism_summary.tsv", summaries, summary_fields)
        manifest.update(status="complete", completed_at_utc=now(), search_hit_count=len(hits),
                        retrieved_structure_count=len(ids), all_complex_rows=len(rows),
                        all_complex_structure_count=len({r["PDB_ID"] for r in rows}),
                        target_contact_structure_count=len({r["PDB_ID"] for r in rows if r["status"] == "target_protein_contact"}),
                        nonredundant_control_count=len(controls), row_status_counts=dict(Counter(r["status"] for r in rows)),
                        mapping_status_counts=dict(Counter(r["UniProt_mapping_status"] for r in rows)),
                        uniprot_source_status_counts=dict(Counter(r["UniProt_source_status"] for r in rows)),
                        processing_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        instance_audit_status_counts=dict(Counter(a["status"] for a in audit)),
                        output_sha256={name: hashlib.sha256((args.out / name).read_bytes()).hexdigest() for name in files})
        atomic_json(manifest_path, manifest)
        atomic_json(args.out / "dataset_manifest.json", dict(manifest, raw_snapshot=str(raw_dir)))
        print(f"Complete: {len(rows)} all-complex rows; {len(controls)} AA x UniProt controls", flush=True)
    except BaseException as exc:
        manifest.update(status="incomplete", last_error=str(exc), last_error_at_utc=now())
        atomic_json(manifest_path, manifest)
        print(f"Run incomplete; raw responses retained. Resume with --resume \"{raw_dir}\" and the same options.", file=sys.stderr)
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
