#!/usr/bin/env python3
"""Select WT/non-metal/single-contact-chain controls and optional monomer controls."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from build_aa_pdb_dataset import ApiClient, DATA_URL, FREE, atomic_json, now, representative_rank, write_tsv
from split_aa_pdb_by_organism import DIRECTORIES, assigned_groups, read_tsv

ALL_FILE = "all_AA_PDB_complexes.tsv"
CONTROL_FILE = "nonredundant_AA_protein_controls.tsv"
STRICT_FILE = "strict_WT_single_protein_AA_controls.tsv"
ULTRA_FILE = "ultra_strict_WT_monomer_AA_controls.tsv"

QUERY = """query StrictMetadata($ids: [String!]!) {
 entries(entry_ids: $ids) {
  rcsb_id
  rcsb_entry_info {
   polymer_entity_count_protein polymer_entity_count
   deposited_polymer_entity_instance_count inter_mol_metalic_bond_count
  }
  polymer_entities {
   rcsb_id
   rcsb_polymer_entity_container_identifiers { entity_id asym_ids uniprot_ids }
   entity_poly {
    rcsb_entity_polymer_type rcsb_mutation_count rcsb_conflict_count
    rcsb_deletion_count rcsb_insertion_count rcsb_artifact_monomer_count
    rcsb_non_std_monomer_count
   }
   rcsb_polymer_entity {
    pdbx_mutation pdbx_fragment rcsb_source_part_count rcsb_source_taxonomy_count
   }
  }
  assemblies {
   rcsb_id
   pdbx_struct_assembly {
    id details oligomeric_count oligomeric_details rcsb_candidate_assembly
   }
   pdbx_struct_assembly_gen { asym_id_list oper_expression }
   rcsb_assembly_info {
    polymer_entity_count polymer_entity_count_protein
    polymer_entity_instance_count polymer_entity_instance_count_protein
    polymer_composition branched_entity_instance_count
   }
  }
  nonpolymer_entities {
   rcsb_id
   nonpolymer_entity_instances {
    rcsb_id
    rcsb_nonpolymer_entity_instance_container_identifiers { asym_id comp_id }
    rcsb_nonpolymer_instance_annotation { comp_id type }
    rcsb_nonpolymer_struct_conn {
     connect_type dist_value
     connect_partner { label_asym_id label_comp_id label_atom_id }
     connect_target { label_asym_id label_comp_id label_atom_id }
    }
    rcsb_target_neighbors { comp_id distance target_asym_id target_entity_id }
   }
  }
 }
}"""

ADDED_FIELDS = [
    "protein_mutations_current", "protein_mutation_count", "protein_sequence_conflict_count",
    "protein_deletion_count", "protein_insertion_count", "protein_artifact_region_count",
    "protein_nonstandard_monomer_count", "protein_fragment", "protein_source_taxonomy_count",
    "protein_UniProt_accession_count", "strict_WT_status", "strict_WT_pass",
    "entry_protein_entity_count", "entry_protein_chain_count", "entry_polymer_chain_count",
    "entry_nonprotein_polymer_entity_count", "entry_protein_label_chain_ids",
    "ligand_contacting_protein_chain_count", "ligand_contacting_protein_chains",
    "entry_metal_link_count", "ligand_has_metal_coordination_current",
    "ligand_metal_connection_count", "ligand_metal_partner_comp_ids", "strict_metal_status",
    "strict_no_metal_pass", "chain_candidate_assembly_ids", "ligand_chain_assembly_ids",
    "monomer_assembly_ids", "nonmonomer_assembly_ids", "unknown_assembly_ids",
    "selected_biological_assembly_id", "assembly_oligomeric_details",
    "assembly_protein_entity_count", "assembly_protein_chain_count", "assembly_polymer_chain_count",
    "assembly_branched_chain_count", "assembly_composition", "assembly_details_json",
    "strict_single_protein_status", "strict_single_protein_pass",
    "strict_control_eligible", "strict_exclusion_reasons",
    "ultra_strict_monomer_status", "ultra_strict_monomer_pass",
    "ultra_strict_control_eligible", "ultra_strict_exclusion_reasons", "strict_metadata_response_file",
]


def fetch_metadata(client, ids):
    data, path = client.post(DATA_URL, {"query": QUERY, "variables": {"ids": ids}}, "metadata")
    entries = ((data or {}).get("data") or {}).get("entries") or []
    if len(entries) != len(ids) or {e["rcsb_id"] for e in entries if e} != set(ids):
        raise RuntimeError(f"Incomplete strict metadata batch: {ids}")
    return entries, path


def operations(expression):
    """Expand mmCIF operation products without guessing malformed expressions."""
    expression = (expression or "").replace(" ", "")
    if not expression:
        return set()
    if "(" in expression:
        groups = re.findall(r"\(([^()]*)\)", expression)
        if "".join(f"({g})" for g in groups) != expression:
            return set()
    else:
        groups = [expression]
    expanded = []
    for group in groups:
        values = []
        for token in group.split(","):
            match = re.fullmatch(r"(\d+)-(\d+)", token)
            if match:
                begin, end = map(int, match.groups())
                if end < begin or end - begin > 10000:
                    return set()
                values.extend(map(str, range(begin, end + 1)))
            elif re.fullmatch(r"[A-Za-z0-9_]+", token):
                values.append(token)
            else:
                return set()
        expanded.append(values)
    if not expanded:
        return set()
    size = 1
    for values in expanded:
        size *= len(values)
    if size > 100000:
        return set()
    return set(itertools.product(*expanded))


def assembly_records(entry, chain, ligand_chain):
    result = []
    for assembly in entry.get("assemblies") or []:
        description = assembly.get("pdbx_struct_assembly") or {}
        if description.get("rcsb_candidate_assembly") != "Y":
            continue  # An asymmetric unit is not evidence of a biological monomer.
        generators = assembly.get("pdbx_struct_assembly_gen") or []
        chain_generators = [g for g in generators if chain in (g.get("asym_id_list") or [])]
        if not chain or not chain_generators:
            continue
        ligand_generators = [g for g in generators if ligand_chain in (g.get("asym_id_list") or [])]
        chain_ops = set().union(*(operations(g.get("oper_expression")) for g in chain_generators))
        ligand_ops = set().union(*(operations(g.get("oper_expression")) for g in ligand_generators))
        info = assembly.get("rcsb_assembly_info") or {}
        counts = [info.get(k) for k in ("polymer_entity_count_protein", "polymer_entity_instance_count_protein", "polymer_entity_instance_count")]
        if any(c is None for c in counts):
            classification = "unknown"
        elif counts == [1, 1, 1] and description.get("oligomeric_count") in (None, 1):
            classification = "monomer"
        else:
            classification = "nonmonomer"
        result.append({"assembly_id": description.get("id") or assembly["rcsb_id"].rsplit("-", 1)[-1],
                       "classification": classification, "protein_entity_count": counts[0],
                       "protein_chain_count": counts[1], "polymer_chain_count": counts[2],
                       "branched_chain_count": info.get("branched_entity_instance_count"),
                       "oligomeric_details": description.get("oligomeric_details"),
                       "composition": info.get("polymer_composition"),
                       "details": description.get("details"),
                       "ligand_present": bool(ligand_generators),
                       "chain_and_ligand_share_operator": bool(chain_ops.intersection(ligand_ops))})
    return sorted(result, key=lambda a: a["assembly_id"])


def no_mutation_text(value):
    return str(value or "").strip().upper() in {"", ".", "?", "NONE", "NO", "NO MUTATION", "NO MUTATIONS", "WILD TYPE", "WILD-TYPE", "WT"}


def annotate(row, entry, raw_path):
    new = dict(row)
    new.update({field: "" for field in ADDED_FIELDS})
    new["strict_metadata_response_file"] = raw_path
    polymers = entry.get("polymer_entities") or []
    proteins = [p for p in polymers if (p.get("entity_poly") or {}).get("rcsb_entity_polymer_type") == "Protein"]
    chains = {chain for p in proteins for chain in (p.get("rcsb_polymer_entity_container_identifiers") or {}).get("asym_ids") or []}
    einfo = entry.get("rcsb_entry_info") or {}
    new.update(entry_protein_entity_count=len(proteins), entry_protein_chain_count=len(chains),
               entry_polymer_chain_count=einfo.get("deposited_polymer_entity_instance_count"),
               entry_nonprotein_polymer_entity_count=len(polymers) - len(proteins),
               entry_protein_label_chain_ids=";".join(sorted(chains)),
               entry_metal_link_count=einfo.get("inter_mol_metalic_bond_count"))
    polymer = next((p for p in proteins if (p.get("rcsb_polymer_entity_container_identifiers") or {}).get("entity_id") == row["protein_entity_id"]), None)
    wt_reasons = []
    if not polymer:
        wt_reasons.append("protein_metadata_unavailable")
    else:
        poly = polymer.get("entity_poly") or {}
        details = polymer.get("rcsb_polymer_entity") or {}
        ids = polymer.get("rcsb_polymer_entity_container_identifiers") or {}
        new.update(protein_mutations_current=details.get("pdbx_mutation") or "", protein_fragment=details.get("pdbx_fragment") or "",
                   protein_source_taxonomy_count=details.get("rcsb_source_taxonomy_count"),
                   protein_UniProt_accession_count=len(ids.get("uniprot_ids") or []))
        mapping = {"protein_mutation_count": "rcsb_mutation_count", "protein_sequence_conflict_count": "rcsb_conflict_count",
                   "protein_deletion_count": "rcsb_deletion_count", "protein_insertion_count": "rcsb_insertion_count",
                   "protein_artifact_region_count": "rcsb_artifact_monomer_count", "protein_nonstandard_monomer_count": "rcsb_non_std_monomer_count"}
        for column, source in mapping.items():
            new[column] = poly.get(source)
        if not no_mutation_text(row.get("protein_mutations")) or not no_mutation_text(details.get("pdbx_mutation")):
            wt_reasons.append("mutation_text_present")
        for name in ("protein_mutation_count", "protein_sequence_conflict_count", "protein_deletion_count", "protein_insertion_count"):
            if new[name] is None:
                wt_reasons.append(name + "_unknown")
            elif new[name] != 0:
                wt_reasons.append(name + "_positive")
        if (details.get("rcsb_source_taxonomy_count") or 0) > 1 or new["protein_UniProt_accession_count"] > 1:
            wt_reasons.append("chimeric_or_multi_accession_construct")
    new["strict_WT_status"] = ";".join(wt_reasons) or "no_reported_mutations_or_sequence_discrepancies"
    new["strict_WT_pass"] = int(not wt_reasons)

    ligand = next((lig for entity in entry.get("nonpolymer_entities") or [] for lig in entity.get("nonpolymer_entity_instances") or []
                   if lig["rcsb_id"] == row["ligand_instance"]), None)
    metal_reasons = []
    free_confirmed = False
    if ligand is None:
        metal_reasons.append("ligand_metadata_unavailable")
    else:
        annotations = [a for a in ligand.get("rcsb_nonpolymer_instance_annotation") or [] if a.get("comp_id") == row["AA"]]
        free_confirmed = any(a.get("type") == FREE for a in annotations)
        links = ligand.get("rcsb_nonpolymer_struct_conn") or []
        metal_links = [link for link in links if link.get("connect_type") == "metal coordination"]
        metal = any(a.get("type") == "HAS_METAL_COORDINATION_LINKAGE" for a in annotations) or bool(metal_links)
        new.update(ligand_has_metal_coordination_current=int(metal), ligand_metal_connection_count=len(metal_links),
                   ligand_metal_partner_comp_ids=";".join(sorted({(link.get("connect_partner") or {}).get("label_comp_id", "") for link in metal_links})))
        if metal or str(row.get("has_metal_coordination")) == "1":
            metal_reasons.append("ligand_metal_coordination")
        if not free_confirmed or any(a.get("type") == "HAS_COVALENT_LINKAGE" for a in annotations) or any(link.get("connect_type") == "covalent bond" for link in links):
            metal_reasons.append("free_noncovalent_ligand_not_confirmed")
    new["strict_no_metal_pass"] = int(not metal_reasons)
    new["strict_metal_status"] = ";".join(metal_reasons) or "no_ligand_metal_linkage_reported"

    contact_chains = {n["target_asym_id"] for n in (ligand or {}).get("rcsb_target_neighbors") or []
                      if n.get("comp_id") == row["AA"] and n.get("target_asym_id") in chains
                      and isinstance(n.get("distance"), (int, float)) and 0 <= n["distance"] <= float(row["contact_cutoff_A"])}
    new["ligand_contacting_protein_chain_count"] = len(contact_chains)
    new["ligand_contacting_protein_chains"] = ";".join(sorted(contact_chains))
    records = assembly_records(entry, row["protein_chain"], row["ligand_label_asym_id"])
    new["assembly_details_json"] = json.dumps(records, separators=(",", ":"))
    new["chain_candidate_assembly_ids"] = ";".join(a["assembly_id"] for a in records)
    new["ligand_chain_assembly_ids"] = ";".join(a["assembly_id"] for a in records if a["chain_and_ligand_share_operator"])
    for classification, column in [("monomer", "monomer_assembly_ids"), ("nonmonomer", "nonmonomer_assembly_ids"), ("unknown", "unknown_assembly_ids")]:
        new[column] = ";".join(a["assembly_id"] for a in records if a["classification"] == classification)
    mono = [a for a in records if a["classification"] == "monomer" and a["chain_and_ligand_share_operator"]]
    assembly_reasons = []
    if not records:
        assembly_reasons.append("no_candidate_biological_assembly_for_chain")
    if new["nonmonomer_assembly_ids"]:
        assembly_reasons.append("nonmonomeric_candidate_assembly")
    if new["unknown_assembly_ids"]:
        assembly_reasons.append("assembly_stoichiometry_unknown")
    if not mono:
        assembly_reasons.append("no_monomer_assembly_with_ligand_and_chain")
    chain_reasons = []
    if contact_chains != {row["protein_chain"]}:
        chain_reasons.append("ligand_contacts_multiple_or_unassigned_protein_chains")
    # Also expose the stoichiometry of rejected complexes in the flat columns.
    # The full set of alternative assemblies remains in assembly_details_json.
    if records:
        selected = (mono or [a for a in records if a["chain_and_ligand_share_operator"]] or records)[0]
        new.update(selected_biological_assembly_id=selected["assembly_id"], assembly_oligomeric_details=selected["oligomeric_details"],
                   assembly_protein_entity_count=selected["protein_entity_count"], assembly_protein_chain_count=selected["protein_chain_count"],
                   assembly_polymer_chain_count=selected["polymer_chain_count"], assembly_branched_chain_count=selected["branched_chain_count"],
                   assembly_composition=selected["composition"])
    new["strict_single_protein_pass"] = int(not chain_reasons)
    new["strict_single_protein_status"] = ";".join(chain_reasons) or "ligand_contacts_exactly_one_protein_chain"
    new["ultra_strict_monomer_pass"] = int(not assembly_reasons)
    new["ultra_strict_monomer_status"] = ";".join(assembly_reasons) or "monomeric_biological_assembly"
    # Biological assembly metadata never gates the main strict set. It is
    # retained for inspection and the optional monomer-only secondary set.
    reasons = wt_reasons + metal_reasons + chain_reasons
    if str(row["control_eligible"]) != "1" or row["status"] != "target_protein_contact":
        reasons.insert(0, "not_eligible_in_original_dataset")
    if row["UniProt_source_status"] != "confirmed_target_source":
        reasons.insert(0, "UniProt_target_source_not_confirmed")
    new["strict_control_eligible"] = int(not reasons)
    new["strict_exclusion_reasons"] = ";".join(reasons)
    new["ultra_strict_control_eligible"] = int(not (reasons + assembly_reasons))
    new["ultra_strict_exclusion_reasons"] = ";".join(reasons + assembly_reasons)
    return new


def row_key(row):
    return tuple(row[k] for k in ("AA", "PDB_ID", "ligand_instance", "protein_chain", "UniProt_accession"))


def select_strict(rows, ultra=False):
    prefix = "ultra_strict" if ultra else "strict"
    groups = defaultdict(list)
    for row in rows:
        if str(row[prefix + "_control_eligible"]) == "1":
            groups[row["AA"], row["UniProt_accession"]].append(row)
    return [dict(min(candidates, key=representative_rank),
                 **{prefix + "_candidate_structure_count": len({r["PDB_ID"] for r in candidates}),
                    prefix + "_candidate_row_count": len(candidates),
                    "control_set": "ultra_strict_monomer" if ultra else "strict"})
            for _, candidates in sorted(groups.items())]


def rebase(rows, source, destination):
    output = []
    for row in rows:
        row = dict(row)
        for field in ("raw_response_file", "strict_metadata_response_file"):
            if row.get(field):
                row[field] = os.path.relpath(source / row[field], destination).replace(os.sep, "/")
        output.append(row)
    return output


def run(args):
    root = args.out.resolve()
    base_fields, original_rows = read_tsv(root / ALL_FILE)
    control_fields, original_controls = read_tsv(root / CONTROL_FILE)
    # Re-running enrichment replaces only these derived annotations.
    base_fields = [f for f in base_fields if f not in ADDED_FIELDS]
    control_fields = [f for f in control_fields if f not in ADDED_FIELDS]
    original_rows = [{f: row[f] for f in base_fields} for row in original_rows]
    original_controls = [{f: row[f] for f in control_fields} for row in original_controls]
    raw = (args.resume.resolve() if args.resume else root / "raw" / ("strict_" + now().replace("-", "").replace(":", "")[:15]))
    manifest_path = raw / "manifest.json"
    input_hash = hashlib.sha256(json.dumps(original_rows, sort_keys=True).encode()).hexdigest()
    if args.resume:
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous["input_rows_sha256"] != input_hash:
            raise ValueError("Resume input rows differ from cached snapshot")
    manifest = {"status": "running", "started_at_utc": now(), "input_rows_sha256": input_hash,
                "filter_version": 2,
                "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "main_control_file": STRICT_FILE, "optional_ultra_strict_control_file": ULTRA_FILE,
                "main_rule": "WT + no ligand metal coordination + ligand contacts exactly the selected protein chain; original control eligibility and confirmed UniProt source required. No biological assembly restriction.",
                "assembly_rule": "Informational only for the main strict set; oligomeric, missing, or conflicting biological assemblies do not exclude a main-set candidate.",
                "ultra_strict_assembly_rule": "Main strict criteria plus all candidate biological assemblies containing the chain must be monomeric; at least one must include the ligand under the same operator. Independent crystal copies allowed.",
                "WT_rule": "No mutation text, mutation/conflict/deletion/insertion counts zero; no multi-source/multi-accession chimeras. Fragments and expression tags are reported but allowed.",
                "metal_rule": "Exclude ligand metal-coordination annotations or explicit metal links. Metals elsewhere in the entry do not by themselves exclude the row."}
    atomic_json(manifest_path, manifest)
    client = ApiClient(raw)
    print(f"Strict metadata snapshot: {raw}", flush=True)
    try:
        ids = sorted({r["PDB_ID"] for r in original_rows})
        metadata, completed = {}, 0
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(fetch_metadata, client, ids[i:i + 10]) for i in range(0, len(ids), 10)]
            for future in as_completed(futures):
                entries, path = future.result()
                for entry in entries:
                    metadata[entry["rcsb_id"]] = (entry, os.path.relpath(path, root).replace(os.sep, "/"))
                completed += len(entries)
                if completed % 100 == 0 or completed == len(ids):
                    print(f"Retrieved assembly/mutation/metal metadata: {completed}/{len(ids)} structures", flush=True)
        rows = [annotate(r, *metadata[r["PDB_ID"]]) for r in original_rows]
        lookup = {row_key(r): r for r in rows}
        controls = [dict(r, **{f: lookup[row_key(r)][f] for f in ADDED_FIELDS}) for r in original_controls]
        strict = select_strict(rows)
        ultra = select_strict(rows, ultra=True)
        original_reps = {(r["AA"], r["UniProt_accession"]): r for r in original_controls}
        for row in strict + ultra:
            previous = original_reps.get((row["AA"], row["UniProt_accession"]), {})
            row["previous_representative_PDB_ID"] = previous.get("PDB_ID", "")
            row["representative_changed_after_strict_filter"] = int(row_key(row) != row_key(previous)) if previous else 1
        selection_fields = ["control_set", "previous_representative_PDB_ID", "representative_changed_after_strict_filter"]
        strict_fields = base_fields + ADDED_FIELDS + ["strict_candidate_structure_count", "strict_candidate_row_count"] + selection_fields
        ultra_fields = base_fields + ADDED_FIELDS + ["ultra_strict_candidate_structure_count", "ultra_strict_candidate_row_count"] + selection_fields
        # Snapshot the source TSVs before adding columns. Raw archive is append-only.
        for name, table, fields in [(ALL_FILE, original_rows, base_fields), (CONTROL_FILE, original_controls, control_fields)]:
            path = raw / ("original_" + name)
            if not path.exists():
                write_tsv(path, table, fields)
        outputs = {ALL_FILE: (rows, base_fields + ADDED_FIELDS), CONTROL_FILE: (controls, control_fields + ADDED_FIELDS),
                   STRICT_FILE: (strict, strict_fields), ULTRA_FILE: (ultra, ultra_fields)}
        summary = []
        for group, directory in DIRECTORIES.items():
            folder = root / directory
            folder.mkdir(exist_ok=True)
            for name, (table, fields) in outputs.items():
                selected = [r for r in table if group in assigned_groups(r, controls=name != ALL_FILE)]
                write_tsv(folder / name, rebase(selected, root, folder), fields)
            group_rows = [r for r in rows if group in assigned_groups(r)]
            group_controls = [r for r in controls if group in assigned_groups(r, controls=True)]
            group_strict = [r for r in strict if group in assigned_groups(r, controls=True)]
            group_ultra = [r for r in ultra if group in assigned_groups(r, controls=True)]
            item = {"organism_group": group, "original_controls": len(group_controls), "strict_controls": len(group_strict),
                    "strict_candidate_rows": sum(r["strict_control_eligible"] for r in group_rows),
                    "strict_structures": len({r["PDB_ID"] for r in group_strict}),
                    "reselected_representatives": sum(r["representative_changed_after_strict_filter"] for r in group_strict),
                    "ultra_strict_controls": len(group_ultra),
                    "ultra_strict_candidate_rows": sum(r["ultra_strict_control_eligible"] for r in group_rows),
                    "ultra_strict_structures": len({r["PDB_ID"] for r in group_ultra})}
            summary.append(item)
            atomic_json(folder / "strict_controls_manifest.json", dict(manifest, status="complete", counts=item,
                        metadata_snapshot=os.path.relpath(raw, folder).replace(os.sep, "/"),
                        output_sha256={name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in outputs}))
            old_manifest_path = folder / "organism_manifest.json"
            if old_manifest_path.exists():
                old_manifest = json.loads(old_manifest_path.read_text(encoding="utf-8"))
                old_manifest.setdefault("output_sha256", {}).update({name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in outputs})
                for name, table in [(STRICT_FILE, group_strict), (ULTRA_FILE, group_ultra)]:
                    old_manifest.setdefault("counts", {})[name] = {"rows": len(table), "structures": len({r["PDB_ID"] for r in table}),
                                                                  "rows_without_confirmed_target_contact": 0}
                old_manifest["strict_enrichment_manifest"] = "strict_controls_manifest.json"
                atomic_json(old_manifest_path, old_manifest)
            print(f"{group}: {len(group_strict)} main strict; {len(group_ultra)} optional monomer controls", flush=True)
        for name, (table, fields) in outputs.items():
            write_tsv(root / name, table, fields)
        for directory in DIRECTORIES.values():
            path = root / directory / "organism_manifest.json"
            if path.exists():
                saved = json.loads(path.read_text(encoding="utf-8"))
                saved.setdefault("pre_strict_enrichment_source_file_sha256", dict(saved.get("source_file_sha256", {})))
                saved.setdefault("source_file_sha256", {}).update({name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in outputs})
                atomic_json(path, saved)
        write_tsv(root / "strict_controls_summary.tsv", summary, list(summary[0]))
        reasons = Counter(reason for row in rows for reason in row["strict_exclusion_reasons"].split(";") if reason)
        ultra_reasons = Counter(reason for row in rows for reason in row["ultra_strict_exclusion_reasons"].split(";") if reason)
        manifest.update(status="complete", completed_at_utc=now(), structures_enriched=len(ids), all_rows=len(rows),
                        strict_eligible_rows=sum(r["strict_control_eligible"] for r in rows), strict_controls=len(strict),
                        ultra_strict_eligible_rows=sum(r["ultra_strict_control_eligible"] for r in rows), ultra_strict_controls=len(ultra),
                        summary=summary, exclusion_reason_counts=dict(reasons),
                        ultra_strict_exclusion_reason_counts=dict(ultra_reasons),
                        output_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in outputs})
        atomic_json(manifest_path, manifest)
        atomic_json(root / "strict_controls_manifest.json", dict(manifest, metadata_snapshot=str(raw)))
        old_manifest_path = root / "dataset_manifest.json"
        old_manifest = json.loads(old_manifest_path.read_text(encoding="utf-8"))
        old_manifest.setdefault("pre_strict_enrichment_output_sha256", dict(old_manifest.get("output_sha256", {})))
        old_manifest["output_sha256"].update(manifest["output_sha256"])
        old_manifest["strict_enrichment_manifest"] = "strict_controls_manifest.json"
        atomic_json(old_manifest_path, old_manifest)
    except BaseException as exc:
        manifest.update(status="incomplete", error=str(exc))
        atomic_json(manifest_path, manifest)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("."))
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 9))
    run(parser.parse_args(argv))


if __name__ == "__main__":
    main()
