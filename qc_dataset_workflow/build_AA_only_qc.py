"""Select WT single-protein + one free-AA assemblies with no other non-water compounds."""
from __future__ import annotations
import collections
import csv
import hashlib
import json
import os
from pathlib import Path

from assess_chain_compatibility import BASE, ROOT, CACHE, COORDS, ORGS, read_categories
from build_strict_aa_controls import operations
from curate_biological_controls import read_table, rebase_paths, rank, row_key, site_decision, strict_biological_veto, write_table
from curate_chain_site_compatibility import classify, key

FOLDERS = {'ecoli': 'E_coli', 'arabidopsis': 'Arabidopsis', 'human': 'Human', 'mouse': 'Mouse', 'yeast': 'Yeast'}
FILES = ('strict_WT_single_protein_AA_controls.tsv', 'nonredundant_AA_protein_controls.tsv')
WATER = {'HOH', 'DOD', 'WAT'}
SIDECHAINS = dict(ALA='CB', ARG='CB CG CD NE CZ NH1 NH2', ASN='CB CG OD1 ND2', ASP='CB CG OD1 OD2',
                 CYS='CB SG', GLN='CB CG CD OE1 NE2', GLU='CB CG CD OE1 OE2', GLY='',
                 HIS='CB CG ND1 CD2 CE1 NE2', ILE='CB CG1 CG2 CD1', LEU='CB CG CD1 CD2',
                 LYS='CB CG CD CE NZ', MET='CB CG SD CE', PHE='CB CG CD1 CD2 CE1 CE2 CZ',
                 PRO='CB CG CD', SER='CB OG', THR='CB OG1 CG2', TRP='CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2',
                 TYR='CB CG CD1 CD2 CE1 CE2 CZ OH', VAL='CB CG1 CG2')
PROTEIN_ATOMS = {aa: set(('N CA C O ' + side).split()) for aa, side in SIDECHAINS.items()}
AUDIT_FIELDS = ['organism', 'AA', 'UniProt_accession', 'PDB_ID', 'ligand_instance', 'protein_chain',
                'source_row_number', 'decision', 'reasons', 'biological_classification', 'chain_compatibility',
                'assembly_composition_json', 'missing_contact_atoms_json', 'selected_as_representative', 'references']


def assembly_inventory(categories, assembly_ids):
    entities = {r['id']: r['type'] for r in categories['_entity']}
    asym_entity = {r['id']: r['entity_id'] for r in categories['_struct_asym']}
    polymers = {r['entity_id']: r['type'] for r in categories['_entity_poly']}
    components = collections.defaultdict(set)
    for row in categories['_pdbx_entity_nonpoly']:
        components[row['entity_id']].add(row['comp_id'])
    for row in categories['_atom_site']:
        components[asym_entity.get(row['label_asym_id'], '')].add(row['label_comp_id'])
    result = []
    for assembly_id in assembly_ids:
        copies = set()
        for generator in categories['_pdbx_struct_assembly_gen']:
            if generator['assembly_id'] != assembly_id:
                continue
            ops = operations(generator['oper_expression'])
            if not ops:
                raise ValueError('Unresolved assembly operators')
            for label in generator['asym_id_list'].split(','):
                copies.update((label.strip(), op) for op in ops)
        if not copies:
            raise ValueError('No assembly composition for ' + assembly_id)
        entries = []
        for label, op in sorted(copies):
            entity = asym_entity.get(label, '')
            kind = entities.get(entity, 'unknown')
            comps = sorted(components.get(entity, set()))
            if kind == 'water' or (kind == 'non-polymer' and comps and set(comps) <= WATER):
                continue
            entries.append(dict(label_asym_id=label, operations=list(op), entity_type=kind,
                                polymer_type=polymers.get(entity, ''), comp_ids=comps))
        labels = {label for label, _ in copies}
        outside_links = []
        for link in categories['_struct_conn']:
            if link.get('conn_type_id', '') not in ('covale', 'metalc', 'disulf'):
                continue
            partners = [link.get('ptnr1_label_asym_id', ''), link.get('ptnr2_label_asym_id', '')]
            if any(label in labels for label in partners):
                for label in partners:
                    if label and label not in labels and entities.get(asym_entity.get(label, '')) != 'water':
                        outside_links.append(dict(label_asym_id=label, connection_type=link['conn_type_id']))
        result.append(dict(assembly_id=assembly_id, non_water_components=entries, outside_linked_components=outside_links))
    return result


def composition_reasons(inventory, row):
    reasons = []
    for assembly in inventory:
        entries = assembly['non_water_components']
        proteins = [e for e in entries if e['entity_type'] == 'polymer' and e['polymer_type'] == 'polypeptide(L)']
        other_polymers = [e for e in entries if e['entity_type'] in ('polymer', 'branched') and e not in proteins]
        ligands = [e for e in entries if e['entity_type'] == 'non-polymer']
        if assembly.get('outside_linked_components'):
            reasons.append('additional_covalently_attached_or_metal_linked_component')
        if len(proteins) != 1 or proteins[0]['label_asym_id'] != row['protein_chain']:
            reasons.append('not_exactly_one_selected_protein_chain')
        if other_polymers:
            reasons.append('additional_polymer_or_glycan')
        if any(e['entity_type'] not in ('polymer', 'non-polymer', 'branched') for e in entries):
            reasons.append('assembly_component_identity_unknown')
        if len(ligands) != 1 or ligands[0]['label_asym_id'] != row['ligand_label_asym_id'] or ligands[0]['comp_ids'] != [row['AA']]:
            reasons.append('not_exactly_one_target_AA_and_no_other_nonwater_compounds')
    return sorted(set(reasons))


def missing_contact_atoms(categories, row, site):
    atoms = collections.defaultdict(set)
    for atom in categories['_atom_site']:
        if atom['label_asym_id'] != row['protein_chain'] or atom.get('pdbx_PDB_model_num', '1') != str(site['selected_model']):
            continue
        if atom['type_symbol'] in ('H', 'D') or float(atom.get('occupancy', '1')) <= 0:
            continue
        atoms[atom['label_seq_id'], atom['label_comp_id']].add(atom['label_atom_id'])
    missing = []
    for contact in site['asymmetric_unit_contacts']:
        if not contact['selected_chain']:
            continue
        for residue in contact['residues']:
            aa = residue['comp_id']
            required = PROTEIN_ATOMS.get(aa)
            present = atoms[str(residue['label_seq_id']), aa]
            absent = sorted(required - present) if required is not None else ['unknown_noncanonical_residue']
            if absent:
                missing.append(dict(label_seq_id=residue['label_seq_id'], comp_id=aa, missing_atoms=absent))
    return missing


def evaluate(row, biological_classification, geometry, categories):
    reasons, composition, missing = [], [], []
    if biological_classification != 'GOLD':
        reasons.append('biological_evidence_not_validated')
    if row['strict_control_eligible'] != '1':
        reasons.append('fails_WT_no_ligand_metal_or_single_contact_chain_rules')
    if row['ultra_strict_monomer_pass'] != '1':
        reasons.append('biological_assembly_not_unambiguously_monomeric')
    if strict_biological_veto(row):
        reasons.append('literature_identifies_ligand_metal_coordination')
    if biological_classification == 'GOLD' and site_decision(row, 'GOLD')[0] != 'GOLD':
        reasons.append('functional_AA_site_not_validated')
    if row.get('protein_nonstandard_monomer_count', '') != '0':
        reasons.append('noncanonical_polymer_modification_or_unknown_count')
    if reasons:
        return reasons, '', composition, missing
    decision = classify(row, geometry)
    if decision['classification'] != 'monomer_compatible':
        reasons.append('binding_site_' + decision['classification'])
    if geometry is None or categories is None:
        reasons.append('occupied_coordinate_audit_unresolved')
        return reasons, decision['classification'], composition, missing
    expected = PROTEIN_ATOMS[row['AA']] | {'OXT'}
    if set(geometry['ligand_atom_names']) != expected:
        reasons.append('free_AA_heavy_atom_model_incomplete_or_unexpected')
    if len([alt for alt in geometry['alt_ids'] if alt not in ('.', '?')]) > 1:
        reasons.append('multiple_ligand_conformers_need_pose_review')
    composition = assembly_inventory(categories, row['chain_candidate_assembly_ids'].split(';'))
    reasons.extend(composition_reasons(composition, row))
    missing = missing_contact_atoms(categories, row, geometry)
    if missing:
        reasons.append('contacting_protein_residues_have_missing_heavy_atoms')
    return sorted(set(reasons)), decision['classification'], composition, missing


def main():
    out = BASE / 'qc_protein_AA_only_controls'
    out.mkdir(exist_ok=True)
    geometry = {key(r): r for r in json.loads((CACHE / 'site_geometry.json').read_text(encoding='utf8'))}
    report = dict(status='PASS', organisms={}, app_executed=False, criteria='WT; one protein chain; one free AA; no other nonwater components; complete modeled pocket')
    for org in ORGS:
        folder, destination = BASE / FOLDERS[org], out / org
        destination.mkdir(exist_ok=True)
        _, bio = read_table(ROOT / org / 'biological_curation.tsv')
        biological = {(r['AA'], r['UniProt_accession']): r for r in bio}
        _, candidates = read_table(folder / 'all_AA_PDB_complexes.tsv')
        headers = {name: read_table(folder / name)[0] for name in FILES}
        _, original_representatives = read_table(folder / FILES[1])
        previous_by_pair = {(r['AA'], r['UniProt_accession']): r for r in original_representatives}
        protected = {name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in ('all_AA_PDB_complexes.tsv', *FILES)}
        groups, audits, coordinate_cache = collections.defaultdict(list), [], {}
        for number, candidate in enumerate(candidates, 2):
            row = dict(candidate, organism_cohort=org)
            bio_record = biological[row['AA'], row['UniProt_accession']]
            site = geometry.get(key(row))
            categories = None
            if site and row['ultra_strict_monomer_pass'] == '1':
                pdb = row['PDB_ID']
                if pdb not in coordinate_cache:
                    coordinate_cache[pdb] = read_categories(COORDS / (pdb + '.cif.gz'))
                categories = coordinate_cache[pdb]
            try:
                reasons, compatibility, composition, missing = evaluate(row, bio_record['classification'], site, categories)
            except Exception as error:
                reasons, compatibility, composition, missing = ['composition_or_pocket_audit_unresolved:' + str(error)], 'uncertain', [], []
            passed = not reasons
            if passed:
                groups[row['AA'], row['UniProt_accession']].append(candidate)
            audits.append(dict(organism=org, AA=row['AA'], UniProt_accession=row['UniProt_accession'], PDB_ID=row['PDB_ID'],
                ligand_instance=row['ligand_instance'], protein_chain=row['protein_chain'], source_row_number=number,
                decision='PASS' if passed else 'EXCLUDE', reasons=';'.join(reasons), biological_classification=bio_record['classification'],
                chain_compatibility=compatibility, assembly_composition_json=json.dumps(composition, separators=(',', ':')),
                missing_contact_atoms_json=json.dumps(missing, separators=(',', ':')), selected_as_representative='0',
                references=bio_record['reference_urls'] + ';https://www.rcsb.org/structure/' + row['PDB_ID']))
        selected = {pair: min(rows, key=rank) for pair, rows in groups.items()}
        selected_ids = {(r['AA'], r['UniProt_accession'], r['PDB_ID'], r['ligand_instance'], r['protein_chain']) for r in selected.values()}
        for audit in audits:
            identity = tuple(audit[field] for field in ('AA', 'UniProt_accession', 'PDB_ID', 'ligand_instance', 'protein_chain'))
            audit['selected_as_representative'] = str(int(identity in selected_ids and audit['decision'] == 'PASS'))
        for name in FILES:
            rows = []
            for pair, candidate in sorted(selected.items()):
                row = dict(candidate)
                if name == FILES[0]:
                    previous = previous_by_pair.get(pair)
                    row.update(strict_candidate_structure_count=str(len({r['PDB_ID'] for r in groups[pair]})),
                               strict_candidate_row_count=str(len(groups[pair])), control_set='strict',
                               previous_representative_PDB_ID=previous['PDB_ID'] if previous else '',
                               representative_changed_after_strict_filter=str(int(previous is None or row_key(previous) != row_key(row))))
                else:
                    row.update(candidate_structure_count=str(len({r['PDB_ID'] for r in groups[pair]})), candidate_row_count=str(len(groups[pair])),
                               representative_selection='lowest_available_resolution_then_most_contacts_then_identifiers')
                row = rebase_paths(row, folder, destination)
                rows.append({column: row[column] for column in headers[name]})
            write_table(destination / name, headers[name], rows)
            actual_headers, actual = read_table(destination / name)
            assert actual_headers == headers[name]
            assert len(actual) == len(selected)
            assert len({(r['AA'], r['UniProt_accession']) for r in actual}) == len(actual)
            for row in actual:
                assert row['strict_control_eligible'] == '1' and row['ultra_strict_monomer_pass'] == '1'
                for field in ('raw_response_file', 'strict_metadata_response_file'):
                    assert (destination / row[field]).resolve().is_file()
        write_table(destination / 'qc_composition_audit.tsv', AUDIT_FIELDS, audits)
        assert all(hashlib.sha256((folder / name).read_bytes()).hexdigest() == expected for name, expected in protected.items())
        report['organisms'][org] = dict(all_candidates_evaluated=len(audits), eligible_alternatives=sum(r['decision'] == 'PASS' for r in audits),
            selected_AA_protein_pairs=len(selected), proteins=len({r['UniProt_accession'] for r in selected.values()}),
            AA_coverage=sorted({r['AA'] for r in selected.values()}), exclusion_reasons=dict(collections.Counter(reason for r in audits for reason in r['reasons'].split(';') if reason)))
        print(f'{org}: {len(selected)} protein + AA-only QC representatives', flush=True)
    (out / 'validation_report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
