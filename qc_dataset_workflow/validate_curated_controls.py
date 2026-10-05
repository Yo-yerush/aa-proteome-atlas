# Adapted from the existing reviewed workflow for fresh datasets.
# Original source SHA256: e77512292cdc1d2783838bd7b711e2d2ab71613e314615638a417d6fe72c78e0
"""Independent TSV validation only. Never imports or executes the app or loader.

Checks the delivered files against original TSV bytes/rows and the curation audit.
No structural computation or evidence fetching occurs here.
"""
import collections
import csv
import datetime
import hashlib
import json
import math
import os
import pathlib
import re

BASE = pathlib.Path(os.environ["AA_QC_DATA_ROOT"]).resolve()
ROOT = BASE / 'curated_controls'
ORGS = {'ecoli': 'E_coli', 'arabidopsis': 'Arabidopsis', 'human': 'Human', 'mouse': 'Mouse', 'yeast': 'Yeast'}
ORGS = {key: value for key, value in ORGS.items() if key in os.environ.get("AA_QC_ORGANISMS", "ecoli,arabidopsis,human,mouse").split(",")}
FILES = ['nonredundant_AA_protein_controls.tsv', 'strict_WT_single_protein_AA_controls.tsv']
PATHS = {'raw_response_file', 'strict_metadata_response_file'}
CANONICAL = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def read(path):
    raw = path.read_bytes()
    assert not raw.startswith(b'\xef\xbb\xbf'), f'Unexpected BOM: {path}'
    assert b'\r\n' in raw and raw.count(b'\n') == raw.count(b'\r\n'), f'Expected original CRLF: {path}'
    with path.open(encoding='utf8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        assert all(None not in r and None not in r.values() for r in rows), f'Malformed TSV: {path}'
    return reader.fieldnames, rows, raw


def key(row):
    return row['AA'], row['UniProt_accession']


def identity(row):
    return '|'.join(row[k] for k in ['AA', 'PDB_ID', 'ligand_instance', 'protein_chain', 'UniProt_accession'])


def ranking(row):
    return (float(row['resolution_A']) if row['resolution_A'] else math.inf,
            -int(row['contact_residue_count'] or 0), row['PDB_ID'], row['ligand_instance'], row['protein_chain'])


def infer(values):
    values = [v for v in values if v != '']
    if not values: return 'empty/string'
    if all(re.fullmatch(r'-?\d+', v) for v in values): return 'integer'
    try:
        for v in values: float(v)
        return 'number'
    except ValueError: pass
    if all(v.startswith(('{', '[')) for v in values):
        try:
            for v in values: json.loads(v)
            return 'JSON/string'
        except ValueError: pass
    return 'string'


def check_types(path, source_headers, source_rows, output_rows):
    profiles = {}
    for column in source_headers:
        dtype = infer([r[column] for r in source_rows])
        profiles[column] = dict(original_inferred_type=dtype, original_blank_allowed=any(r[column] == '' for r in source_rows))
        for row in output_rows:
            value = row[column]
            if value == '':
                assert profiles[column]['original_blank_allowed'], (path, column, 'new missing value')
            elif dtype == 'integer':
                assert re.fullmatch(r'-?\d+', value), (path, column, value)
            elif dtype == 'number':
                assert math.isfinite(float(value)), (path, column, value)
            elif dtype == 'JSON/string':
                json.loads(value)
            assert '\x00' not in value
    return profiles


def main():
    manifest = json.loads((ROOT / 'curation_manifest.json').read_text(encoding='utf8'))
    inventory = json.loads((ROOT / 'candidate_inventory.json').read_text(encoding='utf8'))
    report = dict(validated_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  status='PASS', actual_app_executed=False, app_loader_or_parser_executed=False,
                  method='Independent Python standard-library csv/json/schema/source-provenance checks',
                  originals_unchanged=True, files={}, organisms={}, issues=[])
    for relative, info in manifest['input_files'].items():
        assert hashlib.sha256((BASE / relative).read_bytes()).hexdigest() == info['sha256'], f'Original changed: {relative}'
    total_audit = total_source = 0
    for org, folder in ORGS.items():
        source_dir, output_dir = BASE / folder, ROOT / org
        all_headers, all_rows, _ = read(source_dir / 'all_AA_PDB_complexes.tsv')
        source_index = collections.defaultdict(list)
        for row in all_rows: source_index[identity(row)].append(row)
        _, audit, _ = read(output_dir / 'biological_curation.tsv')
        decisions = {key(r): r for r in audit}
        expected = {(r['AA'], r['accession']) for r in inventory if r['organism'] == org}
        assert len(decisions) == len(audit) and set(decisions) == expected, (org, 'audit coverage')
        assert set(r['classification'] for r in audit) <= {'GOLD', 'SILVER', 'STRUCTURAL_ONLY', 'EXCLUDE', 'MANUAL_REVIEW'}
        covered = []
        alternatives = {}
        for pair, record in decisions.items():
            assert record['reason'] and record['reference_urls'] and record['evidence_summary'], (org, pair, 'missing evidence audit')
            assert pair[0] in CANONICAL
            alts = json.loads(record['structural_alternatives_json'])
            assert len(alts) == int(record['source_all_row_count'])
            alternatives[pair] = alts
            for alt in alts:
                number = alt['source_row_number']
                raw = all_rows[number - 2]
                assert key(raw) == pair
                for column in ('PDB_ID', 'ligand_instance', 'protein_chain', 'resolution_A', 'contact_residue_count',
                               'control_eligible', 'strict_control_eligible', 'protein_mutations',
                               'protein_mutations_current', 'ligand_has_metal_coordination_current',
                               'ligand_contacting_protein_chain_count'):
                    assert alt[column] == raw[column], (org, pair, column, 'alternative differs from original')
                assert (output_dir / alt['source_file']).resolve() == (source_dir / 'all_AA_PDB_complexes.tsv').resolve()
                covered.append(number)
            json.loads(record['PDB_publications_json'])
            json.loads(record['BioLiP2_support_json'])
            json.loads(record['UniProt_functional_annotations_json'])
        assert sorted(covered) == list(range(2, len(all_rows) + 2)), (org, 'lost or double counted raw rows')
        curated_by_file = {}
        for filename in FILES:
            original_headers, original_rows, original_bytes = read(source_dir / filename)
            headers, rows, data = read(output_dir / filename)
            assert headers == original_headers
            assert data.split(b'\r\n', 1)[0] == original_bytes.split(b'\r\n', 1)[0], 'Header bytes/order changed'
            assert len(headers) == (89 if filename == FILES[0] else 91)
            profiles = check_types(output_dir / filename, original_headers, original_rows, rows)
            original_keys = [key(r) for r in original_rows]
            output_keys = [key(r) for r in rows]
            # Inputs inspected separately; both happen to use representatives here.
            # Do not impose this assumption on another file without inspecting it.
            assert len(set(original_keys)) == len(original_keys)
            assert len(set(output_keys)) == len(output_keys)
            assert set(output_keys).issubset(original_keys)
            assert output_keys == [k for k in original_keys if k in output_keys], 'Original ordering lost'
            assert len({identity(r) for r in rows}) == len(rows), 'Duplicate instance rows'
            curated_by_file[filename] = {key(r): r for r in rows}
            expected_keys = set()
            for pair, record in decisions.items():
                flag = 'retained_nonredundant_alternative' if filename == FILES[0] else 'retained_strict_alternative'
                if record['classification'] == 'GOLD' and any(a[flag] for a in alternatives[pair]): expected_keys.add(pair)
            assert set(output_keys) == expected_keys, (org, filename, 'missing GOLD eligible pair')
            source_value_checks = 0
            for row in rows:
                assert decisions[key(row)]['classification'] == 'GOLD' and decisions[key(row)]['target_organism_match'] == '1'
                matches = source_index[identity(row)]
                assert matches, (org, filename, 'new structural row invented')
                matched = False
                for source in matches:
                    unchanged = True
                    for column in all_headers:
                        if column in PATHS and source[column]:
                            assert (output_dir / row[column]).resolve() == (source_dir / source[column]).resolve(), 'Raw path points elsewhere'
                            assert (output_dir / row[column]).resolve().is_file(), 'Broken raw path'
                        elif row[column] != source[column]: unchanged = False
                    if unchanged: matched = True
                assert matched, (org, identity(row), 'Original structural values changed')
                source_value_checks += len(all_headers)
                assert row['control_eligible'] == '1'
                alts = [a for a in alternatives[key(row)] if a['retained_nonredundant_alternative' if filename == FILES[0] else 'retained_strict_alternative']]
                raw_group = [all_rows[a['source_row_number'] - 2] for a in alts]
                assert identity(row) == identity(min(raw_group, key=ranking)), 'Wrong representative ranking'
                count_prefix = 'candidate' if filename == FILES[0] else 'strict_candidate'
                assert int(row[count_prefix + '_structure_count']) == len({a['PDB_ID'] for a in alts})
                assert int(row[count_prefix + '_row_count']) == len(alts)
                if filename == FILES[0]:
                    assert row['representative_selection'] == 'lowest_available_resolution_then_most_contacts_then_identifiers'
                else:
                    assert row['control_set'] == 'strict'
                    for flag in ['strict_control_eligible', 'strict_WT_pass', 'strict_no_metal_pass', 'strict_single_protein_pass']:
                        assert row[flag] == '1', (org, identity(row), flag)
                    assert row['protein_mutations'] == row['protein_mutations_current'] == ''
                    for count in ['protein_mutation_count', 'protein_sequence_conflict_count', 'protein_deletion_count', 'protein_insertion_count', 'ligand_metal_connection_count']:
                        assert int(row[count] or 0) == 0
                    assert row['ligand_has_metal_coordination_current'] == '0'
                    assert int(row['ligand_contacting_protein_chain_count']) == 1
                    assert all(not a['biological_strict_veto'] for a in alts), 'Literature metal veto ignored'
                audit_field = 'corrected_nonredundant_instance' if filename == FILES[0] else 'corrected_strict_instance'
                assert decisions[key(row)][audit_field] == identity(row)
            report['files'][f'{org}/{filename}'] = dict(rows=len(rows), original_rows=len(original_rows),
                columns=len(headers), header_bytes_match=True, TSV_format_and_CRLF_match=True,
                original_structural_values_match=True, original_structural_values_checked=source_value_checks,
                duplicate_pairs=0, duplicate_instance_rows=0, all_pairs_GOLD=True,
                original_selection_semantics_preserved=True, numeric_JSON_and_blank_profiles_match=True,
                raw_provenance_paths_resolve=True, data_type_profiles=profiles,
                sha256=hashlib.sha256(data).hexdigest())
        # previous_representative means the preceding NONRED stage in this cohort.
        for pair, row in curated_by_file[FILES[1]].items():
            previous = curated_by_file[FILES[0]][pair]
            assert row['previous_representative_PDB_ID'] == previous['PDB_ID']
            assert int(row['representative_changed_after_strict_filter']) == int(identity(row) != identity(previous))
        report['organisms'][org] = dict(audit_pairs=len(audit), raw_rows_preserved=len(all_rows),
            GOLD_pairs=sum(r['classification'] == 'GOLD' for r in audit),
            strict_rows=len(curated_by_file[FILES[1]]), nonredundant_rows=len(curated_by_file[FILES[0]]),
            all_candidates_have_classification_reason_and_references=True)
        total_audit += len(audit)
        total_source += len(all_rows)
    assert total_audit == manifest['candidate_pairs_evaluated'] == len(inventory)
    assert total_source == manifest['source_rows_preserved_in_originals_and_audit']
    report.update(candidate_pairs_evaluated=total_audit, original_raw_rows_preserved=total_source)
    (ROOT / 'validation_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf8')
    print(json.dumps({k: v for k, v in report.items() if k != 'files'}, indent=2))


if __name__ == '__main__':
    main()
