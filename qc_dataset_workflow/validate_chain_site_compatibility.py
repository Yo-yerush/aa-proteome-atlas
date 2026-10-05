"""Validate fresh compatibility audits independently, without fixed dataset counts."""
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path

ROOT = Path(os.environ['AA_QC_DATA_ROOT']).resolve() / 'curated_controls'
CACHE = ROOT / 'evidence_cache' / 'chain_compatibility'
ORGS = tuple(os.environ.get('AA_QC_ORGANISMS', 'ecoli,arabidopsis,human,mouse').split(','))
INPUTS = ('strict_WT_single_protein_AA_controls.tsv', 'nonredundant_AA_protein_controls.tsv')
FIELDS = ('AA', 'PDB_ID', 'ligand_instance', 'protein_chain', 'UniProt_accession')
LABELS = {'monomer_compatible', 'complex_dependent', 'uncertain'}


def read(path):
    with path.open(encoding='utf8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        assert len(reader.fieldnames) == len(set(reader.fieldnames)), path
        rows = list(reader)
        assert all(None not in r and None not in r.values() for r in rows), path
        return rows


def key(org, row):
    return (org,) + tuple(row[f] for f in FIELDS)


def main():
    manifest = json.loads((ROOT / 'chain_compatibility_manifest.json').read_text(encoding='utf8'))
    for relative, expected in manifest['protected_input_sha256'].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, relative
    combined = read(ROOT / 'binding_site_chain_compatibility_all.tsv')
    observed = {key(r['organism'], r): r for r in combined}
    assert len(observed) == len(combined) == manifest['unique_sites']
    source_count, expected_pairs, reports = 0, set(), {}
    for org in ORGS:
        bio = {(r['AA'], r['UniProt_accession']): r for r in read(ROOT / org / 'biological_curation.tsv')}
        local = read(ROOT / org / 'binding_site_chain_compatibility.tsv')
        assert local == [r for r in combined if r['organism'] == org]
        indices = {key(org, r): r for r in local}
        joins = 0
        for name in INPUTS:
            for ordinal, row in enumerate(read(ROOT / org / name), 2):
                found = indices[key(org, row)]
                assert name in found['source_files'].split(';')
                assert ordinal in json.loads(found['source_row_numbers_json'])[name]
                assert all(found[f] == row[f] for f in FIELDS)
                assert bio[row['AA'], row['UniProt_accession']]['classification'] == 'GOLD'
                expected_pairs.add((org, row['AA'], row['UniProt_accession']))
                joins += 1
        assert joins == sum(len(numbers) for r in local for numbers in json.loads(r['source_row_numbers_json']).values())
        source_count += joins
        for row in local:
            assert row['chain_compatibility'] in LABELS
            assert row['reason'] and row['compatibility_reference_urls']
            assert row['biological_classification'] == 'GOLD'
            for column, value in row.items():
                if column.endswith('_json'):
                    json.loads(value)
            for column in ('ligand_unique_occupied_heavy_atom_count', 'ligand_expected_heavy_atom_count'):
                assert int(row[column]) >= 0
            if row['other_polymer_minimum_distance_A']:
                assert 0 < float(row['other_polymer_minimum_distance_A']) <= 5
        reports[org] = dict(unique_sites=len(local), source_rows=joins,
                            classifications=dict(collections.Counter(r['chain_compatibility'] for r in local)))
    assert source_count == manifest['source_rows_covered']
    pairs = read(ROOT / 'binding_site_chain_compatibility_pairs.tsv')
    assert len(pairs) == len(expected_pairs) == manifest['unique_GOLD_AA_UniProt_pairs']
    assert {(r['organism'], r['AA'], r['UniProt_accession']) for r in pairs} == expected_pairs
    geometry = json.loads((CACHE / 'site_geometry.json').read_text(encoding='utf8'))
    max_delta = 0.0
    for site in geometry:
        self_asu = [c for c in site['asymmetric_unit_contacts'] if c['selected_chain']]
        assert len(self_asu) == 1
        for assembly in site['biological_assembly_contacts']:
            selected = [c for c in assembly['polymer_contacts'] if c['selected_chain']]
            assert len(selected) == 1
            delta = abs(selected[0]['minimum_distance_A'] - self_asu[0]['minimum_distance_A'])
            assert math.isfinite(delta) and delta < 0.0001, (site['PDB_ID'], delta)
            max_delta = max(max_delta, delta)
    sources = json.loads((CACHE / 'coordinate_sources.json').read_text(encoding='utf8'))
    assert len(sources) == len({r['PDB_ID'] for r in sources})
    for record in sources:
        path = ROOT / record['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']
    report = dict(status='PASS', selected_sites=len(combined), selected_pairs=len(pairs), source_rows=source_count,
                  coordinate_files_verified=len(sources), all_alternative_geometry_sites=len(geometry),
                  maximum_transform_distance_error_A=max_delta, organisms=reports,
                  app_launched=False, app_loader_executed=False)
    (ROOT / 'chain_compatibility_validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf8')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
