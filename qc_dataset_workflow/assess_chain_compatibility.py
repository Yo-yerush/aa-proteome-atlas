# Adapted from the existing reviewed workflow for fresh datasets.
# Original source SHA256: 543509aa63a8a2a74fbfafb3d83350dad358f05f0738421b465952fb5d2f6f3a
"""Assess curated ligand sites against an isolated chain; never imports an app.

Reads existing GOLD tables/evidence. A focused coordinate audit includes biological
symmetry mates and nonprotein polymers, which the original contact count can miss.
No docking or original structural pipeline is run, and app TSVs are not rewritten.
"""
import collections
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import itertools
import json
import os
import pathlib
import re
import time

import numpy as np
import requests
from scipy.spatial import cKDTree

BASE = pathlib.Path(os.environ["AA_QC_DATA_ROOT"]).resolve()
ROOT = BASE / 'curated_controls'
CACHE = ROOT / 'evidence_cache' / 'chain_compatibility'
COORDS = CACHE / 'coordinates'
ORGS = tuple(os.environ.get('AA_QC_ORGANISMS', 'ecoli,arabidopsis,human,mouse').split(','))
INPUTS = ('strict_WT_single_protein_AA_controls.tsv', 'nonredundant_AA_protein_controls.tsv')


def load_controls():
    rows = {}
    for org in ORGS:
        for name in INPUTS:
            with (ROOT / org / name).open(encoding='utf8', newline='') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    identity = (org, row['AA'], row['PDB_ID'], row['ligand_instance'], row['protein_chain'], row['UniProt_accession'])
                    if identity not in rows:
                        rows[identity] = dict(row, organism_cohort=org, source_files=[])
                    rows[identity]['source_files'].append(name)
    if os.environ.get('AA_QC_INCLUDE_ALTERNATIVES') == '1':
        from curate_biological_controls import read_table, site_decision, strict_biological_veto, rebase_paths
        folders = {'ecoli': 'E_coli', 'arabidopsis': 'Arabidopsis', 'human': 'Human', 'mouse': 'Mouse', 'yeast': 'Yeast'}
        for org in ORGS:
            _, audit = read_table(ROOT / org / 'biological_curation.tsv')
            accepted = {(r['AA'], r['UniProt_accession']) for r in audit if r['classification'] == 'GOLD'}
            original = BASE / folders[org]
            _, candidates = read_table(original / 'all_AA_PDB_complexes.tsv')
            for candidate in candidates:
                if (candidate['AA'], candidate['UniProt_accession']) not in accepted:
                    continue
                if candidate['strict_control_eligible'] != '1' or strict_biological_veto(candidate):
                    continue
                if site_decision(candidate, 'GOLD')[0] != 'GOLD':
                    continue
                identity = (org, candidate['AA'], candidate['PDB_ID'], candidate['ligand_instance'], candidate['protein_chain'], candidate['UniProt_accession'])
                if identity not in rows:
                    row = rebase_paths(candidate, original, ROOT / org)
                    rows[identity] = dict(row, organism_cohort=org, source_files=['all_AA_PDB_complexes.tsv'])
    return list(rows.values())


def fetch_coordinates(pdb):
    path = COORDS / (pdb + '.cif.gz')
    url = f'https://files.rcsb.org/download/{pdb}.cif.gz'
    if path.is_file():
        return dict(PDB_ID=pdb, url=url, path=str(path.relative_to(ROOT)), cached=True,
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    for attempt in range(3):
        try:
            response = requests.get(url, timeout=(15, 60))
            response.raise_for_status()
            assert response.content[:2] == b'\x1f\x8b', 'Expected gzip coordinate file'
            with gzip.GzipFile(fileobj=__import__('io').BytesIO(response.content)) as stream:
                assert stream.read(64).startswith(b'data_'), 'Expected mmCIF data block'
            path.write_bytes(response.content)
            return dict(PDB_ID=pdb, url=url, path=str(path.relative_to(ROOT)), cached=False,
                        bytes=len(response.content), last_modified=response.headers.get('Last-Modified', ''),
                        retrieved_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        sha256=hashlib.sha256(response.content).hexdigest())
        except Exception:
            if attempt == 2: raise
            time.sleep(attempt + 1)


def tokens(lines):
    """CIF 1.1 quoting/comments, including semicolon-delimited text fields."""
    iterator = iter(lines)
    for line in iterator:
        if line.startswith(';'):
            content = [line[1:].rstrip('\r\n')]
            for other in iterator:
                if other.startswith(';'): break
                content.append(other.rstrip('\r\n'))
            else: raise ValueError('Unterminated CIF text field')
            yield '\n'.join(content)
            continue
        i, n = 0, len(line)
        while i < n:
            while i < n and line[i].isspace(): i += 1
            if i == n or line[i] == '#': break
            begin = i
            if line[i] in "\"'":
                quote = line[i]
                i += 1
                begin = i
                while i < n:
                    if line[i] == quote and (i + 1 == n or line[i+1].isspace()): break
                    i += 1
                if i == n: raise ValueError('Unterminated CIF quoted field')
                yield line[begin:i]
                i += 1
            else:
                while i < n and not line[i].isspace(): i += 1
                yield line[begin:i]


def read_categories(path):
    needed = {'_atom_site', '_entity_poly', '_entity', '_struct_asym',
              '_pdbx_struct_assembly', '_pdbx_struct_assembly_gen', '_pdbx_struct_oper_list',
              '_pdbx_entity_nonpoly', '_struct_conn'}
    result, pending = collections.defaultdict(list), None
    with gzip.open(path, 'rt', encoding='utf8') as handle:
        stream = iter(tokens(handle))
        while True:
            try: token = pending if pending is not None else next(stream)
            except StopIteration: break
            pending = None
            if token == 'loop_':
                headers = []
                for token in stream:
                    if token.startswith('_'): headers.append(token)
                    else: pending = token; break
                if not headers: raise ValueError('Empty CIF loop')
                category = headers[0].split('.')[0]
                row = []
                while True:
                    try: token = pending if pending is not None else next(stream)
                    except StopIteration:
                        if row: raise ValueError('Incomplete CIF loop row')
                        break
                    pending = None
                    if not row and (token == 'loop_' or token == 'stop_' or token.startswith(('_', 'data_', 'save_'))):
                        pending = token
                        break
                    row.append(token)
                    if len(row) == len(headers):
                        if category in needed:
                            result[category].append(dict(zip([h.split('.', 1)[1] for h in headers], row)))
                        row = []
            elif token.startswith('_'):
                value = next(stream)
                category, name = token.split('.', 1)
                if category in needed:
                    if not result[category]: result[category].append({})
                    result[category][0][name] = value
    return result


def expand_operations(expression):
    expression = expression.replace(' ', '')
    groups = re.findall(r'\(([^()]*)\)', expression) if '(' in expression else [expression]
    if '(' in expression and ''.join('('+g+')' for g in groups) != expression:
        raise ValueError('Malformed operation product: '+expression)
    parts = []
    for group in groups:
        values = []
        for item in group.split(','):
            match = re.fullmatch(r'(\d+)-(\d+)', item)
            if match:
                start, end = map(int, match.groups())
                assert 0 <= end - start <= 10000
                values.extend(str(i) for i in range(start, end+1))
            else:
                assert re.fullmatch(r'[A-Za-z0-9_]+', item), expression
                values.append(item)
        parts.append(values)
    return list(itertools.product(*parts))


def operation_matrices(categories):
    return {op['id']: (np.array([[float(op[f'matrix[{i}][{j}]']) for j in range(1, 4)] for i in range(1, 4)]),
                        np.array([float(op[f'vector[{i}]']) for i in range(1, 4)]))
            for op in categories['_pdbx_struct_oper_list']}


def composed_operation(ids, operations):
    matrix, vector = np.eye(3), np.zeros(3)
    # Product (X)(Y) applies Y first, then X, using column-vector convention.
    for identifier in reversed(ids):
        m, v = operations[identifier]
        matrix, vector = m @ matrix, m @ vector + v
    return matrix, vector


def structured_atoms(categories):
    polymers = {p['entity_id']: p['type'] for p in categories['_entity_poly']}
    atoms = collections.defaultdict(list)
    model_ids = sorted({int(a.get('pdbx_PDB_model_num', '1')) for a in categories['_atom_site']})
    model = model_ids[0]
    for atom in categories['_atom_site']:
        if int(atom.get('pdbx_PDB_model_num', '1')) != model: continue
        if atom['type_symbol'].upper() in ('H', 'D'): continue
        if float(atom.get('occupancy', '1')) <= 0: continue
        if atom.get('label_asym_id') in ('.', '?'): continue
        atoms[atom['label_asym_id']].append(atom)
    chains = {}
    for label, group in atoms.items():
        eid = group[0]['label_entity_id']
        chains[label] = dict(atoms=group, xyz=np.array([[float(a['Cartn_x']), float(a['Cartn_y']), float(a['Cartn_z'])] for a in group]),
                             entity_id=eid, polymer_type=polymers.get(eid, ''), auth_chain=group[0]['auth_asym_id'])
    return chains, model, model_ids


def contacts(ligand_xyz, chain_xyz, atoms, cutoff=5.0):
    # Nearest ligand heavy-atom distance for every polymer heavy atom.
    lower, upper = ligand_xyz.min(0) - cutoff, ligand_xyz.max(0) + cutoff
    candidate = np.flatnonzero(((chain_xyz >= lower) & (chain_xyz <= upper)).all(1))
    if len(candidate) == 0: return None
    distances, lig_indices = cKDTree(ligand_xyz).query(chain_xyz[candidate], k=1)
    close = np.flatnonzero(distances <= cutoff + 1e-8)
    if len(close) == 0: return None
    residues = {}
    for j in close:
        a = atoms[candidate[j]]
        key = a['auth_seq_id']+a.get('pdbx_PDB_ins_code', '').replace('?', '').replace('.', '')
        record = residues.setdefault(key, dict(auth_seq_id=a['auth_seq_id'], label_seq_id=a.get('label_seq_id', ''),
                                               comp_id=a['label_comp_id'], minimum_distance_A=99., contacting_atoms=[], ligand_atom_indices=[]))
        record['minimum_distance_A'] = min(record['minimum_distance_A'], float(distances[j]))
        record['contacting_atoms'].append(a['label_atom_id'])
        record['ligand_atom_indices'].append(int(lig_indices[j]))
    ordered = sorted(residues.values(), key=lambda r: (r['minimum_distance_A'], r['auth_seq_id']))
    return dict(minimum_distance_A=float(distances[close].min()), residue_count_5A=len(residues),
                residue_count_4A=sum(r['minimum_distance_A'] <= 4. for r in ordered),
                residue_count_3_5A=sum(r['minimum_distance_A'] <= 3.5 for r in ordered),
                ligand_atom_indices_contacted=sorted(set(int(lig_indices[j]) for j in close)), residues=ordered)


def audit_site(row, categories, chains, operations, model, model_ids):
    ligand = chains.get(row['ligand_label_asym_id'])
    selected = chains.get(row['protein_chain'])
    if not ligand or not selected: raise ValueError('Selected chain or ligand missing from coordinates')
    ligand_atoms = [i for i, a in enumerate(ligand['atoms']) if a['label_comp_id'] == row['AA'] and a['auth_seq_id'] == row['ligand_auth_seq_id']]
    if not ligand_atoms: raise ValueError('Exact canonical ligand residue not found')
    ligand_xyz = ligand['xyz'][ligand_atoms]
    ligand_atom_names = [ligand['atoms'][i]['label_atom_id'] for i in ligand_atoms]
    asu_contacts = []
    for chain, structure in chains.items():
        if not structure['polymer_type']: continue
        found = contacts(ligand_xyz, structure['xyz'], structure['atoms'])
        if found:
            found.update(label_asym_id=chain, auth_asym_id=structure['auth_chain'], entity_id=structure['entity_id'],
                         polymer_type=structure['polymer_type'], selected_chain=chain == row['protein_chain'],
                         ligand_atoms_contacted=[ligand_atom_names[i] for i in found['ligand_atom_indices_contacted']])
            asu_contacts.append(found)
    assemblies = []
    descriptions = {a['id']: a for a in categories['_pdbx_struct_assembly']}
    generators = categories['_pdbx_struct_assembly_gen']
    for assembly_id, description in descriptions.items():
        generated = collections.defaultdict(set)
        for generator in generators:
            if generator['assembly_id'] != assembly_id: continue
            ops = expand_operations(generator['oper_expression'])
            for chain in generator['asym_id_list'].split(','):
                generated[chain].update(ops)
        anchors = generated[row['protein_chain']] & generated[row['ligand_label_asym_id']]
        if not anchors: continue
        # All copies under the same assembly operators are symmetry-equivalent.
        anchor = min(anchors, key=lambda ids: (not np.allclose(composed_operation(ids, operations)[0], np.eye(3)), ids))
        ligand_m, ligand_v = composed_operation(anchor, operations)
        transformed_ligand = ligand_xyz @ ligand_m.T + ligand_v
        found_contacts = []
        for chain, ids_set in generated.items():
            if chain not in chains or not chains[chain]['polymer_type']: continue
            structure = chains[chain]
            for ids in sorted(ids_set):
                matrix, vector = composed_operation(ids, operations)
                transformed = structure['xyz'] @ matrix.T + vector
                found = contacts(transformed_ligand, transformed, structure['atoms'])
                if found:
                    found.update(label_asym_id=chain, auth_asym_id=structure['auth_chain'], operation_ids=list(ids),
                        entity_id=structure['entity_id'], polymer_type=structure['polymer_type'],
                        selected_chain=chain == row['protein_chain'] and np.allclose(matrix, ligand_m, atol=1e-6) and np.allclose(vector, ligand_v, atol=1e-6),
                        ligand_atoms_contacted=[ligand_atom_names[i] for i in found['ligand_atom_indices_contacted']])
                    found_contacts.append(found)
        # Some deposited assembly ids contain only the AU; retain all for audit,
        # let the classification distinguish author/PISA candidate descriptions.
        assemblies.append(dict(assembly_id=assembly_id, details=description.get('details', ''),
            method_details=description.get('method_details', ''), oligomeric_details=description.get('oligomeric_details', ''),
            ligand_anchor_operation=list(anchor), ligand_atom_count=len(ligand_atoms), polymer_contacts=found_contacts))
    return dict(PDB_ID=row['PDB_ID'], ligand_instance=row['ligand_instance'], AA=row['AA'],
        protein_chain=row['protein_chain'], selected_model=model, available_models=model_ids,
        heavy_atom_count=len(ligand_atoms), ligand_atom_names=ligand_atom_names,
        alt_ids=sorted({ligand['atoms'][i].get('label_alt_id', '.') for i in ligand_atoms}),
        asymmetric_unit_contacts=asu_contacts, biological_assembly_contacts=assemblies)


def cached_api_context(row):
    path = (ROOT / row['organism_cohort'] / row['raw_response_file']).resolve()
    data = json.loads(path.read_text(encoding='utf8'))['response']['data']['entries']
    entry = next(e for e in data if e['rcsb_id'] == row['PDB_ID'])
    lig = next(l for e in entry['nonpolymer_entities'] for l in e['nonpolymer_entity_instances'] if l['rcsb_id'] == row['ligand_instance'])
    polymers = {}
    for entity in entry['polymer_entities']:
        for chain in entity['polymer_entity_instances']:
            ids = chain['rcsb_polymer_entity_instance_container_identifiers']
            polymers[ids['asym_id']] = dict(entity_id=entity['rcsb_id'].rsplit('_', 1)[-1],
                polymer_type=entity['entity_poly']['rcsb_entity_polymer_type'],
                description=entity['rcsb_polymer_entity'].get('pdbx_description', ''),
                UniProt_ids=entity['rcsb_polymer_entity_container_identifiers'].get('uniprot_ids') or [])
    return dict(polymer_chains=polymers,
                API_neighbors=[n for n in lig.get('rcsb_target_neighbors', []) if float(n['distance']) <= 5.0])


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    COORDS.mkdir(parents=True, exist_ok=True)
    rows = load_controls()
    hashes = {f'{org}/{name}': hashlib.sha256((ROOT/org/name).read_bytes()).hexdigest() for org in ORGS for name in INPUTS}
    (CACHE/'input_hashes.json').write_text(json.dumps(hashes, indent=2)+'\n', encoding='utf8')
    pdbs = sorted({r['PDB_ID'] for r in rows})
    status, errors = [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:
        pending = {pool.submit(fetch_coordinates, pdb): pdb for pdb in pdbs}
        for i, future in enumerate(concurrent.futures.as_completed(pending), 1):
            try: status.append(future.result())
            except Exception as error: errors.append(dict(PDB_ID=pending[future], error=str(error)))
            if i % 20 == 0: print('Coordinate audit downloads', i, '/', len(pdbs), 'errors', len(errors), flush=True)
    (CACHE/'coordinate_sources.json').write_text(json.dumps(status, indent=2)+'\n', encoding='utf8')
    output, analysis_errors = [], []
    grouped = collections.defaultdict(list)
    for row in rows: grouped[row['PDB_ID']].append(row)
    for i, (pdb, sites) in enumerate(sorted(grouped.items()), 1):
        try:
            categories = read_categories(COORDS/(pdb+'.cif.gz'))
            chains, model, model_ids = structured_atoms(categories)
            operations = operation_matrices(categories)
            for row in sites:
                site = audit_site(row, categories, chains, operations, model, model_ids)
                site.update(organism_cohort=row['organism_cohort'], UniProt_accession=row['UniProt_accession'], source_files=row['source_files'])
                site['original_API_context'] = cached_api_context(row)
                output.append(site)
        except Exception as error:
            analysis_errors.append(dict(PDB_ID=pdb, error=str(error)))
        if i % 20 == 0: print('Sites/assembly coordinates audited', i, '/',len(pdbs), 'errors', len(analysis_errors), flush=True)
    (CACHE/'site_geometry.json').write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n', encoding='utf8')
    (CACHE/'geometry_errors.json').write_text(json.dumps(errors+analysis_errors, indent=2)+'\n', encoding='utf8')
    assert all(hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == value for path,value in hashes.items()), 'App control TSV changed'
    print('DONE',len(output),'sites;',len(errors+analysis_errors),'errors',flush=True)


if __name__ == '__main__': main()
