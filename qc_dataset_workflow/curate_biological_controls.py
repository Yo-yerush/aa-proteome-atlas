# Adapted from the existing reviewed workflow for fresh datasets.
# Original source SHA256: 9e60ba7b70f49beefbc851dd3285b742eab02ca1ebae0692226b2e09b092bb2d
"""Curate existing AA controls from cached biological evidence; no app/pipeline imports.

This stage filters biochemical identities and adjudicates two site-specific cases.
Coordinates, docking, contacts, mutations, metal and assembly annotations are reused.
"""
import collections
import copy
import csv
import datetime
import hashlib
import json
import math
import os
import pathlib
import re

from biological_decisions import DECISIONS

BASE = pathlib.Path(os.environ["AA_QC_DATA_ROOT"]).resolve()
OUT = BASE / 'curated_controls'
CACHE = OUT / 'evidence_cache'
ORGANISMS = {'ecoli': 'E_coli', 'arabidopsis': 'Arabidopsis',
             'human': 'Human', 'mouse': 'Mouse', 'yeast': 'Yeast'}
ORGANISMS = {key: value for key, value in ORGANISMS.items() if key in os.environ.get("AA_QC_ORGANISMS", "ecoli,arabidopsis,human,mouse").split(",")}
NONRED = 'nonredundant_AA_protein_controls.tsv'
STRICT = 'strict_WT_single_protein_AA_controls.tsv'
ALL = 'all_AA_PDB_complexes.tsv'
PATH_COLUMNS = ('raw_response_file', 'strict_metadata_response_file')
AA_NAMES = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
                   'alanine arginine asparagine aspartate cysteine glutamine glutamate glycine histidine isoleucine leucine lysine methionine phenylalanine proline serine threonine tryptophan tyrosine valine'.split()))
AUDIT_COLUMNS = ['organism', 'AA', 'UniProt_accession', 'protein_name', 'classification',
    'functional_role', 'reason', 'decision_basis', 'source_organism', 'target_organism_match',
    'evidence_summary', 'specific_primary_PMIDs', 'UniProt_function_PMIDs',
    'reference_urls', 'UniProt_functional_annotations_json', 'UniProt_AA_binding_sites_json',
    'EC_numbers', 'KEGG_cached_EC_numbers', 'auxiliary_database_review',
    'PDB_IDs', 'PDB_publications_json', 'PDB_publication_count', 'available_abstract_count',
    'missing_abstract_PMIDs', 'BioLiP2_matched_PDB_IDs', 'BioLiP2_annotation_count',
    'BioLiP2_support_json', 'source_files', 'source_all_row_count',
    'structural_alternatives_json', 'original_nonredundant_instance',
    'corrected_nonredundant_instance', 'original_strict_instance', 'corrected_strict_instance',
    'nonredundant_retained', 'strict_retained', 'representative_changes',
    'evidence_review_scope', 'limitations', 'curation_date_UTC']


def compact(value):
    return json.dumps(value, ensure_ascii=False, separators=(',', ':'))


def read_table(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
        assert all(None not in row and None not in row.values() for row in rows), path
        return reader.fieldnames, rows


def pair(row):
    return row['AA'], row['UniProt_accession']


def row_key(row):
    return row['AA'], row['PDB_ID'], row['ligand_instance'], row['protein_chain'], row['UniProt_accession']


def identifier(row):
    return '' if row is None else '|'.join(row_key(row))


def rank(row):
    # These are the independently inspected ranking rules of both originals.
    resolution = float(row['resolution_A']) if row['resolution_A'] else math.inf
    return resolution, -int(row['contact_residue_count'] or 0), row['PDB_ID'], row['ligand_instance'], row['protein_chain']


def target_match(organism, source):
    if organism == 'ecoli':
        return source.startswith('Escherichia coli')
    if organism == 'arabidopsis':
        return source.startswith('Arabidopsis')
    if organism == 'yeast':
        return source == 'Saccharomyces cerevisiae' or source.startswith('Saccharomyces cerevisiae ')
    return source == {'human': 'Homo sapiens', 'mouse': 'Mus musculus'}[organism]


def functional_comments(uniprot):
    return [c for c in uniprot.get('comments', []) if c['commentType'] in
            ('FUNCTION', 'CATALYTIC ACTIVITY', 'ACTIVITY REGULATION', 'BIOPHYSICOCHEMICAL PROPERTIES')]


def pmids_in(value):
    result = set()
    def walk(v):
        if isinstance(v, dict):
            if v.get('source') == 'PubMed' and v.get('id'):
                result.add(str(v['id']))
            for item in v.values(): walk(item)
        elif isinstance(v, list):
            for item in v: walk(item)
    walk(value)
    return sorted(result)


def aa_features(uniprot, aa):
    name = AA_NAMES[aa]
    alternatives = {name}
    if aa == 'ASP': alternatives.add('aspartic acid')
    if aa == 'GLU': alternatives.add('glutamic acid')
    result = []
    for feature in uniprot.get('features', []):
        if feature['type'] != 'Binding site': continue
        ligand = feature.get('ligand', {}).get('name', '').lower()
        if any(re.fullmatch(r'(?:l-|\(s\)-)?' + re.escape(n), ligand) for n in alternatives):
            result.append(feature)
    return result


def decide(record, uniprot):
    accession, aa = record['accession'], record['AA']
    if not accession:
        return dict(classification='MANUAL_REVIEW', functional_role='unmapped protein identity',
                    reason='The original hit has no mapped UniProt protein identity. A functional free-AA pair cannot be established or safely attributed.',
                    specific_pmids=[], decision_basis='Unresolved original accession mapping')
    source = uniprot.get('organism', {}).get('scientificName', '')
    if not source:
        return dict(classification='MANUAL_REVIEW', functional_role='missing accession evidence', reason='Fresh UniProt evidence is unavailable; do not infer functional binding or organism identity.', specific_pmids=[], decision_basis='Evidence retrieval unresolved')
    if not target_match(record['organism'], source):
        return dict(classification='EXCLUDE', functional_role='foreign organism/fusion component',
                    reason=f'The mapped protein is from {source}, outside this organism cohort. A multispecies/fusion deposition must not transfer ligand function to this accession.',
                    specific_pmids=[], decision_basis='UniProt accession source-organism check')
    if (accession, aa) in DECISIONS:
        result = copy.deepcopy(DECISIONS[accession, aa])
        result['decision_basis'] = 'Included prior accession-AA review, corroborated with freshly retrieved annotations and citation context'
        return result
    # Never infer GOLD from protein names, EC families, PDB presence or BioLiP2.
    function = ' '.join(t.get('value', '') for c in functional_comments(uniprot)
                        for t in c.get('texts', []))
    reaction = '; '.join(c.get('reaction', {}).get('name', '') for c in functional_comments(uniprot)
                         if c['commentType'] == 'CATALYTIC ACTIVITY')
    context = function or reaction
    if context:
        return dict(classification='MANUAL_REVIEW', functional_role='new pair awaiting biological review',
                    reason=f'The documented function of {record["name"]} and the linked structural publication context do not establish {AA_NAMES[aa]} as a native functional FREE canonical-AA ligand. Structural occupancy alone is insufficient; this is not a claim that binding is impossible.',
                    specific_pmids=[], decision_basis='Conservative review of UniProt function/reaction, PDB citations/available abstracts, and BioLiP2 context')
    return dict(classification='MANUAL_REVIEW', functional_role='insufficient functional annotation',
                reason=f'The deposited {AA_NAMES[aa]} association lacks sufficient protein-specific functional annotation or available primary evidence to establish native free-AA recognition.',
                specific_pmids=[], decision_basis='Unresolved function after reviewing available annotations/publication metadata')


def site_decision(row, classification):
    if classification != 'GOLD': return classification, ''
    if row['UniProt_accession'] == 'P27616' and row['AA'] == 'ASP':
        if row['ligand_instance'] == '2CNU.B':
            return 'GOLD', 'Free ASP at the SAICAR substrate/product cleft contacts Lys19, Ser40, Ala41, Tyr42, Asp43, Lys260 and Arg264; these same residues contact the aspartyl portion of SSS in 2CNV. This is a structural inference corroborated by the free-aspartate reaction and primary substrate-binding study PMID 26072057.'
        if row['ligand_instance'] in {'2CNU.C', '2CNU.D', '2CNV.B', '2CNV.C'}:
            return 'STRUCTURAL_ONLY', 'Peripheral ASP at residues 73/76/78/110 or 183-187 remains present alongside bound SAICAR product in 2CNV. Functional free-aspartate recognition at these additional surface sites is not established; retain only the substrate-cleft 2CNU.B pose as a positive control.'
        return 'MANUAL_REVIEW', 'This additional SAICAR-synthetase ASP pose has not been individually assigned to the validated substrate pocket.'
    accession, aa = row['UniProt_accession'], row['AA']
    if accession == 'Q9P2J5' and aa == 'LEU':
        if row['ligand_instance'] == '6KQY.B':
            return 'GOLD', 'Free leucine occupies the functional synthetic/sensing pocket, including UniProt positions 52/54; PMID 33910001.'
        return 'SILVER', 'Free leucine occupies the CP1 editing pocket, not the validated synthetic/sensing site; do not substitute an editing-site product pose for the functional free-Leu control (PMID 33910001).'
    if accession == 'Q9LYU8' and aa == 'LYS':
        if row['ligand_instance'] in ('2CDQ.C', '2CDQ.G'):
            return 'GOLD', 'The native ACT1 allosteric lysine site spans both protein chains (PMID 16731588); nonredundant only, because strict requires exactly one contacting chain.'
        return 'STRUCTURAL_ONLY', 'The lysine instance occupies the Asp substrate pocket in crystals grown without Asp; the demonstrated physiological feedback site is ACT1. A functional role for this catalytic-site lysine occupancy is not established by PMID 16731588.'
    return 'GOLD', 'Functional ligand context supported at the pair level; original contact/mapping information retained.'


def rebase_paths(row, source_dir, output_dir):
    row = dict(row)
    for column in PATH_COLUMNS:
        if row.get(column):
            absolute = (source_dir / row[column]).resolve()
            assert absolute.is_file(), (column, absolute)
            row[column] = pathlib.Path(os.path.relpath(absolute, output_dir)).as_posix()
    return row


def strict_biological_veto(row):
    # Supplement LINK metadata with explicit primary-literature evidence; no
    # geometry recomputation and no exclusion for merely having a remote metal.
    if pair(row) == ('CYS', 'Q16878'):
        return 'CDO directly coordinates cysteine at its catalytic metal site (PMIDs 17135237, 29942080, 30946568); absence of a reported ligand-metal LINK is not proof of metal-independent binding.'
    if pair(row) == ('CYS', 'P21888'):
        return 'CysRS uses a direct cysteine-thiolate/zinc interaction (PMID 12032090); this is not a metal-independent binding control.'
    return ''


def write_table(path, headers, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, delimiter='\t', lineterminator='\r\n', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def indexes():
    entries, pubs, biolip, kegg = {}, {}, collections.defaultdict(list), {}
    for path in (CACHE / 'rcsb').glob('*.json'):
        for entry in json.loads(path.read_text(encoding='utf8')).get('data', {}).get('entries') or []:
            if entry: entries[entry['rcsb_id']] = entry
    for path in (CACHE / 'publications_v2').glob('*.json'):
        for publication in json.loads(path.read_text(encoding='utf8')).get('resultList', {}).get('result', []):
            if publication.get('source') == 'MED': pubs[str(publication['id'])] = publication
    with (CACHE / 'biolip2_matches.tsv').open(encoding='utf8', newline='') as handle:
        for fields in csv.reader(handle, delimiter='\t'):
            assert len(fields) == 21
            biolip[fields[0].upper(), fields[4]].append(fields)
    for path in (CACHE / 'kegg').glob('*.txt'):
        for entry in path.read_text(encoding='utf8').split('///'):
            match = re.search(r'^ENTRY\s+EC (\S+)', entry, re.M)
            if match: kegg[match.group(1)] = path.name
    return entries, pubs, biolip, kegg


def main():
    records = json.loads((OUT / 'candidate_inventory.json').read_text(encoding='utf8'))
    entries, publications, biolip, kegg = indexes()
    uniprot = {p.stem: json.loads(p.read_text(encoding='utf8')) for p in (CACHE / 'uniprot').glob('*.json')}
    manifest = {'curation_date_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'app_executed': False, 'structural_dataset_built_from_APIs': True,
                'organisms': {}, 'input_files': {}, 'output_files': {},
                'policy': 'GOLD-only, exact original app schemas, original independent structural selection rules after biological/site filtering',
                'decision_map_sha256': hashlib.sha256((pathlib.Path(__file__).resolve().parent / 'biological_decisions.py').read_bytes()).hexdigest()}
    for org, directory in ORGANISMS.items():
        source_dir, output_dir = BASE / directory, OUT / org
        output_dir.mkdir(exist_ok=True)
        sources = {name: read_table(source_dir / name) for name in (ALL, NONRED, STRICT)}
        for name, (headers, rows) in sources.items():
            path = source_dir / name
            manifest['input_files'][f'{directory}/{name}'] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                                                              rows=len(rows), columns=len(headers), header=headers)
        all_headers, all_rows = sources[ALL]
        org_records = [r for r in records if r['organism'] == org]
        decisions = { (r['AA'], r['accession']): decide(r, uniprot.get(r['accession'], {})) for r in org_records }
        assert len(decisions) == len(org_records)
        alternatives = collections.defaultdict(list)
        for number, row in enumerate(all_rows, 2):
            assert pair(row) in decisions, (org, pair(row))
            classification, reason = site_decision(row, decisions[pair(row)]['classification'])
            alternatives[pair(row)].append((number, row, classification, reason))
        original = {name: {pair(r): r for r in sources[name][1]} for name in (NONRED, STRICT)}
        # Confirm each original's separately inspected semantics before changing anything.
        for name, flag in ((NONRED, 'control_eligible'), (STRICT, 'strict_control_eligible')):
            eligible_groups = {key: [r for _, r, _, _ in group if r[flag] == '1']
                               for key, group in alternatives.items() if key[1]}
            reproduced = {key: min(group, key=rank) for key, group in eligible_groups.items() if group}
            assert set(reproduced) == set(original[name]), (org, name, 'source semantics changed')
            assert all(row_key(reproduced[key]) == row_key(row) for key, row in original[name].items())
        chosen, selection_groups = {}, {}
        for name, flag in ((NONRED, 'control_eligible'), (STRICT, 'strict_control_eligible')):
            groups = {key: [r for _, r, site_class, _ in group if site_class == 'GOLD' and r[flag] == '1'
                           and (name != STRICT or not strict_biological_veto(r))]
                      for key, group in alternatives.items() if decisions[key]['classification'] == 'GOLD' and key[1]}
            groups = {key: group for key, group in groups.items() if group}
            selection_groups[name] = groups
            chosen[name] = {key: dict(min(group, key=rank)) for key, group in groups.items()}
            assert set(chosen[name]).issubset(original[name])
        outputs, changes = {}, collections.defaultdict(list)
        for name in (NONRED, STRICT):
            output_rows = []
            # Retain original row ordering; each input's selection is independent.
            for old in sources[name][1]:
                key = pair(old)
                if key not in chosen[name]: continue
                row, group = chosen[name][key], selection_groups[name][key]
                if name == NONRED:
                    row.update(candidate_structure_count=str(len({r['PDB_ID'] for r in group})),
                               candidate_row_count=str(len(group)),
                               representative_selection='lowest_available_resolution_then_most_contacts_then_identifiers')
                else:
                    previous = chosen[NONRED].get(key)
                    row.update(strict_candidate_structure_count=str(len({r['PDB_ID'] for r in group})),
                               strict_candidate_row_count=str(len(group)), control_set='strict',
                               previous_representative_PDB_ID=previous['PDB_ID'] if previous else '',
                               representative_changed_after_strict_filter=str(int(previous is not None and row_key(previous) != row_key(row))))
                if row_key(old) != row_key(row):
                    changes[key].append(f'{name}: {identifier(old)} -> {identifier(row)}; existing biological-site alternative selected')
                output_rows.append(rebase_paths(row, source_dir, output_dir))
            outputs[name] = output_rows
            write_table(output_dir / name, sources[name][0], output_rows)
            path = output_dir / name
            manifest['output_files'][f'{org}/{name}'] = dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                rows=len(output_rows), columns=len(sources[name][0]), input_file=f'{directory}/{name}')
        audits = []
        for record in org_records:
            key = record['AA'], record['accession']
            aa, accession = key
            decision, u = decisions[key], uniprot.get(accession, {})
            comments = functional_comments(u)
            func_pmids = pmids_in(comments)
            specific_pmids = list(decision['specific_pmids'])
            if accession == 'Q9P2J5' and aa == 'LEU': specific_pmids.append('33910001')
            if accession == 'Q9LYU8' and aa == 'LYS': specific_pmids.append('16731588')
            if accession == 'Q16878' and aa == 'CYS': specific_pmids.extend(['17135237', '29942080', '30946568'])
            if accession == 'P21888' and aa == 'CYS': specific_pmids.append('12032090')
            pdb_publications, missing_abstracts, pub_seen = [], [], set()
            for pdb in record['pdb_ids']:
                entry = entries.get(pdb, {})
                citation = entry.get('rcsb_primary_citation') or {}
                pmid = str(citation.get('pdbx_database_id_PubMed') or '')
                pub = publications.get(pmid, {})
                abstract = pub.get('abstractText', '')
                aa_terms = {AA_NAMES[aa]}
                if aa == 'ASP': aa_terms.add('aspartic acid')
                if aa == 'GLU': aa_terms.add('glutamic acid')
                abstract_aa_mentioned = any(re.search(r'\b'+re.escape(term)+r'\b', abstract, re.I) for term in aa_terms)
                abstract_context_terms = [term for term in ('free amino acid', 'peptide', 'substrate', 'product', 'activation', 'inhibition', 'transport', 'binding', 'cryoprotectant', 'mutant', 'engineered', 'aminoacyl') if term in abstract.lower()]
                if pmid: pub_seen.add(pmid)
                if pmid and not pub.get('abstractText'): missing_abstracts.append(pmid)
                pdb_publications.append(dict(PDB_ID=pdb, structure_title=entry.get('struct', {}).get('title', ''),
                    PMID=pmid, DOI=citation.get('pdbx_database_id_DOI', ''),
                    title=citation.get('title', ''), abstract_available=bool(pub.get('abstractText')),
                    canonical_AA_named_in_abstract=abstract_aa_mentioned, abstract_context_terms=abstract_context_terms,
                    publication_status=pub.get('pubType', '') or ('unlinked/unpublished' if not pmid else 'metadata available')))
            ec_numbers = sorted({c['reaction']['ecNumber'] for c in comments if c.get('reaction', {}).get('ecNumber')})
            support_by_identity = {}
            for _, row, _, _ in alternatives[key]:
                for fields in biolip[row['PDB_ID'], aa]:
                    receptor_match = fields[1] in (row['protein_chain'], row['protein_auth_chain'])
                    ligand_match = fields[5] in (row['ligand_label_asym_id'], row['ligand_auth_asym_id']) and fields[19].strip() == row['ligand_auth_seq_id'].strip()
                    identity = (fields[0], fields[1], fields[3], fields[4], fields[5], fields[19])
                    if identity in support_by_identity:
                        item = support_by_identity[identity]
                        item['receptor_chain_matches_original'] |= receptor_match
                        item['ligand_auth_chain_and_residue_match'] |= ligand_match
                        if receptor_match and ligand_match:
                            item['matched_original_rows'].append(identifier(row))
                        continue
                    support_by_identity[identity] = dict(PDB_ID=fields[0].upper(), receptor_auth_chain=fields[1],
                        binding_site_id=fields[3], ligand_auth_chain=fields[5], ligand_auth_seq_id=fields[19].strip(),
                        binding_residues=fields[7], annotated_UniProt=fields[17], PMID=fields[18],
                        receptor_chain_matches_original=receptor_match,
                        ligand_auth_chain_and_residue_match=ligand_match,
                        matched_original_rows=[identifier(row)] if receptor_match and ligand_match else [])
            support = list(support_by_identity.values())
            original_nonred, original_strict = original[NONRED].get(key), original[STRICT].get(key)
            corrected_nonred, corrected_strict = chosen[NONRED].get(key), chosen[STRICT].get(key)
            if decision['classification'] == 'GOLD' and original_strict is not None and corrected_strict is None:
                changes[key].append('strict representative removed: no biologically established ligand-site alternative passes the original strict structural flags')
                vetoes = sorted({strict_biological_veto(row) for _, row, _, _ in alternatives[key] if strict_biological_veto(row)})
                changes[key].extend(vetoes)
            urls = []
            if accession: urls.append(f'https://www.uniprot.org/uniprotkb/{accession}/entry')
            urls += [f'https://www.rcsb.org/structure/{pdb}' for pdb in record['pdb_ids']]
            urls += [f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/' for pmid in sorted(set(specific_pmids + func_pmids + list(pub_seen)))]
            urls += [f'https://www.kegg.jp/entry/ec:{ec}' for ec in ec_numbers if ec in kegg]
            urls.append('https://seq2fun.dcmb.med.umich.edu/BioLiP/')
            auxiliary = 'KEGG reaction records checked where EC annotation is present; BioLiP2 is corroboration, never sufficient for GOLD.'
            if accession == 'P0A7B5':
                auxiliary += ' BRENDA EC 2.7.2.11 and MetaCyc GLUTKIN-RXN additionally consulted for the glutamate-kinase substrate/feedback context.'
                urls += ['https://www.brenda-enzymes.org/enzyme.php?ecno=2.7.2.11', 'https://biocyc.org/reaction?orgid=META&id=GLUTKIN-RXN']
            summary = '; '.join([decision['reason']] +
                [f'{c["commentType"]}: ' + ('; '.join(t['value'] for t in c.get('texts', [])) or c.get('reaction', {}).get('name', ''))
                 for c in comments if c['commentType'] in ('FUNCTION', 'CATALYTIC ACTIVITY', 'ACTIVITY REGULATION')])
            alt_json = []
            for number, row, site_class, site_reason in alternatives[key]:
                alt_json.append(dict(source_file=f'../../{directory}/{ALL}', source_row_number=number,
                    PDB_ID=row['PDB_ID'], ligand_instance=row['ligand_instance'], protein_chain=row['protein_chain'],
                    resolution_A=row['resolution_A'], contact_residue_count=row['contact_residue_count'],
                    control_eligible=row['control_eligible'], strict_control_eligible=row['strict_control_eligible'],
                    ligand_contacting_protein_chain_count=row['ligand_contacting_protein_chain_count'],
                    protein_mutations=row['protein_mutations'], protein_mutations_current=row['protein_mutations_current'],
                    ligand_has_metal_coordination_current=row['ligand_has_metal_coordination_current'],
                    site_classification=site_class, site_reason=site_reason,
                    biological_strict_veto=strict_biological_veto(row),
                    retained_nonredundant_alternative=site_class == 'GOLD' and row['control_eligible'] == '1',
                    retained_strict_alternative=site_class == 'GOLD' and row['strict_control_eligible'] == '1' and not strict_biological_veto(row)))
            limits = ['GOLD is conservative for the free canonical AA, not a claim that every withheld pair is nonbinding.',
                      'BioLiP2 mirrored annotation snapshot Last-Modified 2026-03-29; lack of an annotation is not negative evidence.',
                      'Publication abstracts and annotated evidence were reviewed; full text was used for selected adjudications, not claimed for every publication.']
            if missing_abstracts: limits.append('Some linked publications have metadata but no retrieved abstract: '+','.join(sorted(set(missing_abstracts))))
            if not support: limits.append('No BioLiP2 canonical-AA annotation matched this candidate PDB/CCD in the downloaded snapshot.')
            audits.append(dict(zip(AUDIT_COLUMNS, [org, aa, accession, record['name'], decision['classification'],
                decision['functional_role'], decision['reason'], decision['decision_basis'],
                u.get('organism', {}).get('scientificName', ''), str(int(target_match(org, u.get('organism', {}).get('scientificName', '')))) if accession else '',
                summary, ';'.join(sorted(set(specific_pmids))), ';'.join(func_pmids), ';'.join(dict.fromkeys(urls)),
                compact(comments), compact(aa_features(u, aa)), ';'.join(ec_numbers), ';'.join(ec for ec in ec_numbers if ec in kegg), auxiliary,
                ';'.join(record['pdb_ids']), compact(pdb_publications), str(len(pub_seen)),
                str(sum(bool(publications.get(pmid, {}).get('abstractText')) for pmid in pub_seen)), ';'.join(sorted(set(missing_abstracts))),
                ';'.join(sorted({s['PDB_ID'] for s in support})), str(len(support)), compact(support),
                ';'.join(f'../../{directory}/{name}' for name in record['files']), str(len(alternatives[key])), compact(alt_json),
                identifier(original_nonred), identifier(corrected_nonred), identifier(original_strict), identifier(corrected_strict),
                str(int(corrected_nonred is not None)), str(int(corrected_strict is not None)), '; '.join(changes[key]),
                'UniProt functional/reaction/site annotations; all RCSB primary citations; available primary-paper abstracts; BioLiP2 annotations; EC/KEGG context; pair-specific primary-study adjudication where assigned',
                ' '.join(limits), manifest['curation_date_UTC'][:10]])))
        write_table(output_dir / 'biological_curation.tsv', AUDIT_COLUMNS, audits)
        counts = collections.Counter(r['classification'] for r in audits)
        manifest['organisms'][org] = dict(evaluated_pairs=len(audits), source_all_rows=len(all_rows),
            classification_counts=dict(counts), original_nonredundant_rows=len(sources[NONRED][1]),
            corrected_nonredundant_rows=len(outputs[NONRED]), original_strict_rows=len(sources[STRICT][1]),
            corrected_strict_rows=len(outputs[STRICT]), representative_changes=sum(bool(c) for c in changes.values()))
    manifest['candidate_pairs_evaluated'] = len(records)
    manifest['source_rows_preserved_in_originals_and_audit'] = sum(v['source_all_rows'] for v in manifest['organisms'].values())
    manifest['evidence_cache_counts'] = dict(UniProt_entries=len(uniprot), RCSB_entries=len(entries), MED_publications=len(publications), KEGG_EC_records=len(kegg))
    (OUT / 'curation_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    print(json.dumps(manifest['organisms'], indent=2))


if __name__ == '__main__':
    main()
