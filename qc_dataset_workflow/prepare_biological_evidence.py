"""Build candidate inventory and download fresh functional/publication evidence."""
from __future__ import annotations
import concurrent.futures
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import time
import urllib.parse
import requests
from biological_decisions import DECISIONS

BASE = Path(os.environ['AA_QC_DATA_ROOT']).resolve()
ROOT = BASE / 'curated_controls'
CACHE = ROOT / 'evidence_cache'
FOLDERS = {'ecoli': 'E_coli', 'arabidopsis': 'Arabidopsis', 'human': 'Human', 'mouse': 'Mouse', 'yeast': 'Yeast'}
ORGS = tuple(os.environ.get('AA_QC_ORGANISMS', ','.join(FOLDERS)).split(','))
ERRORS, DOWNLOADS = [], []


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf8')
    temporary.replace(path)


def download(url, path, payload=None, text=False):
    if path.is_file():
        return path.read_text(encoding='utf8') if text else json.loads(path.read_text(encoding='utf8'))
    for attempt in range(4):
        try:
            response = (requests.get(url, timeout=(15, 60)) if payload is None else
                        requests.post(url, json=payload, timeout=(15, 60)))
            response.raise_for_status()
            value = response.text if text else response.json()
            if not text and isinstance(value, dict) and value.get('errors'):
                raise ValueError('Partial GraphQL response: ' + str(value['errors']))
            path.parent.mkdir(parents=True, exist_ok=True)
            if text:
                path.write_text(value, encoding='utf8')
            else:
                save(path, value)
            DOWNLOADS.append(dict(url=url, path=str(path.relative_to(ROOT)), request=payload,
                                  retrieved_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                  sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            return value
        except (requests.RequestException, ValueError):
            if attempt == 3:
                raise
            time.sleep(1 + attempt)


def table(path):
    with path.open(encoding='utf8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def make_inventory():
    inventory = []
    for org in ORGS:
        groups = {}
        for name in ('all_AA_PDB_complexes.tsv', 'nonredundant_AA_protein_controls.tsv', 'strict_WT_single_protein_AA_controls.tsv'):
            for row in table(BASE / FOLDERS[org] / name):
                key = row['AA'], row['UniProt_accession']
                record = groups.setdefault(key, dict(organism=org, AA=key[0], accession=key[1],
                                                     name=row['protein_name'], pdb_ids=set(), files=set()))
                record['pdb_ids'].add(row['PDB_ID'])
                record['files'].add(name)
        for _, record in sorted(groups.items()):
            record['pdb_ids'], record['files'] = sorted(record['pdb_ids']), sorted(record['files'])
            inventory.append(record)
    save(ROOT / 'candidate_inventory.json', inventory)
    return inventory


def batch_jobs(jobs, required=False):
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(download, *arguments): key for key, arguments in jobs}
        for number, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                future.result()
            except Exception as error:
                ERRORS.append(dict(source=futures[future], error=str(error)))
                if required:
                    raise
            if number % 20 == 0 or number == len(futures):
                print(f'Biological evidence: {number}/{len(futures)} requests', flush=True)


def find_pmids(value):
    found = set()
    if isinstance(value, dict):
        if value.get('source') == 'PubMed' and str(value.get('id', '')).isdigit():
            found.add(str(value['id']))
        for child in value.values():
            found.update(find_pmids(child))
    elif isinstance(value, list):
        for child in value:
            found.update(find_pmids(child))
    return found


def biolip_matches(inventory):
    destination = CACHE / 'biolip2_matches.tsv'
    if destination.exists():
        return
    url = 'https://seq2fun.dcmb.med.umich.edu/BioLiP/download/BioLiP.txt.gz'
    if os.environ.get('AA_QC_SKIP_BIOLIP') == '1':
        destination.write_text('', encoding='utf8')
        ERRORS.append(dict(source=url, error='BioLiP explicitly skipped; no support inferred.'))
        return
    wanted = {(pdb, r['AA']) for r in inventory for pdb in r['pdb_ids']}
    temporary = destination.with_suffix('.tsv.tmp')
    scanned = matched = 0
    try:
        print('Streaming official BioLiP annotation download...', flush=True)
        with requests.get(url, stream=True, timeout=(15, 60)) as response:
            response.raise_for_status()
            with gzip.GzipFile(fileobj=response.raw) as archive, temporary.open('w', encoding='utf8', newline='') as out:
                for raw in archive:
                    fields = raw.decode('utf8').rstrip('\r\n').split('\t')
                    if len(fields) != 21:
                        raise ValueError(f'Unexpected BioLiP schema: {len(fields)} columns')
                    scanned += 1
                    if (fields[0].upper(), fields[4]) in wanted:
                        out.write('\t'.join(fields) + '\n')
                        matched += 1
                    if scanned % 100000 == 0:
                        print(f'BioLiP: {scanned} annotations scanned; {matched} matching', flush=True)
            DOWNLOADS.append(dict(url=url, last_modified=response.headers.get('Last-Modified', ''),
                                  annotations_scanned=scanned, matching_annotations=matched))
        temporary.replace(destination)
    except Exception as error:
        if temporary.exists():
            temporary.unlink()
        destination.write_text('', encoding='utf8')
        ERRORS.append(dict(source=url, error=str(error)))
        print('BioLiP support unavailable; recorded in evidence-fetch audit.', flush=True)


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    previous_manifest = CACHE / 'evidence_fetch_manifest.json'
    if previous_manifest.exists():
        previous = json.loads(previous_manifest.read_text(encoding='utf8'))
        DOWNLOADS.extend(previous.get('downloads', []))
        ERRORS.extend(previous.get('errors', []))
    for name in ('uniprot', 'rcsb', 'publications_v2', 'kegg'):
        (CACHE / name).mkdir(exist_ok=True)
    inventory = make_inventory()
    accessions = sorted({r['accession'] for r in inventory if r['accession']})
    pdbs = sorted({p for r in inventory for p in r['pdb_ids']})
    batch_jobs([('uniprot:' + a, ('https://rest.uniprot.org/uniprotkb/' + a + '.json', CACHE / 'uniprot' / (a + '.json')))
                for a in accessions])
    jobs = []
    for start in range(0, len(pdbs), 25):
        query = 'query Citations($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} rcsb_primary_citation{title year pdbx_database_id_PubMed pdbx_database_id_DOI journal_abbrev}}}'
        payload = dict(query=query, variables=dict(ids=pdbs[start:start + 25]))
        digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]
        jobs.append(('rcsb:' + str(start), ('https://data.rcsb.org/graphql', CACHE / 'rcsb' / f'{start:05d}_{digest}.json', payload)))
    batch_jobs(jobs, required=True)
    pmids, ecs = set(), set()
    for path in (CACHE / 'rcsb').glob('*.json'):
        for entry in json.loads(path.read_text(encoding='utf8'))['data']['entries']:
            if entry is None:
                raise ValueError('Missing requested RCSB citation entry')
            pmid = str((entry.get('rcsb_primary_citation') or {}).get('pdbx_database_id_PubMed') or '')
            if pmid.isdigit() and int(pmid) > 0:
                pmids.add(pmid)
    for path in (CACHE / 'uniprot').glob('*.json'):
        uniprot = json.loads(path.read_text(encoding='utf8'))
        comments = [c for c in uniprot.get('comments', []) if c.get('commentType') in
                    ('FUNCTION', 'CATALYTIC ACTIVITY', 'ACTIVITY REGULATION', 'BIOPHYSICOCHEMICAL PROPERTIES')]
        pmids.update(find_pmids(comments))
        for comment in comments:
            ec = comment.get('reaction', {}).get('ecNumber')
            if ec and '-' not in ec:
                ecs.add(ec)
    for record in inventory:
        pmids.update(DECISIONS.get((record['accession'], record['AA']), {}).get('specific_pmids', []))
    jobs, pmids, ecs = [], sorted(pmids), sorted(ecs)
    for start in range(0, len(pmids), 25):
        query = '(' + ' OR '.join('EXT_ID:' + p for p in pmids[start:start + 25]) + ') AND SRC:MED'
        url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode(
            dict(query=query, format='json', resultType='core', pageSize=1000))
        digest = hashlib.sha256(url.encode()).hexdigest()[:16]
        jobs.append(('publications:' + str(start), (url, CACHE / 'publications_v2' / f'{start:05d}_{digest}.json')))
    batch_jobs(jobs)
    jobs = []
    for start in range(0, len(ecs), 10):
        url = 'https://rest.kegg.jp/get/' + '+'.join('ec:' + ec for ec in ecs[start:start + 10])
        digest = hashlib.sha256(url.encode()).hexdigest()[:16]
        jobs.append(('kegg:' + str(start), (url, CACHE / 'kegg' / f'{start:05d}_{digest}.txt', None, True)))
    batch_jobs(jobs)
    biolip_matches(inventory)
    save(CACHE / 'evidence_fetch_manifest.json', dict(downloads=DOWNLOADS, errors=ERRORS, candidate_groups=len(inventory),
         UniProt_accessions=len(accessions), PDB_entries=len(pdbs),
         literature_scope='Public metadata, available abstracts and included reviewed decisions; new pairs require review.'))
    print(f'Evidence prepared: {len(inventory)} candidate groups; {len(ERRORS)} retrieval issues logged.', flush=True)


if __name__ == '__main__':
    main()
