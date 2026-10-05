"""Build free-amino-acid QC datasets from scratch using public APIs."""
from __future__ import annotations
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

ORGANISMS = ('ecoli', 'arabidopsis', 'human', 'mouse', 'yeast')
DEFAULT_ORGANISMS = ORGANISMS[:4]
AMINO_ACIDS = 'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parent.parent / 'qc_dataset_output', help='New output directory; defaults to qc_dataset_output beside the scripts folder.')
    parser.add_argument('--organisms', nargs='+', choices=ORGANISMS, default=list(DEFAULT_ORGANISMS))
    parser.add_argument('--aa', nargs='+', choices=AMINO_ACIDS, default=AMINO_ACIDS)
    parser.add_argument('--resume', action='store_true', help='Continue a run owned by this workflow.')
    parser.add_argument('--skip-biolip', action='store_true', help='Omit optional BioLiP corroboration and record the omission.')
    args = parser.parse_args()
    absent = [name for name in ('numpy', 'requests', 'scipy') if importlib.util.find_spec(name) is None]
    if absent:
        parser.error('Install dependencies: python -m pip install ' + ' '.join(absent))
    out = args.out.expanduser().resolve()
    package = Path(__file__).resolve().parent
    if out == package or package in out.parents:
        parser.error('Choose an output directory outside this scripts-only package.')
    # API caches can exceed MAX_PATH in deeply nested Windows project folders.
    if os.name == 'nt' and not str(out).startswith('\\\\?\\'):
        absolute = str(out)
        out = Path('\\\\?\\UNC\\' + absolute[2:] if absolute.startswith('\\\\') else '\\\\?\\' + absolute)
    configuration = dict(organisms=sorted(set(args.organisms)), amino_acids=sorted(set(args.aa)), skip_biolip=args.skip_biolip)
    marker = out / 'workflow_manifest.json'
    current_code = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in package.glob('*.py')}
    if args.resume:
        if not marker.is_file():
            parser.error('--resume requires a workflow_manifest.json from this workflow.')
        manifest = json.loads(marker.read_text(encoding='utf8'))
        if manifest.get('configuration') != configuration or manifest.get('workflow') != 'from_scratch_AA_QC':
            parser.error('Resume configuration differs from this run.')
        if manifest.get('script_sha256') != current_code:
            previous_code = manifest.get('script_sha256', {})
            changed = {name for name in set(previous_code) | set(current_code)
                       if previous_code.get(name) != current_code.get(name)}
            keep = {'search', 'structural'}
            if changed & {'run_workflow.py', 'build_aa_pdb_dataset.py'}:
                keep = set()
            elif changed & {'build_strict_aa_controls.py', 'split_aa_pdb_by_organism.py'}:
                keep = {'search'}
            manifest['completed_stages'] = [s for s in manifest['completed_stages'] if s in keep]
            manifest['script_sha256'] = current_code
    else:
        if out.exists() and any(out.iterdir()):
            parser.error('Output directory is not empty. Choose a new --out or --resume for this workflow.')
        out.mkdir(parents=True, exist_ok=True)
        manifest = dict(workflow='from_scratch_AA_QC', configuration=configuration, status='running',
                        started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(), completed_stages=[],
                        script_sha256=current_code, app_launched=False, app_loader_executed=False)
    env = dict(os.environ, AA_QC_DATA_ROOT=str(out), AA_QC_ORGANISMS=','.join(configuration['organisms']),
               PYTHONDONTWRITEBYTECODE='1', AA_QC_SKIP_BIOLIP='1' if args.skip_biolip else '0')
    env.pop('AA_QC_INCLUDE_ALTERNATIVES', None)
    manifest['status'] = 'running'

    def save_manifest():
        marker.write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf8')

    def stage(name, script, arguments=(), include_alternatives=False):
        if name in manifest['completed_stages']:
            print('Resuming completed stage: ' + name, flush=True)
            return
        print('Stage: ' + name, flush=True)
        stage_env = dict(env)
        if include_alternatives:
            stage_env['AA_QC_INCLUDE_ALTERNATIVES'] = '1'
        result = subprocess.run([sys.executable, '-B', str(package / script), *map(str, arguments)],
                                cwd=out, env=stage_env, check=False)
        if result.returncode:
            manifest.update(status='incomplete', failed_stage=name, exit_code=result.returncode)
            save_manifest()
            raise SystemExit(result.returncode)
        manifest['completed_stages'].append(name)
        save_manifest()

    save_manifest()
    search_args = ['--out', out, '--organisms', *configuration['organisms'], '--aa', *configuration['amino_acids']]
    if args.resume and 'search' not in manifest['completed_stages']:
        cached = sorted(p.parent for p in (out / 'raw').glob('*/manifest.json')
                        if json.loads(p.read_text(encoding='utf8')).get('configuration'))
        if cached:
            search_args += ['--resume', cached[-1]]
    stage('search', 'build_aa_pdb_dataset.py', search_args)
    stage('structural', 'build_strict_aa_controls.py', ['--out', out])
    stage('evidence', 'prepare_biological_evidence.py')
    stage('biological_curation', 'curate_biological_controls.py')
    stage('curated_validation', 'validate_curated_controls.py')
    stage('coordinate_audit', 'assess_chain_compatibility.py', include_alternatives=True)
    stage('chain_compatibility', 'curate_chain_site_compatibility.py')
    stage('chain_validation', 'validate_chain_site_compatibility.py')
    stage('AA_only_QC', 'build_AA_only_qc.py')
    manifest.pop('failed_stage', None)
    manifest.pop('exit_code', None)
    manifest.update(status='complete', completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    save_manifest()
    print('Complete. Validated controls: ' + str(out / 'curated_controls'), flush=True)
    print('Protein + AA-only QC subset: ' + str(out / 'qc_protein_AA_only_controls'), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
