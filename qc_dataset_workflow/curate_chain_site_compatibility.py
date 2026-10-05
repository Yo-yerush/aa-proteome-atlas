# Adapted from the existing reviewed workflow for fresh datasets.
# Original source SHA256: 709ef83665babded78f9744ce9c5518f32502ceb857fdacc5cb6f56a510523ba
"""Create instance-specific, evidence-backed sidecars without modifying app TSVs.

Requires the focused coordinate audit in assess_chain_compatibility.py. This is
not the structural search pipeline, a docking calculation, or an app loader.
"""
import collections
import csv
import datetime
import hashlib
import json
import pathlib
import re

from assess_chain_compatibility import ROOT, CACHE, COORDS, ORGS, INPUTS, load_controls, read_categories

LABELS = ('monomer_compatible', 'complex_dependent', 'uncertain')
SITE_FIELDS = ['organism', 'AA', 'UniProt_accession', 'protein_name', 'PDB_ID', 'ligand_instance', 'ligand_label_asym_id', 'ligand_auth_seq_id', 'protein_chain', 'protein_auth_chain', 'source_files', 'source_row_numbers_json', 'biological_classification', 'chain_compatibility', 'dependency_mechanism', 'confidence', 'frozen_extracted_chain_geometry', 'complex_dependent_conformation', 'reason', 'representative_scope', 'atlas_representation', 'audited_assembly_ids', 'assembly_selection', 'assembly_context_json', 'other_polymer_copy_count_5A', 'other_polymer_residue_count_4A', 'other_polymer_minimum_distance_A', 'symmetry_partner_contributes', 'other_polymer_contacts_json', 'asymmetric_unit_other_polymer_contacts_json', 'selected_model', 'available_models_json', 'ligand_unique_occupied_heavy_atom_count', 'ligand_expected_heavy_atom_count', 'ligand_model_completeness', 'ligand_alt_ids_json', 'protein_mutations_current', 'protein_fragment', 'cofactor_annotations_json', 'UniProt_subunit_domain_annotations_json', 'primary_publication_PMID', 'primary_publication_title', 'primary_publication_abstract_available', 'biological_curation_reference_urls', 'compatibility_reference_urls', 'caveats', 'coordinate_file', 'geometry_evidence_file', 'curation_date_UTC']
PAIR_FIELDS = ['organism', 'AA', 'UniProt_accession', 'protein_name', 'biological_classification', 'pair_summary', 'mixed_representatives', 'observed_site_classifications', 'selected_site_count', 'strict_site_classifications', 'nonredundant_site_classifications', 'site_details_json']
KEYS = ('organism_cohort','AA','PDB_ID','ligand_instance','protein_chain','UniProt_accession')
EXPECTED = dict(ALA=6,ARG=12,ASN=9,ASP=9,CYS=7,GLN=10,GLU=10,GLY=5,HIS=11,
                ILE=9,LEU=9,LYS=10,MET=9,PHE=12,PRO=8,SER=7,THR=8,TRP=15,TYR=13,VAL=8)

# Reviewed exceptions address the observed pocket, rather than oligomeric state.
# Evidence from a different organism/isoform is explicitly identified as inference.
OVERRIDES = {
 'Q03557': dict(classification='complex_dependent', mechanism='partner_scaffolded_glutaminase_site', conformation='required', confidence='moderate',
    reason='The yeast GatFAB primary structure/function study identifies the GatF N-terminal domain as a trans-acting scaffold for the GatA glutaminase active site. Removing that partner lowers glutamine-dependent activity and removes the native site scaffold; an apparently intrachain ligand-contact list does not certify this extracted GatA chain as an independent native binding site.',
    urls=['https://pubmed.ncbi.nlm.nih.gov/24692665/','https://www.rcsb.org/structure/4N0I']),
 'P32178': dict(classification='complex_dependent', mechanism='intersubunit_allosteric_AA_pocket', conformation='required', confidence='high',
    reason='Primary yeast chorismate-mutase studies place the free tryptophan/tyrosine regulatory pocket at the dimer interface. Tyrosine contacts residues of both monomers; tryptophan separates the allosteric domain of one monomer from helix H8 of the other. Extracting one chain removes an experimentally defined part of the native AA pocket.',
    urls=['https://pubmed.ncbi.nlm.nih.gov/8622937/','https://pubmed.ncbi.nlm.nih.gov/7971967/']),
 'P00962': dict(classification='complex_dependent', mechanism='tRNA_induced_site_assembly', conformation='required', confidence='high',
   reason='E. coli GlnRS uses tRNA-mediated induced fit to assemble the selective glutamine pocket. Removing tRNA can preserve a frozen protein conformation but removes the partner responsible for active-site assembly.',
   urls=['https://doi.org/10.1016/S0969-2126(03)00074-1','https://pmc.ncbi.nlm.nih.gov/articles/PMC2516378/','https://pubmed.ncbi.nlm.nih.gov/15845537/']),
 'P04805': dict(classification='complex_dependent', mechanism='tRNA_dependent_AA_recognition', conformation='required_by_homolog_inference', confidence='moderate',
   reason='The local modeled glutamate contacts are intrachain, but selective GluRS recognition is tRNA-dependent. E. coli Glu-Q-RS comparison and homolog structures identify tRNA as part of functional glutamate-site assembly; the dependency is inferred for this E. coli pose.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/18602926/','https://pubmed.ncbi.nlm.nih.gov/17161369/']),
 'O94925': dict(classification='complex_dependent', mechanism='oligomerization_controls_substrate_binding', conformation='required', confidence='high',
   reason='The glutamine contact shell is intrachain. Primary binding measurements and interface mutants nevertheless link kidney-type glutaminase oligomerization to formation of a glutamine-binding conformation; a complete frozen pocket does not establish partner-independent binding.',
   urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC6996896/','https://pmc.ncbi.nlm.nih.gov/articles/PMC5076503/']),
 'Q9UI32': dict(classification='uncertain', mechanism='possible_oligomer_coupled_binding_conformation', conformation='possible_not_established', confidence='moderate',
   reason='GLS2 glutamine contacts are intrachain, but its activation loop couples oligomerization, substrate access and apparent glutamine affinity. The selected K253A filament structure does not establish isolated WT-chain recognition; GLS1 binding requirements cannot simply be transferred to GLS2.',
   urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC8130709/','https://pmc.ncbi.nlm.nih.gov/articles/PMC10770349/','https://www.rcsb.org/structure/8T0Z']),
 'P29477': dict(classification='complex_dependent', mechanism='dimer_supported_active_center_conformation', conformation='required', confidence='high',
   reason='Murine iNOS structural comparisons show interface refolding and recruitment of elements that complete the substrate-bearing active center. Intrachain arginine contacts alone omit this dimer-supported pocket organization.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/9516116/','https://pubmed.ncbi.nlm.nih.gov/9334294/']),
 'P35228': dict(classification='complex_dependent', mechanism='dimer_supported_active_center_conformation', conformation='required_by_homolog_inference', confidence='moderate',
   reason='Human iNOS has an intrachain arginine contact shell, but the conserved inducible-NOS pocket is organized by interface refolding in murine iNOS. Partner dependence for this human conformation is an isoform-matched homolog inference, not an assembly-size rule.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/9516116/','https://www.rcsb.org/structure/1NSI']),
 'P29474': dict(classification='uncertain', mechanism='NOS_interface_coupling_unresolved_for_binding', conformation='possible_not_established', confidence='moderate',
   reason='The modeled eNOS arginine pocket is locally intrachain. NOS interface organization merits a flag, but constitutive-eNOS work differs from inducible-NOS active-site assembly. Dimer-dependent catalysis does not establish that isolated eNOS cannot recognize arginine.',
   urls=['https://doi.org/10.1016/S0092-8674(00)81718-3','https://pubmed.ncbi.nlm.nih.gov/9516116/']),
 'P29475': dict(classification='uncertain', mechanism='NOS_interface_coupling_unresolved_for_binding', conformation='possible_not_established', confidence='moderate',
   reason='The modeled nNOS arginine contacts are intrachain. Conserved NOS interface architecture may stabilize the substrate-bearing conformation, but isoform-specific evidence that another chain is required for arginine recognition is unresolved; catalytic dimer dependence is insufficient.',
   urls=['https://www.nature.com/articles/nsb0399_233','https://pubmed.ncbi.nlm.nih.gov/9516116/']),
 'Q9SCL7': dict(classification='uncertain', mechanism='interfacial_regulatory_helix_coupling', conformation='possible_not_established', confidence='moderate',
   reason='Arginine contacts in 2RD5 are entirely on the selected NAGK chain. Its regulatory helix participates in oligomer contacts and PII alters cleft motions, but inhibition/coupling data do not prove that a partner is necessary for this arginine pocket. Isolated-chain binding-site stability remains unresolved.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/17913711/','https://journals.asm.org/doi/10.1128/jb.01831-07']),
 'Q9SKE2': dict(classification='uncertain', mechanism='partner_stabilized_substrate_gate', conformation='partner_stabilized_observed_conformation;necessity_unresolved', confidence='moderate',
   reason='The 5ECK Ile contacts are intrachain, but FIP1 stabilizes a distinct FIN219 C-terminal orientation that rebuilds the substrate pocket. FIN219 alone can also recruit Ile after JA binds; independent persistence of this particular complex-form pocket after chain extraction is not established.',
   caveat='Primary sequential-binding experiments require prior jasmonate binding for Ile recruitment. A protein-only extraction omitting JA does not reproduce the experimental binding conditions.',
   urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC5347581/','https://pubmed.ncbi.nlm.nih.gov/28223489/']),
 'Q695T7': dict(classification='uncertain', mechanism='ancillary_partner_effect_on_site_unresolved', conformation='possible_not_established', confidence='moderate',
   reason='The B0AT1 amino-acid contact pocket is intrachain in an ACE2-containing experimental complex. ACE2/collectrin supports trafficking and activation, but these observations do not establish whether the observed binding conformation remains independently stable after extraction.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/32132184/','https://doi.org/10.1038/s41421-023-00596-2']),
 'Q9SIE1': dict(classification='complex_dependent', mechanism='partner_completes_AA_site_cofactor_pocket', conformation='cofactor_position_supported_by_partner', confidence='high',
   reason='Free glutamate contacts are intrachain, but adjacent-subunit Tyr132 anchors the PMP/PLP phosphate in the same amino-acid site. Extraction removes a cofactor-pocket contributor even though that residue is more than 5 A from glutamate itself.',
   urls=['https://doi.org/10.1111/tpj.13856']),
 'Q7L266': dict(classification='uncertain', mechanism='mature_proteolytic_segments_share_deposited_chain_label', conformation='not_established', confidence='moderate',
   reason='4OSY is the fully cleaved human asparaginase: mature alpha and beta peptide segments share deposited label A. The glycine neighborhood includes both sequence regions. Retention depends on whether extraction preserves both disconnected segments; one label does not establish one chemical polypeptide.',
   urls=['https://pubmed.ncbi.nlm.nih.gov/23601642/','https://www.rcsb.org/structure/4OSY','https://pmc.ncbi.nlm.nih.gov/articles/PMC4120204/']),
}

CONTEXT = {
 'Q05506': ('Yeast ArgRS binds free L-arginine without tRNA, as directly established by the 1BS2 primary study. tRNA can induce catalytic-center changes and is required for amino acid activation; the occupied-coordinate audit separately checks whether it contributes to this exact observed ligand pocket.', ['https://pubmed.ncbi.nlm.nih.gov/9736621/','https://pubmed.ncbi.nlm.nih.gov/11060012/']),
 'P00960': ('Glycine is buried in the catalytic alpha-chain pocket in 7EIV. The beta chain supports ATP/tRNA functions, without contributing modeled glycine contacts.', ['https://pmc.ncbi.nlm.nih.gov/articles/PMC8464048/']),
 'P27305': ('Glu-Q-RS is experimentally distinguished from GluRS by tRNA-independent glutamate recognition/activation; the GlnRS/GluRS dependency is not assigned to this paralog.', ['https://pubmed.ncbi.nlm.nih.gov/18602926/']),
 'P0A6F1': ('The glutamine pocket is in the small subunit; the large subunit supplies downstream phosphorylation/channeling functions. Whole-complex catalysis is not a binding-site requirement.', ['https://pubmed.ncbi.nlm.nih.gov/10950966/']),
 'P48775': ('The selected 6UD5 ligand occupies the intrachain E105/W208/R211 exosite, rather than the catalytic pocket whose roof includes a neighboring subunit. Classification is site-specific.', ['https://pmc.ncbi.nlm.nih.gov/articles/PMC8892992/']),
 'Q01650': ('The light-chain AA pocket is intrachain. Experiments show LAT1/LAT2 light chains remain functional without 4F2hc, while the heavy chain can modulate affinity/specificity.', ['https://doi.org/10.3390/ijms21207573']),
 'Q9UHI5': ('The light-chain AA pocket is intrachain. LAT2 can function without 4F2hc, although the heavy chain stabilizes transporters and changes affinity/specificity; equivalent physiological affinity is not assumed.', ['https://doi.org/10.3390/ijms21207573']),
 'Q8WTX7': ('CASTOR1 binds arginine between ACT domains of the same polypeptide; its oligomer/signaling interactions do not by themselves make the local site composite.', ['https://pubmed.ncbi.nlm.nih.gov/27487210/']),
 'P58004': ('The experimentally characterized Sestrin2 leucine sensor has a self-contained local pocket; downstream GATOR signaling is distinct from ligand-site construction.', ['https://www.rcsb.org/structure/5DJ4']),
 'P14618': ('The selected PKM2 amino-acid pocket is within one protomer. Tetramer-dependent enzyme regulation is distinguished from a requirement for another chain to supply the modeled AA pocket.', []),
 'P41180': ('The modeled tryptophan pocket is within one VFT chain. Dimeric CaSR activation and calcium-dependent ligand modulation are distinct caveats; they are not proof of an interchain tryptophan site.', []),
 'P60061': ('AdiC recognition is within an individual transporter protomer; native dimerization is not used as a rejection rule.', []),
 'P60063': ('AdiC recognition is within an individual transporter protomer; native dimerization is not used as a rejection rule.', []),
}

COMPOSITE_CONTEXT = {
 'P0A881':'The TrpR corepressor pocket is shared by the intertwined repressor dimer.',
 'P0A9T0':'Serine binds at the PGDH regulatory ACT-domain interface.',
 'P0ACI6':'AsnC effector binding occurs at an intersubunit regulatory interface.',
 'P42738':'The chorismate-mutase regulatory amino-acid pocket is shared between subunits.',
 'P0A786':'The native ATCase catalytic pocket includes a neighboring catalytic-chain residue, rather than depending merely on the regulatory-chain count.',
 'P31153':'Most methionine contacts for the selected chain-A pose come from a symmetry-generated chain-A copy; extracting just the deposited selected chain removes them.',
 'P0A817':'The methionine-adenosyltransferase pocket contains residues from both dimer partners.',
 'P17562':'The methionine-adenosyltransferase pocket contains residues from both dimer partners.',
 'P00439':'This is the phenylalanine-bound regulatory ACT-domain interface, not the intrachain catalytic phenylalanine pocket.',
 'P30047':'The GFRP phenylalanine pocket uses multiple chains of the GCH1/GFRP regulatory complex.',
 'O75311':'The GlyR orthosteric site spans principal and complementary extracellular subunit faces.',
 'P23415':'The GlyR orthosteric site spans principal and complementary extracellular subunit faces.',
 'P23416':'The GlyR orthosteric site spans principal and complementary extracellular subunit faces.',
 'P48167':'The GlyR orthosteric site spans principal and complementary extracellular subunit faces.',
}

def key(row):
    return tuple(row[k] for k in KEYS)

def json_cell(value):
    return json.dumps(value,ensure_ascii=False,separators=(',',':'))

def write_tsv(path, rows, columns=None):
    columns=columns or (list(rows[0]) if rows else (PAIR_FIELDS if path.name.endswith('_pairs.tsv') else SITE_FIELDS))
    with path.open('w',encoding='utf8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)

def provenance():
    entries={e['rcsb_id']:e for f in (ROOT/'evidence_cache/rcsb').glob('*.json')
             for e in (json.loads(f.read_text(encoding='utf8')).get('data',{}).get('entries') or [])}
    publications={e['id']:e for f in (ROOT/'evidence_cache/publications_v2').glob('*.json')
                  for e in json.loads(f.read_text(encoding='utf8'))['resultList']['result']}
    gold={}
    for org in ORGS:
        with (ROOT/org/'biological_curation.tsv').open(encoding='utf8',newline='') as handle:
            for row in csv.DictReader(handle,delimiter='\t'):
                if row['classification']=='GOLD': gold[(org,row['AA'],row['UniProt_accession'])]=row
    return entries,publications,gold

def preferred_assemblies(site):
    assemblies=site['biological_assembly_contacts']
    author=[a for a in assemblies if a['details'] in ('author_defined_assembly','author_and_software_defined_assembly')]
    return author or assemblies

def classify(row,site):
    acc=row['UniProt_accession']
    answer=dict(classification='uncertain',mechanism='unresolved',conformation='not_established',confidence='low',reason='',urls=[])
    if site is None and row['PDB_ID'] != '8WFI':
        answer.update(mechanism='coordinate_audit_unresolved',
          reason='No successful occupied-coordinate audit is available for this exact ligand instance and protein chain.',
          urls=['https://www.rcsb.org/structure/'+row['PDB_ID']])
        return answer
    if site is None:
        answer.update(mechanism='zero_occupancy_ligand',confidence='high',
          reason='All five deposited glycine heavy atoms in 8WFI have occupancy 0.00 and are excluded from the occupied-coordinate contact audit. The publication establishes biological glycine recognition, but this selected pose requires coordinate/experimental reconciliation before certifying compatibility.',
          urls=['https://www.rcsb.org/structure/8WFI','https://pubmed.ncbi.nlm.nih.gov/38513663/'])
        return answer
    assemblies=preferred_assemblies(site)
    others=[c for a in assemblies for c in a['polymer_contacts'] if not c['selected_chain']]
    if row['ligand_instance']=='5UAU.G':
        answer.update(mechanism='functional_site_identity_unresolved',confidence='moderate',
          reason='This proline contacts residues 136-164 in a local intrachain site, whereas the same entry contains a separate native product pocket around Thr238 at a dimer interface. Pair-level GOLD evidence does not establish that this selected secondary pose is the functional control site.',
          urls=['https://pubmed.ncbi.nlm.nih.gov/28258219/','https://www.rcsb.org/structure/5UAU'])
        return answer
    if acc in OVERRIDES:
        return dict(OVERRIDES[acc])
    if others:
        close=[c for c in others if c['residue_count_4A']]
        if close:
            closest=min(c['minimum_distance_A'] for c in close)
            answer.update(classification='complex_dependent',mechanism='native_interpolymer_AA_pocket',confidence='high' if len(close)>1 or sum(c['residue_count_4A'] for c in close)>1 else 'moderate',
                reason=f'The observed occupied ligand has a neighboring polymer copy within {closest:.2f} A in the relevant deposited biological assembly. Extraction removes part of its local contact environment. '+COMPOSITE_CONTEXT.get(acc,'The decision is a structural inference from this native pocket, not from the number of assembly subunits.'))
        else:
            answer.update(mechanism='peripheral_partner_contact_unresolved',confidence='moderate',
              reason='Another biological-assembly polymer enters the 4-5 A shell but supplies no contact at or below 4 A. Whether it forms an essential pocket wall or merely peripheral packing is unresolved; the contact alone is not treated as proof of complex dependence.')
        return answer
    if not assemblies or not any(c['selected_chain'] for a in assemblies for c in a['polymer_contacts']):
        answer.update(mechanism='assembly_or_selected_pocket_unresolved',reason='No applicable assembly containing the selected chain and ligand with a resolved selected-chain pocket was found.')
        return answer
    if row['PDB_ID']=='9NB3':
        answer.update(mechanism='unresolved_ligand_side_chain_atom',confidence='moderate',reason='The selected GLN pose lacks its side-chain amide nitrogen NE2. The remaining modeled contacts are intrachain, but an atom central to glutamine recognition is unresolved, limiting a claim that the full experimental free-AA environment is self-contained.')
        return answer
    note, urls=CONTEXT.get(acc,('',[]))
    name=row['protein_name'].lower()
    if not note:
        if 'glutamate receptor' in name or 'metabotropic glutamate' in name:
            note='The AA-recognition clamshell/VFT domain is within the selected polypeptide; oligomeric channel activation or downstream signaling does not by itself make this pocket interchain.'
        elif 'trna' in name or 'trna' in name.replace('-',''):
            note='The modeled cognate-AA pocket is within the selected catalytic chain. Partner requirements for tRNA binding/aminoacylation are distinct from supplying the observed free-AA pocket.'
        elif 'binding protein' in name or 'binding lipoprotein' in name or 'binding transport' in name or 'leu/ile/val-binding' in name:
            note='The solute-binding protein supplies its local AA-recognition pocket; its transport machinery is a separate downstream function.'
        elif 'ctp synthase' in name:
            note='The glutamine pocket is in one chain\'s glutaminase domain; nucleotide-site oligomerization and filament-dependent activity are distinct from its local AA site.'
            urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC8325340/']
        elif 'transporter' in name or 'antiporter' in name:
            note='The selected transporter chain supplies the local substrate-recognition pocket. Membrane/ion dependence and partner-dependent trafficking are separate caveats.'
        else:
            note='The validated free-AA protein pair and its selected experimental pocket have no modeled polymer contribution outside this chain; partner-independent affinity or whole-enzyme activity has not been measured by this audit.'
    answer.update(classification='monomer_compatible',mechanism='self_contained_observed_AA_pocket',confidence='moderate',reason='All occupied polymer contacts within 5 A belong to the selected chain in the relevant biological assembly. '+note,urls=urls)
    return answer

def main():
    input_hashes=json.loads((CACHE/'input_hashes.json').read_text(encoding='utf8'))
    protected={ROOT/p:h for p,h in input_hashes.items()}
    protected.update({ROOT/org/'biological_curation.tsv':hashlib.sha256((ROOT/org/'biological_curation.tsv').read_bytes()).hexdigest() for org in ORGS})
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==h for p,h in protected.items()),'A protected input changed since coordinate audit'
    controls=load_controls()
    geometry={key(s):s for s in json.loads((CACHE/'site_geometry.json').read_text(encoding='utf8'))}
    entries,publications,gold=provenance()
    membership=collections.defaultdict(lambda:collections.defaultdict(list))
    for org in ORGS:
        for name in INPUTS:
            with (ROOT/org/name).open(encoding='utf8',newline='') as handle:
                for ordinal,row in enumerate(csv.DictReader(handle,delimiter='\t'),2):
                    row['organism_cohort']=org
                    membership[key(row)][name].append(ordinal)
    output=[]
    for row in controls:
        site=geometry.get(key(row))
        decision=classify(row,site)
        biological=gold[(row['organism_cohort'],row['AA'],row['UniProt_accession'])]
        citation=(entries.get(row['PDB_ID'],{}).get('rcsb_primary_citation') or {})
        pmid=str(citation.get('pdbx_database_id_PubMed') or '')
        pub=publications.get(pmid,{})
        uniprot=json.loads((ROOT/'evidence_cache/uniprot'/(row['UniProt_accession']+'.json')).read_text(encoding='utf8'))
        subunits=[c for c in uniprot.get('comments',[]) if c.get('commentType') in ['SUBUNIT','DOMAIN']]
        cofactor=[c for c in uniprot.get('comments',[]) if c.get('commentType')=='COFACTOR']
        assemblies=preferred_assemblies(site) if site else []
        other_contacts=[]
        for a in assemblies:
            for c in a['polymer_contacts']:
                if c['selected_chain']:continue
                other_contacts.append(dict(assembly_id=a['assembly_id'],label_asym_id=c['label_asym_id'],auth_asym_id=c['auth_asym_id'],operation_ids=c['operation_ids'],
                    polymer_type=c['polymer_type'],minimum_distance_A=round(c['minimum_distance_A'],4),
                    residues=[dict(comp_id=r['comp_id'],auth_seq_id=r['auth_seq_id'],label_seq_id=r['label_seq_id'],minimum_distance_A=round(r['minimum_distance_A'],4)) for r in c['residues']],
                    partner_identity=site['original_API_context']['polymer_chains'].get(c['label_asym_id'],{})))
        asu_other=[c for c in (site or {}).get('asymmetric_unit_contacts',[]) if not c['selected_chain']]
        symmetry=any(c['label_asym_id']==row['protein_chain'] or c['operation_ids']!=a['ligand_anchor_operation']
                     for a in assemblies for c in a['polymer_contacts'] if not c['selected_chain'])
        known=set((site or {}).get('ligand_atom_names',[]))
        completeness='no_occupied_ligand_atoms' if site is None else ('complete' if len(known)==EXPECTED[row['AA']] else 'missing_modeled_heavy_atoms')
        caveats=[]
        if decision.get('caveat'):caveats.append(decision['caveat'])
        if completeness!='complete': caveats.append(f'Modeled unique heavy atoms {len(known)}/{EXPECTED[row["AA"]]}; only deposited occupied atoms audited.')
        if site and site['alt_ids']!=['.']: caveats.append('Alternate ligand conformers retained as a union for conservative contact screening; see alt_ids.')
        if cofactor:caveats.append('UniProt lists cofactors: extracting protein atoms alone can remove required small-molecule chemistry; monomer_compatible does not mean cofactor-free.')
        if 'transporter' in row['protein_name'].lower() or 'antiporter' in row['protein_name'].lower():caveats.append('Membrane/ions and whole transport function are outside an isolated pocket docking audit.')
        if row['protein_fragment']:caveats.append('Use the actual deposited construct/domain; compatibility does not transfer automatically to another construct or PDB conformation.')
        frozen='occupied_AA_contact_shell_complete_on_selected_chain' if site and not other_contacts else ('occupied_AA_contact_shell_requires_other_polymer' if site else 'unresolved')
        if row['UniProt_accession']=='Q7L266':frozen='deposited_label_contains_multiple_processed_peptide_segments'
        references=['https://files.rcsb.org/download/'+row['PDB_ID']+'.cif.gz','https://www.rcsb.org/structure/'+row['PDB_ID'],'https://www.uniprot.org/uniprotkb/'+row['UniProt_accession']]
        if pmid:references.append('https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/')
        doi=citation.get('pdbx_database_id_DOI')
        if doi:references.append('https://doi.org/'+doi)
        references=list(dict.fromkeys(references+decision['urls']))
        api_neighbors=(site or {}).get('original_API_context',{}).get('API_neighbors',[])
        output.append(dict(
          organism=row['organism_cohort'],AA=row['AA'],UniProt_accession=row['UniProt_accession'],protein_name=row['protein_name'],PDB_ID=row['PDB_ID'],
          ligand_instance=row['ligand_instance'],ligand_label_asym_id=row['ligand_label_asym_id'],ligand_auth_seq_id=row['ligand_auth_seq_id'],
          protein_chain=row['protein_chain'],protein_auth_chain=row['protein_auth_chain'],source_files=';'.join(row['source_files']),source_row_numbers_json=json_cell(dict(membership[key(row)])),
          biological_classification='GOLD',chain_compatibility=decision['classification'],dependency_mechanism=decision['mechanism'],confidence=decision['confidence'],
          frozen_extracted_chain_geometry=frozen,complex_dependent_conformation=decision['conformation'],reason=decision['reason'],
          representative_scope='exact_PDB_ligand_instance_chain; not_entire_UniProt',atlas_representation='extracted_experimental_PDB_chain',
          audited_assembly_ids=';'.join(a['assembly_id'] for a in assemblies),assembly_selection='author_defined_preferred_else_deposited_software; no_monomer_gate',
          assembly_context_json=json_cell([dict(assembly_id=a['assembly_id'],details=a['details'],oligomeric_details=a['oligomeric_details'],ligand_anchor_operation=a['ligand_anchor_operation']) for a in (site or {}).get('biological_assembly_contacts',[])]),
          other_polymer_copy_count_5A=max((sum(not c['selected_chain'] for c in a['polymer_contacts']) for a in assemblies),default=0) if site else '',
          other_polymer_residue_count_4A=max((sum(c['residue_count_4A'] for c in a['polymer_contacts'] if not c['selected_chain']) for a in assemblies),default=0) if site else '',
          other_polymer_minimum_distance_A=f'{min(c["minimum_distance_A"] for c in other_contacts):.4f}' if other_contacts else '',
          symmetry_partner_contributes=str(symmetry).lower(),other_polymer_contacts_json=json_cell(other_contacts),
          asymmetric_unit_other_polymer_contacts_json=json_cell([dict(label_asym_id=c['label_asym_id'],polymer_type=c['polymer_type'],minimum_distance_A=round(c['minimum_distance_A'],4)) for c in asu_other]),
          selected_model=(site or {}).get('selected_model',''),available_models_json=json_cell((site or {}).get('available_models',[])),
          ligand_unique_occupied_heavy_atom_count=len(known),ligand_expected_heavy_atom_count=EXPECTED[row['AA']],ligand_model_completeness=completeness,
          ligand_alt_ids_json=json_cell((site or {}).get('alt_ids',[])),
          protein_mutations_current=row['protein_mutations_current'],protein_fragment=row['protein_fragment'],
          cofactor_annotations_json=json_cell(cofactor),UniProt_subunit_domain_annotations_json=json_cell(subunits),
          primary_publication_PMID=pmid,primary_publication_title=citation.get('title',''),primary_publication_abstract_available=str(bool(pub.get('abstractText'))).lower(),
          biological_curation_reference_urls=biological['reference_urls'],compatibility_reference_urls=';'.join(references),
          caveats=' '.join(caveats),coordinate_file='../../evidence_cache/chain_compatibility/coordinates/'+row['PDB_ID']+'.cif.gz',
          geometry_evidence_file='../evidence_cache/chain_compatibility/site_geometry.json',curation_date_UTC=datetime.datetime.now(datetime.timezone.utc).date().isoformat()))
    # Coordinate paths are relative to each organism folder, not the workspace.
    for r in output:r['coordinate_file']='../evidence_cache/chain_compatibility/coordinates/'+r['PDB_ID']+'.cif.gz'
    for org in ORGS:
        write_tsv(ROOT/org/'binding_site_chain_compatibility.tsv',[r for r in output if r['organism']==org])
    write_tsv(ROOT/'binding_site_chain_compatibility_all.tsv',output)
    groups=collections.defaultdict(list)
    for r in output:groups[(r['organism'],r['AA'],r['UniProt_accession'])].append(r)
    pairs=[]
    for (org,aa,acc),items in sorted(groups.items()):
        labels=set(r['chain_compatibility'] for r in items)
        # "uncertain" is deliberately used for mixed representatives; retain all
        # classifications rather than treating one favorable pose as universal.
        pairs.append(dict(organism=org,AA=aa,UniProt_accession=acc,protein_name=items[0]['protein_name'],
            biological_classification='GOLD',pair_summary=next(iter(labels)) if len(labels)==1 else 'uncertain',
            mixed_representatives=str(len(labels)>1).lower(),
            observed_site_classifications=';'.join(sorted(labels)),selected_site_count=len(items),
            strict_site_classifications=';'.join(sorted({r['chain_compatibility'] for r in items if INPUTS[0] in r['source_files'].split(';')})),
            nonredundant_site_classifications=';'.join(sorted({r['chain_compatibility'] for r in items if INPUTS[1] in r['source_files'].split(';')})),
            site_details_json=json_cell([{k:r[k] for k in ['PDB_ID','ligand_instance','protein_chain','source_files','chain_compatibility','dependency_mechanism']} for r in items])))
    write_tsv(ROOT/'binding_site_chain_compatibility_pairs.tsv',pairs)
    manifest=dict(created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),atlas='extracted_experimental_PDB_chains',
        unique_sites=len(output),unique_GOLD_AA_UniProt_pairs=len(pairs),unique_PDB_entries=len({r['PDB_ID'] for r in output}),
        source_rows_covered=sum(len(v) for d in membership.values() for v in d.values()),
        site_classification_counts=dict(collections.Counter(r['chain_compatibility'] for r in output)),
        organisms={org:dict(unique_sites=sum(r['organism']==org for r in output),
            classification_counts=dict(collections.Counter(r['chain_compatibility'] for r in output if r['organism']==org)),
            source_file_counts={name:dict(collections.Counter(r['chain_compatibility'] for r in output if r['organism']==org and name in r['source_files'].split(';'))) for name in INPUTS}) for org in ORGS},
        protected_input_sha256={str(p.relative_to(ROOT)).replace('\\','/'):h for p,h in protected.items()},
        app_run=False,app_input_tables_modified=False,structural_search_pipeline_run=False)
    (ROOT/'chain_compatibility_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==h for p,h in protected.items())
    print(json.dumps({k:manifest[k] for k in ['unique_sites','unique_GOLD_AA_UniProt_pairs','source_rows_covered','site_classification_counts','organisms']},indent=2),flush=True)

if __name__=='__main__':main()
