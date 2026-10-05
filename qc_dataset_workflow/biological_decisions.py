# Adapted from the existing reviewed workflow for fresh datasets.
# Original source SHA256: 4f6d6797556ec8a5cfca28550ab196c26c2d9de7dc7d98f658ea02cbd53c15d8
"""Explicit, conservative AA-specific adjudications; no name/keyword GOLD rule.

Entries are (accession, canonical AA) and apply only to the correct source organism.
GOLD requires demonstrated free-AA substrate, transport, or regulatory recognition.
Products alone and proposed functions are held out, rather than promoted from PDB.
"""
DECISIONS = {}

def add(accessions, aas, classification, role, reason, pmids=''):
    for accession in accessions.split():
        for aa in aas.split():
            key=(accession,aa)
            if key in DECISIONS: raise ValueError('Duplicate decision '+str(key))
            DECISIONS[key]={'classification':classification,'functional_role':role,
                'reason':reason,'specific_pmids':pmids.split()}

# Saccharomyces cerevisiae review, 2026-10-04. These decisions cover every mapped
# pair in the fresh yeast candidate inventory. Unmapped identities remain manual.
add('Q05506','ARG','GOLD','cognate aaRS substrate',
    'The yeast ArgRS primary study directly establishes recognition of free L-arginine and demonstrates that tRNA is not required for its binding. tRNA-dependent amino acid activation is a separate requirement; the deposited 1BS2 binary complex is a free-arginine site.',
    '9736621 11060012')
add('P09436','ILE','GOLD','cognate aaRS substrate',
    'Yeast IleRS activates free L-isoleucine for tRNA aminoacylation. The 8WND study and the ScIleRS structures in the reveromycin study establish the cognate free-isoleucine catalytic pocket. Bound tRNA and reveromycin remain structural-composition caveats, not evidence that the observed isoleucine is a peptide residue.',
    '39738040 41839287')
add('Q12109','TRP','GOLD','cognate aaRS substrate',
    'The yeast TrpRS primary study characterizes tryptophan activation and free-tryptophan ligand complexes. Free L-tryptophan is the cognate substrate; sulfate mimics ATP phosphate in the ligand-bound catalytic context.',
    '20123733 9046085')
add('Q01217','ARG','GOLD','feedback inhibitor',
    'Purified yeast acetylglutamate kinase is feedback-inhibited by free L-arginine. The primary structure/function study directly examines arginine inhibition and the regulatory pocket, including the effect of deleting DUF619. The observed 3ZZH construct is the kinase domain and must retain its fragment/assembly caveats.',
    '22529931')
add('P27616','ASP','GOLD','free amino acid substrate',
    'Purified yeast SAICAR synthetase and its curated reaction identify free L-aspartate as the ATP-dependent ligation substrate. In 2CNU, ASP instance B contacts the residues that bind the aspartyl portion of the separate SAICAR product in 2CNV. The additional peripheral ASP sites are not automatically accepted as functional controls.',
    '2684279 1420588 9551557 26072057')
add('Q03557','GLN','GOLD','free amide-donor substrate',
    'Yeast mitochondrial GatFAB uses free glutamine as its ammonia donor. The GatFAB primary study measures glutaminase and transamidase activity with glutamine and identifies GatF as a trans-acting scaffold for the GatA glutaminase active site. This is free glutamine recognition, distinct from glutamyl-tRNA recognition; the site is complex-dependent.',
    '24692665 19417106')
add('P32477','GLU','GOLD','free amino acid substrate',
    'Yeast glutamate-cysteine ligase directly consumes free L-glutamate in glutathione biosynthesis. Its primary structure/function study resolves glutamate-bound enzyme and explains Mg-dependent substrate recognition. Native ligand metal coordination excludes these poses from the no-metal strict set.',
    '19726687 1687097')
add('P40510','SER','GOLD','allosteric feedback inhibitor',
    'The 2026 Ser33 primary study combines ligand-bound structures, enzyme assays, kinetics and biophysical analysis to establish free L-serine as a negative allosteric regulator. The regulatory AA is not the phosphoglycerate substrate or a peptide residue.',
    '42321184')
add('P32178','TRP','GOLD','allosteric activator',
    'Yeast chorismate mutase has directly characterized regulatory binding of free L-tryptophan, which activates the enzyme. Primary structures and functional studies place tryptophan at a native intersubunit regulatory site. Mutant constructs remain subject to the original WT filters.',
    '7971967 9384560 2187528 31498992')
add('P32178','TYR','GOLD','allosteric feedback inhibitor',
    'Free L-tyrosine is the characterized feedback inhibitor of yeast chorismate mutase. The T-state structure and biochemical studies establish regulatory recognition at the same intersubunit site used by activating tryptophan.',
    '8622937 9384560 2187528 31498992')
add('P25605','VAL','GOLD','allosteric feedback inhibitor',
    'Purified yeast acetolactate synthase studies show that the Ilv6 regulatory subunit confers inhibition by free valine. The fungal AHAS structure study establishes the regulatory communication pathway. 6WO1 is a cross-species catalytic/regulatory assembly; the yeast attribution applies to its Ilv6 subunit, not the foreign catalytic subunit.',
    '10213630 8972574 32640464')
add('P32449','PHE','SILVER','noncognate regulatory-site association',
    'This accession is tyrosine-regulated Aro4, not phenylalanine-regulated Aro3. Primary functional experiments show that G226S switches Aro4 to phenylalanine regulation. The WT phenylalanine-bound structure and mutant structure do not establish strong native WT free-Phe control specificity; do not transfer mutant or Aro3 function to this pair.',
    '12540830 15019786')
add('P02994','GLN','EXCLUDE','covalent posttranslational modification',
    'The 5O8W primary study identifies mono-L-glutamine covalently linked through its alpha-amino group to eEF1A Glu45. This glutaminylation is a protein modification, not a standalone noncovalently bound free-glutamine ligand, despite the search annotation.',
    '28801462')
add('P32481','MET','EXCLUDE','aminoacyl-tRNA-associated methionine',
    'The deposited yeast translation-initiation complexes contain eIF2-GTP-Met-tRNAi. Native eIF2 gamma recognizes the aminoacylated initiator tRNA; its methionine moiety is not independent biological recognition of free methionine. Missing or incomplete deposited linkage annotation does not change that ligand context.',
    '26212456')
add('N1P8Q8','GLN','EXCLUDE','N-terminal peptide substrate',
    'The Nta1 primary study characterizes deamidation of N-terminal Gln/Asn residues in N-degron peptides, with peptide-bound structures, peptide affinity and peptide kinetics. This establishes protein/peptide-residue recognition, not native free-glutamine ligand specificity.',
    '27791147')
add('P40051','GLY','EXCLUDE','peptide-cleavage context',
    'Icp55 is a mitochondrial aminopeptidase characterized with protein/peptide substrates and the apstatin inhibitor. Its primary study does not establish a physiological free-glycine substrate or regulatory function. A deposited glycine in this cleavage/inhibitor context cannot be promoted to a free-AA positive control.',
    '30582634 19720832')
add('P00950','ALA','STRUCTURAL_ONLY','phosphoglycerate-mutase context',
    'Gpm1 catalyzes phosphoglycerate interconversion; the 5PGM primary study interprets sulfate sites and the phosphoryl-transfer mechanism. It does not establish alanine as a native functional free-AA ligand.',
    '10064712')
add('P40356','ALA ASP THR','STRUCTURAL_ONLY','transcription-complex context',
    'Med3 is part of the yeast Mediator transcription pre-initiation complex. The 7UIO primary study addresses divergent-promoter transcription and Med-PIC assembly, without establishing functional free-AA recognition by Med3.',
    '36731470')
add('P40969','ARG ASN GLN MET PHE THR','STRUCTURAL_ONLY','kinetochore-complex context',
    'Cep3 recognizes centromeric DNA in the CBF3 kinetochore complex. The 6GYP primary study establishes DNA-sequence recognition and complex architecture, not a biological free-AA ligand function for this accession. Nonpolymer AA annotations alone are insufficient.',
    '30478265')
add('P38629','GLU THR','STRUCTURAL_ONLY','DNA-clamp-loader context',
    'Rfc3 participates in ATP-dependent DNA checkpoint-clamp loading. The 7ST9 primary study resolves Rad24-RFC, the 9-1-1 clamp and DNA; no functional free-AA recognition is established for Rfc3.',
    '35314831')
add('P21524','GLY','STRUCTURAL_ONLY','ribonucleotide-reductase context',
    'Rnr1 recognizes ribonucleotide substrates and nucleotide allosteric effectors. The 1ZYZ primary study addresses dNTP regulation and does not establish free glycine as a native functional ligand.',
    '16537479')
add('P30822','GLU','STRUCTURAL_ONLY','nuclear-export inhibitor context',
    'Crm1 recognizes nuclear-export cargo and is studied in 8HQ6 with aminoratjadone-derived inhibitors. That publication does not establish free glutamate as a physiological functional ligand of yeast Crm1.',
    '37595020')
add('P22515','PRO','STRUCTURAL_ONLY','ubiquitin-activation context',
    'Uba1 activates ubiquitin by recognizing its protein C terminus. The 3CMM primary study concerns the Uba1-ubiquitin complex and does not establish native free-proline binding.',
    '18662542')
add('P49686 Q12315','PRO','STRUCTURAL_ONLY','nuclear-pore protein-interaction context',
    'Nup42 and Gle1 form an mRNA-export regulatory protein complex. The 6B4E primary structural/functional study establishes their interaction and helicase regulation, without demonstrating a free-proline ligand function for either chain.',
    '29899397')
add('Q12125','PRO','STRUCTURAL_ONLY','tail-anchored protein-delivery context',
    'This accession is Get4. The 3LKU primary study characterizes the Get4/Get5 protein complex and its nucleotide-dependent interaction with Get3, not functional recognition of free proline.',
    '20554915')

# Cognate free amino acids are activated before transfer to tRNA. This is distinct
# from elongation factors or transferases that recognize aminoacyl-tRNA donors.
for accessions, aa in [
    ('P11875 P54136','ARG'),('P21888','CYS'),('P00962','GLN'),
    ('P04805','GLU'),('P00960 P41250','GLY'),('P12081','HIS'),
    ('P07813 Q9P2J5','LEU'),('P0A8N3 P0A8N5 Q15046','LYS'),
    ('P00959 P56192','MET'),('P08312 Q9Y285 O95363','PHE'),
    ('P07814','PRO'),('P49591','SER'),('P0A8M3','THR'),
    ('P00954 P67589 P23381','TRP'),('P0AGJ9 A0A0H2UKY9 P54577','TYR')]:
    add(accessions,aa,'GOLD','cognate aaRS substrate',
        'Cognate free amino acid is directly recognized and ATP-activated for aminoacylation; reaction and ligand-bound catalytic context support functional binding.')

# Established metabolic substrates, specified individually to avoid confusing
# free amino acids with peptide residues, phosphorylated analogues, or products.
for accession, aas in [
    ('P00805','ASN'),('P0A962','ASN'),('P0A6E4','ASP'),('P0A786','ASP'),
    ('P08660','ASP'),('P0A6F1','GLN'),('P0A6W9','CYS GLU'),
    ('P0A7B5','GLU'),('P0A7E5','GLN'),('P0A817','MET'),
    ('P22106','GLN'),('P22634','GLU'),('P27305','GLU'),
    ('P0C0U4','GLU'),('P09546','PRO'),('P77444','CYS'),
    ('O15067','GLN'),('O94925','GLN'),('P00439','PHE'),
    ('P00966','ASP'),('P17812','GLN'),('P29474','ARG'),
    ('P29475','ARG'),('P35228','ARG'),('P31153','MET'),
    ('P32929','CYS'),('P34897','GLY'),('P48651','SER'),
    ('P48775','TRP'),('P50440','ARG GLY'),('P14902','TRP'),
    ('Q6ZQW0','TRP'),('Q16878','CYS'),('Q9NRF8','GLN'),
    ('Q9UI32','GLN'),('P17562','MET'),('Q94C74','SER'),
    ('Q9LSQ4','ASP'),('Q9SIE1','GLU'),('Q8RY79','PHE'),
    ('P29477','ARG'),('P70303','GLN'),('P70698','GLN')]:
    add(accession,aas,'GOLD','free-AA metabolic substrate',
        'The canonical free amino acid is a recognized substrate in the documented catalytic reaction, corroborated by protein-specific functional annotations and experimental structural context.')

# Receptors: glutamate belongs to GluN2, glycine to GluN1, not both subunits.
add('O00222 O15399 P39086 P41594 P42262 Q12879 Q13224 Q14416 Q14831 Q14832 Q14833 P19492 Q9Z2W9 Q01097 P35436',
    'GLU','GOLD','physiological receptor agonist',
    'Free L-glutamate is the documented agonist of this receptor/subunit, with ligand-dependent signaling or channel activation evidence.')
add('O75311 P23415 P23416 P48167 Q05586 P35438','GLY','GOLD',
    'physiological receptor agonist/coagonist',
    'Free glycine is the documented ligand at this receptor/subunit; the role is distinct from glycine residues in proteins.')
add('Q9C8E7','CYS GLU GLY MET','GOLD','plant GLR agonist',
    'Protein-specific ligand-binding measurements and amino-acid-elicited calcium-response experiments support free-AA recognition by GLR3.3; canonical ligands occupy its ligand-binding domain.', '31871183 18162597')
add('Q93YT1','GLY MET','GOLD','plant GLR agonist',
    'GLR3.2 ligand-binding structures are supplemented by wild-type channel activation in the presence of CNIH4, supporting glycine agonism and methionine partial agonism.','33027636')
add('Q8GXJ4','GLU SER','GOLD','plant GLR agonist',
    'GLR3.4-specific amino-acid stimulation and ligand-binding-domain studies support free glutamate and serine recognition.','18162597 22447719 34161757')
add('Q8GXJ4','MET','SILVER','plant GLR ligand',
    'Methionine-bound ligand-binding domain is reported, but methionine-specific functional activation is less firmly established than for glutamate/serine in the evidence reviewed.','34161757')

# Transport substrates and periplasmic receptors: only verified specificities.
for accession, aas in [
    ('P02942','SER'),('P07017','ASP'),('P04816','LEU'),
    ('P0AD96','ILE LEU VAL'),('P0AEQ3','GLN'),('P0AEU0','HIS'),
    ('P37902','ASP'),('P28635','MET'),('P60061','ARG'),('P60063','ARG'),
    ('P43003','ASP'),('P43004','GLU'),('P43005','ASP CYS GLU'),
    ('P48067','GLY'),('Q9Y345','GLY'),('P82251','ARG'),
    ('Q01650','LEU PHE TRP TYR'),('Q9UHI5','LEU TRP'),
    ('Q92536','ARG LEU'),('Q8N370','PHE'),('Q9NS82','ALA'),
    ('Q9NP91','PRO'),('Q99884','PRO'),('Q9UPY5','GLU'),
    ('Q695T7','GLN LEU MET PHE TRP'),('Q15758','ALA GLN'),
    ('Q9UM01','ARG LEU LYS'),('Q9FMF7','GLU'),('Q09143','ARG LYS')]:
    add(accession,aas,'GOLD','free-AA transport/chemotaxis substrate',
        'The specific free canonical amino acid is an experimentally supported transported substrate or chemotactic ligand; ligand recognition is part of the native protein function.')
add('Q9UM01','GLN','SILVER','neutral-AA exchange substrate candidate',
    'Neutral-AA exchange is established, but glutamine-specific functional support for the deposited candidate is weaker in the reviewed sources; an unpublished PDB title is insufficient for GOLD.')

for accession, aas in [
    ('P00864','ASP'),('P00888','TYR'),('P04968','ILE'),
    ('P0A6D0','ARG'),('P0A6L2','LYS'),('P08660','LYS'),
    ('P0A7B5','PRO'),('P0A9D4','CYS'),('P0A9T0','SER'),
    ('P0AB91','PHE'),('P0ACI6','ASN'),('P0A881','TRP'),
    ('P0ADF8','VAL'),('P32322','PRO'),('P30047','PHE'),
    ('P58004','LEU'),('Q8WTX7','ARG'),('P41180','TRP'),
    ('Q9SCL7','ARG'),('Q9LYU8','LYS'),('Q9FVC8','LYS'),
    ('Q9LZX6','LYS'),('P42738','PHE TYR'),('Q93YZ7','VAL')]:
    add(accession,aas,'GOLD','free-AA regulation/sensing',
        'The specific free amino acid is a documented feedback inhibitor, activator, or sensing ligand with functional evidence beyond its presence in the crystal.')
add('P14618','ALA PHE PRO SER TRP VAL','GOLD','PKM2 allosteric regulation',
    'Biochemical experiments establish activity modulation by the free amino acid at physiological concentration ranges; the AA pocket mediates allostatic regulation.','29748232 23064226 23530218')
add('P14618','ASN ASP','GOLD','PKM2 allosteric activation',
    'Kinetic, fluorescence-binding, crystallographic and gel-filtration experiments demonstrate activation of PKM2 by free asparagine/aspartate.','32144209')
add('P14618','CYS','GOLD','PKM2 allosteric inhibition',
    'Wild-type PKM2 kinetic and ligand-binding experiments directly establish concentration-dependent cysteine inhibition and coupling to the catalytic site.','31386812')
add('P14618','THR','STRUCTURAL_ONLY','fragment-screen ligand',
    'Threonine is present in a fragment-screen structure; the reviewed amino-acid regulation studies do not establish threonine as a functional PKM2 regulator.','31877353 29748232')
add('Q7L266','GLY','GOLD','autoproteolytic activation',
    'Free glycine specifically accelerates human ASRGL1 autoproteolytic maturation in vitro and in human cells; other small amino acids do not substitute.','23601642')
add('Q71RI9','GLN','GOLD','glutamine transamination',
    'Glutamine is a directly assayed free substrate and competitive ligand of mouse KAT3; independent mouse-kidney work supports its glutamine-transaminase role.','19029248')
add('P02768','TRP','GOLD','plasma transport ligand',
    'Binding of free L-tryptophan to human albumin is measured at physiological pH and temperature; this supports its plasma carrier role, independent of the PDB drug-complex context.','3099525')

# Positive biochemical binding with limited physiological specificity/context.
for accession, aas, reason in [
    ('P00918','HIS PHE TRP','CA activation measured in vitro, but physiological free-AA regulation is not established sufficiently for the primary positive set.'),
    ('P00915','HIS','Histidine activation is measured in vitro; endogenous ligand function is weaker than for cognate substrate/receptor positives.'),
    ('P05187','PHE','Free phenylalanine inhibits placental alkaline phosphatase, but its native regulatory role is proposed and refined models altered ligand placements.'),
    ('P05089','LYS','Lysine is a noncognate arginase ligand/inhibitor, rather than its physiological arginine substrate.'),
    ('Q16773','PHE','Phenylalanine is a competing secondary substrate of KAT1; its physiological role is weaker than the principal kynurenine function.'),
    ('Q9S7N2','PHE','Phenylalanine activity is reported for TAA1 in vitro; the established physiological pathway uses tryptophan for auxin synthesis.'),
    ('Q9SKE2','LEU VAL','JAR1 can conjugate these amino acids in vitro, but the demonstrated native signaling product is JA-isoleucine.'),
    ('Q84VW9','ASP','PEPC-family aspartate feedback is plausible, but the deposited AtPPC3 ligand complex lacks a primary publication and isoform-specific evidence is insufficient.'),
    ('B7M7I1','LYS','Lysine feedback is established for E. coli DHDPS, but this strain-specific unreviewed accession/deposition lacks sufficient independent functional evidence.'),
    ('P0ADF8','ILE','Isoleucine binds IlvN, but the demonstrated physiological feedback specificity is stronger for valine; differential binding alone does not establish equivalent regulatory function.'),
    ('P04816','PHE','Phenylalanine binding is structurally demonstrated for the leucine-specific binding protein, but native phenylalanine uptake is not established here.'),
    ('P05041','TRP','Tryptophan remains tightly bound and may support PabB structural integrity, but the publication finds no effect on catalytic activity and proposes its functional role.'),
    ('P0AGL2','SER','Serine is one of several ligands of TdcF; the proposed functional metabolite is 2-ketobutyrate, not established free-serine recognition.'),
    ('Q8IWU9','PHE','Binding to the TPH2 regulatory domain is demonstrated, but the native phenylalanine-dependent regulatory mechanism remains unresolved.'),
    ('P17643','TYR','Tyrosine binding/activity for human TYRP1 is not sufficient to override its dominant DHICA function; human physiological tyrosine recognition remains uncertain.')]:
    add(accession,aas,'SILVER','secondary/proposed functional ligand',reason)
add('Q9SKE2','ILE','GOLD','JA-Ile biosynthetic substrate',
    'Free L-isoleucine is ATP-activated and conjugated to jasmonate by JAR1 to generate the native JA-Ile signaling molecule.')

# A free product is not automatically an established physiological free ligand.
for accession, aas in [
    ('P00805','ASP GLU'),('P00903','GLU'),('P00963','ASN'),
    ('P06988','HIS'),('P0A962','ASP'),('P17169','GLU'),
    ('O94925','GLU'),('P15104','GLN'),('P19440','GLU'),
    ('P20933','ASP'),('Q06210','GLU'),('P78330','SER'),
    ('O50008','MET'),('Q9LEU8','ARG'),('Q9M0A7','GLU'),
    ('Q9ZVL6','SER'),('D3Z7P3','GLU'),('Q8CHT0','GLU'),
    ('P18956','GLU'),('P37595','ASP'),('Q7L266','ASP')]:
    add(accession,aas,'SILVER','reaction-product binding',
        'The free amino acid is a documented enzymatic product and product-bound structures are mechanistically plausible. Product occupancy alone does not establish an independently functional free-AA ligand for the strict biological positive set.')
for accession in ['O43766','O82392','P36979','Q8CBB9','Q9H9T3']:
    add(accession,'MET','SILVER','radical-SAM reaction byproduct',
        'Methionine is released by cleavage of S-adenosylmethionine; SAM recognition does not itself establish physiological recognition of free methionine as a ligand.')
add('P07813','ILE MET','SILVER','aaRS post-transfer editing product',
    'The noncognate AA occupies an editing domain that hydrolyzes mischarged aminoacyl-tRNA. This differs from the cognate free-AA aminoacylation substrate.','16277600')
add('P0A8M3','SER','SILVER','aaRS post-transfer editing product',
    'The deposited serine is an editing-site product of seryl-tRNA proofreading; direct free-serine misactivation is much weaker than cognate threonine charging.','15525511 10881191')

# Explicit biological context incompatibilities; these are not absence-of-evidence
# EXCLUDE rules based on generic annotation names.
add('A6NHX0','ARG','EXCLUDE','nonbinding paralog',
    'UniProt explicitly states CASTOR2 does not directly bind arginine; CASTOR1 in the heterodimer is the arginine-binding sensor.')
add('P0A8P1','LEU PHE','EXCLUDE','aminoacyl-tRNA donor recognition',
    'The native donor is an aminoacyl-tRNA, not the free amino acid; free-AA analog occupancy does not establish a native free-AA ligand.')
add('P0CE48','ILE PHE','EXCLUDE','aminoacyl-tRNA recognition',
    'EF-Tu recognizes aminoacyl-tRNA, not these free amino acids as physiological ligands.')
add('P41091','MET','EXCLUDE','initiator tRNA recognition',
    'eIF2 recognizes methionyl initiator tRNA; methionine residue recognition is not functional free-methionine binding.')
add('P0AD89 P61175','TRP','EXCLUDE','ribosome-nascent-peptide composite sensor',
    'Free tryptophan sensing is real, but requires the TnaC-ribosome composite pocket. This is not a standalone protein-AA recognition control; do not attribute the complex function to isolated TnaC or L22.','25310980 34504068 34403461')
add('P0A6N4','PRO','EXCLUDE','polyproline translation context',
    'EF-P assists translation of proline-rich peptide motifs; that function does not establish free-proline ligand recognition.')
add('P61964','ARG','EXCLUDE','histone-peptide recognition mimic',
    'Free arginine is used to dissect WDR5 recognition of histone H3 Arg2 methylation states; the native functional ligand is a histone peptide.','32207970')
add('Q13501','ARG','EXCLUDE','N-degron peptide recognition',
    'The primary study measures free-arginine binding but shows that free arginine does not serve as the direct mTORC1 sensor or induce the autophagy response driven by arginylated peptide/protein substrates. Binding alone therefore does not establish a functional free-AA positive control.','30349045')
add('Q16769','GLN LEU TYR','EXCLUDE','peptide-associated cyclase complex',
    'The deposited context is glutaminyl cyclase acting on neurotensin peptide. Native N-terminal peptide cyclization is distinct from free canonical-AA ligand recognition.','32551535')
add('P30419','ASN LYS PRO SER','EXCLUDE','peptide-associated myristoylation context',
    'NMT1 recognizes peptide/protein substrates for myristoylation; the canonical-AA entries do not establish functional binding of these free amino acids.')
add('Q8CHT0','PRO','EXCLUDE','crystallographic cryoprotectant',
    'The primary paper explicitly uses 2-3 M proline as a protein-crystallography cryoprotectant, including active-site occupancy.','22868767')
add('P02768','VAL','EXCLUDE','synthetic metal-complex ligand',
    'Valine is a co-ligand of a synthetic nitrosylruthenium complex bound to albumin; this is not free canonical valine recognition.','38263542')
add('Q15758','ASP','EXCLUDE','chimeric transporter specificity',
    'The candidate includes EAAT/ASCT transporter chimeras. Aspartate recognition cannot be assigned as native free-AA specificity of human ASCT2 from the fusion accession.')
add('P00968','GLN','EXCLUDE','wrong subunit attribution',
    'Carbamoyl-phosphate synthetase free glutamine is bound and hydrolyzed by the small subunit CarA; CarB uses ammonia and must not inherit CarA glutamine specificity.')

# N-terminal peptidase product/analog structures are withheld. Free-AA feedback
# roles would need pair-specific biochemical evidence before promotion to GOLD.
for accession, aas in [
    ('P04825','ALA ARG ASP GLU LEU LYS MET PHE SER TRP TYR'),
    ('P15034','LEU PRO'),('P0AE18','MET'),('P50579','MET'),
    ('P53582','MET'),('Q6UB28','MET'),('P12955','GLY LEU PRO'),
    ('Q03154','GLY'),('Q04609','ASP GLU MET'),('Q07075','ARG ASP GLU'),
    ('Q6P179','LYS'),('Q9UIQ6','GLY LYS'),('Q96IY4','ARG GLY'),
    ('P48052','GLU'),('Q9UI42','ASP VAL'),('Q9Y3Q0','GLU'),
    ('P39377','ASN ASP'),('Q8IYS1','HIS')]:
    add(accession,aas,'STRUCTURAL_ONLY','peptidase product/analog context',
        'The established function recognizes a peptide or modified amino-acid substrate. The free canonical-AA structural entry alone does not establish a physiological free-AA regulatory or substrate-recognition role.')

# Additional accession-specific decisions after checking the remaining inventory.
add('Q8N8M0','HIS','GOLD','histidine acetylation substrate',
    'A 2025 primary study identifies NAT16/HisAT as the human enzyme acetylating free histidine in vitro and in vivo; binding/kinetics and a plasma-associated variant corroborate this function.','40595645')
add('Q8N8M0','ARG','SILVER','noncognate HisAT ligand',
    'Arginine occupies HisAT in a substrate-comparison structure, but the demonstrated native acetylation substrate is histidine; arginine-specific physiological function is not established.','40595645')
add('Q8IWV7 Q8IWV8','TRP','MANUAL_REVIEW','unpublished UBR ligand complex',
    'The new tryptophan-bound UBR deposition has no linked primary publication. Native N-degron recognition and known leucine sensing do not establish free-tryptophan sensing by this protein.')
add('P93836','PHE','STRUCTURAL_ONLY','noncognate HPPD ligand',
    'The established HPPD substrate is 4-hydroxyphenylpyruvate, not free phenylalanine. The ligand deposition has no linked published functional evidence for free Phe.')
add('Q9SKE2','MET','SILVER','secondary JAR1 ligand',
    'Methionine is used in a FIN219/JAR1 substrate-comparison complex; established native jasmonate signaling is through JA-isoleucine, and methionine-specific physiological conjugation is not established.','28223489')
add('P50440','ALA','EXCLUDE','arginine substrate analog',
    'The primary paper explicitly studies alanine as an amidino-donor analog, distinct from the functional free arginine/glycine substrates.','9915841')
add('P75823','SER','SILVER','nonproductive secondary substrate',
    'The primary paper reports free serine in a nonproductive orientation under low-pH conditions; native serine turnover/function is less established than the threonine reaction.','25560296')
add('P00962','GLU','EXCLUDE','noncognate discrimination analog',
    'GlnRS binds glutamate in a nonproductive orientation to discriminate against this noncognate amino acid; it is not the genuine glutamine aminoacylation substrate.','12691748 15845536')
add('P0A9Q9','CYS','EXCLUDE','cysteine-derived inhibitor context',
    'The ASADH deposition concerns a cysteine-derived covalent substrate analog, not the documented physiological free aspartate-semialdehyde substrate.','11724560')
add('P05041','GLU','EXCLUDE','partner glutaminase product context',
    'Glutamate is generated in PabA glutamine hydrolysis. PabB consumes chorismate and transferred ammonia; complex proximity does not assign free-glutamate recognition to PabB.','40365074')
add('P37595','GLY','MANUAL_REVIEW','uncertain glycine activation context',
    'The candidate is an EcAIII mutant structure. Glycine-dependent activation established for human ASRGL1 cannot be transferred to this E. coli enzyme without pair-specific functional evidence.','37095066')
add('P0AD61','GLY','MANUAL_REVIEW','unpublished mutant PK ligand',
    'The deposition concerns a G381A pyruvate-kinase mutant without a linked publication. Glycine-dependent native WT activity is not established by the deposited ligand alone.')
add('Q9UDR5','PRO','STRUCTURAL_ONLY','noncognate lysine-catabolism ligand',
    'AASS uses lysine/saccharopine-pathway substrates. The unpublished proline-bound domain does not establish native free-proline recognition.')
add('Q92600 Q9UKV8','TRP','EXCLUDE','GW-rich protein motif mimic',
    'Tryptophan occupies pockets used to recognize tryptophan-containing GW182/TNRC6 protein motifs; native peptide/protein recognition is distinct from a free-Trp ligand role.','24768540 24768538 22539551 29576456')
add('P41052','ALA','EXCLUDE','muropeptide-associated amino acid',
    'The structure contains peptidoglycan muropeptides with alanine residues; this is not functional recognition of free canonical alanine.')
add('Q9SJQ9','GLY PRO','EXCLUDE','Gly-Pro dipeptide-associated',
    'The primary deposition explicitly concerns the Gly-Pro dipeptide bound to aldolase; activity of a dipeptide does not establish activity of either free amino acid.')
add('A0A8V8TRG9 P62942','GLU','EXCLUDE','fusion-domain accession attribution',
    'The accession refers to a fused mTOR/FKBP domain in a metabotropic-glutamate-receptor construct. Native receptor Glu specificity must not be assigned to the fused partner.')
for accession,aas in [('E9PMV2','GLY'),('P01920','GLY'),('O78189','ARG'),
    ('F6IQS1','GLY LEU MET'),('Q95HB9','GLY'),('P01901','CYS GLY LEU')]:
    add(accession,aas,'EXCLUDE','MHC peptide/dipeptide presentation',
        'The deposited antigen-presentation context uses peptides/dipeptides. Recognition of their amino-acid residues does not establish a native functional free-AA ligand.')
for accession,aas in [('O75112','GLU LEU'),('Q9Y566','ARG LEU'),
    ('O15460','GLY'),('P13674','GLY'),('P06756','GLY'),
    ('Q14974','PHE'),('Q14498','LYS'),('Q96I25','TRP')]:
    add(accession,aas,'EXCLUDE','protein/peptide motif recognition',
        'The ligand context is a recognized protein/peptide motif or peptide hybrid; there is no established free canonical-AA function for the evaluated pair.')
add('P12821','ASP GLN GLU GLY LYS PRO SER VAL','EXCLUDE',
    'peptide/peptidomimetic ACE ligand',
    'The ACE deposits use amyloid-beta, Lys-Pro, or antihypertensive peptidomimetics. Native peptide/dipeptide processing does not establish free-AA specificity.')

# Final review of every previously unassigned functional annotation/citation.
add('A0A140N890','GLU','GOLD','free-L-glutamate racemization substrate',
    'Direct activity assays establish L-glutamate racemization by this E. coli BL21 enzyme, with threefold preference over L-aspartate; the bound substrate structure explains specificity.','26555188')
add('A0A0H3JGH6','ASP','GOLD','free-L-aspartate racemization substrate',
    'The O157-specific primary study reports racemase activity and L-/D-aspartate complexes supporting the mechanism; this is not inferred solely from a homolog.','27001440')
add('P33590','HIS','EXCLUDE','nickel-chelate recognition',
    'NikA recognizes the Ni-(L-His)2 metal chelate. Histidine coordination within a nickel cargo is not recognition of a free uncoordinated amino acid.')
add('P32929','MET','EXCLUDE','engineered nonnative substrate specificity',
    'Methionine degradation and binding are reported for the engineered E59N/R119L/E339V CGL variant. This does not establish free methionine as a native WT human CGL substrate.','28106980')
add('P30793','PHE','EXCLUDE','regulatory-partner attribution',
    'Phenylalanine is recognized by GFRP (P30047) in the stimulatory GCH1-GFRP complex; proximity to GCH1 does not establish a standalone GCH1 free-Phe binding site.','33229582')
add('O75884','PHE','STRUCTURAL_ONLY','aminopeptidase product context',
    'RBBP9 removes aromatic N-terminal amino acids from peptides in human cells. A free phenylalanine product complex does not establish an independent native free-Phe substrate or regulatory function.','35173328')
add('P00441 Q9UK05','CYS','EXCLUDE','cysteinylation modification context',
    'The deposited ligand context explicitly concerns protein cysteinylation, rather than demonstrated reversible functional recognition of free cysteine.')
add('P0AEX9','GLU GLY SER','EXCLUDE','fusion-tag attribution',
    'The maltose-binding protein is a fusion tag in these deposits. Canonical-AA specificity of the fused target or crystallization additive must not be assigned to the tag.')
add('Q8NFU1','GLU','GOLD','intracellular glutamate channel regulation',
    'A primary study identifies the intracellular glutamate-binding site and demonstrates glutamate-dependent regulation of BEST2 together with glutamine synthetase in vivo; this is functional sensing rather than an anion-channel substrate inference.','39737942')
