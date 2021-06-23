##########################################################################################################
#
#       Confidence checks on novel alleles
#
#       This might be better as part of the pipeline, but, at least while the checks are under development,
#       they are run when the alleles view is rendered for the first time, so that they can be modified
#       without re-running the pipeline. To re-run the checks, run the command below from sqlite command line
#       UPDATE database_alleles SET novel=False;
#
##########################################################################################################
import sys
import os.path

from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from sqlalchemy import not_, distinct, or_

from db.vdjbase_model import Sample, Allele, AllelesSample, Gene, GenesDistribution, AllelesPattern, AlleleConfidenceReport, SNP, HaplotypeEvidence, SamplesHaplotype, SeqProtocol

import re
import csv


def check_novel_confidence(ds_dir, session):
    result = ['Running confidence checks']
    novels = gather_novel_data(session)

    gather_haplo_data(novels, ds_dir, session)

    #   Check for deletion on other chromosome, or novel allele present on both chromosomes
    check_for_singleton_infs(novels, ds_dir, session)

    # Mark an allele as low confidence only if the only sightings we have are on both chromosomes
    # Important this runs before we look at any other low-confidence markers

    fix_haplo_confidence(novels, session)

    #   Check for zXXX, XzXX, XXzX, XXXz -> XXXX (where X and z are any base)
    check_XXXX(novels, session)

    #   Check for an amino acid at a location that has not been previously observed in that family
    result.extend(check_novel_aa(novels, session))

    #   Check for RGYW/WRCY hotspot change
    check_hotspot(novels, session)

    #   Check for FR1 SNP in a sample with FR1 primers
    check_snp_in_short(novels, session)

    #   Check for SNP in the same location in >1 novel allele from the same sample
    check_shared_snp(novels, session)

    session.commit()

    # Produce some stats on CNV (not routinely used)
    # allele_copy_stats()

    return result


# Collect info on novel inferences, store in Alleles table
def gather_novel_data(session):
    novels = session.query(Allele).filter(Allele.novel == True).all()
    for novel in novels:
        closest_ref = session.query(Allele).filter(Allele.name == novel.name.split('_')[0]).one_or_none()

        if not closest_ref:
            closest_ref = session.query(Allele).filter(Allele.similar.like('%|' + novel.name.split('_')[0] + '|%')).one_or_none()

        if not closest_ref:
            report_issue(novel, 'No Closest Ref',
                         'The closest reference allele, %s, is not in the database' % (novel.name.split('_')[0]), session)
        elif not closest_ref.seq:
            report_issue(novel, 'No Ref Sequence',
                         'The closest reference allele, %s, has no defined sequence' % (novel.name.split('_')[0]), session)
        else:
            novel.closest_ref = closest_ref
            session.flush()

        snps = novel.name.split('_')

        for snp in snps:
            if snp[0].isalpha() and snp[-1].isalpha():
                s = SNP(
                    from_base=snp[0],
                    pos=int(snp[1:-1]),
                    to_base=snp[-1],
                    allele=novel
                )
                session.add(s)
    return novels

"""
# One-off function, used for a single study, to assess gene duplication. Not used at the moment.
def allele_copy_stats():
    samples = Samples.objects.all()
    HaplotypeEvidence.objects.all().delete()

    allele_report_samples = {}

    for sample in samples:
        if sample.seq_protocol.sequencing_length == 'Full':
            for haplo_file in sample.haplotype.all():
                with open(os.path.join(BASE_DIR, 'uploads', haplo_file.file.name), 'r', newline='') as fi:
                    reader = csv.reader(fi, dialect='excel-tab')
                    header = True
                    for row in reader:
                        if header:
                            header = False
                            continue
                        h1 = set(row[3].split(',')) - {'Unk', 'Del'}
                        h2 = set(row[4].split(',')) - {'Unk', 'Del'}
                        contains_novel = False
                        for haplotyped in (h1, h2):
                            for allele in haplotyped:
                                if '_' in allele:
                                    contains_novel = True
                                    break
                            if not contains_novel and len(haplotyped) > 0:
                                if row[2] not in allele_report_samples:
                                    allele_report_samples[row[2]] = {}
                                if sample.name not in allele_report_samples[row[2]]:
                                    allele_report_samples[row[2]][sample.name] = len(haplotyped)
                                else:
                                    allele_report_samples[row[2]][sample.name] = max(allele_report_samples[row[2]][sample.name], len(haplotyped))

    allele_report = {}
    allele_report_samples3 = {}
    for gene, counts in allele_report_samples.items():
        totals = [0]*5
        for sample, count in counts.items():
            totals[count] += 1
            if count > 2:
                if gene not in allele_report_samples3:
                    allele_report_samples3[gene] = []
                allele_report_samples3[gene].append(sample)
        allele_report[gene] = totals

    with open('allele_report.csv', 'w') as fo:
        for k in sorted(allele_report.keys()):
            fo.write("%s,%s\n" % (k, ', '.join([str(x) for x in allele_report[k]])))

    with open('allele_report_samples.csv', 'w') as fo:
        for k in sorted(allele_report_samples3.keys()):
            fo.write("%s,%s\n" % (k, ', '.join(allele_report_samples3[k])))

"""


# Collect haplotyping info on novel alleles, store in HaplotypeEvidence
def gather_haplo_data(novels, ds_dir, session):
    sample_haplotypes = session.query(SamplesHaplotype).all()
    alleles = session.query(Allele.name, Allele.pipeline_name).all()

    vdjbase_allele = {}

    for allele in alleles:
        if allele.pipeline_name:
            for pa in allele.pipeline_name.split(', '):
                pa = pa.split('*')[0] + '*' + pa.split('*')[1].lower()
                vdjbase_allele[pa] = allele.name

    for sample_haplotype in sample_haplotypes:
        print(sample_haplotype.haplotypes_files.file)
        with open(os.path.join(ds_dir, sample_haplotype.haplotypes_files.file), 'r', newline='') as fi:
            reader = csv.reader(fi, dialect='excel-tab')
            header = True
            row_index = True

            for row in reader:
                if header:
                    header = False
                    # check whether or not rows have an R style index number
                    try:
                        row_index = float(row[0])
                    except:
                        row_index = False
                    continue

                if row_index:
                    row = row[1:]

                haplotyped = set(row[2].split(',')) ^ set(row[3].split(','))
                for allele in haplotyped:
                    if '_' in allele:
                        name = row[1] + '*' + allele.lower()

                        if name in vdjbase_allele:
                            name = vdjbase_allele[name]

                        n = session.query(Allele).filter(Allele.name == name).one_or_none()

                        if n:
                            allele_names = row[4].split(',')
                            allele_counts = [row[7], row[9], row[11], row[13]]
                            counts = []

                            for i in range(len(allele_names)):
                                counts.append('%s (%s)' % (allele_names[i], allele_counts[i]))

                            he = HaplotypeEvidence(
                                hap_gene=sample_haplotype.haplotypes_files.by_gene,
                                sample_id=sample_haplotype.samples_id,
                                allele=n,
                                counts=', '.join(counts)
                            )
                            session.add(he)
                        else:
                            if 'ap01_g291c_t296c_c314t' in name:
                                print('foo')
                            print('foo')
                            print('Allele %s is present in haplotype file %s but is not in the Alleles table' % (name, sample_haplotype.haplotypes_files.file))

    for novel in novels:
        if session.query(HaplotypeEvidence).filter(HaplotypeEvidence.allele_id == novel.id).count():
            detects = []
            for sample_id in set([h.sample_id for h in session.query(HaplotypeEvidence).filter(HaplotypeEvidence.allele_id == novel.id).all()]):

                hap_genes = []
                for h in session.query(HaplotypeEvidence).filter(HaplotypeEvidence.sample_id == sample_id).all():
                    hap_genes.append('%s (%s)' % (h.hap_gene, h.counts))

                detects.append('%s(%s)' % (session.query(Sample).filter(Sample.id == sample_id).one_or_none().name, ', '.join(hap_genes)))

            report_issue(novel, 'Confirmation from Haplotype', 'Inferred allele on one chromosome only in sample(s) %s.' % (', '.join(detects)), session, low_confidence=False)

def check_XXXX(novels, session):
    for novel in novels:
        if novel.closest_ref:
            ref_nt = novel.closest_ref.seq.replace('.', '')
            seq_nt = novel.seq.replace('.', '')

            # walk up each identified repeat of 4nt or more, flag any differences
            seq_qpos = [m.start() for m in re.finditer('(.)\\1+\\1+\\1+', str(seq_nt))]
            q_runs = []

            for p in seq_qpos:
                rep_c = seq_nt[p]
                i = p
                while i < len(seq_nt) and i < len(ref_nt) and seq_nt[i] == rep_c:
                    if ref_nt[i] != rep_c and rep_c != 'n':
                        r_pos = find_gapped_index(i, novel.closest_ref.seq)
                        q_runs.append("%d%s%d" % (r_pos, seq_nt[p:p+4], r_pos+3))
                        break
                    i += 1

            if len(q_runs):
                report_issue(novel, 'Repeated Read', 'Possible repeated read error(s) at %s' % (', '.join(q_runs)), session)


def check_novel_aa(novels, session):
    result = []
    ref_vh_codon_usage, res = get_ref_vh_codon_usage(session)
    result.extend(res)
    for novel in novels:
        try:
            q_codons = []
            ref_aa_gapped = Seq(novel.closest_ref.seq.upper()).translate(gap='.')
            seq_aa_gapped = Seq(novel.seq.upper()).translate(gap='.')

            family = find_family(novel.closest_ref.name)

            for i in range( min(len(ref_aa_gapped), len(seq_aa_gapped))):
                # if ref_aa_gapped[i] != seq_aa_gapped[i] and '*' not in (ref_aa_gapped[i], seq_aa_gapped[i]) and '.' not in (ref_aa_gapped[i], seq_aa_gapped[i]):
                if ref_aa_gapped[i] != seq_aa_gapped[i] and seq_aa_gapped[i] not in ['*', 'X']:
                    if seq_aa_gapped[i] not in ref_vh_codon_usage[family][i+1]:
                        q_codons.append("%s%d" % (seq_aa_gapped[i], i+1))

            if len(q_codons) > 0:
                report_issue(novel, 'Novel AA', 'Amino Acid(s) previously unreported in this family - %s' % ", ".join(q_codons), session, low_confidence=len(q_codons) > 1)

        except TranslationError:
            result.append('Warning: in sequence %s or %s: %s' % (novel.name, novel.closest_ref.name, sys.exc_info()[1]))

    return result


def check_hotspot(novels, session):
    for novel in novels:
        ref_nt = novel.closest_ref.seq.replace('.', '')
        seq_nt = novel.seq.replace('.', '')
        ref_qpos = [m.start() for m in re.finditer('[ag][g][ct][at]', ref_nt)]

        q_hotspots= []

        for p in ref_qpos:
            if len(seq_nt) > p and seq_nt[p+1] == 'c':
                q_hotspots.append("%s%d%s" % (ref_nt[p+1], find_gapped_index(p+1, novel.closest_ref.seq), seq_nt[p+1]))

        ref_qpos = [m.start() for m in re.finditer('[at][ag][c][ct]', ref_nt)]

        for p in ref_qpos:
            if len(seq_nt) > p + 1 and seq_nt[p+2] == 'g':
                q_hotspots.append("%s%d%s" % (ref_nt[p+2], find_gapped_index(p+2, novel.closest_ref.seq), seq_nt[p+2]))

        if len(q_hotspots) > 0:
            report_issue(novel, 'Hotspot SNP', "G/C SNP in RGYW/WRCY hotspot(s) - %s" % ", ".join(q_hotspots), session, low_confidence=False)


def check_snp_in_short(novels, session):
    long_study_alleles = session.query(Allele)\
        .join(AllelesSample)\
        .join(Sample)\
        .join(SeqProtocol)\
        .filter(Allele.novel == True)\
        .filter(SeqProtocol.sequencing_length == 'Full').all()
    print(len(long_study_alleles))
    print(len(novels))
    short_studies_only = list(set(novels) - set(long_study_alleles))

    for novel in short_studies_only:
        min_snp = session.query(SNP).filter(SNP.allele_id == novel.id).order_by(SNP.pos).limit(1).one_or_none()
        if min_snp is not None and min_snp.pos <= 40:
            sample_names = session.query(Sample.name)\
                .join(AllelesSample)\
                .join(Allele)\
                .filter(Allele.id == novel.id).all()
            report_issue(novel, 'Possible SNP in primer region', 'Sequence found only in short-read study/studies %s has SNP %s%d%s in first 40 nucleotides' %
                (', '.join([x[0] for x in sample_names]), min_snp.from_base, min_snp.pos, min_snp.to_base), session)


def check_shared_snp(novels, session):
    positions = {}

    for novel in novels:
        snps = session.query(SNP).filter(SNP.allele_id == novel.id).all()
        positions[novel.name] = set(['%d%s' % (s.pos, s.to_base) for s in snps])

    for novel in novels:
        novel_positions = positions[novel.name]
        issues = []

        for s_assoc in novel.samples:
            sample = s_assoc.sample
            shared_infs = []
            for a_assoc in sample.alleles:
                inf = a_assoc.allele
                if inf.novel and inf != novel:
                    shared = positions[inf.name] & novel_positions
                    if shared:
                        shared_infs.append(inf.name)
            if shared_infs:
                issues.append((', '.join(shared_infs), sample.name))

        if len(issues):
            report_issue(novel, 'Shared SNP', 'SNP(s) shared with inferences in %d sample(s).<br>Example: %s in sample %s' % (len(issues), issues[0][0], issues[0][1]), session)


# Report if an allele is present on both chromosomes, or if there is a deletion
def check_for_singleton_infs(novels, ds_dir, session):
    sample_haplotypes = session.query(SamplesHaplotype).all()
    novel_set = frozenset([n.name for n in novels])
    reported_duplicates = {}
    reported_deletions = {}

    for sample_haplotype in sample_haplotypes:
        with open(os.path.join(ds_dir, sample_haplotype.haplotypes_files.file), 'r', newline='') as fi:
            reader = csv.reader(fi, dialect='excel-tab')
            header = True
            row_index = True

            for row in reader:
                if header:
                    header = False
                    # check whether or not rows have an R style index number
                    try:
                        row_index = float(row[0])
                    except:
                        row_index = False
                    continue

                if row_index:
                    row = row[1:]

                h1 = set([row[1] + '*' + a.lower() for a in row[2].split(',')])
                h2 = set([row[1] + '*' + a.lower() for a in row[3].split(',')])

                haplotyped = h1 | h2
                novel_haplotyped = haplotyped & novel_set
                if len(novel_haplotyped):
                    deleted = False

                    for hap in haplotyped:
                        if 'del' in hap:
                            deleted = True
                            break

                    if deleted:
                        for n in novel_haplotyped:
                            novel = session.query(Allele).filter(Allele.name == n).one_or_none()
                            if novel:
                                novel_nondels = []
                                for a in haplotyped:
                                    if 'del' not in a and 'unk' not in a:
                                        novel_nondels.append(a)

                                if novel.name not in reported_deletions:
                                    reported_deletions[novel.name] = []
                                if sample_haplotype.samples.name not in reported_deletions[novel.name]:
                                    report_issue(novel, 'Deletion on other chromosome', 'In sample %s, allele(s) %s present on one chromosome, deletion on other (%s)'
                                                 % (sample_haplotype.samples.name, ', '.join(novel_nondels), sample_haplotype.haplotypes_files.by_gene), session, low_confidence=False)
                                    reported_deletions[novel.name].append(sample_haplotype.samples.name)

                    else:
                        for n in novel_haplotyped:
                            if n in h1 and n in h2:
                                novel = session.query(Allele).filter(Allele.name == n).one_or_none()
                                if novel:
                                    if novel.name not in reported_duplicates:
                                        reported_duplicates[novel.name] = []
                                    if sample_haplotype.samples.name not in reported_duplicates[novel.name]:
                                        report_issue(novel, 'Inferred allele on both chromosomes', 'In sample %s, allele %s is present on both chromosomes (%s).'
                                                     % (sample_haplotype.samples.name, n, sample_haplotype.haplotypes_files.by_gene), session, low_confidence=True)
                                        reported_duplicates[novel.name].append(sample_haplotype.samples.name)

# An allele should be low_confidence if it is exclusively found on both chromosomes
# If there are some cases in which it's only found on one, we allow it through

def fix_haplo_confidence(novels, session):
    for novel in novels:
        if novel.low_confidence:
            single_events = session.query(AlleleConfidenceReport)\
                .filter(or_(AlleleConfidenceReport.category == 'Deletion on other chromosome', AlleleConfidenceReport.category == 'Inferred allele on both chromosomes'))\
                .filter(AlleleConfidenceReport.allele_id == novel.id)\
                .count()
            novel.low_confidence = (single_events == 0)


def report_issue(novel, issue_category, issue_notes, session, low_confidence=True):
    issue = AlleleConfidenceReport(
        allele=novel,
        category=issue_category,
        notes=issue_notes)
    session.add(issue)
    if low_confidence:
        novel.low_confidence = True
    print("%s: %s: %s: %d: %.2f" % (novel.name, issue_category, issue_notes, novel.appears, novel.max_kdiff))


# find the 1-based index of a nucleotide in a gapped  sequence, given its index in the ungapped sequence
def find_gapped_index(ind_ungapped, gapped_seq):
    ind = 0
    found_nt = 0

    while found_nt < ind_ungapped:
        if gapped_seq[ind] != '.':
            found_nt += 1
        ind += 1

    return ind + 1


# find the family (subgroup) given the name. Assumes IGxxff- type format, where ff is the family
def find_family(gene):
    if '-' not in gene:
        return None

    return gene.split('-')[0][4:]


def get_ref_vh_codon_usage(session):
    result = []
    refs = session.query(Allele).filter(Allele.novel == False).filter(not_(Allele.name.like('%Del%'))).all()
    usage = {}

    for ref in refs:
        if ref.name[3] == 'V':
            family = find_family(ref.name)
            if family not in usage:
                usage[family] = [[],]

            try:
                aa_seq = Seq(ref.seq.upper()).translate(gap='.')

                for i in range(0, len(aa_seq)):
                    if len(usage[family]) <= i+1:
                        usage[family].append([])
                    if aa_seq[i] not in usage[family][i+1]:
                        usage[family][i+1].append(aa_seq[i])
            except TranslationError:
                result.append('Warning: in sequence %s: %s' % (ref.name, sys.exc_info()[1]))

    return usage, result

