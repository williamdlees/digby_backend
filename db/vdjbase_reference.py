import os.path
import csv
import re
import os
from Bio import SeqIO
import importlib.util

from db.vdjbase_model import Gene, Allele, AlleleAliasSets
from db.vdjbase_exceptions import DbCreationError


def read_reference(filename):
    records = {}

    for rec in SeqIO.parse(filename, 'fasta'):
        if '|' in rec.description:
            rd = rec.description.split('|')[1]      # assume IMGT convention
        else:
            rd = rec.description

        records[rd] = rec.seq.lower()

    return records


snp_in_name = re.compile(r'^[acgtACGT]([0-9]+)[acgtACGT-]$')
ins_in_name = re.compile(r'^([0-9]+)[acgtACGT]+([0-9]+)$')


def mark_snp_positions(allele_name, sequence):
    """ Lower-case the SNP and insertion positions named in the allele, as the
        pipeline does. Alleles are matched on sequence, so an all-upper reference
        never matches the same allele arriving from a repertoire. """
    if '_' not in allele_name:
        return sequence

    positions = []
    for tok in allele_name.split('_')[1:]:
        snp = snp_in_name.match(tok)
        ins = ins_in_name.match(tok)
        if snp:
            positions.append(int(snp.group(1)))
        elif ins:
            positions.extend(range(int(ins.group(1)), int(ins.group(2)) + 1))

    if not positions:
        return sequence

    seq = list(sequence)
    for pos in positions:
        if 0 < pos <= len(seq):
            seq[pos - 1] = seq[pos - 1].lower()
    return ''.join(seq)


def read_reference_table(filename):
    """ iuis_allele, asc, sequence, gapped_sequence -> {allele: (asc, gapped_sequence)} """
    recs = {}
    with open(filename, newline='') as fi:
        for row in csv.DictReader(fi, delimiter='\t'):
            name = row['iuis_allele'].strip()
            if name:
                seq = (row['gapped_sequence'] or row['sequence']).strip()
                recs[name] = (row['asc'].strip() or None, mark_snp_positions(name, seq))
    return recs


def asc_of_ambiguous(allele_name, table):
    """ The cluster of an ambiguous name like IGHV1-69*01_12_13, if the alleles it
        spans agree on one. """
    if '*' not in allele_name:
        return None

    gene, rest = allele_name.split('*', 1)
    groups = {table[f'{gene}*{part}'][0] for part in rest.split('_')
              if part.isdigit() and f'{gene}*{part}' in table}

    if len(groups) > 1:
        raise DbCreationError(f'{allele_name} spans alleles in different clusters: {sorted(groups)}')

    return groups.pop() if groups else None


def import_reference_alleles(reference_dir, session, species):
    result = []
    if os.path.isfile(os.path.join(reference_dir, 'gene_order.py')):
        spec = importlib.util.spec_from_file_location("gene_order", os.path.join(reference_dir, 'gene_order.py'))
        gene_order = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gene_order)
    else:
        raise DbCreationError('gene_order.py not found - skipped')

    added_genes = []
    extra_locus = len(gene_order.LOCUS_ORDER)
    extra_alpha = len(gene_order.ALPHA_ORDER)

    # the reference table, then the FASTA for anything it does not name: the FR1
    # ambiguous names (IGHV1-69*01_12_13) are only in the FASTA, and the genotype
    # and ogrdbstats loaders match them as literal strings
    table = {}
    for file in sorted(os.listdir(reference_dir)):
        if file.startswith('reference_table_') and file.endswith('.tsv'):
            table.update(read_reference_table(os.path.join(reference_dir, file)))

    from_table = from_fasta = 0

    for allele, (asc, sequence) in table.items():
        gene_name = allele.split('*')[0] if '*' in allele else allele

        if gene_name not in added_genes:
            (extra_alpha, extra_locus) = add_gene(extra_alpha, extra_locus, gene_name, gene_order, session, species)
            added_genes.append(gene_name)

        save_allele(allele, gene_name, sequence, session, asc=asc, asc_inferred=False)
        from_table += 1

    for file in os.listdir(reference_dir):
        if os.path.splitext(file)[1] == '.fasta':
            recs = read_reference(os.path.join(reference_dir, file))

            for allele, sequence in recs.items():
                if allele in table:
                    continue

                gene_name = allele.split('*')[0] if '*' in allele else allele

                if gene_name not in added_genes:
                    (extra_alpha, extra_locus) = add_gene(extra_alpha, extra_locus, gene_name, gene_order, session, species)
                    added_genes.append(gene_name)

                save_allele(allele, gene_name, sequence, session, asc=asc_of_ambiguous(allele, table))
                from_fasta += 1

        session.commit()

    session.commit()

    if len(added_genes) == 0:
        raise DbCreationError('No genes added from reference set - skipped')

    result.append(f'Reference alleles added: {from_table} from the reference table, '
                  f'{from_fasta} from FASTA, {len(added_genes)} genes')
    return result


def save_allele(allele_name, gene_name, sequence, session, asc=None, asc_inferred=None):
    similar = session.query(Allele).filter(Allele.seq == str(sequence)).one_or_none()

    if similar is not None:
        if similar.similar is None or len(similar.similar) == 0:
            similar.similar = '|%s|' % allele_name
        else:
            similar.similar += ', ' + '|%s|' % allele_name

        # a name collapsed into `similar` still belongs to a cluster
        if asc and not similar.asc:
            similar.asc = asc
            similar.asc_inferred = asc_inferred
    else:
        g = session.query(Gene).filter(Gene.name == gene_name).one_or_none()
        a = Allele(
            name=allele_name,
            seq=str(sequence),
            seq_len=str(len(sequence)),
            appears=0,
            gene_id=g.id,
            is_single_allele=True,
            low_confidence=False,
            novel=False,
            max_kdiff=0,
            similar='',
            pipeline_name='',
            asc=asc,
            asc_inferred=asc_inferred,
        )
        session.add(a)
    session.flush()


def add_gene(extra_alpha, extra_locus, gene, gene_order, session, species):
    if gene in gene_order.LOCUS_ORDER:
        locus_order = gene_order.LOCUS_ORDER.index(gene)
    else:
        locus_order = extra_locus
        extra_locus += 1
    if gene in gene_order.ALPHA_ORDER:
        alpha_order = gene_order.ALPHA_ORDER.index(gene)
    else:
        alpha_order = extra_alpha
        extra_alpha += 1
    g = Gene(
        name=gene,
        type=gene[0:4],
        family=gene.split('-')[0] if '-' in gene else gene,
        species=species,
        locus_order=locus_order,
        alpha_order=alpha_order,
        pseudo_gene=1 if gene in gene_order.PSEUDO_GENES else 0
    )
    session.add(g)
    return (extra_alpha, extra_locus)


# scan subdirs of reference_dir for alias sets
# add an alias set if we find one or more .fasta files in a subdir
def add_alias_sets(reference_dir, session):
    results = []
    alias_num = 1

    alleles = session.query(Allele).all()
    ungapped_seqs = {a.seq.replace('.', ''): a for a in alleles}

    for subdir in os.listdir(reference_dir):
        subdir_path = os.path.join(reference_dir, subdir)
        if os.path.isdir(subdir_path):
            fasta_files = [f for f in os.listdir(subdir_path) if f.endswith('.fasta')]
            if fasta_files:
                alias_set = AlleleAliasSets(
                    set_name=subdir,
                    alias_number=alias_num
                )
                session.add(alias_set)
                alias_name = f'alias_{alias_num}'
                alias_num += 1
                results.append(f'Processing alias set {alias_num}')

                for fasta_file in fasta_files:
                    recs = read_reference(os.path.join(reference_dir, subdir_path, fasta_file))

                    for allele, sequence in recs.items():
                        if sequence in ungapped_seqs:
                            similar = ungapped_seqs[sequence]

                            # add the allele name to the alias_name column of the similar allele
                            sn = getattr(similar, alias_name)
                            if sn is None or len(sn) == 0:
                                setattr(similar, alias_name, allele)
                            else:
                                setattr(similar, alias_name, f'{sn}|{allele}')

                session.commit()

    results.append('Alias sets processed')
    return results
 