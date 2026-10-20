#!/usr/bin/env python3
"""Calculate chromatin potential from RNA-seq and RNA–chromatin contact data."""
import argparse
import csv
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import to_rgba
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

GENES = ('gene_name', 'chr', 'start', 'end', 'gene_type', 'from_source', 'N_counts')
CONTACTS = ('read_id', 'pairtype', 'rna_chr', 'rna_start', 'rna_end', 'rna_strand',
            'gene_name', 'gene_type', 'from_source', 'rna_cigar', 'rna_NM', 'rna_mapq',
            'dna_chr', 'dna_start', 'dna_end', 'dna_strand', 'dna_cigar', 'dna_NM',
            'dna_mapq', 'rna_secondary_alignments', 'dna_secondary_alignments',
            'rna_other_tags', 'dna_other_tags')
ALT = re.compile(r'\(([^,();\s]+),(\d+),([+-])((?:\d+[MIDNSHP=X])+),(\d+|\*)\)')


def fail(path, line, reason):
    """Raise an input error with its file and line number."""
    raise ValueError(f'{path}:{line}: {reason}')


def integer(value, path, line, field, minimum=0):
    """Parse a bounded integer without accepting fractional or padded values."""
    try:
        n = int(value)
        if n >= minimum and str(n) == value:
            return n
    except (ValueError, TypeError):
        pass
    fail(path, line, f'{field} must be an integer >= {minimum}: {value!r}')


def count_integer(value, path, line):
    """Parse a non-negative integer count, accepting decimal notation."""
    try:
        number = Decimal(value)
        if number.is_finite() and number >= 0 and number == number.to_integral_value():
            return int(number)
    except (InvalidOperation, ValueError, TypeError):
        pass

    fail(path, line, f'N_counts must be a non-negative integer: {value!r}')


def read_table(path, columns):
    """Yield validated TSV rows as dictionaries using named columns."""
    with open(path, newline='', encoding='utf-8') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader, None)
        if header is None:
            fail(path, 1, 'missing header')
        if len(header) != len(set(header)):
            fail(path, 1, 'duplicate header columns')
        if not set(columns) <= set(header):
            fail(path, 1, f'missing columns: {sorted(set(columns) - set(header))}')
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                fail(path, line, f'expected {len(header)} fields, got {len(row)}')
            yield line, dict(zip(header, row))


def load_counts(path):
    """Read gene-level counts and reject invalid or duplicate gene keys."""
    genes = {}
    for line, row in read_table(path, GENES):
        key = (row['gene_name'], row['gene_type'])
        if not all(key) or any(not row[k] for k in ('chr', 'from_source')):
            fail(path, line, 'empty gene identifier, biotype, chromosome or source')
        start, end = (integer(row[k], path, line, k, 1) for k in ('start', 'end'))
        if start >= end:
            fail(path, line, 'gene start must be less than end')
        row['N_counts'] = count_integer(row['N_counts'], path, line)
        if key in genes:
            fail(path, line, f'duplicate gene key {key}')
        genes[key] = row
    if not genes:
        fail(path, 2, 'no gene records')
    return genes


def intervals_distance(a, b, x, y):
    """Count unoccupied bases separating two inclusive genomic intervals."""
    """[a, b] — RNA interval"""
    """[x, y] — DNA interval"""
    return max(0, x - b - 1, a - y - 1)


def point_distance(a, b, x):
    """Count bases between an inclusive interval and a genomic position."""
    """[a, b] — RNA interval"""
    """x — DNA start"""
    return max(0, a - x - 1, x - b - 1)


def parse_alternatives(value, path, line):
    """Parse UM DNA alternatives; convert 0-based starts to 1-based."""
    if not value or value == '*':
        fail(path, line, 'UM requires DNA secondary alignments')
    for item in value.split(';'):
        match = ALT.fullmatch(item)
        if not match:
            fail(path, line, f'invalid DNA secondary alignment: {item!r}')
        chrom, start, _, _, _ = match.groups()
        yield chrom, integer(start, path, line, 'secondary start') + 1  # 0-based -> 1-based


def count_contacts(path, expected, dist, genes, counts, stats):
    """Stream voted contacts, validate records and count distant mappings."""
    with open(path, newline='', encoding='utf-8') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader, None)
        if header is None or len(header) < len(CONTACTS):
            fail(path, 1, 'expected a >=23-column Voting header')
        # The Voting stage preserves field positions, including custom upstream header names.
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                fail(path, line, f'expected {len(header)} fields, got {len(row)}')
            
            v = dict(zip(CONTACTS, row[:23]))
            typ = v['pairtype']
            if typ != expected:
                fail(path, line, f'expected pairtype {expected}, got {typ!r}')
            
            if not v['read_id'] or v['read_id'] == '*':
                fail(path, line, 'missing read_id')

            key = (v['gene_name'], v['gene_type'])
            if key not in genes:
                fail(path, line, f'gene not present in RNA-seq annotation: {key}')

            if v['from_source'] != genes[key]['from_source']:
                fail(path, line, f'gene annotation source mismatch: {key}')

            if v['rna_strand'] not in ('+', '-') or v['dna_strand'] not in ('+', '-'):
                fail(path, line, 'invalid RNA/DNA strand')

            if not v['rna_chr'] or v['rna_chr'] == '*' or not v['dna_chr'] or v['dna_chr'] == '*':
                fail(path, line, 'missing RNA/DNA chromosome')

            a, b, x, y = (integer(v[k], path, line, k, 1) for k in ('rna_start', 'rna_end', 'dna_start', 'dna_end'))
            if a > b or x > y:
                fail(path, line, 'alignment start exceeds end')
            
            secondary = list(parse_alternatives(v['dna_secondary_alignments'], path, line)) if typ == 'UM' else []
            
            if typ == 'UU' and v['dna_secondary_alignments'] != '*':
                fail(path, line, 'UU must not have DNA secondary alignments')
                
            keep = dist == 0 or (v['rna_chr'] != v['dna_chr'] or intervals_distance(a, b, x, y) >= dist)
            if keep and typ == 'UM' and dist:
                keep = all(v['rna_chr'] != chrom or point_distance(a, b, pos) >= dist
                           for chrom, pos in secondary)
            stats[f'{expected}_read'] += 1
            stats[f'{expected}_passed' if keep else f'{expected}_excluded'] += 1
            if keep:
                counts[key] += 1


def write_counts(path, genes, values):
    """Write gene-level contact counts, including genes with zero contacts."""
    with open(path, 'w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t', lineterminator='\n')
        writer.writerow(('gene_name', 'gene_type', 'from_source', 'N_counts'))
        for name, biotype in genes:
            n = values.get((name, biotype), 0)
            writer.writerow((name, biotype, genes[name, biotype]['from_source'], n))


def calculate_chp(genes, contacts, N_contacts_min, gene_len_min, fdr_threshold, output):
    """Calculate chromatin potential using a two-sample proportions Z-test."""
    columns = ('gene_name', 'gene_type', 'N_counts_RNAseq_plus_1', 'gene_length', 'N_contacts_plus_1', 'chP', 'pval', 'fdr_bh')
    rows = []

    for (name, biotype), gene in genes.items():
        n_contacts = contacts.get((name, biotype), 0)
        length = int(gene['end']) - int(gene['start']) + 1

        if n_contacts <= N_contacts_min or length <= gene_len_min:
            continue

        rows.append({
            'gene_name': name,
            'gene_type': biotype,
            'N_counts_RNAseq_plus_1': gene['N_counts'] + 1,
            'gene_length': length,
            'N_contacts_plus_1': n_contacts + 1
        })

    if len(rows) >= 2:
        Nc = sum(row['N_contacts_plus_1'] for row in rows)
        Ne = sum(row['N_counts_RNAseq_plus_1'] for row in rows)

        for row in rows:
            z, pval = proportions_ztest(
                [row['N_contacts_plus_1'], row['N_counts_RNAseq_plus_1']],
                [Nc, Ne]
            )
            row['chP'], row['pval'] = float(z), float(pval)

        fdr = multipletests([row['pval'] for row in rows], method='fdr_bh')[1]

        for row, adjusted_pval in zip(rows, fdr):
            row['fdr_bh'] = float(adjusted_pval)

        rows = [row for row in rows if row['fdr_bh'] < fdr_threshold]
        rows.sort(key=lambda row: (row['fdr_bh'], row['gene_name'], row['gene_type']))

    else:
        rows = []

    with open(output, 'w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f'Chromatin Potential: {len(rows)} significant genes. Wrote {output}')
    return rows


def joinplot_chP_Ncontacts(rows, output, args):
    """Plot chromatin potential against contact counts."""
    if not rows:
        return

    test_data = pd.DataFrame(rows)
    test_data['N_contacts'] = test_data['N_contacts_plus_1'] - 1

    # The logarithmic axis requires positive contact counts.
    test_data = test_data[test_data['N_contacts'] > 0].copy()
    if test_data.empty:
        return

    if args.protein_coding_biotype:
        test_data['biotype2'] = test_data['gene_type'].apply(
            lambda value: 'mRNAs' if value == args.protein_coding_biotype else 'ncRNAs'
        )
    else:
        test_data['biotype2'] = 'RNAs'

    counts = test_data['biotype2'].value_counts()
    test_data['biotype'] = test_data['biotype2'].map(
        lambda value: f'{counts[value]} {value}'
    )
    test_data_sorted = test_data.sort_values('biotype2')

    edge_colors = {}
    edge_colors2 = {}
    default_colors = sns.color_palette()[:2]
    biotype2_dict = {}

    for biotype in test_data_sorted['biotype'].unique():
        category = biotype.split(' ', 1)[1]
        biotype2_dict[category] = biotype

        color = default_colors[0] if category in ('mRNAs', 'RNAs') else default_colors[1]
        edge_colors[biotype] = to_rgba(color, 0.3)
        edge_colors2[biotype] = to_rgba(color, 1)

    # Create the jointplot and its marginal distributions.
    g = sns.jointplot(
        data=test_data_sorted,
        x='chP',
        y='N_contacts',
        hue='biotype',
        alpha=0,
        s=50,
        linewidth=0,
        legend=False,
        height=7,
        space=0.2
    )

    g.ax_marg_x.set_yticklabels([])
    g.ax_marg_x.set_ylabel('')
    g.ax_marg_y.set_xticklabels([])
    g.ax_marg_y.set_xlabel('')

    # Logarithmic Y-axis; contact counts remain unchanged.
    g.ax_joint.set_yscale('log')

    y_min = test_data_sorted['N_contacts'].min() * 0.9
    y_max = test_data_sorted['N_contacts'].max() * 1.1
    g.ax_joint.set_ylim(y_min, y_max)

    # Hollow markers with translucent outlines.
    for category in ('mRNAs', 'ncRNAs', 'RNAs'):
        group = test_data_sorted[test_data_sorted['biotype2'] == category]

        if group.empty:
            continue

        g.ax_joint.scatter(
            x=group['chP'],
            y=group['N_contacts'],
            facecolor='none',
            edgecolor=edge_colors[biotype2_dict[category]],
            s=30,
            linewidth=1.8,
            zorder=3 if category == 'ncRNAs' else 2
        )

    # Annotate the gene with the largest number of contacts.
    top = test_data_sorted.loc[test_data_sorted['N_contacts'].idxmax()]

    x_middle = (test_data_sorted['chP'].min() + test_data_sorted['chP'].max()) / 2
    if top['chP'] > x_middle:
        xytext, ha = (-10, 0), 'right'
    else:
        xytext, ha = (10, 0), 'left'

    g.ax_joint.annotate(
        top['gene_name'],
        xy=(top['chP'], top['N_contacts']),
        xytext=xytext,
        textcoords='offset points',
        fontsize=12,
        ha=ha,
        va='center'
    )

    legend_elements = []
    for category in ('mRNAs', 'ncRNAs', 'RNAs'):
        if category not in biotype2_dict:
            continue

        biotype = biotype2_dict[category]
        legend_elements.append(plt.Line2D(
            [0], [0],
            marker='o',
            color='w',
            label=biotype,
            markerfacecolor='none',
            markersize=10,
            markeredgecolor=edge_colors2[biotype],
            markeredgewidth=1.8
        ))

    g.ax_joint.legend(
        handles=legend_elements,
        title='Biotype',
        frameon=True,
        framealpha=0.5,
        fontsize=14,
        title_fontsize=16,
        loc='lower right'
    )

    g.ax_joint.set_xlabel('Chromatin potential', fontsize=16, labelpad=10)
    g.ax_joint.set_ylabel('The total number of contacts', fontsize=16, labelpad=10)
    g.ax_joint.tick_params(axis='x', labelsize=14)
    g.ax_joint.tick_params(axis='y', labelsize=14)

    mode = args.chP_mode
    if mode == 'counts':
        title = 'Chromatin potential (counts)'
    else:
        title = f'Chromatin potential ({mode}, dist = {args.dist:,} nt)'

    subtitle = (
        f'FDR < {args.fdr_threshold} | '
        f'N contacts > {args.N_contacts_min} | '
        f'Gene length > {args.gene_len_min} nt'
    )

    # Reserve space above the marginal distribution.
    g.figure.subplots_adjust(top=0.88)

    g.figure.suptitle(
        f'{title}\n{subtitle}',
        fontsize=14,
        y=0.98
    )

    g.figure.savefig(
        output,
        format='png',
        dpi=360,
        bbox_inches='tight'
    )
    plt.close(g.figure)

    print(f'Wrote {output}')


def arguments():
    """Parse CLI options and check mode-specific file requirements."""
    p = argparse.ArgumentParser(description='Calculate Chromatin Potential from RNA-seq counts and RNA–chromatin contacts.')
    p.add_argument('--chP_mode', required=True, choices=('counts', 'unique_dna', 'all_dna'))
    p.add_argument('--rnaseq-counts', required=True, type=Path)
    p.add_argument('--contact-counts', type=Path)
    p.add_argument('--contacts-uu', type=Path)
    p.add_argument('--contacts-um', type=Path)
    p.add_argument('--output-dir', required=True, type=Path)
    p.add_argument('--dist', type=int, default=0)
    p.add_argument('--N_contacts_min', type=int, default=1000)
    p.add_argument('--gene_len_min', type=int, default=200)
    p.add_argument('--fdr_threshold', type=float, default=0.05)
    p.add_argument('--protein-coding-biotype', type=str, default=None)
    a = p.parse_args()

    if a.dist < 0:
        p.error('--dist must be >= 0')

    if a.N_contacts_min < 0:
        p.error('--N_contacts_min must be >= 0')

    if a.gene_len_min < 0:
        p.error('--gene_len_min must be >= 0')

    if not 0 < a.fdr_threshold <= 1:
        p.error('--fdr_threshold must be in (0, 1]')
    
    if a.chP_mode == 'counts' and a.contact_counts is None:
        p.error('counts requires --contact-counts')
    
    if a.chP_mode == 'unique_dna' and a.contacts_uu is None:
        p.error('unique_dna requires --contacts-uu')
    
    if a.chP_mode == 'all_dna' and (a.dist > 0 or a.contact_counts is None):
        if a.contacts_uu is None or a.contacts_um is None:
            p.error('all_dna requires --contacts-uu and --contacts-um unless dist=0 and --contact-counts is provided')
    
    for path in (a.rnaseq_counts, a.contact_counts, a.contacts_uu, a.contacts_um):
        if path is not None and not path.is_file():
            p.error(f'file does not exist: {path}')
    
    return a


def main():
    """Select contact inputs, calculate chromatin potential, and write results."""
    a = arguments()
    genes = load_counts(a.rnaseq_counts)
    using_ready = a.chP_mode == 'counts' or (a.chP_mode == 'all_dna' and a.dist == 0 and a.contact_counts)
    stats, values = Counter(), Counter()
    if using_ready:
        ready = load_counts(a.contact_counts)
        extra = ready.keys() - genes.keys()
        if extra:
            raise ValueError(f'{a.contact_counts}: genes absent from RNA-seq: {list(sorted(extra))[:5]}')
        values.update({key: row['N_counts'] for key, row in ready.items()})
        print(f'Using prepared counts: {a.contact_counts} ({len(ready)} genes)')
    else:
        count_contacts(a.contacts_uu, 'UU', a.dist, genes, values, stats)
        if a.chP_mode == 'all_dna':
            count_contacts(a.contacts_um, 'UM', a.dist, genes, values, stats)
        a.output_dir.mkdir(parents=True, exist_ok=True)
        result = a.output_dir / f'contact_counts_{a.chP_mode}_dist_{a.dist}.tsv'
        write_counts(result, genes, values)
        with (a.output_dir / f'processing_summary_{a.chP_mode}_dist_{a.dist}.tsv').open('w') as stream:
            stream.write('metric\tvalue\n')
            for key in sorted(stats):
                stream.write(f'{key}\t{stats[key]}\n')
        print(f'Wrote {result}')
        print(f'Counted {sum(values.values())} contacts; {sum(v > 0 for v in values.values())} genes with contacts')
        if stats:
            print('Processing summary: ' + ', '.join(f'{key}={stats[key]}' for key in sorted(stats)))

    a.output_dir.mkdir(parents=True, exist_ok=True)

    name = ('chP.tsv' if a.chP_mode == 'counts'
            else f'chP_{a.chP_mode}_dist_{a.dist}.tsv')

    rows = calculate_chp(
        genes, values, a.N_contacts_min, a.gene_len_min,
        a.fdr_threshold, a.output_dir / name
    )

    plot_name = a.output_dir / name.replace('.tsv', '.png')

    if rows:
        joinplot_chP_Ncontacts(rows, plot_name, a)
    elif plot_name.exists():
        plot_name.unlink()


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as exc:
        raise SystemExit(f'Error: {exc}') from exc
