# 0. Import packages
import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

def collapse_strand(ht):
    """
    Return the deduplicated context by collapsing DNA strands.
    
    Function returns the reverse complement if the reference allele is either 'G' or 'T'.
    
    The reverse_complement_bases function has been made obsolete and should be replaced by `hl.reverse_complement`.
    
    :param ht: Input Table.
    :return: Table with deduplicated context annotation (ref, alt, context, was_flipped).
    """
    collapse_expr = {
        'ref': hl.if_else(((ht.ref == 'G') | (ht.ref == 'T')),
                    hl.reverse_complement(ht.ref), ht.ref),
        'alt': hl.if_else(((ht.ref == 'G') | (ht.ref == 'T')),
                    hl.reverse_complement(ht.alt), ht.alt),
        'context': hl.if_else(((ht.ref == 'G') | (ht.ref == 'T')),
                        hl.reverse_complement(ht.context), ht.context),
        'was_flipped': (ht.ref == 'G') | (ht.ref == 'T')
    }
    return ht.annotate(**collapse_expr) if isinstance(ht, hl.Table) else ht.annotate_rows(**collapse_expr)

def collapse_lof_call(ht):
    ht = ht.annotate(lof_csq_collapsed = hl.if_else(
        hl.is_defined(ht.vep.worst_csq_by_gene_canonical.lof), 
        ht.vep.worst_csq_by_gene_canonical.lof, 
        ht.vep.worst_csq_by_gene_canonical.most_severe_consequence))
    return ht

def main():

    # 2. Import data
    # Import gnomaAD v.3.1.2
    ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')

    # Initial filtering
    ht = gnomad.utils.vep.process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.filter((hl.len(ht.filters) == 0) & (ht.vep.worst_csq_by_gene_canonical.biotype == 'protein_coding') & (ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG')))

    # Import mutation rates from gnomAD paper
    ht_mu = hl.import_table('gs://janucik-dataproc-stage/01_maps/supplementary_dataset_10_mutation_rates.tsv.gz',
                    delimiter='\t', impute=True, force_bgz=True)

    # Convert mu_snp into string
    ht_mu = ht_mu.key_by('context', 'ref', 'alt', 'methylation_level') # has to have a key in order to join using foreign key

    # Import context table from gnomad (https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vep.html?highlight=context#gnomad.utils.vep.get_vep_context)
    context_table = gnomad.utils.vep.get_vep_context("GRCh38").ht()

    context_table_parsed = context_table.select(context=context_table.context[2:5], methyl_mean=hl.float64(context_table.methyl_mean))

    context_table_parsed = context_table_parsed.annotate(methylation_level = hl.case().when(
                hl.is_missing(context_table_parsed.methyl_mean) == 1, 0
            ).when(
                context_table_parsed.methyl_mean > 0.6, 2
            ).when(
                context_table_parsed.methyl_mean > 0.2, 1
            ).default(0))

    # 3. Add context field to main data
    # Split alleles field to ref and alt allele and collapse strands
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1], **context_table_parsed[ht.key])
    ht = collapse_strand(ht) # Collapse strands as some mutation contexts may not be present in mutation table
    # (SNP can be represented as both strands so collapsing makes the schema uniformal) 

    # 3.1 Implement loftee and filter to only synonymous or missense
    ht = collapse_lof_call(ht)
    ht = ht.filter((ht.lof_csq_collapsed == 'HC') | (ht.lof_csq_collapsed == 'LC') | (ht.lof_csq_collapsed == 'OS') | (ht.lof_csq_collapsed == 'missense_variant') | (ht.lof_csq_collapsed == 'synonymous_variant'))
    ht = ht.checkpoint('gs://janucik-dataproc-stage/01_maps/02b_maps_per_variant_array-main_17_Nov_22_v3/maps_per_variant_array_main.ht') # Save the data for further computation

if __name__ == '__main__':
    main()