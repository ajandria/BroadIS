# 0. Import packages
import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

# 1. Define helper functions
## Define data wrangling function for synonymous variants
def train_on_synonymous(ht):
    """
    This function trains the regression model to predict the expected singleton/total variants ratio.
    
    Requirements (things to change as for best parcitces)
    1. Already annotated mutation rates for each tri-nucleotide context in the main hail table
    """
    # Filter only synonymous variants
    ht_syn = ht.filter(ht.vep.most_severe_consequence == "synonymous_variant")

    # Calculate number of variants in each tri-nucleotide context in synonymous variants
    ht_syn_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.methylation_level).aggregate(N_variants = hl.agg.count(), N_singletons = hl.agg.sum(ht_syn.info.singleton)))

    # Merge the N variants and N singletons tables
    ht_syn_ps = ht_syn_variants.annotate(ps = ht_syn_variants.N_singletons/ht_syn_variants.N_variants)
    
    return ht_syn_ps  

## Define regression function
def regress_per_context(ht, ht_syn_lm, ht_mu):
    """
    Calculates MAPS with Standard Error of the Mean for a given consequence.
    
    Requires input in form of Hail table: (so the things that should be changes for easier accessibility):
        1) Mutation rates are to be annotated for each tri-nucleotide context. 
    """
    ht_reg_table = ht.annotate(consequence = ht.most_severe_consequence)

    # Count number of variants and singletons
    ht_reg_table_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.methylation_level, ht_reg_table.consequence).aggregate(N_variants = hl.agg.count(), N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_variants.annotate(ps = ht_reg_table_variants.N_singletons/ht_reg_table_variants.N_variants)
    
    # Get expected number of singletons by applying the model factors
    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_mu[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    
    ht_reg_table_ps = ht_reg_table_ps.filter(hl.is_defined(ht_reg_table_ps.mu_snp))

    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(expected_singletons=(ht_reg_table_ps.mu_snp * ht_syn_lm[1] + ht_syn_lm[0]) * ht_reg_table_ps.N_variants)

    # To aggregate just sum for the context
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("consequence")
                .aggregate(N_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.N_singletons),
                            expected_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.expected_singletons),
                            N_variants=hl.agg.sum(ht_reg_table_ps_lm_cons.N_variants)))

    # Calculate MAPS and aggregated proportions 
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg.annotate(ps_agg=ht_reg_table_ps_lm_cons_agg.N_singletons / ht_reg_table_ps_lm_cons_agg.N_variants,
        maps=(ht_reg_table_ps_lm_cons_agg.N_singletons - ht_reg_table_ps_lm_cons_agg.expected_singletons) / ht_reg_table_ps_lm_cons_agg.N_variants)

    # Add MAPS standard error of the mean (sem)
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg_MAPS.annotate(maps_sem=(ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg * (1 - ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg) / ht_reg_table_ps_lm_cons_agg_MAPS.N_variants) ** 0.5)
    return ht_reg_table_ps_lm_cons_agg_MAPS

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
    ht = ht.annotate(most_severe_consequence = hl.if_else(
        hl.is_defined(ht.vep.transcript_consequences.lof[0]), 
        ht.vep.transcript_consequences.lof[0], 
        ht.vep.most_severe_consequence))
    return ht

def main():

    # 2. Import data
    # Import gnomaAD v.3.1.2
    ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')

    # Initial filtering
    ht = gnomad.utils.vep.process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht.describe()

if __name__ == '__main__':
    main()