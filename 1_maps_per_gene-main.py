# 0. Import packages
import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

# 1. Define helper functions
def train_on_synonymous(ht, ht_mu):
    """
    This function prepares the data for the regression model to predict the expected singleton/total variants ratio.
    """
    # Filter only synonymous variants
    ht_syn = ht.filter(ht.lof_csq_collapsed == "synonymous_variant")
    # Calculate number of variants in each tri-nucleotide context in synonymous variants
    ht_syn_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.methylation_level).aggregate(N_variants = hl.agg.count(), N_singletons = hl.agg.sum(ht_syn.info.singleton)))
    # Merge the N variants and N singletons tables
    ht_syn_ps = ht_syn_variants.annotate(ps = ht_syn_variants.N_singletons/ht_syn_variants.N_variants)
    # Add mutation rates and filter if mu_snp missing after join
    ht_syn_ps = ht_syn_ps.annotate(**ht_mu[ht_syn_ps.context, ht_syn_ps.ref, ht_syn_ps.alt, ht_syn_ps.methylation_level])
    ht_syn_ps = ht_syn_ps.filter(hl.is_defined(ht_syn_ps.mu_snp))
    return ht_syn_ps  

def regress_per_context(ht, ht_syn_lm, ht_mu):
    """
    This function calculates MAPS with Standard Error of the Mean for a given lof_csq_collapsed per gene.
    """
    # Add gene_id info field
    ht_reg_table = ht.annotate(gene_ids = ht.vep.worst_csq_by_gene_canonical.gene_id)
    # Count number of variants and singletons
    ht_reg_table_variants = (ht_reg_table.group_by(ht_reg_table.gene_ids, ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.methylation_level, ht_reg_table.lof_csq_collapsed).aggregate(N_variants = hl.agg.count(), N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))
    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_variants.annotate(ps = ht_reg_table_variants.N_singletons/ht_reg_table_variants.N_variants)
    # Add mu_snp info and filter to non-missing contexts
    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_mu[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    ht_reg_table_ps = ht_reg_table_ps.filter(hl.is_defined(ht_reg_table_ps.mu_snp))
    # Get expected number of singletons by applying the model factors
    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(expected_singletons=(ht_reg_table_ps.mu_snp * ht_syn_lm[1] + ht_syn_lm[0]) * ht_reg_table_ps.N_variants)
    # To aggregate just sum for the context and gene
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("gene_ids", "lof_csq_collapsed")
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
    """
    Collapse lof calls over consequence if present. 
    """
    ht = ht.annotate(lof_csq_collapsed = hl.if_else(
        hl.is_defined(ht.vep.worst_csq_by_gene_canonical.lof), 
        ht.vep.worst_csq_by_gene_canonical.lof, 
        ht.vep.worst_csq_by_gene_canonical.most_severe_consequence))
    return ht

def main():
    # 2. Import data
    # 2.1. Import gnomad v3 genomes
    ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')
    # Initial filtering
    ht = gnomad.utils.vep.process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.filter((hl.len(ht.filters) == 0) & (ht.vep.worst_csq_by_gene_canonical.biotype == 'protein_coding') & (ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG')))

    # 2.2 Import mutation rates from gnomAD paper
    ht_mu = hl.import_table('gs://janucik-dataproc-stage/01_maps/supplementary_dataset_10_mutation_rates.tsv.gz',
                    delimiter='\t', impute=True, force_bgz=True)
    ht_mu = ht_mu.key_by('context', 'ref', 'alt', 'methylation_level') # has to have a key in order to join using foreign key

    # 2.3 Import context table from gnomad
    context_table = gnomad.utils.vep.get_vep_context("GRCh38").ht()
    # Parse according to downstream needs
    context_table_parsed = context_table.select(context=context_table.context[2:5], methyl_mean=hl.float64(context_table.methyl_mean))
    # Recode methylation level into bins
    context_table_parsed = context_table_parsed.annotate(methylation_level = hl.case().when(
                hl.is_missing(context_table_parsed.methyl_mean) == 1, 0
            ).when(
                context_table_parsed.methyl_mean > 0.6, 2
            ).when(
                context_table_parsed.methyl_mean > 0.2, 1
            ).default(0))

    # 3. Annotate with context and methylation level, then export for computation
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1], **context_table_parsed[ht.key])
    # Collapse strands as some mutation contexts may not be present in mutation table
    ht = collapse_strand(ht) #  SNP can be represented as both strands so collapsing makes the schema uniformal
    # Implement loftee and filter to only synonymous or missense
    ht = collapse_lof_call(ht)
    ht = ht.filter((ht.lof_csq_collapsed == 'HC') | (ht.lof_csq_collapsed == 'LC') | (ht.lof_csq_collapsed == 'OS') | (ht.lof_csq_collapsed == 'missense_variant') | (ht.lof_csq_collapsed == 'synonymous_variant'))
    ht = ht.checkpoint('gs://janucik-dataproc-stage/01_maps/maps_per_gene_08_Nov_22_v1/maps_per_gene_main.ht') # Save the data for further computation
    
    # 4. Train linear model on synonymous variants for mutational rate correction
    print("SYNONYMOUS VARIANTS")
    ht_syn_ps = train_on_synonymous(ht, ht_mu)
    # Perform regression
    print("REGRESSING...")
    ht_syn_lm = ht_syn_ps.aggregate(hl.agg.linreg(ht_syn_ps.ps, [1, ht_syn_ps.mu_snp], weight=ht_syn_ps.N_variants).beta)
    # Show intercept and beta
    print("REGRESSION PARAMETERS")
    print(ht_syn_lm)

    # 5. Predict expected number of variants for each context
    maps_table = regress_per_context(ht, ht_syn_lm, ht_mu)

    # 6. Export tables for plotting and exploring
    maps_table.export('gs://janucik-dataproc-stage/01_maps/maps_per_gene_08_Nov_22_v1/02a_f_maps_table.csv', delimiter=',')
    ht_syn_ps.export('gs://janucik-dataproc-stage/01_maps/maps_per_gene_08_Nov_22_v1/ht_syn_ps.csv', delimiter=',')
    
if __name__ == '__main__':
    main()