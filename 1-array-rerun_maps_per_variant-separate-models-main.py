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
    ht_syn = ht.filter(ht.lof_csq_collapsed == "synonymous_variant")

    # ht_syn = ht_syn.annotate(singleton = ht_syn.freq.AC[0] == 1, variant = ht_syn.freq.AC[0] > 0)
    # ht_syn = ht_syn.annotate(is_singleton = ht_syn.aggregate(hl.agg.array_agg(lambda element: hl.agg.sum(element == 1), ht_syn.freq.AC))).show()

    ht_syn = (ht_syn.group_by('locus', 'alleles', 'context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(singleton = ((hl.agg.array_agg(lambda element: hl.agg.sum(element == 1), ht_syn.freq.AC))), variant = ((hl.agg.array_agg(lambda element: hl.agg.sum(element > 0), ht_syn.freq.AC)))))

                    

    # Calculate number of variants in each tri-nucleotide context in synonymous variants
    ht_syn_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.methylation_level).aggregate(N_variants = hl.agg.array_sum(ht_syn.variant), N_singletons = hl.agg.array_sum(ht_syn.singleton)))

    # Merge the N variants and N singletons tables
    ht_syn_ps = ht_syn_variants.annotate(ps = ht_syn_variants.N_singletons/ht_syn_variants.N_variants)
    
    return ht_syn_ps  

## Define regression function
def regress_per_context(ht, ht_syn_lm, ht_mu):
    """
    Calculates MAPS with Standard Error of the Mean for a given lof_csq_collapsed.
    
    Requires input in form of Hail table: (so the things that should be changes for easier accessibility):
        1) Mutation rates are to be annotated for each tri-nucleotide context. 
    """
    ht_reg_table = ht

    # Count number of variants and singletons
    #ht_reg_table_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.methylation_level, ht_reg_table.lof_csq_collapsed).aggregate(N_variants = hl.agg.array_sum(ht_reg_table.freq.AC), N_singletons = ht_reg_table.aggregate(hl.agg.array_agg(lambda AC: hl.agg.count_where(AC == 1), ht_reg_table.freq.AC))))

    ht_reg_table_variants = (ht_reg_table.group_by('locus', 'alleles', 'context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(singletons = ((hl.agg.array_agg(lambda x: hl.agg.sum(x[0] == 1), ht_reg_table.freq))), variants = ((hl.agg.array_agg(lambda x: hl.agg.sum(x[0] > 0), ht_reg_table.freq)))))

    ht_reg_table_variants = (ht_reg_table_variants.group_by('context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(N_variants = hl.agg.array_sum(ht_reg_table_variants.variants), N_singletons = hl.agg.array_sum(ht_reg_table_variants.singletons)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_variants.annotate(ps = ht_reg_table_variants.N_singletons/ht_reg_table_variants.N_variants)
    
    # Get expected number of singletons by applying the model factors
    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_mu[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    
    ht_reg_table_ps = ht_reg_table_ps.filter(hl.is_defined(ht_reg_table_ps.mu_snp))

         
    ht_syn_lm = ht_syn_lm.key_by()

    ht_syn_lm_sub = ht_syn_lm.select('context', 'ref', 'alt', 'methylation_level', 'betas')

    ht_syn_lm_sub = ht_syn_lm_sub.key_by('context', 'ref', 'alt', 'methylation_level')

    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_syn_lm_sub[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    
    ht_reg_table_ps.show()
    
    ht_reg_table_ps = ht_reg_table_ps.annotate(intercept = ht_reg_table_ps.betas[0], beta1 = ht_reg_table_ps.betas[1])

    ht_reg_table_ps = ht_reg_table_ps.annotate(betas_variants = hl.zip(ht_reg_table_ps.intercept, ht_reg_table_ps.beta1, ht_reg_table_ps.N_variants))

    ht_reg_table_ps.show()

    # This adds the expected singleton info
    ht_reg_table_ps_lm_cons = (ht_reg_table_ps.group_by('lof_csq_collapsed', 'context', 'ref', 'alt', 'methylation_level').aggregate(expected_singletons = hl.agg.array_agg(lambda element: hl.agg.sum((ht_reg_table_ps.mu_snp * element[1] + element[0]) * element[2]), ht_reg_table_ps.betas_variants)))

    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(**ht_reg_table_ps_lm_cons[ht_reg_table_ps.lof_csq_collapsed, ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])

    # To aggregate just sum for the context
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("lof_csq_collapsed")
                .aggregate(N_singletons=hl.agg.array_sum(ht_reg_table_ps_lm_cons.N_singletons),
                            expected_singletons=hl.agg.array_sum(ht_reg_table_ps_lm_cons.expected_singletons),
                            N_variants=hl.agg.array_sum(ht_reg_table_ps_lm_cons.N_variants)))

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
    ht = ht.annotate(lof_csq_collapsed = hl.if_else(
        hl.is_defined(ht.vep.worst_csq_by_gene_canonical.lof), 
        ht.vep.worst_csq_by_gene_canonical.lof, 
        ht.vep.worst_csq_by_gene_canonical.most_severe_consequence))
    return ht

def main():
    """
    # 2. Import data
    # Import gnomaAD v.3.1.2
    ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')

    # Initial filtering
    ht = gnomad.utils.vep.process_consequences(ht)
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    ht = ht.filter((hl.len(ht.filters) == 0) & (ht.vep.worst_csq_by_gene_canonical.biotype == 'protein_coding') & (ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG')))
    """
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
    """ 
    # 3. Add context field to main data
    # Split alleles field to ref and alt allele and collapse strands
    ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1], **context_table_parsed[ht.key])
    ht = collapse_strand(ht) # Collapse strands as some mutation contexts may not be present in mutation table
    # (SNP can be represented as both strands so collapsing makes the schema uniformal) 

    # 3.1 Implement loftee and filter to only synonymous or missense
    ht = collapse_lof_call(ht)
    ht = ht.filter((ht.lof_csq_collapsed == 'HC') | (ht.lof_csq_collapsed == 'LC') | (ht.lof_csq_collapsed == 'OS') | (ht.lof_csq_collapsed == 'missense_variant') | (ht.lof_csq_collapsed == 'synonymous_variant'))
    """
    # Checkpoint or read here
    ht = hl.read_table('gs://janucik-dataproc-stage/01_maps/02b_maps_per_variant_array-main_17_Nov_22_v2/maps_per_variant_array_main.ht')
    

    # 4. Train linear model on synonymous variants for mutational class correction
    print("SYNONYMOUS VARIANTS")
    ht_synonymous_count = ht.filter(ht.lof_csq_collapsed == "synonymous_variant")

    ht_syn_ps = train_on_synonymous(ht)
    ht_syn_ps = ht_syn_ps.annotate(**ht_mu[ht_syn_ps.context, ht_syn_ps.ref, ht_syn_ps.alt, ht_syn_ps.methylation_level])
    ht_syn_ps = ht_syn_ps.filter(hl.is_defined(ht_syn_ps.mu_snp))

    # Perform regression
    ht_syn_ps = ht_syn_ps.annotate(ps_nVariants = hl.zip(ht_syn_ps.ps, ht_syn_ps.N_variants))
    ht_syn_lm = ht_syn_ps.annotate(betas = ht_syn_ps.aggregate(hl.agg.array_agg(lambda element: hl.agg.linreg(element[0], [1, ht_syn_ps.mu_snp], weight=element[1]).beta, ht_syn_ps.ps_nVariants)))

    # Show intercept and beta
    print("REGRESSION PARAMETERS")
    #ht_syn_lm.show()

    
    # 5. Predict expected number of variants for each context
    maps_table = regress_per_context(ht, ht_syn_lm, ht_mu)

    # Add meta fields
    #meta_fields = ht.freq_meta.collect()
    #maps_table = maps_table.annotate_globals(freq_metas = meta_fields)

    maps_table.show()

    #maps_table = maps_table.annotate_globals(ht.freq_meta.collect())

    maps_table.write('gs://janucik-dataproc-stage/01_maps/1_array_rerun_maps_per_variant_main_21_Nov_23_v2/02a_f_maps_table.ht')

    # 6. Export tables for plotting and exploring
    maps_table.export('gs://janucik-dataproc-stage/01_maps/1_array_rerun_maps_per_variant_main_21_Nov_23_v2/02a_f_maps_table.csv', delimiter=',')
    ht_syn_ps.export('gs://janucik-dataproc-stage/01_maps/1_array_rerun_maps_per_variant_main_21_Nov_23_v2/ht_syn_ps.csv', delimiter=',')
    
if __name__ == '__main__':
    main()