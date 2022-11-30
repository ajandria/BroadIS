# 0. Import packages
import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

# 1. Define helper functions
def train_on_synonymous(ht):
    """
    Trains the regression model to predict the expected singleton/total variants ratio.
    """
    # Filter only synonymous variants.
    ht_syn = ht.filter(ht.lof_csq_collapsed == "synonymous_variant")

    # Get number of singletons and number of variants in total for each variant type. 
    ht_syn = (ht_syn.group_by('locus', 'alleles', 'context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(
        singleton = ((hl.agg.array_agg(lambda element: hl.agg.sum(element == 1), ht_syn.freq.AC))), 
        variant = ((hl.agg.array_agg(lambda element: hl.agg.sum(element > 0), ht_syn.freq.AC)))))

    # Calculate number of variants in each tri-nucleotide context and methylation level in synonymous variants.
    ht_syn_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.methylation_level).aggregate(
        N_variants = hl.agg.array_sum(ht_syn.variant),
        N_singletons = hl.agg.array_sum(ht_syn.singleton)))

    # Compute singleton to total number of variants ratio.
    ht_syn_ps = ht_syn_variants.annotate(
        ps = ht_syn_variants.N_singletons/ht_syn_variants.N_variants)
    
    return ht_syn_ps  

## Define regression function
def regress_per_context(ht, ht_syn_lm, ht_mu):
    """
    Calculates MAPS with Standard Error of the Mean for a given lof_csq_collapsed.
    """
    # Was based on *_reg_table.
    ht_reg_table = ht

    # Annotate gene ids to the main ht
    ht_reg_table = ht_reg_table.annotate(gene_ids = ht_reg_tableht.vep.worst_csq_by_gene_canonical.gene_id)

    # Get number of singletons and number of variants in total for each variant type. 
    ht_reg_table_variants = (ht_reg_table.group_by('locus', 'alleles', 'gene_ids', 'context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(
        singletons = ((hl.agg.array_agg(lambda x: hl.agg.sum(x[0] == 1), ht_reg_table.freq))), 
        variants = ((hl.agg.array_agg(lambda x: hl.agg.sum(x[0] > 0), ht_reg_table.freq)))))

    # Collapse the singleton and variant frequencies into per tri-nucletide context and methylation level.
    ht_reg_table_variants = (ht_reg_table_variants.group_by('gene_ids', 'context', 'ref', 'alt', 'methylation_level', 'lof_csq_collapsed').aggregate(
        N_variants = hl.agg.array_sum(ht_reg_table_variants.variants), 
        N_singletons = hl.agg.array_sum(ht_reg_table_variants.singletons)))

    # Compute singleton to total number of variants ratio (this is actually not needed as it gets recomputed later, but can be usefull for debugging).
    ht_reg_table_ps = ht_reg_table_variants.annotate(
        ps = ht_reg_table_variants.N_singletons/ht_reg_table_variants.N_variants)
    
    # Annotate mutational rates of each tri-nucleotide context.
    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_mu[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    
    # Filter to contexts with only defined mutational rates.
    ht_reg_table_ps = ht_reg_table_ps.filter(hl.is_defined(ht_reg_table_ps.mu_snp))

    # Remove keys prior to annotation by foreign key.
    ht_syn_lm = ht_syn_lm.key_by()

    # Pick a subset for easier debugging.
    ht_syn_lm_sub = ht_syn_lm.select('context', 'ref', 'alt', 'methylation_level', 'intercept', 'beta1')

    # Add keys to the table with N_singletons, N_variants and PS prior to annotation with synonymous model factors. 
    ht_syn_lm_sub = ht_syn_lm_sub.key_by('context', 'ref', 'alt', 'methylation_level')

    # Add regression model factors from synonymous variants into the main hail table. 
    ht_reg_table_ps = ht_reg_table_ps.annotate(**ht_syn_lm_sub[ht_reg_table_ps.context, ht_reg_table_ps.ref, ht_reg_table_ps.alt, ht_reg_table_ps.methylation_level])
    
    # This adds the expected singleton info adjusting for the synonymous variants stability accross the whole genome.
    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(
        expected_singletons = (ht_reg_table_ps.mu_snp * ht_reg_table_ps.beta1 + ht_reg_table_ps.intercept) * ht_reg_table_ps.N_variants)

    # This aggregates the N_singletons, expected_singletons and N_variants into the SNVs class. 
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("gene_ids", "lof_csq_collapsed")
                .aggregate(N_singletons=hl.agg.array_sum(ht_reg_table_ps_lm_cons.N_singletons),
                            expected_singletons=hl.agg.array_sum(ht_reg_table_ps_lm_cons.expected_singletons),
                            N_variants=hl.agg.array_sum(ht_reg_table_ps_lm_cons.N_variants)))

    # Recompute singleton to total number of variants proportion per SNVs class and MAPS scores using this
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg.annotate(
        ps_agg=ht_reg_table_ps_lm_cons_agg.N_singletons / ht_reg_table_ps_lm_cons_agg.N_variants,
        maps=(ht_reg_table_ps_lm_cons_agg.N_singletons - ht_reg_table_ps_lm_cons_agg.expected_singletons) / ht_reg_table_ps_lm_cons_agg.N_variants)

    # Add standard error of the mean (sem) to MAPS
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg_MAPS.annotate(
        maps_sem=(ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg * (1 - ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg) / ht_reg_table_ps_lm_cons_agg_MAPS.N_variants) ** 0.5)

    return ht_reg_table_ps_lm_cons_agg_MAPS

def collapse_strand(ht):
    """
    Return the deduplicated context by collapsing DNA strands.
    Function returns the reverse complement if the reference allele is either 'G' or 'T'.
    The reverse_complement_bases function has been made obsolete and should be replaced by `hl.reverse_complement`.
    """
    # Combined if_else() statements to form a switch case situation.
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
    Overwrites consequence with LOF calls.
    """
    # Collapse LOF call on top of SNVs classes if exist. 
    ht = ht.annotate(lof_csq_collapsed = hl.if_else(
        hl.is_defined(ht.vep.worst_csq_by_gene_canonical.lof), 
        ht.vep.worst_csq_by_gene_canonical.lof, 
        ht.vep.worst_csq_by_gene_canonical.most_severe_consequence))

    return ht

def main():
    # 2. Import data
    # Import gnomaAD v.3.1.2 Genomes
    ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')

    # Get consequences by canonical transcript for genes
    ht = gnomad.utils.vep.process_consequences(ht)
    # Include overlapping genes
    ht = ht.explode(ht.vep.worst_csq_by_gene_canonical)
    # Keep SNPs passing filters, related to protein_coding ENSG genes. 
    ht = ht.filter((hl.len(ht.filters) == 0) & (ht.vep.worst_csq_by_gene_canonical.biotype == 'protein_coding') & (ht.vep.worst_csq_by_gene_canonical.gene_id.startswith('ENSG')))
    
    # Import mutation rates from gnomAD paper.
    ht_mu = hl.import_table('gs://janucik-dataproc-stage/01_maps/supplementary_dataset_10_mutation_rates.tsv.gz',
                    delimiter='\t', impute=True, force_bgz=True)

    # Add keys prior to using foreign join.
    ht_mu = ht_mu.key_by('context', 'ref', 'alt', 'methylation_level') 

    # Import context table from gnomad (https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vep.html?highlight=context#gnomad.utils.vep.get_vep_context).
    context_table = gnomad.utils.vep.get_vep_context("GRCh38").ht()

    # Parse contexts and methylation info.
    context_table_parsed = context_table.select(
        context = context_table.context[2:5], 
        methyl_mean = hl.float64(context_table.methyl_mean))

    # Recode methylation info into binned values.
    context_table_parsed = context_table_parsed.annotate(methylation_level = hl.case().when(
                hl.is_missing(context_table_parsed.methyl_mean) == 1, 0
            ).when(
                context_table_parsed.methyl_mean > 0.6, 2
            ).when(
                context_table_parsed.methyl_mean > 0.2, 1
            ).default(0))
    
    # 3. Merge context info into the main ht.
    # Split alleles into ref, alt and add contexts by matching locus and alleles fields.
    ht = ht.annotate(
        ref = ht.alleles[0], 
        alt = ht.alleles[1], 
        **context_table_parsed[ht.key])

    # Collapse strands as some mutation contexts may not be present in mutation table.
    # (SNP can be represented as both strands so collapsing makes the schema uniformal).
    ht = collapse_strand(ht) 
    
    # 4. Implement loftee and filter to only synonymous, missense or LOF calls.
    ht = collapse_lof_call(ht)
    ht = ht.filter((ht.lof_csq_collapsed == 'HC') | (ht.lof_csq_collapsed == 'LC') | (ht.lof_csq_collapsed == 'OS') | (ht.lof_csq_collapsed == 'missense_variant') | (ht.lof_csq_collapsed == 'synonymous_variant'))
    
    # Checkpoint or read here if debug needed.
    #ht = hl.read_table('gs://janucik-dataproc-stage/01_maps/1-array-rerun_maps_per_variant-separate-models-main_2022_Nov_24_v1/main_ht_upstream.ht')
    ht = ht.checkpoint('gs://janucik-dataproc-stage/01_maps/v1_30_Nov_2022_MAPS_per_gene_using_arrays/main_ht_upstream.ht')
    
    # 5. Train linear model on synonymous variants for mutational class correction.
    print("SYNONYMOUS VARIANTS")
    ht_synonymous_count = ht.filter(ht.lof_csq_collapsed == "synonymous_variant")

    ht_syn_ps = train_on_synonymous(ht)
    # Annotate mutational rate just before computation to decrease the risk of malformed values.
    ht_syn_ps = ht_syn_ps.annotate(**ht_mu[ht_syn_ps.context, ht_syn_ps.ref, ht_syn_ps.alt, ht_syn_ps.methylation_level])
    # Filter to non missing mutation rates.
    ht_syn_ps = ht_syn_ps.filter(hl.is_defined(ht_syn_ps.mu_snp))

    # Perform regression.
    # Zip PS and total number of variants as required by the model.
    ht_syn_ps = ht_syn_ps.annotate(ps_nVariants = hl.zip(ht_syn_ps.ps, ht_syn_ps.N_variants))
    # Compute intercept.
    ht_syn_lm = ht_syn_ps.annotate(intercept = ht_syn_ps.aggregate(hl.agg.array_agg(
        lambda element: hl.agg.linreg(element[0], [1, ht_syn_ps.mu_snp], weight=element[1]).beta[0], ht_syn_ps.ps_nVariants)))
    # Compute slope for mutational rate.
    ht_syn_lm = ht_syn_lm.annotate(beta1 = ht_syn_lm.aggregate(hl.agg.array_agg(
        lambda element: hl.agg.linreg(element[0], [1, ht_syn_lm.mu_snp], weight=element[1]).beta[1], ht_syn_lm.ps_nVariants)))
    
    # 6. Predict number of singletons based on synonymous variation accross the whole genome and compute MAPS. 
    maps_table = regress_per_context(ht, ht_syn_lm, ht_mu)

    # 7. Export the data.
    # Final hail table.
    maps_table.write('gs://janucik-dataproc-stage/01_maps/v1_30_Nov_2022_MAPS_per_gene_using_arrays/maps_downsampling_multi_model.ht')
    # Export final hail table as csv for debugging/plotting.
    maps_table.export('gs://janucik-dataproc-stage/01_maps/v1_30_Nov_2022_MAPS_per_gene_using_arrays/maps_downsampling_multi_model.csv', delimiter=',')
    # Export MAPS table as csv for debugging/plotting.
    ht_syn_ps.export('gs://janucik-dataproc-stage/01_maps/v1_30_Nov_2022_MAPS_per_gene_using_arrays/maps_downsampling_multi_model_syn_ps.csv', delimiter=',')

# Run the pipeline from this script.     
if __name__ == '__main__':
    main()