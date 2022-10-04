## 0. Define helper functions

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
    ht_syn_N_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.mu_snp, ht_syn.methylation_level).aggregate(N_variants = hl.agg.count()))

    # Calculate number of singletons for each tri-nucleotide context in synonymous variants
    ht_syn_singletons = ht_syn.filter(ht_syn.info.singleton == 1)
    ht_syn_N_singletons = (ht_syn_singletons.group_by(ht_syn_singletons.context, ht_syn_singletons.ref, ht_syn_singletons.alt, ht_syn_singletons.mu_snp, ht_syn_singletons.methylation_level).aggregate(N_singletons = hl.agg.count()))

    # Merge the N variants and N singletons tables
    ht_syn_ps  = ht_syn_N_variants.join(ht_syn_N_singletons, how = 'outer') # outer as all will match and we want all
    ht_syn_ps = ht_syn_ps.annotate(ps = ht_syn_ps.N_singletons/ht_syn_ps.N_variants)
    
    return ht_syn_ps 

## Define regression function
def regress_per_context(ht, ht_syn_lm):
    """
    Calculates MAPS with Standard Error of the Mean for a given consequence.
    
    Requires input in form of Hail table: (so the things that should be changes for easier accessibility):
        1) Mutation rates are to be annotated for each tri-nucleotide context. 
    """
    ht_reg_table = ht.annotate(consequence = ht.vep.most_severe_consequence)

    # Count number of variants and singletons
    ht_reg_table_N_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.methylation_level, ht_reg_table.consequence).aggregate(N_variants = hl.agg.count()))
    ht_reg_table_N_singletons = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.methylation_level, ht_reg_table.consequence).aggregate(N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_N_variants.annotate(**ht_reg_table_N_singletons[ht_reg_table_N_variants.context, ht_reg_table_N_variants.ref, ht_reg_table_N_variants.alt, ht_reg_table_N_variants.mu_snp, ht_reg_table_N_variants.methylation_level, ht_reg_table_N_variants.consequence])
    ht_reg_table_ps = ht_reg_table_ps.annotate(ps = ht_reg_table_ps.N_singletons/ht_reg_table_ps.N_variants)
    
    # Get expected number of singletons by applying the model factors
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
    
## 1. Import packages
import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

## 2. Import data
# Import gnomaAD v.3.1.2
ht = hl.read_table('gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht')
ht = ht.head(100000) # Subset the data

# Import mutation rates from gnomAD paper
ht_mu = hl.import_table('data/supplementary_dataset_10_mutation_rates.tsv.gz',
                delimiter='\t', impute=True, force_bgz=True)
ht_mu = ht_mu.key_by('context', 'ref', 'alt', 'methylation_level') # has to have a key in order to join using foreign key

# Import context table from gnomad (https://broadinstitute.github.io/gnomad_methods/api_reference/utils/vep.html?highlight=context#gnomad.utils.vep.get_vep_context)
context_table = gnomad.utils.vep.get_vep_context("GRCh38").ht()
context_table = context_table.filter(hl.is_defined(context_table.methyl_mean))

context_table_parsed = context_table.select(context_table.context, context_table.methyl_mean)
context_table_parsed = context_table_parsed.transmute(context = context_table_parsed.context[2:5])

# Bin methylation mean to methylation levels
context_table_parsed = context_table_parsed.annotate(methyl_mean = hl.float64(context_table_parsed.methyl_mean))
context_table_parsed = context_table_parsed.annotate(methylation_level = hl.if_else(context_table_parsed.methyl_mean <= 0.2,0,
                                                            hl.if_else(context_table_parsed.methyl_mean <= 0.8,1,2)))

### Show the data structure
ht.show(3)
# Table with methylation level and mutational rate in the trinucleotide context
ht_mu.show(3)
# This table contains already precalculated nucleotides -3/+3 from mutation site 
context_table_parsed.show(3)

## 3. Add context field to main data
# Before joining the tri-nucleotide context of mutation
ht.count()
# Join only matching rows from context to ht table.
#ht = ht.key_by('locus', 'alleles').join(context_table_parsed.key_by('locus', 'alleles'), how = 'left') # resorts the data making it slow
ht = ht.annotate(**context_table_parsed[ht.locus, ht.alleles])
# After
ht.count()

# Subset the data for local computation
ht = ht.filter(hl.is_defined(ht.methylation_level)) # makes sure there is methylation_level variables to be laters used in the grouping
ht = ht.head(10000) # Subset the data
ht.show(3)


## 4. Add mutation rates for added contexts
# Split alleles field to ref and alt allele
ht = ht.annotate(ref=ht.alleles[0], alt=ht.alleles[1])

# Add mutation rates according to the context, but also ref and alt allele for this context
#ht = ht.key_by("context", "ref", "alt").join(ht_mu.key_by("context", "ref", "alt"), how = 'left') # resorts the data making it slow
ht = ht.annotate(**ht_mu[ht.context, ht.ref, ht.alt, ht.methylation_level])

## 5. Train linear model on synonymous variants for mutational class correction
ht_syn_ps = train_on_synonymous(ht)
ht_syn_ps.show(3)

# Perform regression
ht_syn_lm = ht_syn_ps.aggregate(hl.agg.linreg(ht_syn_ps.ps, [1, ht_syn_ps.mu_snp], weight=ht_syn_ps.N_variants).beta)
# Show intercept and beta
ht_syn_lm

## 6. Predict expected number of variants for each context
### Function for regression eventually will be made starting here and put in `/utils/utils.py` script
maps_table = regress_per_context(ht, ht_syn_lm)
maps_table.show(maps_table.count())

### Visualise
ggplot(maps_table, aes(x=maps_table.consequence, y = maps_table.maps)) + geom_col(aes(fill=maps_table.consequence))


