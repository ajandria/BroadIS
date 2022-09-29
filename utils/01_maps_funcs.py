
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
    ht_syn_N_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.mu_snp).aggregate(N_variants = hl.agg.count()))

    # Calculate number of singletons for each tri-nucleotide context in synonymous variants
    ht_syn_singletons = ht_syn.filter(ht_syn.info.singleton == 1)
    ht_syn_N_singletons = (ht_syn_singletons.group_by(ht_syn_singletons.context, ht_syn_singletons.ref, ht_syn_singletons.alt, ht_syn_singletons.mu_snp).aggregate(N_singletons = hl.agg.count()))

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
    ht_reg_table_N_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence).aggregate(N_variants = hl.agg.count()))
    ht_reg_table_N_singletons = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence).aggregate(N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_N_variants.annotate(**ht_reg_table_N_singletons[ht_reg_table_N_variants.context, ht_reg_table_N_variants.ref, ht_reg_table_N_variants.alt, ht_reg_table_N_variants.mu_snp, ht_reg_table_N_variants.consequence])
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



## Define data wrangling function for synonymous variants
def train_on_synonymous_per_gene(ht):
    """
    This function trains the regression model to predict the expected singleton/total variants ratio.
    
    Requirements (things to change as for best parcitces)
    1. Already annotated mutation rates for each tri-nucleotide context in the main hail table
    """
    ht = ht.annotate(gene_id = ht.vep.transcript_consequences.gene_id)

    # Filter only synonymous variants
    ht_syn = ht.filter(ht.vep.most_severe_consequence == "synonymous_variant")

    # Calculate number of variants in each tri-nucleotide context in synonymous variants
    ht_syn_N_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.mu_snp, ht_syn.gene_id).aggregate(N_variants = hl.agg.count()))

    # Calculate number of singletons for each tri-nucleotide context in synonymous variants
    ht_syn_singletons = ht_syn.filter(ht_syn.info.singleton == 1)
    ht_syn_N_singletons = (ht_syn_singletons.group_by(ht_syn_singletons.context, ht_syn_singletons.ref, ht_syn_singletons.alt, ht_syn_singletons.mu_snp, ht_syn_singletons.gene_id).aggregate(N_singletons = hl.agg.count()))

    # Merge the N variants and N singletons tables
    ht_syn_ps  = ht_syn_N_variants.join(ht_syn_N_singletons, how = 'outer') # outer as all will match and we want all
    ht_syn_ps = ht_syn_ps.annotate(ps = ht_syn_ps.N_singletons/ht_syn_ps.N_variants)
    
    return ht_syn_ps 

## Define regression function
def regress_per_context_per_gene(ht, ht_syn_lm):
    """
    Calculates MAPS with Standard Error of the Mean for a given consequence.
    
    Requires input in form of Hail table: (so the things that should be changes for easier accessibility):
        1) Mutation rates are to be annotated for each tri-nucleotide context. 
    """
    ht = ht.annotate(gene_id = ht.vep.transcript_consequences.gene_id)
    ht_reg_table = ht.annotate(consequence = ht.vep.most_severe_consequence)

    # Count number of variants and singletons
    ht_reg_table_N_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence, ht_reg_table.gene_id).aggregate(N_variants = hl.agg.count()))
    ht_reg_table_N_singletons = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence, ht_reg_table.gene_id).aggregate(N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_N_variants.annotate(**ht_reg_table_N_singletons[ht_reg_table_N_variants.context, ht_reg_table_N_variants.ref, ht_reg_table_N_variants.alt, ht_reg_table_N_variants.mu_snp, ht_reg_table_N_variants.consequence, ht_reg_table_N_variants.gene_id])
    ht_reg_table_ps = ht_reg_table_ps.annotate(ps = ht_reg_table_ps.N_singletons/ht_reg_table_ps.N_variants)
    
    # Get expected number of singletons by applying the model factors
    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(expected_singletons=(ht_reg_table_ps.mu_snp * ht_syn_lm[1] + ht_syn_lm[0]) * ht_reg_table_ps.N_variants)

    # To aggregate just sum for the context
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("consequence", "gene_id")
                  .aggregate(N_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.N_singletons),
                             expected_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.expected_singletons),
                             N_variants=hl.agg.sum(ht_reg_table_ps_lm_cons.N_variants)))

    # Calculate MAPS and aggregated proportions 
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg.annotate(ps_agg=ht_reg_table_ps_lm_cons_agg.N_singletons / ht_reg_table_ps_lm_cons_agg.N_variants,
        maps=(ht_reg_table_ps_lm_cons_agg.N_singletons - ht_reg_table_ps_lm_cons_agg.expected_singletons) / ht_reg_table_ps_lm_cons_agg.N_variants)

    # Add MAPS standard error of the mean (sem)
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg_MAPS.annotate(maps_sem=(ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg * (1 - ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg) / ht_reg_table_ps_lm_cons_agg_MAPS.N_variants) ** 0.5)
    return ht_reg_table_ps_lm_cons_agg_MAPS

## Define data wrangling function for synonymous variants
def train_on_synonymous_per_gene_03b(ht):
    """
    This function trains the regression model to predict the expected singleton/total variants ratio.
    
    Requirements (things to change as for best parcitces)
    1. Already annotated mutation rates for each tri-nucleotide context in the main hail table
    """
    ht = ht.annotate(gene_id = ht.vep.worst_csq_for_variant.gene_id)

    # Filter only synonymous variants
    ht_syn = ht.filter(ht.vep.most_severe_consequence == "synonymous_variant")

    # Calculate number of variants in each tri-nucleotide context in synonymous variants
    ht_syn_N_variants = (ht_syn.group_by(ht_syn.context, ht_syn.ref, ht_syn.alt, ht_syn.mu_snp, ht_syn.gene_id).aggregate(N_variants = hl.agg.count()))

    # Calculate number of singletons for each tri-nucleotide context in synonymous variants
    ht_syn_singletons = ht_syn.filter(ht_syn.info.singleton == 1)
    ht_syn_N_singletons = (ht_syn_singletons.group_by(ht_syn_singletons.context, ht_syn_singletons.ref, ht_syn_singletons.alt, ht_syn_singletons.mu_snp, ht_syn_singletons.gene_id).aggregate(N_singletons = hl.agg.count()))

    # Merge the N variants and N singletons tables
    ht_syn_ps  = ht_syn_N_variants.join(ht_syn_N_singletons, how = 'outer') # outer as all will match and we want all
    ht_syn_ps = ht_syn_ps.annotate(ps = ht_syn_ps.N_singletons/ht_syn_ps.N_variants)
    
    return ht_syn_ps 

## Define regression function
def regress_per_context_per_gene_03b(ht, ht_syn_lm):
    """
    Calculates MAPS with Standard Error of the Mean for a given consequence.
    
    Requires input in form of Hail table: (so the things that should be changes for easier accessibility):
        1) Mutation rates are to be annotated for each tri-nucleotide context. 
    """
    ht = ht.annotate(gene_id = ht.vep.worst_csq_for_variant.gene_id)
    ht_reg_table = ht.annotate(consequence = ht.vep.most_severe_consequence)

    # Count number of variants and singletons
    ht_reg_table_N_variants = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence, ht_reg_table.gene_id).aggregate(N_variants = hl.agg.count()))
    ht_reg_table_N_singletons = (ht_reg_table.group_by(ht_reg_table.context, ht_reg_table.ref, ht_reg_table.alt, ht_reg_table.mu_snp, ht_reg_table.consequence, ht_reg_table.gene_id).aggregate(N_singletons = hl.agg.sum(ht_reg_table.info.singleton)))

    # Merge the tables to obtain proportions (ps)
    ht_reg_table_ps = ht_reg_table_N_variants.annotate(**ht_reg_table_N_singletons[ht_reg_table_N_variants.context, ht_reg_table_N_variants.ref, ht_reg_table_N_variants.alt, ht_reg_table_N_variants.mu_snp, ht_reg_table_N_variants.consequence, ht_reg_table_N_variants.gene_id])
    ht_reg_table_ps = ht_reg_table_ps.annotate(ps = ht_reg_table_ps.N_singletons/ht_reg_table_ps.N_variants)
    
    # Get expected number of singletons by applying the model factors
    ht_reg_table_ps_lm_cons = ht_reg_table_ps.annotate(expected_singletons=(ht_reg_table_ps.mu_snp * ht_syn_lm[1] + ht_syn_lm[0]) * ht_reg_table_ps.N_variants)

    # To aggregate just sum for the context
    ht_reg_table_ps_lm_cons_agg = (ht_reg_table_ps_lm_cons.group_by("consequence", "gene_id")
                  .aggregate(N_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.N_singletons),
                             expected_singletons=hl.agg.sum(ht_reg_table_ps_lm_cons.expected_singletons),
                             N_variants=hl.agg.sum(ht_reg_table_ps_lm_cons.N_variants)))

    # Calculate MAPS and aggregated proportions 
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg.annotate(ps_agg=ht_reg_table_ps_lm_cons_agg.N_singletons / ht_reg_table_ps_lm_cons_agg.N_variants,
        maps=(ht_reg_table_ps_lm_cons_agg.N_singletons - ht_reg_table_ps_lm_cons_agg.expected_singletons) / ht_reg_table_ps_lm_cons_agg.N_variants)

    # Add MAPS standard error of the mean (sem)
    ht_reg_table_ps_lm_cons_agg_MAPS = ht_reg_table_ps_lm_cons_agg_MAPS.annotate(maps_sem=(ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg * (1 - ht_reg_table_ps_lm_cons_agg_MAPS.ps_agg) / ht_reg_table_ps_lm_cons_agg_MAPS.N_variants) ** 0.5)
    return ht_reg_table_ps_lm_cons_agg_MAPS