import hail as hl
from bokeh.io import output_notebook,show
import gnomad.utils.vep
from hail.ggplot import *
import plotly
import plotly.io as pio
pio.renderers.default='iframe'

def main():
    ht_final = hl.read_table('gs://janucik-dataproc-stage/01_maps/data_full_ht/02a_f_ht_final.ht')

    # Import mutation rates from gnomAD paper
    ht_mu = hl.import_table('gs://janucik-dataproc-stage/01_maps/supplementary_dataset_10_mutation_rates.tsv.gz',
                    delimiter='\t', impute=True, force_bgz=True)

    # Fix the keys
    ht_mu = ht_mu.key_by('context', 'ref', 'alt', 'methylation_level') # has to have a key in order to join using foreign key

    # Annotate mu_snp in the final ht that is fed into main functions
    ht_final = ht_final.annotate(**ht_mu[ht_final.context, ht_final.ref, ht_final.alt, ht_final.methylation_level])

    # Subset into specific variant types
    ht_f_inframe_del = ht_final.filter(ht_final.vep.most_severe_consequence == "inframe_deletion")
    ht_f_inframe_del = ht_f_inframe_del.select(ht_f_inframe_del.vep.most_severe_consequence, ht_f_inframe_del.methyl_mean, ht_f_inframe_del.methylation_level, ht_f_inframe_del.mu_snp, ht_f_inframe_del.context, ht_f_inframe_del.ref, ht_f_inframe_del.alt)
    print('inframe deletion')
    ht_f_inframe_del.show()

    ht_f_inframe_ins = ht_final.filter(ht_final.vep.most_severe_consequence == "inframe_insertion")
    ht_f_inframe_ins = ht_f_inframe_ins.select(ht_f_inframe_ins.vep.most_severe_consequence, ht_f_inframe_ins.methyl_mean, ht_f_inframe_ins.methylation_level, ht_f_inframe_ins.mu_snp, ht_f_inframe_ins.context, ht_f_inframe_ins.ref, ht_f_inframe_ins.alt)
    print('inframe insertion')
    ht_f_inframe_ins.show()

if __name__ == '__main__':
    main()