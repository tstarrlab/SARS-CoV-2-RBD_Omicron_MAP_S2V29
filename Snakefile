"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_rulegraph,
            make_summary,
            get_mut_bind_expr,
            get_early2020_mut_bind_expr,
            save_pinned_env

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Global variables extracted from config --------------------------------------
# Information on samples and barcode runs for Wuhan_Hu_1 ----------------------
barcode_runs_Wuhan_Hu_1 = pd.read_csv(config['barcode_runs_Wuhan_Hu_1'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Wuhan_Hu_1.groupby(['library', 'sample'])) == len(barcode_runs_Wuhan_Hu_1)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Wuhan_Hu_1 = (
    barcode_runs_Wuhan_Hu_1
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Wuhan_Hu_1['equal'].all(), sample_vs_expect_Wuhan_Hu_1.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Wuhan_Hu_1_expandR1 = (
    barcode_runs_Wuhan_Hu_1
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Wuhan_Hu_1_expandR1['sample_lib'].nunique() == len(barcode_runs_Wuhan_Hu_1_expandR1)
if any(barcode_runs_Wuhan_Hu_1_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Wuhan_Hu_1_expandR1.query('n_R1 < 1')}")
    

# Information on samples and barcode runs for Omicron_BA2 ----------------------
barcode_runs_Omicron_BA2 = pd.read_csv(config['barcode_runs_Omicron_BA2'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_BA2.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_BA2)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_BA2 = (
    barcode_runs_Omicron_BA2
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_BA2['equal'].all(), sample_vs_expect_Omicron_BA2.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_BA2_expandR1 = (
    barcode_runs_Omicron_BA2
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_BA2_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_BA2_expandR1)
if any(barcode_runs_Omicron_BA2_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_BA2_expandR1.query('n_R1 < 1')}")

# Information on samples and barcode runs for Omicron_BQ11 ----------------------
barcode_runs_Omicron_BQ11 = pd.read_csv(config['barcode_runs_Omicron_BQ11'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_BQ11.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_BQ11)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_BQ11 = (
    barcode_runs_Omicron_BQ11
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_BQ11['equal'].all(), sample_vs_expect_Omicron_BQ11.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_BQ11_expandR1 = (
    barcode_runs_Omicron_BQ11
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_BQ11_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_BQ11_expandR1)
if any(barcode_runs_Omicron_BQ11_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_BQ11_expandR1.query('n_R1 < 1')}")

# Information on samples and barcode runs for Omicron_XBB15 ----------------------
barcode_runs_Omicron_XBB15 = pd.read_csv(config['barcode_runs_Omicron_XBB15'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_XBB15.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_XBB15)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_XBB15 = (
    barcode_runs_Omicron_XBB15
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_XBB15['equal'].all(), sample_vs_expect_Omicron_XBB15.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_XBB15_expandR1 = (
    barcode_runs_Omicron_XBB15
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_XBB15_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_XBB15_expandR1)
if any(barcode_runs_Omicron_XBB15_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_XBB15_expandR1.query('n_R1 < 1')}")

# Information on samples and barcode runs for Omicron_EG5 ----------------------
barcode_runs_Omicron_EG5 = pd.read_csv(config['barcode_runs_Omicron_EG5'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_EG5.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_EG5)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_EG5 = (
    barcode_runs_Omicron_EG5
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_EG5['equal'].all(), sample_vs_expect_Omicron_EG5.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_EG5_expandR1 = (
    barcode_runs_Omicron_EG5
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_EG5_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_EG5_expandR1)
if any(barcode_runs_Omicron_EG5_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_EG5_expandR1.query('n_R1 < 1')}")

# Information on samples and barcode runs for Omicron_FLip ----------------------
barcode_runs_Omicron_FLip = pd.read_csv(config['barcode_runs_Omicron_FLip'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_FLip.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_FLip)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_FLip = (
    barcode_runs_Omicron_FLip
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_FLip['equal'].all(), sample_vs_expect_Omicron_FLip.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_FLip_expandR1 = (
    barcode_runs_Omicron_FLip
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_FLip_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_FLip_expandR1)
if any(barcode_runs_Omicron_FLip_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_FLip_expandR1.query('n_R1 < 1')}")

# Information on samples and barcode runs for Omicron_BA286 ----------------------
barcode_runs_Omicron_BA286 = pd.read_csv(config['barcode_runs_Omicron_BA286'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs_Omicron_BA286.groupby(['library', 'sample'])) == len(barcode_runs_Omicron_BA286)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *sort_bin*.
sample_vs_expect_Omicron_BA286 = (
    barcode_runs_Omicron_BA286
    .query('experiment_type=="ab_selection"')
    .assign(concentration=lambda x: x['concentration'].astype(int),
            expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'sort_bin']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect_Omicron_BA286['equal'].all(), sample_vs_expect_Omicron_BA286.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_Omicron_BA286_expandR1 = (
    barcode_runs_Omicron_BA286
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )

assert barcode_runs_Omicron_BA286_expandR1['sample_lib'].nunique() == len(barcode_runs_Omicron_BA286_expandR1)
if any(barcode_runs_Omicron_BA286_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_Omicron_BA286_expandR1.query('n_R1 < 1')}")


# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        rulegraph=os.path.join(config['summary_dir'], 'rulegraph.svg'),
        env='environment_pinned.yml',
        get_mut_bind_expr=config['mut_bind_expr'],
        get_early2020_mut_bind_expr=config['early2020_mut_bind_expr'],
        
        bind_expr_filters_Wuhan_Hu_1=nb_markdown('bind_expr_filters_Wuhan_Hu_1.ipynb'),
        bind_expr_filters_Omicron_BA2=nb_markdown('bind_expr_filters_Omicron_BA2.ipynb'),
        aggregate_variant_counts_Wuhan_Hu_1=nb_markdown('aggregate_variant_counts_Wuhan_Hu_1.ipynb'),
        aggregate_variant_counts_Omicron_BA2=nb_markdown('aggregate_variant_counts_Omicron_BA2.ipynb'),
        variant_counts_Wuhan_Hu_1=config['variant_counts_Wuhan_Hu_1'],
        variant_counts_Omicron_BA2=config['variant_counts_Omicron_BA2'],
        counts_to_cells_ratio_Wuhan_Hu_1=nb_markdown('counts_to_cells_ratio_Wuhan_Hu_1.ipynb'),
        counts_to_cells_ratio_Omicron_BA2=nb_markdown('counts_to_cells_ratio_Omicron_BA2.ipynb'),
        counts_to_cells_csv_Wuhan_Hu_1=config['counts_to_cells_csv_Wuhan_Hu_1'],
        counts_to_cells_csv_Omicron_BA2=config['counts_to_cells_csv_Omicron_BA2'],
        counts_to_scores_Wuhan_Hu_1=nb_markdown('counts_to_scores_Wuhan_Hu_1.ipynb'),
        counts_to_scores_Omicron_BA2=nb_markdown('counts_to_scores_Omicron_BA2.ipynb'),
        escape_scores_Wuhan_Hu_1=config['escape_scores_Wuhan_Hu_1'],
        escape_fracs_Wuhan_Hu_1=config['escape_fracs_Wuhan_Hu_1'],
        escape_scores_Omicron_BA2=config['escape_scores_Omicron_BA2'],
        escape_fracs_Omicron_BA2=config['escape_fracs_Omicron_BA2'],
        call_strong_escape_sites_Wuhan_Hu_1=nb_markdown('call_strong_escape_sites_Wuhan_Hu_1.ipynb'),
        strong_escape_sites_Wuhan_Hu_1=config['strong_escape_sites_Wuhan_Hu_1'],
        call_strong_escape_sites_Omicron_BA2=nb_markdown('call_strong_escape_sites_Omicron_BA2.ipynb'),
        strong_escape_sites_Omicron_BA2=config['strong_escape_sites_Omicron_BA2'],
        escape_profiles_Wuhan_Hu_1=nb_markdown('escape_profiles_Wuhan_Hu_1.ipynb'),
        escape_profiles_Omicron_BA2=nb_markdown('escape_profiles_Omicron_BA2.ipynb'),
        output_pdbs_Wuhan_Hu_1=nb_markdown('output_pdbs_Wuhan_Hu_1.ipynb'),
        output_pdbs_Omicron_BA2=nb_markdown('output_pdbs_Omicron_BA2.ipynb'),
        make_supp_data_Wuhan_Hu_1=nb_markdown('make_supp_data_Wuhan_Hu_1.ipynb'),
        make_supp_data_Omicron_BA2=nb_markdown('make_supp_data_Omicron_BA2.ipynb'),
        
        bind_expr_filters_Omicron_BQ11=nb_markdown('bind_expr_filters_Omicron_BQ11.ipynb'),
        bind_expr_filters_Omicron_XBB15=nb_markdown('bind_expr_filters_Omicron_XBB15.ipynb'),
        aggregate_variant_counts_Omicron_BQ11=nb_markdown('aggregate_variant_counts_Omicron_BQ11.ipynb'),
        aggregate_variant_counts_Omicron_XBB15=nb_markdown('aggregate_variant_counts_Omicron_XBB15.ipynb'),
        variant_counts_Omicron_BQ11=config['variant_counts_Omicron_BQ11'],
        variant_counts_Omicron_XBB15=config['variant_counts_Omicron_XBB15'],
        counts_to_cells_ratio_Omicron_BQ11=nb_markdown('counts_to_cells_ratio_Omicron_BQ11.ipynb'),
        counts_to_cells_ratio_Omicron_XBB15=nb_markdown('counts_to_cells_ratio_Omicron_XBB15.ipynb'),
        counts_to_cells_csv_Omicron_BQ11=config['counts_to_cells_csv_Omicron_BQ11'],
        counts_to_cells_csv_Omicron_XBB15=config['counts_to_cells_csv_Omicron_XBB15'],
        counts_to_scores_Omicron_BQ11=nb_markdown('counts_to_scores_Omicron_BQ11.ipynb'),
        counts_to_scores_Omicron_XBB15=nb_markdown('counts_to_scores_Omicron_XBB15.ipynb'),
        escape_scores_Omicron_BQ11=config['escape_scores_Omicron_BQ11'],
        escape_fracs_Omicron_BQ11=config['escape_fracs_Omicron_BQ11'],
        escape_scores_Omicron_XBB15=config['escape_scores_Omicron_XBB15'],
        escape_fracs_Omicron_XBB15=config['escape_fracs_Omicron_XBB15'],
        call_strong_escape_sites_Omicron_BQ11=nb_markdown('call_strong_escape_sites_Omicron_BQ11.ipynb'),
        strong_escape_sites_Omicron_BQ11=config['strong_escape_sites_Omicron_BQ11'],
        call_strong_escape_sites_Omicron_XBB15=nb_markdown('call_strong_escape_sites_Omicron_XBB15.ipynb'),
        strong_escape_sites_Omicron_XBB15=config['strong_escape_sites_Omicron_XBB15'],
        escape_profiles_Omicron_BQ11=nb_markdown('escape_profiles_Omicron_BQ11.ipynb'),
        escape_profiles_Omicron_XBB15=nb_markdown('escape_profiles_Omicron_XBB15.ipynb'),
        output_pdbs_Omicron_BQ11=nb_markdown('output_pdbs_Omicron_BQ11.ipynb'),
        output_pdbs_Omicron_XBB15=nb_markdown('output_pdbs_Omicron_XBB15.ipynb'),
        make_supp_data_Omicron_BQ11=nb_markdown('make_supp_data_Omicron_BQ11.ipynb'),
        make_supp_data_Omicron_XBB15=nb_markdown('make_supp_data_Omicron_XBB15.ipynb'),
        
        bind_expr_filters_Omicron_EG5=nb_markdown('bind_expr_filters_Omicron_EG5.ipynb'),
        bind_expr_filters_Omicron_FLip=nb_markdown('bind_expr_filters_Omicron_FLip.ipynb'),
        aggregate_variant_counts_Omicron_EG5=nb_markdown('aggregate_variant_counts_Omicron_EG5.ipynb'),
        aggregate_variant_counts_Omicron_FLip=nb_markdown('aggregate_variant_counts_Omicron_FLip.ipynb'),
        variant_counts_Omicron_EG5=config['variant_counts_Omicron_EG5'],
        variant_counts_Omicron_FLip=config['variant_counts_Omicron_FLip'],
        counts_to_cells_ratio_Omicron_EG5=nb_markdown('counts_to_cells_ratio_Omicron_EG5.ipynb'),
        counts_to_cells_ratio_Omicron_FLip=nb_markdown('counts_to_cells_ratio_Omicron_FLip.ipynb'),
        counts_to_cells_csv_Omicron_EG5=config['counts_to_cells_csv_Omicron_EG5'],
        counts_to_cells_csv_Omicron_FLip=config['counts_to_cells_csv_Omicron_FLip'],
        counts_to_scores_Omicron_EG5=nb_markdown('counts_to_scores_Omicron_EG5.ipynb'),
        counts_to_scores_Omicron_FLip=nb_markdown('counts_to_scores_Omicron_FLip.ipynb'),
        escape_scores_Omicron_EG5=config['escape_scores_Omicron_EG5'],
        escape_fracs_Omicron_EG5=config['escape_fracs_Omicron_EG5'],
        escape_scores_Omicron_FLip=config['escape_scores_Omicron_FLip'],
        escape_fracs_Omicron_FLip=config['escape_fracs_Omicron_FLip'],
        call_strong_escape_sites_Omicron_EG5=nb_markdown('call_strong_escape_sites_Omicron_EG5.ipynb'),
        strong_escape_sites_Omicron_EG5=config['strong_escape_sites_Omicron_EG5'],
        call_strong_escape_sites_Omicron_FLip=nb_markdown('call_strong_escape_sites_Omicron_FLip.ipynb'),
        strong_escape_sites_Omicron_FLip=config['strong_escape_sites_Omicron_FLip'],
        escape_profiles_Omicron_EG5=nb_markdown('escape_profiles_Omicron_EG5.ipynb'),
        escape_profiles_Omicron_FLip=nb_markdown('escape_profiles_Omicron_FLip.ipynb'),
        output_pdbs_Omicron_EG5=nb_markdown('output_pdbs_Omicron_EG5.ipynb'),
        output_pdbs_Omicron_FLip=nb_markdown('output_pdbs_Omicron_FLip.ipynb'),
        make_supp_data_Omicron_EG5=nb_markdown('make_supp_data_Omicron_EG5.ipynb'),
        make_supp_data_Omicron_FLip=nb_markdown('make_supp_data_Omicron_FLip.ipynb'),
        
        bind_expr_filters_Omicron_BA286=nb_markdown('bind_expr_filters_Omicron_BA286.ipynb'),
        aggregate_variant_counts_Omicron_BA286=nb_markdown('aggregate_variant_counts_Omicron_BA286.ipynb'),
        variant_counts_Omicron_BA286=config['variant_counts_Omicron_BA286'],
        counts_to_cells_ratio_Omicron_BA286=nb_markdown('counts_to_cells_ratio_Omicron_BA286.ipynb'),
        counts_to_cells_csv_Omicron_BA286=config['counts_to_cells_csv_Omicron_BA286'],
        counts_to_scores_Omicron_BA286=nb_markdown('counts_to_scores_Omicron_BA286.ipynb'),
        escape_scores_Omicron_BA286=config['escape_scores_Omicron_BA286'],
        escape_fracs_Omicron_BA286=config['escape_fracs_Omicron_BA286'],
        call_strong_escape_sites_Omicron_BA286=nb_markdown('call_strong_escape_sites_Omicron_BA286.ipynb'),
        strong_escape_sites_Omicron_BA286=config['strong_escape_sites_Omicron_BA286'],
        escape_profiles_Omicron_BA286=nb_markdown('escape_profiles_Omicron_BA286.ipynb'),
        output_pdbs_Omicron_BA286=nb_markdown('output_pdbs_Omicron_BA286.ipynb'),
        make_supp_data_Omicron_BA286=nb_markdown('make_supp_data_Omicron_BA286.ipynb'),
        
        gisaid_muts=nb_markdown('gisaid_rbd_mutations.ipynb'),
        custom_plots='results/summary/custom-plots.md',
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the rule graph of the computational workflow:
            ![{path(input.rulegraph)}]({path(input.rulegraph)})

            Here is the Markdown output of each notebook in the workflow:
            
            1. Get prior RBD DMS mutation-level binding and expression measurements and barcode-variant lookup table from the [SARS-CoV-2-RBD_DMS_Omicron repository](https://github.com/tstarrlab/SARS-CoV-2-RBD_DMS_Omicron-XBB-BQ) and the original DMS library for SARS-CoV-2 (PCR-based mutagenesis) [here](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS). 

            2. Count variants and then aggregate counts for
               [Wuhan_Hu_1]({path(input.aggregate_variant_counts_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.aggregate_variant_counts_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.aggregate_variant_counts_Omicron_BQ11)}),
               [Omicron_XBB15]({path(input.aggregate_variant_counts_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.aggregate_variant_counts_Omicron_EG5)}),
               [Omicron_FLip]({path(input.aggregate_variant_counts_Omicron_FLip)}),
               [Omicron_BA286]({path(input.aggregate_variant_counts_Omicron_BA286)})
               to create variant counts files for [Wuhan_Hu_1]({path(input.variant_counts_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.variant_counts_Omicron_BA2)}), 
               [Omicron_BQ11]({path(input.variant_counts_Omicron_BQ11)}),
               [Omicron_XBB15]({path(input.variant_counts_Omicron_XBB15)}), 
               [Omicron_EG5]({path(input.variant_counts_Omicron_EG5)}), 
               [Omicron_FLip]({path(input.variant_counts_Omicron_FLip)}), and
               [Omicron_BA286]({path(input.variant_counts_Omicron_BA286)}).

            3. Analyze sequencing counts to cells ratio for [Wuhan_Hu_1]({path(input.counts_to_cells_ratio_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.counts_to_cells_ratio_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.counts_to_cells_ratio_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.counts_to_cells_ratio_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.counts_to_cells_ratio_Omicron_EG5)}),
               [Omicron_FLip]({path(input.counts_to_cells_ratio_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.counts_to_cells_ratio_Omicron_BA286)}).
               this prints a list of any samples where this ratio too low. Also
               creates a CSV for [Wuhan_Hu_1]({path(input.counts_to_cells_csv_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.counts_to_cells_csv_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.counts_to_cells_csv_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.counts_to_cells_csv_Omicron_XBB15)}), 
               [Omicron_EG5]({path(input.counts_to_cells_csv_Omicron_EG5)}),
               [Omicron_FLip]({path(input.counts_to_cells_csv_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.counts_to_cells_csv_Omicron_BA286)}) with the
               sequencing counts, number of sorted cells, and ratios for
               all samples.

            4. Calculate escape scores from variant counts for [Wuhan_Hu_1]({path(input.counts_to_scores_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.counts_to_scores_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.counts_to_scores_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.counts_to_scores_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.counts_to_scores_Omicron_EG5)}),
               [Omicron_FLip]({path(input.counts_to_scores_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.counts_to_scores_Omicron_BA286)}).

            5. Call sites of strong escape for [Wuhan_Hu_1]({path(input.call_strong_escape_sites_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.call_strong_escape_sites_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.call_strong_escape_sites_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.call_strong_escape_sites_Omicron_XBB15)}), 
               [Omicron_EG5]({path(input.call_strong_escape_sites_Omicron_EG5)}),
               [Omicron_FLip]({path(input.call_strong_escape_sites_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.call_strong_escape_sites_Omicron_BA286)}).

            6. Plot escape profiles for [Wuhan_Hu_1]({path(input.escape_profiles_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.escape_profiles_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.escape_profiles_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.escape_profiles_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.escape_profiles_Omicron_EG5)}),
               [Omicron_FLip]({path(input.escape_profiles_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.escape_profiles_Omicron_BA286)}).

            7. Map escape profiles to ``*.pdb`` files using notebooks here for 
               [Wuhan_Hu_1]({path(input.output_pdbs_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.output_pdbs_Omicron_BA2)}),
               [Omicron_BQ11]({path(input.output_pdbs_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.output_pdbs_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.output_pdbs_Omicron_EG5)}),
               [Omicron_FLip]({path(input.output_pdbs_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.output_pdbs_Omicron_BA286)}).
            
            8. Compare DMS escape to circulating variants from GISAID. First requires manually running
               the query_GISAID.Rmd script to download sequences from GISAID (doesn't run on cluster), followed by
               [this notebook]({path(input.gisaid_muts)}) to align and call mutations and [this notebook]({path(input.custom_plots)})
               to make plots. 
            
           9. Make supplementary data files for [Wuhan_Hu_1]({path(input.make_supp_data_Wuhan_Hu_1)}), 
               [Omicron_BA2]({path(input.make_supp_data_Omicron_BA2)}), 
               [Omicron_BQ11]({path(input.make_supp_data_Omicron_BQ11)}), 
               [Omicron_XBB15]({path(input.make_supp_data_Omicron_XBB15)}),
               [Omicron_EG5]({path(input.make_supp_data_Omicron_EG5)}), 
               [Omicron_FLip]({path(input.make_supp_data_Omicron_FLip)}) and 
               [Omicron_BA286]({path(input.make_supp_data_Omicron_BA286)})
               which are here for [Wuhan_Hu_1]({path(config['supp_data_dir_Wuhan_Hu_1'])}), 
               [Omicron_BA2]({path(config['supp_data_dir_Omicron_BA2'])}), 
               [Omicron_BQ11]({path(config['supp_data_dir_Omicron_BQ11'])}), 
               [Omicron_XBB15]({path(config['supp_data_dir_Omicron_XBB15'])}), 
               [Omicron_EG5]({path(config['supp_data_dir_Omicron_EG5'])}), 
               [Omicron_FLip]({path(config['supp_data_dir_Omicron_FLip'])}) and 
               [Omicron_BA286]({path(config['supp_data_dir_Omicron_BA286'])}). These include
               `dms-view` input files.
               
            


            """
            ).strip())


rule make_rulegraph:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'rulegraph.svg')
    shell:
        "snakemake --forceall --rulegraph | dot -Tsvg > {output}"

rule save_pinned_env:
    input:
        dag=os.path.join(config['summary_dir'], 'rulegraph.svg'),
    log:
    	"environment_pinned.yml"
    shell:
        """
        conda env export > {log}
        """

rule make_supp_data_Wuhan_Hu_1:
    input:
        config['escape_profiles_config_Wuhan_Hu_1'],
        config['output_pdbs_config'],
        config['escape_fracs_Wuhan_Hu_1'],
        config['escape_profiles_dms_colors_Wuhan_Hu_1']
    output:
        nb_markdown=nb_markdown('make_supp_data_Wuhan_Hu_1.ipynb'),
        outdir=directory(config['supp_data_dir_Wuhan_Hu_1']),
    params:
        nb='make_supp_data_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule make_supp_data_Omicron_BA2:
    input:
        config['escape_profiles_config_Omicron_BA2'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_BA2'],
        config['escape_profiles_dms_colors_Omicron_BA2']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_BA2.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_BA2']),
    params:
        nb='make_supp_data_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule make_supp_data_Omicron_BQ11:
    input:
        config['escape_profiles_config_Omicron_BQ11'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_BQ11'],
        config['escape_profiles_dms_colors_Omicron_BQ11']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_BQ11.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_BQ11']),
    params:
        nb='make_supp_data_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule make_supp_data_Omicron_XBB15:
    input:
        config['escape_profiles_config_Omicron_XBB15'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_XBB15'],
        config['escape_profiles_dms_colors_Omicron_XBB15']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_XBB15.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_XBB15']),
    params:
        nb='make_supp_data_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule make_supp_data_Omicron_EG5:
    input:
        config['escape_profiles_config_Omicron_EG5'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_EG5'],
        config['escape_profiles_dms_colors_Omicron_EG5']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_EG5.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_EG5']),
    params:
        nb='make_supp_data_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

rule make_supp_data_Omicron_FLip:
    input:
        config['escape_profiles_config_Omicron_FLip'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_FLip'],
        config['escape_profiles_dms_colors_Omicron_FLip']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_FLip.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_FLip']),
    params:
        nb='make_supp_data_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule make_supp_data_Omicron_BA286:
    input:
        config['escape_profiles_config_Omicron_BA286'],
        config['output_pdbs_config'],
        config['escape_fracs_Omicron_BA286'],
        config['escape_profiles_dms_colors_Omicron_BA286']
    output:
        nb_markdown=nb_markdown('make_supp_data_Omicron_BA286.ipynb'),
        outdir=directory(config['supp_data_dir_Omicron_BA286']),
    params:
        nb='make_supp_data_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"        
        
        

rule custom_plots:
    input:
        config['escape_fracs_Omicron_BA2'],
        config['gisaid_mutation_counts'],
        'custom-plots.Rmd'
    output:
        md='results/summary/custom-plots.md',
        md_files=directory('results/summary/custom-plots_files'),
    envmodules:
        'R/4.1.3'
    params:
        nb='custom-plots.Rmd',
        md='custom-plots.md',
        md_files='custom-plots_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule gisaid_rbd_mutations:
    input:
        config['spike_sequences'],
        config['RBD_sequence_Omicron_BA2'],
    output:
        config['gisaid_mutation_counts'],
        nb_markdown=nb_markdown('gisaid_rbd_mutations.ipynb'),
    params:
        nb='gisaid_rbd_mutations.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

#rule download_gisaid:
#    input:
#        'query_GISAID.Rmd'
#    output:
#        md='results/summary/query_GISAID.md',
#        fasta=config['spike_sequences'],
#    envmodules:
#        'R/4.1.3'
#    params:
#        nb='query_GISAID.Rmd',
#        md='query_GISAID.md',
#    shell:
#        """
#        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
#        mv {params.md} {output.md};
#        """





rule output_pdbs_Wuhan_Hu_1:
    input:
        config['escape_fracs_Wuhan_Hu_1'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Wuhan_Hu_1.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Wuhan_Hu_1']),
    params:
        nb='output_pdbs_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
                
rule output_pdbs_Omicron_BA2:
    input:
        config['escape_fracs_Omicron_BA2'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_BA2.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_BA2']),
    params:
        nb='output_pdbs_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule output_pdbs_Omicron_BQ11:
    input:
        config['escape_fracs_Omicron_BQ11'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_BQ11.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_BQ11']),
    params:
        nb='output_pdbs_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
                
rule output_pdbs_Omicron_XBB15:
    input:
        config['escape_fracs_Omicron_XBB15'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_XBB15.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_XBB15']),
    params:
        nb='output_pdbs_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule output_pdbs_Omicron_EG5:
    input:
        config['escape_fracs_Omicron_EG5'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_EG5.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_EG5']),
    params:
        nb='output_pdbs_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule output_pdbs_Omicron_FLip:
    input:
        config['escape_fracs_Omicron_FLip'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_FLip.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_FLip']),
    params:
        nb='output_pdbs_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
                
rule output_pdbs_Omicron_BA286:
    input:
        config['escape_fracs_Omicron_BA286'],
        config['output_pdbs_config'],
    output:
        nb_markdown=nb_markdown('output_pdbs_Omicron_BA286.ipynb'),
        outdir=directory(config['pdb_outputs_dir_Omicron_BA286']),
    params:
        nb='output_pdbs_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
          
        

rule escape_profiles_Wuhan_Hu_1:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Wuhan_Hu_1'],
        escape_profiles_config=config['escape_profiles_config_Wuhan_Hu_1'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Wuhan_Hu_1'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Wuhan_Hu_1'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Wuhan_Hu_1.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Wuhan_Hu_1'],
    params:
        nb='escape_profiles_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule escape_profiles_Omicron_BA2:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA2'],
        escape_profiles_config=config['escape_profiles_config_Omicron_BA2'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_BA2'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_BA2.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_BA2'],
    params:
        nb='escape_profiles_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule escape_profiles_Omicron_BQ11:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BQ11'],
        escape_profiles_config=config['escape_profiles_config_Omicron_BQ11'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_BQ11'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_BQ11'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_BQ11.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_BQ11'],
    params:
        nb='escape_profiles_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule escape_profiles_Omicron_XBB15:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_XBB15'],
        escape_profiles_config=config['escape_profiles_config_Omicron_XBB15'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_XBB15'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_XBB15'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_XBB15.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_XBB15'],
    params:
        nb='escape_profiles_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule escape_profiles_Omicron_EG5:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_EG5'],
        escape_profiles_config=config['escape_profiles_config_Omicron_EG5'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_EG5'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_EG5'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_EG5.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_EG5'],
    params:
        nb='escape_profiles_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule escape_profiles_Omicron_FLip:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_FLip'],
        escape_profiles_config=config['escape_profiles_config_Omicron_FLip'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_FLip'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_FLip'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_FLip.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_FLip'],
    params:
        nb='escape_profiles_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule escape_profiles_Omicron_BA286:
    """Make stacked logo plots of antibody escape profiles."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA286'],
        escape_profiles_config=config['escape_profiles_config_Omicron_BA286'],
        site_color_schemes=config['site_color_schemes'],
        wildtype_sequence=config['wildtype_sequence_Omicron_BA286'],
        mut_bind_expr=config['mut_bind_expr'],
        strong_escape_sites=config['strong_escape_sites_Omicron_BA286'],
    output:
        nb_markdown=nb_markdown('escape_profiles_Omicron_BA286.ipynb'),
        escape_profiles_dms_colors=config['escape_profiles_dms_colors_Omicron_BA286'],
    params:
        nb='escape_profiles_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
          
        
        
        
        
        

rule call_strong_escape_sites_Wuhan_Hu_1:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Wuhan_Hu_1'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Wuhan_Hu_1.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Wuhan_Hu_1'],
    params:
        nb='call_strong_escape_sites_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule call_strong_escape_sites_Omicron_BA2:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_BA2.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_BA2'],
    params:
        nb='call_strong_escape_sites_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule call_strong_escape_sites_Omicron_BQ11:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BQ11'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_BQ11.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_BQ11'],
    params:
        nb='call_strong_escape_sites_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule call_strong_escape_sites_Omicron_XBB15:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_XBB15'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_XBB15.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_XBB15'],
    params:
        nb='call_strong_escape_sites_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule call_strong_escape_sites_Omicron_EG5:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_EG5'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_EG5.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_EG5'],
    params:
        nb='call_strong_escape_sites_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule call_strong_escape_sites_Omicron_FLip:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_FLip'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_FLip.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_FLip'],
    params:
        nb='call_strong_escape_sites_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule call_strong_escape_sites_Omicron_BA286:
    """Call sites of strong escape."""
    input:
        escape_fracs=config['escape_fracs_Omicron_BA286'],
    output:
        nb_markdown=nb_markdown('call_strong_escape_sites_Omicron_BA286.ipynb'),
        strong_escape_sites=config['strong_escape_sites_Omicron_BA286'],
    params:
        nb='call_strong_escape_sites_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"     
        
        
        

rule counts_to_scores_Wuhan_Hu_1:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Wuhan_Hu_1'],
        config['wildtype_sequence_Wuhan_Hu_1'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Wuhan_Hu_1.ipynb'),
        escape_scores=config['escape_scores_Wuhan_Hu_1'],
        escape_fracs=config['escape_fracs_Wuhan_Hu_1'],
        escape_score_samples=config['escape_score_samples_Wuhan_Hu_1'],
    params:
        nb='counts_to_scores_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule counts_to_scores_Omicron_BA2:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_BA2'],
        config['wildtype_sequence_Omicron_BA2'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_BA2.ipynb'),
        escape_scores=config['escape_scores_Omicron_BA2'],
        escape_fracs=config['escape_fracs_Omicron_BA2'],
        escape_score_samples=config['escape_score_samples_Omicron_BA2'],
    params:
        nb='counts_to_scores_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_scores_Omicron_BQ11:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_BQ11'],
        config['wildtype_sequence_Omicron_BQ11'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_BQ11.ipynb'),
        escape_scores=config['escape_scores_Omicron_BQ11'],
        escape_fracs=config['escape_fracs_Omicron_BQ11'],
        escape_score_samples=config['escape_score_samples_Omicron_BQ11'],
    params:
        nb='counts_to_scores_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule counts_to_scores_Omicron_XBB15:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_XBB15'],
        config['wildtype_sequence_Omicron_XBB15'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_XBB15.ipynb'),
        escape_scores=config['escape_scores_Omicron_XBB15'],
        escape_fracs=config['escape_fracs_Omicron_XBB15'],
        escape_score_samples=config['escape_score_samples_Omicron_XBB15'],
    params:
        nb='counts_to_scores_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
        

rule counts_to_scores_Omicron_EG5:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_EG5'],
        config['wildtype_sequence_Omicron_EG5'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_EG5.ipynb'),
        escape_scores=config['escape_scores_Omicron_EG5'],
        escape_fracs=config['escape_fracs_Omicron_EG5'],
        escape_score_samples=config['escape_score_samples_Omicron_EG5'],
    params:
        nb='counts_to_scores_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_scores_Omicron_FLip:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_FLip'],
        config['wildtype_sequence_Omicron_FLip'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_FLip.ipynb'),
        escape_scores=config['escape_scores_Omicron_FLip'],
        escape_fracs=config['escape_fracs_Omicron_FLip'],
        escape_score_samples=config['escape_score_samples_Omicron_FLip'],
    params:
        nb='counts_to_scores_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"


rule counts_to_scores_Omicron_BA286:
    """Analyze variant counts to compute escape scores."""
    input:
        config['variant_counts_Omicron_BA286'],
        config['wildtype_sequence_Omicron_BA286'],
        # config['mut_bind_expr'],
        # config['variant_expr'],
        # config['variant_bind'],
    output:
        nb_markdown=nb_markdown('counts_to_scores_Omicron_BA286.ipynb'),
        escape_scores=config['escape_scores_Omicron_BA286'],
        escape_fracs=config['escape_fracs_Omicron_BA286'],
        escape_score_samples=config['escape_score_samples_Omicron_BA286'],
    params:
        nb='counts_to_scores_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
        
        


rule counts_to_cells_ratio_Wuhan_Hu_1:
    input:
        config['variant_counts_Wuhan_Hu_1'],
        config['barcode_runs_Wuhan_Hu_1'],
        config['wildtype_sequence_Wuhan_Hu_1'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Wuhan_Hu_1.ipynb'),
        counts_to_cells_csv_Wuhan_Hu_1=config['counts_to_cells_csv_Wuhan_Hu_1'],
    params:
        nb='counts_to_cells_ratio_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_cells_ratio_Omicron_BA2:
    input:
        config['variant_counts_Omicron_BA2'],
        config['barcode_runs_Omicron_BA2'],
        config['wildtype_sequence_Omicron_BA2'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_BA2.ipynb'),
        counts_to_cells_csv_Omicron_BA2=config['counts_to_cells_csv_Omicron_BA2'],
    params:
        nb='counts_to_cells_ratio_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule counts_to_cells_ratio_Omicron_BQ11:
    input:
        config['variant_counts_Omicron_BQ11'],
        config['barcode_runs_Omicron_BQ11'],
        config['wildtype_sequence_Omicron_BQ11'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_BQ11.ipynb'),
        counts_to_cells_csv_Omicron_BQ11=config['counts_to_cells_csv_Omicron_BQ11'],
    params:
        nb='counts_to_cells_ratio_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_cells_ratio_Omicron_XBB15:
    input:
        config['variant_counts_Omicron_XBB15'],
        config['barcode_runs_Omicron_XBB15'],
        config['wildtype_sequence_Omicron_XBB15'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_XBB15.ipynb'),
        counts_to_cells_csv_Omicron_XBB15=config['counts_to_cells_csv_Omicron_XBB15'],
    params:
        nb='counts_to_cells_ratio_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_cells_ratio_Omicron_EG5:
    input:
        config['variant_counts_Omicron_EG5'],
        config['barcode_runs_Omicron_EG5'],
        config['wildtype_sequence_Omicron_EG5'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_EG5.ipynb'),
        counts_to_cells_csv_Omicron_EG5=config['counts_to_cells_csv_Omicron_EG5'],
    params:
        nb='counts_to_cells_ratio_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule counts_to_cells_ratio_Omicron_FLip:
    input:
        config['variant_counts_Omicron_FLip'],
        config['barcode_runs_Omicron_FLip'],
        config['wildtype_sequence_Omicron_FLip'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_FLip.ipynb'),
        counts_to_cells_csv_Omicron_FLip=config['counts_to_cells_csv_Omicron_FLip'],
    params:
        nb='counts_to_cells_ratio_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule counts_to_cells_ratio_Omicron_BA286:
    input:
        config['variant_counts_Omicron_BA286'],
        config['barcode_runs_Omicron_BA286'],
        config['wildtype_sequence_Omicron_BA286'],
    output:
        nb_markdown=nb_markdown('counts_to_cells_ratio_Omicron_BA286.ipynb'),
        counts_to_cells_csv_Omicron_BA286=config['counts_to_cells_csv_Omicron_BA286'],
    params:
        nb='counts_to_cells_ratio_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"     
        
        
        
        
        

rule aggregate_variant_counts_Wuhan_Hu_1:
    input:
        counts=expand(os.path.join(config['counts_dir_Wuhan_Hu_1'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Wuhan_Hu_1_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Wuhan_Hu_1'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Wuhan_Hu_1_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Wuhan_Hu_1'],
        wt_seq=config['wildtype_sequence_Wuhan_Hu_1'],
        barcode_runs=config['barcode_runs_Wuhan_Hu_1'],
    output:
        config['variant_counts_Wuhan_Hu_1'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Wuhan_Hu_1.ipynb')
    params:
        nb='aggregate_variant_counts_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule aggregate_variant_counts_Omicron_BA2:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_BA2'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_BA2_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_BA2'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_BA2_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_BA2'],
        wt_seq=config['wildtype_sequence_Omicron_BA2'],
        barcode_runs=config['barcode_runs_Omicron_BA2'],
    output:
        config['variant_counts_Omicron_BA2'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_BA2.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule aggregate_variant_counts_Omicron_BQ11:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_BQ11'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_BQ11_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_BQ11'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_BQ11_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_BQ11'],
        wt_seq=config['wildtype_sequence_Omicron_BQ11'],
        barcode_runs=config['barcode_runs_Omicron_BQ11'],
    output:
        config['variant_counts_Omicron_BQ11'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_BQ11.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule aggregate_variant_counts_Omicron_XBB15:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_XBB15'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_XBB15_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_XBB15'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_XBB15_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_XBB15'],
        wt_seq=config['wildtype_sequence_Omicron_XBB15'],
        barcode_runs=config['barcode_runs_Omicron_XBB15'],
    output:
        config['variant_counts_Omicron_XBB15'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_XBB15.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
       
rule aggregate_variant_counts_Omicron_EG5:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_EG5'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_EG5_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_EG5'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_EG5_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_EG5'],
        wt_seq=config['wildtype_sequence_Omicron_EG5'],
        barcode_runs=config['barcode_runs_Omicron_EG5'],
    output:
        config['variant_counts_Omicron_EG5'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_EG5.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
rule aggregate_variant_counts_Omicron_FLip:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_FLip'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_FLip_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_FLip'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_FLip_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_FLip'],
        wt_seq=config['wildtype_sequence_Omicron_FLip'],
        barcode_runs=config['barcode_runs_Omicron_FLip'],
    output:
        config['variant_counts_Omicron_FLip'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_FLip.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule aggregate_variant_counts_Omicron_BA286:
    input:
        counts=expand(os.path.join(config['counts_dir_Omicron_BA286'],
                                   "{sample_lib}_counts.csv"),
                      sample_lib=barcode_runs_Omicron_BA286_expandR1['sample_lib']),
        fates=expand(os.path.join(config['counts_dir_Omicron_BA286'],
                                  "{sample_lib}_fates.csv"),
                     sample_lib=barcode_runs_Omicron_BA286_expandR1['sample_lib']),
        variant_table=config['bc_variant_lookup_Omicron_BA286'],
        wt_seq=config['wildtype_sequence_Omicron_BA286'],
        barcode_runs=config['barcode_runs_Omicron_BA286'],
    output:
        config['variant_counts_Omicron_BA286'],
        nb_markdown=nb_markdown('aggregate_variant_counts_Omicron_BA286.ipynb')
    params:
        nb='aggregate_variant_counts_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
          
        
        
        

rule count_variants_Wuhan_Hu_1:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Wuhan_Hu_1'],
        wt_seq=config['wildtype_sequence_Wuhan_Hu_1'],
        r1s=lambda wildcards: (barcode_runs_Wuhan_Hu_1_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Wuhan_Hu_1'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Wuhan_Hu_1'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Wuhan_Hu_1_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Wuhan_Hu_1_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule count_variants_Omicron_BA2:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_BA2'],
        wt_seq=config['wildtype_sequence_Omicron_BA2'],
        r1s=lambda wildcards: (barcode_runs_Omicron_BA2_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_BA2'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_BA2'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_BA2_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_BA2_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule count_variants_Omicron_BQ11:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_BQ11'],
        wt_seq=config['wildtype_sequence_Omicron_BQ11'],
        r1s=lambda wildcards: (barcode_runs_Omicron_BQ11_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_BQ11'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_BQ11'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_BQ11_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_BQ11_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target,
                    allowgaps=True)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule count_variants_Omicron_XBB15:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_XBB15'],
        wt_seq=config['wildtype_sequence_Omicron_XBB15'],
        r1s=lambda wildcards: (barcode_runs_Omicron_XBB15_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_XBB15'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_XBB15'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_XBB15_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_XBB15_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target,
                    allowgaps=True)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule count_variants_Omicron_BA286:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_BA286'],
        wt_seq=config['wildtype_sequence_Omicron_BA286'],
        r1s=lambda wildcards: (barcode_runs_Omicron_BA286_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_BA286'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_BA286'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_BA286_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_BA286_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target,
                    allowgaps=True)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)
        
rule count_variants_Omicron_EG5:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_EG5'],
        wt_seq=config['wildtype_sequence_Omicron_EG5'],
        r1s=lambda wildcards: (barcode_runs_Omicron_EG5_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_EG5'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_EG5'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_EG5_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_EG5_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target,
                    allowgaps=True)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)        
        
rule count_variants_Omicron_FLip:
    """Count variants for a specific sample."""
    input:
        variant_table=config['bc_variant_lookup_Omicron_FLip'],
        wt_seq=config['wildtype_sequence_Omicron_FLip'],
        r1s=lambda wildcards: (barcode_runs_Omicron_FLip_expandR1
                               .set_index('sample_lib')
                               .at[wildcards.sample_lib, 'R1']
                               ),
    output:
        counts=os.path.join(config['counts_dir_Omicron_FLip'], "{sample_lib}_counts.csv"),
        fates=os.path.join(config['counts_dir_Omicron_FLip'], "{sample_lib}_fates.csv"),
    params:
        sample_lib="{sample_lib}"
    run:
        # parse sample and library from `sample_lib` wildcard
        lib = params.sample_lib.split('_')[-1]
        sample = params.sample_lib[: -len(lib) - 1]
        assert sample == (barcode_runs_Omicron_FLip_expandR1
                          .set_index('sample_lib')
                          .at[params.sample_lib, 'sample']
                          )
        assert lib == (barcode_runs_Omicron_FLip_expandR1
                       .set_index('sample_lib')
                       .at[params.sample_lib, 'library']
                       )
        # initialize `CodonVariantTable` (used to get valid barcodes)
        wt_seqrecord = Bio.SeqIO.read(input.wt_seq, 'fasta')
        geneseq = str(wt_seqrecord.seq)
        primary_target = wt_seqrecord.name
        variants=dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=input.variant_table,
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target,
                    allowgaps=True)
        # initialize `IlluminaBarcodeParser`
        parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    valid_barcodes=variants.valid_barcodes(lib),
                    **config['illumina_barcode_parser_params'])
        # parse barcodes
        counts, fates = parser.parse(input.r1s,
                                     add_cols={'library': lib,
                                               'sample': sample})
        # write files
        counts.to_csv(output.counts, index=False)
        fates.to_csv(output.fates, index=False)  
        
        
        

rule bind_expr_filters_Wuhan_Hu_1:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Wuhan_Hu_1.ipynb')
    params:
        nb='bind_expr_filters_Wuhan_Hu_1.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule bind_expr_filters_Omicron_BA2:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_BA2.ipynb')
    params:
        nb='bind_expr_filters_Omicron_BA2.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule bind_expr_filters_Omicron_BQ11:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_BQ11.ipynb')
    params:
        nb='bind_expr_filters_Omicron_BQ11.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule bind_expr_filters_Omicron_XBB15:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_XBB15.ipynb')
    params:
        nb='bind_expr_filters_Omicron_XBB15.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule bind_expr_filters_Omicron_EG5:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_EG5.ipynb')
    params:
        nb='bind_expr_filters_Omicron_EG5.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
        
        
rule bind_expr_filters_Omicron_FLip:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_FLip.ipynb')
    params:
        nb='bind_expr_filters_Omicron_FLip.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"

        
rule bind_expr_filters_Omicron_BA286:
    """QC checks on bind & expression filters from DMS data.
    """
    input:
    	file=config['mut_bind_expr']
    output:
        nb_markdown=nb_markdown('bind_expr_filters_Omicron_BA286.ipynb')
    params:
        nb='bind_expr_filters_Omicron_BA286.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"   
        
        
        
        


rule get_mut_bind_expr:
    """Download SARS-CoV-2 mutation ACE2-binding and expression from URL for Wuhan-Hu-1, BA2, BQ11, and XBB15 site-saturation mutagenesis library."""
    output:
        file=config['mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['mut_bind_expr_url'], output.file)
        
rule get_early2020_mut_bind_expr:
    """Download SARS-CoV-2 Wuhan-1 mutation ACE2-binding and expression from URL."""
    output:
        file=config['early2020_mut_bind_expr']
    run:
        urllib.request.urlretrieve(config['early2020_mut_bind_expr_url'], output.file)
