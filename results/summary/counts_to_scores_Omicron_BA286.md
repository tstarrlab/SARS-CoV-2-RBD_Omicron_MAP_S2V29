# Analyze counts and compute escape scores
This Python Jupyter notebook analyzes the variant counts and looks at mutation coverage and jackpotting.
It then computes an "escape scores" for each variant after grouping by barcode or substitutions as specified in the configuration.

## Set up analysis

This notebook primarily makes use of the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package, and uses [plotnine](https://github.com/has2k1/plotnine) for ggplot2-like plotting syntax:


```python
import collections
import math
import os
import warnings

import Bio.SeqIO

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import numpy

import pandas as pd

from plotnine import *

import seaborn

import yaml

%matplotlib inline
```

Set [plotnine](https://github.com/has2k1/plotnine) theme to the gray-grid one defined in [dms_variants](https://jbloomlab.github.io/dms_variants):


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 1.4.3


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['escape_scores_dir'], exist_ok=True)
```

Read information about the samples:


```python
samples_df = pd.read_csv(config['barcode_runs_Omicron_BA286'])
```

## Initialize codon-variant table
Initialize [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) from wildtype gene sequence and the variant counts CSV file.
We will then use the plotting functions of this variant table to analyze the counts per sample:


```python
wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence_Omicron_BA286'], 'fasta')
geneseq = str(wt_seqrecord.seq)
primary_target = wt_seqrecord.name
print(f"Read sequence of {len(geneseq)} nt for {primary_target} from {config['wildtype_sequence_Omicron_BA286']}")
      
print(f"Initializing CodonVariantTable from gene sequence and {config['variant_counts_Omicron_BA286']}")
      
variants = dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df(
                geneseq=geneseq,
                variant_count_df_file=config['variant_counts_Omicron_BA286'],
                primary_target=primary_target,
                allowgaps=True)
      
print('Done initializing CodonVariantTable.')
```

    Read sequence of 600 nt for Omicron_BA286 from data/wildtype_sequence_Omicron_BA286.fasta
    Initializing CodonVariantTable from gene sequence and results/counts/Omicron_BA286/variant_counts.csv.gz
    Done initializing CodonVariantTable.


## Sequencing counts per sample
Average counts per variant for each sample.
Note that these are the **sequencing** counts, in some cases they may outstrip the actual number of sorted cells:


```python
p = variants.plotAvgCountsPerVariant(libraries=variants.libraries,
                                     by_target=False,
                                     orientation='v')
p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_19_0.png)
    


And the numerical values plotted above:


```python
display(HTML(
 variants.avgCountsPerVariant(libraries=variants.libraries,
                               by_target=False)
 .pivot_table(index='sample',
              columns='library',
              values='avg_counts_per_variant')
 .round(1)
 .to_html()
 ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>lib92</th>
      <th>lib93</th>
    </tr>
    <tr>
      <th>sample</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>exptREF-none-0-ref</th>
      <td>740.2</td>
      <td>771.0</td>
    </tr>
    <tr>
      <th>expt7-S2V29_v37_2-48-abneg</th>
      <td>24.7</td>
      <td>17.6</td>
    </tr>
    <tr>
      <th>expt8-S2V29-36-abneg</th>
      <td>21.8</td>
      <td>17.4</td>
    </tr>
  </tbody>
</table>


## Mutations per variant
Average number of mutations per gene among all variants of the primary target, separately for each date:


```python
#this plotting is very slow when lots of samples, so for now plots are commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotNumCodonMutsByType(variant_type='all',
                                       orientation='v',
                                       libraries=variants.libraries,
                                       samples=date_df['sample'].unique().tolist(),
                                       widthscale=2)
   p = p + theme(panel_grid_major_x=element_blank())  # no vertical grid lines
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_23_0.png)
    


Now similar plots but showing mutation frequency across the gene:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for date, date_df in samples_df.groupby('date', sort=False):
   p = variants.plotMutFreqs(variant_type='all',
                             mut_type='codon',
                             orientation='v',
                             libraries=variants.libraries,
                             samples=date_df['sample'].unique().tolist(),
                             widthscale=1.5)
   fig = p.draw()
   display(fig)
   plt.close(fig)
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_25_0.png)
    


## Jackpotting and mutation coverage in pre-selection libraries
We look at the distribution of counts in the "reference" (pre-selection) libraries to see if they seem jackpotted (a few variants at very high frequency):


```python
pre_samples_df = samples_df.query('selection == "reference"')
```

Distribution of mutations along the gene for the pre-selection samples; big spikes may indicate jackpotting:


```python
# this plotting is very slow when lots of samples, so for now code commented out

p = variants.plotMutFreqs(variant_type='all',
                         mut_type='codon',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.5)
_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_29_0.png)
    


How many mutations are observed frequently in pre-selection libraries?
Note that the libraries have been pre-selected for ACE2 binding, so we expect stop variants to mostly be missing.
Make the plot both for all variants and just single-mutant variants:


```python
# this plotting is very slow when lots of samples, so for now code commented out

for variant_type in ['all', 'single']:
   p = variants.plotCumulMutCoverage(
                         variant_type=variant_type,
                         mut_type='aa',
                         orientation='v',
                         libraries=variants.libraries,
                         samples=pre_samples_df['sample'].unique().tolist(),
                         widthscale=1.8,
                         heightscale=1.2)
   _ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_31_0.png)
    



    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_31_1.png)
    


Now make a plot showing what number and fraction of counts are for each variant in each pre-selection sample / library.
If some variants constitute a very high fraction, then that indicates extensive bottlenecking:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_33_0.png)
    



    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_33_1.png)
    


Now make the same plot breaking down by variant class, which enables determination of which types of variants are at high and low frequencies.
For this plot (unlike one above not classified by category) we only show variants of the primary target (not the homologs), and also group synonymous with wildtype in order to reduce number of plotted categories to make more interpretable:


```python
for ystat in ['frac_counts', 'count']:
    p = variants.plotCountsPerVariant(ystat=ystat,
                                      libraries=variants.libraries,
                                      samples=pre_samples_df['sample'].unique().tolist(),
                                      orientation='v',
                                      widthscale=1.75,
                                      by_variant_class=True,
                                      classifyVariants_kwargs={'syn_as_wt': True},
                                      primary_target_only=True,
                                      )
    _ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_35_0.png)
    



    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_35_1.png)
    


We also directly look to see what is the variant in each reference library / sample with the highest fraction counts.
Knowing if the highest frequency variant is shared helps determine **where** in the experiment the jackpotting happened:


```python
frac_counts_per_variant = (
        variants.add_frac_counts(variants.variant_count_df)
        .query(f"sample in {pre_samples_df['sample'].unique().tolist()}")
        )

display(HTML(
    frac_counts_per_variant
    .sort_values('frac_counts', ascending=False)
    .groupby(['library', 'sample'])
    .head(n=1)
    .sort_values(['library', 'sample'])
    .set_index(['library', 'sample'])
    [['frac_counts', 'target', 'barcode', 'aa_substitutions', 'codon_substitutions']]
    .round(4)
    .to_html()
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>frac_counts</th>
      <th>target</th>
      <th>barcode</th>
      <th>aa_substitutions</th>
      <th>codon_substitutions</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>lib92</th>
      <th>exptREF-none-0-ref</th>
      <td>0.0003</td>
      <td>Omicron_BA286</td>
      <td>AACCGTTAAATGCAAA</td>
      <td>K73H</td>
      <td>AAA73CAT</td>
    </tr>
    <tr>
      <th>lib93</th>
      <th>exptREF-none-0-ref</th>
      <td>0.0002</td>
      <td>Omicron_BA286</td>
      <td>CCAATGCGATACGGGT</td>
      <td>T100A</td>
      <td>ACC100GCT</td>
    </tr>
  </tbody>
</table>


To further where the jackpotting relative to the generation of the reference samples, we plot the correlation among the fraction of counts for the different reference samples.
If the fractions are highly correlated, that indicates that the jackpotting occurred in some upstream step common to the reference samples:


```python
# this code makes a full matrix of scatter plots, but is REALLY SLOW. So for now,
# it is commented out in favor of code that just makes correlation matrix.
for lib, lib_df in frac_counts_per_variant.groupby('library'):
   wide_lib_df = lib_df.pivot_table(index=['target', 'barcode'],
                                    columns='sample',
                                    values='frac_counts')
   g = seaborn.pairplot(wide_lib_df, corner=True, plot_kws={'alpha': 0.5}, diag_kind='kde')
   _ = g.fig.suptitle(lib, size=18)
   plt.show()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_39_0.png)
    



    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_39_1.png)
    


## Examine counts for wildtype variants
The type of score we use to quantify escape depends on how well represented wildtype is in the selected libraries.
If wildtype is still well represented, we can use a more conventional functional score that gives differential selection relative to wildtype.
If wildtype is not well represented, then we need an alternative score that does not involve normalizing frequencies to wildtype.

First get average fraction of counts per variant for each variant class:


```python
counts_by_class = (
    variants.variant_count_df
    .pipe(variants.add_frac_counts)
    .pipe(variants.classifyVariants,
          primary_target=variants.primary_target,
          non_primary_target_class='homolog',
          class_as_categorical=True)
    .groupby(['library', 'sample', 'variant_class'])
    .aggregate(avg_frac_counts=pd.NamedAgg('frac_counts', 'mean'))
    .reset_index()
    .merge(samples_df[['sample', 'library', 'date', 'antibody', 'concentration', 'selection']],
           on=['sample', 'library'], validate='many_to_one')
    )

display(HTML(counts_by_class.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>sample</th>
      <th>variant_class</th>
      <th>avg_frac_counts</th>
      <th>date</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>selection</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>wildtype</td>
      <td>0.000012</td>
      <td>240305</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>synonymous</td>
      <td>0.000014</td>
      <td>240305</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>1 nonsynonymous</td>
      <td>0.000012</td>
      <td>240305</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>&gt;1 nonsynonymous</td>
      <td>0.000012</td>
      <td>240305</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
    <tr>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>stop</td>
      <td>0.000014</td>
      <td>240305</td>
      <td>none</td>
      <td>0</td>
      <td>reference</td>
    </tr>
  </tbody>
</table>


Plot average fraction of all counts per variant for each variant class.
If the values for wildtype are low for the non-reference samples (such as more similar to stop the nonsynonymous), then normalizing by wildtype in calculating scores will probably not work well as wildtype is too depleted:


```python
min_frac = 1e-7  # plot values < this as this

p = (ggplot(counts_by_class
            .assign(avg_frac_counts=lambda x: numpy.clip(x['avg_frac_counts'], min_frac, None))
            ) +
     aes('avg_frac_counts', 'sample', color='selection') +
     geom_point(size=2) +
     scale_color_manual(values=CBPALETTE[1:]) +
     facet_grid('library ~ variant_class') +
     scale_x_log10() +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.5 * counts_by_class['variant_class'].nunique(),
                        0.2 * counts_by_class['library'].nunique() * 
                        counts_by_class['sample'].nunique())
           ) +
     geom_vline(xintercept=min_frac, linetype='dotted', color=CBPALETTE[3])
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_43_0.png)
    


## Compute escape scores
We use the escape score metric, which does **not** involve normalizing to wildtype and so isn't strongly affected by low wildtype counts.
We compute the scores using the method [dms_variants.codonvarianttable.CodonVariantTable.escape_scores](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores).

First, define what samples to compare for each calculation, matching each post-selection (escape) to the pre-selection (reference) sample on the same date:


```python
score_sample_df = (
    samples_df
    .query('selection == "escape"')
    .rename(columns={'sample': 'post_sample',
                     'number_cells': 'pre_cells_sorted'})
    .merge(samples_df
           .query('selection == "reference"')
           [['sample', 'library', 'date', 'number_cells']]
           .rename(columns={'sample': 'pre_sample',
                            'number_cells': 'post_cells_sorted'}),
           how='left', on=['date', 'library'], validate='many_to_one',
           )
    .assign(name=lambda x: x['antibody'] + '_' + x['concentration'].astype(str))
    # add dates to names where needed to make unique
    .assign(n_libs=lambda x: x.groupby(['name', 'date'])['pre_sample'].transform('count'))
    .sort_values(['name', 'date', 'n_libs'], ascending=False)
    .assign(i_name=lambda x: x.groupby(['library', 'name'], sort=False)['pre_sample'].cumcount(),
            name=lambda x: x.apply(lambda r: r['name'] + '_' + str(r['date']) if r['i_name'] > 0 else r['name'],
                                   axis=1),
            )
    .sort_values(['antibody', 'concentration', 'library', 'i_name'])
    # get columns of interest
    [['name', 'library', 'antibody', 'concentration', 'date',
      'pre_sample', 'post_sample', 'frac_escape', 'pre_cells_sorted', 'post_cells_sorted']]
    )

assert len(score_sample_df.groupby(['name', 'library'])) == len(score_sample_df)

print(f"Writing samples used to compute escape scores to {config['escape_score_samples_Omicron_BA286']}\n")
score_sample_df.to_csv(config['escape_score_samples_Omicron_BA286'], index=False)

display(HTML(score_sample_df.to_html(index=False)))
```

    Writing samples used to compute escape scores to results/escape_scores/samples_Omicron_BA286.csv
    



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>date</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>frac_escape</th>
      <th>pre_cells_sorted</th>
      <th>post_cells_sorted</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>S2V29</td>
      <td>36</td>
      <td>240305</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>0.101</td>
      <td>404000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>lib93</td>
      <td>S2V29</td>
      <td>36</td>
      <td>240305</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>0.086</td>
      <td>344000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>S2V29_v37_2_48</td>
      <td>lib92</td>
      <td>S2V29_v37_2</td>
      <td>48</td>
      <td>240305</td>
      <td>exptREF-none-0-ref</td>
      <td>expt7-S2V29_v37_2-48-abneg</td>
      <td>0.093</td>
      <td>372000.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>S2V29_v37_2_48</td>
      <td>lib93</td>
      <td>S2V29_v37_2</td>
      <td>48</td>
      <td>240305</td>
      <td>exptREF-none-0-ref</td>
      <td>expt7-S2V29_v37_2-48-abneg</td>
      <td>0.081</td>
      <td>324000.0</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>


Compute the escape scores for variants of the primary target and classify the variants:


```python
print(f"Computing escape scores for {primary_target} variants using {config['escape_score_type']} "
      f"score type with a pseudocount of {config['escape_score_pseudocount']} and "
      f"an escape fraction floor {config['escape_score_floor_E']}, an escape fraction ceiling "
      f"{config['escape_score_ceil_E']}, and grouping variants by {config['escape_score_group_by']}.")

escape_scores = (variants.escape_scores(score_sample_df,
                                        score_type=config['escape_score_type'],
                                        pseudocount=config['escape_score_pseudocount'],
                                        floor_E=config['escape_score_floor_E'],
                                        ceil_E=config['escape_score_ceil_E'],
                                        by=config['escape_score_group_by'],
                                        )
                 .query('target == @primary_target')
                 .pipe(variants.classifyVariants,
                       primary_target=variants.primary_target,
                       syn_as_wt=(config['escape_score_group_by'] == 'aa_substitutions'),
                       )
                 )
print('Here are the first few lines of the resulting escape scores:')
display(HTML(escape_scores.head().to_html(index=False)))
```

    Computing escape scores for Omicron_BA286 variants using frac_escape score type with a pseudocount of 0.5 and an escape fraction floor 0, an escape fraction ceiling 1, and grouping variants by barcode.
    Here are the first few lines of the resulting escape scores:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AACCGTTAAATGCAAA</td>
      <td>0.009930</td>
      <td>1.582215e-06</td>
      <td>21129</td>
      <td>62</td>
      <td>AAA73CAT</td>
      <td>1</td>
      <td>K73H</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AAAAAATCAAACTGAA</td>
      <td>0.549101</td>
      <td>1.248995e-04</td>
      <td>17127</td>
      <td>2801</td>
      <td>TAC177GAT</td>
      <td>1</td>
      <td>Y177D</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CAATATCACCCAAGTG</td>
      <td>0.002073</td>
      <td>4.093863e-07</td>
      <td>17006</td>
      <td>10</td>
      <td>GGT165GAA</td>
      <td>1</td>
      <td>G165E</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AGACGAAATAGATGAC</td>
      <td>0.007114</td>
      <td>1.428507e-06</td>
      <td>16752</td>
      <td>35</td>
      <td>ACC3---</td>
      <td>1</td>
      <td>T3-</td>
      <td>1</td>
      <td>deletion</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CTAACTGACCCACTAA</td>
      <td>0.013908</td>
      <td>3.011216e-06</td>
      <td>15568</td>
      <td>64</td>
      <td>TCC69TTG</td>
      <td>1</td>
      <td>S69L</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
    </tr>
  </tbody>
</table>


## Apply pre-selection count filter to variant escape scores
Now determine a pre-selection count filter in order to flag for removal variants with counts that are so low that the estimated score is probably noise.
We know that stop codons should be largely purged pre-selection, and so the counts for them are a good indication of the "noise" threshold.
We therefore set the filter using the number of pre-selection counts for the stop codons.

To do this, we first compute the number of pre-selection counts for stop-codon variants at various quantiles and look at these.
We then take the number of pre-selection counts at the specified quantile as the filter cutoff, and filter scores for all variants with pre-selection counts less than this filter cutoff:

NOTE, since we didn't pre-sort this library, this is not the appropriate method. This has historically beena round 100, so we'll just apply a filter of 100.


```python
filter_precount = config['escape_score_precount_filter']

print(f"\nSetting the pre-count filter cutoff to the absolute value of {config['escape_score_precount_filter']}")
```

    
    Setting the pre-count filter cutoff to the absolute value of 100


Apply the filter to the escape scores, so that scores that fail the pre-selection count filter are now marked with `pass_pre_count_filter` of `False`:


```python
escape_scores = (
    escape_scores
    .assign(pass_pre_count_filter=lambda x: x['pre_count'] >= config['escape_score_precount_filter'])
    )
escape_scores.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pass_pre_count_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AACCGTTAAATGCAAA</td>
      <td>0.009930</td>
      <td>1.582215e-06</td>
      <td>21129</td>
      <td>62</td>
      <td>AAA73CAT</td>
      <td>1</td>
      <td>K73H</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AAAAAATCAAACTGAA</td>
      <td>0.549101</td>
      <td>1.248995e-04</td>
      <td>17127</td>
      <td>2801</td>
      <td>TAC177GAT</td>
      <td>1</td>
      <td>Y177D</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CAATATCACCCAAGTG</td>
      <td>0.002073</td>
      <td>4.093863e-07</td>
      <td>17006</td>
      <td>10</td>
      <td>GGT165GAA</td>
      <td>1</td>
      <td>G165E</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AGACGAAATAGATGAC</td>
      <td>0.007114</td>
      <td>1.428507e-06</td>
      <td>16752</td>
      <td>35</td>
      <td>ACC3---</td>
      <td>1</td>
      <td>T3-</td>
      <td>1</td>
      <td>deletion</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CTAACTGACCCACTAA</td>
      <td>0.013908</td>
      <td>3.011216e-06</td>
      <td>15568</td>
      <td>64</td>
      <td>TCC69TTG</td>
      <td>1</td>
      <td>S69L</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
</div>



## Apply ACE2-binding / expression filter to variant mutations
We also used deep mutational scanning to estimate how each mutation affected ACE2 binding and expression in the B.1.351 background.
Here we flag for removal any variants of the primary target that have (or have mutations) that were measured to decrease ACE2-binding or expression beyond a minimal threshold, in order to avoid these variants muddying the signal as spurious escape mutants.

To do this, we first determine all mutations that do / do-not having binding that exceeds the thresholds.

Note that because we are working on this serum-mapping project at the same time as we are working on the ACE2-binding / RBD-expression project, the scores will be preliminary until all final analyses have been done on the DMS project end. So, we will allow either preliminary or "final" measurements to be used. 


```python
mut_bind_expr_file = config['mut_bind_expr']
    
print(f"Reading ACE2-binding and expression for mutations from {mut_bind_expr_file}, "
      f"and filtering for variants that have single mutations that "
      f"only have mutations with binding >={config['escape_score_min_bind_mut_Omicron_BA286']} and "
      f"expression >={config['escape_score_min_expr_mut_Omicron_BA286']}.")

mut_bind_expr = (pd.read_csv(mut_bind_expr_file)
                 .query("target == 'Omicron_BA286' and position != 483")  # Remove rows where position == 483
                 .assign(RBD_site=lambda x: numpy.where(x['position'] < 484,
                                                        x['position'] - config['site_number_offset'],
                                                        x['position'] - config['site_number_offset'] - 1),
                         RBD_mutation=lambda x: x['wildtype'] + x['RBD_site'].astype(str) + x['mutant']
                        )
                )

print('Here is what that dataframe looks like:')

display(HTML(mut_bind_expr.query("position > 480 & position < 490").to_html(index=False)))
```

    Reading ACE2-binding and expression for mutations from results/prior_DMS_data/mutant_ACE2binding_expression.csv, and filtering for variants that have single mutations that only have mutations with binding >=-3.0 and expression >=-0.75.
    Here is what that dataframe looks like:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>wildtype</th>
      <th>position</th>
      <th>mutant</th>
      <th>mutation</th>
      <th>bind</th>
      <th>delta_bind</th>
      <th>n_bc_bind</th>
      <th>n_libs_bind</th>
      <th>bind_rep1</th>
      <th>bind_rep2</th>
      <th>bind_rep3</th>
      <th>expr</th>
      <th>delta_expr</th>
      <th>n_bc_expr</th>
      <th>n_libs_expr</th>
      <th>expr_rep1</th>
      <th>expr_rep2</th>
      <th>RBD_site</th>
      <th>RBD_mutation</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>-</td>
      <td>K481-</td>
      <td>8.38</td>
      <td>-1.69</td>
      <td>40</td>
      <td>2</td>
      <td>8.70</td>
      <td>8.06</td>
      <td>NaN</td>
      <td>6.98</td>
      <td>-0.61</td>
      <td>45</td>
      <td>2</td>
      <td>6.96</td>
      <td>7.01</td>
      <td>151</td>
      <td>K151-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>A</td>
      <td>K481A</td>
      <td>10.08</td>
      <td>0.01</td>
      <td>26</td>
      <td>2</td>
      <td>10.31</td>
      <td>9.85</td>
      <td>NaN</td>
      <td>7.79</td>
      <td>0.19</td>
      <td>35</td>
      <td>2</td>
      <td>7.69</td>
      <td>7.88</td>
      <td>151</td>
      <td>K151A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>C</td>
      <td>K481C</td>
      <td>9.66</td>
      <td>-0.41</td>
      <td>25</td>
      <td>2</td>
      <td>9.90</td>
      <td>9.41</td>
      <td>NaN</td>
      <td>7.36</td>
      <td>-0.24</td>
      <td>28</td>
      <td>2</td>
      <td>7.34</td>
      <td>7.38</td>
      <td>151</td>
      <td>K151C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>D</td>
      <td>K481D</td>
      <td>10.12</td>
      <td>0.05</td>
      <td>26</td>
      <td>2</td>
      <td>10.29</td>
      <td>9.94</td>
      <td>NaN</td>
      <td>7.86</td>
      <td>0.27</td>
      <td>37</td>
      <td>2</td>
      <td>7.74</td>
      <td>7.99</td>
      <td>151</td>
      <td>K151D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>E</td>
      <td>K481E</td>
      <td>10.02</td>
      <td>-0.05</td>
      <td>19</td>
      <td>2</td>
      <td>10.31</td>
      <td>9.74</td>
      <td>NaN</td>
      <td>7.87</td>
      <td>0.28</td>
      <td>31</td>
      <td>2</td>
      <td>7.77</td>
      <td>7.97</td>
      <td>151</td>
      <td>K151E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>F</td>
      <td>K481F</td>
      <td>10.06</td>
      <td>-0.01</td>
      <td>17</td>
      <td>2</td>
      <td>10.28</td>
      <td>9.84</td>
      <td>NaN</td>
      <td>7.64</td>
      <td>0.05</td>
      <td>19</td>
      <td>2</td>
      <td>7.67</td>
      <td>7.62</td>
      <td>151</td>
      <td>K151F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>G</td>
      <td>K481G</td>
      <td>9.89</td>
      <td>-0.18</td>
      <td>48</td>
      <td>2</td>
      <td>10.14</td>
      <td>9.63</td>
      <td>NaN</td>
      <td>7.56</td>
      <td>-0.03</td>
      <td>55</td>
      <td>2</td>
      <td>7.51</td>
      <td>7.61</td>
      <td>151</td>
      <td>K151G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>H</td>
      <td>K481H</td>
      <td>10.04</td>
      <td>-0.03</td>
      <td>22</td>
      <td>2</td>
      <td>10.31</td>
      <td>9.77</td>
      <td>NaN</td>
      <td>7.70</td>
      <td>0.11</td>
      <td>25</td>
      <td>2</td>
      <td>7.59</td>
      <td>7.81</td>
      <td>151</td>
      <td>K151H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>I</td>
      <td>K481I</td>
      <td>10.10</td>
      <td>0.03</td>
      <td>11</td>
      <td>2</td>
      <td>10.19</td>
      <td>10.00</td>
      <td>NaN</td>
      <td>7.83</td>
      <td>0.23</td>
      <td>17</td>
      <td>2</td>
      <td>7.75</td>
      <td>7.90</td>
      <td>151</td>
      <td>K151I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>K</td>
      <td>K481K</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>151</td>
      <td>K151K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>L</td>
      <td>K481L</td>
      <td>10.05</td>
      <td>-0.02</td>
      <td>30</td>
      <td>2</td>
      <td>10.26</td>
      <td>9.84</td>
      <td>NaN</td>
      <td>7.78</td>
      <td>0.18</td>
      <td>37</td>
      <td>2</td>
      <td>7.73</td>
      <td>7.82</td>
      <td>151</td>
      <td>K151L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>M</td>
      <td>K481M</td>
      <td>10.06</td>
      <td>-0.01</td>
      <td>27</td>
      <td>2</td>
      <td>10.29</td>
      <td>9.84</td>
      <td>NaN</td>
      <td>7.77</td>
      <td>0.17</td>
      <td>35</td>
      <td>2</td>
      <td>7.68</td>
      <td>7.85</td>
      <td>151</td>
      <td>K151M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>N</td>
      <td>K481N</td>
      <td>10.08</td>
      <td>0.01</td>
      <td>26</td>
      <td>2</td>
      <td>10.41</td>
      <td>9.75</td>
      <td>NaN</td>
      <td>7.87</td>
      <td>0.27</td>
      <td>35</td>
      <td>2</td>
      <td>7.78</td>
      <td>7.95</td>
      <td>151</td>
      <td>K151N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>P</td>
      <td>K481P</td>
      <td>9.91</td>
      <td>-0.16</td>
      <td>19</td>
      <td>2</td>
      <td>10.24</td>
      <td>9.57</td>
      <td>NaN</td>
      <td>7.73</td>
      <td>0.13</td>
      <td>23</td>
      <td>2</td>
      <td>7.77</td>
      <td>7.68</td>
      <td>151</td>
      <td>K151P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>Q</td>
      <td>K481Q</td>
      <td>10.23</td>
      <td>0.16</td>
      <td>27</td>
      <td>2</td>
      <td>10.36</td>
      <td>10.09</td>
      <td>NaN</td>
      <td>7.75</td>
      <td>0.16</td>
      <td>33</td>
      <td>2</td>
      <td>7.52</td>
      <td>7.99</td>
      <td>151</td>
      <td>K151Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>R</td>
      <td>K481R</td>
      <td>9.92</td>
      <td>-0.15</td>
      <td>28</td>
      <td>2</td>
      <td>10.36</td>
      <td>9.48</td>
      <td>NaN</td>
      <td>7.55</td>
      <td>-0.04</td>
      <td>34</td>
      <td>2</td>
      <td>7.55</td>
      <td>7.54</td>
      <td>151</td>
      <td>K151R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>S</td>
      <td>K481S</td>
      <td>10.05</td>
      <td>-0.02</td>
      <td>20</td>
      <td>2</td>
      <td>10.12</td>
      <td>9.98</td>
      <td>NaN</td>
      <td>7.84</td>
      <td>0.25</td>
      <td>31</td>
      <td>2</td>
      <td>7.78</td>
      <td>7.90</td>
      <td>151</td>
      <td>K151S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>T</td>
      <td>K481T</td>
      <td>10.13</td>
      <td>0.06</td>
      <td>20</td>
      <td>2</td>
      <td>10.38</td>
      <td>9.88</td>
      <td>NaN</td>
      <td>8.04</td>
      <td>0.44</td>
      <td>25</td>
      <td>2</td>
      <td>7.84</td>
      <td>8.24</td>
      <td>151</td>
      <td>K151T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>V</td>
      <td>K481V</td>
      <td>10.13</td>
      <td>0.06</td>
      <td>27</td>
      <td>2</td>
      <td>10.39</td>
      <td>9.87</td>
      <td>NaN</td>
      <td>7.82</td>
      <td>0.22</td>
      <td>33</td>
      <td>2</td>
      <td>7.65</td>
      <td>7.98</td>
      <td>151</td>
      <td>K151V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>W</td>
      <td>K481W</td>
      <td>9.88</td>
      <td>-0.19</td>
      <td>34</td>
      <td>2</td>
      <td>10.22</td>
      <td>9.54</td>
      <td>NaN</td>
      <td>7.34</td>
      <td>-0.25</td>
      <td>39</td>
      <td>2</td>
      <td>7.27</td>
      <td>7.42</td>
      <td>151</td>
      <td>K151W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>481</td>
      <td>Y</td>
      <td>K481Y</td>
      <td>10.23</td>
      <td>0.15</td>
      <td>17</td>
      <td>2</td>
      <td>10.46</td>
      <td>9.99</td>
      <td>NaN</td>
      <td>7.78</td>
      <td>0.18</td>
      <td>20</td>
      <td>2</td>
      <td>7.72</td>
      <td>7.83</td>
      <td>151</td>
      <td>K151Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>-</td>
      <td>G482-</td>
      <td>7.16</td>
      <td>-2.91</td>
      <td>26</td>
      <td>2</td>
      <td>7.27</td>
      <td>7.04</td>
      <td>NaN</td>
      <td>6.81</td>
      <td>-0.78</td>
      <td>31</td>
      <td>2</td>
      <td>6.69</td>
      <td>6.94</td>
      <td>152</td>
      <td>G152-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>A</td>
      <td>G482A</td>
      <td>10.05</td>
      <td>-0.02</td>
      <td>25</td>
      <td>2</td>
      <td>10.24</td>
      <td>9.86</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>30</td>
      <td>2</td>
      <td>7.52</td>
      <td>7.66</td>
      <td>152</td>
      <td>G152A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>C</td>
      <td>G482C</td>
      <td>9.03</td>
      <td>-1.04</td>
      <td>23</td>
      <td>2</td>
      <td>9.38</td>
      <td>8.68</td>
      <td>NaN</td>
      <td>7.14</td>
      <td>-0.45</td>
      <td>28</td>
      <td>2</td>
      <td>7.08</td>
      <td>7.20</td>
      <td>152</td>
      <td>G152C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>D</td>
      <td>G482D</td>
      <td>10.02</td>
      <td>-0.05</td>
      <td>14</td>
      <td>2</td>
      <td>10.43</td>
      <td>9.62</td>
      <td>NaN</td>
      <td>7.85</td>
      <td>0.25</td>
      <td>23</td>
      <td>2</td>
      <td>7.71</td>
      <td>7.98</td>
      <td>152</td>
      <td>G152D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>E</td>
      <td>G482E</td>
      <td>10.06</td>
      <td>-0.01</td>
      <td>16</td>
      <td>2</td>
      <td>10.29</td>
      <td>9.84</td>
      <td>NaN</td>
      <td>7.68</td>
      <td>0.09</td>
      <td>17</td>
      <td>2</td>
      <td>7.58</td>
      <td>7.78</td>
      <td>152</td>
      <td>G152E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>F</td>
      <td>G482F</td>
      <td>9.78</td>
      <td>-0.29</td>
      <td>21</td>
      <td>2</td>
      <td>10.05</td>
      <td>9.50</td>
      <td>NaN</td>
      <td>7.43</td>
      <td>-0.16</td>
      <td>25</td>
      <td>2</td>
      <td>7.39</td>
      <td>7.48</td>
      <td>152</td>
      <td>G152F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>G</td>
      <td>G482G</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>152</td>
      <td>G152G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>H</td>
      <td>G482H</td>
      <td>9.74</td>
      <td>-0.33</td>
      <td>25</td>
      <td>2</td>
      <td>10.07</td>
      <td>9.41</td>
      <td>NaN</td>
      <td>7.20</td>
      <td>-0.39</td>
      <td>29</td>
      <td>2</td>
      <td>7.00</td>
      <td>7.40</td>
      <td>152</td>
      <td>G152H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>I</td>
      <td>G482I</td>
      <td>9.81</td>
      <td>-0.26</td>
      <td>16</td>
      <td>2</td>
      <td>10.05</td>
      <td>9.58</td>
      <td>NaN</td>
      <td>7.47</td>
      <td>-0.12</td>
      <td>19</td>
      <td>2</td>
      <td>7.36</td>
      <td>7.58</td>
      <td>152</td>
      <td>G152I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>K</td>
      <td>G482K</td>
      <td>9.93</td>
      <td>-0.14</td>
      <td>26</td>
      <td>2</td>
      <td>10.23</td>
      <td>9.64</td>
      <td>NaN</td>
      <td>7.08</td>
      <td>-0.51</td>
      <td>30</td>
      <td>2</td>
      <td>6.86</td>
      <td>7.31</td>
      <td>152</td>
      <td>G152K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>L</td>
      <td>G482L</td>
      <td>9.68</td>
      <td>-0.39</td>
      <td>24</td>
      <td>2</td>
      <td>9.95</td>
      <td>9.40</td>
      <td>NaN</td>
      <td>7.45</td>
      <td>-0.14</td>
      <td>30</td>
      <td>2</td>
      <td>7.27</td>
      <td>7.62</td>
      <td>152</td>
      <td>G152L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>M</td>
      <td>G482M</td>
      <td>9.93</td>
      <td>-0.14</td>
      <td>26</td>
      <td>2</td>
      <td>10.28</td>
      <td>9.58</td>
      <td>NaN</td>
      <td>7.47</td>
      <td>-0.12</td>
      <td>28</td>
      <td>2</td>
      <td>7.34</td>
      <td>7.60</td>
      <td>152</td>
      <td>G152M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>N</td>
      <td>G482N</td>
      <td>9.90</td>
      <td>-0.17</td>
      <td>31</td>
      <td>2</td>
      <td>10.13</td>
      <td>9.67</td>
      <td>NaN</td>
      <td>7.50</td>
      <td>-0.09</td>
      <td>34</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.51</td>
      <td>152</td>
      <td>G152N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>P</td>
      <td>G482P</td>
      <td>8.95</td>
      <td>-1.12</td>
      <td>26</td>
      <td>2</td>
      <td>9.30</td>
      <td>8.61</td>
      <td>NaN</td>
      <td>6.99</td>
      <td>-0.61</td>
      <td>32</td>
      <td>2</td>
      <td>6.90</td>
      <td>7.07</td>
      <td>152</td>
      <td>G152P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>Q</td>
      <td>G482Q</td>
      <td>9.92</td>
      <td>-0.16</td>
      <td>24</td>
      <td>2</td>
      <td>10.14</td>
      <td>9.69</td>
      <td>NaN</td>
      <td>7.65</td>
      <td>0.05</td>
      <td>32</td>
      <td>2</td>
      <td>7.60</td>
      <td>7.69</td>
      <td>152</td>
      <td>G152Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>R</td>
      <td>G482R</td>
      <td>9.77</td>
      <td>-0.30</td>
      <td>30</td>
      <td>2</td>
      <td>10.14</td>
      <td>9.40</td>
      <td>NaN</td>
      <td>7.31</td>
      <td>-0.28</td>
      <td>36</td>
      <td>2</td>
      <td>7.17</td>
      <td>7.44</td>
      <td>152</td>
      <td>G152R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>S</td>
      <td>G482S</td>
      <td>10.00</td>
      <td>-0.07</td>
      <td>28</td>
      <td>2</td>
      <td>10.26</td>
      <td>9.74</td>
      <td>NaN</td>
      <td>7.52</td>
      <td>-0.07</td>
      <td>34</td>
      <td>2</td>
      <td>7.54</td>
      <td>7.50</td>
      <td>152</td>
      <td>G152S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>T</td>
      <td>G482T</td>
      <td>10.10</td>
      <td>0.03</td>
      <td>21</td>
      <td>2</td>
      <td>10.35</td>
      <td>9.86</td>
      <td>NaN</td>
      <td>7.70</td>
      <td>0.11</td>
      <td>24</td>
      <td>2</td>
      <td>7.63</td>
      <td>7.77</td>
      <td>152</td>
      <td>G152T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>V</td>
      <td>G482V</td>
      <td>9.92</td>
      <td>-0.15</td>
      <td>20</td>
      <td>2</td>
      <td>10.19</td>
      <td>9.65</td>
      <td>NaN</td>
      <td>7.69</td>
      <td>0.09</td>
      <td>31</td>
      <td>2</td>
      <td>7.56</td>
      <td>7.81</td>
      <td>152</td>
      <td>G152V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>W</td>
      <td>G482W</td>
      <td>9.67</td>
      <td>-0.40</td>
      <td>30</td>
      <td>2</td>
      <td>9.80</td>
      <td>9.55</td>
      <td>NaN</td>
      <td>7.13</td>
      <td>-0.47</td>
      <td>37</td>
      <td>2</td>
      <td>6.86</td>
      <td>7.39</td>
      <td>152</td>
      <td>G152W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>482</td>
      <td>Y</td>
      <td>G482Y</td>
      <td>9.72</td>
      <td>-0.35</td>
      <td>22</td>
      <td>2</td>
      <td>10.02</td>
      <td>9.42</td>
      <td>NaN</td>
      <td>7.33</td>
      <td>-0.27</td>
      <td>26</td>
      <td>2</td>
      <td>7.16</td>
      <td>7.49</td>
      <td>152</td>
      <td>G152Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>-</td>
      <td>K484-</td>
      <td>7.68</td>
      <td>-2.39</td>
      <td>29</td>
      <td>2</td>
      <td>7.73</td>
      <td>7.63</td>
      <td>NaN</td>
      <td>7.18</td>
      <td>-0.41</td>
      <td>33</td>
      <td>2</td>
      <td>7.06</td>
      <td>7.30</td>
      <td>153</td>
      <td>K153-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>A</td>
      <td>K484A</td>
      <td>9.91</td>
      <td>-0.16</td>
      <td>22</td>
      <td>2</td>
      <td>10.19</td>
      <td>9.63</td>
      <td>NaN</td>
      <td>7.79</td>
      <td>0.20</td>
      <td>28</td>
      <td>2</td>
      <td>7.77</td>
      <td>7.81</td>
      <td>153</td>
      <td>K153A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>C</td>
      <td>K484C</td>
      <td>9.06</td>
      <td>-1.01</td>
      <td>37</td>
      <td>2</td>
      <td>9.54</td>
      <td>8.58</td>
      <td>NaN</td>
      <td>7.48</td>
      <td>-0.12</td>
      <td>42</td>
      <td>2</td>
      <td>7.38</td>
      <td>7.58</td>
      <td>153</td>
      <td>K153C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>D</td>
      <td>K484D</td>
      <td>9.74</td>
      <td>-0.33</td>
      <td>34</td>
      <td>2</td>
      <td>10.06</td>
      <td>9.43</td>
      <td>NaN</td>
      <td>7.70</td>
      <td>0.11</td>
      <td>36</td>
      <td>2</td>
      <td>7.61</td>
      <td>7.80</td>
      <td>153</td>
      <td>K153D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>E</td>
      <td>K484E</td>
      <td>10.03</td>
      <td>-0.05</td>
      <td>15</td>
      <td>2</td>
      <td>10.39</td>
      <td>9.66</td>
      <td>NaN</td>
      <td>8.11</td>
      <td>0.52</td>
      <td>27</td>
      <td>2</td>
      <td>8.02</td>
      <td>8.21</td>
      <td>153</td>
      <td>K153E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>F</td>
      <td>K484F</td>
      <td>8.37</td>
      <td>-1.70</td>
      <td>27</td>
      <td>2</td>
      <td>8.64</td>
      <td>8.09</td>
      <td>NaN</td>
      <td>7.31</td>
      <td>-0.28</td>
      <td>28</td>
      <td>2</td>
      <td>7.21</td>
      <td>7.41</td>
      <td>153</td>
      <td>K153F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>G</td>
      <td>K484G</td>
      <td>9.55</td>
      <td>-0.52</td>
      <td>29</td>
      <td>2</td>
      <td>9.85</td>
      <td>9.26</td>
      <td>NaN</td>
      <td>7.63</td>
      <td>0.04</td>
      <td>35</td>
      <td>2</td>
      <td>7.59</td>
      <td>7.66</td>
      <td>153</td>
      <td>K153G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>H</td>
      <td>K484H</td>
      <td>9.04</td>
      <td>-1.03</td>
      <td>29</td>
      <td>2</td>
      <td>9.47</td>
      <td>8.62</td>
      <td>NaN</td>
      <td>7.37</td>
      <td>-0.22</td>
      <td>29</td>
      <td>2</td>
      <td>7.40</td>
      <td>7.35</td>
      <td>153</td>
      <td>K153H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>I</td>
      <td>K484I</td>
      <td>8.89</td>
      <td>-1.18</td>
      <td>24</td>
      <td>2</td>
      <td>9.24</td>
      <td>8.53</td>
      <td>NaN</td>
      <td>7.58</td>
      <td>-0.01</td>
      <td>28</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.67</td>
      <td>153</td>
      <td>K153I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>K</td>
      <td>K484K</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>153</td>
      <td>K153K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>L</td>
      <td>K484L</td>
      <td>9.09</td>
      <td>-0.98</td>
      <td>30</td>
      <td>2</td>
      <td>9.57</td>
      <td>8.60</td>
      <td>NaN</td>
      <td>7.67</td>
      <td>0.08</td>
      <td>34</td>
      <td>2</td>
      <td>7.35</td>
      <td>7.99</td>
      <td>153</td>
      <td>K153L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>M</td>
      <td>K484M</td>
      <td>9.31</td>
      <td>-0.76</td>
      <td>32</td>
      <td>2</td>
      <td>9.60</td>
      <td>9.01</td>
      <td>NaN</td>
      <td>7.69</td>
      <td>0.10</td>
      <td>34</td>
      <td>2</td>
      <td>7.53</td>
      <td>7.86</td>
      <td>153</td>
      <td>K153M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>N</td>
      <td>K484N</td>
      <td>9.73</td>
      <td>-0.34</td>
      <td>23</td>
      <td>2</td>
      <td>10.01</td>
      <td>9.45</td>
      <td>NaN</td>
      <td>7.85</td>
      <td>0.25</td>
      <td>31</td>
      <td>2</td>
      <td>7.75</td>
      <td>7.95</td>
      <td>153</td>
      <td>K153N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>P</td>
      <td>K484P</td>
      <td>9.36</td>
      <td>-0.72</td>
      <td>29</td>
      <td>2</td>
      <td>9.65</td>
      <td>9.06</td>
      <td>NaN</td>
      <td>7.68</td>
      <td>0.09</td>
      <td>34</td>
      <td>2</td>
      <td>7.45</td>
      <td>7.91</td>
      <td>153</td>
      <td>K153P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>Q</td>
      <td>K484Q</td>
      <td>9.99</td>
      <td>-0.08</td>
      <td>20</td>
      <td>2</td>
      <td>10.29</td>
      <td>9.68</td>
      <td>NaN</td>
      <td>7.92</td>
      <td>0.33</td>
      <td>21</td>
      <td>2</td>
      <td>7.68</td>
      <td>8.17</td>
      <td>153</td>
      <td>K153Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>R</td>
      <td>K484R</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>34</td>
      <td>2</td>
      <td>10.37</td>
      <td>9.77</td>
      <td>NaN</td>
      <td>7.53</td>
      <td>-0.07</td>
      <td>40</td>
      <td>2</td>
      <td>7.50</td>
      <td>7.55</td>
      <td>153</td>
      <td>K153R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>S</td>
      <td>K484S</td>
      <td>9.95</td>
      <td>-0.12</td>
      <td>32</td>
      <td>2</td>
      <td>10.24</td>
      <td>9.66</td>
      <td>NaN</td>
      <td>7.84</td>
      <td>0.25</td>
      <td>41</td>
      <td>2</td>
      <td>7.74</td>
      <td>7.94</td>
      <td>153</td>
      <td>K153S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>T</td>
      <td>K484T</td>
      <td>9.87</td>
      <td>-0.20</td>
      <td>36</td>
      <td>2</td>
      <td>10.17</td>
      <td>9.56</td>
      <td>NaN</td>
      <td>7.93</td>
      <td>0.34</td>
      <td>41</td>
      <td>2</td>
      <td>7.78</td>
      <td>8.08</td>
      <td>153</td>
      <td>K153T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>V</td>
      <td>K484V</td>
      <td>9.04</td>
      <td>-1.03</td>
      <td>33</td>
      <td>2</td>
      <td>9.37</td>
      <td>8.71</td>
      <td>NaN</td>
      <td>7.70</td>
      <td>0.11</td>
      <td>40</td>
      <td>2</td>
      <td>7.47</td>
      <td>7.93</td>
      <td>153</td>
      <td>K153V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>W</td>
      <td>K484W</td>
      <td>8.86</td>
      <td>-1.21</td>
      <td>32</td>
      <td>2</td>
      <td>9.28</td>
      <td>8.44</td>
      <td>NaN</td>
      <td>7.30</td>
      <td>-0.29</td>
      <td>34</td>
      <td>2</td>
      <td>7.16</td>
      <td>7.43</td>
      <td>153</td>
      <td>K153W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>K</td>
      <td>484</td>
      <td>Y</td>
      <td>K484Y</td>
      <td>8.87</td>
      <td>-1.20</td>
      <td>28</td>
      <td>2</td>
      <td>9.26</td>
      <td>8.47</td>
      <td>NaN</td>
      <td>7.37</td>
      <td>-0.22</td>
      <td>28</td>
      <td>2</td>
      <td>7.28</td>
      <td>7.47</td>
      <td>153</td>
      <td>K153Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>-</td>
      <td>G485-</td>
      <td>6.95</td>
      <td>-3.12</td>
      <td>34</td>
      <td>2</td>
      <td>7.31</td>
      <td>6.60</td>
      <td>NaN</td>
      <td>6.96</td>
      <td>-0.64</td>
      <td>39</td>
      <td>2</td>
      <td>6.81</td>
      <td>7.11</td>
      <td>154</td>
      <td>G154-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>A</td>
      <td>G485A</td>
      <td>9.58</td>
      <td>-0.49</td>
      <td>30</td>
      <td>2</td>
      <td>9.84</td>
      <td>9.32</td>
      <td>NaN</td>
      <td>7.26</td>
      <td>-0.33</td>
      <td>33</td>
      <td>2</td>
      <td>7.13</td>
      <td>7.40</td>
      <td>154</td>
      <td>G154A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>C</td>
      <td>G485C</td>
      <td>7.46</td>
      <td>-2.61</td>
      <td>43</td>
      <td>2</td>
      <td>7.41</td>
      <td>7.50</td>
      <td>NaN</td>
      <td>6.62</td>
      <td>-0.97</td>
      <td>45</td>
      <td>2</td>
      <td>6.47</td>
      <td>6.77</td>
      <td>154</td>
      <td>G154C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>D</td>
      <td>G485D</td>
      <td>9.10</td>
      <td>-0.97</td>
      <td>29</td>
      <td>2</td>
      <td>9.53</td>
      <td>8.68</td>
      <td>NaN</td>
      <td>7.10</td>
      <td>-0.49</td>
      <td>37</td>
      <td>2</td>
      <td>7.06</td>
      <td>7.15</td>
      <td>154</td>
      <td>G154D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>E</td>
      <td>G485E</td>
      <td>9.72</td>
      <td>-0.35</td>
      <td>24</td>
      <td>2</td>
      <td>9.96</td>
      <td>9.48</td>
      <td>NaN</td>
      <td>7.17</td>
      <td>-0.42</td>
      <td>27</td>
      <td>2</td>
      <td>7.07</td>
      <td>7.27</td>
      <td>154</td>
      <td>G154E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>F</td>
      <td>G485F</td>
      <td>8.55</td>
      <td>-1.52</td>
      <td>24</td>
      <td>2</td>
      <td>8.80</td>
      <td>8.30</td>
      <td>NaN</td>
      <td>6.91</td>
      <td>-0.69</td>
      <td>25</td>
      <td>2</td>
      <td>6.74</td>
      <td>7.08</td>
      <td>154</td>
      <td>G154F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>G</td>
      <td>G485G</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>154</td>
      <td>G154G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>H</td>
      <td>G485H</td>
      <td>8.62</td>
      <td>-1.45</td>
      <td>28</td>
      <td>2</td>
      <td>8.91</td>
      <td>8.33</td>
      <td>NaN</td>
      <td>6.72</td>
      <td>-0.88</td>
      <td>30</td>
      <td>2</td>
      <td>6.48</td>
      <td>6.95</td>
      <td>154</td>
      <td>G154H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>I</td>
      <td>G485I</td>
      <td>9.01</td>
      <td>-1.06</td>
      <td>29</td>
      <td>2</td>
      <td>9.29</td>
      <td>8.74</td>
      <td>NaN</td>
      <td>6.92</td>
      <td>-0.67</td>
      <td>34</td>
      <td>2</td>
      <td>6.76</td>
      <td>7.08</td>
      <td>154</td>
      <td>G154I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>K</td>
      <td>G485K</td>
      <td>8.76</td>
      <td>-1.31</td>
      <td>34</td>
      <td>2</td>
      <td>9.04</td>
      <td>8.48</td>
      <td>NaN</td>
      <td>6.72</td>
      <td>-0.87</td>
      <td>33</td>
      <td>2</td>
      <td>6.63</td>
      <td>6.82</td>
      <td>154</td>
      <td>G154K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>L</td>
      <td>G485L</td>
      <td>8.81</td>
      <td>-1.26</td>
      <td>27</td>
      <td>2</td>
      <td>9.10</td>
      <td>8.52</td>
      <td>NaN</td>
      <td>6.78</td>
      <td>-0.81</td>
      <td>28</td>
      <td>2</td>
      <td>6.65</td>
      <td>6.91</td>
      <td>154</td>
      <td>G154L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>M</td>
      <td>G485M</td>
      <td>9.59</td>
      <td>-0.48</td>
      <td>22</td>
      <td>2</td>
      <td>9.86</td>
      <td>9.32</td>
      <td>NaN</td>
      <td>7.01</td>
      <td>-0.59</td>
      <td>24</td>
      <td>2</td>
      <td>6.94</td>
      <td>7.07</td>
      <td>154</td>
      <td>G154M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>N</td>
      <td>G485N</td>
      <td>8.73</td>
      <td>-1.34</td>
      <td>27</td>
      <td>2</td>
      <td>9.02</td>
      <td>8.44</td>
      <td>NaN</td>
      <td>7.00</td>
      <td>-0.59</td>
      <td>28</td>
      <td>2</td>
      <td>6.94</td>
      <td>7.07</td>
      <td>154</td>
      <td>G154N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>P</td>
      <td>G485P</td>
      <td>7.87</td>
      <td>-2.20</td>
      <td>33</td>
      <td>2</td>
      <td>8.22</td>
      <td>7.53</td>
      <td>NaN</td>
      <td>6.80</td>
      <td>-0.80</td>
      <td>33</td>
      <td>2</td>
      <td>6.67</td>
      <td>6.93</td>
      <td>154</td>
      <td>G154P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>Q</td>
      <td>G485Q</td>
      <td>9.44</td>
      <td>-0.63</td>
      <td>36</td>
      <td>2</td>
      <td>9.71</td>
      <td>9.16</td>
      <td>NaN</td>
      <td>7.06</td>
      <td>-0.54</td>
      <td>40</td>
      <td>2</td>
      <td>6.96</td>
      <td>7.15</td>
      <td>154</td>
      <td>G154Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>R</td>
      <td>G485R</td>
      <td>8.80</td>
      <td>-1.27</td>
      <td>29</td>
      <td>2</td>
      <td>9.06</td>
      <td>8.54</td>
      <td>NaN</td>
      <td>6.71</td>
      <td>-0.89</td>
      <td>30</td>
      <td>2</td>
      <td>6.58</td>
      <td>6.83</td>
      <td>154</td>
      <td>G154R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>S</td>
      <td>G485S</td>
      <td>9.53</td>
      <td>-0.54</td>
      <td>33</td>
      <td>2</td>
      <td>9.77</td>
      <td>9.29</td>
      <td>NaN</td>
      <td>7.18</td>
      <td>-0.41</td>
      <td>39</td>
      <td>2</td>
      <td>7.06</td>
      <td>7.30</td>
      <td>154</td>
      <td>G154S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>T</td>
      <td>G485T</td>
      <td>9.31</td>
      <td>-0.76</td>
      <td>26</td>
      <td>2</td>
      <td>9.70</td>
      <td>8.93</td>
      <td>NaN</td>
      <td>7.09</td>
      <td>-0.50</td>
      <td>30</td>
      <td>2</td>
      <td>7.00</td>
      <td>7.19</td>
      <td>154</td>
      <td>G154T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>V</td>
      <td>G485V</td>
      <td>9.10</td>
      <td>-0.97</td>
      <td>38</td>
      <td>2</td>
      <td>9.53</td>
      <td>8.67</td>
      <td>NaN</td>
      <td>6.95</td>
      <td>-0.64</td>
      <td>42</td>
      <td>2</td>
      <td>6.86</td>
      <td>7.04</td>
      <td>154</td>
      <td>G154V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>W</td>
      <td>G485W</td>
      <td>8.30</td>
      <td>-1.77</td>
      <td>38</td>
      <td>2</td>
      <td>8.49</td>
      <td>8.12</td>
      <td>NaN</td>
      <td>6.78</td>
      <td>-0.81</td>
      <td>42</td>
      <td>2</td>
      <td>6.69</td>
      <td>6.87</td>
      <td>154</td>
      <td>G154W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>G</td>
      <td>485</td>
      <td>Y</td>
      <td>G485Y</td>
      <td>8.46</td>
      <td>-1.61</td>
      <td>21</td>
      <td>2</td>
      <td>8.56</td>
      <td>8.36</td>
      <td>NaN</td>
      <td>6.64</td>
      <td>-0.96</td>
      <td>25</td>
      <td>2</td>
      <td>6.40</td>
      <td>6.87</td>
      <td>154</td>
      <td>G154Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>-</td>
      <td>P486-</td>
      <td>6.73</td>
      <td>-3.34</td>
      <td>28</td>
      <td>2</td>
      <td>7.02</td>
      <td>6.44</td>
      <td>NaN</td>
      <td>7.04</td>
      <td>-0.56</td>
      <td>29</td>
      <td>2</td>
      <td>6.91</td>
      <td>7.16</td>
      <td>155</td>
      <td>P155-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>A</td>
      <td>P486A</td>
      <td>9.80</td>
      <td>-0.27</td>
      <td>33</td>
      <td>2</td>
      <td>10.10</td>
      <td>9.50</td>
      <td>NaN</td>
      <td>7.37</td>
      <td>-0.23</td>
      <td>37</td>
      <td>2</td>
      <td>7.25</td>
      <td>7.48</td>
      <td>155</td>
      <td>P155A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>C</td>
      <td>P486C</td>
      <td>7.69</td>
      <td>-2.38</td>
      <td>40</td>
      <td>2</td>
      <td>7.85</td>
      <td>7.53</td>
      <td>NaN</td>
      <td>6.78</td>
      <td>-0.81</td>
      <td>45</td>
      <td>2</td>
      <td>6.64</td>
      <td>6.92</td>
      <td>155</td>
      <td>P155C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>D</td>
      <td>P486D</td>
      <td>9.04</td>
      <td>-1.03</td>
      <td>44</td>
      <td>2</td>
      <td>9.38</td>
      <td>8.71</td>
      <td>NaN</td>
      <td>7.50</td>
      <td>-0.09</td>
      <td>48</td>
      <td>2</td>
      <td>7.34</td>
      <td>7.66</td>
      <td>155</td>
      <td>P155D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>E</td>
      <td>P486E</td>
      <td>9.05</td>
      <td>-1.02</td>
      <td>28</td>
      <td>2</td>
      <td>9.39</td>
      <td>8.70</td>
      <td>NaN</td>
      <td>7.76</td>
      <td>0.17</td>
      <td>29</td>
      <td>2</td>
      <td>7.75</td>
      <td>7.78</td>
      <td>155</td>
      <td>P155E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>F</td>
      <td>P486F</td>
      <td>10.36</td>
      <td>0.29</td>
      <td>25</td>
      <td>2</td>
      <td>10.47</td>
      <td>10.26</td>
      <td>NaN</td>
      <td>7.43</td>
      <td>-0.16</td>
      <td>41</td>
      <td>2</td>
      <td>7.36</td>
      <td>7.50</td>
      <td>155</td>
      <td>P155F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>G</td>
      <td>P486G</td>
      <td>8.85</td>
      <td>-1.22</td>
      <td>44</td>
      <td>2</td>
      <td>9.17</td>
      <td>8.53</td>
      <td>NaN</td>
      <td>7.10</td>
      <td>-0.49</td>
      <td>47</td>
      <td>2</td>
      <td>7.00</td>
      <td>7.20</td>
      <td>155</td>
      <td>P155G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>H</td>
      <td>P486H</td>
      <td>9.61</td>
      <td>-0.46</td>
      <td>19</td>
      <td>2</td>
      <td>9.87</td>
      <td>9.34</td>
      <td>NaN</td>
      <td>6.91</td>
      <td>-0.68</td>
      <td>23</td>
      <td>2</td>
      <td>6.85</td>
      <td>6.98</td>
      <td>155</td>
      <td>P155H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>I</td>
      <td>P486I</td>
      <td>9.91</td>
      <td>-0.16</td>
      <td>19</td>
      <td>2</td>
      <td>10.17</td>
      <td>9.65</td>
      <td>NaN</td>
      <td>7.44</td>
      <td>-0.16</td>
      <td>29</td>
      <td>2</td>
      <td>7.39</td>
      <td>7.48</td>
      <td>155</td>
      <td>P155I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>K</td>
      <td>P486K</td>
      <td>9.26</td>
      <td>-0.81</td>
      <td>26</td>
      <td>2</td>
      <td>9.64</td>
      <td>8.89</td>
      <td>NaN</td>
      <td>6.97</td>
      <td>-0.62</td>
      <td>27</td>
      <td>2</td>
      <td>6.87</td>
      <td>7.08</td>
      <td>155</td>
      <td>P155K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>L</td>
      <td>P486L</td>
      <td>9.67</td>
      <td>-0.40</td>
      <td>40</td>
      <td>2</td>
      <td>9.97</td>
      <td>9.37</td>
      <td>NaN</td>
      <td>7.52</td>
      <td>-0.07</td>
      <td>42</td>
      <td>2</td>
      <td>7.52</td>
      <td>7.52</td>
      <td>155</td>
      <td>P155L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>M</td>
      <td>P486M</td>
      <td>9.44</td>
      <td>-0.63</td>
      <td>31</td>
      <td>2</td>
      <td>9.68</td>
      <td>9.21</td>
      <td>NaN</td>
      <td>7.34</td>
      <td>-0.25</td>
      <td>36</td>
      <td>2</td>
      <td>7.22</td>
      <td>7.46</td>
      <td>155</td>
      <td>P155M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>N</td>
      <td>P486N</td>
      <td>9.17</td>
      <td>-0.90</td>
      <td>29</td>
      <td>2</td>
      <td>9.54</td>
      <td>8.81</td>
      <td>NaN</td>
      <td>7.33</td>
      <td>-0.26</td>
      <td>33</td>
      <td>2</td>
      <td>7.27</td>
      <td>7.40</td>
      <td>155</td>
      <td>P155N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>P</td>
      <td>P486P</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>155</td>
      <td>P155P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>Q</td>
      <td>P486Q</td>
      <td>9.42</td>
      <td>-0.65</td>
      <td>40</td>
      <td>2</td>
      <td>9.66</td>
      <td>9.18</td>
      <td>NaN</td>
      <td>7.39</td>
      <td>-0.21</td>
      <td>46</td>
      <td>2</td>
      <td>7.31</td>
      <td>7.46</td>
      <td>155</td>
      <td>P155Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>R</td>
      <td>P486R</td>
      <td>9.38</td>
      <td>-0.69</td>
      <td>44</td>
      <td>2</td>
      <td>9.72</td>
      <td>9.04</td>
      <td>NaN</td>
      <td>6.98</td>
      <td>-0.61</td>
      <td>46</td>
      <td>2</td>
      <td>6.85</td>
      <td>7.11</td>
      <td>155</td>
      <td>P155R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>S</td>
      <td>P486S</td>
      <td>9.41</td>
      <td>-0.66</td>
      <td>43</td>
      <td>2</td>
      <td>9.74</td>
      <td>9.09</td>
      <td>NaN</td>
      <td>7.31</td>
      <td>-0.29</td>
      <td>50</td>
      <td>2</td>
      <td>7.21</td>
      <td>7.40</td>
      <td>155</td>
      <td>P155S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>T</td>
      <td>P486T</td>
      <td>9.37</td>
      <td>-0.70</td>
      <td>29</td>
      <td>2</td>
      <td>9.65</td>
      <td>9.10</td>
      <td>NaN</td>
      <td>7.41</td>
      <td>-0.19</td>
      <td>32</td>
      <td>2</td>
      <td>7.33</td>
      <td>7.48</td>
      <td>155</td>
      <td>P155T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>V</td>
      <td>P486V</td>
      <td>9.87</td>
      <td>-0.20</td>
      <td>36</td>
      <td>2</td>
      <td>10.17</td>
      <td>9.58</td>
      <td>NaN</td>
      <td>7.51</td>
      <td>-0.08</td>
      <td>45</td>
      <td>2</td>
      <td>7.39</td>
      <td>7.64</td>
      <td>155</td>
      <td>P155V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>W</td>
      <td>P486W</td>
      <td>9.76</td>
      <td>-0.31</td>
      <td>39</td>
      <td>2</td>
      <td>10.01</td>
      <td>9.50</td>
      <td>NaN</td>
      <td>7.28</td>
      <td>-0.31</td>
      <td>44</td>
      <td>2</td>
      <td>7.25</td>
      <td>7.32</td>
      <td>155</td>
      <td>P155W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>P</td>
      <td>486</td>
      <td>Y</td>
      <td>P486Y</td>
      <td>9.64</td>
      <td>-0.43</td>
      <td>30</td>
      <td>2</td>
      <td>9.86</td>
      <td>9.43</td>
      <td>NaN</td>
      <td>7.19</td>
      <td>-0.40</td>
      <td>35</td>
      <td>2</td>
      <td>6.95</td>
      <td>7.43</td>
      <td>155</td>
      <td>P155Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>-</td>
      <td>N487-</td>
      <td>6.07</td>
      <td>-4.00</td>
      <td>23</td>
      <td>2</td>
      <td>6.41</td>
      <td>5.73</td>
      <td>NaN</td>
      <td>6.83</td>
      <td>-0.76</td>
      <td>26</td>
      <td>2</td>
      <td>6.83</td>
      <td>6.84</td>
      <td>156</td>
      <td>N156-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>A</td>
      <td>N487A</td>
      <td>8.69</td>
      <td>-1.38</td>
      <td>39</td>
      <td>2</td>
      <td>9.14</td>
      <td>8.23</td>
      <td>NaN</td>
      <td>7.42</td>
      <td>-0.18</td>
      <td>41</td>
      <td>2</td>
      <td>7.24</td>
      <td>7.59</td>
      <td>156</td>
      <td>N156A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>C</td>
      <td>N487C</td>
      <td>7.82</td>
      <td>-2.25</td>
      <td>23</td>
      <td>2</td>
      <td>8.08</td>
      <td>7.55</td>
      <td>NaN</td>
      <td>7.23</td>
      <td>-0.36</td>
      <td>25</td>
      <td>2</td>
      <td>7.09</td>
      <td>7.37</td>
      <td>156</td>
      <td>N156C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>D</td>
      <td>N487D</td>
      <td>9.36</td>
      <td>-0.71</td>
      <td>28</td>
      <td>2</td>
      <td>9.69</td>
      <td>9.03</td>
      <td>NaN</td>
      <td>7.70</td>
      <td>0.11</td>
      <td>30</td>
      <td>2</td>
      <td>7.71</td>
      <td>7.69</td>
      <td>156</td>
      <td>N156D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>E</td>
      <td>N487E</td>
      <td>8.94</td>
      <td>-1.13</td>
      <td>29</td>
      <td>2</td>
      <td>9.28</td>
      <td>8.60</td>
      <td>NaN</td>
      <td>7.82</td>
      <td>0.23</td>
      <td>35</td>
      <td>2</td>
      <td>7.78</td>
      <td>7.85</td>
      <td>156</td>
      <td>N156E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>F</td>
      <td>N487F</td>
      <td>8.52</td>
      <td>-1.55</td>
      <td>29</td>
      <td>2</td>
      <td>9.05</td>
      <td>8.00</td>
      <td>NaN</td>
      <td>7.56</td>
      <td>-0.03</td>
      <td>32</td>
      <td>2</td>
      <td>7.45</td>
      <td>7.67</td>
      <td>156</td>
      <td>N156F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>G</td>
      <td>N487G</td>
      <td>9.04</td>
      <td>-1.03</td>
      <td>34</td>
      <td>2</td>
      <td>9.42</td>
      <td>8.66</td>
      <td>NaN</td>
      <td>7.69</td>
      <td>0.10</td>
      <td>37</td>
      <td>2</td>
      <td>7.47</td>
      <td>7.90</td>
      <td>156</td>
      <td>N156G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>H</td>
      <td>N487H</td>
      <td>8.82</td>
      <td>-1.25</td>
      <td>16</td>
      <td>2</td>
      <td>9.02</td>
      <td>8.63</td>
      <td>NaN</td>
      <td>7.12</td>
      <td>-0.47</td>
      <td>19</td>
      <td>2</td>
      <td>6.98</td>
      <td>7.26</td>
      <td>156</td>
      <td>N156H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>I</td>
      <td>N487I</td>
      <td>7.86</td>
      <td>-2.21</td>
      <td>27</td>
      <td>2</td>
      <td>7.92</td>
      <td>7.79</td>
      <td>NaN</td>
      <td>6.86</td>
      <td>-0.73</td>
      <td>31</td>
      <td>2</td>
      <td>6.81</td>
      <td>6.91</td>
      <td>156</td>
      <td>N156I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>K</td>
      <td>N487K</td>
      <td>8.96</td>
      <td>-1.11</td>
      <td>22</td>
      <td>2</td>
      <td>9.37</td>
      <td>8.56</td>
      <td>NaN</td>
      <td>7.02</td>
      <td>-0.57</td>
      <td>28</td>
      <td>2</td>
      <td>6.94</td>
      <td>7.10</td>
      <td>156</td>
      <td>N156K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>L</td>
      <td>N487L</td>
      <td>8.70</td>
      <td>-1.37</td>
      <td>29</td>
      <td>2</td>
      <td>9.13</td>
      <td>8.27</td>
      <td>NaN</td>
      <td>7.66</td>
      <td>0.07</td>
      <td>32</td>
      <td>2</td>
      <td>7.62</td>
      <td>7.70</td>
      <td>156</td>
      <td>N156L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>M</td>
      <td>N487M</td>
      <td>8.90</td>
      <td>-1.17</td>
      <td>28</td>
      <td>2</td>
      <td>9.39</td>
      <td>8.42</td>
      <td>NaN</td>
      <td>7.48</td>
      <td>-0.11</td>
      <td>34</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.48</td>
      <td>156</td>
      <td>N156M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>N</td>
      <td>N487N</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>156</td>
      <td>N156N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>P</td>
      <td>N487P</td>
      <td>6.62</td>
      <td>-3.45</td>
      <td>17</td>
      <td>2</td>
      <td>7.08</td>
      <td>6.15</td>
      <td>NaN</td>
      <td>6.97</td>
      <td>-0.62</td>
      <td>21</td>
      <td>2</td>
      <td>6.86</td>
      <td>7.08</td>
      <td>156</td>
      <td>N156P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>Q</td>
      <td>N487Q</td>
      <td>9.24</td>
      <td>-0.83</td>
      <td>21</td>
      <td>2</td>
      <td>9.61</td>
      <td>8.88</td>
      <td>NaN</td>
      <td>7.63</td>
      <td>0.04</td>
      <td>22</td>
      <td>2</td>
      <td>7.48</td>
      <td>7.78</td>
      <td>156</td>
      <td>N156Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>R</td>
      <td>N487R</td>
      <td>8.50</td>
      <td>-1.57</td>
      <td>25</td>
      <td>2</td>
      <td>8.74</td>
      <td>8.26</td>
      <td>NaN</td>
      <td>7.04</td>
      <td>-0.55</td>
      <td>29</td>
      <td>2</td>
      <td>6.82</td>
      <td>7.27</td>
      <td>156</td>
      <td>N156R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>S</td>
      <td>N487S</td>
      <td>9.05</td>
      <td>-1.02</td>
      <td>28</td>
      <td>2</td>
      <td>9.45</td>
      <td>8.65</td>
      <td>NaN</td>
      <td>7.69</td>
      <td>0.10</td>
      <td>30</td>
      <td>2</td>
      <td>7.61</td>
      <td>7.77</td>
      <td>156</td>
      <td>N156S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>T</td>
      <td>N487T</td>
      <td>8.96</td>
      <td>-1.11</td>
      <td>20</td>
      <td>2</td>
      <td>9.29</td>
      <td>8.63</td>
      <td>NaN</td>
      <td>7.45</td>
      <td>-0.14</td>
      <td>25</td>
      <td>2</td>
      <td>7.29</td>
      <td>7.61</td>
      <td>156</td>
      <td>N156T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>V</td>
      <td>N487V</td>
      <td>8.01</td>
      <td>-2.06</td>
      <td>27</td>
      <td>2</td>
      <td>8.21</td>
      <td>7.82</td>
      <td>NaN</td>
      <td>6.97</td>
      <td>-0.62</td>
      <td>27</td>
      <td>2</td>
      <td>6.92</td>
      <td>7.02</td>
      <td>156</td>
      <td>N156V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>W</td>
      <td>N487W</td>
      <td>8.03</td>
      <td>-2.04</td>
      <td>30</td>
      <td>2</td>
      <td>8.38</td>
      <td>7.67</td>
      <td>NaN</td>
      <td>7.29</td>
      <td>-0.30</td>
      <td>32</td>
      <td>2</td>
      <td>7.23</td>
      <td>7.36</td>
      <td>156</td>
      <td>N156W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>N</td>
      <td>487</td>
      <td>Y</td>
      <td>N487Y</td>
      <td>8.16</td>
      <td>-1.91</td>
      <td>27</td>
      <td>2</td>
      <td>8.36</td>
      <td>7.95</td>
      <td>NaN</td>
      <td>7.42</td>
      <td>-0.17</td>
      <td>31</td>
      <td>2</td>
      <td>7.25</td>
      <td>7.60</td>
      <td>156</td>
      <td>N156Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>-</td>
      <td>C488-</td>
      <td>5.64</td>
      <td>-4.43</td>
      <td>36</td>
      <td>2</td>
      <td>6.09</td>
      <td>5.19</td>
      <td>NaN</td>
      <td>6.88</td>
      <td>-0.72</td>
      <td>38</td>
      <td>2</td>
      <td>6.79</td>
      <td>6.97</td>
      <td>157</td>
      <td>C157-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>A</td>
      <td>C488A</td>
      <td>5.83</td>
      <td>-4.24</td>
      <td>34</td>
      <td>2</td>
      <td>6.30</td>
      <td>5.37</td>
      <td>NaN</td>
      <td>6.78</td>
      <td>-0.81</td>
      <td>37</td>
      <td>2</td>
      <td>6.57</td>
      <td>7.00</td>
      <td>157</td>
      <td>C157A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>C</td>
      <td>C488C</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>157</td>
      <td>C157C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>D</td>
      <td>C488D</td>
      <td>5.56</td>
      <td>-4.52</td>
      <td>24</td>
      <td>2</td>
      <td>5.81</td>
      <td>5.30</td>
      <td>NaN</td>
      <td>6.85</td>
      <td>-0.75</td>
      <td>27</td>
      <td>2</td>
      <td>6.64</td>
      <td>7.06</td>
      <td>157</td>
      <td>C157D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>E</td>
      <td>C488E</td>
      <td>5.97</td>
      <td>-4.10</td>
      <td>25</td>
      <td>2</td>
      <td>6.16</td>
      <td>5.78</td>
      <td>NaN</td>
      <td>6.74</td>
      <td>-0.86</td>
      <td>27</td>
      <td>2</td>
      <td>6.61</td>
      <td>6.86</td>
      <td>157</td>
      <td>C157E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>F</td>
      <td>C488F</td>
      <td>5.92</td>
      <td>-4.15</td>
      <td>26</td>
      <td>2</td>
      <td>6.61</td>
      <td>5.23</td>
      <td>NaN</td>
      <td>6.68</td>
      <td>-0.91</td>
      <td>29</td>
      <td>2</td>
      <td>6.52</td>
      <td>6.85</td>
      <td>157</td>
      <td>C157F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>G</td>
      <td>C488G</td>
      <td>6.18</td>
      <td>-3.89</td>
      <td>34</td>
      <td>2</td>
      <td>6.45</td>
      <td>5.91</td>
      <td>NaN</td>
      <td>6.75</td>
      <td>-0.85</td>
      <td>40</td>
      <td>2</td>
      <td>6.75</td>
      <td>6.74</td>
      <td>157</td>
      <td>C157G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>H</td>
      <td>C488H</td>
      <td>5.73</td>
      <td>-4.34</td>
      <td>28</td>
      <td>2</td>
      <td>6.05</td>
      <td>5.42</td>
      <td>NaN</td>
      <td>6.74</td>
      <td>-0.85</td>
      <td>35</td>
      <td>2</td>
      <td>6.59</td>
      <td>6.88</td>
      <td>157</td>
      <td>C157H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>I</td>
      <td>C488I</td>
      <td>5.82</td>
      <td>-4.25</td>
      <td>25</td>
      <td>2</td>
      <td>6.26</td>
      <td>5.37</td>
      <td>NaN</td>
      <td>6.58</td>
      <td>-1.01</td>
      <td>31</td>
      <td>2</td>
      <td>6.46</td>
      <td>6.71</td>
      <td>157</td>
      <td>C157I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>K</td>
      <td>C488K</td>
      <td>5.71</td>
      <td>-4.36</td>
      <td>31</td>
      <td>2</td>
      <td>5.97</td>
      <td>5.46</td>
      <td>NaN</td>
      <td>6.82</td>
      <td>-0.77</td>
      <td>33</td>
      <td>2</td>
      <td>6.68</td>
      <td>6.96</td>
      <td>157</td>
      <td>C157K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>L</td>
      <td>C488L</td>
      <td>5.60</td>
      <td>-4.47</td>
      <td>29</td>
      <td>2</td>
      <td>5.81</td>
      <td>5.38</td>
      <td>NaN</td>
      <td>6.66</td>
      <td>-0.93</td>
      <td>30</td>
      <td>2</td>
      <td>6.48</td>
      <td>6.85</td>
      <td>157</td>
      <td>C157L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>M</td>
      <td>C488M</td>
      <td>5.94</td>
      <td>-4.13</td>
      <td>39</td>
      <td>2</td>
      <td>6.23</td>
      <td>5.65</td>
      <td>NaN</td>
      <td>6.78</td>
      <td>-0.81</td>
      <td>41</td>
      <td>2</td>
      <td>6.64</td>
      <td>6.92</td>
      <td>157</td>
      <td>C157M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>N</td>
      <td>C488N</td>
      <td>6.06</td>
      <td>-4.01</td>
      <td>17</td>
      <td>2</td>
      <td>6.54</td>
      <td>5.58</td>
      <td>NaN</td>
      <td>6.75</td>
      <td>-0.85</td>
      <td>18</td>
      <td>2</td>
      <td>6.64</td>
      <td>6.85</td>
      <td>157</td>
      <td>C157N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>P</td>
      <td>C488P</td>
      <td>5.73</td>
      <td>-4.34</td>
      <td>43</td>
      <td>2</td>
      <td>5.85</td>
      <td>5.61</td>
      <td>NaN</td>
      <td>6.79</td>
      <td>-0.80</td>
      <td>43</td>
      <td>2</td>
      <td>6.56</td>
      <td>7.03</td>
      <td>157</td>
      <td>C157P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>Q</td>
      <td>C488Q</td>
      <td>5.76</td>
      <td>-4.31</td>
      <td>30</td>
      <td>2</td>
      <td>5.98</td>
      <td>5.54</td>
      <td>NaN</td>
      <td>6.88</td>
      <td>-0.71</td>
      <td>34</td>
      <td>2</td>
      <td>6.77</td>
      <td>7.00</td>
      <td>157</td>
      <td>C157Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>R</td>
      <td>C488R</td>
      <td>5.72</td>
      <td>-4.35</td>
      <td>36</td>
      <td>2</td>
      <td>6.07</td>
      <td>5.38</td>
      <td>NaN</td>
      <td>6.72</td>
      <td>-0.87</td>
      <td>40</td>
      <td>2</td>
      <td>6.67</td>
      <td>6.77</td>
      <td>157</td>
      <td>C157R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>S</td>
      <td>C488S</td>
      <td>5.98</td>
      <td>-4.09</td>
      <td>40</td>
      <td>2</td>
      <td>6.30</td>
      <td>5.65</td>
      <td>NaN</td>
      <td>6.90</td>
      <td>-0.69</td>
      <td>45</td>
      <td>2</td>
      <td>6.72</td>
      <td>7.08</td>
      <td>157</td>
      <td>C157S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>T</td>
      <td>C488T</td>
      <td>5.90</td>
      <td>-4.17</td>
      <td>46</td>
      <td>2</td>
      <td>6.41</td>
      <td>5.39</td>
      <td>NaN</td>
      <td>6.82</td>
      <td>-0.77</td>
      <td>49</td>
      <td>2</td>
      <td>6.71</td>
      <td>6.93</td>
      <td>157</td>
      <td>C157T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>V</td>
      <td>C488V</td>
      <td>5.79</td>
      <td>-4.28</td>
      <td>22</td>
      <td>2</td>
      <td>6.37</td>
      <td>5.20</td>
      <td>NaN</td>
      <td>6.73</td>
      <td>-0.86</td>
      <td>24</td>
      <td>2</td>
      <td>6.61</td>
      <td>6.86</td>
      <td>157</td>
      <td>C157V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>W</td>
      <td>C488W</td>
      <td>5.56</td>
      <td>-4.51</td>
      <td>54</td>
      <td>2</td>
      <td>5.82</td>
      <td>5.29</td>
      <td>NaN</td>
      <td>6.65</td>
      <td>-0.94</td>
      <td>57</td>
      <td>2</td>
      <td>6.55</td>
      <td>6.75</td>
      <td>157</td>
      <td>C157W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>C</td>
      <td>488</td>
      <td>Y</td>
      <td>C488Y</td>
      <td>5.52</td>
      <td>-4.55</td>
      <td>28</td>
      <td>2</td>
      <td>5.83</td>
      <td>5.20</td>
      <td>NaN</td>
      <td>6.59</td>
      <td>-1.00</td>
      <td>28</td>
      <td>2</td>
      <td>6.46</td>
      <td>6.73</td>
      <td>157</td>
      <td>C157Y</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>-</td>
      <td>Y489-</td>
      <td>5.70</td>
      <td>-4.37</td>
      <td>42</td>
      <td>2</td>
      <td>5.99</td>
      <td>5.40</td>
      <td>NaN</td>
      <td>7.12</td>
      <td>-0.47</td>
      <td>44</td>
      <td>2</td>
      <td>6.99</td>
      <td>7.26</td>
      <td>158</td>
      <td>Y158-</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>A</td>
      <td>Y489A</td>
      <td>6.45</td>
      <td>-3.62</td>
      <td>24</td>
      <td>2</td>
      <td>6.51</td>
      <td>6.40</td>
      <td>NaN</td>
      <td>7.08</td>
      <td>-0.51</td>
      <td>27</td>
      <td>2</td>
      <td>6.91</td>
      <td>7.26</td>
      <td>158</td>
      <td>Y158A</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>C</td>
      <td>Y489C</td>
      <td>5.76</td>
      <td>-4.31</td>
      <td>27</td>
      <td>2</td>
      <td>6.17</td>
      <td>5.35</td>
      <td>NaN</td>
      <td>7.17</td>
      <td>-0.42</td>
      <td>30</td>
      <td>2</td>
      <td>7.11</td>
      <td>7.23</td>
      <td>158</td>
      <td>Y158C</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>D</td>
      <td>Y489D</td>
      <td>5.77</td>
      <td>-4.30</td>
      <td>31</td>
      <td>2</td>
      <td>6.16</td>
      <td>5.39</td>
      <td>NaN</td>
      <td>7.11</td>
      <td>-0.48</td>
      <td>33</td>
      <td>2</td>
      <td>6.91</td>
      <td>7.31</td>
      <td>158</td>
      <td>Y158D</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>E</td>
      <td>Y489E</td>
      <td>6.42</td>
      <td>-3.65</td>
      <td>31</td>
      <td>2</td>
      <td>6.99</td>
      <td>5.85</td>
      <td>NaN</td>
      <td>7.24</td>
      <td>-0.35</td>
      <td>35</td>
      <td>2</td>
      <td>7.13</td>
      <td>7.35</td>
      <td>158</td>
      <td>Y158E</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>F</td>
      <td>Y489F</td>
      <td>9.05</td>
      <td>-1.02</td>
      <td>33</td>
      <td>2</td>
      <td>9.55</td>
      <td>8.55</td>
      <td>NaN</td>
      <td>7.78</td>
      <td>0.19</td>
      <td>35</td>
      <td>2</td>
      <td>7.65</td>
      <td>7.91</td>
      <td>158</td>
      <td>Y158F</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>G</td>
      <td>Y489G</td>
      <td>5.68</td>
      <td>-4.39</td>
      <td>24</td>
      <td>2</td>
      <td>5.95</td>
      <td>5.40</td>
      <td>NaN</td>
      <td>7.07</td>
      <td>-0.52</td>
      <td>26</td>
      <td>2</td>
      <td>6.94</td>
      <td>7.21</td>
      <td>158</td>
      <td>Y158G</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>H</td>
      <td>Y489H</td>
      <td>8.91</td>
      <td>-1.16</td>
      <td>32</td>
      <td>2</td>
      <td>9.32</td>
      <td>8.49</td>
      <td>NaN</td>
      <td>7.14</td>
      <td>-0.45</td>
      <td>36</td>
      <td>2</td>
      <td>6.92</td>
      <td>7.36</td>
      <td>158</td>
      <td>Y158H</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>I</td>
      <td>Y489I</td>
      <td>7.88</td>
      <td>-2.19</td>
      <td>21</td>
      <td>2</td>
      <td>7.98</td>
      <td>7.78</td>
      <td>NaN</td>
      <td>7.56</td>
      <td>-0.03</td>
      <td>22</td>
      <td>2</td>
      <td>7.52</td>
      <td>7.60</td>
      <td>158</td>
      <td>Y158I</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>K</td>
      <td>Y489K</td>
      <td>6.08</td>
      <td>-3.99</td>
      <td>40</td>
      <td>2</td>
      <td>6.22</td>
      <td>5.94</td>
      <td>NaN</td>
      <td>6.88</td>
      <td>-0.71</td>
      <td>43</td>
      <td>2</td>
      <td>6.79</td>
      <td>6.97</td>
      <td>158</td>
      <td>Y158K</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>L</td>
      <td>Y489L</td>
      <td>7.84</td>
      <td>-2.23</td>
      <td>27</td>
      <td>2</td>
      <td>8.25</td>
      <td>7.42</td>
      <td>NaN</td>
      <td>7.61</td>
      <td>0.01</td>
      <td>32</td>
      <td>2</td>
      <td>7.46</td>
      <td>7.75</td>
      <td>158</td>
      <td>Y158L</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>M</td>
      <td>Y489M</td>
      <td>7.88</td>
      <td>-2.19</td>
      <td>35</td>
      <td>2</td>
      <td>8.27</td>
      <td>7.49</td>
      <td>NaN</td>
      <td>7.45</td>
      <td>-0.14</td>
      <td>41</td>
      <td>2</td>
      <td>7.38</td>
      <td>7.53</td>
      <td>158</td>
      <td>Y158M</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>N</td>
      <td>Y489N</td>
      <td>6.65</td>
      <td>-3.42</td>
      <td>32</td>
      <td>2</td>
      <td>6.99</td>
      <td>6.30</td>
      <td>NaN</td>
      <td>7.12</td>
      <td>-0.47</td>
      <td>33</td>
      <td>2</td>
      <td>7.05</td>
      <td>7.19</td>
      <td>158</td>
      <td>Y158N</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>P</td>
      <td>Y489P</td>
      <td>6.89</td>
      <td>-3.18</td>
      <td>24</td>
      <td>2</td>
      <td>7.25</td>
      <td>6.54</td>
      <td>NaN</td>
      <td>7.14</td>
      <td>-0.45</td>
      <td>29</td>
      <td>2</td>
      <td>7.01</td>
      <td>7.28</td>
      <td>158</td>
      <td>Y158P</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>Q</td>
      <td>Y489Q</td>
      <td>7.01</td>
      <td>-3.06</td>
      <td>32</td>
      <td>2</td>
      <td>7.16</td>
      <td>6.85</td>
      <td>NaN</td>
      <td>7.12</td>
      <td>-0.47</td>
      <td>33</td>
      <td>2</td>
      <td>6.99</td>
      <td>7.25</td>
      <td>158</td>
      <td>Y158Q</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>R</td>
      <td>Y489R</td>
      <td>5.98</td>
      <td>-4.09</td>
      <td>32</td>
      <td>2</td>
      <td>6.29</td>
      <td>5.67</td>
      <td>NaN</td>
      <td>7.01</td>
      <td>-0.58</td>
      <td>36</td>
      <td>2</td>
      <td>6.85</td>
      <td>7.17</td>
      <td>158</td>
      <td>Y158R</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>S</td>
      <td>Y489S</td>
      <td>5.87</td>
      <td>-4.20</td>
      <td>28</td>
      <td>2</td>
      <td>6.27</td>
      <td>5.47</td>
      <td>NaN</td>
      <td>7.72</td>
      <td>0.13</td>
      <td>32</td>
      <td>2</td>
      <td>7.62</td>
      <td>7.83</td>
      <td>158</td>
      <td>Y158S</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>T</td>
      <td>Y489T</td>
      <td>5.64</td>
      <td>-4.43</td>
      <td>20</td>
      <td>2</td>
      <td>5.96</td>
      <td>5.31</td>
      <td>NaN</td>
      <td>8.02</td>
      <td>0.42</td>
      <td>22</td>
      <td>2</td>
      <td>7.83</td>
      <td>8.21</td>
      <td>158</td>
      <td>Y158T</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>V</td>
      <td>Y489V</td>
      <td>7.44</td>
      <td>-2.63</td>
      <td>32</td>
      <td>2</td>
      <td>8.03</td>
      <td>6.85</td>
      <td>NaN</td>
      <td>7.54</td>
      <td>-0.05</td>
      <td>35</td>
      <td>2</td>
      <td>7.36</td>
      <td>7.71</td>
      <td>158</td>
      <td>Y158V</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>W</td>
      <td>Y489W</td>
      <td>9.27</td>
      <td>-0.81</td>
      <td>38</td>
      <td>2</td>
      <td>9.62</td>
      <td>8.91</td>
      <td>NaN</td>
      <td>7.52</td>
      <td>-0.08</td>
      <td>41</td>
      <td>2</td>
      <td>7.31</td>
      <td>7.72</td>
      <td>158</td>
      <td>Y158W</td>
    </tr>
    <tr>
      <td>Omicron_BA286</td>
      <td>Y</td>
      <td>489</td>
      <td>Y</td>
      <td>Y489Y</td>
      <td>10.07</td>
      <td>0.00</td>
      <td>4048</td>
      <td>2</td>
      <td>10.33</td>
      <td>9.81</td>
      <td>NaN</td>
      <td>7.59</td>
      <td>0.00</td>
      <td>5260</td>
      <td>2</td>
      <td>7.49</td>
      <td>7.70</td>
      <td>158</td>
      <td>Y158Y</td>
    </tr>
  </tbody>
</table>



```python
assert mut_bind_expr['RBD_mutation'].nunique() == len(mut_bind_expr)
for prop in ['bind', 'expr']:
    muts_adequate = set(mut_bind_expr
                        .query(f"delta_{prop} >= {config[f'escape_score_min_{prop}_mut_Omicron_BA286']}")
                        ['RBD_mutation']
                        )
    print(f"{len(muts_adequate)} of {len(mut_bind_expr)} mutations have adequate {prop}.")
    escape_scores[f"muts_pass_{prop}_filter"] = (
        escape_scores
        ['aa_substitutions']
        .map(lambda s: set(s.split()).issubset(muts_adequate))
        )

# annotate as passing overall filter if passes all mutation and binding filters:
escape_scores['pass_ACE2bind_expr_filter'] = (
        escape_scores['muts_pass_bind_filter'] &
        escape_scores['muts_pass_expr_filter'] 
        )

display(HTML(escape_scores.query('not pass_ACE2bind_expr_filter & variant_class != "stop"').head().to_html(index=False)))
```

    3314 of 4200 mutations have adequate bind.
    2437 of 4200 mutations have adequate expr.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pass_pre_count_filter</th>
      <th>muts_pass_bind_filter</th>
      <th>muts_pass_expr_filter</th>
      <th>pass_ACE2bind_expr_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AACCGTTAAATGCAAA</td>
      <td>0.009930</td>
      <td>1.582215e-06</td>
      <td>21129</td>
      <td>62</td>
      <td>AAA73CAT</td>
      <td>1</td>
      <td>K73H</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AAAAAATCAAACTGAA</td>
      <td>0.549101</td>
      <td>1.248995e-04</td>
      <td>17127</td>
      <td>2801</td>
      <td>TAC177GAT</td>
      <td>1</td>
      <td>Y177D</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CAATATCACCCAAGTG</td>
      <td>0.002073</td>
      <td>4.093863e-07</td>
      <td>17006</td>
      <td>10</td>
      <td>GGT165GAA</td>
      <td>1</td>
      <td>G165E</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CTAACTGACCCACTAA</td>
      <td>0.013908</td>
      <td>3.011216e-06</td>
      <td>15568</td>
      <td>64</td>
      <td>TCC69TTG</td>
      <td>1</td>
      <td>S69L</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>GTATAAGACAGAAGAG</td>
      <td>0.006827</td>
      <td>1.531027e-06</td>
      <td>14998</td>
      <td>30</td>
      <td>GGT146AAA</td>
      <td>1</td>
      <td>G146K</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Print the number of mutations that pass RBD bind, RBD expression, and are not to sites that are disulfide bonds (if specified in config) 


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
if config['exclude_cysteines']:
    print("Here are the number of mutations that pass the bind, express, and disulfide filters:")
    print(len(mut_bind_expr
              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut_Wuhan_Hu_1']} & \
                       delta_expr >= {config[f'escape_score_min_expr_mut_Wuhan_Hu_1']} & \
                       pass_cysteine_filter")
             ))
    print("There are these many possible mutations (excluding wildtype and disulfides!):")
    print(mut_bind_expr.query('wildtype!="C"')['position'].nunique()*19
         )

else:
    print("Here are the number of mutations that pass the bind and express filters:")
    print(len(mut_bind_expr
              .assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C")
              .query(f"delta_bind >= {config[f'escape_score_min_bind_mut_Wuhan_Hu_1']} & \
                       delta_expr >= {config[f'escape_score_min_expr_mut_Wuhan_Hu_1']}")
             ))
    print("There are these many possible mutations (excluding wildtype!):")
    print(mut_bind_expr['position'].nunique()*19
         )
```

    Here are the number of mutations that pass the bind, express, and disulfide filters:
    3002
    There are these many possible mutations (excluding wildtype and disulfides!):
    3648



```python
print('These are the sites that are involved in disulfide bonds:')
print(mut_bind_expr.query('wildtype=="C"')['position'].unique())
```

    These are the sites that are involved in disulfide bonds:
    [336 361 379 391 432 480 488 525]



```python
frac_ACE2bind_expr_pass_filter = (
    escape_scores
    .query('pass_pre_count_filter == True')
    [['pre_sample', 'library', 'target', config['escape_score_group_by'],
      'pre_count', 'pass_ACE2bind_expr_filter', 'variant_class']]
    .drop_duplicates()
    .groupby(['pre_sample', 'library', 'variant_class'], observed=True)
    .aggregate(n_variants=pd.NamedAgg('pass_ACE2bind_expr_filter', 'count'),
               n_pass_filter=pd.NamedAgg('pass_ACE2bind_expr_filter', 'sum')
               )
    .reset_index()
    .assign(frac_pass_filter=lambda x: x['n_pass_filter'] / x['n_variants'],
            pre_sample=lambda x: pd.Categorical(x['pre_sample'], x['pre_sample'].unique(), ordered=True).remove_unused_categories())
    )

p = (ggplot(frac_ACE2bind_expr_pass_filter) +
     aes('variant_class', 'frac_pass_filter', fill='variant_class') +
     geom_bar(stat='identity') +
     facet_grid('library ~ pre_sample') +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(3.3 * frac_ACE2bind_expr_pass_filter['pre_sample'].nunique(),
                        2 * frac_ACE2bind_expr_pass_filter['library'].nunique()),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_58_0.png)
    


## Examine and write escape scores
Plot the distribution of escape scores across variants of different classes **among those that pass both the pre-selection count filter and the ACE2-binding / expression filter**.
If things are working correctly, we don't expect escape in wildtype (or synonymous variants), but do expect escape for some small fraction of nonsynymous variants.
Also, we do not plot the scores for the stop codon variant class, as most stop-codon variants have already been filtered out so this category is largely noise:


```python
nfacets = len(escape_scores.groupby(['library', 'name']).nunique())
ncol = min(8, nfacets)
nrow = math.ceil(nfacets / ncol)

df = (escape_scores
      .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
      .query('variant_class != "stop"')
      )
     
p = (ggplot(df) +
     aes('variant_class', 'score', color='variant_class') +
     geom_boxplot(outlier_size=1.5, outlier_alpha=0.1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 3 * nrow),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:])
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_60_0.png)
    


Also, we want to see how much the high escape scores are correlated with simple coverage.
To do this, we plot the correlation between escape score and pre-selection count just for the nonsynonymous variants (which are the ones that we expect to have true escape).
The plots below have a lot of overplotting, but are still sufficient to test of the score is simply correlated with the pre-selection counts or not.
The hoped for result is that the escape score doesn't appear to be strongly correlated with pre-selection counts:


```python
p = (ggplot(escape_scores
            .query('pass_pre_count_filter == True')
            .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
            .query('variant_class=="1 nonsynonymous"')
            ) +
     aes('pre_count', 'score') +
     geom_point(alpha=0.1, size=1) +
     facet_wrap('~ name + library', ncol=ncol) +
     theme(axis_text_x=element_text(angle=90),
           figure_size=(2.35 * ncol, 2.35 * nrow),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     scale_color_manual(values=CBPALETTE[1:]) +
     scale_x_log10()
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_62_0.png)
    


Write the escape scores to a file:


```python
print(f"Writing escape scores for Omicron_BA286 SSM to {config['escape_scores_Omicron_BA286']}")
escape_scores.to_csv(config['escape_scores_Omicron_BA286'], index=False, float_format='%.4g')
```

    Writing escape scores for Omicron_BA286 SSM to results/escape_scores/scores_Omicron_BA286.csv


### Now we will also remove anything that did not pass all the filters above. 


```python
escape_scores_primary = (escape_scores
                         .query('(pass_pre_count_filter == True) & (pass_ACE2bind_expr_filter == True)')
                        )

display(HTML(escape_scores_primary.head().to_html()))
print(f"Read {len(escape_scores_primary)} scores.")
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>name</th>
      <th>target</th>
      <th>library</th>
      <th>pre_sample</th>
      <th>post_sample</th>
      <th>barcode</th>
      <th>score</th>
      <th>score_var</th>
      <th>pre_count</th>
      <th>post_count</th>
      <th>codon_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_aa_substitutions</th>
      <th>variant_class</th>
      <th>pass_pre_count_filter</th>
      <th>muts_pass_bind_filter</th>
      <th>muts_pass_expr_filter</th>
      <th>pass_ACE2bind_expr_filter</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>3</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AGACGAAATAGATGAC</td>
      <td>0.007114</td>
      <td>1.428507e-06</td>
      <td>16752</td>
      <td>35</td>
      <td>ACC3---</td>
      <td>1</td>
      <td>T3-</td>
      <td>1</td>
      <td>deletion</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>6</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>TTACTTTAGCGACTCA</td>
      <td>0.005857</td>
      <td>1.347439e-06</td>
      <td>14616</td>
      <td>25</td>
      <td>TAC119AGA</td>
      <td>1</td>
      <td>Y119R</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>8</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>CTAAAATCCACCTCTA</td>
      <td>0.003830</td>
      <td>8.901341e-07</td>
      <td>14461</td>
      <td>16</td>
      <td>ACT200TTT</td>
      <td>1</td>
      <td>T200F</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>11</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>AAGGAATCTAGAGCGG</td>
      <td>0.002515</td>
      <td>6.026694e-07</td>
      <td>14017</td>
      <td>10</td>
      <td>TTG95ATT</td>
      <td>1</td>
      <td>L95I</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
    <tr>
      <th>12</th>
      <td>S2V29_36</td>
      <td>Omicron_BA286</td>
      <td>lib92</td>
      <td>exptREF-none-0-ref</td>
      <td>expt8-S2V29-36-abneg</td>
      <td>TTGCGTTCGGTGAAGG</td>
      <td>0.002998</td>
      <td>7.197201e-07</td>
      <td>13996</td>
      <td>12</td>
      <td>TTC159TGT</td>
      <td>1</td>
      <td>F159C</td>
      <td>1</td>
      <td>1 nonsynonymous</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
      <td>True</td>
    </tr>
  </tbody>
</table>


    Read 159636 scores.


### Count number of barcodes per mutation and remove variants with >1 amino acid substitution
Also add the number of barocdes per mutation to the `escape_scores` dataframe and plot this. 
But first see how many variants there are with >1 mutation, and query the dataframe to look at them qualitatively. 


```python
p = (ggplot(escape_scores_primary) +
     aes('n_aa_substitutions', fill='variant_class') +
     geom_bar() +
     facet_wrap('~library + pre_sample', ncol=5) +
     theme(
#            figure_size=(3.3 * escape_scores_primary['pre_sample'].nunique(),
#                         2 * escape_scores_primary['library'].nunique()),
         figure_size=(12, 4),
           panel_grid_major_x=element_blank(),
           ) +
     scale_fill_manual(values=CBPALETTE[1:]) +
     expand_limits(y=(0, 1))
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_68_0.png)
    


### Filter dataframe on single mutations that are present in at least `n` number of variants (specified in `config.yaml` file)
Now see how many `n_single_mut_measurements` there are for each variant:


```python
print(f'Remove anything with fewer than {config["escape_frac_min_single_mut_measurements"]} single mutant measurements (barcodes)')

raw_avg_single_mut_scores = (
    escape_scores_primary
    .query('n_aa_substitutions == 1')
    .rename(columns={'name': 'selection',
                     'aa_substitutions': 'mutation'})
    .groupby(['selection', 'library', 'mutation'])
    .aggregate(raw_single_mut_score=pd.NamedAgg('score', 'mean'),
               n_single_mut_measurements=pd.NamedAgg('barcode', 'count')
              )
    .assign(sufficient_measurements=lambda x: (
                (x['n_single_mut_measurements'] >= config['escape_frac_min_single_mut_measurements'])))
    .reset_index()
    )

# remove mutations with insufficient measurements
effects_df = (raw_avg_single_mut_scores
              .query('sufficient_measurements == True')
              .drop(columns='sufficient_measurements')
              )

# some error checks
assert len(effects_df) == len(effects_df.drop_duplicates()), 'duplicate rows in `effects_df`'
assert all(effects_df['raw_single_mut_score'].notnull() | (effects_df['n_single_mut_measurements'] == 0))

display(HTML(effects_df.head().to_html()))
```

    Remove anything with fewer than 2 single mutant measurements (barcodes)



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>selection</th>
      <th>library</th>
      <th>mutation</th>
      <th>raw_single_mut_score</th>
      <th>n_single_mut_measurements</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>A105I</td>
      <td>0.007985</td>
      <td>19</td>
    </tr>
    <tr>
      <th>1</th>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>A105S</td>
      <td>0.005138</td>
      <td>20</td>
    </tr>
    <tr>
      <th>2</th>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>A105T</td>
      <td>0.021114</td>
      <td>15</td>
    </tr>
    <tr>
      <th>3</th>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>A105V</td>
      <td>0.008595</td>
      <td>19</td>
    </tr>
    <tr>
      <th>4</th>
      <td>S2V29_36</td>
      <td>lib92</td>
      <td>A105W</td>
      <td>0.004693</td>
      <td>27</td>
    </tr>
  </tbody>
</table>


### Should we exclude mutations for which the wildtype identity is a cysteine?
These would be mutations that break disulfide bonds (unless there is an unpaired cysteine in the protein that does not form a disulfide bond, of course). 
Note that this approach would not work well for something like the SARS-CoV-2 NTD, where it has been documented (by the Veesler lab) that mutations in the B.1.427/429 (epislon) variant can rearrange disulfide bonds, leading to major structural rearrangements in the NTD, and yet an apparently fully functional spike. 

If we are excluding cysteines, do that now:


```python
# if we are excluding all cysteines to remove spurious mutations that break disulfide bonds:
if config['exclude_cysteines']:
    print(f'Excluding mutations where the wildtype identity is a cysteine')
    effects_df = effects_df.assign(pass_cysteine_filter=lambda x: x['mutation'].str[0] != "C",
                                  )
    disulfides_to_drop = effects_df.query('pass_cysteine_filter==False')['mutation'].unique()
    print(f'Specifically, excluding: {disulfides_to_drop}')
    effects_df=effects_df.query('pass_cysteine_filter').drop(columns='pass_cysteine_filter')

else:
    print(f'Retaining mutations where the wildtype identity is a cysteine')
```

    Excluding mutations where the wildtype identity is a cysteine
    Specifically, excluding: ['C194-' 'C194A' 'C194D' 'C194E' 'C194F' 'C194G' 'C194H' 'C194I' 'C194L'
     'C194M' 'C194N' 'C194P' 'C194Q' 'C194S' 'C194T' 'C194V' 'C194W' 'C194Y'
     'C31-' 'C31D' 'C31F' 'C31P' 'C49W' 'C61A' 'C61D' 'C61E' 'C61F' 'C61G'
     'C61H' 'C61I' 'C61K' 'C61L' 'C61M' 'C61N' 'C61P' 'C61Q' 'C61R' 'C61S'
     'C61T' 'C61V' 'C61Y' 'C6F' 'C6N']


We need to compute the escape scores (calculated as [here](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html?highlight=escape_scores#dms_variants.codonvarianttable.CodonVariantTable.escape_scores)) back to escape fractions. We define a function to do this depending on the score type:


```python
def score_to_frac(score):
    """Convert escape score to escape fraction."""
    if pd.isnull(score):
        return pd.NA
    floor = config['escape_score_floor_E']
    ceil = config['escape_score_ceil_E']
    if config['escape_score_type'] == 'frac_escape':
        return min(ceil, max(floor, score))  # just the score after applying ceiling and floor
    elif config['escape_score_type'] == 'log_escape':
        # convert score back to fraction, apply ceiling, then re-scale so floor is 0
        frac = 2**score
        frac = min(ceil, max(floor, frac))
        frac = (frac - floor) / (ceil - floor)
        return frac
    else:
        raise ValueError(f"invalid `escape_score_type` of {config['escape_score_type']}")

effects_df = (
    effects_df
    .assign(
            mut_escape_frac_single_mut=lambda x: x['raw_single_mut_score'].map(score_to_frac),
            )
    )
```

### Average the escape score across all barcodes of the same mutation, for each library, and the average of both libraries. 
Add rows that are the average of the two libraries for the fraction escape for all mutations that are present in both libraries (and if in just one library, the value in that or purge depending on config values printed here), the number of libraries each mutation is measured in, and the sum of the statistics giving the number of measurements:


```python
min_libs = config['escape_frac_avg_min_libraries']
min_single = config['escape_frac_avg_min_single']
print(f"Only taking average of mutations with escape fractions in >={min_libs} libraries "
      f"or with >={min_single} single-mutant measurements total.")

effects_df = (
    effects_df
    .query('library != "average"')  # git rid of averages if already there
    .assign(nlibs=1)
    .append(effects_df
            .query('library != "average"')
            .groupby(['selection', 'mutation'])
            .aggregate(nlibs=pd.NamedAgg('library', 'count'),
                       mut_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                              lambda s: s.mean(skipna=True)),
                       n_single_mut_measurements=pd.NamedAgg('n_single_mut_measurements', 'sum'),
                       )
            .query('(nlibs >= @min_libs) or (n_single_mut_measurements >= @min_single)')
            .assign(library="average")
            .reset_index(),
            ignore_index=True, sort=False,
            )
    )

print(len(effects_df.query('nlibs>1')))
print(len(effects_df.query('nlibs==1')))
```

    Only taking average of mutations with escape fractions in >=2 libraries or with >=2 single-mutant measurements total.
    4228
    8456


Plot the correlations of the escape fractions among the two libraries for all selections performed on both libraries. 


```python
libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

# wide data frame with each library's score in a different column
effects_df_wide = (
    effects_df
    .query('library != "average"')
    .query(f"n_single_mut_measurements >= 1")
    # just get selections with 2 libraries
    .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
    .query('nlibs == 2')
    # now make columns for each library, only keep mutants with scores for both libs
    [['selection', 'mutation', 'library', 'mut_escape_frac_single_mut']]
    .pivot_table(index=['selection', 'mutation'],
                 columns='library',
                 values='mut_escape_frac_single_mut',
                 aggfunc='first')
    .reset_index()
    .dropna(axis=0)
    )

# correlations between libraries
corrs = (
    effects_df_wide
    .groupby('selection')
    [libraries]
    .corr(method='pearson')
    .reset_index()
    .query('library == @libraries[0]')
    .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
    [['selection', 'correlation']]
    # add number of mutations measured
    .merge(effects_df_wide
           .groupby('selection')
           .size()
           .rename('n')
           .reset_index()
           )
    .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
    )

# plot correlations
nfacets = effects_df_wide['selection'].nunique()
ncol = min(nfacets, 5)
nrow = math.ceil(nfacets / ncol)
xmin = effects_df_wide[libraries[0]].min()
xspan = effects_df_wide[libraries[0]].max() - xmin
ymin = effects_df_wide[libraries[1]].min()
yspan = effects_df_wide[libraries[1]].max() - ymin
p = (ggplot(effects_df_wide) +
     aes(libraries[0], libraries[1]) +
     geom_point(alpha=0.2) +
     geom_text(mapping=aes(label='correlation'),
               data=corrs,
               x=0.01 * xspan + xmin,
               y=0.99 * yspan + ymin,
               size=10,
               ha='left',
               va='top',
               ) +
     facet_wrap('~ selection', ncol=ncol) +
     theme(figure_size=(2.5 * ncol, 2.5 * nrow),
           plot_title=element_text(size=14)) +
     ggtitle('Mutation-level escape fractions')
     )

_ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_78_0.png)
    


### Escape at site level
The above analysis estimates the effects of mutations. We also compute escape statistics at the site level. First, add sites to the data frame of mutational effects:


```python
effects_df = (
    effects_df
    .assign(site=lambda x: x['mutation'].str[1: -1].astype(int),
            wildtype=lambda x: x['mutation'].str[0],
            mutant=lambda x: x['mutation'].str[-1],
            )
    )
```

Now compute some site-level metrics. These are the average and total escape fraction at each site over all mutations at the site:


```python
site_effects_df = (
    effects_df
    .groupby(['selection', 'library', 'site'])
    .aggregate(
        site_avg_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut',
                                                    lambda s: s.mean(skipna=True)),
        site_total_escape_frac_single_mut=pd.NamedAgg('mut_escape_frac_single_mut', 'sum'),
        )
    .reset_index()
    )

display(HTML(site_effects_df.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>site</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>1</td>
      <td>0.007245</td>
      <td>0.115927</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>2</td>
      <td>0.006506</td>
      <td>0.130121</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>3</td>
      <td>0.008774</td>
      <td>0.131614</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>4</td>
      <td>0.011663</td>
      <td>0.233251</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>5</td>
      <td>0.008945</td>
      <td>0.178904</td>
    </tr>
  </tbody>
</table>


Plot correlations between libraries of the same selection for the site-level statistics:


```python
libraries = [lib for lib in effects_df['library'].unique() if lib != "average"]
assert len(libraries) == 2, 'plot only makes sense if 2 libraries'

for val in ['site_avg_escape_frac_single_mut', 'site_total_escape_frac_single_mut']:

    # wide data frame with each library's score in a different column
    site_effects_df_wide = (
        site_effects_df
        .query('library != "average"')
        # just get selections with 2 libraries
        .assign(nlibs=lambda x: x.groupby('selection')['library'].transform('nunique'))
        .query('nlibs == 2')
        # now make columns for each library, only keep sites with scores for both libs
        .pivot_table(index=['selection', 'site'],
                     columns='library',
                     values=val)
        .reset_index()
        .dropna(axis=0)
        )

    # correlations between libraries
    corrs = (
        site_effects_df_wide
        .groupby('selection')
        [libraries]
        .corr(method='pearson')
        .reset_index()
        .query('library == @libraries[0]')
        .assign(correlation=lambda x: 'R=' + x[libraries[1]].round(2).astype(str))
        [['selection', 'correlation']]
        # add number of mutations measured
        .merge(site_effects_df_wide
               .groupby('selection')
               .size()
               .rename('n')
               .reset_index()
               )
        .assign(correlation=lambda x: x['correlation'] + ', N=' + x['n'].astype(str))
        )

    # plot correlations
    nfacets = site_effects_df_wide['selection'].nunique()
    ncol = min(nfacets, 5)
    nrow = math.ceil(nfacets / ncol)
    xmin = site_effects_df_wide[libraries[0]].min()
    xspan = site_effects_df_wide[libraries[0]].max() - xmin
    ymin = site_effects_df_wide[libraries[1]].min()
    yspan = site_effects_df_wide[libraries[1]].max() - ymin
    p = (ggplot(site_effects_df_wide) +
         aes(libraries[0], libraries[1]) +
         geom_point(alpha=0.2) +
         geom_text(mapping=aes(label='correlation'),
                   data=corrs,
                   x=0.01 * xspan + xmin,
                   y=0.99 * yspan + ymin,
                   size=10,
                   ha='left',
                   va='top',
                   ) +
         facet_wrap('~ selection', ncol=ncol) +
         theme(figure_size=(2.5 * ncol, 2.5 * nrow),
               plot_title=element_text(size=14)) +
         ggtitle(val)
         )

    _ = p.draw()
```


    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_84_0.png)
    



    
![png](counts_to_scores_Omicron_BA286_files/counts_to_scores_Omicron_BA286_84_1.png)
    


## Write file with escape fractions at mutation and site levels
We write a files that has the mutation- and site-level escape fractions. This file has both the separate measurements for each library plus the average across libraries for all mutations measured in both libraries. We name the columns in such a way that this file can be used as [dms-view data file](https://dms-view.github.io/docs/dataupload):


```python
escape_fracs_to_write = (
    effects_df
    .merge(site_effects_df,
           how='left',
           validate='many_to_one',
           on=['selection', 'library', 'site'])
    .assign(protein_chain=config['escape_frac_protein_chain'],
            protein_site=lambda x: numpy.where(x['site'] >= 153,
                                               x['site'] + config['site_number_offset'] + 1,
                                               x['site'] + config['site_number_offset']),
            label_site=lambda x: x['protein_site'],
            condition=lambda x: x['selection'].where(x['library'] == "average", x['selection'] + '_' + x['library']),
            mutation=lambda x: x['mutant'],  # dms-view uses mutation to refer to mutant amino acid
            )
    [['selection', 'library', 'condition', 'site', 'label_site', 'wildtype', 'mutation',
      'protein_chain', 'protein_site', 'mut_escape_frac_single_mut', 'site_total_escape_frac_single_mut',
      'site_avg_escape_frac_single_mut', 'nlibs', 'n_single_mut_measurements',
      ]]
    .sort_values(['library', 'selection', 'site', 'mutation'])
    )

print('Here are the first few lines that will be written to the escape-fraction file:')
display(HTML(escape_fracs_to_write.head().to_html(index=False)))

print(f"\nWriting to {config['escape_fracs_Omicron_BA286']}")
escape_fracs_to_write.to_csv(config['escape_fracs_Omicron_BA286'], index=False, float_format='%.4g')

```

    Here are the first few lines that will be written to the escape-fraction file:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>selection</th>
      <th>library</th>
      <th>condition</th>
      <th>site</th>
      <th>label_site</th>
      <th>wildtype</th>
      <th>mutation</th>
      <th>protein_chain</th>
      <th>protein_site</th>
      <th>mut_escape_frac_single_mut</th>
      <th>site_total_escape_frac_single_mut</th>
      <th>site_avg_escape_frac_single_mut</th>
      <th>nlibs</th>
      <th>n_single_mut_measurements</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>S2V29_36</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>-</td>
      <td>E</td>
      <td>331</td>
      <td>0.004016</td>
      <td>0.115927</td>
      <td>0.007245</td>
      <td>2</td>
      <td>31</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>S2V29_36</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>A</td>
      <td>E</td>
      <td>331</td>
      <td>0.011084</td>
      <td>0.115927</td>
      <td>0.007245</td>
      <td>2</td>
      <td>20</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>S2V29_36</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>D</td>
      <td>E</td>
      <td>331</td>
      <td>0.004846</td>
      <td>0.115927</td>
      <td>0.007245</td>
      <td>2</td>
      <td>18</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>S2V29_36</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>E</td>
      <td>E</td>
      <td>331</td>
      <td>0.004962</td>
      <td>0.115927</td>
      <td>0.007245</td>
      <td>2</td>
      <td>29</td>
    </tr>
    <tr>
      <td>S2V29_36</td>
      <td>average</td>
      <td>S2V29_36</td>
      <td>1</td>
      <td>331</td>
      <td>N</td>
      <td>G</td>
      <td>E</td>
      <td>331</td>
      <td>0.008418</td>
      <td>0.115927</td>
      <td>0.007245</td>
      <td>2</td>
      <td>27</td>
    </tr>
  </tbody>
</table>


    
    Writing to results/escape_scores/escape_fracs_Omicron_BA286.csv



```python

```
