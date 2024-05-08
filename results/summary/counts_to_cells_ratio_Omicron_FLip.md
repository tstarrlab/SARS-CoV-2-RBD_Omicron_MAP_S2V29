# Ratio of sequencing counts to cells sorted
This Python Jupyter notebook looks at the ratio of sequencing counts to cells sorted for each sample, and flags any where this ratio seems too low.

First, import Python modules:


```python
import Bio.SeqIO

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml

%matplotlib inline
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Read information about the samples:


```python
samples_df = pd.read_csv(config['barcode_runs_Omicron_FLip'])
```

Read the variant counts:


```python
variant_counts = pd.read_csv(config['variant_counts_Omicron_FLip'])
```

Get the "primary target"


```python
wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence_Omicron_FLip'], 'fasta')
primary_target = wt_seqrecord.name
print(f"Primary target is {primary_target}")
```

    Primary target is Omicron_FLip


Merge the sample information and aggregated variant counts data frames into a new data frame that has the number of cells sorted and the total variant counts for each sample, as well as the ratio of the variant counts to cells.
Importantly, counts and ratios are computed using counts for all variants (SARS-CoV-2) and counts for just variants of the primary target (SARS-CoV-2).
The cell counts aren't meaningful for reference samples, so set to `NA`.


```python
counts_cells = (
    samples_df
    .assign(sample_lib=lambda x: x['sample'] + ', ' + x['library'],
            number_cells=lambda x: x['number_cells'].where(x['sort_bin'] != 'ref', pd.NA))
    [['sample_lib', 'date', 'experiment', 'antibody', 'concentration', 'sort_bin', 'number_cells']]
    .merge(variant_counts
           .assign(sample_lib=lambda x: x['sample'] + ', ' + x['library'])
           .groupby('sample_lib')
           .aggregate(counts=pd.NamedAgg('count', 'sum'))
           .astype(float)
           .reset_index(),
           on='sample_lib', how='left', validate='one_to_one',
           )
    .assign(sample_lib=lambda x: pd.Categorical(x['sample_lib'], reversed(x['sample_lib'].unique()), ordered=True),
            counts_to_cells=lambda x: x['counts'] / x['number_cells'],
            is_reference=lambda x: x['sort_bin'].isin(['trans', 'ref']),
           )
    .rename(columns={'counts_to_cells': 'counts:cells',
                     })
    )

print(f"First few lines of the data frame, writing entirety to {config['counts_to_cells_csv_Omicron_FLip']}:")
display(HTML(counts_cells.head().to_html(index=False, float_format='{:.2g}'.format)))
counts_cells.to_csv(config['counts_to_cells_csv_Omicron_FLip'], index=False, float_format='%.3g')
```

    First few lines of the data frame, writing entirety to results/counts/Omicron_FLip/counts_to_cells_csv.csv:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>sample_lib</th>
      <th>date</th>
      <th>experiment</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>sort_bin</th>
      <th>number_cells</th>
      <th>counts</th>
      <th>counts:cells</th>
      <th>is_reference</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>exptREF-none-0-ref, lib96</td>
      <td>230615</td>
      <td>exptREF</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>NaN</td>
      <td>6.6e+07</td>
      <td>NaN</td>
      <td>True</td>
    </tr>
    <tr>
      <td>exptREF-none-0-ref, lib97</td>
      <td>230615</td>
      <td>exptREF</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>NaN</td>
      <td>6.5e+07</td>
      <td>NaN</td>
      <td>True</td>
    </tr>
    <tr>
      <td>expt4-S2V29_v37_2-89-abneg, lib96</td>
      <td>230615</td>
      <td>expt4</td>
      <td>S2V29_v37_2</td>
      <td>89</td>
      <td>abneg</td>
      <td>4.4e+05</td>
      <td>2.4e+06</td>
      <td>5.5</td>
      <td>False</td>
    </tr>
    <tr>
      <td>expt4-S2V29_v37_2-89-abneg, lib97</td>
      <td>230615</td>
      <td>expt4</td>
      <td>S2V29_v37_2</td>
      <td>89</td>
      <td>abneg</td>
      <td>3.8e+05</td>
      <td>2.2e+06</td>
      <td>5.8</td>
      <td>False</td>
    </tr>
    <tr>
      <td>expt5-S2V29-290-abneg, lib96</td>
      <td>230615</td>
      <td>expt5</td>
      <td>S2V29</td>
      <td>290</td>
      <td>abneg</td>
      <td>9.2e+05</td>
      <td>5.2e+06</td>
      <td>5.7</td>
      <td>False</td>
    </tr>
  </tbody>
</table>


Now we plot the number of counts or the counts:cells ratio for all samples.
We stratify by reference samples and escape samples, and only plot the counts:cells ratio for escape samples as cells are not sorted for reference samples.
We plot the counts / ratios for all variants, and primary target (SARS-CoV-2) variants.

In order to have the axis limits be reasonable, we clip very high / low values and draw dotted black lines to indicate the clipping.
For reference sample counts, and escape sample counts:cells ratio, we have minimum desired values.
We draw dashed green lines at these minimum values, and keep track of samples that don't achieve the minimum.


```python
min_fail_sample_libs = set([])  # sample_lib that fail desired min

for var, is_reference, lower_clip, upper_clip, desired_min in [
            ('counts', True, 1e5, None, config['reference_min_counts']),
            ('counts', False, 1e5, None, None),
            ('counts:cells', False, 0.1, 10, config['min_counts_to_cells_ratio']),
            ]:

    # get tidy data frame with samples of interest
    tidy_df = (
        counts_cells
        .query('is_reference == @is_reference')
        )
    tidy_df[var] = tidy_df[var].clip(lower=lower_clip, upper=upper_clip).astype(float)

    # make plot
    p = (ggplot(tidy_df) +
         aes(var, 'sample_lib') +
         geom_point(size=2, alpha=0.7) +
         theme(figure_size=(4, 0.25 * tidy_df['sample_lib'].nunique())) +
         ggtitle(f"{var} for reference=={is_reference} samples")
         )
    if var == 'counts':  # plot counts on log scale
        p = p + scale_x_log10()
        
    # add dotted vertical lines if clipping if data
    if (lower_clip is not None) and (lower_clip >= tidy_df[var].min()):
        p = p + geom_vline(xintercept=lower_clip, color='black', linetype='dotted')
    if (upper_clip is not None) and (upper_clip <= tidy_df[var].max()):
        p = p + geom_vline(xintercept=upper_clip, color='black', linetype='dotted')
        
    # draw line at desired minimum value, and identify any samples that fail minimum
    if desired_min is not None:
        p = p + geom_vline(xintercept=desired_min, color='green', linetype='dashed')
        min_fail_sample_libs.update(set(
            tidy_df
            .assign(fail_min=lambda x: x[var] < desired_min)
            .query('fail_min')
            ['sample_lib']
            ))
    
    # draw figure
    fig = p.draw()
    display(fig)
    plt.close(fig)
```

    /scratch/local/u6042467/13245058/ipykernel_2717971/1580948032.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy



    
![png](counts_to_cells_ratio_Omicron_FLip_files/counts_to_cells_ratio_Omicron_FLip_13_1.png)
    


    /scratch/local/u6042467/13245058/ipykernel_2717971/1580948032.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy



    
![png](counts_to_cells_ratio_Omicron_FLip_files/counts_to_cells_ratio_Omicron_FLip_13_3.png)
    


    /scratch/local/u6042467/13245058/ipykernel_2717971/1580948032.py:14: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy



    
![png](counts_to_cells_ratio_Omicron_FLip_files/counts_to_cells_ratio_Omicron_FLip_13_5.png)
    


Now list any samples that fail the minimum counts or counts:cell ratio:


```python
print(f"Reference samples with < {config['reference_min_counts']} counts, "
      f"or escape samples with a counts:cell ratio < {config['min_counts_to_cells_ratio']}.")

display(HTML(
    counts_cells
    .query('sample_lib in @min_fail_sample_libs')
    .to_html(index=False, float_format='{:.2g}'.format)
    ))
```

    Reference samples with < 25000000.0 counts, or escape samples with a counts:cell ratio < 2.5.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>sample_lib</th>
      <th>date</th>
      <th>experiment</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>sort_bin</th>
      <th>number_cells</th>
      <th>counts</th>
      <th>counts:cells</th>
      <th>is_reference</th>
    </tr>
  </thead>
  <tbody>
  </tbody>
</table>



```python

```
