import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

def get_formated_count_dataframes(df_, dicc_variants, time_period, omit_other=False):
    ### get all linages list
    all_linages = []
    for value in dicc_variants.values():
        all_linages += value

    ### get counts of grouped linages
    df_variants = pd.DataFrame({})
    for variant in dicc_variants.keys():
        mask = df_['pangolin_lineage'].str.contains('|'.join(dicc_variants[variant]))
        df_partial = df_.loc[mask][['count_']].resample(time_period, origin=df_.index[0]).sum()
        df_partial['Variant'] = variant
        df_variants = pd.concat([df_variants, df_partial])

    ### get counts of excluded linages
    if not omit_other:
        mask_ = df_['pangolin_lineage'].str.contains('|'.join(all_linages))
        df_partial = df_.loc[~mask_][['count_']].resample(time_period, origin=df_.index[0]).sum()
        df_partial['Variant'] = 'Otros'
        df_variants = pd.concat([df_variants, df_partial])

    ## return final df
    df_counts = df_variants.pivot(columns='Variant', values='count_').fillna(0).astype(int)
    
    return df_counts

def get_GISAID_bar_plot(df_counts, time_period = '14D', periods = 50, bar_width = 10, fontsize = 12, fontsize_dates = 10, title="Prevalencia por variantes - Bogotá"):
    ################################################
    df_prop = df_counts.div(df_counts.T.sum(), axis=0).fillna(0)
    ################################################
    date_list = df_prop.index
    variant_list = df_prop.columns
    ################################################
    plt.rcParams.update({'font.size': fontsize})
    fig, ax = plt.subplots()
    fig.set_size_inches(w=18, h=7)
    category_colors = plt.get_cmap('RdYlGn')(np.linspace(0.15, 0.85, len(variant_list)))
    P = []
    legend_label = []
    ################################################

    bottom = np.zeros(len(date_list))
    for color, variant in zip(category_colors, variant_list):
        props = df_prop[variant].to_numpy()
        counts = df_counts[variant].astype(str).str.replace('0','').to_numpy()

        p = ax.bar(date_list, props, bar_width, bottom=bottom, linewidth=1, edgecolor='black', color=color)
        P.append((p[0],))
        legend_label.append(variant)
        bottom += props

        for rect, count in zip(p, counts):
            x = rect.get_x() + rect.get_width() / 2 + 1
            y = rect.get_y() + rect.get_height() / 2
            ax.annotate(count, (x, y), xytext=(0, 0), textcoords="offset points", ha='center', va='center', fontsize=fontsize_dates, rotation=90)

    ax.set_title(title, fontsize=fontsize, loc='left')

    ### Plot details
    lower_date = date_list[-periods] - dt.timedelta(days=4)
    upper_date = date_list[-1] + dt.timedelta(days=4)

    time_labels = pd.date_range(date_list[-periods], periods=periods, freq=time_period)
    plt.xticks(time_labels, rotation=90) #, fontsize=fontsize_dates
    ax.set_xlim(lower_date, upper_date)

    ax.yaxis.set_ticks(np.arange(0,1.1,0.1)) #, fontsize=fontsize_dates
    ax.set_ylabel('Prevalencia por variante / Linaje viral')

    ## Legend
    ax.legend(P, legend_label, loc='upper center', bbox_to_anchor=(0.5, -0.2),
              fancybox=True, shadow=True, ncol=10, fontsize = 13.5)