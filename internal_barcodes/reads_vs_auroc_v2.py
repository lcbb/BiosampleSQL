import os
from pathlib import Path
import pickle
import random

import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

output_dir = f'tmp_{random.randint(0,10**6-1)}'
os.makedirs(output_dir)

plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.size'] = 12

sequencing_depth_values = [1000, 10000, 100000, 1000000]
log_sequencing_depth_values = [np.log10(depth) for depth in sequencing_depth_values]

def plot_auc_by_reads(replicate_count_paths, ground_truth_path):
    ground_truth_data = pd.read_csv(ground_truth_path)
    ground_truth_labels = ground_truth_data['Truth Value'].values
    
    auroc_scores_by_depth = [[] for _ in range(len(sequencing_depth_values))]
    
    for replicate_paths in replicate_count_paths:
        for depth_index, (data_file_path, current_sequencing_depth) in enumerate(zip(replicate_paths, sequencing_depth_values)):
            count_data = pd.read_csv(data_file_path)
            match_values = count_data['Match Count'].values
            
            false_positive_rates, true_positive_rates, _ = roc_curve(ground_truth_labels, match_values)
            area_under_roc_curve = auc(false_positive_rates, true_positive_rates)
            auroc_scores_by_depth[depth_index].append(area_under_roc_curve)
    
    plotting_data_list = []
    for depth_index, current_sequencing_depth in enumerate(sequencing_depth_values):
        for _, auroc_value in enumerate(auroc_scores_by_depth[depth_index]):
            plotting_data_list.append([np.log10(current_sequencing_depth), auroc_value])
    
    results_dataframe = pd.DataFrame(plotting_data_list, columns=['log_reads', 'auc'])
    
    figure_height = 3
    figure_width = figure_height
    plt.figure(figsize=(figure_width, figure_height))
    
    sns.barplot(x='log_reads', y='auc', data=results_dataframe, 
                errorbar=None, width=0.6, color='#E5E5E5')
    
    sns.swarmplot(x='log_reads', y='auc', data=results_dataframe, 
                  color='#4D4D4D', size=5)
    
    plt.xlabel('log$_{10}$(Number of reads)')
    plt.ylabel('AUROC')
    plt.ylim([0.7, 1.05]) 
    sns.despine()
    
    return results_dataframe

def plot_roc_curves(replicate_count_paths, ground_truth_path, dataset_key):
    ground_truth_data = pd.read_csv(ground_truth_path)
    ground_truth_labels = ground_truth_data['Truth Value'].values
    
    figure_height = 3
    figure_width = figure_height
    plt.figure(figsize=(figure_width, figure_height))
    
    curve_colors = plt.cm.viridis(np.linspace(0, 1, len(sequencing_depth_values)))
    
    for depth_index, current_sequencing_depth in enumerate(sequencing_depth_values):
        all_true_positive_rates = []
        reference_false_positive_rates = np.linspace(0, 1, 100)
        auroc_values_list = []
        
        for replicate_id, replicate_group_paths in enumerate(replicate_count_paths):
            data_file_path = replicate_group_paths[depth_index]
            count_data = pd.read_csv(data_file_path)
            match_values = count_data['Match Count'].values
            
            false_positive_rates, true_positive_rates, _ = roc_curve(ground_truth_labels, match_values)
            area_under_roc_curve = auc(false_positive_rates, true_positive_rates)
            auroc_values_list.append(area_under_roc_curve)
            
            interpolated_true_positive_rate = np.interp(reference_false_positive_rates, false_positive_rates, true_positive_rates)
            interpolated_true_positive_rate[0] = 0.0
            all_true_positive_rates.append(interpolated_true_positive_rate)
        
        average_true_positive_rate = np.mean(all_true_positive_rates, axis=0)
        average_true_positive_rate[-1] = 1.0
        std_dev_true_positive_rate = np.std(all_true_positive_rates, axis=0)
        average_auroc = np.mean(auroc_values_list)
        
        plt.plot(reference_false_positive_rates, average_true_positive_rate, color=curve_colors[depth_index],
                 label=f'{current_sequencing_depth:,} reads (AUC={average_auroc:.3f})')
        
        upper_confidence_tpr = np.minimum(average_true_positive_rate + std_dev_true_positive_rate, 1)
        lower_confidence_tpr = np.maximum(average_true_positive_rate - std_dev_true_positive_rate, 0)
        plt.fill_between(reference_false_positive_rates, lower_confidence_tpr, upper_confidence_tpr, color=curve_colors[depth_index], alpha=0.2)
    
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.8, label='Random')
    
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.xlim([0, 1])
    plt.ylim([0, 1.05])
    plt.legend(loc='lower right')
    
    sns.despine()
    
    plt.savefig(os.path.join(output_dir,f'roc_curves_{dataset_key}.pdf'))
    plt.close()
    
    return auroc_values_list

def calc_confusion_matrix(replicate_count_paths, ground_truth_path, dataset_key):
    ground_truth_data = pd.read_csv(ground_truth_path)
    ground_truth_labels = ground_truth_data['Truth Value'].values
    
    processed_ground_truth_labels = ground_truth_labels.copy()
    if dataset_key == 'not_symptomatic':
        processed_ground_truth_labels = 1 - ground_truth_labels
    
    actual_positives_count = np.sum(processed_ground_truth_labels)
    actual_negatives_count = len(processed_ground_truth_labels) - actual_positives_count
    total_sample_count = len(processed_ground_truth_labels)
    
    confusion_matrix_data_list = []
    
    for depth_index, current_sequencing_depth in enumerate(sequencing_depth_values):
        for replicate_id, replicate_group_paths in enumerate(replicate_count_paths):
            data_file_path = replicate_group_paths[depth_index]
            count_data = pd.read_csv(data_file_path)
            match_values = count_data['Match Count'].values
            
            classification_thresholds = np.sort(np.unique(match_values))
            optimal_f1_score = 0
            optimal_threshold_value = None
            optimal_performance_metrics = None
            
            for current_threshold in classification_thresholds:
                if dataset_key == 'not_symptomatic':
                    predicted_labels = (match_values <= current_threshold).astype(int)
                else:
                    predicted_labels = (match_values >= current_threshold).astype(int)
                
                true_positives = np.sum((predicted_labels == 1) & (processed_ground_truth_labels == 1))
                false_positives = np.sum((predicted_labels == 1) & (processed_ground_truth_labels == 0))
                true_negatives = np.sum((predicted_labels == 0) & (processed_ground_truth_labels == 0))
                false_negatives = np.sum((predicted_labels == 0) & (processed_ground_truth_labels == 1))
                
                precision_score = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
                recall_score = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
                f1_score_value = 2 * precision_score * recall_score / (precision_score + recall_score) if (precision_score + recall_score) > 0 else 0
                
                if f1_score_value > optimal_f1_score:
                    optimal_f1_score = f1_score_value
                    optimal_threshold_value = current_threshold
                    optimal_performance_metrics = (true_positives, false_positives, true_negatives, false_negatives, precision_score, recall_score, f1_score_value)
            
            true_positives, false_positives, true_negatives, false_negatives, precision_score, recall_score, f1_score_value = optimal_performance_metrics
            
            false_positive_rate_value = false_positives / actual_negatives_count if actual_negatives_count > 0 else 0
            accuracy_score_val = (true_positives + true_negatives) / len(processed_ground_truth_labels)
            
            confusion_matrix_data_list.append({
                'log_reads': log_sequencing_depth_values[depth_index],
                'read_count': current_sequencing_depth,
                'replicate': replicate_id,
                'threshold': optimal_threshold_value,
                'tp': true_positives,
                'fp': false_positives,
                'tn': true_negatives,
                'fn': false_negatives,
                'precision': precision_score,
                'recall': recall_score,
                'f1': f1_score_value,
                'fpr': false_positive_rate_value,
                'accuracy': accuracy_score_val,
                'tp_pct': true_positives / total_sample_count * 100,
                'fp_pct': false_positives / total_sample_count * 100,
                'tn_pct': true_negatives / total_sample_count * 100,
                'fn_pct': false_negatives / total_sample_count * 100
            })
    
    confusion_matrix_dataframe = pd.DataFrame(confusion_matrix_data_list)
    
    fig_cm, (axis_cm, axis_metrics) = plt.subplots(1, 2, figsize=(10, 5))
    
    max_depth_log_index = log_sequencing_depth_values.index(max(log_sequencing_depth_values))
    data_at_max_depth = confusion_matrix_dataframe[confusion_matrix_dataframe['log_reads'] == log_sequencing_depth_values[max_depth_log_index]]
    
    average_tp_percentage = data_at_max_depth['tp_pct'].mean()
    average_fp_percentage = data_at_max_depth['fp_pct'].mean()
    average_fn_percentage = data_at_max_depth['fn_pct'].mean()
    average_tn_percentage = data_at_max_depth['tn_pct'].mean()
    
    percentage_confusion_matrix = np.array([
        [average_tp_percentage, average_fn_percentage],
        [average_fp_percentage, average_tn_percentage]
    ])
    
    average_true_positives = data_at_max_depth['tp'].mean()
    average_false_positives = data_at_max_depth['fp'].mean()
    average_false_negatives = data_at_max_depth['fn'].mean()
    average_true_negatives = data_at_max_depth['tn'].mean()
    
    formatted_matrix_annotations = np.array([
        [f"{average_tp_percentage:.1f}%\n({average_true_positives:.0f})", f"{average_fn_percentage:.1f}%\n({average_false_negatives:.0f})"],
        [f"{average_fp_percentage:.1f}%\n({average_false_positives:.0f})", f"{average_tn_percentage:.1f}%\n({average_true_negatives:.0f})"]
    ])
    
    sns.heatmap(percentage_confusion_matrix, annot=formatted_matrix_annotations, fmt="", cmap='Blues', 
                         cbar=False, ax=axis_cm, xticklabels=['Target\npresent', 'Target\nabsent'], 
                         yticklabels=['Predicted\npresent', 'Predicted\nabsent'])
    
    if dataset_key == 'not_symptomatic':
        target_condition_label = "'NOT symptomatic'"
    elif dataset_key == 'immune':
        target_condition_label = "'Immunodeficient'"
    elif 'age' in dataset_key:
        if dataset_key == 'age1':
            target_condition_label = "'Age = 0'"
        elif dataset_key == 'age2':
            target_condition_label = "'15 ≤ Age ≤ 19'"
        elif dataset_key == 'age3':
            target_condition_label = "'50 ≤ Age ≤ 74'"
        else:
            target_condition_label = f"'{dataset_key}'"
    else:
        target_condition_label = f"'{dataset_key}'"
    
    axis_cm.set_title(f'Target: {target_condition_label}\n(at {sequencing_depth_values[max_depth_log_index]:,} reads)')
    
    mean_precision_score = data_at_max_depth['precision'].mean()
    mean_recall_score = data_at_max_depth['recall'].mean()
    mean_f1_score = data_at_max_depth['f1'].mean()
    mean_accuracy_score = data_at_max_depth['accuracy'].mean()
    
    performance_summary_text = (f"Precision: {mean_precision_score:.3f}\n"
                               f"Recall: {mean_recall_score:.3f}\n"
                               f"F1: {mean_f1_score:.3f}\n"
                               f"Accuracy: {mean_accuracy_score:.3f}")
    
    textbox_properties = dict(boxstyle='round', facecolor='white', alpha=0.5)
    axis_cm.text(0.5, -0.25, performance_summary_text, transform=axis_cm.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='center', bbox=textbox_properties)
    
    metrics_for_plotting = ['precision', 'recall', 'f1', 'accuracy']
    
    average_metrics_by_depth = confusion_matrix_dataframe.groupby('log_reads')[metrics_for_plotting].mean().reset_index()
    
    long_format_metrics_data = pd.melt(average_metrics_by_depth, id_vars='log_reads', 
                             value_vars=metrics_for_plotting, 
                             var_name='Metric', value_name='Value')
    
    sns.lineplot(x='log_reads', y='Value', hue='Metric', data=long_format_metrics_data, marker='o', ax=axis_metrics)
    axis_metrics.set_xlabel('log$_{10}$(Number of reads)')
    axis_metrics.set_ylabel('Metric value')
    axis_metrics.set_title('Performance metrics')
    axis_metrics.set_ylim([0, 1.05])
    axis_metrics.grid(True, alpha=0.3)
    
    target_class_prevalence = actual_positives_count / total_sample_count
    print(f"Target prevalence: {target_class_prevalence:.1%} ({actual_positives_count}/{total_sample_count})")
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'confusion_matrix_{dataset_key}.pdf'), dpi=300)
    plt.close()
    
    fig_table, axis_table = plt.subplots(figsize=(10, 3))
    axis_table.axis('tight')
    axis_table.axis('off')
    
    summary_table_data = confusion_matrix_dataframe.groupby('log_reads')[['tp', 'fp', 'tn', 'fn', 'precision', 'recall', 'f1', 'accuracy']].mean().reset_index()
    
    summary_table_data['tp_pct'] = (summary_table_data['tp'] / total_sample_count * 100).round(1)
    summary_table_data['fp_pct'] = (summary_table_data['fp'] / total_sample_count * 100).round(1)
    summary_table_data['tn_pct'] = (summary_table_data['tn'] / total_sample_count * 100).round(1)
    summary_table_data['fn_pct'] = (summary_table_data['fn'] / total_sample_count * 100).round(1)
    
    for col_name in ['precision', 'recall', 'f1', 'accuracy']:
        summary_table_data[col_name] = summary_table_data[col_name].round(3)
    
    summary_table_data['Reads'] = [f"{sequencing_depth_values[int(log_val-3)]}" for log_val in summary_table_data['log_reads']]
    
    table_column_order = ['Reads', 
                   'tp', 'tp_pct', 
                   'fn', 'fn_pct', 
                   'fp', 'fp_pct', 
                   'tn', 'tn_pct', 
                   'precision', 'recall', 'f1', 'accuracy']
    summary_table_data = summary_table_data[table_column_order]
    
    table_header_labels = ['Reads', 
                    'TP', 'TP %', 
                    'FN', 'FN %', 
                    'FP', 'FP %', 
                    'TN', 'TN %', 
                    'Precision', 'Recall', 'F1', 'Accuracy']
    
    matplotlib_table_object = axis_table.table(cellText=summary_table_data.values,
                     colLabels=table_header_labels,
                     cellLoc='center',
                     loc='center')
    
    matplotlib_table_object.auto_set_font_size(False)
    matplotlib_table_object.set_fontsize(9)
    matplotlib_table_object.scale(1, 1.5)
    
    plt.title(f"Retrieval Performance for {target_condition_label}", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'performance_table_{dataset_key}.pdf'), dpi=300)
    plt.close()
    
    return confusion_matrix_dataframe

def plot_confusion_counts_by_depth(confusion_matrix_dataframe, dataset_key, target_condition_label):
    figure_width = 3
    figure_height = 3 * 1.1
    fig_counts, axis_counts = plt.subplots(figsize=(figure_width, figure_height))
    
    first_row_sample = confusion_matrix_dataframe.iloc[0]
    actual_positives_count = first_row_sample['tp'] + first_row_sample['fn']
    actual_negatives_count = first_row_sample['tn'] + first_row_sample['fp']
    
    plotting_dataframe = confusion_matrix_dataframe.copy()
    plotting_dataframe['tp_pct'] = plotting_dataframe['tp'] / actual_positives_count * 100
    plotting_dataframe['fp_pct'] = plotting_dataframe['fp'] / actual_negatives_count * 100
    plotting_dataframe['fn_pct'] = plotting_dataframe['fn'] / actual_positives_count * 100
    plotting_dataframe['tn_pct'] = plotting_dataframe['tn'] / actual_negatives_count * 100
    
    metric_color_map = {'tp_pct': '#4daf4a', 'fp_pct': '#e41a1c', 'fn_pct': '#377eb8', 'tn_pct': '#984ea3'}
    replicate_marker_styles = ['o', 's', '^']
    metric_legend_labels = {'tp_pct': 'True positives', 
              'fp_pct': 'False positives', 
              'fn_pct': 'False negatives',
              'tn_pct': 'True negatives'}
    
    legend_handles_for_metrics = {}
    
    for metric_name in ['tp_pct', 'fp_pct', 'fn_pct', 'tn_pct']:
        for replicate_index_val in range(3):
            replicate_specific_data = plotting_dataframe[plotting_dataframe['replicate'] == replicate_index_val]
            plot_line_object, = axis_counts.plot(replicate_specific_data['log_reads'], replicate_specific_data[metric_name], marker=replicate_marker_styles[replicate_index_val], 
                    linestyle='-', color=metric_color_map[metric_name], alpha=0.7)
            
            if replicate_index_val == 0:
                legend_handles_for_metrics[metric_name] = plot_line_object
    
    legend_handles_for_replicates = []
    for replicate_index_val in range(3):
        plot_line_object, = axis_counts.plot([], [], marker=replicate_marker_styles[replicate_index_val], linestyle='', color='black')
        legend_handles_for_replicates.append(plot_line_object)
    
    combined_legend_handles = []
    combined_legend_labels = []
    
    for metric_name in ['tp_pct', 'fp_pct', 'fn_pct', 'tn_pct']:
        if metric_name in legend_handles_for_metrics:
            combined_legend_handles.append(legend_handles_for_metrics[metric_name])
            combined_legend_labels.append(metric_legend_labels[metric_name])
    
    combined_legend_handles.append(plt.Line2D([], [], linestyle='none'))
    combined_legend_labels.append('_' * 20)
    
    for replicate_index_val in range(3):
        combined_legend_handles.append(legend_handles_for_replicates[replicate_index_val])
        combined_legend_labels.append(f'Replicate {replicate_index_val+1}')
    
    axis_counts.set_xticks(log_sequencing_depth_values)
    axis_counts.set_xticklabels(log_sequencing_depth_values) 
    axis_counts.set_ylim([0.5, 105])
    
    axis_counts.spines['top'].set_visible(False)
    axis_counts.spines['right'].set_visible(False)
    
    axis_counts.set_xlabel('log$_{10}$(Number of reads)')
    axis_counts.set_ylabel('Percentage [%]')
    
    axis_counts.legend(combined_legend_handles, combined_legend_labels, loc='best')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,f'confusion_percentages_{dataset_key}.pdf'), dpi=300)
    plt.close()

def plot_performance_metrics_summary(confusion_matrix_dataframe, dataset_key, target_condition_label):
    confusion_matrix_dataframe['specificity'] = confusion_matrix_dataframe['tn'] / (confusion_matrix_dataframe['tn'] + confusion_matrix_dataframe['fp'])
    
    figure_height = 3
    figure_width = figure_height * 1.5
    
    fig_metrics, axis_metrics_plot = plt.subplots(figsize=(figure_width, figure_height))
    
    performance_metric_keys = ['precision', 'recall', 'accuracy', 'specificity']
    metric_bar_colors = ['#4daf4a', '#377eb8', '#e41a1c', '#984ea3']
    
    plotting_data_list = []
    for log_depth_value in log_sequencing_depth_values:
        data_for_current_depth = confusion_matrix_dataframe[confusion_matrix_dataframe['log_reads'] == log_depth_value]
        
        for metric_key in performance_metric_keys:
            for replicate_index_val in range(3): 
                replicate_specific_data = data_for_current_depth[data_for_current_depth['replicate'] == replicate_index_val]
                if len(replicate_specific_data) == 0:
                    continue
                
                plotting_data_list.append({
                    'log_reads': log_depth_value,
                    'metric': metric_key.capitalize(),
                    'value': replicate_specific_data[metric_key].values[0],
                    'replicate': replicate_index_val + 1
                })
    
    plotting_dataframe = pd.DataFrame(plotting_data_list)
    
    sns.barplot(x='log_reads', y='value', hue='metric', data=plotting_dataframe, 
                ax=axis_metrics_plot, alpha=0.7, palette=metric_bar_colors, errorbar=None)
    
    sns.swarmplot(x='log_reads', y='value', hue='metric', data=plotting_dataframe, 
                 dodge=True, ax=axis_metrics_plot, palette=metric_bar_colors, size=3, edgecolor='black')
    
    axis_metrics_plot.set_xlabel('log$_{10}$(Number of reads)')
    axis_metrics_plot.set_ylabel('Value')
    axis_metrics_plot.set_ylim([0, 1.05])
    
    axis_metrics_plot.spines['top'].set_visible(False)
    axis_metrics_plot.spines['right'].set_visible(False)
    
    legend_handles_from_plot, legend_labels_from_plot = axis_metrics_plot.get_legend_handles_labels()
    axis_metrics_plot.legend(legend_handles_from_plot[:len(performance_metric_keys)], legend_labels_from_plot[:len(performance_metric_keys)], loc='lower right')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir,f'performance_metrics_{dataset_key}.pdf'), dpi=300)
    plt.close()

def print_detailed_metrics(all_confusion_matrix_dataframes):
    performance_metric_keys_list = ['precision', 'recall', 'accuracy', 'specificity']
    confusion_matrix_count_keys = ['tp', 'fp', 'tn', 'fn']
    normalized_metric_keys = ['tp_pct', 'fp_pct', 'tn_pct', 'fn_pct']
    
    for dataset_id, current_confusion_matrix_df in all_confusion_matrix_dataframes.items():
        current_confusion_matrix_df['specificity'] = current_confusion_matrix_df['tn'] / (current_confusion_matrix_df['tn'] + current_confusion_matrix_df['fp'])
        
        print(f"\n{'='*140}")
        print(f"GROUP: {dataset_id}")
        print(f"{'='*140}")
        
        first_row_sample_data = current_confusion_matrix_df.iloc[0]
        actual_positives_count = first_row_sample_data['tp'] + first_row_sample_data['fn']
        actual_negatives_count = first_row_sample_data['tn'] + first_row_sample_data['fp']
        
        for depth_idx, current_sequencing_depth in enumerate(sequencing_depth_values):
            print(f"\n{'-'*140}")
            print(f"Read Depth: {current_sequencing_depth:,}")
            print(f"{'-'*140}")
            
            table_header_string = f"{'Replicate':<10}"
            for metric_key_header in confusion_matrix_count_keys:
                table_header_string += f"{metric_key_header.upper():<10}"
            for metric_key_header in normalized_metric_keys:
                table_header_string += f"{metric_key_header.upper().replace('_PCT', '%'):<10}"
            for metric_key_header in performance_metric_keys_list:
                table_header_string += f"{metric_key_header.capitalize():<15}"
            print(table_header_string)
            print('-' * len(table_header_string))
            
            data_for_current_depth = current_confusion_matrix_df[current_confusion_matrix_df['read_count'] == current_sequencing_depth]
            
            for replicate_index_val in range(3): 
                replicate_specific_data = data_for_current_depth[data_for_current_depth['replicate'] == replicate_index_val]
                if len(replicate_specific_data) == 0:
                    continue
                
                true_positives_val = replicate_specific_data['tp'].values[0]
                false_positives_val = replicate_specific_data['fp'].values[0]
                true_negatives_val = replicate_specific_data['tn'].values[0]
                false_negatives_val = replicate_specific_data['fn'].values[0]
                
                tp_percentage = (true_positives_val / actual_positives_count) * 100 if actual_positives_count > 0 else 0
                fp_percentage = (false_positives_val / actual_negatives_count) * 100 if actual_negatives_count > 0 else 0
                fn_percentage = (false_negatives_val / actual_positives_count) * 100 if actual_positives_count > 0 else 0
                tn_percentage = (true_negatives_val / actual_negatives_count) * 100 if actual_negatives_count > 0 else 0
                
                table_row_string = f"{replicate_index_val+1:<10}"
                
                table_row_string += f"{int(true_positives_val):<10}"
                table_row_string += f"{int(false_positives_val):<10}"
                table_row_string += f"{int(true_negatives_val):<10}"
                table_row_string += f"{int(false_negatives_val):<10}"
                
                table_row_string += f"{tp_percentage:.2f}%    "
                table_row_string += f"{fp_percentage:.2f}%    "
                table_row_string += f"{tn_percentage:.2f}%    "
                table_row_string += f"{fn_percentage:.2f}%    "
                
                for metric_key_print in performance_metric_keys_list:
                    table_row_string += f"{replicate_specific_data[metric_key_print].values[0]:.4f}       "
                print(table_row_string)

sequencing_data_paths_map = {
    'immune': [
    [
        r'match_analysis/20230808_immune_1_S3.assembled_match_counts_1k.csv',
        r'match_analysis/20230808_immune_1_S3.assembled_match_counts_10k.csv',
        r'match_analysis/20230808_immune_1_S3.assembled_match_counts.csv',
        r'match_analysis/20230808_immune_1_S3.assembled_match_counts_1M.csv'
    ],
    [
        r'match_analysis/20230808_immune_2_S4.assembled_match_counts_1k.csv',
        r'match_analysis/20230808_immune_2_S4.assembled_match_counts_10k.csv',
        r'match_analysis/20230808_immune_2_S4.assembled_match_counts.csv',
        r'match_analysis/20230808_immune_2_S4.assembled_match_counts_1M.csv'
    ],
    [
        r'match_analysis/20230808_immune_4_S6.assembled_match_counts_1k.csv',
        r'match_analysis/20230808_immune_4_S6.assembled_match_counts_10k.csv',
        r'match_analysis/20230808_immune_4_S6.assembled_match_counts.csv',
        r'match_analysis/20230808_immune_4_S6.assembled_match_counts_1M.csv'
    ]
    ],
    'not_symptomatic': [
    [
        r'match_analysis/2023_08_02_symptomatic_S2.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_02_symptomatic_S2.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_02_symptomatic_S2.assembled_match_counts.csv',
        r'match_analysis/2023_08_02_symptomatic_S2.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_symptomatic_1_S2.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_symptomatic_1_S2.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_symptomatic_1_S2.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_symptomatic_1_S2.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_symptomatic_2_S6.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_symptomatic_2_S6.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_symptomatic_2_S6.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_symptomatic_2_S6.assembled_match_counts_1M.csv',
    ]
    ],
    'age1': [
    [
        r'match_analysis/2023_08_02_age1_S3.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_02_age1_S3.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_02_age1_S3.assembled_match_counts.csv',
        r'match_analysis/2023_08_02_age1_S3.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age1_1_S3.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age1_1_S3.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age1_1_S3.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age1_1_S3.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age1_2_S7.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age1_2_S7.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age1_2_S7.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age1_2_S7.assembled_match_counts_1M.csv',
    ]
    ],
    'age2': [
    [
        r'match_analysis/2023_08_02_age2_S4.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_02_age2_S4.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_02_age2_S4.assembled_match_counts.csv',
        r'match_analysis/2023_08_02_age2_S4.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age2_1_S4.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age2_1_S4.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age2_1_S4.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age2_1_S4.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age2_2_S8.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age2_2_S8.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age2_2_S8.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age2_2_S8.assembled_match_counts_1M.csv',
    ]
    ],
    'age3': [
    [
        r'match_analysis/2023_08_02_age3_S5.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_02_age3_S5.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_02_age3_S5.assembled_match_counts.csv',
        r'match_analysis/2023_08_02_age3_S5.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age3_1_S5.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age3_1_S5.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age3_1_S5.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age3_1_S5.assembled_match_counts_1M.csv',
    ],
    [
        r'match_analysis/2023_08_03_age3_2_S9.assembled_match_counts_1k.csv',
        r'match_analysis/2023_08_03_age3_2_S9.assembled_match_counts_10k.csv',
        r'match_analysis/2023_08_03_age3_2_S9.assembled_match_counts.csv',
        r'match_analysis/2023_08_03_age3_2_S9.assembled_match_counts_1M.csv',
    ]
    ],
}

ground_truth_paths_map = {
  'immune': r'immune_truthtable.csv',
  'not_symptomatic': r'symp_neg_truthtable.csv',
  'age1': r'age1_truthtable.csv',
  'age2': r'age2_truthtable.csv',
  'age3': r'age3_truthtable.csv',
}

auroc_results_collection = {}
confusion_matrix_results_collection = {}  
roc_curve_data_collection = {}

for dataset_key_main, replicate_group_paths_main in sequencing_data_paths_map.items():
    print(f"Processing {dataset_key_main}...")
    
    auroc_dataframe_for_group = plot_auc_by_reads(replicate_group_paths_main, ground_truth_paths_map[dataset_key_main])
    auroc_results_collection[dataset_key_main] = auroc_dataframe_for_group
    
    plt.savefig(os.path.join(output_dir,f'auc_{dataset_key_main}.pdf'))
    plt.close()
    
    roc_auroc_scores_for_group = plot_roc_curves(replicate_group_paths_main, ground_truth_paths_map[dataset_key_main], dataset_key_main)
    roc_curve_data_collection[dataset_key_main] = roc_auroc_scores_for_group
    
    confusion_matrix_dataframe_for_group = calc_confusion_matrix(replicate_group_paths_main, ground_truth_paths_map[dataset_key_main], dataset_key_main)
    confusion_matrix_results_collection[dataset_key_main] = confusion_matrix_dataframe_for_group
    
    if dataset_key_main == 'not_symptomatic':
        target_condition_display_name = "'NOT symptomatic'"
    elif dataset_key_main == 'immune':
        target_condition_display_name = "'Immunodeficient'"
    elif 'age' in dataset_key_main:
        if dataset_key_main == 'age1':
            target_condition_display_name = "'Age = 0'"
        elif dataset_key_main == 'age2':
            target_condition_display_name = "'15 ≤ Age ≤ 19'"
        elif dataset_key_main == 'age3':
            target_condition_display_name = "'50 ≤ Age ≤ 74'"
        else:
            target_condition_display_name = f"'{dataset_key_main}'"
    else:
        target_condition_display_name = f"'{dataset_key_main}'"
    
    plot_confusion_counts_by_depth(confusion_matrix_dataframe_for_group, dataset_key_main, target_condition_display_name)
    
    plot_performance_metrics_summary(confusion_matrix_dataframe_for_group, dataset_key_main, target_condition_display_name)

print_detailed_metrics(confusion_matrix_results_collection)
