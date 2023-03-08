import glob
import math
import os
import random
from collections import defaultdict
from pathlib import Path
from statistics import fmean
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import pm4py
from pm4py import ProcessTree
from pm4py.util.lp import solver as lp_solver
from pm4py.objects.conversion.process_tree.converter import apply as pt_to_petri_net
from pm4py.objects.log.obj import EventLog, Trace, Event
from pm4py.objects.log.importer.xes.importer import apply as xes_import
from cortado_core.lca_approach import add_trace_to_pt_language
from cortado_core.models.infix_type import InfixType
from cortado_core.process_tree_utils.reduction import apply_reduction_rules
from cortado_core.utils.trace import TypedTrace

sns.set(style="ticks", rc={'figure.figsize': (5, 3.5)})

LOG_FILE = os.getenv('LOG_FILE', 'BPI_Challenge_2012_short.xes')
EVENT_LOG_DIRECTORY = os.getenv('EVENT_LOG_DIRECTORY', 'C:\\sources\\arbeit\\cortado\\event-logs')
REL_PERCENTAGE_INITIAL_MODEL = float(os.getenv('REL_SUPP_INIT_MODEL', '0.01'))
RESULTS_DIRECTORY = os.getenv('RESULT_DIRECTORY', 'C:\\sources\\arbeit\\cortado\\trace_fragments')
AVG_TRACE_LENGTH_TO_CUT = float(os.getenv('AVG_TL_TO_CUT', "0.2"))

COLUMNS = ['processed_variants', 'f_measure', 'fitness', 'precision', 'no_model_change', 'variant']
COLUMNS_CLASSICAL = ['processed_variants', 'f_measure', 'fitness', 'precision', 'no_model_change']


def run_experiments(log_filename: str, log_directory: str):
    log = xes_import(os.path.join(log_directory, log_filename))
    log = filter_log(log)
    log_path = os.path.join(RESULTS_DIRECTORY, f"{LOG_FILE}_{AVG_TRACE_LENGTH_TO_CUT}", 'logs')
    Path(log_path).mkdir(parents=True, exist_ok=True)
    pm4py.write_xes(log, os.path.join(log_path, 'filtered_log.xes'))
    traces_to_add = artificially_alter_log(log)
    __sort_typed_traces(traces_to_add)
    __export_typed_log(traces_to_add, os.path.join(log_path, 'preprocessed_log_with_fragments.xes'))
    variants_to_add, variants_dict = __get_sorted_variants(traces_to_add)
    initial_model, remaining_variants_to_add, already_added = __get_initial_model(variants_to_add,
                                                                                  REL_PERCENTAGE_INITIAL_MODEL)
    print('Starting experiments for trace fragment supporting ipd')
    __add_to_initial_model(initial_model, remaining_variants_to_add, already_added,
                           eval=lambda model: __eval_function(log, model), is_trace_fragments_experiment=True)
    print('Starting experiments for classical ipd')
    __add_to_initial_model(initial_model, remaining_variants_to_add, already_added,
                           eval=lambda model: __eval_function(log, model), is_trace_fragments_experiment=False)

    print('Starting experiments for classical discovery techniques')
    __run_experiments_for_classical_discovery_techniques(already_added, remaining_variants_to_add, variants_dict,
                                                         os.path.join(log_path, 'filtered_log.xes'))


def __add_to_initial_model(initial_model, variants_to_add, already_added, eval, is_trace_fragments_experiment):
    if is_trace_fragments_experiment:
        to_add = [v for v in variants_to_add]
    else:
        to_add = [TypedTrace(v.trace, InfixType.NOT_AN_INFIX) for v in variants_to_add]

    added = [v for v in already_added]
    evaluation_result = eval(initial_model)
    results = [[0] + evaluation_result]
    model = initial_model
    previous_evaluation = evaluation_result
    for idx, variant_to_add in enumerate(to_add):
        previous_model = model
        model = add_trace_to_pt_language(model, added, variant_to_add, try_pulling_lca_down=True,
                                         add_artificial_start_end=True)
        added.append(variant_to_add)

        if model == previous_model:
            results.append([idx + 1] + previous_evaluation + [True, __typed_trace_to_tuple(variant_to_add)])
        else:
            evaluation_result = eval(model)
            results.append([idx + 1] + evaluation_result + [False, __typed_trace_to_tuple(variant_to_add)])
            previous_evaluation = evaluation_result

    results_df = pd.DataFrame(results, columns=COLUMNS)
    results_path = os.path.join(RESULTS_DIRECTORY, f"{LOG_FILE}_{AVG_TRACE_LENGTH_TO_CUT}",
                                f'results_{"infix" if is_trace_fragments_experiment else "classical"}_preprocessed_logs_{REL_PERCENTAGE_INITIAL_MODEL}')
    Path(results_path).mkdir(parents=True, exist_ok=True)
    results_df.to_csv(os.path.join(results_path, 'data.csv'))
    create_plots(results_path, is_trace_fragments_experiment)


def filter_log(log: EventLog) -> EventLog:
    start_timestamp_to_filter, end_timestamp_to_filter = __get_start_end_timestamp_to_filter(log, 0.1)
    assert start_timestamp_to_filter < end_timestamp_to_filter
    log = __filter_log(log, lambda trace: __filter_by_start_timestamp(trace, start_timestamp_to_filter))
    log = __filter_log(log, lambda trace: __filter_by_end_timestamp(trace, end_timestamp_to_filter))

    return log


def __get_start_end_timestamp_to_filter(log: EventLog, percentage: float):
    timestamps = []
    for trace in log:
        for event in trace:
            timestamp = event['time:timestamp']
            timestamps.append(timestamp)

    timestamps.sort()
    abs_threshold = round(len(timestamps) * percentage)

    return timestamps[abs_threshold], timestamps[len(timestamps) - abs_threshold]


def __filter_log(log: EventLog, filter_func) -> EventLog:
    return EventLog([trace for trace in log if not filter_func(trace)], attributes=log.attributes,
                    properties=log.properties, extensions=log.extensions)


def __filter_by_start_timestamp(trace: Trace, start_timestamp) -> bool:
    for event in trace:
        if event['time:timestamp'] <= start_timestamp:
            return True

    return False


def __filter_by_end_timestamp(trace: Trace, end_timestamp) -> bool:
    for event in trace:
        if event['time:timestamp'] >= end_timestamp:
            return True

    return False


def artificially_alter_log(log: EventLog):
    traces = []
    avg_trace_length = fmean([len(trace) for trace in log])
    n_activities_to_remove = round(avg_trace_length * AVG_TRACE_LENGTH_TO_CUT)
    n_activities_to_remove = max(1, n_activities_to_remove)

    for trace in log:
        should_change_trace = random.random() >= 0.5
        if not should_change_trace:
            traces.append(TypedTrace(trace, InfixType.NOT_AN_INFIX))
            continue

        rnd = random.random()
        if rnd < 0.3333:
            new_trace = __trace_to_prefix(trace, n_activities_to_remove)
            if len(new_trace) > 0:
                traces.append(TypedTrace(new_trace, InfixType.PREFIX))
        elif 0.3333 <= rnd < 0.6666:
            new_trace = __trace_to_postfix(trace, n_activities_to_remove)
            if len(new_trace) > 0:
                traces.append(TypedTrace(new_trace, InfixType.POSTFIX))
        else:
            new_trace = __trace_to_infix(trace, n_activities_to_remove)
            if len(new_trace) > 0:
                traces.append(TypedTrace(new_trace, InfixType.PROPER_INFIX))

    assert len(traces) <= len(log)

    return traces


def __trace_to_prefix(trace: Trace, n_activities_to_remove: int) -> Trace:
    if len(trace) <= n_activities_to_remove:
        return Trace()

    return Trace(trace[:-n_activities_to_remove])


def __trace_to_postfix(trace: Trace, n_activities_to_remove: int) -> Trace:
    if len(trace) <= n_activities_to_remove:
        return Trace()

    return Trace(trace[n_activities_to_remove:])


def __trace_to_infix(trace: Trace, n_activities_to_remove: int) -> Trace:
    trace = __trace_to_prefix(trace, n_activities_to_remove)

    return __trace_to_postfix(trace, n_activities_to_remove)


def __export_typed_log(traces: list[TypedTrace], file: str):
    log = EventLog()

    for trace in traces:
        t = trace.trace
        for e in t:
            e['infix_type'] = trace.infix_type.value
        log.append(t)

    pm4py.write_xes(log, file)


def __sort_typed_traces(traces: list[TypedTrace]):
    variants = __get_variants(traces)

    traces.sort(key=lambda t: len(variants[__typed_trace_to_tuple(t)]), reverse=True)


def __typed_trace_to_tuple(trace: TypedTrace):
    return tuple([e['concept:name'] for e in trace.trace]), trace.infix_type


def __get_variants(traces: list[TypedTrace]):
    variants = defaultdict(list)

    for trace in traces:
        variants[__typed_trace_to_tuple(trace)].append(trace)

    return variants


def __get_sorted_variants(traces: list[TypedTrace]):
    variants_dict = __get_variants(traces)
    variants = sorted(variants_dict.keys(), key=lambda x: len(variants_dict[x]), reverse=True)
    sorted_variants = []

    for variant_tuple, infix_type in variants:
        t = Trace()
        for act in variant_tuple:
            e = Event()
            e['concept:name'] = act
            e['infix_type'] = infix_type
            t.append(e)

        sorted_variants.append(TypedTrace(t, infix_type))

    return sorted_variants, variants_dict


def __get_initial_model(variants: list[TypedTrace], rel_percentage: float):
    absolute_n_variants = math.ceil(
        rel_percentage * len([v for v in variants if v.infix_type == InfixType.NOT_AN_INFIX]))
    initial_variants = []
    for v in variants:
        if v.infix_type != InfixType.NOT_AN_INFIX:
            continue
        initial_variants.append(v)
        if len(initial_variants) >= absolute_n_variants:
            break

    for init_variant in initial_variants:
        variants.remove(init_variant)

    log = EventLog([v.trace for v in initial_variants])
    tree = pm4py.discover_process_tree_inductive(log, noise_threshold=0)
    apply_reduction_rules(tree)

    return tree, variants, initial_variants


def __run_experiments_for_classical_discovery_techniques(initial_model_variants, inc_added_variants, variants_dict,
                                                         eval_log_file):
    eval_log = xes_import(eval_log_file)
    partitioned_filenames = __partition_variants(initial_model_variants, inc_added_variants, variants_dict)
    for noise_threshold in [0, 0.5, 0.7, 0.9]:
        __discover_models_im(noise_threshold, partitioned_filenames,
                             res_path=os.path.join(RESULTS_DIRECTORY, f"{LOG_FILE}_{AVG_TRACE_LENGTH_TO_CUT}",
                                                   f'classical'))

    __evaluate_models(eval_log,
                      results_path=os.path.join(RESULTS_DIRECTORY, f"{LOG_FILE}_{AVG_TRACE_LENGTH_TO_CUT}",
                                                f'classical'))


def __partition_variants(initial_model_variants: list[TypedTrace], inc_added_variants: list[TypedTrace], variants_dict):
    partition_boundaries = [0.2, 0.4, 0.6, 0.8, 1]
    n_variants = len(inc_added_variants)
    filenames = []

    initial_model_traces = []
    for variant in initial_model_variants:
        initial_model_traces += variants_dict[__typed_trace_to_tuple(variant)]

    for boundary in partition_boundaries:
        threshold = math.ceil(n_variants * boundary)
        inc_variants = inc_added_variants[:threshold]

        variants_to_add = []
        for variant in inc_variants:
            variants_to_add += variants_dict[__typed_trace_to_tuple(variant)]

        filename = os.path.join(RESULTS_DIRECTORY, f"{LOG_FILE}_{AVG_TRACE_LENGTH_TO_CUT}", "logs", f"{boundary}.xes")
        filenames.append(filename)
        __export_typed_log(initial_model_traces + variants_to_add, filename)

    return filenames


def __evaluate_models(log, results_path=os.path.join(RESULTS_DIRECTORY, LOG_FILE, f'classical')):
    for folder, _, _ in os.walk(results_path):
        if folder == results_path:
            continue

        results = []
        for model_file in glob.glob(os.path.join(folder, "*.ptml")):
            percentage_added_traces = Path(model_file).stem
            model = pm4py.read_ptml(model_file)
            results.append([percentage_added_traces] + __eval_function(log, model) + [None])

        results_df = pd.DataFrame(results, columns=COLUMNS_CLASSICAL)
        results_path = os.path.join(folder, 'results.csv')
        results_df.to_csv(results_path)


def __discover_models_im(noise_threshold: float, logs: list[str],
                         res_path=os.path.join(RESULTS_DIRECTORY, LOG_FILE, f'classical')):
    for log_name in logs:
        log = xes_import(log_name)
        model = pm4py.discover_process_tree_inductive(log, noise_threshold=noise_threshold)
        filename = Path(log_name).stem
        results_path = os.path.join(res_path, f'im_{noise_threshold}')
        Path(results_path).mkdir(parents=True, exist_ok=True)
        pm4py.write_ptml(model, os.path.join(results_path, f'{filename}.ptml'))


def __eval_function(log, model):
    f, fitness, precision = __calculate_f_measure(model, log)
    return [f, fitness, precision]


def __calculate_f_measure(process_tree: ProcessTree, log: EventLog) -> tuple[float, float, float]:
    net, im, fm = pt_to_petri_net(process_tree)

    fitness = pm4py.fitness_alignments(log, net, im, fm)['averageFitness']
    precision = pm4py.precision_alignments(log, net, im, fm)

    return 2 * ((fitness * precision) / (fitness + precision)), fitness, precision


def create_plots(data_path: str, is_fragments_plots: bool = True):
    results_df = pd.read_csv(os.path.join(data_path, 'data.csv'))
    create_f_measure_fitness_precision_plots(data_path, results_df, is_fragments_plots)


def create_f_measure_fitness_precision_plots(filename: str, df: pd.DataFrame, is_fragments_plots: bool):
    Path(os.path.join(filename, 'plots')).mkdir(parents=True, exist_ok=True)
    df = df.drop_duplicates(subset='processed_variants', keep='last')
    create_measurement_plot(filename, 'f_measure', df)
    create_measurement_plot(filename, 'fitness', df)
    create_measurement_plot(filename, 'precision', df)
    create_combined_plot(filename, df, is_fragments_plots)


def create_measurement_plot(filename: str, measure: str, df: pd.DataFrame):
    sns.lineplot(data=df, x='processed_variants', y=measure)
    plt.ylim(0, 1.05)
    plt.savefig(os.path.join(filename, 'plots', 'variants_' + measure + '.pdf'),
                bbox_inches='tight')
    plt.close()


def create_combined_plot(filename: str, df: pd.DataFrame, is_fragments_plots: bool):
    df = df.rename(columns={"fitness": "Fitness", "precision": "Precision", 'f_measure': 'F-Measure'})
    plot_df = df.melt(id_vars=['processed_variants'], value_vars=['Fitness', 'F-Measure', 'Precision'])
    plot_df = plot_df.rename(columns={"variable": "Metric"})
    sns.lineplot(data=plot_df, x='processed_variants', y='value', hue='Metric')
    if is_fragments_plots:
        plt.xlabel('Added trace (fragment) variants')
    else:
        plt.xlabel('Added trace (fragment) variants')
    plt.ylabel('Metric')
    plt.ylim(0, 1.05)
    plt.savefig(os.path.join(filename, 'plots', 'combined_plot' + '.pdf'),
                bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    print('Experimental setup', 'new')
    print('LP Solver:', lp_solver.DEFAULT_LP_SOLVER_VARIANT)
    print('rel n activities to artificially cut', AVG_TRACE_LENGTH_TO_CUT)
    run_experiments(LOG_FILE, EVENT_LOG_DIRECTORY)
