from pathlib import Path
import re
from itertools import groupby
import sys

from scipy import stats

SOLVE_DURATION_PREFIX = "SOLVER   Solved in"
SOLUTION_PREFIX = "SOLVER   Optimal objective"


def analyse_performance(paths_to_logs, path_to_output):
    with open(path_to_output, "w") as f_output:
        sys.stdout = f_output
        for parameter_setting, logs_for_parameter_setting in groupby(paths_to_logs, _parameter_setting):
            logs_for_parameter_setting = list(logs_for_parameter_setting)
            print(f"Parameter setting {parameter_setting}:")
            print(f"{len(logs_for_parameter_setting)} runs.")
            _analye_parameter_setting(logs_for_parameter_setting)
            print("")


def _parameter_setting(path_to_log):
    return Path(path_to_log).parent.name


def _analye_parameter_setting(paths_to_logs):
    results = [_preprocess_log(path) for path in paths_to_logs]
    successfull_results = list(filter(_found_solution, results))
    if len(successfull_results) == 0:
        print("All runs crashed.")
    elif len(successfull_results) != len(results):
        number_fails = len(results) - len(successfull_results)
        print(f"{number_fails} run(s) failed.")
        _print_stats(successfull_results)
    else:
        _print_stats(successfull_results)


def _print_stats(successfull_results):
    durations = [parse_solution_time(lines) for lines in successfull_results]
    solutions = [parse_solution(lines) for lines in successfull_results]
    iterations = [parse_iterations(lines) for lines in successfull_results]

    mean, std = stats.norm.fit(solutions)
    print(f"Gurobi found the following optimal value on average: {mean} (std {std}).")
    mean, std = stats.norm.fit(iterations)
    print(f"Gurobi needed on average {mean:.0f} (std {std:.0f}) iterations to find the optimal value.")
    mean, std = stats.norm.fit(durations)
    print(f"Gurobi needed on average {mean:.0f}s (std {std:.0f}s) to find the optimal value.")


def _preprocess_log(path_to_log):
    path_to_log = Path(path_to_log)
    with path_to_log.open("r") as log_file:
        lines = log_file.readlines()
    return [line[22:] for line in lines]


def _found_solution(lines):
    duration_found = any([SOLVE_DURATION_PREFIX in line for line in lines])
    solution_found = any([SOLUTION_PREFIX in line for line in lines])
    return duration_found & solution_found


def parse_solution_time(lines):
    duration_line = list(filter(lambda line: line.startswith(SOLVE_DURATION_PREFIX), lines))[0]
    return float(re.findall('\d+.?\d+', duration_line)[1])


def parse_iterations(lines):
    iterations_line = list(filter(lambda line: line.startswith(SOLVE_DURATION_PREFIX), lines))[0]
    return float(re.findall('\d+.?\d+', iterations_line)[0])


def parse_solution(lines):
    solution_line = list(filter(lambda line: line.startswith(SOLUTION_PREFIX), lines))[0]
    return float(re.findall('\d+.?\d+e?\+?\d+', solution_line)[0])


if __name__ == "__main__":
    analyse_performance(
        paths_to_logs=snakemake.input.logs,
        path_to_output=snakemake.output[0]
    )
