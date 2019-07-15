import importlib.util
from pathlib import Path

from doit.action import CmdAction

DATA_DIR = Path.cwd() / "outputs"
RESULT_DIR = Path.cwd() / "results"


def load_benchmark(benchmark_path):
    file_path = benchmark_path / "benchmark.py"
    module_name = benchmark_path.name + ".benchmark"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BENCHMARKS = [
    folder
    for folder in Path(".").iterdir()
    if folder.is_dir() and (folder / "benchmark.py").exists()
]


def task_build():
    for workdir in BENCHMARKS:
        benchmark = workdir.name
        benchmark_file = load_benchmark(workdir)
        if not hasattr(benchmark_file, "BUILD_ACTIONS"):
            continue

        actions = benchmark_file.BUILD_ACTIONS

        yield {
            "name": benchmark,
            "actions": [CmdAction(a, cwd=workdir) for a in actions],
        }


def wrap_run(result, cmd, data, poetry=False):
    return 'hyperfine -w1 -r3 --export-json {} "{} {}"'.format(result, cmd, data)
    if poetry:
        cmd = "poetry run " + cmd
    return cmd


def task_time():
    for workdir in BENCHMARKS:
        benchmark = workdir.name
        data_file = (DATA_DIR / benchmark).with_suffix(".csv")
        result_file = (RESULT_DIR / benchmark).with_suffix(".json")

        benchmark_file = load_benchmark(workdir)
        if not hasattr(benchmark_file, "TIME_ACTION"):
            continue

        poetry = hasattr(benchmark_file, "USE_POETRY") and benchmark_file.USE_POETRY

        action = benchmark_file.TIME_ACTION

        yield {
            "name": benchmark,
            "actions": [
                CmdAction(wrap_run(result_file, action, data_file, poetry), cwd=workdir)
            ],
            "task_dep": ["build"],
            "targets": [data_file],
        }
