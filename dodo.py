import importlib.util
from pathlib import Path

from doit.action import CmdAction


def load_benchmark(benchmark_path):
    file_path = benchmark_path / "benchmark.py"
    module_name = benchmark_path.name + ".benchmark"
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BENCHMARKS = ["cpp", "rust"]


def task_build():
    for benchmark in BENCHMARKS:
        workdir = Path.cwd() / benchmark
        benchmark_file = load_benchmark(workdir)
        if not hasattr(benchmark_file, "BUILD_ACTIONS"):
            continue

        actions = benchmark_file.BUILD_ACTIONS

        yield {
            "name": benchmark,
            "actions": [CmdAction(a, cwd=workdir) for a in actions],
        }


def wrap_run(cmd, poetry=False):
    cmd = "hyperfine " + cmd
    if poetry:
        cmd = 'poetry run "{}"'.format(cmd)
    return cmd


def task_time():
    for benchmark in BENCHMARKS:
        workdir = Path.cwd() / benchmark
        benchmark_file = load_benchmark(workdir)
        if not hasattr(benchmark_file, "TIME_ACTIONS"):
            continue

        poetry = hasattr(benchmark_file, "USE_POETRY") and benchmark_file.USE_POETRY

        actions = benchmark_file.TIME_ACTIONS

        yield {
            "name": benchmark,
            "actions": [CmdAction(wrap_run(a, poetry), cwd=workdir) for a in actions],
            "task_dep": ["build:" + benchmark],
        }
