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


BENCHMARKS = ["cpp"]


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
