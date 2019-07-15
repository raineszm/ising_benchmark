BUILD_ACTIONS = [
    "cmake -H. -B_builds -G 'Ninja' -DCMAKE_BUILD_TYPE=Release",
    "cmake --build _builds",
]

TIME_ACTIONS = ["_builds/ising"]
