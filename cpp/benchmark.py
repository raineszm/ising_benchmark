BUILD_ACTIONS = [
    "cmake -H. -B_builds -G 'Ninja' -DCMAKE_BUILD_TYPE=Release",
    "cmake --build _builds",
]

TIME_ACTION = "_builds/ising"
