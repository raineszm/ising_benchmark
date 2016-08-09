project('ising_cpp', 'cpp', 'c',
        default_options : ['cpp_std=c++11',
        'buildtype=release'])
threaddep = dependency('threads')
executable('ising', 'src/ising.cpp', dependencies : threaddep)