project('ising_cpp', 'cpp', 'c',
        default_options : ['cpp_std=c++11',
        'buildtype=release'])
threaddep = dependency('threads')
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
executable('ising', 'src/ising.cpp', dependencies : [m_dep, threaddep])
