from libc.math cimport log10

# https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
cdef union struct_union:
    double d
    int x[2]

cdef inline double fastPow(double a, double b):
    cdef struct_union u

    u = struct_union(a)
    u.x[1] = (int)(b * (u.x[1] - 1072632447) + 1072632447)
    u.x[0] = 0
    return u.d

def compute_average_quality(list phred_list, long l):
    cdef double q = 0.0
    cdef double prob = 0.0
    cdef double average_phred = 0.0

    for q in phred_list:
        prob += fastPow(10, q / -10)

    average_phred = -10 * log10(prob / l)

    return average_phred

# python setup.py build_ext --inplace
