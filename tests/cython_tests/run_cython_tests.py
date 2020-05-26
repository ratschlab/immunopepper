
#! python tests/cython_tests/setup.py build_ext --inplace
import timeit


cy = timeit.timeit('''cpython_test.translate_dna_to_peptide('ATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGG')''',setup='import cpython_test',number=100)
py = timeit.timeit('''python_test.translate_dna_to_peptide('ATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGG')''',setup='import python_test', number=100)

print(cy, py)
print('Cyttests/cython_tests/run_cython_tests.py:7hotests/cython_tests/run_cython_tests.py:7n is {}x faster'.format(py/cy))




import pyximport; pyximport.install()
import cpython_test
import python_test
start_time = timeit.default_timer()
cpython_test.translate_dna_to_peptide('ATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGG')
print('cpython single test', timeit.default_timer() - start_time)

start_time = timeit.default_timer()
python_test.translate_dna_to_peptide('ATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGGATCATTCTTCCACATCAGGTACGG')
print('python single test', timeit.default_timer() - start_time)