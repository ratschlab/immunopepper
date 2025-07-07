import immunopepper.utils as utils

def test_get_all_comb():
    assert [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)] == utils.get_all_comb([1,2,3])
    assert [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3)] == utils.get_all_comb([1,2,3], 2)


