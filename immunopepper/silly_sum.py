""" a silly code example """

import sys

import click
import numpy as np


def silly_sum(arr, step=2):
    """
    Taking the sum of every n-th number in an array

    :param arr: a numerical array or array like (supporting striding)
    :param step: positive step size (striding size)
    :return: the sum
    """
    if step <= 0:
        raise ValueError(
            'Step parameter cannot be {}. Should be strictly positive.'.format(
                step))

    return arr[::step].sum()


def silly_sum_from_file(path, step):
    """
    Taking the sum every n-th number in a text file

    Convenience function.

    :param path: path to file containing one number per line
    :param step: positive step size
    :return: the sum
    """
    a = _load_array_from_file(path)
    return silly_sum(a, step)


def _load_array_from_file(file):
    return np.loadtxt(file)


@click.command()
@click.argument("file", type=click.Path())
@click.option('--step', type=int)
def main(file, step=2):
    result = silly_sum_from_file(file, step)
    print(result)
    return 0


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    sys.exit(main())
