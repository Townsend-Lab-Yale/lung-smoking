import numpy as np
from itertools import product
from scipy.special import binom


# Numbers of possible positive lambdas (I doubt M would ever be >20)
numbers_positive_lambdas = [np.sum([binom(n, i)*i
                                    for i in range(n+1)]).astype(int)
                            for n in range(0, 20+1)]
"""List with number of possible positive lambdas per mutation number
M.

"""



def build_S_with_tuples(M):
    """Build entire space of mutations combinations S.

    It will represented by a list of size 2^M with items being tuples
    of size M with 0s and 1s.

    :type M: int
    :param M: Number of mutations.

    :rtype: list
    :return: S as an ordered list of tuples.

    """
    return [x for x in product([0, 1], repeat=M)]


def build_S_with_tuples_alt(M):
    """Build entire space of mutations combinations S.

    It will represented by a list of size 2^M with items being tuples
    of size M with 0s and 1s.

    :type M: int
    :param M: Number of mutations.

    :rtype: list
    :return: S as an ordered list of tuples.

    """
    if M == 1:
        return [(0,), (1,)]
    else:
        S_M_minus_1 = build_S_with_tuples_alt(M-1)
        S = ([x + (0,) for x in S_M_minus_1]
             + [x + (1,) for x in S_M_minus_1])
        return S


def build_S_as_array(M):
    """Build entire space of mutations combinations S.

    It will represented by a numpy array of shape (2^M, M).

    :type M: int
    :param M: Number of mutations.

    :rtype S: numpy.ndarray
    :return: S as a numpy array of 0s and 1s.

    """
    return np.array(build_S_with_tuples(M))


def obtain_pos_lambdas_indices(S):
    """Obtain the lambda indices that could be positive.

    If Q is a matrix with entries its i-th row, j-th column being
    lambda[S[j], S[i]], then this function would produce a matrix with
    i-th row, j-th column being True if it is posible for lambda[S[j],
    S[i]] > 0, just by the definition of the lambdas, in other words
    if S[i] is equal S[j] plus one more mutation.

    :type S: numpy.ndarray
    :param S: S as a numpy array of 0s and 1s.

    :rtype: numpy.ndarray
    :return: A numpy array of booleans.

    """

    M = np.shape(S)[1]

    pos_lambdas_indices = np.array(
        [[True
          if np.sum(S[j] >= S[i]) == M
          and np.sum(S[j] - S[i]) == 1
          and i < len(S)-1 # for last x=(1,...,1) all are 0
          else False
          for i in range(len(S))]
         for j in range(len(S))])

    return pos_lambdas_indices


def order_pos_lambdas(S):
    """Order the tuples in S that subscript the positive lambdas.

    This function can be used to translate a flatten array of positive
    lambdas (such as the one returned by the MCMC or MAP finding
    routine) to a dictionary indexed by two elements of S.
    If pos_lambdas_ordered represents the ordered list returned by
    this function, lambda[x, y] is the flux from x to y, and Q is a
    matrix with i-th row j-th column being lambda[S[j], S[i]] then:

        Q[pos_lambdas_ordered]

    is the same as

        [lambdas[build_S_with_tuples(M).index(j),
                 build_S_with_tuples(M).index(i)]
         for i, j in positive_lambdas_ordered]


    :type S: numpy.ndarray
    :param S: S as a numpy array of 0s and 1s.

    :rtype: list
    :return: A list with the lambda tuples numpy array of booleans.

    """

    pos_lambdas_indices = obtain_pos_lambdas_indices(S)
    pos_lambdas_ordered = [(tuple(S[j]), tuple(S[i]))
                           for i in range(len(S))
                           for j in range(len(S))
                           if pos_lambdas_indices[i, j]]
    return pos_lambdas_ordered


def human_order_single_jumps(n, first_by_destiny=False):
    if not first_by_destiny:
        pos_lambdas = order_pos_lambdas(build_S_as_array(n))
        ordered = sorted(
            pos_lambdas,
            key=lambda xy: (np.sum(xy[0])*10**(2*n) +
                            np.dot(xy[0], [10**i for i in range(n, 2*n)]) +
                            np.dot(xy[1], [10**i for i in range(n)])))
    else:
        ordered = []
        for i in range(n):
            ordered = (
                ordered
                + [xy for xy in human_order_single_jumps(n, False)
                   if tuple(np.array(xy[1])-np.array(xy[0])).index(1) == i])
    return ordered


def human_order_paths(n):
    all_paths = []
    for path in product(*[[x for x in build_S_with_tuples(n)
                           if np.sum(x) == i]
                          for i in range(n+1)]):
        valid_path = True
        k = 1
        while valid_path and k < n:
            y_minus_x = np.array(path[k+1])-np.array(path[k])
            if not np.sum(np.power(y_minus_x, 2)) == 1:
                valid_path = False
            k += 1
        if valid_path:
            all_paths.append(path)
    return all_paths[::-1]


def generate_paths(target):
    """For a final `target' somatic genotype destination give all
    possible paths of somatic genotypes that lead there.

    :type target: tuple
    :param target: A destination somatic genotype as a tuple.

    :rtype: list
    :return: A list with the paths. Each path is represented by a list
        that has pairs representing each jump in path, for example,
        [(x1,x2), (x2, x3), (x3, x4)] would represent a path to the
        target x4, where the first jump is from x1 to x2, then to x3
        and finally to x4. The first entry of the first tuple in each
        path (x1, in the example) is always the normal genotype.

    """
    M = len(target)
    initial_state = tuple([0] * M)

    # Function to generate next states from a given state
    def next_states(state):
        return [state[:i] + (1,) + state[i+1:]
                for i in range(M)
                if state[i] == 0]

    # Recursive function to generate paths
    def recursive_paths(state):
        if state == target:
            return [[state]]
        paths = []
        for next_state in next_states(state):
            for path in recursive_paths(next_state):
                paths.append([state] + path)
        return paths

    all_paths = recursive_paths(initial_state)

    # Organize paths as pairs
    organized_paths = []
    for path in all_paths:
        organized_paths.append(list(zip(path[:-1], path[1:])))

    return organized_paths
