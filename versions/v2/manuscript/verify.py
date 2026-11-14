import argparse
import json
from math import comb
from pathlib import Path


def partitions(weight, ceiling=None):
    if weight == 0:
        yield ()
        return
    ceiling = weight if ceiling is None else min(ceiling, weight)
    for first in range(ceiling, 0, -1):
        for tail in partitions(weight - first, first):
            yield (first,) + tail


def reverse_sum(partition):
    endpoint = 0
    for part in partition:
        endpoint = part - endpoint
    return endpoint


def alternating_formula(partition):
    return sum((-1) ** (len(partition) - index - 1) * part
               for index, part in enumerate(partition))


def paired(partition):
    return len(partition) % 2 == 0 and all(
        partition[index] == partition[index + 1]
        for index in range(0, len(partition), 2))


def raise_partition(partition):
    if len(partition) % 2:
        return (partition[0] + 1,) + partition[1:]
    return partition + (1,)


def negative_raise(partition):
    return (partition[0] + 1,) + partition[1:]


def conjugate(partition):
    if not partition:
        return ()
    return tuple(sum(part >= height for part in partition)
                 for height in range(1, partition[0] + 1))


def parity_counts(maximum, distinct=False):
    counts = [[0, 0] for weight in range(maximum + 1)]
    counts[0][0] = 1
    for part in range(1, maximum + 1):
        weights = range(maximum, part - 1, -1) if distinct else range(part, maximum + 1)
        for weight in weights:
            even, odd = counts[weight - part]
            counts[weight][1] += even
            counts[weight][0] += odd
    return counts


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def compositions(weight):
    if weight == 0:
        yield ()
        return
    for first in range(1, weight + 1):
        for tail in compositions(weight - first):
            yield (first,) + tail


def ordered_checks(maximum=16):
    expected_positive = [1, 1, 2, 3, 7, 11, 26, 42, 99, 163, 382, 638, 1486,
                         2510, 5812, 9908, 22819]
    expected_negative_shifted = [0, 0, 1, 1, 5, 6, 22, 29, 93, 130, 386, 562,
                                 1586, 2380, 6476, 9949]
    checked = 0
    table = []
    for weight in range(maximum + 1):
        positive_count = 0
        distribution = {}
        for composition in compositions(weight):
            checked += 1
            endpoint = reverse_sum(composition)
            forward = sum((-1) ** index * part for index, part in enumerate(composition))
            require(endpoint == sum((-1) ** index * part for index, part in enumerate(reversed(composition))),
                    'reversal exchanges the two alternating sums')
            positive_count += endpoint >= 0
            distribution[forward] = distribution.get(forward, 0) + 1
        formula = 1 if weight == 0 else sum(comb(weight - 1, index)
                                           for index in range(weight // 2 + 1))
        require(positive_count == formula == expected_positive[weight], 'ordered nonnegative formula')
        total = 1 if weight == 0 else 2 ** (weight - 1)
        negative_count = total - positive_count
        if weight:
            require(negative_count == expected_negative_shifted[weight - 1], 'A294175 shift')
            for endpoint in range(-weight, weight + 1):
                lower_twice = weight + endpoint - 2
                coefficient = (comb(weight - 1, lower_twice // 2)
                               if lower_twice % 2 == 0 and 0 <= lower_twice // 2 < weight else 0)
                require(distribution.get(endpoint, 0) == coefficient, 'bivariate coefficient')
        if weight and weight % 2 == 0:
            half = weight // 2
            require(positive_count == 4 ** (half - 1) + comb(2 * half, half) // 2,
                    'A114121 bisection')
            require(negative_count == 4 ** (half - 1) - comb(2 * half - 1, half - 1),
                    'A008549 shifted bisection')
        table.append({'weight': weight, 'nonnegative': positive_count, 'negative': negative_count})
    values = [1] + [sum(comb(weight - 1, index) for index in range(weight // 2 + 1))
                    for weight in range(1, 101)]
    for weight in range(3, 101):
        recurrence = (weight * (weight - 3) * values[weight]
                      + 2 * (-weight ** 2 + 4 * weight - 2) * values[weight - 1]
                      - 4 * (weight - 2) ** 2 * values[weight - 2]
                      + 8 * (weight - 2) * (weight - 3) * values[weight - 3])
        require(recurrence == 0, 'Mathar recurrence')
    nonnegative_paths = 0
    for encoding in range(16):
        height = 0
        minimum = 0
        for bit in range(4):
            height += 1 if encoding & (1 << bit) else -1
            minimum = min(minimum, height)
        nonnegative_paths += minimum >= 0
    require(nonnegative_paths == 6 and values[4] == 7, 'endpoint versus prefix path condition')
    return {'weights': [0, maximum], 'checked_compositions': checked, 'counts': table,
            'recurrence_indices': [3, 100],
            'path_condition_counterexample': {'weight': 4, 'H': 7, 'nonnegative_paths': 6},
            'complement_counterexample': {'weight': 4, 'H': 7, 'unshifted_expression': 3},
            'even_complement_counterexample': {'index': 2, 'a114121': 7, 'unshifted_expression': 2}}


def run(maximum_n):
    maximum_weight = 2 * maximum_n + 2
    ordinary = parity_counts(maximum_weight)
    distinct = parity_counts(maximum_weight, distinct=True)
    counts = []
    endpoint_distributions = []
    distinct_distributions = []
    checked_partitions = 0
    checked_nonnegative_maps = 0
    checked_negative_maps = 0
    for weight in range(maximum_weight + 1):
        row = {'weight': weight, 'partitions': 0, 'nonnegative': 0,
               'negative': 0, 'odd_length': 0, 'even_length': 0}
        nonnegative_images = set()
        negative_images = set()
        strict_positive = 0
        strict_negative = 0
        endpoint_distribution = {}
        distinct_distribution = {}
        for partition in partitions(weight):
            checked_partitions += 1
            row['partitions'] += 1
            odd = len(partition) % 2
            row['odd_length' if odd else 'even_length'] += 1
            endpoint = reverse_sum(partition)
            endpoint_distribution[endpoint] = endpoint_distribution.get(endpoint, 0) + 1
            require(endpoint == alternating_formula(partition), 'endpoint recurrence')
            require((endpoint >= 0) == (bool(odd) or paired(partition)), 'sign characterization')
            transposed = conjugate(partition)
            if partition:
                require(endpoint == (-1) ** (len(partition) - 1) * sum(part % 2 for part in transposed),
                        'exact conjugate odd-part formula')
            if len(set(partition)) == len(partition):
                strict_positive += endpoint > 0
                strict_negative += endpoint < 0
                distinct_distribution[endpoint] = distinct_distribution.get(endpoint, 0) + 1
            conjugate_condition = (bool(transposed) and transposed[0] % 2 == 1) or all(
                part % 2 == 0 for part in transposed)
            require((endpoint >= 0) == conjugate_condition, 'conjugation characterization')
            negative_conjugate = (bool(transposed) and transposed[0] % 2 == 0
                                  and any(part % 2 == 1 for part in transposed))
            require((endpoint < 0) == negative_conjugate, 'negative conjugation characterization')
            row['nonnegative' if endpoint >= 0 else 'negative'] += 1
            if weight == maximum_weight:
                continue
            if endpoint >= 0:
                image = raise_partition(partition)
                require(sum(image) == weight + 1 and len(image) % 2 == 1, 'target weight and length')
                require(all(image[index] >= image[index + 1]
                            for index in range(len(image) - 1)), 'target order')
                require(image not in nonnegative_images, 'nonnegative injection collision')
                nonnegative_images.add(image)
                checked_nonnegative_maps += 1
            else:
                image = negative_raise(partition)
                require(reverse_sum(image) < 0 and sum(image) == weight + 1, 'negative map')
                require(image not in negative_images, 'negative injection collision')
                negative_images.add(image)
                checked_negative_maps += 1
        require(row['even_length'] == ordinary[weight][0], 'even-length dynamic program')
        require(row['odd_length'] == ordinary[weight][1], 'odd-length dynamic program')
        correction = sum(ordinary[weight // 2]) if weight % 2 == 0 else 0
        require(row['nonnegative'] == ordinary[weight][1] + correction, 'count identity C')
        require(row['negative'] == ordinary[weight][0] - correction, 'count identity N')
        require(strict_positive == distinct[weight][1], 'strict positive odd length')
        require(strict_negative == distinct[weight][0] - (weight == 0), 'strict negative even length')
        if weight < maximum_weight:
            require(len(nonnegative_images) == row['nonnegative'], 'nonnegative image size')
            require(len(negative_images) == row['negative'], 'negative image size')
            require(row['nonnegative'] <= ordinary[weight + 1][1], 'all-weight comparison')
            require((row['nonnegative'] == ordinary[weight + 1][1])
                    == (weight in {0, 1, 2, 3, 4, 6}), 'general-weight equality set')
            if weight >= 8 and weight % 2 == 0:
                index_n = weight // 2
                for index_j in range((index_n - 4) // 2 + 1):
                    witness = (3, 3, 3) + (2,) * (2 * index_j) + (1,) * (2 * index_n - 8 - 4 * index_j)
                    require(sum(witness) == weight + 1 and len(witness) % 2 == 1, 'strictness witness weight')
                    require(reverse_sum(witness) == 3, 'strictness witness sum')
                    require(witness not in nonnegative_images, 'strictness missing image')
            if weight >= 9 and weight % 2 == 1:
                index_n = (weight - 1) // 2
                for index_j in range((index_n - 4) // 2 + 1):
                    witness = (3, 3, 3) + (2,) * (2 * index_j) + (1,) * (2 * index_n - 7 - 4 * index_j)
                    require(sum(witness) == weight + 1 and reverse_sum(witness) == -2, 'companion witness')
                    require(witness not in negative_images, 'companion missing image')
        counts.append(row)
        endpoint_distributions.append(endpoint_distribution)
        distinct_distributions.append(distinct_distribution)
    signed_coefficients = 0
    for weight in range(maximum_weight + 1):
        for endpoint in range(-weight, weight + 1):
            convolution = sum(sum(ordinary[index])
                              * distinct_distributions[weight - 2 * index].get(endpoint, 0)
                              for index in range(weight // 2 + 1))
            require(convolution == endpoint_distributions[weight].get(endpoint, 0),
                    'all-weight signed-endpoint convolution')
            signed_coefficients += 1
    table = []
    for index_n in range(maximum_n + 1):
        source = counts[2 * index_n]['nonnegative']
        target = ordinary[2 * index_n + 1][1]
        deficit = target - source
        require((deficit == 0) == (index_n <= 3), 'equality characterization')
        if index_n >= 4:
            require(deficit >= (index_n - 4) // 2 + 1, 'quantitative deficit')
        residual_counts = [1] + [distinct[2 * index_r][1] for index_r in range(1, index_n + 1)]
        convolution = sum(residual_counts[index_r] * sum(ordinary[index_n - index_r])
                          for index_r in range(index_n + 1))
        require(source == convolution, 'distinct residual convolution')
        companion = counts[2 * index_n + 2]['negative'] - ordinary[2 * index_n + 1][0]
        require(companion >= 0, 'shifted companion comparison')
        if index_n >= 4:
            require(companion >= (index_n - 4) // 2 + 1, 'quantitative shifted companion')
        require((companion == 0) == (index_n <= 3), 'shifted companion equality set')
        table.append({'n': index_n, 'a344611': source, 'a160786': target,
                      'deficit': deficit, 'shifted_companion_deficit': companion})
    for weight in range(maximum_weight):
        require(counts[weight]['nonnegative'] <= counts[weight + 1]['nonnegative'], 'C monotonicity')
        require((counts[weight]['nonnegative'] == counts[weight + 1]['nonnegative'])
                == (weight in {0, 2, 4, 6}), 'C plateaus')
        require(counts[weight]['negative'] <= counts[weight + 1]['negative'], 'N monotonicity')
        require((counts[weight]['negative'] == counts[weight + 1]['negative'])
                == (weight in {0, 1, 3, 5, 7}), 'N plateaus')
    expected_a = [1, 2, 4, 8, 15, 27, 48, 81, 135, 220, 352, 553, 859, 1313,
                  1986, 2969, 4394, 6439, 9357, 13479, 19273]
    expected_target = [1, 2, 4, 8, 16, 29, 52, 90, 151, 248, 400, 632, 985,
                       1512, 2291, 3431, 5084, 7456, 10836, 15613, 22316]
    expected_deficit = [target - source for source, target in zip(expected_a, expected_target)]
    require([row['a344611'] for row in table[:21]] == expected_a[:len(table[:21])], 'OEIS A344611 values')
    require([row['a160786'] for row in table[:21]] == expected_target[:len(table[:21])], 'OEIS A160786 values')
    require([row['deficit'] for row in table[:21]] == expected_deficit[:len(table[:21])], 'OEIS deficit values')
    require(counts[2]['negative'] == 0 and ordinary[3][0] == 1, 'unshifted counterexample')
    require(distinct[5] == [2, 1], 'distinct negative cross-reference counterexample')
    return {
        'status': 'ALL_CHECKS_OK',
        'finite_scope': {'n': [0, maximum_n], 'partition_weights': [0, maximum_weight],
                         'map_source_weights': [0, maximum_weight - 1]},
        'checked_partitions': checked_partitions,
        'checked_nonnegative_maps': checked_nonnegative_maps,
        'checked_negative_maps': checked_negative_maps,
        'qualifying_even_weight_sources': sum(row['a344611'] for row in table),
        'common_index_reading_counterexample': {'n': 1, 'a344743': 0, 'a236914': 1,
                                                'source_status': 'right-hand argument unspecified'},
        'strict_negative_counterexample': {'weight': 5, 'negative': 2, 'a067659': 1},
        'signed_convolution_coefficients': signed_coefficients,
        'ordered_compositions': ordered_checks(),
        'counts': table,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-n', type=int, default=20)
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    if not 3 <= arguments.max_n <= 20:
        parser.error('--max-n must be between 3 and 20')
    result = run(arguments.max_n)
    encoded = json.dumps(result, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(encoded)
    print(encoded, end='')


if __name__ == '__main__':
    main()
