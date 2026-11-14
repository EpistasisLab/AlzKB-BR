import random
import pandas as pd

INPUT = 'alzkb_br_triples_memgraph_v2.csv'
OUTPUT = 'alzkb_br_permuted.csv'
OUTPUT_PKL = 'alzkb_br_permuted.pkl'
SEED = 42
SWAPS_PER_EDGE = 10
RELATIONS = [
    'BRAINREGIONOVEREXPRESSESGENE',
    'BRAINREGIONUNDEREXPRESSESGENE',
]


def degree_preserving_rewire(edge_df, rng, forbidden_pairs=None, swaps_per_edge=10):
    forbidden_pairs = forbidden_pairs or set()
    edges = [(str(s), str(t)) for s, t in edge_df[['source', 'target']].itertuples(index=False, name=None)]
    edge_set = set(edges)
    if len(edge_set) != len(edges):
        raise ValueError('Duplicate input pairs detected.')

    n = len(edges)
    target_swaps = swaps_per_edge * n
    success = 0
    attempts = 0
    max_attempts = target_swaps * 50

    while success < target_swaps and attempts < max_attempts:
        attempts += 1
        i, j = rng.sample(range(n), 2)
        a, b = edges[i]
        c, d = edges[j]
        if a == c or b == d:
            continue
        new1, new2 = (a, d), (c, b)
        if new1 in edge_set or new2 in edge_set:
            continue
        if new1 in forbidden_pairs or new2 in forbidden_pairs:
            continue
        edge_set.remove((a, b)); edge_set.remove((c, d))
        edge_set.add(new1); edge_set.add(new2)
        edges[i], edges[j] = new1, new2
        success += 1

    out = edge_df.copy()
    out['source'] = [e[0] for e in edges]
    out['target'] = [e[1] for e in edges]
    return out, success, attempts


# The uploaded .csv is tab-delimited.
df = pd.read_csv(INPUT, sep='\t', dtype=str)
original = df.copy()
rng = random.Random(SEED)
summary = []


def pairs_for(relation):
    """Current (BrainRegion, Gene) pairs held by `relation` in `df`."""
    return set(map(tuple, df.loc[df['type'].eq(relation), ['source', 'target']].to_numpy()))


# A BrainRegion must never both over- and under-express the same gene. The real
# graph contains no such pair, so the permuted graph must not create one either.
#
# Relations are rewired one at a time, so each pass needs to know which pairs are
# off-limits. It is not enough to block only the pairs of relations already
# rewired: `degree_preserving_rewire` applies `forbidden_pairs` to edges it is
# about to ADD, never to the edges the relation already holds. So if the first
# relation were rewired unconstrained, it could land on a pair the second
# relation is already sitting on, and the second pass would inherit that clash
# with no way to detect it -- it would only get cleared if that edge happened to
# be drawn for a swap, which is luck, not a guarantee.
#
# Blocking BOTH directions removes the problem at the source:
#   - pairs of relations already rewired  -> their final state
#   - pairs of relations not yet rewired  -> their original state
# The second relation then starts from a state that is already clash-free, and
# the loop's own check keeps it that way. Correct by construction, at any
# SWAPS_PER_EDGE value.
order = sorted(RELATIONS, key=lambda r: (df['type'] == r).sum(), reverse=True)

for position, rel in enumerate(order):
    mask = df['type'].eq(rel)
    part = df.loc[mask, ['source', 'target', 'type']].copy()
    head_before = part.groupby('source').size().sort_index()
    tail_before = part.groupby('target').size().sort_index()
    original_pairs = set(map(tuple, part[['source','target']].to_numpy()))

    # Every other relation's pairs: already-rewired ones are final in `df`,
    # not-yet-rewired ones are still their originals. Reading from `df` picks
    # up whichever applies without tracking the two cases separately.
    forbidden = set()
    for other in order:
        if other != rel:
            forbidden |= pairs_for(other)

    # No clash can be present before rewiring even begins.
    assert not (original_pairs & forbidden), f'{rel}: started with a clash'

    rewired, success, attempts = degree_preserving_rewire(part, rng, forbidden, SWAPS_PER_EDGE)

    # The rewire loop stops on EITHER the swap target or the attempt cap. Nothing
    # below can tell the difference: degrees, row counts and duplicate checks all
    # pass on a half-shuffled graph just as well as on a fully shuffled one. A
    # truncated shuffle leaves real biological signal in the ablation graph, making it a
    # weak null, so fail loudly rather than print "validation passed".
    assert success == SWAPS_PER_EDGE * len(part), (
        f'{rel}: only {success:,} of {SWAPS_PER_EDGE * len(part):,} swaps completed '
        f'in {attempts:,} attempts - raise max_attempts'
    )

    df.loc[mask, 'source'] = rewired['source'].to_numpy()
    df.loc[mask, 'target'] = rewired['target'].to_numpy()

    head_after = rewired.groupby('source').size().sort_index()
    tail_after = rewired.groupby('target').size().sort_index()
    final_pairs = set(map(tuple, rewired[['source','target']].to_numpy()))

    assert head_before.equals(head_after)
    assert tail_before.equals(tail_after)
    assert len(final_pairs) == len(rewired)
    assert not (final_pairs & forbidden)

    unchanged = len(original_pairs & final_pairs)
    summary.append((rel, len(part), success, attempts, unchanged, 1 - unchanged/len(part)))
    # No bookkeeping needed here: `df` was updated above, so the next pass picks
    # up this relation's final pairs via pairs_for().

assert len(df) == len(original)
assert original['type'].value_counts().sort_index().equals(df['type'].value_counts().sort_index())
untouched = ~original['type'].isin(RELATIONS)
assert original.loc[untouched].reset_index(drop=True).equals(df.loc[untouched].reset_index(drop=True))

# Final check on the whole graph: no (BrainRegion, Gene) pair carries more than
# one of the permuted relations. Holds in the original graph, must hold here.
# Note this constrains only the SAME region-gene pair. A gene over-expressed in
# one region and under-expressed in a DIFFERENT region is real and common, and
# is deliberately left untouched.
clashes = set.intersection(*(pairs_for(r) for r in RELATIONS))
assert not clashes, f'{len(clashes):,} region-gene pairs carry both relations'

# Save as a true comma-separated CSV.
df.to_csv(OUTPUT, index=False)
df.to_pickle(OUTPUT_PKL)

print(f'Saved: {OUTPUT}')
print(f'Saved: {OUTPUT_PKL}')
print(f'Total triples: {len(df):,}')
for rel, n, success, attempts, unchanged, changed_fraction in summary:
    print(f'{rel}: edges={n:,}, successful_swaps={success:,}, unchanged_original_edges={unchanged:,}, fraction_changed={changed_fraction:.4f}')
print('Validation passed: row counts, relation counts, BrainRegion degrees, Gene degrees, and untouched relations are preserved.')
