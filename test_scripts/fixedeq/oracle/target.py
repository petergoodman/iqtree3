"""IQ-TREE's empirical amino-acid frequency estimator, ported from alignment/alignment.cpp and
model/partitionmodel.cpp at 63c330d9, and the target parser of decision 006."""
import re

import numpy as np

from .likelihood import AMBIGUOUS, AMINO, UNKNOWN

NUM_STATES, STATE_UNKNOWN = 20, 23   # alignment/alignment.cpp:1051
MIN_STATE_FREQ = 1e-4                # alignment/alignment.h:21; utils/tools.cpp:7263
NUM_TIME = 8                         # alignment/alignment.cpp:4898
AMBI_AA = (4 + 8, 32 + 64, 512 + 1024)  # getAppearance, alignment.cpp:5205


def state_code(c):
    c = c.upper()
    if c in AMINO:
        return AMINO.index(c)
    if c in AMBIGUOUS:
        return NUM_STATES + "BZJ".index(c)
    if c in UNKNOWN:
        return STATE_UNKNOWN
    raise ValueError(f"{c!r} is not a protein character")


def appearance(code):
    """getAppearance (5190-5221)."""
    if code == STATE_UNKNOWN:
        return np.ones(NUM_STATES)
    app = np.zeros(NUM_STATES)
    if code < NUM_STATES:
        app[code] = 1.0
    else:
        for i in range(11):
            if AMBI_AA[code - NUM_STATES] & (1 << i):
                app[i] = 1.0
    return app


def remove_gappy(alignment):
    """removeGappySeq (723-745): drop all-unknown sequences, but keep at least three."""
    names = list(alignment)
    gap_only = {n: all(state_code(c) == STATE_UNKNOWN for c in alignment[n]) for n in names}
    kept = [n for n in names if not gap_only[n]]
    if len(kept) == len(names):
        return alignment
    if len(names) >= 3:
        for n in names:
            if len(kept) >= 3:
                break
            if gap_only[n]:
                kept.append(n)
    return {n: alignment[n] for n in names if n in kept}


def count_states(alignment, num_unknown=0):
    """countStates (4845-4884): every cell counted once, unknown cells padded by num_unknown."""
    counts = np.zeros(STATE_UNKNOWN + 1, dtype=np.int64)
    counts[STATE_UNKNOWN] = num_unknown
    for seq in alignment.values():
        for c in seq:
            counts[state_code(c)] += 1
    return counts


def convfreq(freq, min_state_freq=MIN_STATE_FREQ):
    """convfreq (5882-5909) past its keep_zero_freq return: raise entries below the minimum and
    take the deficit from the largest original entry."""
    maxi, maxfreq, total = 0, 0.0, 0.0
    for i in range(NUM_STATES):
        f = freq[i]
        if f < min_state_freq:
            freq[i] = min_state_freq
        if f > maxfreq:
            maxfreq, maxi = f, i
        total += freq[i]
    freq[maxi] += 1.0 - total
    return freq


def convert_count_to_freq(counts, keep_zero_freq=True, min_state_freq=MIN_STATE_FREQ):
    """convertCountToFreq (4886-4933): 8 fixed rounds that share ambiguous and unknown counts out
    in proportion to the current frequencies. keep_zero_freq defaults to true (tools.cpp:7262),
    so the floor applies only under --inc-zero-freq."""
    apps = [appearance(i) for i in range(STATE_UNKNOWN + 1)]
    freq = [1.0 / NUM_STATES] * NUM_STATES
    for _ in range(NUM_TIME):
        new = [0.0] * NUM_STATES
        for i in range(STATE_UNKNOWN + 1):
            if counts[i] == 0:
                continue
            nf = [freq[j] * apps[i][j] for j in range(NUM_STATES)]
            s = 0.0
            for v in nf:
                s += v
            s = 1.0 / s
            for j in range(NUM_STATES):
                new[j] += nf[j] * s * float(counts[i])
        total = 0.0
        for v in new:
            total += v
        if total == 0.0:
            break
        total = 1.0 / total
        freq = [v * total for v in new]
    freq = np.array(freq)
    return freq if keep_zero_freq else convfreq(freq, min_state_freq)


def pooled_frequencies(alignment, charsets, partition_type="p", remove_empty_seq=True):
    """The linked +F/+FO pooling of model/partitionmodel.cpp:110-160. Each partition loses its
    all-unknown sequences when it is built (alignment/superalignment.cpp:495-499); except under
    -S (TOPO_UNLINKED), its missing taxa are padded as unknown cells (132-133)."""
    ntaxa = len(alignment)
    total = np.zeros(STATE_UNKNOWN + 1, dtype=np.int64)
    for sites in charsets.values():
        part = {n: "".join(s[i] for i in sites) for n, s in alignment.items()}
        if remove_empty_seq:
            part = remove_gappy(part)
        unknown = len(sites) * (ntaxa - len(part)) if partition_type != "S" else 0
        total += count_states(part, unknown)
    return convert_count_to_freq(total)


def parse_target(text, min_state_freq=MIN_STATE_FREQ, tol=1e-6):
    """Decision 006: exactly 20 finite entries, each > 0 and >= min_state_freq, separated by
    commas, spaces or slashes; a sum within tol of 1 is normalized once. Returns (pi, raw sum)."""
    tokens = [t for t in re.split(r"[,\s/]+", text.strip()) if t]
    if len(tokens) != NUM_STATES:
        raise ValueError(f"{len(tokens)} entries; the target needs {NUM_STATES}")
    try:
        values = np.array([float(t) for t in tokens])
    except ValueError as e:
        raise ValueError(f"non-numeric entry: {e}") from None
    if not np.all(np.isfinite(values)) or np.any(values <= 0) or np.any(values < min_state_freq):
        raise ValueError(f"every entry must be finite, positive and at least {min_state_freq}")
    raw = float(values.sum())
    if abs(raw - 1.0) > tol:
        raise ValueError(f"the entries sum to {raw!r}, not within {tol} of 1")
    return values / raw, raw
