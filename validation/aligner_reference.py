"""
Writes test/reference/aligner.tsv and test/reference/kmers.tsv with results of
cutadapt's Aligner.locate and of its k-mer heuristic on random cases.

Usage (Python with cutadapt installed, e.g. in WSL):
    python validation/aligner_reference.py test/reference [ncases] [seed]
"""
import random
import sys

import cutadapt
from cutadapt._align import Aligner
from cutadapt._kmer_finder import KmerFinder
from cutadapt.kmer_heuristic import create_positions_and_kmers

outdir = sys.argv[1]
ncases = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 7)

FLAGS = [7, 11, 8, 2, 9, 6, 15, 0, 1, 4, 3, 12]
IUPAC = "ACGTRYSWKMBDHVN"


def mutate(s, rate):
    out = []
    for c in s:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice("ACGT"))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(c)
            out.append(rng.choice("ACGTN"))
        else:
            out.append(c)
    return "".join(out)


def concrete(s):
    table = dict(R="AG", Y="CT", S="CG", W="AT", K="GT", M="AC", B="CGT", D="AGT", H="ACT", V="ACG", N="ACGT")
    return "".join(rng.choice(table[c]) if c in table else c for c in s)


def rseq(n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


with open(f"{outdir}/aligner.tsv", "w") as out:
    out.write(f"# cutadapt {cutadapt.__version__}\n")
    out.write("reference\tquery\tflags\twildcard_ref\twildcard_query\tindel_cost\t"
              "max_error_rate\tmin_overlap\tresult\n")
    for _ in range(ncases):
        m = rng.randint(1, 30)
        wildcard_ref = rng.random() < 0.3
        reference = rseq(m, IUPAC if wildcard_ref else "ACGT")
        if wildcard_ref and set(reference) == {"N"}:
            reference = "A" + reference[1:]
        piece = concrete(reference)
        mode = rng.random()
        if mode < 0.3:
            piece = piece[rng.randint(0, len(piece) - 1):]
        elif mode < 0.6:
            piece = piece[: rng.randint(1, len(piece))]
        piece = mutate(piece, rng.choice([0, 0.05, 0.1, 0.2, 0.4]))
        r = rng.random()
        if r < 0.25:
            query = piece + rseq(rng.randint(0, 40))
        elif r < 0.5:
            query = rseq(rng.randint(0, 40)) + piece
        elif r < 0.9:
            query = rseq(rng.randint(0, 30)) + piece + rseq(rng.randint(0, 30))
        else:
            query = rseq(rng.randint(0, 60))
        if rng.random() < 0.1:
            query = query.lower()
        flags = rng.choice(FLAGS)
        wildcard_query = rng.random() < 0.2
        indel_cost = 1 if rng.random() < 0.8 else 100000
        rate = rng.choice([0.0, 0.1, 0.15, 0.2, 0.3, 0.5])
        min_overlap = rng.randint(1, 6)
        aligner = Aligner(reference, rate, flags=flags, wildcard_ref=wildcard_ref,
                          wildcard_query=wildcard_query, indel_cost=indel_cost,
                          min_overlap=min_overlap)
        result = aligner.locate(query)
        result = "none" if result is None else ",".join(map(str, result))
        out.write(f"{reference}\t{query}\t{flags}\t{int(wildcard_ref)}\t{int(wildcard_query)}\t"
                  f"{indel_cost}\t{rate}\t{min_overlap}\t{result}\n")

with open(f"{outdir}/kmers.tsv", "w") as out:
    out.write(f"# cutadapt {cutadapt.__version__}\n")
    out.write("adapter\tmin_overlap\terror_rate\tback\tfront\tinternal\twildcard_ref\t"
              "wildcard_query\tsequence\tpresent\n")
    for _ in range(ncases):
        m = rng.randint(3, 40)
        wildcard_ref = rng.random() < 0.3
        adapter = rseq(m, IUPAC if wildcard_ref else "ACGT")
        min_overlap = rng.randint(1, m)
        rate = rng.choice([0.0, 0.1, 0.15, 0.2, 0.3])
        back, front = rng.choice([(True, False), (False, True), (True, True)])
        internal = rng.random() < 0.7
        wildcard_query = rng.random() < 0.2
        pk = create_positions_and_kmers(adapter, min_overlap, rate, back, front, internal)
        try:
            finder = KmerFinder(pk, wildcard_ref, wildcard_query)
        except ValueError:
            finder = None
        piece = mutate(concrete(adapter), rng.choice([0, 0.1, 0.3]))
        if rng.random() < 0.5:
            piece = piece[: rng.randint(0, len(piece))]
        seq = rseq(rng.randint(0, 30)) + piece + rseq(rng.randint(0, 30))
        if rng.random() < 0.3:
            seq = rseq(rng.randint(0, 50), "ACGTN")
        present = True if finder is None else finder.kmers_present(seq)
        out.write(f"{adapter}\t{min_overlap}\t{rate}\t{int(back)}\t{int(front)}\t{int(internal)}\t"
                  f"{int(wildcard_ref)}\t{int(wildcard_query)}\t{seq}\t{int(present)}\n")
