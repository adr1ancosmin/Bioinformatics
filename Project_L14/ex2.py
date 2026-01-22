#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Exercise 2 — Plagiarism detection using Markov models + Log-Likelihood Ratio
(DIACRITICS REMOVED)

We compare:
- Mihai Eminescu poem
- Nichita Stanescu poem

Using:
- bigram transition probabilities
- log-likelihood ratio matrix
- sliding window classification
"""

import math
import random
import re
from collections import Counter, defaultdict


# -------------------------------------------------
# INPUT TEXTS (DIACRITICS REMOVED)
# -------------------------------------------------

EMINESCU = """
Ce te legeni?

- Ce te legeni, codrule,
Fara ploaie, fara vant,
Cu crengile la pamant?
- De ce nu m-as legana,
Daca trece vremea mea!
Ziua scade, noaptea creste
Si frunzisul mi-l rareste.
Bate vantul frunza-n dunga -
Cantaretii mi-i alunga;
Bate vantul dintr-o parte -
Iarna-i ici, vara-i departe.
Si de ce sa nu ma plec,
Daca pasarile trec!
Peste varf de ramurele
Trec in stoluri randurele,
Ducand gandurile mele
Si norocul meu cu ele.
Si se duc pe rand, pe rand,
Zarea lumii-intunecand,
Si se duc ca clipele,
Scuturand aripele,
Si ma lasa pustiit,
Vestegit si amortit
Si cu doru-mi singurel,
De ma-ngan numai cu el!
""".strip()


NICHITA = """
Cu mana stanga ti-am intors spre mine chipul,
sub cortul adormitilor gutui
si de-as putea sa-mi rup din ochii tai privirea,
vazduhul serii mi-ar parea caprui.

Mi s-ar parea ca deslusessc, prin crenge,
zvelti vanatori, in arcuitii lei
din goana calului, cum isi subtie arcul.
O, tinde-ti mana stanga catre ei

si stinge tu conturul lor de lemn subtire
pe care ramurile i-au aprins,
suind sub luna-n seve caii repezi
ce-au ratacit cu timpul, pe intins.

Eu te privesc in ochi si-n jur sa sterg copacii
In ochii tai cu luna ma rasfrang
si ai putea, uitand, sa ne strivesti in gene
dar chipul ti-l intorn, pe bratul stang.
""".strip()


# -------------------------------------------------
# CONFIG
# -------------------------------------------------

ALPHA = 1.0
WINDOW_SIZE = 18
THRESHOLD = 0.25
RANDOM_SEED = 7

GENERATE_TOKENS = 120
MIX_WEIGHT = 0.50


# -------------------------------------------------
# TOKENIZATION
# -------------------------------------------------

TOKEN_RE = re.compile(r"[a-z]+")

def tokenize(text: str):
    text = text.lower()
    return TOKEN_RE.findall(text)


# -------------------------------------------------
# BIGRAM MODEL
# -------------------------------------------------

def build_bigram_counts(tokens):
    unigram = Counter()
    bigram = defaultdict(Counter)
    vocab = set(tokens)

    for i in range(len(tokens) - 1):
        w1, w2 = tokens[i], tokens[i + 1]
        unigram[w1] += 1
        bigram[w1][w2] += 1

    if tokens:
        unigram[tokens[-1]] += 1

    return unigram, bigram, vocab


def transition_prob(w1, w2, unigram, bigram, vocab, alpha=1.0):
    V = len(vocab)
    c_w1 = unigram.get(w1, 0)

    if c_w1 == 0:
        return 1.0 / V

    return (bigram[w1].get(w2, 0) + alpha) / (c_w1 + alpha * V)


def build_llr_matrix(vocab, uni_E, bi_E, uni_N, bi_N, alpha=1.0):
    llr = defaultdict(dict)

    for w1 in vocab:
        for w2 in vocab:
            pE = transition_prob(w1, w2, uni_E, bi_E, vocab, alpha)
            pN = transition_prob(w1, w2, uni_N, bi_N, vocab, alpha)
            llr[w1][w2] = math.log(pN / pE)

    return llr


# -------------------------------------------------
# SLIDING WINDOW SCORING
# -------------------------------------------------

def score_window(tokens, start, window, llr):
    end = min(len(tokens), start + window)
    if end - start < 2:
        return 0.0, 0.0

    s = 0.0
    n = 0

    for i in range(start, end - 1):
        s += llr.get(tokens[i], {}).get(tokens[i + 1], 0.0)
        n += 1

    return s, s / n


def classify(avg_llr):
    if avg_llr > THRESHOLD:
        return "N"   # Nichita
    if avg_llr < -THRESHOLD:
        return "E"   # Eminescu
    return "?"


def label_tokens(tokens, llr):
    votes = [Counter() for _ in tokens]

    for i in range(len(tokens)):
        total, avg = score_window(tokens, i, WINDOW_SIZE, llr)
        lab = classify(avg)
        for j in range(i, min(len(tokens), i + WINDOW_SIZE)):
            votes[j][lab] += 1

    labels = []
    for v in votes:
        if not v:
            labels.append("?")
        else:
            labels.append(v.most_common(1)[0][0])

    return labels


def annotate(tokens, labels):
    return " ".join(f"[{l}]{w}" for w, l in zip(tokens, labels))


# -------------------------------------------------
# TEXT GENERATION (MIXTURE MODEL)
# -------------------------------------------------

def sample_next(w1, unigram, bigram, vocab, alpha, rng):
    words = list(vocab)
    probs = [transition_prob(w1, w2, unigram, bigram, vocab, alpha) for w2 in words]

    r = rng.random() * sum(probs)
    acc = 0.0
    for w2, p in zip(words, probs):
        acc += p
        if acc >= r:
            return w2
    return words[-1]


def generate_text(start, vocab, uni_E, bi_E, uni_N, bi_N):
    rng = random.Random(RANDOM_SEED)
    tokens = [start]
    w = start

    for _ in range(GENERATE_TOKENS - 1):
        if rng.random() < MIX_WEIGHT:
            w = sample_next(w, uni_N, bi_N, vocab, ALPHA, rng)
        else:
            w = sample_next(w, uni_E, bi_E, vocab, ALPHA, rng)
        tokens.append(w)

    return tokens


# -------------------------------------------------
# MAIN
# -------------------------------------------------

def main():
    tok_E = tokenize(EMINESCU)
    tok_N = tokenize(NICHITA)

    uni_E, bi_E, vE = build_bigram_counts(tok_E)
    uni_N, bi_N, vN = build_bigram_counts(tok_N)

    vocab = vE | vN
    llr = build_llr_matrix(vocab, uni_E, bi_E, uni_N, bi_N)

    labels_E = label_tokens(tok_E, llr)

    print("\n=== ANNOTATED EMINESCU TEXT ===")
    print(annotate(tok_E, labels_E))

    start_word = random.choice(list(vocab))
    gen = generate_text(start_word, vocab, uni_E, bi_E, uni_N, bi_N)
    labels_gen = label_tokens(gen, llr)

    print("\n=== GENERATED MIXED TEXT ===")
    print(" ".join(gen))

    print("\n=== ANNOTATED GENERATED TEXT ===")
    print(annotate(gen, labels_gen))


if __name__ == "__main__":
    main()
