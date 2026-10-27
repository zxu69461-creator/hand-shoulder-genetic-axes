"""Multivariate LD score regression (observed scale) for many traits at once.

For traits i, j with per-SNP u = Z / sqrt(N):  E[u_i u_j] = S_ij * l / M + c_ij, with S the SNP-heritability / genetic
covariance matrix on the observed scale (the same S as GenomicSEM::ldsc(stand = FALSE) without prevalence conversion).
  * missing SNPs and the per-trait large-effect filter (chi2 > max(80, 0.001 N)) are 0/1 masks, so each pair uses the
    SNPs valid in both traits;
  * regression weights are separable, w_s * a_i(s) * a_j(s), with w_s = 1 / max(l_w, 1) and the heteroscedasticity
    factor a_i = 1 / (1 + N_i h2_i l / M) from a first-step h2;
  * standard errors of the genetic covariance and of rg by delete-one-block jackknife over 200 contiguous blocks (rg is
    recomputed inside every replicate).
Traits are processed in tiles so that hundreds of traits fit in memory.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_ld(ld_dir: Path):
    parts = [pd.read_csv(Path(ld_dir) / f"{c}.l2.ldscore.gz", sep="\t") for c in range(1, 23)]
    ld = pd.concat(parts, ignore_index=True).sort_values(["CHR", "BP"]).reset_index(drop=True)
    M = sum(float(open(Path(ld_dir) / f"{c}.l2.M_5_50").read().split()[0]) for c in range(1, 23))
    return ld, M


def load_traits(files, ld, hm3):
    """Z and N aligned to the HapMap 3 A1 allele; NaN where missing."""
    ref = ld[["SNP"]].merge(hm3, on="SNP", how="left")
    k = len(files); n = len(ld)
    Z = np.full((n, k), np.nan); N = np.full((n, k), np.nan)
    pos = pd.Series(np.arange(n), index=ld.SNP)
    for t, f in enumerate(files):
        d = pd.read_csv(f, sep="\t", usecols=["SNP", "A1", "A2", "Z", "N"]).dropna()
        d = d[d.SNP.isin(pos.index)].drop_duplicates("SNP")
        r = ref.set_index("SNP").loc[d.SNP]
        same = (d.A1.str.upper().values == r.A1.values) & (d.A2.str.upper().values == r.A2.values)
        flip = (d.A1.str.upper().values == r.A2.values) & (d.A2.str.upper().values == r.A1.values)
        keep = same | flip
        idx = pos.loc[d.SNP.values[keep]].values
        Z[idx, t] = np.where(flip[keep], -1.0, 1.0) * d.Z.values[keep]
        N[idx, t] = d.N.values[keep]
    return Z, N


def _sums(Ua, Ma, Ub, Mb, l, w, blocks):
    """Per-block sums for a tile pair: arrays (B, ta, tb) x 5 = C0, C1, C2, Y0, Y1."""
    B = len(blocks)
    ta, tb = Ua.shape[1], Ub.shape[1]
    out = [np.empty((B, ta, tb)) for _ in range(5)]
    for k, (b0, b1) in enumerate(blocks):
        ma, mb, ua, ub = Ma[b0:b1], Mb[b0:b1], Ua[b0:b1], Ub[b0:b1]
        ll, ww = l[b0:b1, None], w[b0:b1, None]
        wm, wu = ma * ww, ua * ww
        out[0][k] = wm.T @ mb
        out[1][k] = (wm * ll).T @ mb
        out[2][k] = (wm * ll ** 2).T @ mb
        out[3][k] = wu.T @ ub
        out[4][k] = (wu * ll).T @ ub
    return out


def _solve(C0, C1, C2, Y0, Y1):
    det = C0 * C2 - C1 ** 2
    return (C0 * Y1 - C1 * Y0) / det, (C2 * Y0 - C1 * Y1) / det


def prepare(Z, N, ld_score, w_ld, M):
    """Masks, u = Z / sqrt(N) and the step-1 heteroscedasticity factors A."""
    valid = np.isfinite(Z) & np.isfinite(N)
    thr = np.maximum(80.0, 0.001 * np.nanmax(N, axis=0))
    valid &= ~(Z ** 2 > thr[None, :])
    U = np.where(valid, Z / np.sqrt(np.where(valid, N, 1.0)), 0.0)
    Mk = valid.astype(np.float64)
    w0 = 1.0 / np.maximum(w_ld, 1.0)
    wm = Mk * w0[:, None]
    C0 = (wm * Mk).sum(0); C1 = (wm * ld_score[:, None]).sum(0); C2 = (wm * ld_score[:, None] ** 2).sum(0)
    Y0 = (wm * U ** 2).sum(0); Y1 = (wm * ld_score[:, None] * U ** 2).sum(0)
    s1, _ = _solve(C0, C1, C2, Y0, Y1)
    h2 = np.clip(s1 * M, 1e-4, None)
    Nm = np.nanmedian(np.where(valid, N, np.nan), axis=0)
    A = 1.0 / (1.0 + np.outer(ld_score, Nm * h2) / M)
    return (U * A), (Mk * A), w0, h2, valid.sum(0)


def run(Z, N, ld_score, w_ld, M, n_blocks=200, tile=200, log=print):
    """k x k arrays: gcov, gcov_se, intercept (conventional LDSC scale), rg, rg_se; h2 and h2_se on the diagonal."""
    n, k = Z.shape
    U, Mk, w0, h2_step1, nsnp = prepare(Z, N, ld_score, w_ld, M)
    Nm = np.nanmedian(np.where(np.isfinite(Z) & np.isfinite(N), N, np.nan), axis=0)
    del Z, N
    edges = np.linspace(0, n, n_blocks + 1).astype(int)
    blocks = list(zip(edges[:-1], edges[1:]))
    B = len(blocks)
    tiles = [(s, min(s + tile, k)) for s in range(0, k, tile)]
    S = np.full((k, k), np.nan); icpt = np.full((k, k), np.nan)
    Sse = np.full((k, k), np.nan); RG = np.full((k, k), np.nan); RGse = np.full((k, k), np.nan)
    diag_blk = np.empty((B, k))
    tot_diag = np.empty(k)

    def jack(parts):
        tot = [p.sum(0) for p in parts]
        sl, ic = _solve(*tot)
        return tot, sl * M, ic

    for (a0, a1) in tiles:                                  # diagonal tiles: per-block heritabilities
        parts = _sums(U[:, a0:a1], Mk[:, a0:a1], U[:, a0:a1], Mk[:, a0:a1], ld_score, w0, blocks)
        tot, Sab, ic = jack(parts)
        for b in range(B):
            sl, _ = _solve(*[t - p[b] for t, p in zip(tot, parts)])
            diag_blk[b, a0:a1] = np.diagonal(sl) * M
        tot_diag[a0:a1] = np.diagonal(Sab)
    for ia, (a0, a1) in enumerate(tiles):                   # all tile pairs
        for (b0, b1) in tiles[ia:]:
            parts = _sums(U[:, a0:a1], Mk[:, a0:a1], U[:, b0:b1], Mk[:, b0:b1], ld_score, w0, blocks)
            tot, Sab, ic = jack(parts)
            gm = np.zeros(Sab.shape); gv = np.zeros(Sab.shape)
            rm = np.zeros(Sab.shape); rv = np.zeros(Sab.shape)
            for b in range(B):
                sl, _ = _solve(*[t - p[b] for t, p in zip(tot, parts)])
                sb = sl * M
                with np.errstate(invalid='ignore', divide='ignore'):
                    rb = sb / np.sqrt(np.outer(diag_blk[b, a0:a1], diag_blk[b, b0:b1]))
                for acc_m, acc_v, x in ((gm, gv, sb), (rm, rv, rb)):
                    d = x - acc_m
                    acc_m += d / (b + 1)
                    acc_v += d * (x - acc_m)
            gse = np.sqrt((B - 1) / B * gv)
            rse = np.sqrt((B - 1) / B * rv)
            with np.errstate(invalid='ignore', divide='ignore'):
                r_full = Sab / np.sqrt(np.outer(tot_diag[a0:a1], tot_diag[b0:b1]))
            S[a0:a1, b0:b1] = Sab; Sse[a0:a1, b0:b1] = gse; icpt[a0:a1, b0:b1] = ic
            RG[a0:a1, b0:b1] = r_full; RGse[a0:a1, b0:b1] = rse
            if a0 != b0:
                S[b0:b1, a0:a1] = Sab.T; Sse[b0:b1, a0:a1] = gse.T; icpt[b0:b1, a0:a1] = ic.T
                RG[b0:b1, a0:a1] = r_full.T; RGse[b0:b1, a0:a1] = rse.T
    icpt *= np.sqrt(np.outer(Nm, Nm))
    return dict(gcov=S, gcov_se=Sse, intercept=icpt, rg=RG, rg_se=RGse, h2=np.diag(S).copy(), h2_se=np.diag(Sse).copy(),
                n_snp=nsnp)
